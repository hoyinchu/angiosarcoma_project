library(maftools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(forcats)
library(ggsci)


source("../util_scripts/project_palettes.R")
# Load processed MAF file and clinical metadata
asc_maf_path = "../data/processed/tumor_WES/ASC_mutations.maf"
asc_maf_metadata_path = "../data/processed/tumor_WES/ASC_mutations_metadata.tsv"

asc_maf_metadata = read.csv(asc_maf_metadata_path,sep="\t")
asc_maf_metadata

## Added copy-number calls
asc_cn_call_path = "../data/processed/tumor_WES/combined_cnvkit_calls.csv"
asc_cn_df = read.csv(asc_cn_call_path)
## Keep only deep deletion or amp
asc_cn_df_filtered = asc_cn_df %>% 
  filter(CN == "Amp" | CN == "DeepDel") 

asc_maf = read.maf(maf=asc_maf_path,clinicalData=asc_maf_metadata_path,cnTable = asc_cn_df_filtered)
asc_maf@clinical.data$`PrimarySite` = asc_maf@clinical.data$`Primary_Site_(Recombined)`
asc_maf@clinical.data$`IsCutaneous` = asc_maf@clinical.data$`CUTANEOUS_AS_(EHR_EXTRACTED)`
asc_maf@clinical.data$`IsCutaneous`[asc_maf@clinical.data$`IsCutaneous` == 1] = "Cutaneous AS"
asc_maf@clinical.data$`IsCutaneous`[asc_maf@clinical.data$`IsCutaneous` == 0] = "Non-cutaneous AS"
asc_maf@clinical.data$`IsCutaneous`[is.na(asc_maf@clinical.data$`Is Cutaneous`)] = "Unknown"
asc_maf@clinical.data$`Subtype` = asc_maf@clinical.data$RAAS_LAAS_Class
asc_maf@clinical.data$`TMB` = asc_maf@clinical.data$`TMB_all_mutations`

quantile(as.numeric(asc_maf@clinical.data$TMB))

primary_site_palette_temp = c(
  "Other Visceral Organs"=pal_npg("nrc")(9)[1],
  "Hepatobiliary"=pal_npg("nrc")(9)[2],
  "Breast (Cutaneous)"=pal_npg("nrc")(9)[3],
  "HNFS"=pal_npg("nrc")(9)[4],
  "Extremities"=pal_npg("nrc")(9)[5],
  "Heart"=pal_npg("nrc")(9)[6],
  "Musculoskeletal"=pal_npg("nrc")(9)[7],
  "Breast (Parenchymal)"=pal_npg("nrc")(9)[8],
  "Other Rare Sites"=pal_npg("nrc")(9)[9]
)

## Calculate ti/tv per sample
asc_maf_titv = titv(asc_maf,useSyn = TRUE, plot = FALSE)
asc_maf_titv_proportions = asc_maf_titv$fraction.contribution
asc_maf_titv_proportions_merged = asc_maf_titv_proportions %>%
  left_join(asc_maf@clinical.data, by = "Tumor_Sample_Barcode")
ct_prop_plot = ggplot(asc_maf_titv_proportions_merged,aes(x=`C>T`,y=fct_reorder(`Primary_Site_(Recombined)`,`C>T`,.fun = mean))) + 
  geom_boxplot(width=0.5,outlier.shape = NA) +
  geom_jitter(width = 0.2,size=0.25,aes(color=`Primary_Site_(Recombined)`)) +
  pretty_plot() + L_border() +
  scale_color_manual(values=primary_site_palette_temp) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text = element_text(size = 5),legend.position = "none")
ct_prop_plot
ggsave("./outputs/plots/ct_proportion_plot.pdf",ct_prop_plot,dpi=300,width=2,height=1.6)



# Function to parse GTF specifically for gene coordinates
get_gtf_coords <- function(gtf_path) {
  # Reading only the necessary columns to save memory
  # Filter for 'gene' in the 3rd column
  gtf <- read_tsv(gtf_path, comment = "#", col_names = FALSE, 
                  col_types = "c-cdd-c-c") %>%
    filter(X3 == "gene") %>%
    dplyr::select(chr = X1, start = X4, end = X5, info = X9) %>%
    # Extract gene_name using regex from the info column
    mutate(gene_name = str_extract(info, 'gene_name "[^"]+"') %>% 
             str_replace('gene_name "', "") %>% 
             str_replace('"', "")) %>%
    dplyr::select(gene_name, chr, start) %>%
    distinct(gene_name, .keep_all = TRUE)
  return(gtf)
}

# Run this once
gene_ref = get_gtf_coords("../data/public/gencode.v19.annotation.gtf")

## Also load oncokb anotated genes
oncokb_genes = read_tsv("../data/public/oncokb_cancer_gene_list.tsv") %>% pull(`Hugo Symbol`)

library(tidyverse)
annotation_color = list(
  `PrimarySite`=primary_site_palette_maftools,
  `IsCutaneous`=cutaneous_palette,
  `Subtype`=RAAS_class_palette,
  `RAAS_LAAS_Class`=RAAS_class_palette
)
variant_palette
make_oncoplot_priority_filter = function(maf, gene_ref, priority_genes = c(), save_path="", top_n=15, top_by="total", raw_palette=FALSE) {
  
  # 1. Get Gene Summary from MAF
  gene_sum <- getGeneSummary(maf)
  gene_sum$mut_amp_del = gene_sum$MutatedSamples + gene_sum$Amp + gene_sum$DeepDel

  
  # 2. Join with GTF and define Priority
  ranked_pool <- gene_sum %>%
    inner_join(gene_ref, by = c("Hugo_Symbol" = "gene_name")) %>%
    # Create a priority flag: 0 for priority genes, 1 for others
    # This ensures priority genes come first in the arrange()
    mutate(is_priority = ifelse(Hugo_Symbol %in% priority_genes, 0, 1)) %>%
    # Primary sort by priority, secondary sort by mutation frequency
    arrange(is_priority, desc(.data[[top_by]]))
  
  print("ranked pool")
  print(head(ranked_pool,20))
  
  # 3. Greedy Filtering for 1MB proximity
  final_selection <- c()
  pool <- ranked_pool
  
  while(nrow(pool) > 0 && length(final_selection) < top_n) {
    # Select the current top (respecting priority first)
    top_gene <- pool[1, ]
    final_selection <- c(final_selection, top_gene$Hugo_Symbol)
    
    # Define the 1MB exclusion zone
    # We remove all genes within 1MB of the gene we just picked
    pool <- pool %>%
      filter(!(chr == top_gene$chr &
                 abs(start - top_gene$start) <= 1000000))
  }
  
  # 4. Generate Plot
  if (save_path != "") {
    pdf(file = save_path, height = 4, width = 6)
  }
  
  if(raw_palette) {
    oncoplot(maf = maf,
             genes = final_selection,
             clinicalFeatures = c("IsCutaneous","Primary_Site_(Wagner2024)","Subtype"),
             topBarData = "TMB",
             draw_titv = TRUE,
             sortByAnnotation = TRUE,
             #annotationColor = annotation_color,
             color = variant_palette,
             titleFontSize=0.8,
             legendFontSize=0.8,
             annotationFontSize=0.8,
             fontSize = 0.6)
  }
  else {
    oncoplot(maf = maf,
             genes = final_selection,
             #clinicalFeatures = c("IsCutaneous","PrimarySite"),
             clinicalFeatures = c("IsCutaneous","RAAS_LAAS_Class"),
             topBarData = "TMB",
             draw_titv = TRUE,
             sortByAnnotation = TRUE,
             annotationColor = annotation_color,
             color = variant_palette,
             titleFontSize=0.8,
             legendFontSize=0.8,
             annotationFontSize=0.8,
             fontSize = 0.6)
  }

  if (save_path != "") {
    dev.off()
  }
  
  return(final_selection)
}

oncoplot_all_samples_with_cn = make_oncoplot_priority_filter(
  asc_maf,gene_ref,
  priority_genes=oncokb_genes,
  save_path = "./outputs/plots/04_oncoplot_all_with_CN_new.pdf",top_by="mut_amp_del",
  top=20
)

oncoplot_all_samples_with_cn_no_priority = make_oncoplot_priority_filter(
  asc_maf,gene_ref,
  priority_genes=c(),
  save_path = "./outputs/plots/04_oncoplot_all_with_CN_new_no_priority.pdf",top_by="mut_amp_del",
  top=20
)

oncoplot_all_samples_with_cn
#make_oncoplot_priority_filter(asc_maf,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_all_with_CN.pdf",top_by="total")
#make_oncoplot_priority_filter(asc_maf,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_all_with_CN_new.pdf",top_by="total")
#make_oncoplot_gtf_filter(asc_maf,gene_ref,save_path = "./outputs/plots/04_oncoplot_all_with_CN.pdf")

## Subset by sites

## Cutaneous
asc_maf_cut_tsb = asc_maf@clinical.data[asc_maf@clinical.data$`IsCutaneous` == "Cutaneous AS"]$Tumor_Sample_Barcode
asc_maf_cut = subsetMaf(asc_maf,tsb = asc_maf_cut_tsb)
#make_oncoplot(asc_maf_cut,save_path = "./outputs/plots/04_oncoplot_cutaneous.pdf",top_n=30)
#make_oncoplot(asc_maf_cut,save_path = "./outputs/plots/04_oncoplot_cutaneous_with_CN.pdf",top_n=30)
#make_oncoplot_gtf_filter(asc_maf_cut,gene_ref,save_path = "./outputs/plots/04_oncoplot_cutaneous_with_CN.pdf",top_n=30)
make_oncoplot_priority_filter(asc_maf_cut,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_cutaneous_with_CN.pdf",top_by="mut_amp_del")
make_oncoplot_priority_filter(asc_maf_cut,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_cutaneous_with_CN_raw.pdf",top_by="mut_amp_del",raw_palette = TRUE)
make_oncoplot_priority_filter(asc_maf_cut,gene_ref,priority_genes=c(),save_path = "./outputs/plots/04_oncoplot_cutaneous_with_CN_no_priority.pdf",top_by="total",top_n = 20)

somaticInteractions(maf = asc_maf_cut, top = 25, pvalue = c(0.05, 0.1))

## Non-Cutaneous
asc_maf_noncut_tsb = asc_maf@clinical.data[asc_maf@clinical.data$`IsCutaneous` == "Non-cutaneous AS"]$Tumor_Sample_Barcode
asc_maf_noncut = subsetMaf(asc_maf,tsb = asc_maf_noncut_tsb)
#make_oncoplot(asc_maf_noncut,save_path = "./outputs/plots/04_oncoplot_non_cutaneous.pdf",top_n=30)
#make_oncoplot(asc_maf_noncut,save_path = "./outputs/plots/04_oncoplot_non_cutaneous_with_CN.pdf",top_n=30)
#make_oncoplot_gtf_filter(asc_maf_noncut,gene_ref,save_path = "./outputs/plots/04_oncoplot_non_cutaneous_with_CN.pdf",top_n=30)
make_oncoplot_priority_filter(asc_maf_noncut,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_non_cutaneous_with_CN.pdf",top_by="mut_amp_del")

somaticInteractions(maf = asc_maf_noncut, top = 25, pvalue = c(0.05, 0.1))

## Cutaneous Breast
asc_maf_cut_breast_tsb = asc_maf@clinical.data[asc_maf@clinical.data$`Primary_Site_(Recombined)` == "Breast (Cutaneous)"]$Tumor_Sample_Barcode
asc_maf_cut_breast = subsetMaf(asc_maf,tsb = asc_maf_cut_breast_tsb)
make_oncoplot_priority_filter(asc_maf_cut_breast,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_cutaneous_breast_with_CN.pdf",top_by="mut_amp_del")
make_oncoplot_priority_filter(asc_maf_cut_breast,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_cutaneous_breast_with_CN_with_subtype.pdf",top_by="mut_amp_del")

## NonCutaneous Breast
asc_maf_noncut_breast_tsb = asc_maf@clinical.data[asc_maf@clinical.data$`Primary_Site_(Recombined)` == "Breast (Parenchymal)"]$Tumor_Sample_Barcode
asc_maf_noncut_breast = subsetMaf(asc_maf,tsb = asc_maf_noncut_breast_tsb)
make_oncoplot_priority_filter(asc_maf_noncut_breast,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_non_cutaneous_breast_with_CN.pdf",top_by="mut_amp_del")
make_oncoplot_priority_filter(asc_maf_noncut_breast,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_non_cutaneous_breast_with_CN_with_subtype.pdf",top_by="mut_amp_del")

## HNFS
asc_maf_hnfs_tsb = asc_maf@clinical.data[asc_maf@clinical.data$`Primary_Site_(Recombined)` == "HNFS"]$Tumor_Sample_Barcode
asc_maf_hnfs = subsetMaf(asc_maf,tsb = asc_maf_hnfs_tsb)
make_oncoplot_priority_filter(asc_maf_hnfs,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_HNFS_with_CN.pdf",top_by="mut_amp_del")
make_oncoplot_priority_filter(asc_maf_hnfs,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_HNFS_with_CN_raw_palette.pdf",top_by="mut_amp_del",raw_palette = TRUE)
make_oncoplot_priority_filter(asc_maf_hnfs,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_HNFS_with_CN_with_subtype.pdf",top_by="mut_amp_del",raw_palette = FALSE)

## Heart
asc_maf_heart_tsb = asc_maf@clinical.data[asc_maf@clinical.data$`Primary_Site_(Recombined)` == "Heart"]$Tumor_Sample_Barcode
asc_maf_heart = subsetMaf(asc_maf,tsb = asc_maf_heart_tsb)
make_oncoplot_priority_filter(asc_maf_heart,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_Heart_with_CN.pdf",top_by="mut_amp_del")

## Extremities
asc_maf_extremities_tsb = asc_maf@clinical.data[asc_maf@clinical.data$`Primary_Site_(Recombined)` == "Extremities"]$Tumor_Sample_Barcode
asc_maf_extremities = subsetMaf(asc_maf,tsb = asc_maf_extremities_tsb)
make_oncoplot_priority_filter(asc_maf_extremities,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_Extremities_with_CN.pdf",top_by="mut_amp_del")
make_oncoplot_priority_filter(asc_maf_extremities,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_Extremities_with_CN_raw_palette.pdf",top_by="mut_amp_del",raw_palette = TRUE)

## Hepatobiliary
asc_maf_hepatobiliary_tsb = asc_maf@clinical.data[asc_maf@clinical.data$`Primary_Site_(Recombined)` == "Hepatobiliary"]$Tumor_Sample_Barcode
asc_maf_hepatobiliary = subsetMaf(asc_maf,tsb = asc_maf_hepatobiliary_tsb)
make_oncoplot_priority_filter(asc_maf_hepatobiliary,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_Hepatobiliary_with_CN.pdf",top_by="mut_amp_del")

## Musculoskeletal
asc_maf_musculoskeletal_tsb = asc_maf@clinical.data[asc_maf@clinical.data$`Primary_Site_(Recombined)` == "Musculoskeletal"]$Tumor_Sample_Barcode
asc_maf_musculoskeletal = subsetMaf(asc_maf,tsb = asc_maf_musculoskeletal_tsb)
make_oncoplot_priority_filter(asc_maf_musculoskeletal,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_Musculoskeletal_with_CN.pdf",top_by="mut_amp_del")

## Others
asc_maf_others_tsb = asc_maf@clinical.data[asc_maf@clinical.data$`Primary_Site_(Recombined)` %in% c("Musculoskeletal","Hepatobiliary","Other Visceral Organs","Other Rare Sites")]$Tumor_Sample_Barcode
asc_maf_others = subsetMaf(asc_maf,tsb = asc_maf_others_tsb)
make_oncoplot_priority_filter(asc_maf_others,gene_ref,priority_genes=oncokb_genes,save_path = "./outputs/plots/04_oncoplot_Others_with_CN.pdf",top_by="mut_amp_del",raw_palette=TRUE)

## MYC-Amp samples
myc_amp_samples = asc_cn_df %>% filter(Gene == "MYC") %>% pull(Sample_name)
myc_amp_asc_maf = subsetMaf(asc_maf,tsb = myc_amp_samples)
make_oncoplot(myc_amp_asc_maf,save_path = "./outputs/plots/04_oncoplot_MYC_amp_samples.pdf",top_n=30)

## Highlight Mutsig Significant genes
mutsig_output_path = "../data/processed/tumor_WES/mutsig2_ASC_20260202.txt"
mutsig_gene = read.csv(mutsig_output_path,sep="\t")

## Add number of mutations
gene_mutsig_merged = merge(asc_maf@gene.summary,mutsig_gene,by.x="Hugo_Symbol",by.y="gene",all.x=TRUE)
qval_threshold = 0.1
gene_mutsig_merged[["neg_log_p"]]=-log(gene_mutsig_merged$p)
gene_mutsig_merged[["q_pass_threshold"]]=gene_mutsig_merged$q < qval_threshold
mutsig_gene_sub = gene_mutsig_merged[gene_mutsig_merged$MutatedSamples>=2,] %>%
  drop_na(q_pass_threshold) %>%
  mutate(q_pass_threshold = ifelse(q_pass_threshold, "q < 0.1", "q >= 0.1"))

## add num mutated samples percentage
mutsig_gene_sub$mutated_percentage = mutsig_gene_sub$MutatedSamples / length(unique(asc_maf@data$Tumor_Sample_Barcode))

# Plot mutsig significant genes
mutsig_plot = ggplot(mutsig_gene_sub,aes(x=mutated_percentage,y=neg_log_p,label=Hugo_Symbol,color=q_pass_threshold)) +
  geom_text_repel(data = subset(mutsig_gene_sub, (p < 0.05 & MutatedSamples>=5) | (MutatedSamples>=10) | (q < 0.5)),size=2) +
  geom_point(size=1) +
  theme_minimal() +
  geom_hline(yintercept=-log(0.05),linetype="dashed",color="gray") +
  #geom_text(aes(0,-log(0.05),label = "p = 0.05", vjust = -1)) +
  scale_x_continuous(labels = scales::percent) + pretty_plot() + L_border() +
  scale_color_manual(values=c("q < 0.1"="dodgerblue3", "q >= 0.1"="gray"))
  #labs(color='Mutsig q-value', y = "-log(Mutsig p-value)", x="% of Samples with Mutation")
mutsig_plot
mutsig_plot_clean = mutsig_plot + theme(legend.position = "none", axis.title = element_blank())
mutsig_plot_clean
cowplot::ggsave2(mutsig_plot_clean,file="./outputs/plots/04_MutSig_num_samples_by_significance_plot_clean.pdf",width=1.5,height=1.5)
cowplot::ggsave2(mutsig_plot_clean,file="./outputs/plots/04_MutSig_num_samples_by_significance_plot.pdf",width=1.5,height=1.5)

## Identify recurrently mutated genes with nominal MutSig significance
asc_maf_onehot = table(asc_maf@data$Tumor_Sample_Barcode,asc_maf@data$Hugo_Symbol)
asc_maf_onehot = (asc_maf_onehot>0)*1
recurrent_genes = gene_mutsig_merged %>% filter(MutatedSamples >= 5, p < 0.05)

## Do another more curated overall plot 
## Exclude genes that do not occur at least twice in samples that do not have a known driver
mutsig_genes = c("TP53","KDR","POT1","PLCG1","ASXL1","BRAF")
literature_genes = c("PTPRB","PIK3CA","FLT4","MUC16")
known_mut_genes = c(mutsig_genes,literature_genes)
known_gene_samples = subsetMaf(asc_maf,genes = known_mut_genes)@data$Tumor_Sample_Barcode
unknown_gene_samples = setdiff(asc_maf@data$Tumor_Sample_Barcode,known_gene_samples)
unknown_gene_maf = subsetMaf(asc_maf,tsb = unknown_gene_samples)
unknown_gene_maf_summary = unknown_gene_maf@gene.summary
recurrent_unknown_genes = unknown_gene_maf_summary[unknown_gene_maf_summary$AlteredSamples >= 2]
recurrent_unknown_genes = recurrent_unknown_genes[order(-recurrent_unknown_genes$AlteredSamples),]$Hugo_Symbol
viz_genes = c(known_mut_genes,recurrent_unknown_genes)

asc_maf@clinical.data$`Primary Site` = asc_maf@clinical.data$`Primary_Site_(Recombined)`
asc_maf@clinical.data$`Is Cutaneous` = asc_maf@clinical.data$`CUTANEOUS_AS_(EHR_EXTRACTED)`
asc_maf@clinical.data$`Is Cutaneous`[asc_maf@clinical.data$`Is Cutaneous` == 1] = "Cutaneous"
asc_maf@clinical.data$`Is Cutaneous`[asc_maf@clinical.data$`Is Cutaneous` == 0] = "Non-cutaneous"
asc_maf@clinical.data$`TMB` = asc_maf@clinical.data$`TMB_all_mutations`


pdf("./outputs/plots/04_oncoplot_combined.pdf",width=10,height=8)
overall_mut_plot = oncoplot(
  asc_maf,
  genes = viz_genes,
  clinicalFeatures=c("Is Cutaneous","Primary Site"),#,"LOCAL_RECURRENCE_(EHR_EXTRACTED)"),
  #clinicalFeatures=c("PRIMARY_SITE_(Combined)","CUTANEOUS_AS_(EHR_EXTRACTED)"),#,"LOCAL_RECURRENCE_(EHR_EXTRACTED)"),
  topBarData="TMB",
  draw_titv = TRUE,
  sortByAnnotation = TRUE,
  fontSize = 0.8,
)
dev.off()


viz_genes

asc_maf_cut_noncut_barcodes = getClinicalData(asc_maf) %>%
  filter(!is.na(`IsCutaneous`)) %>%
  pull(Tumor_Sample_Barcode)
asc_maf_cut_noncut_filtered = subsetMaf(maf = asc_maf, tsb = asc_maf_cut_noncut_barcodes)
cut_noncut_enrichment = clinicalEnrichment(maf = asc_maf_cut_noncut_filtered, clinicalFeature = "IsCutaneous",minMut=5)

## Two cateogry comparison is equivalent to one versus rest
cut_noncut_enrichment_output = cut_noncut_enrichment$groupwise_comparision
cut_noncut_enrichment_output$padj = p.adjust(cut_noncut_enrichment_output$p_value,method="BH")
cut_noncut_enrichment_output = cut_noncut_enrichment_output %>%
  separate(n_mutated_group1, into = c("mut1", "total1"), sep = " of ", convert = TRUE) %>%
  separate(n_mutated_group2, into = c("mut2", "total2"), sep = " of ", convert = TRUE) %>%
  filter(Group1 == "Non-cutaneous AS") %>%
  mutate(
    pct_non_cutaneous = (mut1 / total1) * 100,
    pct_cutaneous = (mut2 / total2) * 100,
    significance = case_when(
      fdr < 0.05 ~ "FDR < 0.05",
      TRUE ~ "Not Significant"
    )
  )

cut_noncut_enrichment_output_plot = ggplot(cut_noncut_enrichment_output, aes(x = pct_non_cutaneous, y = pct_cutaneous, color = significance)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey80") +
  geom_point(size = 2) +
  geom_text_repel(
    data = cut_noncut_enrichment_output %>% filter(((padj < 0.05) & (pct_cutaneous >= 15)) | (pct_non_cutaneous >= 10) | ((pct_non_cutaneous > 7)&(pct_cutaneous > 20))), 
    aes(label = Hugo_Symbol),
    size = 2,
    max.overlaps = Inf,
    box.padding = 0.1
  ) +
  scale_color_manual(values = c(
    "FDR < 0.05" = "dodgerblue3",          # Bright Red
    "Not Significant" = "grey80"
  )) +
  pretty_plot() + L_border() +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  scale_y_continuous(labels = function(y) paste0(y, "%")) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none")

cut_noncut_enrichment_output_plot
ggsave("./outputs/plots/cut_noncut_enrichment_output_plot.pdf",cut_noncut_enrichment_output_plot,dpi=300,width=5,height=1.3)


## Plot gene by clinical attribute (Whether one gene mutation is enriched in cutaneous vs. non-cut)
asc_maf_ct = table(asc_maf@data$Tumor_Sample_Barcode,asc_maf@data$Hugo_Symbol)
asc_maf_ct_subset = as.data.frame((asc_maf_ct[,viz_genes] >= 1)) %>% add_rownames(var = "Tumor_Sample_Barcode")  #mutate_if(,as.numeric)
asc_maf_ct_subset = asc_maf_ct_subset %>% mutate(across(-Tumor_Sample_Barcode, as.integer))
asc_maf_ct_subset_merged = merge(asc_maf_ct_subset,asc_maf@clinical.data,by="Tumor_Sample_Barcode",all.x=TRUE) %>% drop_na("Is Cutaneous")

## Calculate how many percent of Cutaneous ssamples carried mutation in TP53
cut_tp53_mut_total = sum(asc_maf_ct_subset_merged[asc_maf_ct_subset_merged$`Is Cutaneous`=="Cutaneous",]$TP53)
cut_total = dim(asc_maf_ct_subset_merged[asc_maf_ct_subset_merged$`Is Cutaneous`=="Cutaneous",])[1]
binom.test(cut_tp53_mut_total,cut_total)

## Calculate how many percent of non-cutaneou ssamples carried mutation in KDR,PLCG1,TP53,PTPRB,ASXL1, PIK3CA, and POT1
asc_maf_ct_subset_merged$has_mut_in_literature_gene = as.integer(rowSums(asc_maf_ct_subset_merged[,c("KDR","PLCG1","TP53","PTPRB","ASXL1","PIK3CA","POT1")])>0)
non_cut_lit_mut_total = sum(asc_maf_ct_subset_merged[asc_maf_ct_subset_merged$`Is Cutaneous`=="Non-cutaneous",]$has_mut_in_literature_gene)
non_cut_total = dim(asc_maf_ct_subset_merged[asc_maf_ct_subset_merged$`Is Cutaneous`=="Non-cutaneous",])[1]
binom.test(non_cut_lit_mut_total,non_cut_total)



## Count cutaneous mutations
asc_maf_ct_by_cut = asc_maf_ct_subset_merged %>%
  dplyr::select(c(viz_genes,"Tumor_Sample_Barcode","Is Cutaneous")) %>%      # Exclude sample names
  group_by(`Is Cutaneous`) %>%              # Group by sample type
  summarise(across(viz_genes, sum)) %>%  # Sum mutations within each group
  ungroup() %>%
  pivot_longer(cols = -`Is Cutaneous`, names_to = "Gene", values_to = "Count")
## Make it a wide table and consider wildtype counts
asc_maf_ct_by_cut_fisher = asc_maf_ct_by_cut %>% group_by(Gene) %>%
  summarise(
    Cutaneous_Mut = sum(Count[`Is Cutaneous`=="Cutaneous"]),
    NonCutaneous_Mut = sum(Count[`Is Cutaneous`=="Non-cutaneous"]),
    Cutaneous_Wt = nrow(asc_maf_ct_subset_merged[asc_maf_ct_subset_merged$`Is Cutaneous` == "Cutaneous",]) - sum(Count[`Is Cutaneous`=="Cutaneous"]),
    NonCutaneous_Wt = nrow(asc_maf_ct_subset_merged[asc_maf_ct_subset_merged$`Is Cutaneous` == "Non-cutaneous",]) - sum(Count[`Is Cutaneous`=="Non-cutaneous"]),
  )
## Compute fisher test statistics
asc_maf_ct_by_cut_fisher = asc_maf_ct_by_cut_fisher %>% rowwise() %>%
  mutate(
    cutaneous_total = Cutaneous_Mut+Cutaneous_Wt,
    noncutaneous_total = NonCutaneous_Mut+NonCutaneous_Wt,
    cutaneous_mut_prop = Cutaneous_Mut/cutaneous_total,
    noncutaneous_mut_prop = NonCutaneous_Mut/noncutaneous_total,
    fisher_pval = {
      contingency_table <- matrix(
        c(Cutaneous_Mut, NonCutaneous_Mut, Cutaneous_Wt, NonCutaneous_Wt),
        nrow = 2
      )
      fisher.test(contingency_table)$p.value
    }
  )
## Convert to plot friendly format
fisher_plot_data = asc_maf_ct_by_cut_fisher %>%
  dplyr::select(Gene, cutaneous_mut_prop, noncutaneous_mut_prop, fisher_pval) %>%
  pivot_longer(cols = c(cutaneous_mut_prop, noncutaneous_mut_prop), 
               names_to = "Group", values_to = "Proportion") %>%
  mutate(Group = ifelse(Group == "cutaneous_mut_prop", "Cutaneous", "Non-cutaneous")) %>%
  arrange(-Proportion) %>%
  mutate(Percentage = Proportion * 100)

fisher_plot_data = fisher_plot_data %>%
  group_by(Gene) %>%
  mutate(Max_Percentage = max(Percentage)) %>%
  ungroup() %>%
  arrange(desc(Max_Percentage))

fisher_plot_segment_data = fisher_plot_data %>%
  distinct(Gene, fisher_pval, Max_Percentage) 

## Make the plot
driver_fisher_plot = ggplot(fisher_plot_data, aes(x = fct_reorder(Gene, -Max_Percentage), y = Percentage, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  
  # Add segments above the bars (only once per gene)
  geom_segment(data = fisher_plot_segment_data,
               aes(x = as.numeric(fct_reorder(Gene, -Max_Percentage)) - 0.2,
                   xend = as.numeric(fct_reorder(Gene, -Max_Percentage)) + 0.2,
                   y = Max_Percentage + 2, yend = Max_Percentage + 2),
               inherit.aes = FALSE) +
  
  # Add p-values above the segments (only once per gene)
  geom_text(data = fisher_plot_segment_data,
            aes(x = fct_reorder(Gene, -Max_Percentage),
                y = Max_Percentage + 3,
                label = paste0("p = ", signif(fisher_pval, 2))),
            inherit.aes = FALSE) +
  
  labs(#title = "Mutation Proportion by Gene in Cutaneous vs NonCutaneous Samples",
       x = "Gene", y = "Percentage Mutated %") +
  scale_fill_manual(values = c("Cutaneous" = pal_npg("nrc")(5)[1], "Non-cutaneous" = pal_npg("nrc")(5)[2])) +
  theme_minimal() +
  theme(legend.title = element_blank())  +
  pretty_plot() + L_border() + scale_y_continuous(expand = c(0,Inf))

driver_fisher_plot
ggsave(filename = "03_Tumor_WES_Analysis/outputs/plots/04_driver_gene_prop_by_cutaenous_fisher.png",driver_fisher_plot,dpi=300,width=16,height=6)
ggsave(filename = "03_Tumor_WES_Analysis/outputs/plots/04_driver_gene_prop_by_cutaenous_fisher.pdf",driver_fisher_plot,dpi=300,width=16,height=6)


# Plot genes that are enriched by a clinical attributes
make_clin_enrichment_table = function(maf_df,clin_feat,clin_feat_label) {
  clin_enrich = clinicalEnrichment(
    maf = asc_maf,
    clinicalFeature = c(clin_feat),
    minMut = 5,
  )
  clin_enrich_group_wise = clin_enrich$groupwise_comparision
  clin_enrich_group_wise_renamed = clin_enrich_group_wise
  clin_enrich_group_wise_renamed$p_bin = cut(
    clin_enrich_group_wise_renamed$p_value,
    breaks = c(0,0.01,0.05,0.1,1.1),
    labels = c("p<0.01","p<0.05","p<0.1","p>=0.1")
  )
  clin_enrich_group_wise_sig_genes = clin_enrich_group_wise_renamed %>% filter(p_value < 0.1 & OR > 1)
  return(clin_enrich_group_wise_sig_genes)
}

make_clin_enrichment_plot = function(enrichment_table,clin_feat_label) {
  enrich_dotplot = ggplot(enrichment_table,aes(x=Group1,y=Hugo_Symbol,size=OR,color=p_bin)) +
    geom_point() +
    #geom_text(size=2) +
    theme_minimal() + 
    labs(color="Fisher p-value", x=clin_feat_label, y = "Genes", title="Gene Mutation Frequency Enrichment (one vs. rest)")
  return(enrich_dotplot)
}

# Genes enriched by primary sites
#primary_site_enrichment_table = make_clin_enrichment_table(asc_maf,"CUTANEOUS_AS_(EHR_EXTRACTED)")
primary_site_enrichment_table = make_clin_enrichment_table(asc_maf,"Primary_Site_(Recombined)")

primary_site_enrichment_table_subset = primary_site_enrichment_table[primary_site_enrichment_table$Hugo_Symbol %in% viz_genes,]
#primary_site_enrichment_table_subset = primary_site_enrichment_table[primary_site_enrichment_table$Hugo_Symbol %in% recurrent_genes$Hugo_Symbol,]
primary_site_enrichment_table_subset = primary_site_enrichment_table
#primary_site_enrichment_plot = make_clin_enrichment_plot(primary_site_enrichment_table_subset,"CUTANEOUS_AS_(EHR_EXTRACTED)")
primary_site_enrichment_plot = make_clin_enrichment_plot(primary_site_enrichment_table_subset,"Primary_Site_(Recombined)")

primary_site_enrichment_plot
ggsave(primary_site_enrichment_plot,file="03_Tumor_WES_Analysis/outputs/plots/04_gene_mutaton_frequency_enrichment_by_primary_sites.png",width=9)
ggsave(primary_site_enrichment_plot,file="03_Tumor_WES_Analysis/outputs/plots/04_gene_mutaton_frequency_enrichment_by_primary_sites.png",width=9)


