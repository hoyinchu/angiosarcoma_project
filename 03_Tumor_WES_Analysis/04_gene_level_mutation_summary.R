library(maftools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(forcats)
library(ggsci)

setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts")

source("util_scripts/project_palettes.R")
# Load processed MAF file and clinical metadata
asc_maf_path = "data/processed/tumor_WES/ASC_mutations.maf"
asc_maf_metadata_path = "data/processed/tumor_WES/ASC_mutations_metadata.tsv"
asc_maf = read.maf(maf=asc_maf_path,clinicalData=asc_maf_metadata_path)

asc_maf@clinical.data$`Primary Site` = asc_maf@clinical.data$`Primary_Site_(Recombined)`
asc_maf@clinical.data$`Is Cutaneous` = asc_maf@clinical.data$`CUTANEOUS_AS_(EHR_EXTRACTED)`
asc_maf@clinical.data$`Is Cutaneous`[asc_maf@clinical.data$`Is Cutaneous` == 1] = "Cutaneous AS"
asc_maf@clinical.data$`Is Cutaneous`[asc_maf@clinical.data$`Is Cutaneous` == 0] = "Non-cutaneous AS"
asc_maf@clinical.data$`Is Cutaneous`[is.na(asc_maf@clinical.data$`Is Cutaneous`)] = "Unknown"
asc_maf@clinical.data$`TMB` = asc_maf@clinical.data$`TMB_all_mutations`

primary_site_palette_temp = c(
  "Other Visceral Organs"=pal_npg("nrc")(9)[1],
  "Hepatobiliary"=pal_npg("nrc")(9)[2],
  "Breast (Cutaneous)"=pal_npg("nrc")(9)[3],
  "HNFS"=pal_npg("nrc")(9)[4],
  "Extremities"=pal_npg("nrc")(9)[5],
  "Heart"=pal_npg("nrc")(9)[6],
  "Musculoskeletal"=pal_npg("nrc")(9)[7],
  "Breast (Parenchymal)"=pal_npg("nrc")(9)[8],
  #"NA"=pal_npg("nrc")(9)[9],
  "Other Rare Sites"=pal_npg("nrc")(9)[9]
)

make_oncoplot = function(maf,save_path="",top_n=15) {
  if (save_path!="") {
    pdf(file=save_path,height=8,width=12)
    #pdf(file = "Tumor_WES_Analysis/outputs/01_oncoplot_all_by_pathways.pdf",height=6)
  }
  oncoplot(maf = maf,
           top = top_n,
           clinicalFeatures=c("Is Cutaneous","Primary Site"),#,"LOCAL_RECURRENCE_(EHR_EXTRACTED)"),
           #clinicalFeatures=c("PRIMARY_SITE_(Combined)","CUTANEOUS_AS_(EHR_EXTRACTED)"),#,"LOCAL_RECURRENCE_(EHR_EXTRACTED)"),
           topBarData="TMB",
           draw_titv = TRUE,
           sortByAnnotation = TRUE,
           fontSize = 0.8,
           annotationColor = list(`Is Cutaneous`=cutaneous_palette)#,`Primary Site`=primary_site_palette_temp)
  )
  if (save_path!="") {
    dev.off()
  }
}
make_oncoplot(asc_maf,save_path = "03_Tumor_WES_Analysis/outputs/plots/04_oncoplot_all.pdf")


## Subset to specific samples
## Non-HNFS
asc_maf_non_hnfs_tsb = asc_maf@clinical.data[asc_maf@clinical.data$`Primary_Site_(Recombined)` != "HNFS"]$Tumor_Sample_Barcode
asc_maf_non_hnfs = subsetMaf(asc_maf,tsb = asc_maf_non_hnfs_tsb)
make_oncoplot(asc_maf_non_hnfs,save_path = "03_Tumor_WES_Analysis/outputs/plots/04_oncoplot_non_HNFS.pdf")
somaticInteractions(maf = asc_maf_non_hnfs, top = 25, pvalue = c(0.05, 0.1))

## HNFS
asc_maf_hnfs_tsb = asc_maf@clinical.data[asc_maf@clinical.data$`Primary_Site_(Recombined)` == "HNFS"]$Tumor_Sample_Barcode
asc_maf_hnfs = subsetMaf(asc_maf,tsb = asc_maf_hnfs_tsb)
make_oncoplot(asc_maf_hnfs,save_path = "03_Tumor_WES_Analysis/outputs/plots/04_oncoplot_HNFS.pdf")
somaticInteractions(maf = asc_maf_hnfs, top = 25, pvalue = c(0.05, 0.1))

## Cutaneous
asc_maf_cut_tsb = asc_maf@clinical.data[asc_maf@clinical.data$`Is Cutaneous` == "Cutaneous"]$Tumor_Sample_Barcode
asc_maf_cut = subsetMaf(asc_maf,tsb = asc_maf_cut_tsb)
make_oncoplot(asc_maf_cut,save_path = "03_Tumor_WES_Analysis/outputs/plots/04_oncoplot_cutaneous.pdf",top_n=30)
somaticInteractions(maf = asc_maf_cut, top = 25, pvalue = c(0.05, 0.1))

## Non-Cutaneous
asc_maf_noncut_tsb = asc_maf@clinical.data[asc_maf@clinical.data$`Is Cutaneous` == "Non-cutaneous"]$Tumor_Sample_Barcode
asc_maf_noncut = subsetMaf(asc_maf,tsb = asc_maf_noncut_tsb)
make_oncoplot(asc_maf_noncut,save_path = "03_Tumor_WES_Analysis/outputs/plots/04_oncoplot_non_cutaneous.pdf",top_n=30)
somaticInteractions(maf = asc_maf_noncut, top = 25, pvalue = c(0.05, 0.1))

## Low-TMB
asc_maf_low_tmb_tsb = asc_maf@clinical.data[asc_maf@clinical.data$TMB <= 10]$Tumor_Sample_Barcode
asc_maf_low_tmb = subsetMaf(asc_maf,tsb = asc_maf_low_tmb_tsb)
make_oncoplot(asc_maf_low_tmb,save_path = "03_Tumor_WES_Analysis/outputs/plots/04_oncoplot_sub10_tmb.pdf",top_n=30)
somaticInteractions(maf = asc_maf_low_tmb, top = 25, pvalue = c(0.05, 0.1))


## Highlight Mutsig Significant genes
mutsig_output_path = "data/processed/tumor_WES/mutsig/Apr16_2024_sig_genes.txt"
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
  geom_text_repel(data = subset(mutsig_gene_sub, (p < 0.05 & MutatedSamples>=5) | (MutatedSamples>=10) | (q < 0.5))) +
  geom_point(size=3) +
  theme_minimal() +
  geom_hline(yintercept=-log(0.05),linetype="dashed",color="gray") +
  geom_text(aes(0,-log(0.05),label = "p = 0.05", vjust = -1)) +
  scale_x_continuous(labels = scales::percent) +
  labs(color='Mutsig q-value', y = "-log(Mutsig p-value)", x="% of Samples with Mutation")
mutsig_plot
ggsave(mutsig_plot,file="03_Tumor_WES_Analysis/outputs/plots/04_MutSig_num_samples_by_significance_plot.png",width=10,height=6)
ggsave(mutsig_plot,file="03_Tumor_WES_Analysis/outputs/plots/04_MutSig_num_samples_by_significance_plot.pdf",width=10,height=6)

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


pdf("03_Tumor_WES_Analysis/outputs/plots/04_oncoplot_combined.pdf",width=10,height=8)
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
  select(c(viz_genes,"Tumor_Sample_Barcode","Is Cutaneous")) %>%      # Exclude sample names
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
  select(Gene, cutaneous_mut_prop, noncutaneous_mut_prop, fisher_pval) %>%
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
  theme(legend.title = element_blank()) 


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
primary_site_enrichment_table = make_clin_enrichment_table(asc_maf,"CUTANEOUS_AS_(EHR_EXTRACTED)")
primary_site_enrichment_table_subset = primary_site_enrichment_table[primary_site_enrichment_table$Hugo_Symbol %in% viz_genes,]
#primary_site_enrichment_table_subset = primary_site_enrichment_table[primary_site_enrichment_table$Hugo_Symbol %in% recurrent_genes$Hugo_Symbol,]
primary_site_enrichment_table_subset = primary_site_enrichment_table
primary_site_enrichment_plot = make_clin_enrichment_plot(primary_site_enrichment_table_subset,"CUTANEOUS_AS_(EHR_EXTRACTED)")
primary_site_enrichment_plot
ggsave(primary_site_enrichment_plot,file="03_Tumor_WES_Analysis/outputs/plots/04_gene_mutaton_frequency_enrichment_by_primary_sites.png",width=9)

## Perform ORA on these recurrent mutated genes (GO BP and KEGG Pathways)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

mutsig_bg_genes = mutsig_gene$gene
mutsig_recurrent_sig_genes = recurrent_genes$Hugo_Symbol
mutsig_enrich_obj = enrichGO(
  gene = mutsig_recurrent_sig_genes,
  universe = mutsig_bg_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont="BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)
mutsig_sig_genes_uniprot = bitr(mutsig_recurrent_sig_genes,fromType="SYMBOL",toType="UNIPROT",OrgDb = org.Hs.eg.db)$UNIPROT
mutsig_sig_genes_uniprot = mutsig_sig_genes_uniprot[!duplicated(mutsig_sig_genes_uniprot)]
mutsig_bg_genes_uniprot = bitr(mutsig_bg_genes,fromType="SYMBOL",toType="UNIPROT",OrgDb = org.Hs.eg.db)$UNIPROT
mutsig_bg_genes_uniprot = mutsig_bg_genes_uniprot[!duplicated(mutsig_bg_genes_uniprot)]
kegg_gson = gson_KEGG("hsa",keyType="uniprot")
mutsig_kegg_obj = enricher(
  gene=mutsig_sig_genes_uniprot,
  universe=mutsig_bg_genes_uniprot,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  gson = kegg_gson,
  minGSSize = 10,
  maxGSSize = 500
)
# MSigDB oncogenic signature enrichment (upregulated gene sets only)
msigdb_c6_gene_sets = msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::filter(grepl("up-regulated",gs_description)) %>%
  dplyr::select(gs_name, gene_symbol)
mutsig_msig_c6_obj = enricher(
  gene=mutsig_recurrent_sig_genes,
  universe=mutsig_bg_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  TERM2GENE=msigdb_c6_gene_sets,
  minGSSize = 10,
  maxGSSize = 500
)

mutsig_enrich_res = list(
  "go"=mutsig_enrich_obj,
  "kegg"=mutsig_kegg_obj,
  "msigdb_c6"=mutsig_msig_c6_obj
)

draw_enriched_bars = function(enrichment_res,enrich_type,title) {
  plot = barplot(
    enrichment_res[[enrich_type]], 
    drop = TRUE, 
    showCategory = 15, 
    title = title,
    font.size = 8
  )
  return(plot)
}

mutsig_go_bar = draw_enriched_bars(mutsig_enrich_res,"go","MutSig Recurrently Mutated Genes:\nGO-BP")
mutsig_kegg_bar = draw_enriched_bars(mutsig_enrich_res,"kegg", "MutSig Recurrently Mutated Genes:\nKEGG")
mutsig_msigdb_bar = draw_enriched_bars(mutsig_enrich_res,"msigdb_c6", "MutSig Recurrently Mutated Genes: Oncogenic Signatures")

#mutsig_bars_combined =  cowplot::plot_grid(mutsig_go_bar,mutsig_kegg_bar,mutsig_msigdb_bar,labels="AUTO",ncol=3)
# Plot the enrichment results
mutsig_bars_combined =  cowplot::plot_grid(mutsig_go_bar,mutsig_kegg_bar,labels="AUTO",ncol=2)
ggsave(mutsig_bars_combined,file="03_Tumor_WES_Analysis/outputs/plots/04_MutSig_recurrent_gene_enrichment_analysis.png",width=12)
ggsave(mutsig_bars_combined,file="03_Tumor_WES_Analysis/outputs/plots/04_MutSig_recurrent_gene_enrichment_analysis.pdf",width=12)

## Make a onehot matrix for samples carrying recurrent mutations
asc_maf_onehot = asc_maf_onehot[,recurrent_genes$Hugo_Symbol]
asc_maf_onehot = as.data.frame.matrix(asc_maf_onehot)
asc_maf_onehot = cbind(Tumor_Sample_Barcode=rownames(asc_maf_onehot),asc_maf_onehot)
write.table(asc_maf_onehot,file="data/processed/tumor_WES/ASC_mutations_MutSig_Recurrent_Gene_One_Hot.tsv",sep="\t",row.names=FALSE,quote = FALSE)

## Jaccard similarity?
jaccard_dist_matrix = as.matrix(proxy::dist(asc_maf_onehot, method = "Jaccard",convert_similarities = FALSE))
ComplexHeatmap::Heatmap(jaccard_dist_matrix)


# asc_maf_onehot$Tumor_Sample_Barcode = rownames(asc_maf_onehot)
# #asc_maf_onehot = cbind(Tumor_Sample_Barcode=rownames(asc_maf_onehot),asc_maf_onehot)
# asc_maf_onehot_merged = merge(asc_maf_onehot,asc_maf@clinical.data,by="Tumor_Sample_Barcode",all.x=TRUE)
# 
# 
# 
# asc_maf_onehot_merged$TP53 = factor(asc_maf_onehot_merged$TP53)
# asc_maf_onehot_merged$TMB_nonsyn = as.numeric(asc_maf_onehot_merged$TMB_nonsyn)
# 
# asc_maf_onehot_merged
# 
# wilcox.test(TMB_nonsyn ~ TP53, data = asc_maf_onehot_merged)
# 
# tp53_tmb_plot = ggplot(asc_maf_onehot_merged,aes(x=TP53,y=TMB_nonsyn)) +
#   geom_boxplot()
# ggsave(tp53_tmb_plot,file="03_Tumor_WES_Analysis/outputs/plots/04_tp53_tmb_plot.png",height=12,width=12)
# 
# tmb_sbs7_plot = ggplot(asc_maf_onehot_merged,aes(x=TMB_nonsyn,y=SBS7)) +
#   geom_point()
# ggsave(tmb_sbs7_plot,file="03_Tumor_WES_Analysis/outputs/plots/04_tmb_sbs7_plot.png",height=12,width=12)
# 
# asc_maf_onehot_merged$`Age_(Combined)` = as.numeric(asc_maf_onehot_merged$`Age_(Combined)`)
# age_tmb_plot = ggplot(asc_maf_onehot_merged,aes(x=`Age_(Combined)`,y=TMB_nonsyn,color=`PRIMARY_SITE_(Combined)`)) +
#   geom_point()
# ggsave(age_tmb_plot,file="03_Tumor_WES_Analysis/outputs/plots/04_age_tmb_plot.png",height=12,width=12)



# # See if certian mutations are associated with clinical attributes
# clin_enrich = clinicalEnrichment(
#   maf = asc_maf,
#   clinicalFeature = c("PRIMARY_SITE_(Combined)"),
#   minMut = 5,
# )
# clin_enrich_group_wise = clin_enrich$groupwise_comparision
# clin_enrich_group_wise_renamed = clin_enrich_group_wise
# clin_enrich_group_wise_renamed = clin_enrich_group_wise %>%
#   mutate(Group1 = case_when(
#     Group1==0 ~ "cluster 0\n(Parenchymal Breast)",
#     Group1==1 ~ "cluster 1\n(Cutaneous Breast)",
#     Group1==2 ~ "cluster 2",
#     Group1==3 ~ "cluster 3\n(HFNS)",
#     Group1==4 ~ "cluster 4",
#   ))
# clin_enrich_group_wise_renamed$p_bin = cut(
#   clin_enrich_group_wise_renamed$p_value,
#   breaks = c(0,0.01,0.05,0.1,1.1),
#   labels = c("p<0.01","p<0.05","p<0.1","p>=0.1")
#   )
# clin_enrich_group_wise_sig_genes = clin_enrich_group_wise_renamed %>% filter(p_value < 0.1 & OR > 1)
# enrich_dotplot = ggplot(clin_enrich_group_wise_sig_genes,aes(x=Group1,y=Hugo_Symbol,size=OR,color=p_bin)) +
#   geom_point() +
#   #geom_text(size=2) +
#   theme_minimal() + 
#   labs(color="Fisher p-value", x="Primary Site", y = "Genes", title="Gene Mutation Frequency Enrichment (one vs. rest)")
#   #scale_size(range=c(2,10))
# enrich_dotplot
# ggsave(enrich_dotplot,file="Tumor_WES_Analysis/outputs/plots/01_gene_mutaton_frequency_enrichment_by_seurat_clusters.png",width=12)


# # Add genes that are cancer implicated
# cosmic_table = read.csv("data/public/cosmic_cancer_gene_census.csv",sep=",",header=TRUE)
# oncokb_table = read.csv("data/public/oncokb_cancer_gene_list.tsv",sep="\t",header=TRUE,check.names = FALSE)
# oncokb_table_filtered = oncokb_table %>% filter(`OncoKB Annotated` == "Yes" & `COSMIC CGC (v99)`=="Yes")
# oncokb_filtered_gene_list = oncokb_table_filtered$`Hugo Symbol`
# reported_gene_list <- c("PTPRB","CIC","TRF1","TRF2","RAP1","TIN2","TPP1") # Add shelterin genes
# angiosarcoma_gene_list <- c(oncokb_filtered_gene_list,reported_gene_list)
# asc_oncokb_maf = subsetMaf(maf = asc_maf, genes = angiosarcoma_gene_list)
# make_oncoplot(asc_oncokb_maf,save_path = "Tumor_WES_Analysis/outputs/plots/01_oncoplot_oncoKB_only.pdf")

# Plot tmb by gene alteration frequency
# gene_mutation_by_tmb =  asc_maf@clinical.data
