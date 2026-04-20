library(Seurat)
library(cowplot)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)
library(EnhancedVolcano)
library(readr)
library(BuenColors)

library(msigdbr)
library(fgsea)

## Load MAF 
asc_maf_path = "../data/processed/tumor_WES/ASC_mutations.maf"
asc_maf_metadata_path = "../data/processed/tumor_WES/ASC_mutations_metadata.tsv"
asc_cn_call_path = "../data/processed/tumor_WES/combined_cnvkit_calls.csv"

asc_cn_df = read.csv(asc_cn_call_path)
asc_cn_df_filtered = asc_cn_df %>% filter(CN == "Amp" | CN == "DeepDel")
asc_maf = read.maf(maf=asc_maf_path,clinicalData=asc_maf_metadata_path,cnTable = asc_cn_df_filtered)

## Load SeuratObject
so = LoadSeuratRds("../data/processed/rna/ASCSeuratObj.rds")

## Look for intersecting samples
intersect_samples = intersect(so@meta.data$sample_alias,unique(asc_maf@clinical.data$sample_alias))

# 2. Create a metadata slice from MAF containing the Tumor_Sample_Barcode and the sample alias
so@meta.data$seurat_id = row.names(so@meta.data)
so_meta_short = so@meta.data[,c("seurat_id","sample_alias_cleaned")] %>% rename("sample_alias_cleaned" = "sample_alias")
maf_metadata = asc_maf@clinical.data[, .(Tumor_Sample_Barcode, sample_alias)]
mapping_df = merge(so_meta_short, maf_metadata, by = "sample_alias", all.x = FALSE)

## Subset to intersecting samples
so_subset = subset(so, subset = sample_alias %in% intersect_samples)

asc_maf_subset_tsb = unique(asc_maf@clinical.data[asc_maf@clinical.data$sample_alias %in% intersect_samples,] %>% pull(Tumor_Sample_Barcode))
asc_maf_subset = subsetMaf(asc_maf,tsb=asc_maf_subset_tsb)

## Get recurrently mutated genes
asc_subset_gene_summary = getGeneSummary(asc_maf_subset)
asc_subset_gene_summary

high_freq_genes = asc_subset_gene_summary %>% filter(AlteredSamples >= 5) %>% pull(Hugo_Symbol)

# get_gtf_coords <- function(gtf_path) {
#   # Reading only the necessary columns to save memory
#   # Filter for 'gene' in the 3rd column
#   gtf <- read_tsv(gtf_path, comment = "#", col_names = FALSE, 
#                   col_types = "c-cdd-c-c") %>%
#     filter(X3 == "gene") %>%
#     dplyr::select(chr = X1, start = X4, end = X5, info = X9) %>%
#     # Extract gene_name using regex from the info column
#     mutate(gene_name = str_extract(info, 'gene_name "[^"]+"') %>% 
#              str_replace('gene_name "', "") %>% 
#              str_replace('"', "")) %>%
#     dplyr::select(gene_name, chr, start) %>%
#     distinct(gene_name, .keep_all = TRUE)
#   return(gtf)
# }
# 
# gene_ref = get_gtf_coords("../data/public/gencode.v19.annotation.gtf")
# 
## Also load oncokb anotated genes
oncokb_genes = read_tsv("../data/public/oncokb_cancer_gene_list.tsv") %>% pull(`Hugo Symbol`)
# 
# ## Subset to genes in oncoKB
# high_freq_genes_filtered = intersect(high_freq_genes,oncokb_genes)


## Quick inspection of top mutated genes
oncoplot(maf = asc_maf_subset,
         genes = c("TIAM1","POT1","TMPRSS15"),
         #genes = high_freq_genes,
         #genes = high_freq_genes_filtered,
         clinicalFeatures = c("Primary_Site_(Recombined)", "CUTANEOUS_AS_(EHR_EXTRACTED)"),
         #topBarData = "TMB",
         draw_titv = TRUE,
         sortByAnnotation = TRUE,
         fontSize = 0.8)

## TODO: For every gene mutated at least 5 times, do mut vs wt, concat results, then check which
## genes have significant differences

# genes_to_check = c("POT1","TP53","KDR","MYC","PLCG1")
# samples_with_muts = genesToBarcodes(maf = asc_maf_subset, genes = genes_to_check, justNames = TRUE)
# mut_sample_barcodes = unique(unname(unlist(samples_with_muts)))
# asc_maf_non_driver_subset = subsetMaf(
#   maf = asc_maf_subset, 
#   tsb = setdiff(unique(asc_maf_subset@clinical.data$Tumor_Sample_Barcode), mut_sample_barcodes)
# )

# oncoplot(maf = asc_maf_non_driver_subset,
#          #genes = high_freq_genes,
#          clinicalFeatures = c("Primary_Site_(Recombined)", "CUTANEOUS_AS_(EHR_EXTRACTED)"),
#          #topBarData = "TMB",
#          draw_titv = TRUE,
#          sortByAnnotation = TRUE,
#          fontSize = 0.8)

## Add this as metadata to the seurat object
gene_mut_matrix = genesToBarcodes(maf = asc_maf_subset, genes = high_freq_genes, justNames = TRUE)
#gene_mut_matrix = genesToBarcodes(maf = asc_maf_subset, genes = high_freq_genes_filtered, justNames = TRUE)

for(gene in high_freq_genes){
  mut_samples <- gene_mut_matrix[[gene]]
  mut_seurat_ids <- mapping_df$seurat_id[mapping_df$Tumor_Sample_Barcode %in% mut_samples]
  so[[gene]] <- ifelse(colnames(so) %in% mut_seurat_ids, "Mutant", "WT")
}


library(pbapply)
run_tracked_mutation_signatures = function(seurat_obj, 
                                            high_freq_genes, 
                                            output_dir = "./signature_results",
                                            latent_var = "ESTIMATE_purity") {
  # 1. Create directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 2. Use pblapply to track progress
  message(paste("Starting analysis for", length(high_freq_genes), "genes..."))
  
  results_list <- pblapply(high_freq_genes, function(gene_name) {
    file_path <- file.path(output_dir, paste0(gene_name, "_markers.csv"))
    if (file.exists(file_path)) {
      return(read.csv(file_path))
    }
    markers <- tryCatch({
      FindMarkers(
        object = seurat_obj,
        group.by = gene_name,
        ident.1 = "Mutant",
        ident.2 = "WT",
        test.use = "LR",
        latent.vars = latent_var,
        logfc.threshold = 0,
        min.pct = 0,
        verbose = FALSE
      )
    }, error = function(e) {
      warning(paste("Error in gene", gene_name, ":", e$message))
      return(NULL)
    })
    if(!is.null(markers) && nrow(markers) > 0) {
      markers$gene_symbol <- rownames(markers)
      markers$driver_mutation <- gene_name
      write.csv(markers, file = file_path, row.names = FALSE)
      return(markers)
    }
    return(NULL)
  })
  
  # 3. Combine all non-NULL results into one table
  final_table <- bind_rows(results_list)
  # Save the master table
  write.csv(final_table, file = file.path(output_dir, "00_FINAL_combined_results.csv"), row.names = FALSE)
  return(final_table)
}

## This takes a few hours
run_tracked_mutation_signatures(so,high_freq_genes,output_dir = "./outputs/tables/by_gene")

## Read the combined csv and identify the best hit from each gene
rna_wes_combined_table = read.csv("./outputs/tables/by_gene/00_FINAL_combined_results.csv")

rna_wes_combined_table$padj_overall = p.adjust(rna_wes_combined_table$p_val,method = "BH")
rna_wes_combined_table$is_padj_overall_sig = rna_wes_combined_table$padj_overall < 0.05

rna_wes_combined_table$pair = paste0(rna_wes_combined_table$driver_mutation,":",rna_wes_combined_table$gene_symbol)
rna_wes_combined_table$gene_in_oncokb = rna_wes_combined_table$gene_symbol %in% oncokb_genes

rna_wes_combined_table_best_only = rna_wes_combined_table %>% arrange(p_val) %>% distinct(driver_mutation,.keep_all = TRUE)
rna_wes_combined_table_best_only$rank = as.numeric(rownames(rna_wes_combined_table_best_only))

## Filter to genes in oncoKB
rna_wes_combined_table_filtered = rna_wes_combined_table %>% filter(gene_in_oncokb)

rna_wes_combined_oncokb_volcano = ggplot(rna_wes_combined_table_filtered,aes(x=avg_log2FC,y=-log10(padj_overall),color=is_padj_overall_sig,label=pair)) +
  geom_point() + pretty_plot() + L_border() + scale_color_manual(values=c("TRUE"="firebrick","FALSE"="black")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept=-log10(0.05),linetype="dashed") +
  geom_vline(xintercept=0) +
  # geom_text_repel(data=rna_wes_combined_table_filtered %>% 
  #                   filter(gene_in_oncokb, is_padj_overall_sig, -log10(padj_overall) > 2.2 | (avg_log2FC > 4)),
  #                 box.padding = 0.1, size=2,max.overlaps = Inf) +
  theme(legend.position = "none",axis.title = element_blank())

rna_wes_combined_oncokb_volcano
#ggsave("./outputs/plots/rna_wes_scan_volcano.png",rna_wes_combined_oncokb_volcano,dpi=600,width=3.8,height=1.6)
ggsave("./outputs/plots/rna_wes_scan_volcano_no_text.png",rna_wes_combined_oncokb_volcano,dpi=600,width=3.8,height=1.6)


# qqnorm(rna_wes_combined_table_filtered$p_val, frame = FALSE)
# qqline(rna_wes_combined_table_filtered$p_val, col = "steelblue", lwd = 2)


rna_wes_combined_table_pot1 = rna_wes_combined_table %>% filter(driver_mutation == "POT1")

## Make volcano plot
rna_wes_combined_table_pot1$padj_significant = rna_wes_combined_table_pot1$p_val_adj < 0.05
rna_wes_combined_table_pot1$is_pot1 = rna_wes_combined_table_pot1$gene_symbol == "POT1"
rna_wes_combined_table_pot1$color_cat = ifelse(rna_wes_combined_table_pot1$is_pot1, "POT1", 
                                               ifelse(rna_wes_combined_table_pot1$padj_significant, "TRUE", "FALSE"))
rna_wes_combined_table_pot1 = rna_wes_combined_table_pot1[order(rna_wes_combined_table_pot1$color_cat == "POT1"), ]

#rna_wes_combined_table_pot1$gene = rownames(pot1_signature)
pot1_sig_volcano = ggplot(rna_wes_combined_table_pot1,aes(x=avg_log2FC,y=-log10(p_val),color=color_cat)) +
  geom_point() + pretty_plot() + L_border() + scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="black","POT1"="red")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_vline(xintercept = 0) +
  # geom_text_repel(data=rna_wes_combined_table_pot1 %>% filter(padj_significant | (-log10(p_val)>5.3)),aes(label=gene_symbol),
  #                 max.overlaps = Inf, size=2, box.padding = 0.2) +
  geom_hline(yintercept=-log10(0.05),color="gray",linetype="dashed") +
  theme(legend.position = "none",axis.title = element_blank())

pot1_sig_volcano
cowplot::ggsave2("./outputs/plots/POT1_volcano.png",pot1_sig_volcano,dpi=600,width=1.6,height=1.3)
cowplot::ggsave2("./outputs/plots/POT1_volcano_no_text.png",pot1_sig_volcano,dpi=600,width=1.6,height=1.3)

write.csv(pot1_signature,"./outputs/tables/POT1_volcano.csv")

pot1_signature

# Do pathway enrichment analysis
fgsea_hallmark_set = msigdbr(species = "Homo sapiens", category = "H") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_kegg_set = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c1_set = msigdbr(species = "Homo sapiens", category = "C1") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c5_set = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c6_set = msigdbr(species = "Homo sapiens", category = "C6") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c6_set_up_only = fgsea_c6_set[!grepl("_DN$", names(fgsea_c6_set))]
fgsea_c8_set = msigdbr(species = "Homo sapiens", category = "C8") %>% split(x = .$gene_symbol, f = .$gs_name)

## Helper function to add enrichment ratio to output
add_enriched_or = function(df, n_degs, universe_size) {
  df %>%
    mutate(
      a = overlap + 0.5, # a: In pathway, in DEG list (overlap)
      b = (n_degs - overlap) + 0.5, # b: Not in pathway, in DEG list
      c = (size - overlap) + 0.5, # c: In pathway, not in DEG list
      d = (universe_size - size - n_degs + overlap) + 0.5, # d: Not in pathway, not in DEG list
      odds_ratio = (a * d) / (b * c) # Calculate corrected Odds Ratio
    ) %>%
    dplyr::select(-a, -b, -c, -d) # Remove the temporary calculation columns
}

## Helper function to save the DEG into a csv
save_deg_csv = function(df, file_path) {
  df_to_save = as.data.frame(df)
  df_to_save[] = lapply(df_to_save, function(x) {
    if (is.list(x)) {
      return(sapply(x, function(y) paste(y, collapse = ", ")))
    } else {
      return(x)
    }
  })
  write.csv(df_to_save, file = file_path, row.names = FALSE)
}

## Declare the gene universe
fora_gene_universe = rownames(so@assays$RNA$counts)
pot1_mut_pos_degs = pot1_signature %>% filter(p_val < 0.05, avg_log2FC > 0) %>% pull(gene)

pot1_mut_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,pot1_mut_pos_degs,fora_gene_universe) %>% as.data.frame()

pot1_mut_pos_degs_fora_hallmark$padj_sig = pot1_mut_pos_degs_fora_hallmark$padj < 0.05
pot1_mut_pos_degs_fora_hallmark = add_enriched_or(pot1_mut_pos_degs_fora_hallmark, n_degs = length(pot1_mut_pos_degs), universe_size = length(fora_gene_universe))
pot1_mut_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",pot1_mut_pos_degs_fora_hallmark$pathway)
pot1_mut_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",pot1_mut_pos_degs_fora_hallmark$pretty_pathway)

pot1_mut_pos_degs_fora_hallmark_plot = ggplot(pot1_mut_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(pval),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(pot1_mut_pos_degs_fora_hallmark, padj_sig),
    aes(label = pretty_pathway), size=2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) +
  geom_hline(yintercept=-log10(0.05),color="gray",linetype="dashed")
pot1_mut_pos_degs_fora_hallmark_plot

ggsave("./outputs/plots/POT1_pos_degs_hallmark.pdf",pot1_mut_pos_degs_fora_hallmark_plot,dpi=300,width=1.6,height=1.3)

## Perform similar analysis for other highly mutated genes

## TP53
rna_wes_combined_table_tp53 = rna_wes_combined_table %>% filter(driver_mutation == "TP53")
## Make volcano plot
rna_wes_combined_table_tp53$padj_significant = rna_wes_combined_table_tp53$p_val_adj < 0.05
rna_wes_combined_table_tp53$is_tp53 = rna_wes_combined_table_tp53$gene_symbol == "TP53"
rna_wes_combined_table_tp53$color_cat = ifelse(rna_wes_combined_table_tp53$is_tp53, "TP53", 
                                              ifelse(rna_wes_combined_table_tp53$padj_significant, "TRUE", "FALSE"))
rna_wes_combined_table_tp53 = rna_wes_combined_table_tp53[order(rna_wes_combined_table_tp53$color_cat == "TP53"), ]


tp53_sig_volcano = ggplot(rna_wes_combined_table_tp53,aes(x=avg_log2FC,y=-log10(p_val),color=color_cat)) +
  geom_point() + pretty_plot() + L_border() + scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="black","TP53"="red")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_vline(xintercept = 0) +
  # geom_text_repel(data=rna_wes_combined_table_tp53 %>% filter(padj_significant,(-log10(p_val)>6), avg_log2FC > 4 | avg_log2FC < -4.5 | -log10(p_val)>7.7),aes(label=gene_symbol),
  #                 max.overlaps = Inf, size=2, box.padding = 0.2) +
  geom_hline(yintercept=-log10(0.05),color="gray",linetype="dashed") +
  theme(legend.position = "none",axis.title = element_blank())

tp53_sig_volcano
ggsave("./outputs/plots/tp53_sig_volcano.png",tp53_sig_volcano,dpi=600,width=3,height=1.3)
ggsave("./outputs/plots/tp53_sig_volcano_no_text.png",tp53_sig_volcano,dpi=600,width=3,height=1.3)


tp53_mut_pos_degs = rna_wes_combined_table_tp53 %>% filter(p_val < 0.05, avg_log2FC > 0) %>% pull(gene_symbol)
tp53_mut_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,tp53_mut_pos_degs,fora_gene_universe) %>% as.data.frame()

tp53_mut_pos_degs_fora_hallmark$padj_sig = tp53_mut_pos_degs_fora_hallmark$padj < 0.05
tp53_mut_pos_degs_fora_hallmark = add_enriched_or(tp53_mut_pos_degs_fora_hallmark, n_degs = length(tp53_mut_pos_degs), universe_size = length(fora_gene_universe))
tp53_mut_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",tp53_mut_pos_degs_fora_hallmark$pathway)
tp53_mut_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",tp53_mut_pos_degs_fora_hallmark$pretty_pathway)

tp53_mut_pos_degs_fora_hallmark_plot = ggplot(tp53_mut_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(pval),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = tp53_mut_pos_degs_fora_hallmark %>% filter(padj_sig, odds_ratio > 3.1),
    aes(label = pretty_pathway), size=2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) +
  geom_hline(yintercept=-log10(0.05),color="gray",linetype="dashed")
tp53_mut_pos_degs_fora_hallmark_plot
ggsave("./outputs/plots/tp53_mut_pos_degs_fora_hallmark_plot.pdf",tp53_mut_pos_degs_fora_hallmark_plot,dpi=300,width=3,height=1.3)


## KDR
rna_wes_combined_table_kdr = rna_wes_combined_table %>% filter(driver_mutation == "KDR")
## Make volcano plot
rna_wes_combined_table_kdr$padj_significant = rna_wes_combined_table_kdr$p_val_adj < 0.05
rna_wes_combined_table_kdr$is_kdr = rna_wes_combined_table_kdr$gene_symbol == "KDR"
rna_wes_combined_table_kdr$color_cat = ifelse(rna_wes_combined_table_kdr$is_kdr, "KDR", 
                                                ifelse(rna_wes_combined_table_kdr$padj_significant, "TRUE", "FALSE"))
rna_wes_combined_table_kdr = rna_wes_combined_table_kdr[order(rna_wes_combined_table_kdr$color_cat == "KDR"), ]

kdr_sig_volcano = ggplot(rna_wes_combined_table_kdr,aes(x=avg_log2FC,y=-log10(p_val),color=color_cat)) +
  geom_point() + pretty_plot() + L_border() + scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="black","KDR"="red")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_vline(xintercept = 0) +
  # geom_text_repel(data=rna_wes_combined_table_kdr %>% filter(padj_significant),aes(label=gene_symbol),
  #                 max.overlaps = Inf, size=2, box.padding = 0.2) +
  geom_hline(yintercept=-log10(0.05),color="gray",linetype="dashed") +
  theme(legend.position = "none",axis.title = element_blank())

kdr_sig_volcano
ggsave("./outputs/plots/kdr_sig_volcano.png",kdr_sig_volcano,dpi=600,width=3,height=1.3)
ggsave("./outputs/plots/kdr_sig_volcano_no_text.png",kdr_sig_volcano,dpi=600,width=3,height=1.3)


kdr_mut_pos_degs = rna_wes_combined_table_kdr %>% filter(p_val < 0.05, avg_log2FC > 0) %>% pull(gene_symbol)
kdr_mut_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,kdr_mut_pos_degs,fora_gene_universe) %>% as.data.frame()

kdr_mut_pos_degs_fora_hallmark$padj_sig = kdr_mut_pos_degs_fora_hallmark$padj < 0.05
kdr_mut_pos_degs_fora_hallmark = add_enriched_or(kdr_mut_pos_degs_fora_hallmark, n_degs = length(kdr_mut_pos_degs), universe_size = length(fora_gene_universe))
kdr_mut_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",kdr_mut_pos_degs_fora_hallmark$pathway)
kdr_mut_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",kdr_mut_pos_degs_fora_hallmark$pretty_pathway)

kdr_mut_pos_degs_fora_hallmark_plot = ggplot(kdr_mut_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(pval),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = kdr_mut_pos_degs_fora_hallmark %>% filter(padj_sig, odds_ratio > 3.2),
    aes(label = pretty_pathway), size=2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) +
  geom_hline(yintercept=-log10(0.05),color="gray",linetype="dashed")
kdr_mut_pos_degs_fora_hallmark_plot
ggsave("./outputs/plots/kdr_mut_pos_degs_fora_hallmark_plot.pdf",kdr_mut_pos_degs_fora_hallmark_plot,dpi=300,width=3,height=1.3)

## DO it for PLCG1
## DO it for LRP1B
rna_wes_combined_table_plcg1 = rna_wes_combined_table %>% filter(driver_mutation == "PLCG1")
## Make volcano plot
rna_wes_combined_table_plcg1$padj_significant = rna_wes_combined_table_plcg1$p_val_adj < 0.05
rna_wes_combined_table_plcg1$is_plcg1 = rna_wes_combined_table_plcg1$gene_symbol == "PLCG1"
rna_wes_combined_table_plcg1$color_cat = ifelse(rna_wes_combined_table_plcg1$is_plcg1, "PLCG1", 
                                                ifelse(rna_wes_combined_table_plcg1$padj_significant, "TRUE", "FALSE"))
rna_wes_combined_table_plcg1 = rna_wes_combined_table_plcg1[order(rna_wes_combined_table_plcg1$color_cat == "PLCG1"), ]

plcg1_sig_volcano = ggplot(rna_wes_combined_table_plcg1,aes(x=avg_log2FC,y=-log10(p_val),color=color_cat)) +
  geom_point() + pretty_plot() + L_border() + scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="black", "PLCG1"="red")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_vline(xintercept = 0) +
  # geom_text_repel(data=rna_wes_combined_table_lrp1b %>% filter(padj_significant, avg_log2FC > 4 | avg_log2FC < -4.4 | -log10(p_val) > 7.5),aes(label=gene_symbol),
  #                 max.overlaps = Inf, size=2, box.padding = 0.2) +
  geom_hline(yintercept=-log10(0.05),color="gray",linetype="dashed") +
  theme(legend.position = "none",axis.title = element_blank())

plcg1_sig_volcano
ggsave("./outputs/plots/plcg1_sig_volcano.png",plcg1_sig_volcano,dpi=600,width=3,height=1.3)
ggsave("./outputs/plots/plcg1_sig_volcano_no_text.png",plcg1_sig_volcano,dpi=600,width=3,height=1.3)


lrp1b_mut_pos_degs = rna_wes_combined_table_lrp1b %>% filter(p_val < 0.05, avg_log2FC > 0) %>% pull(gene_symbol)
lrp1b_mut_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,lrp1b_mut_pos_degs,fora_gene_universe) %>% as.data.frame()

lrp1b_mut_pos_degs_fora_hallmark$padj_sig = lrp1b_mut_pos_degs_fora_hallmark$padj < 0.05
lrp1b_mut_pos_degs_fora_hallmark = add_enriched_or(lrp1b_mut_pos_degs_fora_hallmark, n_degs = length(lrp1b_mut_pos_degs), universe_size = length(fora_gene_universe))
lrp1b_mut_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",lrp1b_mut_pos_degs_fora_hallmark$pathway)
lrp1b_mut_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",lrp1b_mut_pos_degs_fora_hallmark$pretty_pathway)

lrp1b_mut_pos_degs_fora_hallmark_plot = ggplot(lrp1b_mut_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(pval),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = lrp1b_mut_pos_degs_fora_hallmark %>% filter(padj_sig, odds_ratio > 2.6),
    aes(label = pretty_pathway), size=2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) +
  geom_hline(yintercept=-log10(0.05),color="gray",linetype="dashed")
lrp1b_mut_pos_degs_fora_hallmark_plot
ggsave("./outputs/plots/lrp1b_mut_pos_degs_fora_hallmark_plot.pdf",lrp1b_mut_pos_degs_fora_hallmark_plot,dpi=300,width=3,height=1.3)



## DO it for LRP1B
rna_wes_combined_table_lrp1b = rna_wes_combined_table %>% filter(driver_mutation == "LRP1B")
## Make volcano plot
rna_wes_combined_table_lrp1b$padj_significant = rna_wes_combined_table_lrp1b$p_val_adj < 0.05
rna_wes_combined_table_lrp1b$is_lrp1b = rna_wes_combined_table_lrp1b$gene_symbol == "LRP1B"
rna_wes_combined_table_lrp1b$color_cat = ifelse(rna_wes_combined_table_lrp1b$is_lrp1b, "LRP1B", 
                                              ifelse(rna_wes_combined_table_lrp1b$padj_significant, "TRUE", "FALSE"))
rna_wes_combined_table_lrp1b = rna_wes_combined_table_lrp1b[order(rna_wes_combined_table_lrp1b$color_cat == "LRP1B"), ]

lrp1b_sig_volcano = ggplot(rna_wes_combined_table_lrp1b,aes(x=avg_log2FC,y=-log10(p_val),color=color_cat)) +
  geom_point() + pretty_plot() + L_border() + scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="black", "LRP1B"="red")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_vline(xintercept = 0) +
  # geom_text_repel(data=rna_wes_combined_table_lrp1b %>% filter(padj_significant, avg_log2FC > 4 | avg_log2FC < -4.4 | -log10(p_val) > 7.5),aes(label=gene_symbol),
  #                 max.overlaps = Inf, size=2, box.padding = 0.2) +
  geom_hline(yintercept=-log10(0.05),color="gray",linetype="dashed") +
  theme(legend.position = "none",axis.title = element_blank())

lrp1b_sig_volcano
ggsave("./outputs/plots/lrp1b_sig_volcano.png",lrp1b_sig_volcano,dpi=600,width=3,height=1.3)
ggsave("./outputs/plots/lrp1b_sig_volcano_no_text.png",lrp1b_sig_volcano,dpi=600,width=3,height=1.3)


lrp1b_mut_pos_degs = rna_wes_combined_table_lrp1b %>% filter(p_val < 0.05, avg_log2FC > 0) %>% pull(gene_symbol)
lrp1b_mut_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,lrp1b_mut_pos_degs,fora_gene_universe) %>% as.data.frame()

lrp1b_mut_pos_degs_fora_hallmark$padj_sig = lrp1b_mut_pos_degs_fora_hallmark$padj < 0.05
lrp1b_mut_pos_degs_fora_hallmark = add_enriched_or(lrp1b_mut_pos_degs_fora_hallmark, n_degs = length(lrp1b_mut_pos_degs), universe_size = length(fora_gene_universe))
lrp1b_mut_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",lrp1b_mut_pos_degs_fora_hallmark$pathway)
lrp1b_mut_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",lrp1b_mut_pos_degs_fora_hallmark$pretty_pathway)

lrp1b_mut_pos_degs_fora_hallmark_plot = ggplot(lrp1b_mut_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(pval),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = lrp1b_mut_pos_degs_fora_hallmark %>% filter(padj_sig, odds_ratio > 2.6),
    aes(label = pretty_pathway), size=2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) +
  geom_hline(yintercept=-log10(0.05),color="gray",linetype="dashed")
lrp1b_mut_pos_degs_fora_hallmark_plot
ggsave("./outputs/plots/lrp1b_mut_pos_degs_fora_hallmark_plot.pdf",lrp1b_mut_pos_degs_fora_hallmark_plot,dpi=300,width=3,height=1.3)

## DO it for MYC
rna_wes_combined_table_myc = rna_wes_combined_table %>% filter(driver_mutation == "MYC")
## Make volcano plot
rna_wes_combined_table_myc$padj_significant = rna_wes_combined_table_myc$p_val_adj < 0.05
rna_wes_combined_table_myc$is_myc = rna_wes_combined_table_myc$gene_symbol == "MYC"
rna_wes_combined_table_myc$color_cat = ifelse(rna_wes_combined_table_myc$is_myc, "MYC", 
                                              ifelse(rna_wes_combined_table_myc$padj_significant, "TRUE", "FALSE"))
rna_wes_combined_table_myc = rna_wes_combined_table_myc[order(rna_wes_combined_table_myc$color_cat == "MYC"), ]
myc_sig_volcano = ggplot(rna_wes_combined_table_myc,aes(x=avg_log2FC,y=-log10(p_val),color=color_cat)) +
  geom_point() + pretty_plot() + L_border() + scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="black","MYC"="red")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_vline(xintercept = 0) +
  # geom_text_repel(data=rna_wes_combined_table_myc %>% filter(-log10(p_val)> 3.5),aes(label=gene_symbol),
  #                 max.overlaps = Inf, size=2, box.padding = 0.2) +
  geom_hline(yintercept=-log10(0.05),color="gray",linetype="dashed") +
  theme(legend.position = "none",axis.title = element_blank())

myc_sig_volcano
ggsave("./outputs/plots/myc_sig_volcano.png",myc_sig_volcano,dpi=600,width=3,height=1.3)
ggsave("./outputs/plots/myc_sig_volcano_no_text.png",myc_sig_volcano,dpi=600,width=3,height=1.3)


myc_mut_pos_degs = rna_wes_combined_table_myc %>% filter(p_val < 0.05, avg_log2FC > 0) %>% pull(gene_symbol)
myc_mut_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,myc_mut_pos_degs,fora_gene_universe) %>% as.data.frame()

myc_mut_pos_degs_fora_hallmark$padj_sig = myc_mut_pos_degs_fora_hallmark$padj < 0.05
myc_mut_pos_degs_fora_hallmark = add_enriched_or(myc_mut_pos_degs_fora_hallmark, n_degs = length(myc_mut_pos_degs), universe_size = length(fora_gene_universe))
myc_mut_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",myc_mut_pos_degs_fora_hallmark$pathway)
myc_mut_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",myc_mut_pos_degs_fora_hallmark$pretty_pathway)

myc_mut_pos_degs_fora_hallmark_plot = ggplot(myc_mut_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(pval),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = myc_mut_pos_degs_fora_hallmark %>% filter(padj_sig, odds_ratio > 5),
    aes(label = pretty_pathway), size=2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) +
  geom_hline(yintercept=-log10(0.05),color="gray",linetype="dashed")
myc_mut_pos_degs_fora_hallmark_plot
ggsave("./outputs/plots/myc_mut_pos_degs_fora_hallmark_plot.pdf",myc_mut_pos_degs_fora_hallmark_plot,dpi=300,width=3,height=1.3)

### Save the outputs

rna_wes_combined_table_subset = rna_wes_combined_table %>% filter(
  driver_mutation %in% c("POT1","TP53","LRP1B","KDR","PLCG1","MYC")
)
write.csv(rna_wes_combined_table_subset,"./outputs/tables/rna_wes_combined_table_subset.csv",row.names = FALSE)

key_enriched_pathways_fora = rbind(
  pot1_mut_pos_degs_fora_hallmark %>% mutate(driver="POT1"),
  kdr_mut_pos_degs_fora_hallmark %>% mutate(driver="KDR"),
  tp53_mut_pos_degs_fora_hallmark %>% mutate(driver="TP53"),
  plcg1_mut_pos_degs_fora_hallmark %>% mutate(driver="PLCG1"),
  lrp1b_mut_pos_degs_fora_hallmark %>% mutate(driver="LRP1B"),
  myc_mut_pos_degs_fora_hallmark %>% mutate(driver="MYC")
)

save_deg_csv(key_enriched_pathways_fora,"./outputs/tables/key_enriched_pathways_fora.csv")




# 
# 
# 
# 
# # ## Adding metadata
# # mutation_cols = c("Total_Mutations","TMB_all_mutations","Total_Nonsyn_Mutations","TMB_nonsyn")
# # sbs_cols = names(tumor_WES_meta)[grepl("SBS",names(tumor_WES_meta))]
# # one_hot_cols = names(tumor_WES_onehot_data)
# # mutation_subset_cols_to_keep = c("Tumor_Sample_Barcode","sample_alias",mutation_cols,sbs_cols)
# # tumor_WES_meta_subset = tumor_WES_meta[,mutation_subset_cols_to_keep]
# 
# ## Add one-hot somatic mutations
# tumor_WES_onehot_data = read.csv("../data/processed/tumor_WES/ASC_mutations_MutSig_Recurrent_Gene_One_Hot.tsv",sep="\t",check.names = FALSE)
# ## Change column names to disambiguate
# colnames(tumor_WES_onehot_data) = paste(colnames(tumor_WES_onehot_data),"mutation", sep = "_")
# names(tumor_WES_onehot_data)[names(tumor_WES_onehot_data) == 'Tumor_Sample_Barcode_mutation'] = 'Tumor_Sample_Barcode'
# tumor_WES_meta_subset_merged = merge(tumor_WES_meta_subset,tumor_WES_onehot_data,on="Tumor_Sample_Barcode",all.x=TRUE)
# 
# ## Add to Seurat Metadata
# seurat_meta_expanded = merge(so_subset@meta.data,tumor_WES_meta_subset_merged,by="sample_alias",all.x=TRUE)
# seurat_meta_expanded = seurat_meta_expanded %>% select("og_id", everything())
# so_subset@meta.data = seurat_meta_expanded
# rownames(so_subset@meta.data) = so_subset@meta.data$og_id
# 
# ## Re-scale the data w.r.t. to the rest of the subset samples
# vst_rescaled = so_subset@assays$RNA$counts
# vst_rescaled = as.data.frame(t(scale(t(vst_rescaled))),check.names=FALSE)
# so_subset@assays$RNA$vst_rescaled = vst_rescaled
# 
# FeaturePlot(so,features = c("rna_PLCG1"))
# 
# ## Check SBS by cluster
# VlnPlot(so_subset,features = c("SBS6","SBS7","SBS15","SBS87"))
# 
# ## Check Shelterin-related gene expression in POT1-mutant samples
# telomere_related_geneset = c(
#   "POT1","TERT","ATRX","DAXX",
#   "TERF1","TERF2","TINF2","RPA",
#   "TPP1","WRN","ATM","ATR",
#   "BRCA1","BRCA2","RAD51","BLM","BARD1",
# )
# so_subset = SetIdent(so_subset, value = so_subset@meta.data$POT1_mutation)
# pot1_mutant_heatmap = DoHeatmap(so_subset, features = telomere_related_geneset, 
#                                 disp.max = 3.5,disp.min = -3.5,slot="vst_rescaled")
# pot1_mutant_heatmap
# ggsave(pot1_mutant_heatmap,filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_POT1_mutant_shelterin_expression_plot.png",dpi=300)
# 
# ## Check telomere-related gene expression in POT1-mutant samples in parenchymal breast samples only
# so_site_subset = subset(so_subset, idents = "BREAST (PARENCHYMAL)")
# so_site_subset = SetIdent(so_site_subset, value = so_site_subset@meta.data$POT1_mutation)
# pot1_mutant_parenchymal_breast_heatmap = DoHeatmap(so_site_subset, features = telomere_related_geneset,disp.max = 3.5,disp.min = -3.5,slot="vst_rescaled")
# pot1_mutant_parenchymal_breast_heatmap
# ggsave(pot1_mutant_parenchymal_breast_heatmap,filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_POT1_mutant_telomere_related_expression_in_parenchymal_breast_plot.png",dpi=300)
# 
# 
# ## Check Shelterin-related gene expression in KDR-mutant samples
# so_subset = SetIdent(so_subset, value = so_subset@meta.data$KDR_mutation)
# kdr_mutant_heatmap = DoHeatmap(so_subset, features = telomere_related_geneset,
#                                disp.max = 3.5,disp.min = -3.5,slot="vst_rescaled")
# kdr_mutant_heatmap
# ggsave(kdr_mutant_heatmap,filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_KDR_mutant_shelterin_expression_plot.png",dpi=300)
# 
# ## Check general gene markers per primary sites
# so_subset = SetIdent(so_subset, value = so_subset@meta.data$`PRIMARY SITE (Combined)`)
# primary_site_telomere_genes_heatmap = DoHeatmap(so_subset, features = telomere_related_geneset,
#                                disp.max = 3.5,disp.min = -3.5,slot="vst_rescaled"
# )
# primary_site_telomere_genes_heatmap
# ggsave(primary_site_telomere_genes_heatmap,
#        filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_telomere_related_genes_expression_plot_by_primary_sites.png",dpi=300,height=10,width=20)
# 
# so_subset@meta.data[so_subset@meta.data$FLT4_mutation==1,]
# 
# # so_subset = SetIdent(so_subset, value = so_subset@meta.data$`PRIMARY SITE (Combined)`)
# # breast_kdr_diff_expr = FindMarkers(so_subset, ident.1 = 1, group.by = 'POT1_mutation', subset.ident = "BREAST (PARENCHYMAL)")
# # so_site_subset = subset(so_subset, idents = "BREAST (PARENCHYMAL)")
# # so_site_subset = SetIdent(so_site_subset, value = so_site_subset@meta.data$POT1_mutation)
# # pos_enrich_kdr_genes = breast_kdr_diff_expr#[(breast_kdr_diff_expr$p_val<0.05) & (breast_kdr_diff_expr$avg_log2FC>0),]
# #DoHeatmap(so_site_subset, features = rownames(pos_enrich_kdr_genes)[1:50],disp.max = 3.5,disp.min = -3.5,slot="vst_rescaled")
# #DoHeatmap(so_site_subset, features = c("TERT","ATRX","DAXX"),disp.max = 3.5,disp.min = -3.5,slot="vst_rescaled")
# 
# pot1_subset = subset(so_subset, subset = POT1 < 1,slot="vst_rescaled")
# 
# rna_mutation_overlay = DimPlot(so_subset, reduction = "pca",group.by=c(
#   "seurat_clusters","PRIMARY SITE (Combined)",
#   "TP53_mutation","KDR_mutation","POT1_mutation",
#   "FLT4_mutation","ASXL1_mutation","ATRX_mutation",
#   "PTPRO_mutation","PLCG1_mutation","PTPRB_mutation",
#   "SETBP1_mutation","ARID1A_mutation","LRP2_mutation"
#   ),ncol=5,pt.size=4)
# 
# rna_mutation_overlay
# 
# ## Check the correlation between gene expression and TMB
# rna_vst_expression = so_subset@assays$RNA$counts
# tmb_nonsyn = so_subset@meta.data$TMB_nonsyn
# names(tmb_nonsyn) = rownames(so_subset@meta.data)
# 
# calc_correlation = function(expr_mat,feature_vec) {
#   gene_feature_correlations = as.data.frame(t(apply(expr_mat, 1, function(gene_data) {
#     test_res = cor.test(gene_data, feature_vec, use = "complete.obs", method = "pearson")
#     parital_res = c(
#       gene = rownames(gene_data)[[1]],
#       correlation = test_res$estimate,
#       pvalue = test_res$p.value,
#       conf_low = test_res$conf.int[1],
#       conf_high = test_res$conf.int[2]
#     )
#     return(parital_res)
#   })),check.names=FALSE)
#   gene_feature_correlations_ordered = gene_feature_correlations[order(gene_feature_correlations$pvalue),]
#   gene_feature_correlations$neg_logp = -log(gene_feature_correlations$pvalue)
#   gene_feature_correlations$neg_logp_capped = pmin(gene_feature_correlations$neg_logp,20)
#   gene_feature_correlations$gene = rownames(gene_feature_correlations)
#   return(gene_feature_correlations)
# }
# 
# plot_correlation_volcano = function(corr_df) {
#   correlation_volcano = ggplot(corr_df,aes(x=correlation.cor,y=neg_logp_capped,label=gene)) +
#     geom_text_repel(data = subset(corr_df, pvalue < 0.05), max.overlaps=50) +
#     geom_point() +
#     geom_hline(yintercept = -log(0.05)) +
#     geom_vline(xintercept = 0) +
#     theme_minimal() +
#     labs(x="TMB vs. Gene Expression Pearson Correlation",y="-log(p-value)")
#   return(correlation_volcano)
# }
# 
# tmb_gene_correlations = calc_correlation(rna_vst_expression,tmb_nonsyn)
# tmb_correlation_volcano = plot_correlation_volcano(tmb_gene_correlations)
# tmb_correlation_volcano
# ggsave(correlation_volcano,filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_tmb_rna_expression_correlation_volcano.png",dpi=300)
# 
# ## Remove those with TMB > 10 and redo
# tmb_nonsyn_no_hyper = tmb_nonsyn[tmb_nonsyn<10]
# rna_vst_expression_no_hyper = rna_vst_expression[,colnames(rna_vst_expression) %in% names(tmb_nonsyn_no_hyper)]
# tmb_gene_correlations_no_hyper = calc_correlation(rna_vst_expression_no_hyper,tmb_nonsyn_no_hyper)
# correlation_no_hyper_volcano = plot_correlation_volcano(tmb_gene_correlations_no_hyper)
# 
# ## Restrict to same site only
# hnfs_only = rownames(so_subset@meta.data[so_subset@meta.data$`PRIMARY SITE (Combined)` == "HNFS",])
# tmb_nonsyn_hfns_only = tmb_nonsyn[hnfs_only]
# rna_vst_expression_hfns_only = rna_vst_expression[,colnames(rna_vst_expression) %in% names(tmb_nonsyn_hfns_only)]
# tmb_gene_correlations_hfns_only = calc_correlation(rna_vst_expression_hfns_only,tmb_nonsyn_hfns_only)
# correlation_hfns_volcano = plot_correlation_volcano(tmb_gene_correlations_hfns_only)
# correlation_hfns_volcano
# ggsave(correlation_hfns_volcano,filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_tmb_rna_expression_correlation_hfns_only_volcano.png",dpi=300)
# 
# ## Make POT1-TPP1 plots
# so_subset = SetIdent(so_subset, value = so_subset@meta.data$POT1_mutation)
# pot1_tpp1_scatter_by_mut = FeatureScatter(so_subset,feature1 = "POT1",feature2 = "TPP1")
# ggsave(pot1_tpp1_scatter_by_mut,filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_pot1_vs_tpp1_scatter_by_pot1_mutation.png",dpi=300)
# 
# pot1_tpp1_violin = VlnPlot(so_subset,features = c("TPP1"))
# ggsave(pot1_tpp1_violin,filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_pot1_vs_tpp1_violin.png",dpi=300)
# 
# so_subset = SetIdent(so_subset, value = so_subset@meta.data$seurat_clusters)
# pot1_tpp1_scatter_by_cluster = FeatureScatter(so_subset,feature1 = "POT1",feature2 = "TPP1")
# ggsave(pot1_tpp1_scatter_by_cluster,filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_pot1_vs_tpp1_scatter_by_cluster.png",dpi=300)
# 
# ## Check correlation between SBS7 and gene expression?
# sbs7_expr = so_subset@meta.data$SBS1
# sbs7_gene_correlations = calc_correlation(rna_vst_expression,sbs7_expr)
# sbs7_gene_correlations_volcano = plot_correlation_volcano(sbs7_gene_correlations)
# sbs7_gene_correlations_volcano

##
# ir_genes = read.csv("data/public/ionizing_radiation_up_genes.csv")
# so_subset = AddModuleScore(so_subset,features = list(c(ir_genes$gene)),name="ir_signature")
# FeaturePlot(so_subset,reduction = "pca",features = "ir_signature1",pt.size=4)
# so_subset = SetIdent(so_subset, value = so_subset@meta.data$seurat_clusters)
# VlnPlot(so_subset,features = "ir_signature1")
# FeatureScatter(so_subset,feature1 = "FLT4",feature2 = "FLT1")


# ## Check if Cluster 2&3 samples are separated by keratin expression
# so_subset = UpdateSlots(so_subset)
# colnames(so_subset) = Cells(so_subset[["RNA_snn"]])
# so_subset = SetIdent(so_subset, value = so_subset@meta.data$seurat_clusters)
# so_subset_cluster_2_and_3 = subset(so_subset, subset = seurat_clusters == 2)
# 
# 
# 
# VlnPlot(so_subset_cluster_2_and_3,features = c("KRT1"))
# FeaturePlot(so_subset_cluster_2_and_3,reduction = "pca",features = c("TMB_nonsyn"))
# FeatureScatter(so_subset_cluster_2_and_3,feature1 = "TMB_nonsyn",feature2 = "KRT4")

