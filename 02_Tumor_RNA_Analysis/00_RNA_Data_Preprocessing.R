library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(DESeq2)

## Load data
setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts")


count_data = read.csv("data/raw/rna/All_Oct2023samples_159bxs.gene_reads.gct",sep="\t",skip=2,check.names = FALSE)
tpm_data = read.csv("data/raw/rna/All_Oct2023samples_159bxs.gene_tpm.gct",sep="\t",skip=2,check.names = FALSE)

## Perform filtering based on sequencing metrics
## Script based on QC RNA metrics provided by Jorge
sample_qc_summary = function(metrics_df,count_data) {
  # Flag samples with contamination ≥ 5%:
  metrics_df = metrics_df %>% mutate(CONTAM_HARD_FLAG = ifelse(pct_contamination >= 0.05, 'FAIL', 'PASS'))
  
  # Samples that meet both of these criteria should also be flagged:
  # - Reads aligned in pairs < 25 million
  # - % mRNA bases < 60%
  metrics_df = metrics_df %>% mutate(RAP_mRNA_HARD_FLAG = ifelse(reads_aligned_in_pairs < 25e6 & rna_pct_mrna_bases < 0.6, 'FAIL', 'PASS'))

  #Also flag outliers in chimeras %:
  metrics_df = metrics_df %>% mutate(CHIMERAS_HARD_FLAG = ifelse(pct_chimerism > 0.05, 'FAIL', 'PASS'))
  
  # Flag samples with < 8000 genes mapped
  genes_mapped = colSums(count_data > 0)
  valid_detected_samples = names(genes_mapped[genes_mapped >= 8000])
  metrics_df$GENES_DETECTED_FLAG = ifelse(metrics_df$`entity:sample_id` %in% valid_detected_samples, "PASS","FAIL")

  flag_cols = c("CONTAM_HARD_FLAG","RAP_mRNA_HARD_FLAG","CHIMERAS_HARD_FLAG","GENES_DETECTED_FLAG")
  is_failed_qc = rowSums(metrics_df[,flag_cols]=="FAIL") > 0
  metrics_df$QC_FLAG = ifelse(is_failed_qc, "FAILED", "PASS")
  
  metrics_df_subset = metrics_df[,c("entity:sample_id","LC-SET","PDO",flag_cols,"QC_FLAG")]
  metrics_df_subset
}

filter_transcripts = function(count_data) {
  # Return names of the transcripts to keep
  tx_counts = rowSums(count_data[,!names(count_data) %in% c("Name","Description")])
  tx_counts_df = as.data.frame(tx_counts)
  tx_counts_df$name = count_data$Name
  tx_counts_df$log_tx_counts = log(tx_counts_df$tx_counts+1)
  tx_counts_df_filtered = tx_counts_df[tx_counts_df$tx_counts>=1000,]
  
  # Compare the distributions of tx counts before and after filtering
  prefilter_tx_histo = ggplot(tx_counts_df,aes(log_tx_counts)) +
    geom_histogram() + labs(x="log(tx count + 1)",title="Pre-filter Tx counts")
  postfilter_tx_histo = ggplot(tx_counts_df_filtered,aes(log_tx_counts)) +
    geom_histogram() + labs(x="log(tx count + 1)",title="Post-filter Tx counts (Tx <= 10 removed)")
  tx_plots = plot_grid(prefilter_tx_histo,postfilter_tx_histo,labels="AUTO")
  #ggsave("reference_data/RNA_Seq/outputs/plots/00_tx_counts_pre_post_filtering.png",tx_plots,width=12)
  
  return(tx_counts_df_filtered$name)
}

collapse_count_data = function(count_data) {
  # Given count data, collapse transcripts associated with the same gene into one entry
  count_data$Name = NULL
  collapased_counts = count_data %>% 
    group_by(Description) %>%
    summarise(across(everything(), sum))
  collapased_counts = as.data.frame(collapased_counts)
  rownames(collapased_counts) = collapased_counts$Description
  collapased_counts$Description = NULL
  return(collapased_counts)
}


## We first filter out samples that failed QC
metrics_df = read.csv("data/raw/sample_tables/RNA_sample_Feb19_2024.tsv",sep="\t",check.names = FALSE)
sample_qc_summary_df = sample_qc_summary(metrics_df,count_data) 
sample_qc_summary_df["BATCH"] = sapply(strsplit(sample_qc_summary_df[["LC-SET"]], ","), `[`,1)
sample_qc_summary_df[["individual_alias"]] = paste0("ASCProject_",sapply(strsplit(sample_qc_summary_df$`entity:sample_id`, "_"), `[`, 3))
sample_qc_summary_df[["sample_alias"]] = paste0(sample_qc_summary_df[["individual_alias"]],"_",sapply(strsplit(sample_qc_summary_df$`entity:sample_id`, "_"), `[`, 4))

## We then filter out samples that do not have clinical metadata
#clin_data = read.csv("data/processed/sample_clin_data.tsv",sep="\t",check.names = FALSE)
clin_data = read.csv("data/processed/sample_clin_data.tsv",sep="\t",check.names = FALSE)

sample_qc_summary_df[["HAS_METADATA"]] = sample_qc_summary_df$individual_alias %in% clin_data$individual_alias
sample_qc_summary_df[["ALL_FILTERS_PASSED"]] = (sample_qc_summary_df[["QC_FLAG"]]=="PASS") & (sample_qc_summary_df[["HAS_METADATA"]]==TRUE)
sample_qc_summary_df[["SAMPLE_DATA_TYPE"]] = "tumor_RNA"
## For record keeping, write a list outlining the pass / fail status of each sample and the reason
write.table(sample_qc_summary_df,file="data/processed/qc_tables/rna_sample_filter_status.tsv",sep="\t",row.names = FALSE)


# Save filtered samples
pass_qc_samples = sample_qc_summary_df[sample_qc_summary_df$ALL_FILTERS_PASSED==TRUE,"entity:sample_id"]
count_data_filtered = count_data[,c("Name","Description",pass_qc_samples)]

sample_meta_data = read.csv("data/processed/sample_clin_data.tsv",sep="\t",check.names = FALSE)
# sample_meta_data_merged = merge(sample_qc_summary_df,clin_data,by.x="individual_alias",by.y="individual_alias",all.x=TRUE)

common_ids = intersect(sample_meta_data$`entity:sample_id`, colnames(count_data_filtered))
sample_meta_data_filtered = sample_meta_data[sample_meta_data$`entity:sample_id` %in% common_ids,]

# Remove tx with 100 or less total counts across all samples
keep_tx = filter_transcripts(count_data_filtered)
count_data_filtered = count_data_filtered[count_data_filtered$Name %in% keep_tx,]

## Collapse transcript counts into gene counts
collapsed_filtered_counts = collapse_count_data(count_data_filtered)

### 2025-05-04
## Check if we do a z-score normalization against GTEX breast tissues whether that would change the outcome
gtex_breast_rna_counts = read.csv("data/raw/rna/GTEx/gene_reads_v10_breast_mammary_tissue.gct",sep="\t",skip=2,check.names = FALSE)
#gtex_breast_rna_counts = read.csv("data/raw/rna/GTEx/gene_reads_v10_skin_sun_exposed_lower_leg.gct",sep="\t",skip=2,check.names = FALSE)

gtex_breast_rna_counts_collapsed = collapse_count_data(gtex_breast_rna_counts)
gtex_breast_rna_shared_name = intersect(gtex_breast_rna_counts$Description,count_data_filtered$Description)


asc_breast_samples = sample_meta_data %>% filter(`Primary Site (Recombined)`=="Breast (Cutaneous)" | `Primary Site (Recombined)`=="Breast (Parenchymal)") %>% pull(`entity:sample_id`)
asc_breast_samples_with_rna = intersect(asc_breast_samples,colnames(collapsed_filtered_counts))
collapsed_filtered_counts_breast = collapsed_filtered_counts[,asc_breast_samples_with_rna]
gtex_asc_shared_genes = intersect(rownames(gtex_breast_rna_counts_collapsed),rownames(collapsed_filtered_counts_breast))
gtex_asc_breast_merged = cbind(collpased_filtered_counts_breast[gtex_asc_shared_genes,],gtex_breast_rna_counts_collapsed[gtex_asc_shared_genes,])

## Write the merged version for purity calculation
library(estimate)
asc_gtex_collapsed_count_path="/Users/hoyin/Downloads/gene_tpm_v10_bladder_processed.tsv"
#asc_gtex_collapsed_count_path="data/processed/rna/asc_merged_with_gtex.txt"

write.table(gtex_asc_breast_merged,asc_gtex_collapsed_count_path,sep="\t",quote = FALSE)
write.table(t(gtex_asc_breast_merged),asc_gtex_collapsed_count_path,sep="\t",quote = FALSE)

filterCommonGenes(input.f=asc_gtex_collapsed_count_path, output.f="data/processed/rna/asc_merged_with_gtex_collapsed_counts.gct", id="GeneSymbol")
estimateScore(input.ds = "data/processed/rna/asc_merged_with_gtex_collapsed_counts.gct",output.ds = "data/processed/rna/asc_merged_with_gtex_estimate_score.gct")#,platform = "illumina"
asc_tex_merged_estimate_score_table = read.csv("data/processed/rna/asc_merged_with_gtex_estimate_score.gct",sep="\t",skip=2,check.names = FALSE)
asc_tex_merged_estimate_scores = asc_tex_merged_estimate_score_table[asc_tex_merged_estimate_score_table$NAME == "ESTIMATEScore", -c(1,2)]
## Using the equation indicated in the original publication
asc_tex_merged_purity_estimates = cos(0.6049872018+0.0001467884*asc_tex_merged_estimate_scores)
asc_tex_merged_purity_estimates_table = as.data.frame(t(asc_tex_merged_purity_estimates))
colnames(asc_tex_merged_purity_estimates_table) = c("ESTIMATE_purity")
rownames(asc_tex_merged_purity_estimates_table) <- gsub("\\.", "-", rownames(asc_tex_merged_purity_estimates_table))
asc_tex_merged_purity_estimates_table = cbind("entity:sample_id" = rownames(asc_tex_merged_purity_estimates_table),asc_tex_merged_purity_estimates_table)
write.table(asc_tex_merged_purity_estimates_table,  # Add row names as a new column
            file = "data/processed/rna/asc_merged_with_gtex_estimate_score.tsv",  # Change to .csv if needed
            sep = "\t",  # Use "," for CSV files
            quote = FALSE,
            row.names = FALSE  # Prevent duplicate row names
)


#
# # 6. Create metadata
# breast_tumor_ids = colnames(collapsed_filtered_counts_breast)
# breast_gtex_ids = colnames(gtex_breast_rna_counts_collapsed)
# breast_condition = c(rep("tumor", length(breast_tumor_ids)),
#               rep("normal", length(breast_gtex_ids)))
# breast_sample_info = data.frame(
#   row.names = c(breast_tumor_ids, breast_gtex_ids),
#   condition = factor(breast_condition, levels = c("normal", "tumor"))
# )
# 
# # 7. DESeq2 setup
# breast_dds = DESeqDataSetFromMatrix(
#   countData = gtex_asc_breast_merged,
#   colData = breast_sample_info,
#   design = ~ condition
# )
# breast_dds = estimateSizeFactors(breast_dds)
# breast_norm_counts = counts(breast_dds, normalized = TRUE)
# 
# library(sva)
# # 5. Run svaseq to estimate surrogate variables
# mod <- model.matrix(~ condition, data = colData(breast_dds))
# mod0 <- model.matrix(~ 1, data = colData(breast_dds))  # null model
# 
# svobj <- svaseq(breast_norm_counts, mod, mod0)
# 
# # 6. Add surrogate variables to sample_info
# for (i in seq_len(ncol(svobj$sv))) {
#   breast_sample_info[[paste0("SV", i)]] <- svobj$sv[, i]
# }
# 
# # 7. Redefine DESeq2 object with surrogate variables
# breast_design_formula <- as.formula(paste("~", paste0("SV", seq_len(ncol(svobj$sv)), collapse = " + "), "+ condition"))
# breast_dds = DESeqDataSetFromMatrix(
#   countData = gtex_asc_breast_merged,
#   colData = breast_sample_info,
#   design = breast_design_formula
# )
# 
# breast_dds = DESeq(breast_dds)
# breast_dds_res = results(breast_dds, contrast = c("condition", "tumor", "normal"))
# 
# # 10. Order and filter for overexpressed genes
# breast_dds_res_ordered = breast_dds_res[order(breast_dds_res$pvalue), ]
# breast_overexpressed_genes = breast_dds_res_ordered[breast_dds_res_ordered$log2FoldChange > 1 & breast_dds_res_ordered$padj < 0.05, ]
# write.table(breast_dds_res_ordered,file = "02_Tumor_RNA_Analysis/outputs/DEGs/2025_05_04_GTEx_DEGs/02_GTEx_breast_overexpression_with_sva.tsv", sep = "\t", quote = FALSE)
# 
# # 1. Perform variance-stabilizing transformation
# vsd <- vst(breast_dds, blind = TRUE)  # blind = TRUE avoids using condition info for transformation
# 
# # 2. Extract PCA data
# pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# 
# # 3. Quick PCA plot
# ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
#   geom_point(size = 3) +
#   labs(
#     title = "PCA of Tumor vs. GTEx RNA-seq (Breast)",
#     x = paste0("PC1: ", percentVar[1], "% variance"),
#     y = paste0("PC2: ", percentVar[2], "% variance")
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 14, face = "bold"),
#     axis.title = element_text(size = 12),
#     legend.title = element_blank(),
#     legend.position = "top"
#   )
# 
# #Heatmap(gtex_asc_breast_merged)


## Estimate Tumor Purity using ESTIMATE (Yoshihara K., 2013)
library(estimate)
collapsed_count_path="data/processed/rna/collapsed_counts.txt"
write.table(collapsed_filtered_counts,collapsed_count_path,sep="\t",quote = FALSE)
filterCommonGenes(input.f=collapsed_count_path, output.f="data/processed/rna/collapsed_counts.gct", id="GeneSymbol")
estimateScore(input.ds = "data/processed/rna/collapsed_counts.gct",output.ds = "data/processed/rna/estimate_score.gct",platform = "illumina")
estimate_score_table = read.csv("data/processed/rna/estimate_score.gct",sep="\t",skip=2,check.names = FALSE)
estimate_scores = estimate_score_table[estimate_score_table$NAME == "ESTIMATEScore", -c(1,2)]
## Using the equation indicated in the original publication
purity_estimates = cos(0.6049872018+0.0001467884*estimate_scores)
purity_estimates_table = as.data.frame(t(purity_estimates))
colnames(purity_estimates_table) = c("ESTIMATE_purity")
rownames(purity_estimates_table) <- gsub("\\.", "-", rownames(purity_estimates_table))
purity_estimates_table = cbind("entity:sample_id" = rownames(purity_estimates_table),purity_estimates_table)
write.table(purity_estimates_table,  # Add row names as a new column
  file = "data/processed/rna/estimate_score.tsv",  # Change to .csv if needed
  sep = "\t",  # Use "," for CSV files
  quote = FALSE,
  row.names = FALSE  # Prevent duplicate row names
)

## Add estimated tumor purity to metadata
sample_meta_data_filtered = merge(sample_meta_data_filtered,purity_estimates_table,by="entity:sample_id",all.x=TRUE)

## Check per sample expression
expression_by_sample_long = as.data.frame(collapsed_filtered_counts,check.names=FALSE)
expression_by_sample_long$gene = rownames(collapsed_filtered_counts)
expression_by_sample_long = pivot_longer(
  expression_by_sample_long,
  cols = -gene, 
  names_to = "Sample",
  values_to = "raw_count"
)

expression_by_sample_plot_prefilter = ggplot(expression_by_sample_long, aes(x = Sample, y = log2(raw_count+1))) +
  geom_boxplot() +
  labs(x = "Sample", y = "log2(Raw Count+1)") +
  theme_minimal() +
  theme(axis.text.x=element_blank())
expression_by_sample_plot_prefilter
ggsave("02_Tumor_RNA_Analysis/outputs/plots/00_qc_plots/00_expression_by_sample_plot_raw_counts.png",expression_by_sample_plot_prefilter,width=10,height=3)


## Applying similar filtering and transformation to tpm data
tpm_data_filtered = tpm_data[tpm_data$Name %in% keep_tx , c("Name","Description",sample_meta_data_filtered$`entity:sample_id`)]
collpased_tpm_data = collapse_count_data(tpm_data_filtered)

## Apply VST transform and save vst-normalized count data
dds = DESeqDataSetFromMatrix(
  countData = collapsed_filtered_counts, # the counts values for all samples in our dataset
  #colData = meta_data_filtered, 
  colData = sample_meta_data_filtered, # annotation data for the samples in the counts data frame
  design = ~1 # Here we are not specifying a model
)
dds_norm = vst(dds,blind=FALSE)
collpased_filtered_vst = assays(dds_norm)[[1]]

## Visualize the post-vst count distribution by sample
vst_data_wide = as.data.frame(collpased_filtered_vst,check.names=FALSE)
vst_data_wide$gene = rownames(collpased_filtered_vst)
long_df = pivot_longer(
  vst_data_wide,
  cols = -gene, 
  names_to = "Sample",
  values_to = "VSTcount"
)

expression_by_sample_plot = ggplot(long_df, aes(x = Sample, y = VSTcount)) +
  geom_boxplot() +
  labs(x = "Sample", y = "VST-transformed Expression") +
  theme_minimal() +
  theme(axis.text.x=element_blank())
expression_by_sample_plot
ggsave("02_Tumor_RNA_Analysis/outputs/plots/00_qc_plots/00_expression_by_sample_plot_post_vst.png",expression_by_sample_plot,width=10,height=3)

## Add library size to meta data
total_library_sizes = data.frame("library_size"=colSums(count_data[,-c(1,2)]))
total_library_sizes[["entity:sample_id"]] = rownames(total_library_sizes)
sample_meta_data_filtered = merge(sample_meta_data_filtered,total_library_sizes,by.x="entity:sample_id",all.x=TRUE)

## Write the data
if(FALSE){
  write.csv(collapsed_filtered_counts,file="data/processed/rna/00_filtered_gene_counts.csv")
  write.csv(collpased_filtered_vst,file="data/processed/rna/00_filtered_gene_counts.vst.csv")
  write.csv(collpased_tpm_data, file="data/processed/rna/00_filtered_gene_tpm.csv")
  #write.csv(meta_data_filtered_merged,file="reference_data/RNA_Seq/outputs/filtered_metadata.csv")
  write.csv(sample_meta_data_filtered,file="data/processed/rna/00_filtered_sample_metadata.csv")
}



# 
# ## Process external data
# ## Loh (ASC)
# loh_gene_count_data = read.csv("reference_data/RNA_Seq/external/raw/GSE163359_Raw_gene_counts_matrix_angio.txt",sep="\t",check.names = FALSE)
# collapse_count_data_external = function(count_data) {
#   # Given count data, collapse transcripts associated with the same gene into one entry
#   count_data$Name = NULL
#   collapased_counts = count_data %>% 
#     group_by(GeneName) %>%
#     summarise(across(everything(), sum))
#   collapased_counts = as.data.frame(collapased_counts)
#   rownames(collapased_counts) = collapased_counts$GeneName
#   collapased_counts$GeneName = NULL
#   collapased_counts
# }
# loh_gene_count_data_collpased = collapse_count_data_external(loh_gene_count_data)
# loh_gene_count_data_collpased_filtered = loh_gene_count_data_collpased[rowSums(loh_gene_count_data_collpased) > 100,]
# loh_gene_count_data_collpased_filtered
# 
# ## GTEX (Healthy Breast Tissue)
# gtex_breast_count_data = read.csv("reference_data/RNA_Seq/external/raw/bulk-gex_v8_rna-seq_counts-by-tissue_gene_reads_2017-06-05_v8_breast_mammary_tissue.gct",sep="\t",skip=2,check.names = FALSE,row.names = 1)
# keep_tx_gtex_breast = filter_transcripts(gtex_breast_count_data)
# gtex_breast_count_data_filtered = gtex_breast_count_data[gtex_breast_count_data$Name %in% keep_tx_gtex_breast,]
# ## Collapse transcript counts into gene counts
# gtex_breast_count_data_filtered_collpased = collapse_count_data(gtex_breast_count_data_filtered)
# 
# if(FALSE) {
#   write.csv(loh_gene_count_data_collpased_filtered,file="reference_data/RNA_Seq/external/processed/GSE163359_Raw_gene_counts_matrix_angio_filtered.csv")
#   write.csv(gtex_breast_count_data_filtered_collpased,file="reference_data/RNA_Seq/external/processed/GTEX_v8_breast_mammary_tissue_counts_filtered.csv")
# }
