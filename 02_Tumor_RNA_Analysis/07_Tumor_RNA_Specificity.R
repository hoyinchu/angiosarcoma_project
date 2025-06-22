library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(DESeq2)
library(Seurat)
library(fgsea)
library(ggrepel)
#library(scSigR)

setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts")

## Load Seurat object from 10x format and save
## Takes a quite a bit of time so only do this once
if (FALSE) {
  so_raw = Seurat::Read10X(data.dir="./data/public/GTEx/GTEx_8_tissues_snRNAseq_atlas_071421_matrices")
  so_metadata = read.csv("./data/public/GTEx/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad.metadata.csv")
  gtex_so = CreateSeuratObject(counts = so_raw, meta.data = so_metadata)
  SaveSeuratRds(gtex_so,file = "./data/public/GTEx/GTExSeuratObj2025.rds")
}

## Load Seurat Object
gtex_so = readRDS("./data/public/GTEx/GTExSeuratObj2025.rds")
#gtex_so@meta.data$Cell_Type = gtex_so@meta.data$Cell.types.level.3

## Subset into tissue of interest
gtex_breast_so = subset(gtex_so,tissue == "breast")
gtex_skin_so = subset(gtex_so,tissue == "skin")
gtex_heart_so = subset(gtex_so, tissue == "heart")
gtex_extremities_so = subset(gtex_so, tissue == "skeletalmuscle")

## Calculate signatureusing scSigR
## Temporarily switch directory to scSigR and source all the functions
## Import the RunSigR function from scSigR
scSigR_scripts_path = "/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scSigR/R"
setwd(scSigR_scripts_path)
# Source all R scripts in the directory
r_files <- list.files(pattern = "\\.R$")
for (file in r_files) {
  print(file)
  source(file)
}

gtex_breast_signature = RunSigR(gtex_breast_so,cell_types=c("Stromal", "Epithelial", "Immune"),celltype_column="Cell.types.level.3")
gtex_skin_signature = RunSigR(gtex_skin_so,cell_types=c("Stromal", "Epithelial", "Immune","Other"),celltype_column="Cell.types.level.3")
gtex_heart_signature = RunSigR(gtex_heart_so,cell_types=c("Glia", "Stromal", "Immune"),celltype_column="Cell.types.level.3")
gtex_extremities_signature = RunSigR(gtex_extremities_so,cell_types=c("Glia", "Immune", "Stromal"),celltype_column="Cell.types.level.3")

gtex_breast_signature_out = tibble::rownames_to_column(gtex_breast_signature, var = "Gene_name")
gtex_skin_signature_out = tibble::rownames_to_column(gtex_skin_signature, var = "Gene_name")
gtex_heart_signature_out = tibble::rownames_to_column(gtex_heart_signature, var = "Gene_name")
gtex_extremities_signature_out = tibble::rownames_to_column(gtex_extremities_signature, var = "Gene_name")

write.table(gtex_breast_signature_out, file = "/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts/data/public/GTEx/GTEx_8_tissues_snRNAseq_atlas_071421_sig_matrices/breast_signature.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gtex_skin_signature_out, file = "/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts/data/public/GTEx/GTEx_8_tissues_snRNAseq_atlas_071421_sig_matrices/skin_signature.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gtex_heart_signature_out, file = "/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts/data/public/GTEx/GTEx_8_tissues_snRNAseq_atlas_071421_sig_matrices/heart_signature.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gtex_extremities_signature_out, file = "/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts/data/public/GTEx/GTEx_8_tissues_snRNAseq_atlas_071421_sig_matrices/extremities_signature.txt", sep = "\t", quote = FALSE, row.names = FALSE)

## Set the working directory back to the project directory
setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts")

gtex_so


unique(gtex_so@meta.data$tissue)
max(gtex_so@meta.data$PercentMito)

# ## Load seurat object from GTEX
# library(SeuratDisk)
# SeuratDisk::Convert(
#   source="./data/public/GTEx/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs_converted.h5ad",
#   dest="h5seurat",
#   overwrite = TRUE
# )
# gtex_so = SeuratDisk::LoadH5Seurat(
#   "./data/public/GTEx/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs_converted.h5seurat",
#   assays = "RNA"
# )
# 
# library(reticulate)
# library(scater)
# library(SeuratDisk)
# use_python("/Users/hoyin/miniforge3/bin/python")
# ad = import("anndata", convert = FALSE)
# #gtex_ad = ad$read_h5ad("./data/public/GTEx/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs_converted.h5ad")
# gtex_ad = ad$read_h5ad("./data/public/GTEx/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad")
# gtex_so = Convert(gtex_ad, to = "seurat")
# 
# #seuratObject = SeuratDisk::LoadH5Seurat("example_dir/example_ad.h5Seurat")
# 
# ## Conversion didnt go well. Need to make changes in scanpy before conversion
# 
# # This creates a copy of this .h5ad object reformatted into .h5seurat inside the example_dir directory
# 
# # This .d5seurat object can then be read in manually
# seuratObject <- LoadH5Seurat("example_dir/example_ad.h5Seurat")


## Check if the DEG results are tissue-specific

## Load GTEX result
gtex_breast = read.csv("./data/public/GTEx/gene_reads_v10_breast_mammary_tissue.gct",sep="\t",skip=2,check.names = FALSE)
gtex_skin = read.csv("./data/public/GTEx/gene_reads_v10_skin_sun_exposed_lower_leg.gct",sep="\t",skip=2,check.names = FALSE)
gtex_muscle = read.csv("./data/public/GTEx/gene_reads_v10_muscle_skeletal.gct",sep="\t",skip=2,check.names = FALSE)
gtex_heart = read.csv("./data/public/GTEx/gene_reads_v10_heart_left_ventricle.gct",sep="\t",skip=2,check.names = FALSE)

## Load the DEG results
hclust_deg_table = read.csv("./02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_hclust_expr_sites_combined_results_all.tsv",sep="\t")
sites_deg_table = read.csv("./02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_sites_combined_results_all.tsv",sep="\t")

## Load the processed count data
collapsed_filtered_counts = read.csv("data/processed/rna/00_filtered_gene_counts.csv",check.names = FALSE,row.names = 1)

## Subset each gtex tables to the genes that are overlapped
gtex_breast_subset = gtex_breast[gtex_breast$Description %in% rownames(collapsed_filtered_counts),]
gtex_skin_subset = gtex_skin[gtex_skin$Description %in% rownames(collapsed_filtered_counts),]
gtex_muscle_subset = gtex_muscle[gtex_muscle$Description %in% rownames(collapsed_filtered_counts),]
gtex_heart_subset = gtex_heart[gtex_heart$Description %in% rownames(collapsed_filtered_counts),]

## Collapse the same way it was done to preprocessed count data
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
gtex_breast_subset = collapse_count_data(gtex_breast_subset)
gtex_skin_subset = collapse_count_data(gtex_skin_subset)
gtex_muscle_subset = collapse_count_data(gtex_muscle_subset)
gtex_heart_subset = collapse_count_data(gtex_heart_subset)

## Combined the collapsed subset tables
gtex_collapsed_counts_combined = cbind(gtex_breast_subset,gtex_skin_subset,gtex_muscle_subset,gtex_heart_subset)

## Create sample metadata
sample_ids = colnames(gtex_collapsed_counts_combined)
tissue_labels = rep(c("Breast", "Skin", "Muscle", "Heart"),
                    times = c(ncol(gtex_breast_subset),
                              ncol(gtex_skin_subset),
                              ncol(gtex_muscle_subset),
                              ncol(gtex_heart_subset)))
sample_metadata = data.frame(
  sample_id = sample_ids,
  tissue = factor(tissue_labels)
)
rownames(sample_metadata) = sample_metadata$sample_id

## Set one hot features
sample_metadata = sample_metadata %>% mutate(
  site.breast = ifelse(tissue == "Breast","Breast","Rest"),
  site.heart = ifelse(tissue == "Heart","Heart","Rest"),
  site.muscle = ifelse(tissue == "Muscle","Muscle","Rest"),
  site.skin = ifelse(tissue == "Skin","Skin","Rest"),
)

one_versus_rest_deseq2 = function(raw_counts,deseq2_metadata,formula_string,constrast,outpath) {
  deseq2_dds = DESeqDataSetFromMatrix(
    countData = raw_counts,
    colData = deseq2_metadata,
    design = as.formula(formula_string)  # Adjusting for covariates
  )
  deseq2_estimates = DESeq(deseq2_dds)
  deseq2_results = results(deseq2_estimates,contrast = constrast)
  deseq2_results$gene = rownames(deseq2_results)
  write.table(deseq2_results,file = outpath, sep = "\t", quote = FALSE)
  return(deseq2_results)
}



gtex_breast_deseq_res = one_versus_rest_deseq2(gtex_collapsed_counts_combined,sample_metadata,
                                                 "~ site.breast",c("site.breast","Breast","Rest"),
                                                 "02_Tumor_RNA_Analysis/outputs/DEGs/2025_05_04_GTEx_DEGs/02_GTEX_breast_vs_rest.tsv")

gtex_heart_deseq_res = one_versus_rest_deseq2(gtex_collapsed_counts_combined,sample_metadata,
                                               "~ site.heart",c("site.heart","Heart","Rest"),
                                               "02_Tumor_RNA_Analysis/outputs/DEGs/2025_05_04_GTEx_DEGs/02_GTEX_heart_vs_rest.tsv")

gtex_muscle_deseq_res = one_versus_rest_deseq2(gtex_collapsed_counts_combined,sample_metadata,
                                               "~ site.muscle",c("site.muscle","Muscle","Rest"),
                                               "02_Tumor_RNA_Analysis/outputs/DEGs/2025_05_04_GTEx_DEGs/02_GTEX_muscle_vs_rest.tsv")

gtex_skin_deseq_res = one_versus_rest_deseq2(gtex_collapsed_counts_combined,sample_metadata,
                                               "~ site.skin",c("site.skin","Skin","Rest"),
                                               "02_Tumor_RNA_Analysis/outputs/DEGs/2025_05_04_GTEx_DEGs/02_GTEX_skin_vs_rest.tsv")

## Combine all results into one table
#deseq_combined_results = bind_rows(deseq_results_list)



# Calculate per-tissue specificity score
so = readRDS("data/processed/rna/ASCSeuratObj2025.rds")

gtex_df = read.csv("./data/public/rna_tissue_gtex.tsv",sep="\t")

calculate_tissue_specificity = function(df) {
  df %>%
    # Step 1: Tissue specificity score and rank (per gene)
    group_by(Gene, Gene.name) %>%
    mutate(
      sum_nTPM = sum(nTPM),
      tissue_specificity = ifelse(sum_nTPM == 0, 0, nTPM / sum_nTPM),
      tissue_specificity_rank = rank(-tissue_specificity, ties.method = "min")
    ) %>%
    ungroup() %>%
    
    # Step 2: Rank and z-score within each tissue (across genes)
    group_by(Tissue) %>%
    mutate(
      nTPM_rank_within_tissue = rank(-nTPM, ties.method = "min"),
      nTPM_mean = mean(nTPM, na.rm = TRUE),
      nTPM_sd = sd(nTPM, na.rm = TRUE),
      nTPM_zscore = ifelse(nTPM_sd == 0, NA, (nTPM - nTPM_mean) / nTPM_sd)
    ) %>%
    ungroup() %>%
    
    # Final columns
    select(
      Gene, Gene.name, Tissue, nTPM, tissue_specificity, tissue_specificity_rank,
      nTPM_rank_within_tissue, nTPM_zscore
    )
}

gtex_specificity_df = calculate_tissue_specificity(gtex_df)
gtex_specificity_df_dedupped = gtex_specificity_df %>%
  arrange(-tissue_specificity) %>%
  distinct(Gene.name,Tissue,.keep_all = TRUE)

gtex_specificity_table = gtex_specificity_df_dedupped %>%
  select(Gene.name, Tissue, tissue_specificity) %>%
  pivot_wider(
    names_from = Tissue,
    values_from = tissue_specificity
  )


## Plot specificty index distribution

# ggplot(gtex_specificity_df, aes(x = tissue_specificity_rank, fill = Tissue, color = Tissue)) +
#   geom_density(alpha = 0.4) +
#   theme_minimal() +
#   labs(
#     title = "Tissue-Specificity Density per Tissue",
#     x = "Tissue Specificity",
#     y = "Density"
#   ) +
#   theme(legend.position = "right")

#heatmap(gtex_specificity_mat)

## Load the DEG results
hclust_deg_table = read.csv("./02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_hclust_expr_sites_combined_results_all.tsv",sep="\t")
sites_deg_table = read.csv("./02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_sites_combined_results_all.tsv",sep="\t")

hclust_deg_table_merged = hclust_deg_table %>% left_join(gtex_specificity_table,by = c("gene"="Gene.name"))
sites_deg_table_merged = sites_deg_table %>% left_join(gtex_specificity_table,by = c("gene"="Gene.name"))

## Helper function to run fast-ORA
run_fora = function(degs,fgsea_sets,pval_col="p_val") {
  #genes = degs[degs$p_val_adj < 0.1,]$gene
  genes = degs[degs[[pval_col]] < 0.05,]$gene
  #print(genes)
  universe = rownames(so@assays$RNA$counts)
  fora_res = fora(fgsea_sets, genes, universe, minSize = 5, maxSize = 500)
  return(fora_res)
}

quantile(sites_deg_table_merged[sites_deg_table_merged$cluster=="ParenchymalBreast",]$breast,na.rm=TRUE)
quantile(sites_deg_table_merged[sites_deg_table_merged$cluster=="CutaneousBreast",]$breast,na.rm=TRUE)
quantile(sites_deg_table_merged[sites_deg_table_merged$cluster=="HNFS",]$skin,na.rm=TRUE)
quantile(sites_deg_table_merged[sites_deg_table_merged$cluster=="Extremities",]$`skeletal muscle`,na.rm=TRUE)
quantile(sites_deg_table_merged[sites_deg_table_merged$cluster=="Heart",]$`heart muscle`,na.rm=TRUE)

## Show gene set selected for ORA is biased towards tissues we claimed
## Show this bias can be removed by thresholding at different specificty values
## Do a gradient of specificity (showing removing specific genes do not impact pathway significance)

# Step 0: Compute neg log p-values
hclust_deg_table_merged$neg_log_pval <- -log10(hclust_deg_table_merged$pvalue)
sites_deg_table_merged$neg_log_pval <- -log10(sites_deg_table_merged$pvalue)

# Step 1: Filter for upregulated genes
sites_deg_table_merged <- sites_deg_table_merged %>%
  filter(log2FoldChange > 0)

# Step 2: Subset by tissue/cluster
sites_deg_table_merged_parenchymal <- sites_deg_table_merged %>% filter(cluster == "ParenchymalBreast")
sites_deg_table_merged_cutaneous <- sites_deg_table_merged %>% filter(cluster == "CutaneousBreast")
sites_deg_table_merged_hnfs <- sites_deg_table_merged %>% filter(cluster == "HNFS")
sites_deg_table_merged_heart <- sites_deg_table_merged %>% filter(cluster == "Heart")
sites_deg_table_merged_extremities <- sites_deg_table_merged %>% filter(cluster == "Extremities")

# Step 3: ORA set assignment (all log2FC > 0 already by filtering above)
sites_deg_table_merged_parenchymal$ora_set <- sites_deg_table_merged_parenchymal$pvalue < 0.05
sites_deg_table_merged_cutaneous$ora_set <- sites_deg_table_merged_cutaneous$pvalue < 0.05
sites_deg_table_merged_hnfs$ora_set <- sites_deg_table_merged_hnfs$pvalue < 0.05
sites_deg_table_merged_heart$ora_set <- sites_deg_table_merged_heart$pvalue < 0.05
sites_deg_table_merged_extremities$ora_set <- sites_deg_table_merged_extremities$pvalue < 0.05

# Step 4: Compute tissue-specific thresholds (based on filtered data)
specificity_thresholds <- tibble(
  Tissue = c("Parenchymal (breast)", "Cutaneous (breast)", "HNFS (skin)", "Heart", "Extremities"),
  cluster = c("ParenchymalBreast", "CutaneousBreast", "HNFS", "Heart", "Extremities"),
  specificity_col = c("breast", "breast", "skin", "heart muscle", "skeletal muscle")
) %>%
  rowwise() %>%
  mutate(specificity_threshold = quantile(
    sites_deg_table_merged[sites_deg_table_merged$cluster == cluster, ][[specificity_col]],
    probs = 0.5,
    na.rm = TRUE
  )) %>%
  ungroup()

# Step 5: Build plotting dataframe
combined_sites_deg_df <- bind_rows(
  sites_deg_table_merged_parenchymal %>%
    transmute(neg_log_pval, specificity = breast, Tissue = "Parenchymal (breast)"),
  sites_deg_table_merged_cutaneous %>%
    transmute(neg_log_pval, specificity = breast, Tissue = "Cutaneous (breast)"),
  sites_deg_table_merged_cutaneous %>%
    transmute(neg_log_pval, specificity = skin, Tissue = "HNFS (skin)"),
  sites_deg_table_merged_heart %>%
    transmute(neg_log_pval, specificity = `heart muscle`, Tissue = "Heart"),
  sites_deg_table_merged_extremities %>%
    transmute(neg_log_pval, specificity = `skeletal muscle`, Tissue = "Extremities")
) %>%
  left_join(specificity_thresholds %>% select(Tissue, specificity_threshold), by = "Tissue") %>%
  mutate(
    Significant = neg_log_pval > -log10(0.05),
    Tissue_Specific = specificity > specificity_threshold,
    Group = case_when(
      Significant & Tissue_Specific ~ "Significant + Tissue-specific",
      Significant ~ "Significant only",
      Tissue_Specific ~ "Tissue-specific only",
      TRUE ~ "Not significant"
    )
  )

# Step 6: Median line per tissue
median_df <- combined_sites_deg_df %>%
  group_by(Tissue) %>%
  summarize(median_specificity = median(specificity, na.rm = TRUE), .groups = "drop")

# Step 7: Scatter plot
combined_sites_deg_df_plot <- ggplot(combined_sites_deg_df, aes(x = neg_log_pval, y = specificity, color = Group)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_hline(data = median_df, aes(yintercept = median_specificity), linetype = "dotted", color = "gray30") +
  facet_wrap(~ Tissue, ncol = 5) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_color_manual(
    values = c(
      "Significant + Tissue-specific" = "firebrick",
      "Significant only" = "steelblue",
      "Tissue-specific only" = "darkgreen",
      "Not significant" = "gray"
    )
  ) +
  labs(
    x = "Gene -log10(DEG p-values)",
    y = "Gene Tissue Specificity Score",
    color = "Annotation"
  )

ggsave(
  "./02_Tumor_RNA_Analysis/outputs/plots/07_specificity_plots/combined_specificity_annotated_plot.png",
  combined_sites_deg_df_plot,
  dpi = 300,
  width = 10,
  height = 4
)

# Step 8: Barplot of tissue specificity among significant vs not
tissue_specificity_summary_v2 <- combined_sites_deg_df %>%
  mutate(Significance = ifelse(Significant, "DEG Significant (p < 0.05)", "DEG Non-Significant (p >= 0.05)")) %>%
  group_by(Tissue, Significance) %>%
  summarize(
    total = n(),
    tissue_specific = sum(Tissue_Specific, na.rm = TRUE),
    proportion_tissue_specific = tissue_specific / total,
    .groups = "drop"
  )

tissue_specificity_barplot_v2 <- ggplot(tissue_specificity_summary_v2, aes(x = Tissue, y = proportion_tissue_specific, fill = Significance)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("DEG Significant (p < 0.05)" = "firebrick", "DEG Non-Significant (p >= 0.05)" = "gray60")) +
  labs(
    x = "Tissue",
    y = "Proportion of Tissue-Specific Genes\n(Tissue-specificity score quantile > 50%)",
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    #axis.text.x = element_text(hjust = 1),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  )

ggsave(
  "./02_Tumor_RNA_Analysis/outputs/plots/07_specificity_plots/tissue_specificity_proportions_sig_vs_nonsig.png",
  tissue_specificity_barplot_v2,
  dpi = 300,
  width = 10,
  height = 4
)

## Among gene used to perform ORA, what are the proportion of those that are specific to the tissue?
library(msigdbr)
## Load gene set libraries
fgsea_hallmark_set = msigdbr(species = "Homo sapiens", category = "H") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_kegg_set = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c1_set = msigdbr(species = "Homo sapiens", category = "C1") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c5_set = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c6_set = msigdbr(species = "Homo sapiens", category = "C6") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c6_set_up_only = fgsea_c6_set[!grepl("_DN$", names(fgsea_c6_set))]
fgsea_c8_set = msigdbr(species = "Homo sapiens", category = "C8") %>% split(x = .$gene_symbol, f = .$gs_name)

# Extract thresholds
parenchymal_thresh <- specificity_thresholds %>% filter(cluster == "ParenchymalBreast") %>% pull(specificity_threshold)
cutaneous_thresh <- specificity_thresholds %>% filter(cluster == "CutaneousBreast") %>% pull(specificity_threshold)
hnfs_thresh <- specificity_thresholds %>% filter(cluster == "HNFS") %>% pull(specificity_threshold)
extremities_thresh <- specificity_thresholds %>% filter(cluster == "Extremities") %>% pull(specificity_threshold)
heart_thresh <- specificity_thresholds %>% filter(cluster == "Heart") %>% pull(specificity_threshold)

# Filter genes below median specificity for each tissue
parenchymal_ora_set_new <- sites_deg_table_merged %>%
  filter(pvalue < 0.05, breast < parenchymal_thresh, cluster == "ParenchymalBreast")

cutaneous_ora_set_new <- sites_deg_table_merged %>%
  filter(pvalue < 0.05, breast < cutaneous_thresh, cluster == "CutaneousBreast")

hnfs_ora_set_new <- sites_deg_table_merged %>%
  filter(pvalue < 0.05, skin < hnfs_thresh, cluster == "HNFS")

extremities_ora_set_new <- sites_deg_table_merged %>%
  filter(pvalue < 0.05, `skeletal muscle` < extremities_thresh, cluster == "Extremities")

heart_ora_set_new <- sites_deg_table_merged %>%
  filter(pvalue < 0.05, `heart muscle` < heart_thresh, cluster == "Heart")


parenchymal_fora = run_fora(parenchymal_ora_set_new,fgsea_hallmark_set,pval_col="pvalue")
cutaneous_fora = run_fora(cutaneous_ora_set_new,fgsea_hallmark_set,pval_col="pvalue")
hnfs_fora = run_fora(hnfs_ora_set_new,fgsea_hallmark_set,pval_col="pvalue")
extremities_fora = run_fora(extremities_ora_set_new,fgsea_hallmark_set,pval_col="pvalue")
heart_fora = run_fora(heart_ora_set_new,fgsea_hallmark_set,pval_col="pvalue")

parenchymal_c6_fora = run_fora(parenchymal_ora_set_new,fgsea_c6_set_up_only,pval_col="pvalue")
cutaneous_c6_fora = run_fora(cutaneous_ora_set_new,fgsea_c6_set_up_only,pval_col="pvalue")
hnfs_c6_fora = run_fora(hnfs_ora_set_new,fgsea_c6_set_up_only,pval_col="pvalue")
extremities_c6_fora = run_fora(extremities_ora_set_new,fgsea_c6_set_up_only,pval_col="pvalue")
heart_c6_fora = run_fora(heart_ora_set_new,fgsea_c6_set_up_only,pval_col="pvalue")

# DEGs (p < 0.05, log2FC > 0) — without specificity filter
parenchymal_ora_set_all <- sites_deg_table_merged %>%
  filter(pvalue < 0.05, cluster == "ParenchymalBreast")

cutaneous_ora_set_all <- sites_deg_table_merged %>%
  filter(pvalue < 0.05, cluster == "CutaneousBreast")

hnfs_ora_set_all <- sites_deg_table_merged %>%
  filter(pvalue < 0.05, cluster == "HNFS")

extremities_ora_set_all <- sites_deg_table_merged %>%
  filter(pvalue < 0.05, cluster == "Extremities")

heart_ora_set_all <- sites_deg_table_merged %>%
  filter(pvalue < 0.05, cluster == "Heart")

# Run FORA on unfiltered DEG sets (Hallmark)
parenchymal_fora_all <- run_fora(parenchymal_ora_set_all, fgsea_hallmark_set, pval_col = "pvalue")
cutaneous_fora_all <- run_fora(cutaneous_ora_set_all, fgsea_hallmark_set, pval_col = "pvalue")
hnfs_fora_all <- run_fora(hnfs_ora_set_all, fgsea_hallmark_set, pval_col = "pvalue")
extremities_fora_all <- run_fora(extremities_ora_set_all, fgsea_hallmark_set, pval_col = "pvalue")
heart_fora_all <- run_fora(heart_ora_set_all, fgsea_hallmark_set, pval_col = "pvalue")

# Run FORA on unfiltered DEG sets (C6 oncogenic)
parenchymal_c6_fora_all <- run_fora(parenchymal_ora_set_all, fgsea_c6_set_up_only, pval_col = "pvalue")
cutaneous_c6_fora_all <- run_fora(cutaneous_ora_set_all, fgsea_c6_set_up_only, pval_col = "pvalue")
hnfs_c6_fora_all <- run_fora(hnfs_ora_set_all, fgsea_c6_set_up_only, pval_col = "pvalue")
extremities_c6_fora_all <- run_fora(extremities_ora_set_all, fgsea_c6_set_up_only, pval_col = "pvalue")
heart_c6_fora_all <- run_fora(heart_ora_set_all, fgsea_c6_set_up_only, pval_col = "pvalue")



plot_fora_comparison <- function(
  fora_all_list,
  fora_filtered_list,
  gene_set_label,
  output_path,
  csv_path = NULL  # Optional CSV export
) {
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(purrr)
  
  # Function to clean pathway labels
  clean_pathway_label <- function(x) {
    x <- gsub("^HALLMARK_", "", x)
    x <- gsub("^C6_", "", x)
    x <- gsub("(_V[0-9]+)?_?(UP|DN)$", "", x)  # remove suffixes like _UP, _DN, _V1_UP
    x
  }
  
  # Compare FORA results for a single tissue
  compare_fora_sets <- function(fora_all, fora_filtered, tissue_name) {
    df <- full_join(
      fora_all %>% select(pathway, pval_all = pval),
      fora_filtered %>% select(pathway, pval_filtered = pval),
      by = "pathway"
    ) %>%
      mutate(
        log10_pval_all = -log10(pval_all),
        log10_pval_filtered = -log10(pval_filtered),
        Tissue = tissue_name,
        sig_all = !is.na(pval_all) & pval_all < 0.05,
        sig_filtered = !is.na(pval_filtered) & pval_filtered < 0.05,
        Enrichment = case_when(
          sig_filtered & sig_all ~ "Significant in Both",
          sig_filtered ~ "Significant using Specificity-Filtered DEGs",
          sig_all ~ "Significant using All DEGs",
          TRUE ~ "Not Significant"
        )
      )
    
    # Shared
    shared <- df %>%
      filter(sig_filtered & sig_all) %>%
      mutate(Annotate = clean_pathway_label(pathway))
    
    # Top 5 unique to all DEGs
    all_only <- df %>%
      filter(sig_all & !sig_filtered & !pathway %in% shared$pathway) %>%
      arrange(pval_all) %>%
      slice_head(n = 5) %>%
      mutate(Annotate = clean_pathway_label(pathway))
    
    # Top 5 unique to filtered
    filtered_only <- df %>%
      filter(sig_filtered & !sig_all & !pathway %in% shared$pathway) %>%
      arrange(pval_filtered) %>%
      slice_head(n = 5) %>%
      mutate(Annotate = clean_pathway_label(pathway))
    
    # Merge
    annotated <- bind_rows(shared, all_only, filtered_only) %>%
      select(pathway, Annotate)
    
    df %>% left_join(annotated, by = "pathway")
  }
  
  # Process all tissues
  comparison_df <- pmap_dfr(
    list(fora_all_list, fora_filtered_list, names(fora_all_list)),
    ~ compare_fora_sets(..1, ..2, tissue_name = ..3)
  )
  
  # Write annotated CSV if requested
  if (!is.null(csv_path)) {
    write.csv(comparison_df, csv_path, row.names = FALSE)
  }
  
  # Colors
  enrichment_colors <- c(
    "Significant in Both" = "purple",
    "Significant using Specificity-Filtered DEGs" = "firebrick",
    "Significant using All DEGs" = "steelblue",
    "Not Significant" = "gray80"
  )
  
  # Plot
  p <- ggplot(comparison_df, aes(x = log10_pval_all, y = log10_pval_filtered, color = Enrichment)) +
    geom_point(alpha = 0.8) +
    geom_text_repel(aes(label = Annotate), size = 2.5, max.overlaps = 100, na.rm = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
    geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "black") +
    facet_wrap(~ Tissue, ncol = 5) +
    scale_color_manual(values = enrichment_colors) +
    theme_minimal() +
    theme(
      legend.position = "top",
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      x = "-log10(p-value) from All DEGs",
      y = "-log10(p-value) from Specificity-Filtered DEGs",
      color = "Enrichment Class"
    )
  
  ggsave(output_path, p, width = 16, height = 6, dpi = 300)
}

# Named lists of results by tissue
hallmark_all <- list(
  "Parenchymal (breast)" = parenchymal_fora_all,
  "Cutaneous (breast)"   = cutaneous_fora_all,
  "HNFS (skin)"          = hnfs_fora_all,
  "Extremities"          = extremities_fora_all,
  "Heart"                = heart_fora_all
)

hallmark_filtered <- list(
  "Parenchymal (breast)" = parenchymal_fora,
  "Cutaneous (breast)"   = cutaneous_fora,
  "HNFS (skin)"          = hnfs_fora,
  "Extremities"          = extremities_fora,
  "Heart"                = heart_fora
)

c6_all <- list(
  "Parenchymal (breast)" = parenchymal_c6_fora_all,
  "Cutaneous (breast)"   = cutaneous_c6_fora_all,
  "HNFS (skin)"          = hnfs_c6_fora_all,
  "Extremities"          = extremities_c6_fora_all,
  "Heart"                = heart_c6_fora_all
)

c6_filtered <- list(
  "Parenchymal (breast)" = parenchymal_c6_fora,
  "Cutaneous (breast)"   = cutaneous_c6_fora,
  "HNFS (skin)"          = hnfs_c6_fora,
  "Extremities"          = extremities_c6_fora,
  "Heart"                = heart_c6_fora
)

# Run the function
library(purrr)
plot_fora_comparison(
  fora_all_list = hallmark_all,
  fora_filtered_list = hallmark_filtered,
  gene_set_label = "Hallmark",
  output_path = "./02_Tumor_RNA_Analysis/outputs/plots/07_specificity_plots/hallmark_fora_logpval_comparison_all_tissues.png"
)

plot_fora_comparison(
  fora_all_list = c6_all,
  fora_filtered_list = c6_filtered,
  gene_set_label = "C6",
  output_path = "./02_Tumor_RNA_Analysis/outputs/plots/07_specificity_plots/c6_fora_logpval_comparison_all_tissues.png"
)

## Specificity should be evaluated by percentile / cutoff


## Seems like overall ORA selected genes are more tissue-specific
## Route 1: Recognize this and say it as a limitation
## Route 2: Try deconvolution, repeat analysis

ggplot(sites_deg_table_merged_extremities,aes(x=breast,color=ora_set)) + geom_boxplot()

# Wilcoxon test: compare 'breast' specificity between significant vs non-significant
wilcox_test = wilcox.test(
  skin ~ ora_set,
  data = sites_deg_table_merged_hnfs,
  alternative = "less"
)

# Print test statistic
print(wilcox_test)

t_test_result <- t.test(
  breast ~ ora_set,
  data = sites_deg_table_merged_parenchymal,
  alternative = "less"  # tests if mean in nom_sig == TRUE is < mean in nom_sig == FALSE
)

print(t_test_result)

# Density plot
ggplot(sites_deg_table_merged_heart, aes(x = skin, color = ora_set, fill = ora_set)) +
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(
    title = "Tissue specificity of ORA-selected genes",
    subtitle = paste0("Wilcoxon p = ", signif(wilcox_test$p.value, 3)),
    x = "Tissue specificity (breast)",
    fill = "Nominally significant",
    color = "Nominally significant"
  )
