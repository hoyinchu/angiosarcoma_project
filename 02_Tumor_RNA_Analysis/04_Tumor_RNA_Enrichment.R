library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(EnhancedVolcano)
library(msigdbr)
library(fgsea)
library(forcats)
library(DESeq2)
library(enrichR)
library(cowplot)

setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts")

## Load the previously created seurat data
so = readRDS("data/processed/rna/ASCSeuratObj2025.rds")

## Load the palettes
source("util_scripts/project_palettes.R")

## Make one-versus-rest columns
so@meta.data = so@meta.data %>% mutate(
  expr.clust1 = ifelse(seurat_clusters_by_expression_str == "Cluster 1","Cluster1","Rest"),
  expr.clust2 = ifelse(seurat_clusters_by_expression_str == "Cluster 2","Cluster2","Rest"),
  expr.clust3 = ifelse(seurat_clusters_by_expression_str == "Cluster 3","Cluster3","Rest"),
  expr.clust4 = ifelse(seurat_clusters_by_expression_str == "Cluster 4","Cluster4","Rest"),
  expr.clust5 = ifelse(seurat_clusters_by_expression_str == "Cluster 5","Cluster5","Rest"),
  expr.hclust1 = ifelse(hcluster_by_expr_str == "Cluster 1","Cluster1","Rest"),
  expr.hclust2 = ifelse(hcluster_by_expr_str == "Cluster 2","Cluster2","Rest"),
  expr.hclust3 = ifelse(hcluster_by_expr_str == "Cluster 3","Cluster3","Rest"),
  expr.hclust4 = ifelse(hcluster_by_expr_str == "Cluster 4","Cluster4","Rest"),
  expr.hclust5 = ifelse(hcluster_by_expr_str == "Cluster 5","Cluster5","Rest"),
  site.breast.noncut = ifelse(seurat_clusters_by_site_str == "Breast (Parenchymal)","ParenchymalBreast","Rest"),
  site.breast.cut = ifelse(seurat_clusters_by_site_str == "Breast (Cutaneous)","CutaneousBreast","Rest"),
  site.hnfs = ifelse(seurat_clusters_by_site_str == "HNFS","HNFS","Rest"),
  site.extremities = ifelse(seurat_clusters_by_site_str == "Extremities","Extremities","Rest"),
  site.heart = ifelse(seurat_clusters_by_site_str == "Heart","Heart","Rest"),
  site.others = ifelse(seurat_clusters_by_site_str == "Others","Others","Rest"),
)

## DESEQ2-based marker identification
# Convert Seurat object to count matrix and metadata
raw_counts = as.matrix(GetAssayData(so, slot = "counts"))
deseq2_metadata = so@meta.data
deseq2_metadata = deseq2_metadata[colnames(raw_counts),]


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

## Make the tables by site
noncut_breast_deseq_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                    "~ ESTIMATE_purity + rna_inferred_female + site.breast.noncut",c("site.breast.noncut","ParenchymalBreast","Rest"),
                                    "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_noncut_breast_vs_rest.tsv")
cut_breast_deseq_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                                 "~ ESTIMATE_purity + rna_inferred_female + site.breast.cut",c("site.breast.cut","CutaneousBreast","Rest"),
                                                 "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_cut_breast_vs_rest.tsv")
hnfs_deseq_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                                 "~ ESTIMATE_purity + rna_inferred_female + site.hnfs",c("site.hnfs","HNFS","Rest"),
                                                 "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_hnfs_vs_rest.tsv")
heart_deseq_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                                 "~ ESTIMATE_purity + rna_inferred_female + site.heart",c("site.heart","Heart","Rest"),
                                                 "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_heart_vs_rest.tsv")
extremities_deseq_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                                 "~ ESTIMATE_purity + rna_inferred_female + site.extremities",c("site.extremities","Extremities","Rest"),
                                                 "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_extremities_vs_rest.tsv")
others_deseq_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                                 "~ ESTIMATE_purity + rna_inferred_female + site.others",c("site.others","Others","Rest"),
                                                 "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_others_vs_rest.tsv")

## All results regardless of whether they are positive
noncut_breast_deseq_res_all = noncut_breast_deseq_res%>% as.data.frame() %>% mutate(cluster="ParenchymalBreast")
cut_breast_deseq_res_all = cut_breast_deseq_res %>% as.data.frame() %>% mutate(cluster="CutaneousBreast")
hnfs_deseq_res_all = hnfs_deseq_res %>% as.data.frame() %>% mutate(cluster="HNFS")
heart_deseq_res_all = heart_deseq_res %>% as.data.frame() %>% mutate(cluster="Heart")
extremities_deseq_res_all = extremities_deseq_res %>% as.data.frame() %>% mutate(cluster="Extremities")
others_deseq_res_all = others_deseq_res %>% as.data.frame() %>% mutate(cluster="Others")
all_deseq_res_all = rbind(noncut_breast_deseq_res_all,cut_breast_deseq_res_all,hnfs_deseq_res_all,
                          heart_deseq_res_all,extremities_deseq_res_all,others_deseq_res_all)
write.table(all_deseq_res_all,file = "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_sites_combined_results_all.tsv", sep = "\t", quote = FALSE)


## Subset DESEQ2 markers to positives only
noncut_breast_deseq_res_pos = noncut_breast_deseq_res[noncut_breast_deseq_res$log2FoldChange >= 0,] %>% as.data.frame() %>% mutate(cluster="ParenchymalBreast")
cut_breast_deseq_res_pos = cut_breast_deseq_res[cut_breast_deseq_res$log2FoldChange >= 0,]  %>% as.data.frame() %>% mutate(cluster="CutaneousBreast")
hnfs_deseq_res_pos = hnfs_deseq_res[hnfs_deseq_res$log2FoldChange >= 0,] %>% as.data.frame() %>% mutate(cluster="HNFS")
heart_deseq_res_pos = heart_deseq_res[heart_deseq_res$log2FoldChange >= 0,] %>% as.data.frame() %>% mutate(cluster="Heart")
extremities_deseq_res_pos = extremities_deseq_res[extremities_deseq_res$log2FoldChange >= 0,] %>% as.data.frame() %>% mutate(cluster="Extremities")
others_deseq_res_pos = others_deseq_res[others_deseq_res$log2FoldChange >= 0,] %>% as.data.frame() %>% mutate(cluster="Others")
## Combine into one
all_deseq_res_pos = rbind(noncut_breast_deseq_res_pos,cut_breast_deseq_res_pos,hnfs_deseq_res_pos,
                          heart_deseq_res_pos,extremities_deseq_res_pos,others_deseq_res_pos)

## Write the combined results
write.table(all_deseq_res_pos,file = "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_sites_combined_results_pos.tsv", sep = "\t", quote = FALSE)

## Repeat the analysis above but with HClust derived clusters instead
## Make the tables by site
hclust1_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                                 "~ ESTIMATE_purity + rna_inferred_female + expr.hclust1",c("expr.hclust1","Cluster1","Rest"),
                                                 "02_Tumor_RNA_Analysis/outputs/DEGs/02_hclust_Cluster1_vs_rest.tsv")
hclust2_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                              "~ ESTIMATE_purity + rna_inferred_female + expr.hclust2",c("expr.hclust2","Cluster2","Rest"),
                                              "02_Tumor_RNA_Analysis/outputs/DEGs/02_hclust_Cluster2_vs_rest.tsv")
hclust3_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                        "~ ESTIMATE_purity + rna_inferred_female + expr.hclust3",c("expr.hclust3","Cluster3","Rest"),
                                        "02_Tumor_RNA_Analysis/outputs/DEGs/02_hclust_Cluster3_vs_rest.tsv")
hclust4_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                         "~ ESTIMATE_purity + rna_inferred_female + expr.hclust4",c("expr.hclust4","Cluster4","Rest"),
                                         "02_Tumor_RNA_Analysis/outputs/DEGs/02_hclust_Cluster4_vs_rest.tsv")
hclust5_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                               "~ ESTIMATE_purity + rna_inferred_female + expr.hclust5",c("expr.hclust5","Cluster5","Rest"),
                                               "02_Tumor_RNA_Analysis/outputs/DEGs/02_hclust_Cluster5_vs_rest.tsv")

hclust1_res_all = hclust1_res %>% as.data.frame() %>% mutate(cluster="Cluster 1")
hclust2_res_all = hclust2_res %>% as.data.frame() %>% mutate(cluster="Cluster 2")
hclust3_res_all = hclust3_res %>% as.data.frame() %>% mutate(cluster="Cluster 3")
hclust4_res_all = hclust4_res %>% as.data.frame() %>% mutate(cluster="Cluster 4")
hclust5_res_all = hclust5_res %>% as.data.frame() %>% mutate(cluster="Cluster 5")
all_hclust_deseq_res_all = rbind(hclust1_res_all,hclust2_res_all,hclust3_res_all,hclust4_res_all,hclust5_res_all)
write.table(all_hclust_deseq_res_all,file = "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_hclust_expr_sites_combined_results_all.tsv", sep = "\t", quote = FALSE)

## Subset DESEQ2 markers to positives only
hclust1_res_pos = hclust1_res[hclust1_res$log2FoldChange >= 0,] %>% as.data.frame() %>% mutate(cluster="Cluster 1")
hclust2_res_pos = hclust2_res[hclust2_res$log2FoldChange >= 0,]  %>% as.data.frame() %>% mutate(cluster="Cluster 2")
hclust3_res_pos = hclust3_res[hclust3_res$log2FoldChange >= 0,] %>% as.data.frame() %>% mutate(cluster="Cluster 3")
hclust4_res_pos = hclust4_res[hclust4_res$log2FoldChange >= 0,] %>% as.data.frame() %>% mutate(cluster="Cluster 4")
hclust5_res_pos = hclust5_res[hclust5_res$log2FoldChange >= 0,] %>% as.data.frame() %>% mutate(cluster="Cluster 5")
## Combine into one
all_hclust_deseq_res_pos = rbind(hclust1_res_pos,hclust2_res_pos,hclust3_res_pos,hclust4_res_pos,hclust5_res_pos)
## Write the combined results
write.table(all_hclust_deseq_res_pos,file = "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_hclust_expr_sites_combined_results_pos.tsv", sep = "\t", quote = FALSE)


## Alternatively try using Seurat's built in wilcoxon test
## Find markers associated with each cluster
Idents(so) = so@meta.data$seurat_clusters_by_expression_str

## Marker by unsupervised clusters
all_markers_by_expr_cluster = FindAllMarkers(
  so,
  #slot="counts",
  slots="scale.data",
  test.use = "wilcox",
  #fc.slot = "counts",
  only.pos = TRUE,
  #min.pct = 0,
  #min.cells.group=0,
  #min.cells.feature = 0,
  logfc.threshold = 0,
  return.thresh = 1
)
## Marker by clusters marked by primary sites
Idents(so) = so@meta.data$seurat_clusters_by_site_str
all_markers_by_site_cluster = FindAllMarkers(
  so,
  #slot="counts",
  slots="scale.data",
  test.use = "wilcox",
  #fc.slot = "counts",
  only.pos = TRUE,
  #min.pct = 0,
  #min.cells.group=0,
  #min.cells.feature = 0,
  logfc.threshold = 0,
  return.thresh = 1
)

## Make volcano plot based on the marker dataframe provided
make_volcano = function(
  marker_df,
  highlight_genes,
  fc_cutoff=0.5,
  log2FC_col = "avg_log2FC",
  pval_col="p_val_adj",
  pval_cutoff=5e-2,
  title="",
  subtitle="",
  draw_bonferroni=FALSE
) {
  bonferroni_thresh = 0.05/dim(marker_df)[1]
  marker_df_sig = marker_df %>% filter(!!sym(log2FC_col) >= fc_cutoff,!!sym(pval_col) < pval_cutoff)
  top_fc_genes = marker_df_sig %>% arrange(-!!sym(log2FC_col)) %>% head(10) %>% pull(gene)
  top_pval_genes = marker_df_sig %>% arrange(!!sym(pval_col)) %>% head(10) %>% pull(gene)
  top_bonferrino_sig_fc_genes = marker_df_sig %>% filter(!!sym(pval_col) < bonferroni_thresh) %>% arrange(-!!sym(log2FC_col)) %>% head(10) %>% pull(gene)
  
  genes_to_highlight = c(top_fc_genes,top_pval_genes,top_bonferrino_sig_fc_genes,highlight_genes)
  vol_plot = EnhancedVolcano(
    marker_df,
    lab=marker_df$gene,
    x=log2FC_col,
    y=pval_col,
    selectLab = genes_to_highlight,
    pCutoff = pval_cutoff,
    FCcutoff = fc_cutoff,
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    title = title,
    subtitle = subtitle
  )
  if (draw_bonferroni) {
    vol_plot = vol_plot + ggplot2::geom_hline(yintercept=-log10(bonferroni_thresh),linetype="dashed",color="purple")
  }
  return(vol_plot)
}

## Make volcanoes for site clusters (DESEQ2)
deseq_s1_volcano = make_volcano(data.frame(noncut_breast_deseq_res_pos),c(),title="Overexpressed Genes",subtitle = "Breast (Parenchymal) vs. Rest",pval_col = "pvalue", log2FC_col = "log2FoldChange", draw_bonferroni=TRUE)
deseq_s2_volcano = make_volcano(data.frame(cut_breast_deseq_res_pos),c(),title="Overexpressed Genes",subtitle = "Breast (Cutaneous) vs. Rest",pval_col = "pvalue", log2FC_col = "log2FoldChange", draw_bonferroni=TRUE)
deseq_s3_volcano = make_volcano(data.frame(hnfs_deseq_res),c(),title="Overexpressed Genes",subtitle = "HNFS vs. Rest",pval_col = "pvalue", log2FC_col = "log2FoldChange", draw_bonferroni=TRUE)
deseq_s4_volcano = make_volcano(data.frame(heart_deseq_res),c(),title="Overexpressed Genes",subtitle = "Heart vs. Rest",pval_col = "pvalue", log2FC_col = "log2FoldChange", draw_bonferroni=TRUE)
deseq_s5_volcano = make_volcano(data.frame(extremities_deseq_res),c(),title="Overexpressed Genes",subtitle = "Extremities vs. Rest",pval_col = "pvalue", log2FC_col = "log2FoldChange", draw_bonferroni=TRUE)
deseq_s6_volcano = make_volcano(data.frame(others_deseq_res_pos),c(),title="Overexpressed Genes",subtitle = "Others vs. Rest",pval_col = "pvalue", log2FC_col = "log2FoldChange", draw_bonferroni=TRUE)

deseq_volcanoes_combined = plot_grid(
  deseq_s1_volcano,deseq_s2_volcano,deseq_s3_volcano,
  deseq_s4_volcano,deseq_s5_volcano,deseq_s6_volcano,ncol = 6
)
ggsave("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_cluster_DESEQ2_DEGs_nominal_site_cluster_volcanoes_combined.png",deseq_volcanoes_combined,dpi=300,width=32,height = 12)
ggsave("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_cluster_DESEQ2_DEGs_nominal_site_cluster_volcanoes_combined.pdf",deseq_volcanoes_combined,dpi=300,width=32,height = 12)


## Make marker and volcanoes for the unsupervised clusters
c1_markers = all_markers_by_expr_cluster %>% filter(cluster=="Cluster 1")
c2_markers = all_markers_by_expr_cluster %>% filter(cluster=="Cluster 2")
c3_markers = all_markers_by_expr_cluster %>% filter(cluster=="Cluster 3")
c4_markers = all_markers_by_expr_cluster %>% filter(cluster=="Cluster 4")
c5_markers = all_markers_by_expr_cluster %>% filter(cluster=="Cluster 5")

c1_volcano = make_volcano(c1_markers,c(),title="Overexpressed Genes",subtitle = "Cluster 1",pval_col = "p_val", draw_bonferroni=TRUE)
c2_volcano = make_volcano(c2_markers,c(),title="Overexpressed Genes",subtitle = "Cluster 2",pval_col = "p_val", draw_bonferroni=TRUE)
c3_volcano = make_volcano(c3_markers,c(),title="Overexpressed Genes",subtitle = "Cluster 3",pval_col = "p_val", draw_bonferroni=TRUE)
c4_volcano = make_volcano(c4_markers,c(),title="Overexpressed Genes",subtitle = "Cluster 4",pval_col = "p_val", draw_bonferroni=TRUE)
c5_volcano = make_volcano(c5_markers,c(),title="Overexpressed Genes",subtitle = "Cluster 5",pval_col = "p_val", draw_bonferroni=TRUE)

# c1_volcano = make_volcano(c1_markers,rownames(angiogenesis_df),title="Overexpressed Genes",subtitle = "Cluster 1")
# c2_volcano = make_volcano(c2_markers,rownames(lymphoangiogensis_df),title="Overexpressed Genes",subtitle = "Cluster 2")
# c3_volcano = make_volcano(c3_markers,rownames(flt3_related_df),title="Overexpressed Genes",subtitle = "Cluster 3")
# c4_volcano = make_volcano(c4_markers,rownames(skin_epithelial_df),title="Overexpressed Genes",subtitle = "Cluster 4")
# c5_volcano = make_volcano(c5_markers,rownames(other_markers_df),title="Overexpressed Genes",subtitle = "Cluster 5")

expr_cluster_volcanoes_combined = plot_grid(
  c1_volcano,c2_volcano,c3_volcano,c4_volcano,c5_volcano,ncol = 5
)
ggsave("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_cluster_DEGs_nominal_expr_cluster_volcanoes_combined.png",expr_cluster_volcanoes_combined,dpi=300,width=32,height = 12)
ggsave("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_cluster_DEGs_nominal_expr_cluster_volcanoes_combined.pdf",expr_cluster_volcanoes_combined,dpi=300,width=32,height = 12)


## Make marker and volcanoes for the primary site clusters
## Make marker and volcanoes for the unsupervised clusters
s1_markers = all_markers_by_site_cluster %>% filter(cluster=="Breast (Parenchymal)")
s2_markers = all_markers_by_site_cluster %>% filter(cluster=="Breast (Cutaneous)")
s3_markers = all_markers_by_site_cluster %>% filter(cluster=="HNFS")
s4_markers = all_markers_by_site_cluster %>% filter(cluster=="Extremities")
s5_markers = all_markers_by_site_cluster %>% filter(cluster=="Heart")
s6_markers = all_markers_by_site_cluster %>% filter(cluster=="Others")

s1_volcano = make_volcano(s1_markers,c(),title="Overexpressed Genes",subtitle = "Breast (Parenchymal)",pval_col = "p_val", draw_bonferroni=TRUE)
s2_volcano = make_volcano(s2_markers,c(),title="Overexpressed Genes",subtitle = "Breast (Cutaneous)",pval_col = "p_val", draw_bonferroni=TRUE)
s3_volcano = make_volcano(s3_markers,c(),title="Overexpressed Genes",subtitle = "HNFS",pval_col = "p_val", draw_bonferroni=TRUE)
s4_volcano = make_volcano(s4_markers,c(),title="Overexpressed Genes",subtitle = "Extremities",pval_col = "p_val", draw_bonferroni=TRUE)
s5_volcano = make_volcano(s5_markers,c(),title="Overexpressed Genes",subtitle = "Heart",pval_col = "p_val", draw_bonferroni=TRUE)
s6_volcano = make_volcano(s6_markers,c(),title="Overexpressed Genes",subtitle = "Others",pval_col = "p_val", draw_bonferroni=TRUE)

# c1_volcano = make_volcano(c1_markers,rownames(angiogenesis_df),title="Overexpressed Genes",subtitle = "Cluster 1")
# c2_volcano = make_volcano(c2_markers,rownames(lymphoangiogensis_df),title="Overexpressed Genes",subtitle = "Cluster 2")
# c3_volcano = make_volcano(c3_markers,rownames(flt3_related_df),title="Overexpressed Genes",subtitle = "Cluster 3")
# c4_volcano = make_volcano(c4_markers,rownames(skin_epithelial_df),title="Overexpressed Genes",subtitle = "Cluster 4")
# c5_volcano = make_volcano(c5_markers,rownames(other_markers_df),title="Overexpressed Genes",subtitle = "Cluster 5")

site_cluster_volcanoes_combined = plot_grid(
  s1_volcano,s2_volcano,s3_volcano,s4_volcano,s5_volcano,s6_volcano,ncol = 6
)
ggsave("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_cluster_nominal_DEGs_site_cluster_volcanoes_combined.png",site_cluster_volcanoes_combined,dpi=300,width=32,height=12)
ggsave("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_cluster_nominal_DEGs_site_cluster_volcanoes_combined.pdf",site_cluster_volcanoes_combined,dpi=300,width=32,height=12)


## Make marker and volcanoes for the unsupervised hierarchical clusters
hc1_markers = all_hclust_deseq_res_pos %>% filter(cluster=="Cluster 1")
hc2_markers = all_hclust_deseq_res_pos %>% filter(cluster=="Cluster 2")
hc3_markers = all_hclust_deseq_res_pos %>% filter(cluster=="Cluster 3")
hc4_markers = all_hclust_deseq_res_pos %>% filter(cluster=="Cluster 4")
hc5_markers = all_hclust_deseq_res_pos %>% filter(cluster=="Cluster 5")

hc1_volcano = make_volcano(hc1_markers,c(),title="Overexpressed Genes",subtitle = "Cluster 1",log2FC_col = "log2FoldChange",pval_col = "pvalue", draw_bonferroni=TRUE)
hc2_volcano = make_volcano(hc2_markers,c(),title="Overexpressed Genes",subtitle = "Cluster 2",log2FC_col = "log2FoldChange",pval_col = "pvalue", draw_bonferroni=TRUE)
hc3_volcano = make_volcano(hc3_markers,c(),title="Overexpressed Genes",subtitle = "Cluster 3",log2FC_col = "log2FoldChange",pval_col = "pvalue", draw_bonferroni=TRUE)
hc4_volcano = make_volcano(hc4_markers,c(),title="Overexpressed Genes",subtitle = "Cluster 4",log2FC_col = "log2FoldChange",pval_col = "pvalue", draw_bonferroni=TRUE)
hc5_volcano = make_volcano(hc5_markers,c(),title="Overexpressed Genes",subtitle = "Cluster 5",log2FC_col = "log2FoldChange",pval_col = "pvalue", draw_bonferroni=TRUE)


hclust_expr_cluster_volcanoes_combined = plot_grid(
  hc1_volcano,hc2_volcano,hc3_volcano,hc4_volcano,hc5_volcano,ncol = 5
)
ggsave("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_hclust_DEGs_nominal_expr_cluster_volcanoes_combined.png",hclust_expr_cluster_volcanoes_combined,dpi=300,width=32,height = 12)
ggsave("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_hclust_DEGs_nominal_expr_cluster_volcanoes_combined.pdf",hclust_expr_cluster_volcanoes_combined,dpi=300,width=32,height = 12)



## Check if the endogenous virus expression is correlated with IFN-stimulated genes?

 ## Load gene set libraries
fgsea_hallmark_set = msigdbr(species = "Homo sapiens", category = "H") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_kegg_set = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c1_set = msigdbr(species = "Homo sapiens", category = "C1") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c5_set = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c6_set = msigdbr(species = "Homo sapiens", category = "C6") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c6_set_up_only = fgsea_c6_set[!grepl("_DN$", names(fgsea_c6_set))]
fgsea_c8_set = msigdbr(species = "Homo sapiens", category = "C8") %>% split(x = .$gene_symbol, f = .$gs_name)

# ## Helper function to run fast-gsea
# run_fgsea = function(degs,fgsea_sets) {
#   deg_genes = degs %>%
#     arrange(desc(avg_log2FC)) %>% 
#     dplyr::select(gene, avg_log2FC)
#   vec = deg_genes$avg_log2FC; names(vec) = deg_genes$gene
#   fgseaRes = fgseaMultilevel(fgsea_sets, stats = vec)
#   return(fgseaRes)
# }
# 
# ## Helper function to plot fgsea results
# plot_fgsea = function(fgsea_res,title) {
#   fgsea_res_sorted = fgsea_res %>% arrange((pval))
#   p1 = ggplot(fgsea_res_sorted, aes(x = ES, y = -log10(padj), label=pathway)) + 
#     geom_text_repel(data = fgsea_res_sorted[fgsea_res_sorted$pval < 0.05,]) +
#     geom_hline(yintercept = -log10(0.05),linetype = "dashed") +
#     geom_vline(xintercept = 0,linetype = "dashed") +
#     geom_point() + labs(x = "Enrichment score", y = "-log10 padj", title=title)# +
#   #pretty_plot(fontsize = 8) + L_border()
#   return(p1)
# }

# Enrichr Analysis
dbs = listEnrichrDbs()
selected_dbs = c("MSigDB_Hallmark_2020","KEGG_2021_Human","NCI-Nature_2016","MSigDB_Oncogenic_Signatures")
background_set = rownames(so@assays$RNA$counts)

run_enrichr_analysis = function(input_genes,selected_dbs,prefix) {
  enriched_results = enrichr(input_genes, selected_dbs,include_overlap = TRUE)#, background = background_set)
  printEnrich(enriched_results,outFile = "excel",prefix=prefix)
}

run_enrichr_analysis(all_deseq_res_pos[(all_deseq_res_pos$cluster=="CutaneousBreast")&(all_deseq_res_pos$padj<0.1),]$gene,selected_dbs,
                     "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_cut_breast_vs_rest_Enrichr")
run_enrichr_analysis(all_deseq_res_pos[(all_deseq_res_pos$cluster=="ParenchymalBreast")&(all_deseq_res_pos$padj<0.1),]$gene,selected_dbs,
                     "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_noncut_breast_vs_rest_Enrichr")
run_enrichr_analysis(all_deseq_res_pos[(all_deseq_res_pos$cluster=="HNFS")&(all_deseq_res_pos$padj<0.1),]$gene,selected_dbs,
                     "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_hnfs_vs_rest_Enrichr")
run_enrichr_analysis(all_deseq_res_pos[(all_deseq_res_pos$cluster=="Extremities")&(all_deseq_res_pos$padj<0.1),]$gene,selected_dbs,
                     "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_extremities_vs_rest_Enrichr")
run_enrichr_analysis(all_deseq_res_pos[(all_deseq_res_pos$cluster=="Heart")&(all_deseq_res_pos$padj<0.1),]$gene,selected_dbs,
                     "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_heart_vs_rest_Enrichr")
run_enrichr_analysis(all_deseq_res_pos[(all_deseq_res_pos$cluster=="Others")&(all_deseq_res_pos$padj<0.1),]$gene,selected_dbs,
                     "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_others_vs_rest_Enrichr")

run_enrichr_analysis(all_hclust_deseq_res_pos[(all_hclust_deseq_res_pos$cluster=="Cluster 1")&(all_hclust_deseq_res_pos$padj<0.1),]$gene,selected_dbs,
                     "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_hclust1_vs_rest_Enrichr")
run_enrichr_analysis(all_hclust_deseq_res_pos[(all_hclust_deseq_res_pos$cluster=="Cluster 2")&(all_hclust_deseq_res_pos$padj<0.1),]$gene,selected_dbs,
                     "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_hclust2_vs_rest_Enrichr")
run_enrichr_analysis(all_hclust_deseq_res_pos[(all_hclust_deseq_res_pos$cluster=="Cluster 3")&(all_hclust_deseq_res_pos$padj<0.1),]$gene,selected_dbs,
                     "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_hclust3_vs_rest_Enrichr")
run_enrichr_analysis(all_hclust_deseq_res_pos[(all_hclust_deseq_res_pos$cluster=="Cluster 4")&(all_hclust_deseq_res_pos$padj<0.1),]$gene,selected_dbs,
                     "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_hclust4_vs_rest_Enrichr")
run_enrichr_analysis(all_hclust_deseq_res_pos[(all_hclust_deseq_res_pos$cluster=="Cluster 5")&(all_hclust_deseq_res_pos$padj<0.1),]$gene,selected_dbs,
                     "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_hclust5_vs_rest_Enrichr")




## Helper function to run fast-ORA
run_fora = function(degs,fgsea_sets,pval_col="p_val") {
  #genes = degs[degs$p_val_adj < 0.1,]$gene
  genes = degs[degs[[pval_col]] < 0.05,]$gene
  #print(genes)
  universe = rownames(so@assays$RNA$counts)
  fora_res = fora(fgsea_sets, genes, universe, minSize = 5, maxSize = 500)
  return(fora_res)
}

## Helper function to prettify fast-ORA results
process_fora_res = function(raw_res,title) {
  raw_res$cluster = title
  raw_res$pathway_minimized = gsub("GOBP_","",raw_res$pathway)
  raw_res$pathway_minimized = gsub("_"," ",raw_res$pathway_minimized)
  raw_res$neglog_pval = -log(raw_res$pval)
  raw_res$neglog_padj = -log(raw_res$padj)
  raw_res$enrichment_ratio = raw_res$overlap/raw_res$size
  return(raw_res)
}

## Plot one with bars instead of scatter
plot_fora_bar = function(res) {
  fora_barplot = ggplot(res,aes(x=neglog_pval,y=fct_reorder(pathway_minimized,neglog_pval))) + 
    geom_bar(stat="identity") +
    theme_minimal() +
    labs(x="-log(p-value)",y="GO-BP Term") 
  return(fora_barplot)
}

## Helper functionfor running cluster enrichments
run_fora_analysis = function(markers_df, fgsea_set, ylabel, filename, n_top = 10, pval_col="p_val") {
  unique_clusters = unique(markers_df[["cluster"]])
  plot_list = list()
  result_list = list()
  for (idx in 1:length(unique_clusters)) {
    cluster = unique_clusters[[idx]]
    # Subset markers for the current cluster
    cluster_markers = markers_df[markers_df[["cluster"]] == cluster,]
    fora_res = run_fora(cluster_markers, fgsea_set,pval_col=pval_col)
    processed_res = process_fora_res(fora_res, cluster)
    processed_res$term_set = ylabel
    processed_res$cluster = cluster
    processed_res_to_show = processed_res %>% top_n(n_top, wt = neglog_pval)
    result_list[[idx]] = as.data.frame(processed_res)
    # Plot results
    plot = plot_fora_bar(processed_res_to_show) + 
      labs(title =cluster, y = ylabel)
    plot_list[[idx]] = plot
  }
  # Combine all plots into a single grid
  combined_plot = plot_grid(plotlist = plot_list, nrow = length(unique_clusters))
  ggsave(filename,combined_plot,dpi=300,width=10,height=24)
  return(bind_rows(result_list))
}

#new_fgsea_hallmark_set = fgsea_hallmark_set
#new_fgsea_hallmark_set$EndothelialMesenchymalTransition=c("TGFB1","TGFB2","AIFM2","SNAI1","SNAI2","ZEB1","ZEB2","ACTA2")
## Maybe instead of FORA we just do module scores instead

## Run ORA analysis for a bunch of gene sets (DESEQ2 Sites)
all_deseq_res_pos = read.csv("02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_sites_combined_results_pos.tsv", sep = "\t")

## Run FORA analysis for a bunch of gene sets (DESEQ2 Sites)
deseq2_site_hallmark_fora = run_fora_analysis(all_deseq_res_pos,fgsea_hallmark_set,"Hallmarks","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_hallmark_fora_bar_by_DESEQ2_site_clusters.png",pval_col="pvalue")
deseq2_site_kegg_fora = run_fora_analysis(all_deseq_res_pos,fgsea_kegg_set,"KEGG","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_KEGG_fora_bar_by_DESEQ2_site_clusters.png",pval_col="pvalue")
deseq2_site_c1_fora = run_fora_analysis(all_deseq_res_pos,fgsea_c1_set,"C1","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C1_fora_bar_by_DESEQ2_site_clusters.png",pval_col="pvalue")
deseq2_site_c5_fora = run_fora_analysis(all_deseq_res_pos,fgsea_c5_set,"C5","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C5_fora_bar_by_DESEQ2_site_clusters.png",pval_col="pvalue")
deseq2_site_c6_fora = run_fora_analysis(all_deseq_res_pos,fgsea_c6_set_up_only,"C6","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C6_fora_bar_by_DESEQ2_site_clusters.png",pval_col="pvalue")
deseq2_site_c8_fora = run_fora_analysis(all_deseq_res_pos,fgsea_c8_set,"C8","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C8_fora_bar_by_DESEQ2_site_clusters.png",pval_col="pvalue")
deseq2_combined_fora = rbind(deseq2_site_hallmark_fora,deseq2_site_kegg_fora,deseq2_site_c1_fora,deseq2_site_c5_fora,deseq2_site_c6_fora,deseq2_site_c8_fora)
deseq2_combined_fora$overlapGenes = sapply(deseq2_combined_fora$overlapGenes,function(x) paste0(x,collapse=","))
## Write down FORA results
write.table(deseq2_combined_fora,file = "02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_sites_combined_FORA_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

## Run FORA analysis for a bunch of gene sets (hclust clusters)
hclust_res_pos = read.csv("02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_hclust_expr_sites_combined_results_pos.tsv", sep = "\t")
hclust_hallmark_fora = run_fora_analysis(hclust_res_pos,fgsea_hallmark_set,"Hallmarks","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_hallmark_fora_bar_by_expr_hclust.png",pval_col="pvalue")
hclust_kegg_fora = run_fora_analysis(hclust_res_pos,fgsea_kegg_set,"KEGG","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_KEGG_fora_bar_by_expr_hclust.png",pval_col="pvalue")
hclust_c1_fora = run_fora_analysis(hclust_res_pos,fgsea_c1_set,"C1","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C1_fora_bar_by_expr_hclust.png",pval_col="pvalue")
hclust_c5_fora = run_fora_analysis(hclust_res_pos,fgsea_c5_set,"C5","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C5_fora_bar_by_expr_hclust.png",pval_col="pvalue")
hclust_c6_fora = run_fora_analysis(hclust_res_pos,fgsea_c6_set_up_only,"C6","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C6_fora_bar_by_expr_hclust.png",pval_col="pvalue")
hclust_c8_fora = run_fora_analysis(hclust_res_pos,fgsea_c8_set,"C8","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C8_fora_bar_by_expr_hclust.png",pval_col="pvalue")
hclust_combined_fora = rbind(hclust_hallmark_fora,hclust_kegg_fora,hclust_c1_fora,hclust_c5_fora,hclust_c6_fora,hclust_c8_fora)
hclust_combined_fora$overlapGenes = sapply(hclust_combined_fora$overlapGenes,function(x) paste0(x,collapse=","))
## Write down FORA results
write.table(hclust_combined_fora,file = "02_Tumor_RNA_Analysis/outputs/DEGs/02_expr_hclust_combined_FORA_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

## (Update May 12 2025)
all_deseq_res_pos


plot_combined_fora = function(fora_df,term_set_filter,pathway_col="pathway_minimized") {
  top_hallmark_foras = fora_df %>% filter(term_set==term_set_filter) %>%
    group_by(cluster) %>% 
    arrange(pval, .by_group = TRUE) %>%
    filter(pval < 0.05) %>%
    slice_head(n = 5) %>%
    ungroup() %>%
    mutate(unique_pathway = paste(cluster,pathway_minimized, sep = ": ")) %>%
    group_by(cluster) %>%
    mutate(unique_pathway = factor(unique_pathway, levels = unique_pathway[order(neglog_pval)])) %>%
    ungroup() %>%
    mutate(significance=ifelse(padj < 0.05, "adjusted p < 0.05","nominal p < 0.05"))# %>%
  significance_colors = c("adjusted p < 0.05" = "blue", "nominal p < 0.05" = "lightblue")
  fora_barplot = ggplot(top_hallmark_foras,aes(x=neglog_pval,y=unique_pathway,fill = significance)) + 
    geom_bar(stat="identity") +
    scale_fill_manual(values = significance_colors) +
    theme_minimal() +
    labs(x="-log(p-value)",fill="Significance", y=term_set_filter) +
    facet_wrap(~ cluster, ncol = 1,scales="free_y") +
    theme(
      strip.text = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom"
    )
  return(fora_barplot)
}

#hclust_combined_fora$pretty_pathway = gsub("\\HALLMARK","",hclust_combined_fora$pathway_minimized)

hclust_hallmark_combined_fora = plot_combined_fora(hclust_combined_fora,"Hallmarks") + labs(y="MsigDB Hallmarks") + coord_cartesian(xlim = c(0, 50))
hclust_c6_combined_fora = plot_combined_fora(hclust_combined_fora,"C6") + labs(y="MsigDB Oncogenic Signature (C6)") + coord_cartesian(xlim = c(0, 50))
ggsave("02_Tumor_RNA_Analysis/outputs/plots/04_enrichment_plots/02_hclust_pathways_hallmarks_combined.pdf",hclust_hallmark_combined_fora,width=10,height=8)
ggsave("02_Tumor_RNA_Analysis/outputs/plots/04_enrichment_plots/02_hclust_pathways_c6_combined.pdf",hclust_c6_combined_fora,width=10,height=8)

deseq2_sites_hallmark_combined_fora = plot_combined_fora(deseq2_combined_fora,"Hallmarks") + labs(y="MsigDB Hallmarks") + coord_cartesian(xlim = c(0, 50))# + labs(y="C6")
deseq2_sites_c6_combined_fora = plot_combined_fora(deseq2_combined_fora,"C6") + labs(y="MsigDB Oncogenic Signature (C6)") + coord_cartesian(xlim = c(0, 50))# + labs(y="C6")
ggsave("02_Tumor_RNA_Analysis/outputs/plots/04_enrichment_plots/02_deseq2_sites_pathways_hallmarks_combined.pdf",deseq2_sites_hallmark_combined_fora,width=10,height=8)
ggsave("02_Tumor_RNA_Analysis/outputs/plots/04_enrichment_plots/02_deseq2_sites_pathways_c6_combined.pdf",deseq2_sites_c6_combined_fora,width=10,height=8)


## Test Showing enrichment of all pathways at the same time
## Dot-plot of Most commonly enriched pathways
deseq2_combined_fora_C6 = deseq2_combined_fora[deseq2_combined_fora$term_set=="Hallmarks",]
deseq2_combined_fora_C6$padj_sig = deseq2_combined_fora_C6$padj < 0.1
test_plot = ggplot(deseq2_combined_fora_C6, aes(x = cluster, y = pathway, size = neglog_pval, color=padj_sig)) +
  geom_point(alpha = 0.7) +  # Blue dots with transparency
  scale_size(range = c(2, 10)) +  # Adjust dot size range
  scale_color_npg() +
  theme_minimal() +
  labs(
    title = "Pathway Enrichment Across Clusters",
    x = "Cluster",
    y = "Pathway",
    size = "-log10(pval)",
    color = "Significant (padj<0.1)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_fora_Hallmarks_pathway_enrichment_dot_plot.png",test_plot,dpi=300,height=10,width = 10)
ggsave("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_fora_Hallmarks_pathway_enrichment_dot_plot.pdf",test_plot,dpi=300,height=10,width = 10)


## Run FORA analysis for a bunch of gene sets (expr clusters)
# run_fora_analysis(all_markers_by_expr_cluster,fgsea_hallmark_set,"Hallmarks","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_hallmark_fora_bar_by_expr_clusters.png")
# run_fora_analysis(all_markers_by_expr_cluster,fgsea_kegg_set,"KEGG","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_KEGG_fora_bar_by_expr_clusters.png")
# run_fora_analysis(all_markers_by_expr_cluster,fgsea_c1_set,"C1","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C1_fora_bar_by_expr_clusters.png")
# run_fora_analysis(all_markers_by_expr_cluster,fgsea_c5_set,"C5","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C5_fora_bar_by_expr_clusters.png")
# run_fora_analysis(all_markers_by_expr_cluster,fgsea_c6_set_up_only,"C6","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C6_fora_bar_by_expr_clusters.png")
# run_fora_analysis(all_markers_by_expr_cluster,fgsea_c8_set,"C8","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C8_fora_bar_by_expr_clusters.png")

run_fora_analysis(all_markers_by_expr_cluster,fgsea_hallmark_set,"Hallmarks","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_hallmark_fora_bar_by_expr_clusters.pdf")
run_fora_analysis(all_markers_by_expr_cluster,fgsea_kegg_set,"KEGG","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_KEGG_fora_bar_by_expr_clusters.pdf")
run_fora_analysis(all_markers_by_expr_cluster,fgsea_c1_set,"C1","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C1_fora_bar_by_expr_clusters.pdf")
run_fora_analysis(all_markers_by_expr_cluster,fgsea_c5_set,"C5","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C5_fora_bar_by_expr_clusters.pdf")
run_fora_analysis(all_markers_by_expr_cluster,fgsea_c6_set_up_only,"C6","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C6_fora_bar_by_expr_clusters.pdf")
run_fora_analysis(all_markers_by_expr_cluster,fgsea_c8_set,"C8","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C8_fora_bar_by_expr_clusters.pdf")

## Run FORA analysis for a bunch of gene sets (site clusters)
# run_fora_analysis(all_markers_by_site_cluster,fgsea_hallmark_set,"Hallmarks","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_hallmark_fora_bar_by_site_clusters.png")
# run_fora_analysis(all_markers_by_site_cluster,fgsea_kegg_set,"KEGG","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_KEGG_fora_bar_by_site_clusters.png")
# run_fora_analysis(all_markers_by_site_cluster,fgsea_c1_set,"C1","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C1_fora_bar_by_site_clusters.png")
# run_fora_analysis(all_markers_by_site_cluster,fgsea_c5_set,"C5","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C5_fora_bar_by_site_clusters.png")
# run_fora_analysis(all_markers_by_site_cluster,fgsea_c6_set_up_only,"C6","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C6_fora_bar_by_site_clusters.png")
# run_fora_analysis(all_markers_by_site_cluster,fgsea_c8_set,"C8","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C8_fora_bar_by_site_clusters.png")

run_fora_analysis(all_markers_by_site_cluster,fgsea_hallmark_set,"Hallmarks","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_hallmark_fora_bar_by_site_clusters.pdf")
run_fora_analysis(all_markers_by_site_cluster,fgsea_kegg_set,"KEGG","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_KEGG_fora_bar_by_site_clusters.pdf")
run_fora_analysis(all_markers_by_site_cluster,fgsea_c1_set,"C1","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C1_fora_bar_by_site_clusters.pdf")
run_fora_analysis(all_markers_by_site_cluster,fgsea_c5_set,"C5","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C5_fora_bar_by_site_clusters.pdf")
run_fora_analysis(all_markers_by_site_cluster,fgsea_c6_set_up_only,"C6","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C6_fora_bar_by_site_clusters.pdf")
run_fora_analysis(all_markers_by_site_cluster,fgsea_c8_set,"C8","02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_C8_fora_bar_by_site_clusters.pdf")


## For heatmap (Fig 2.) Highlight genes from top enriched sets
## Something to include (TGFB Up, PTEN DN, MYC UP, PGF Up, Cyclin D1 Up,) 
all_markers_by_site_cluster_nomsig = all_deseq_res_pos[all_deseq_res_pos$padj < 0.05,]
test1 = all_markers_by_site_cluster_nomsig[all_markers_by_site_cluster_nomsig$cluster=="ParenchymalBreast",]$gene
test2 = fgsea_c6_set$PTEN_DN.V1_UP
intersect(test1,test2)

