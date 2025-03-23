library(Seurat)
library(cowplot)
library(ggplot2)
library(tidyr)
library(dplyr)
#library(fgsea)
#library(msigdbr)
library(ggsci)
library(dendextend)
#library(dynamicTreeCut)
#library(harmony)
#library(ComplexHeatmap)
#library(ggrepel)
#library(circlize)
#library(RColorBrewer)
#library(dendextend)
#library(forcats)
#library(ggpubr)
#library(ggbeeswarm)

setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts")

## Load count data and reformat IDs for downstream compatibility
count_data = read.csv("data/processed/rna/00_filtered_gene_counts.csv",row.names = 1,check.names = FALSE)
vst_data = read.csv("data/processed/rna/00_filtered_gene_counts.vst.csv",row.names = 1,check.names = FALSE)
tpm_data = read.csv("data/processed/rna/00_filtered_gene_tpm.csv",row.names = 1,check.names = FALSE)
  
colnames(count_data) = gsub("-", ".", colnames(count_data))
colnames(count_data) = gsub("_", ".", colnames(count_data))
colnames(vst_data) = gsub("-", ".", colnames(vst_data))
colnames(vst_data) = gsub("_", ".", colnames(vst_data))
colnames(tpm_data) = gsub("-", ".", colnames(tpm_data))
colnames(tpm_data) = gsub("_", ".", colnames(tpm_data))

# Store a scaled version of the data row-wise (gene) for visualization
vst_data_scaled = as.data.frame(t(scale(t(vst_data))),check.names=FALSE)
vst_data_unscaled = as.data.frame((vst_data),check.names=FALSE)

#vst_data_scaled = as.data.frame(t(scale(t(vst_data))),check.names=FALSE)

#vst_data_scaled_pretty_names = gsub("\\.","_",colnames(vst_data_scaled))
#colnames(vst_data_scaled) =  vst_data_scaled_pretty_names

## Apply same formatting to metadata and reorder it
meta_data = read.csv("data/processed/rna/00_filtered_sample_metadata.csv",check.names = FALSE)
meta_data$seurat_id = meta_data$`entity:sample_id`
meta_data$seurat_id = gsub("_", ".", meta_data$seurat_id)
meta_data$seurat_id = gsub("-", ".", meta_data$seurat_id)
meta_data[meta_data == ""] = NA
rownames(meta_data) = meta_data$seurat_id

# Create a seurat object and add meta data
#so = CreateSeuratObject(counts = vst_data,min.cells = 0, min.features = 0, meta.data = meta_data)
so = CreateSeuratObject(counts = count_data,min.cells = 0, min.features = 0, meta.data = meta_data)
#so = CreateSeuratObject(counts = tpm_data,min.cells = 0, min.features = 0, meta.data = meta_data)


so@meta.data[["og_id"]] = rownames(so@meta.data)
so@meta.data = merge(so@meta.data,meta_data,by.x="og_id",by.y="seurat_id",all.x=TRUE)
rownames(so@meta.data) = so@meta.data[["og_id"]]

so@assays$RNA$vst_scaled = vst_data_scaled
so@assays$RNA$vst_data = vst_data_unscaled

so = SetIdent(so, value = rownames(so@meta.data))

# Adapted based on rocanja's recommendation in handling bulk data in seurat
# reference: https://github.com/satijalab/seurat/issues/7496
#so = SetAssayData(object=so, layer="data", new.data = vst_data)
so = SetAssayData(object=so, layer="data", new.data = log1p(tpm_data))
so = FindVariableFeatures(so,selection.method = "vst", nfeatures = 2000)
#so = ScaleData(so)
## Regress out estimated purity
so@assays$RNA$data = as.matrix(so@assays$RNA$data)
so = ScaleData(so,vars.to.regress = "ESTIMATE_purity",features = rownames(so),assay = "RNA")
## Dimension reduction
so = RunPCA(so, features = VariableFeatures(object = so))
## Run unsuperivsed clustering
so = FindNeighbors(so, dims = 1:10)
so = FindClusters(so, resolution = 1.5)

## For easier reading start the unsupervised cluster index at 1 instead of 0
so@meta.data$seurat_clusters_by_expression = as.integer(so@meta.data$seurat_clusters)
so@meta.data$seurat_clusters_by_expression_str = paste0("Cluster ",so@meta.data$seurat_clusters_by_expression)

## Also define "clusters" by site information only
so@meta.data = so@meta.data %>% mutate(
  seurat_clusters_by_site = case_when(
    `Primary Site (Recombined)` == "Breast (Parenchymal)" ~ 1,
    `Primary Site (Recombined)` == "Breast (Cutaneous)" ~ 2,
    `Primary Site (Recombined)` == "HNFS" ~ 3,
    `Primary Site (Recombined)` == "Heart" ~ 4,
    `Primary Site (Recombined)` == "Extremities" ~ 5,
    TRUE ~ 6
  ),
  seurat_clusters_by_site_str = case_when(
    `Primary Site (Recombined)` == "Breast (Parenchymal)" ~ "Breast (Parenchymal)",
    `Primary Site (Recombined)` == "Breast (Cutaneous)" ~ "Breast (Cutaneous)",
    `Primary Site (Recombined)` == "HNFS" ~ "HNFS",
    `Primary Site (Recombined)` == "Heart" ~ "Heart",
    `Primary Site (Recombined)` == "Extremities" ~ "Extremities",
    TRUE ~ "Others"
  ),
)

## Also compute hierarchical clusters based on vst data (Using variable genes)
genes_to_use_for_hclust = so@assays$RNA$vst_data[VariableFeatures(so,nfeatures = 2000),]
#genes_to_use_for_hclust = so@assays$RNA$vst_data#[VariableFeatures(so),]
hclust_dist_matrix = dist(t(genes_to_use_for_hclust), method = "euclidean")
hclust_res = hclust(hclust_dist_matrix, method = "ward.D2")

## Assign clusters via DynamicCut
sample_clusters = cutreeDynamic(dendro = hclust_res,minClusterSize = 10)
## Asign clusters via cut height
sample_clusters = cutree(hclust_res, k = NULL, h = 200)

# Get sample labels and primary site information
sample_labels = colnames(genes_to_use_for_hclust)
primary_sites = so@meta.data[sample_labels, ]$seurat_clusters_by_site_str#`Primary Site (Recombined)`
unique_clusters = unique(sample_clusters)

# Convert hclust object to dendrogram
dend = as.dendrogram(hclust_res)
dend = color_branches(dend, k = length(unique_clusters), groupLabels = as.character(unique_clusters))

# Color the leaves by primary site
# Match primary site colors to the corresponding sample labels
dend_order = order.dendrogram(dend)
ordered_labels = primary_sites[dend_order]
leaf_colors = cluster_by_primary_site_palette[ordered_labels]
names(leaf_colors) = sample_labels
dend = set(dend, "labels_col", leaf_colors)
dend = set(dend, "labels_cex", 0.8)
dend = set(dend, "labels", ordered_labels)  # Label by primary site

# Plot colored dendrogram
pdf("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_gene_expression_cluster_dendrogram.pdf",width=18,height=14)
plot(dend, main = "Sample Clustering Dendrogram by Primary Site", ylab = "Height")
# Add legend for primary sites
legend("topright", legend = unique(primary_sites), 
       fill = cluster_by_primary_site_palette, 
       title = "Primary Site",
       cex = 0.8)
dev.off()

# Add cluster assignments to metadata
so@meta.data$hcluster_by_expr = as.factor(sample_clusters)
so@meta.data$hcluster_by_expr_str = paste0("Cluster ",so@meta.data$hcluster_by_expr)

## Plot cluster membership by site
cluster_membership_table = table(so@meta.data$hcluster_by_expr_str,so@meta.data$seurat_clusters_by_site_str)
#cluster_membership_table_hmap = pheatmap(cluster_membership_table, display_numbers = T)
pdf("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_cluster_gene_set_heatmap_dendro_table.pdf",width=10,height=10)
plot(cluster_membership_table_hmap)
dev.off()



## Infer sample sex by XIST expression
#VlnPlot(so,c("XIST"),group.by = "SEX (EHR_EXTRACTED)",slot = "count")
xist_expression = FetchData(so,vars="XIST",slot = "counts")
so@meta.data$xist_inferred_female = ifelse(xist_expression >= 1000, 1, 0)
so = AddMetaData(so, metadata = so$xist_inferred_female, col.name = "xist_inferred_female")
so@meta.data = so@meta.data %>% mutate(rna_inferred_female = ifelse(!is.na(`SEX (EHR_EXTRACTED)`), ifelse(`SEX (EHR_EXTRACTED)` == "F", 1, 0), xist_inferred_female))
## For cluster names start index at 1 instead of 0 for easier reading
#so@meta.data$seurat_clusters_renamed = as.integer(so@meta.data$seurat_clusters)
#so@meta.data$seurat_clusters_renamed = as.integer(so@meta.data$seurat_clusters_by_site)
#so@meta.data$seurat_clusters_renamed_str = paste0("Cluster ",so@meta.data$seurat_clusters_renamed)

### This section adds the module scores ###
### Does not work at the current moment ###
so = JoinLayers(so)
## Add radioresistance score
radioresistance_genes_df = read.csv("data/public/Marcone2021_Radioresistance_Geneset.csv",check.names = FALSE)
radioresistance_genes = radioresistance_genes_df[radioresistance_genes_df$`Fold change` > 1,]$`Gene names`
so = AddModuleScore(so, list(radioresistance_genes),name="radioresistance_score",slot="data")

## Add radiation exposure score
radiation_exposure_genes_df = read.csv("data/public/Paul2013_RadiationExposure_Geneset.csv",check.names = FALSE)
radiation_exposure_genes = radiation_exposure_genes_df$Symbol
so = AddModuleScore(so, list(radiation_exposure_genes),name="radiation_exposure_score",slot="data")

## Add angiogenesis score
## TODO: retrospectively add gene set here
fgsea_hallmark_set = msigdbr(species = "Homo sapiens", category = "H") %>% split(x = .$gene_symbol, f = .$gs_name)
so = AddModuleScore(so, list(fgsea_hallmark_set$HALLMARK_ANGIOGENESIS),name="angiogenesis_score",slot="data")

## Add Lympoangiogenesis score
fgsea_c2_set = msigdbr(species = "Homo sapiens", category = "C2")  %>% split(x = .$gene_symbol, f = .$gs_name)
lymphoangiogenesis_set = fgsea_c2_set$KEGG_CHRONIC_MYELOID_LEUKEMIA #c("IL4", "CSF2", "PROX1", "TEK")
so = AddModuleScore(so, list(lymphoangiogenesis_set),name="lymphangiogenesis_score",slot="data")

## Add EndMT score
## From Choi K. et al: Endothelial-to-mesenchymal transition in anticancer therapy and normal tissue damage
EndMT_genes = c("TGFB1","TGFB2","AIFM2","SNAI1","SNAI2","ZEB1","ZEB2","ACTA2")
so = AddModuleScore(so, list(EndMT_genes),name="endMT_score",slot="data")

## Add fibroblast score
fibroblast_genes = c("SCN1A","DPT","HSD11B1","CD34","CELF2","ENTPD2","COL5A1","GSN","SLC43A3","S100A16","S100A10","OLFM3","LPAR1","HTRA3","UGDH","PDGFRA","C1QTNF3-AMACR","MEDAG","COL1A2","FBLN2","MFAP5","MGST1","LSP1","MMP2","DPEP1","LOXL1","PCOLCE2","BICC1","DCN","LUM","GFPT2","ADAMTS2","MFAP4","SERPINF1","COL1A1","ABCA8","RNASE")
so = AddModuleScore(so, list(fibroblast_genes),name="fibroblast_score",slot="data")

## Add EMT score
emt_score = read.csv("data/public/Chakravarthy2018_TGFB_ECM.tsv",sep="\t")
emt_genes = emt_score$ID
so = AddModuleScore(so, list(emt_genes),name="EMT_score",slot="data")

# ## Add fibrosis score
# fibrosis_genes = c("LGALS3")
# temp_plot = VlnPlot(so,features = c("MUC5B","GPR",""))
# ggsave(filename="02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_MUC5B_expr.png",temp_plot)

#tgf_beta_induced_genes = c("GSN","MYLK2","SOX9")

## Save Seurat Object
if (FALSE) {
  SaveSeuratRds(so,file = "data/processed/rna/ASCSeuratObj2025.rds")
  #SaveSeuratRds(so,file = "data/processed/rna/ASCSeuratObj2025_Purity.rds")
}

## Load Seurat Object for basic sumary statistics

so = readRDS("data/processed/rna/ASCSeuratObj2025.rds")

## Number of parencymal breast samples by site
site_by_cluster_table = table(so@meta.data$seurat_clusters_by_site_str,so@meta.data$hcluster_by_expr_str)
site_totals = rowSums(site_by_cluster_table)
cluster_totals = colSums(site_by_cluster_table)


## Plot scores by clusters
score_subset = so@meta.data[,c("radioresistance_score1","radiation_exposure_score1",
                               "angiogenesis_score1","endMT_score1","fibroblast_score1",
                               "EMT_score1"
                               )]
score_subset$og_id = rownames(score_subset)
score_subset_long = pivot_longer(score_subset,cols = -c("og_id"))
score_subset_merged = merge(score_subset_long,so@meta.data,by="og_id",all.x=TRUE)

score_labeller = labeller(name=c("radioresistance_score1"="Radioresistance (Marcone 2021)",
                                 "radiation_exposure_score1"="Radiation Exposure (Paul 2011)",
                                 "angiogenesis_score1"="Angiogenesis (MSigDB Hallmark)",
                                 "endMT_score1"="EndMT Score (Choi 2020)",
                                 "fibroblast_score1"="Fibroblast Score (Muhl 2020)",
                                 "EMT_score1"="ECM Score (Chakravarthy 2018)"
                                 ))
## Wrap as two rows
scores_boxplot = ggplot(score_subset_merged, 
                        aes(x = `seurat_clusters_renamed_str`, y = value)) +
  geom_violin(width = 1.2) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  facet_wrap(facets = vars(name),nrow=2,labeller=score_labeller) +
  #facet_grid(cols = vars(name),labeller=score_labeller) +
  # stat_compare_means(comparisons = list(c("Cluster 1", "Cluster 2"),
  #                                       c("Cluster 1", "Cluster 3"),
  #                                       c("Cluster 2", "Cluster 3"),
  #                                       c("Cluster 2", "Cluster 4")
  # ),
  # na.rm = TRUE, label = "p.format") +
  stat_compare_means(method = "anova", label.y = 2.2) +
  geom_quasirandom(aes(color = `Primary Site (Recombined)`)) +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0),text = element_text(size=22)) + # Adjust facet label orientation
  labs(x = "Cluster", y = "Normalized Expression Level", color = "Cluster") +
  guides(color=guide_legend(title="Primary Site"))
ggsave(scores_boxplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_scores_boxplot_with_endMT_fibroblast_ECM.png",dpi=300,width=20,height=16)


#VlnPlot(so,features = c("angiogenesis_score1","lymphangiogenesis_score1","radioresistance_score1","radiation_exposure_score1"), pt.size =1,group.by = "seurat_clusters_renamed")

## Find markers associated with each cluster
Idents(so) = so@meta.data$seurat_clusters_renamed
all_markers = FindAllMarkers(
  so,
  #slot="counts",
  test.use = "wilcox",
  #fc.slot = NULL,
  only.pos = TRUE,
  #min.pct = 0,
  #min.cells.group=0,
  #min.cells.feature = 0,
  logfc.threshold = 0,
  return.thresh = 1.1
)

# all_markers_LR = FindAllMarkers(
#   so,
#   slot="counts",
#   test.use = "LR",
#   only.pos = TRUE
# )

# all_markers_incl_neg = FindAllMarkers(
#   so,
#   slot="counts",
#   test.use = "wilcox",
#   only.pos = FALSE
# )

## Save the output
# write.csv(all_markers,"02_Tumor_RNA_Analysis/outputs/DEGs/02_Seurat_all_markers.csv")
#write.csv(all_markers_LR,"02_Tumor_RNA_Analysis/outputs/DEGs/02_Seurat_all_markers_LR.csv")



## Check the top genes in each cluster
top_cluster_markers = all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val < 0.01) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  #dplyr::filter(!grepl("KRT",gene)) %>%
  slice_head(n = 10) %>%
  ungroup()

top_cluster_markers_heatmap = DoHeatmap(so, features = top_cluster_markers$gene,slot="vst.scaled",disp.max = 2.5,disp.min = -2.5)
top_cluster_markers_heatmap
ggsave(top_cluster_markers_heatmap,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_top_cluster_marker_heatmap.png",dpi=300,width=24,height=12)

## overlap the top 5 in the PCs
top_5_cluster_markers = all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.1) %>%
  dplyr::filter(avg_log2FC > 2) %>%
  #dplyr::filter(!grepl("KRT",gene)) %>%
  slice_head(n = 10) %>%
  ungroup()

## The same content but on PCs
top_markers_dimplot = FeaturePlot(so, features = top_5_cluster_markers$gene,pt.size = 2,ncol = 5,reduction = "pca")
top_markers_dimplot
ggsave(top_markers_dimplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_top_5_markers_per_cluster_on_pca_plot.png",dpi=300,width=24,height=12)

## Save the Seurat Object
if (FALSE) {
  saveRDS(so, "data/processed/rna/ASCSeuratObj.rds")
}

## Check genes that are shared between cluster 3 and 4
clust3_samples = rownames(so@meta.data[so@meta.data$seurat_clusters_renamed==3,])
clust4_samples = rownames(so@meta.data[so@meta.data$seurat_clusters_renamed==4,])
clust34_subset = subset(x = so, subset = seurat_clusters_renamed == 3 | seurat_clusters_renamed == 4)
clust34_markers = FindMarkers(
  clust34_subset,ident.1 = 3,ident.2 = 4,
  slot="counts",
  test.use = "wilcox"
)
clust34_marker_df_sig = clust34_markers %>% filter(p_val_adj < 0.1)
top_clust34_fc_pos_genes = clust34_marker_df_sig %>% arrange(-avg_log2FC) %>% head(10) %>% rownames()
top_clust34_fc_neg_genes = clust34_marker_df_sig %>% arrange(avg_log2FC) %>% head(10) %>% rownames()
top_clust34_pval_genes = clust34_marker_df_sig %>% arrange(p_val_adj) %>% head(10) %>% rownames()
clust34_genes_to_highlight = c(top_clust34_fc_pos_genes,top_clust34_fc_neg_genes,top_clust34_pval_genes,c("FLT3","MYCL","KRT17"))
clust34_vol_plot = EnhancedVolcano(
  clust34_markers,
  lab=rownames(clust34_markers),
  x="avg_log2FC",
  y="p_val_adj",
  selectLab = clust34_genes_to_highlight,
  pCutoff = 10e-2,
  FCcutoff = 0.5,
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  title = "Cluster 3 vs. Cluster 4",
  subtitle = "Differential Expression"
) + xlab("log2(Cluster 3/Cluster 4)")
ggsave("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_clust3_vs_clust4_volcano.png",clust34_vol_plot,width = 12,height = 8)

  # 
# all_markers %>% filter(cluster == 0, grepl("CD", gene))
# m3 = all_markers[all_markers$cluster==3,] %>% filter(p_val < 0.05)
# m4 = all_markers[all_markers$cluster==4,] %>% filter(p_val < 0.05)
# inm3_notin_m4 = m3[!(m3$gene %in% m4$gene),] %>%
#   arrange(p_val_adj) %>%
#   filter(p_val < 0.05)
# inm4_notin_m3 = m4[!(m4$gene %in% m3$gene),] %>%
#   arrange(p_val_adj) %>%
#   filter(p_val < 0.05)
# inm3_notin_m4

# ## Plot the ASC markers combining clinical and new features
# angiogenesis_markers = c(
#   #"VWF","ERG","FLI1",
#   #"CD36","CD93"
#   #"ITGA5","FIGF"
#   #"TGFBR1","PIGF",
#   #"VEGFC"
#   #"KDR","PIK3CA","PGF"
#   #"STAT1"
#   "FLT1","NRP1","KDR","APLN","APLNR",#"PDGFA","PDGFRA",
#   "PDGFB","PDGFRB","FGF1","FGFR1","CD34"#"TBX2","TBX3","NOTCH4", #"HFE",
#   #,"PECAM1","PECAM-1"
# )
# 
# angiogenesis_markers = head(all_markers[all_markers$cluster==1,]$gene,10)
# 
# angiogenesis_df = data.frame(row.names = angiogenesis_markers)
# angiogenesis_df$`Gene Set` = "Angiogenesis"
# 
# # stem_cell_proliferation_markers = c(
# #   "CD34","FERMT2"
# #   #"TGFB1","TGFBR1","CCNE1","SIX2"
# #   #"TERT"
# #   #"FGF2","PTPRC","NANOG","WNT1","WNT3"
# # )
# # stem_cell_df = data.frame(row.names = stem_cell_proliferation_markers)
# # stem_cell_df$`Gene Set` = "Stem Cell Proliferation"
# 
# lymphoangiogensis_markers = c(
#   #"CD320",
#   #"LYVE1",
#   "FLT4","NRP2","PROX1","PDPN","MYC", #"PECAM1",
#   "TBX1","TPX2","TIE1","MAP4K2","CTLA4"
#   #"CD274"
#   #"CLEC14A",#"ANGPT2",
#   #"CCBE1","EPHA2","VEGFC","VASH1"
# )
# 
# lymphoangiogensis_markers = head(all_markers[all_markers$cluster==2,]$gene,10)
# 
# lymphoangiogensis_df = data.frame(row.names = lymphoangiogensis_markers)
# lymphoangiogensis_df$`Gene Set` = "Lymphangiogenesis"
# 
# flt3_related_markers = c(
#   "FGFR2","FGFR3","FGF22","FGFBP1",
#   "FLT3","FLT3LG","MYCL",
#   "IL37","IL18",#"CD44","CDH1",
#   "LMO1","WNT7B","FBXW7","CXCL12","DKK2",
#   #"KDR",
#   #"VCAM1",
#   #"FGF10"
#   #"KRT1","FLG",
#   #"KRT5","SERPINB5","GATA3", #"CASP14","TP73","CDH3",
#   #"CLEC2A","KLRF2"
#   #,"KRT10"
#   #"KRT1"#"KRT5"#"KRT17","KRT24",
#   "KRT1","KRT2","KRT5","KRT10","FLG"
#   #"KRT14","KRT15","KRT77","KRT78","KRT80"
# )
# 
# flt3_related_markers = head(all_markers[all_markers$cluster==3,]$gene,10)
# 
# flt3_related_df = data.frame(row.names = flt3_related_markers)
# flt3_related_df$`Gene Set` = "FGFR Enriched"
# 
# skin_epithelial_markers = c(
#   #"KRT5",#"KRT14","PEX","GLI2",
#   "KRT17","KRT6A","KRT6B","KRT6C",
#   "EGFR",#"BIRC5","TP63",
#   "KRT16","GSDMA","CTSV","ERBB3","PTPRF",
#   #"IL36RN","FLG","IL36B","SOX9","FOXE1","VDR","IL36G",
#   #"RUNX2","CARM1","TGFB1","ANGIO","AAMP","ANGPT1","ANGPT2","CALR","CXCL9","CXCL10",
#   #"EPO","FGF1","FGF2","HOXB4","PGF","SERPINB5","TIMP1",
#   "ERBB2",#"IL37",
#   "CD82",#"FOXC1",
#   "CST6",
#   "SOX9",
#   "DPP4",
#   "SMAD3",
#   "IGFL1","IGFL2"
#   #,"MUC1"
#   #"SOX9"
#   #,"VIL1"
#   #"KRT9","KRT16","KRT17"
# )
# 
# skin_epithelial_markers = head(all_markers[all_markers$cluster==4,]$gene,10)
# 
# skin_epithelial_df = data.frame(row.names = skin_epithelial_markers)
# skin_epithelial_df$`Gene Set` = "17q Amplified"
# 
# other_markers = c(
#   "SELP","IL6","IL8",#"IL6R","IL10",#"CTNNB1","CD163",#"TNFRSF1A","DLL4",
#   "STAT3",#"SOCS3",
#   "BCL3",#"IL21",
#   "PTPRD","OLR1",
#   "ST6GALNAC5","CXCL1",
#   "CXCL2","CXCL3",
#   #"EPCAM","CD163",
#   "HIF1A",
#   #"KRT8","KRT18","ICAM1",
#   "THBS1","PTPRO","EGR1",
#   "VEGFA"
#   #,"KDR"
#   
#   #"MMP2","MMP9","CTNNB1","CDH1","ITGA","S100A4"
#   #,"VEGFB"#,"VEGFC"
#   #"CD302"
#   #,"KRT19","KRT7"
# )
# 
# # other_markers = all_markers[grepl("KRT",all_markers$gene),]$gene
# # other_markers = other_markers[!duplicated(other_markers)]
# 
# other_markers = head(all_markers[all_markers$cluster==5,]$gene,10)
# 
# other_markers_df = data.frame(row.names = other_markers)
# other_markers_df$`Gene Set` = "Radioresistant"
# 
# #gene_set_df = rbind(angiogenesis_df,stem_cell_df,lymphoangiogensis_df,epitheliod_df,diff_epitheliod_df,other_markers_df)
# gene_set_df = rbind(
#   angiogenesis_df,
#   lymphoangiogensis_df,
#   flt3_related_df,
#   skin_epithelial_df,
#   other_markers_df
# )
# 
# endMT_markers = c(
#   "XAF1","LRP1",
#   "SOX9",
#   "PYCR1",
#   "TWIST","SMAD3","ZEB2","SNAI1","SNAI2",
#   "ACTA2","VWF",
#   "TGFB1I1","TGFBRAP1","TGFA"
# )
# 
# caf_markers = c(
#   "ACTA2","FAP",
#   "THY1", # General fibroblast marker
#   "DLL4",
#   "IGFBP3","IGF2"
# )
# 
# LEC_markers = c(
#   "ACKR4"
# )
# 
# gene_set_df$`Gene Set` = factor(
#   gene_set_df$`Gene Set`,
#   levels = c("Angiogenesis",
#              "Lymphangiogenesis",
#              "FGFR Enriched",
#              #"FLT3 Enriched",
#              "17q Amplified",
#              "Radioresistant"
#              ))
# gene_set_df
# 
# # cd_markers = paste0("CD",1:400)
# # il_markers = paste0("IL",1:400)
# 
# asc_new_markers = rownames(gene_set_df)
# Idents(so) = so@meta.data$seurat_clusters_renamed
# asc_new_marker_heatmap = DoHeatmap(so, features = asc_new_markers, 
#                                         disp.max = 3.5,disp.min = -3.5,slot="vst_scaled",label=TRUE)
# 
# ## prettify the results for figure
# # First sort by cluster ID to make sure they are in the same order
# asc_cluster_idents = so@meta.data %>% arrange(seurat_clusters_renamed,`entity:sample_id`) %>% dplyr::select(seurat_clusters_renamed)  %>% mutate(seurat_clusters_renamed = paste0("Cluster ",seurat_clusters_renamed))
# asc_primary_sites =  so@meta.data %>% arrange(seurat_clusters_renamed,`entity:sample_id`) %>% dplyr::select(`Primary Site (Recombined)`) 
# asc_new_marker_heatmap_data = GetAssayData(so,slot="vst_scaled")[asc_new_markers,rownames(asc_cluster_idents)]
# 
# gene_set_df
# all_markers_dedupped = all_markers %>% arrange(-avg_log2FC) %>% distinct(gene, .keep_all = TRUE)
# all_new_markers_FC = all_markers_dedupped %>% filter(gene %in% rownames(gene_set_df))
# all_new_markers_FC = all_new_markers_FC %>% mutate(
#   significance = case_when(
#     p_val_adj < 0.01 ~ "adj. p < 0.01",
#     p_val_adj < 0.1 ~ "adj. p < 0.1",
#     p_val < 0.05 ~ "nom. p < 0.05",
#     TRUE ~ "nom. p >= 0.05"
#   )
# )
# 
# all_new_markers_FC
# sig_high_fc = all_markers %>% filter(p_val_adj < 0.05)
# 
# cluster_color_map = seurart_cluster_palette
# names(cluster_color_map) = paste0("Cluster ",names(seurart_cluster_palette))
# 
# ## Load heatmap data
# asc_new_marker_heatmap_matrix = asc_new_marker_heatmap_data
# 
# ## Make sure the palette is consistent to the palette used for clustering
# gene_set_palette = c(
#   "Angiogenesis"=pal_npg("nrc")(5)[1],
#   "Lymphangiogenesis"=pal_npg("nrc")(5)[2],
#   "FGFR Enriched"=pal_npg("nrc")(5)[3],
#   "17q Amplified"=pal_npg("nrc")(5)[4],
#   "Radioresistant"=pal_npg("nrc")(5)[5]
# )
# # Annotate the overall category of genes each highlighted genes represent
# asc_gene_set_annotations = rowAnnotation(
#   ` `=gene_set_df$`Gene Set`,
#   annotation_legend_param = list(` ` = list(title = "Gene Set")),
#   col=list(` `=gene_set_palette)
# )
# 
# FC_col_annot = colorRamp2(c(0, 5), c("white", "red"))
# #logp_col_annot = colorRamp2(c(0, 5), c("white", "purple"))
# purples_cols = brewer.pal(4,"Purples")
# significance_pal = c(
#   "adj. p < 0.01" = purples_cols[4],
#   "adj. p < 0.1" = purples_cols[3],
#   "nom. p < 0.05" = purples_cols[2],
#   "nom. p >= 0.05" = purples_cols[1]
# )
# # Annotate the gene's FC in its highest expressed cluster relative to other clusters
# asc_FC_annotations = rowAnnotation(
#   `Gene log2FC`=all_new_markers_FC[rownames(gene_set_df),"avg_log2FC"],
#   `Gene Significance`=all_new_markers_FC[rownames(gene_set_df),"significance"],
#   #`Gene -log10(p-value)`=-log10(all_new_markers_FC[rownames(gene_set_df),"p_val_adj"]),
#   annotation_name_rot=90,
#   col=list(`Gene log2FC`=FC_col_annot,`Gene Significance`=significance_pal)
# )
# 
# # Annotate each sample by their most representative somatic, germline, or both mutations
# # This annotation require inputs from running all other scripts so if it is your first time 
# # Running this block you should come back after running all other scripts first
# rep_mut_df = read.csv("03_Tumor_WES_Analysis/outputs/tables/representative_mutation_table.csv",check.names = FALSE)
# rep_mut_df = rep_mut_df[,c("sample_alias_cleaned","tumor_mutation_id","tumor_mutation_id_short","Hugo_Symbol")] %>%
#   dplyr::rename(
#     rep_tumor_mutation_id = "tumor_mutation_id",
#     rep_tumor_mutation_id_short = "tumor_mutation_id_short",
#     rep_Hugo_Symbol = "Hugo_Symbol"
#   )
# representative_somatic_muts_merged = merge(so@meta.data,rep_mut_df,by="sample_alias_cleaned",all.x=TRUE)
# 
# ## Also add in sample TMB
# tmb_df = read.csv("03_Tumor_WES_Analysis/outputs/tables/total_mutations_per_sample.csv")
# tmb_df_smol = tmb_df %>%
#   dplyr::select(sample_alias,TMB_all_mutations)
# representative_somatic_muts_merged_with_tmb = merge(
#   representative_somatic_muts_merged,tmb_df_smol,by.x="sample_alias_cleaned",
#   by.y="sample_alias",all.x=TRUE
# )
# # assign tmb categories
# representative_somatic_muts_merged_with_tmb = representative_somatic_muts_merged_with_tmb %>%
#   mutate(
#     tmb_cat = case_when(
#       TMB_all_mutations <= 1 ~ "<= 1 mut/MB",
#       TMB_all_mutations <= 3 ~ "<= 3 mut/MB",
#       TMB_all_mutations <= 10 ~ "<= 10 mut/MB",
#       TMB_all_mutations > 10 ~ "> 10 mut/MB",
#       TRUE ~ "NA"
#     )
#   )
# # TMB palette
# orange_cols = brewer.pal(4,"Oranges")
# tmb_palette = c(
#   "<= 1 mut/MB" = orange_cols[1],
#   "<= 3 mut/MB" = orange_cols[2],
#   "<= 10 mut/MB" = orange_cols[3],
#   "> 10 mut/MB" = orange_cols[4],
#   "NA" = "gray"
# )
# 
# germline_rep_mut_df = read.csv("05_Germline_WES_Analysis/outputs/tables/representative_germline_pv_table.csv",check.names=FALSE)
# germline_cols_to_keep = c("Sample","inferred_ancestry_PCA","inferred_sex","Hugo_Symbol","STUDY ID","germline_gene_priority_score")
# germline_rep_mut_df_smol = germline_rep_mut_df[,germline_cols_to_keep]
# germline_rep_mut_df_smol = germline_rep_mut_df_smol %>%
#   dplyr::rename(
#     "germline_inferred_sex" = "inferred_sex",
#     "germline_inferred_ancestry_PCA" = "inferred_ancestry_PCA",
#     "germline_Hugo_Symbol" = "Hugo_Symbol"
#   )
# # Making sure only one sample per study id
# length(unique(germline_rep_mut_df_smol$`STUDY ID`)) == length(unique(germline_rep_mut_df_smol$`Sample`))
# representative_somatic_muts_merged_with_germline = merge(
#   representative_somatic_muts_merged_with_tmb,
#   germline_rep_mut_df_smol,
#   by="STUDY ID",
#   all.x=TRUE
# )
# 
# representative_somatic_muts_merged = representative_somatic_muts_merged_with_germline %>%
#   arrange(seurat_clusters_renamed,`entity:sample_id`) %>% 
#   replace_na(list(rep_Hugo_Symbol="Not available"))
# 
# representative_germline_vars_merged = representative_somatic_muts_merged_with_germline %>%
#   arrange(seurat_clusters_renamed,`entity:sample_id`) %>% 
#   #select(`germline_Hugo_Symbol`) %>% 
#   replace_na(list(germline_Hugo_Symbol="No PV Detected"))
# 
# ## Calculate the enrichment of KDR mutations 
# tumor_dna_avail_samples = representative_somatic_muts_merged[representative_somatic_muts_merged$highlight_som_rep != "Not available",]
# both_tumor_dna_rna_avail_samples = so@meta.data[so@meta.data$`entity:sample_id` %in% tumor_dna_avail_samples$`entity:sample_id`,]
# kdr_mutated_samples = tumor_dna_avail_samples[tumor_dna_avail_samples$highlight_som_rep=="KDR",]
# both_tumor_dna_rna_avail_samples$has_kdr_mutation = both_tumor_dna_rna_avail_samples$`entity:sample_id` %in% kdr_mutated_samples$`entity:sample_id`
# kdr_count_table = table(both_tumor_dna_rna_avail_samples$has_kdr_mutation,both_tumor_dna_rna_avail_samples$seurat_clusters_renamed)
# kdr_clust1_enrichment = fisher.test(matrix(c(8, 1, 10, 30),nrow = 2))
# kdr_clust1_enrichment
# 
# ## Only highlight genes that are interesting (somatic)
# rep_som_mut_palette = c(
#   "POT1"=pal_npg("nrc")(10)[1],
#   "TP53"=pal_npg("nrc")(10)[2],
#   "CFTR"=pal_npg("nrc")(10)[3],
#   "KDR"=pal_npg("nrc")(10)[4],
#   "PLCG1"=pal_npg("nrc")(10)[5],
#   #"FLG"=pal_npg("nrc")(10)[3],
#   #"FLT1"=pal_npg("nrc")(10)[1],
#   #"FLT3"=pal_npg("nrc")(10)[3],
#   #"FLT4"=pal_npg("nrc")(10)[4],
#   #"FLG"=pal_npg("nrc")(10)[8],
#   #"BRAF"=pal_npg("nrc")(10)[9],
#   "Others"="black",
#   "Not available"="gray"
# )
# 
# rep_germ_var_palette = c(
#   "POT1"=pal_npg("nrc")(10)[1],
#   "TP53"=pal_npg("nrc")(10)[2],
#   "CFTR"=pal_npg("nrc")(10)[3],
#   "BRCA2"=pal_npg("nrc")(10)[6],
#   "CHEK2"=pal_npg("nrc")(10)[7],
#   #"FLG"=pal_npg("nrc")(10)[3],
#   #"BRCA1"=pal_npg("nrc")(10)[7],
#   #"MUTYH"=pal_npg("nrc")(10)[10],
#   #"PKHD1"="red",
#   #"USH2A"="pink",
#   #"PAH"=pal_npg("nrc")(10)[11],
#   #"GJB2"=pal_npg("nrc")(10)[12],
#   "Others"="black",
#   "No PV Detected"="lightblue"
# )
# 
# rep_som_genes_to_highlight = names(rep_som_mut_palette)
# rep_germ_genes_to_highlight = names(rep_germ_var_palette)
# 
# representative_somatic_muts_merged = representative_somatic_muts_merged %>%
#   mutate(highlight_som_rep=case_when(
#     rep_Hugo_Symbol %in% rep_som_genes_to_highlight ~ rep_Hugo_Symbol,
#     TRUE ~ "Others"
#   )) %>%
#   mutate(`SEX (EHR_EXTRACTED)`=case_when(
#     `SEX (EHR_EXTRACTED)`=="F" ~ "Female",
#     `SEX (EHR_EXTRACTED)`=="M" ~ "Male",
#     TRUE ~ NA
#   ))
# representative_germline_vars_merged = representative_germline_vars_merged %>%
#   mutate(highlight_germ_rep=case_when(
#     germline_Hugo_Symbol %in% rep_germ_genes_to_highlight ~ germline_Hugo_Symbol,
#     TRUE ~ "Others"
#   ))
# representative_germline_vars_merged
# 
# bx_palette = c(
#   "YES" = "red",
#   "NO" = "lightpink"
# )
# representative_somatic_muts_merged = representative_somatic_muts_merged %>%
#   mutate(is_epithelioid=ifelse(as.character(BX_EPITHELIOID)=="YES (FOCAL)","YES",as.character(BX_EPITHELIOID)))
# representative_somatic_muts_merged = representative_somatic_muts_merged %>%
#   mutate(is_spindle=ifelse(as.character(BX_SPINDLE_CELL)=="YES (FOCAL)","YES",as.character(BX_SPINDLE_CELL)))
# 
# nuclear_grade_palette = c(
#   "HIGH" = "orange",
#   "LOW" = "lightpink"
# )
# 
# representative_somatic_muts_merged = representative_somatic_muts_merged %>%
#   mutate(has_mets = case_when(
#     `HAS METS AT DX (EHR_EXTRACTED)`==1 ~ "YES",
#     `HAS METS AT DX (EHR_EXTRACTED)`==0 ~ "NO",
#     TRUE ~ NA
#   ))
# 
# mets_dx_palette = c(
#   "YES" = "darkblue",
#   "NO" = "lightblue"
# )
# 
# representative_somatic_muts_merged = representative_somatic_muts_merged %>% mutate(
#   cutaneous_viz = case_when(
#     `CUTANEOUS AS (EHR_EXTRACTED)` == 0 ~ "Non-cutaneous AS",
#     `CUTANEOUS AS (EHR_EXTRACTED)` == 1 ~ "Cutaneous AS",
#     TRUE ~ NA
#   )
# )
# 
# ## Top Annotation
# age_col_annot = colorRamp2(c(20, 80), c("white", "purple"))
# asc_top_annotations = HeatmapAnnotation(
#   `Primary Site`= asc_primary_sites$`Primary Site (Recombined)`,
#   `Cutaneous` = representative_somatic_muts_merged$`cutaneous_viz`,
#   "RAAS/LAAS" = representative_somatic_muts_merged$RAAS_LAAS_Class,
#   `Age`= representative_somatic_muts_merged$`Age (Combined)`,
#   "Sex" = representative_somatic_muts_merged$`SEX (EHR_EXTRACTED)`,
#   "Epithelioid" = representative_somatic_muts_merged$is_epithelioid,
#   "Spindle Cell" = representative_somatic_muts_merged$is_spindle,
#   "Nuclear Grade" = representative_somatic_muts_merged$BX_NUCLEAR_GRADE,
#   #"Mets at Dx" = representative_somatic_muts_merged$has_mets,
#   #"Vasoformative" = representative_somatic_muts_merged$BX_VASOFORMATIVE,
#   col=list(
#     `Primary Site`=primary_site_palette,
#     `Cutaneous` = cutaneous_palette,
#     "RAAS/LAAS" = RAAS_class_palette,
#     "Sex" = sex_clin_palette,
#     "Age" = age_col_annot,
#     "Epithelioid" = bx_palette,
#     "Spindle Cell"=bx_palette,
#     "Nuclear Grade" = nuclear_grade_palette
#     #"Mets at Dx" = mets_dx_palette
#     )
# )
# 
# ## Bottom Annotation
# asc_rep_mut_annnotation = HeatmapAnnotation(
#   `Repr. Somatic Mut.` = representative_somatic_muts_merged$highlight_som_rep,
#   `Repr. Germline Var.` = representative_germline_vars_merged$highlight_germ_rep,
#   `TMB` = representative_germline_vars_merged$tmb_cat,
#   col = list(
#     `Repr. Somatic Mut.`=rep_som_mut_palette,
#     `Repr. Germline Var.`=rep_germ_var_palette,
#     `TMB` = tmb_palette
#     )
# )
# 
# ## Create a column dendrogram based on highlighted genes
# so_small_embeddings= Embeddings(so,reduction = "pca")[,1:2]
# #so_small_embedding = FetchData(so,vars=rownames(gene_set_df),layer = "vst_scaled")
# so_small_hclust = hclust(dist(so_small_embeddings))
# so_small_dendro = as.dendrogram(so_small_hclust)
# 
# # Make the actual heatmap plot
# complex_heatmap_fig = ComplexHeatmap::Heatmap(
#   asc_new_marker_heatmap_matrix,
#   cluster_rows=FALSE,
#   cluster_columns=TRUE,
#   #cluster_column_slices = T,
#   #cluster_columns = so_small_dendro,
#   #show_row_names=FALSE,
#   show_column_names=FALSE,
#   top_annotation=asc_top_annotations,
#   bottom_annotation = asc_rep_mut_annnotation,
#   right_annotation=asc_gene_set_annotations,
#   left_annotation=asc_FC_annotations,
#   column_split = asc_cluster_idents,
#   cluster_column_slices = FALSE,
#   column_dend_reorder = F,
#   #column_split = FALSE,
#   row_split = gene_set_df,
#   heatmap_legend_param = list(
#     title = "Scaled Expression", at = c(-3, 0, 3) 
#     #labels = c("neg_two", "zero", "pos_two")
#   )
# )
# #pdf("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_cluster_gene_set_heatmap_dendro.pdf",width=24,height=18)
# pdf("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_cluster_gene_set_heatmap_dendro_cluster_by_site.pdf",width=24,height=18)
# draw(complex_heatmap_fig, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",merge_legend = TRUE)
# dev.off()
# 
# 
# ## How many of the total parenchymal breast samples available were classified into cluster 1?
# parenchymal_breast_df = so@meta.data %>% filter(`Primary Site (Recombined)`=="Breast (Parenchymal)")
# parenchymal_breast_in_clust1_df = total_parenchymal_breast %>% filter(seurat_clusters_renamed==1)
# total_parenchymal_breast = dim(parenchymal_breast_df)[1]
# total_parenchymal_breast_in_clust1 = dim(parenchymal_breast_in_clust0_df)[1]
# print(paste0(total_parenchymal_breast_in_clust1,"/",total_parenchymal_breast,"=",total_parenchymal_breast_in_clust0/total_parenchymal_breast))
# 
# ## How many of the total cutaneous breast samplse available were classified into either cluster 1 or 2?
# cutaneous_breast_df = so@meta.data %>% filter(`Primary Site (Recombined)`=="Breast (Cutaneous)")
# cutaneous_breast_df_in_clust2_df = cutaneous_breast_df %>% filter(seurat_clusters_renamed==2)
# cutaneous_breast_df_in_clust3_df = cutaneous_breast_df %>% filter(seurat_clusters_renamed==3)
# total_cutaneous_breast = dim(cutaneous_breast_df)[1]
# total_cutaneous_breast_in_clust2 = dim(cutaneous_breast_df_in_clust2_df)[1]
# total_cutaneous_breast_in_clust3 = dim(cutaneous_breast_df_in_clust3_df)[1]
# total_cutaneous_breast_in_clust23 = total_cutaneous_breast_in_clust2 + total_cutaneous_breast_in_clust3
# total_cutaneous_breast_in_clust2/total_cutaneous_breast
# total_cutaneous_breast_in_clust3/total_cutaneous_breast
# total_cutaneous_breast_in_clust23/total_cutaneous_breast
# 
# ## How many of the KDR-containing tumors fall into cluster 0?
# tumors_with_kdr_mut = representative_somatic_muts_merged %>% filter(highlight_som_rep=="KDR")
# 
# ## How many of the TP53 mutation-containing tumors fall into cluster 0?
# #tumors_with_tp53_mut = representative_somatic_muts_merged %>% filter(highlight_som_rep=="TP53")
# 
# 
# ## Check whether clusters with enriched genomic locations show WES segment enrichment as well
# ## Export metadata for downstream analysis
# write.csv(so@meta.data,"02_Tumor_RNA_Analysis/outputs/post_analysis_metadata.csv")
# 
# 
# 
# ## Find the gene sets each gene belongs to
# msigdb_c2 = msigdbr(species = "Homo sapiens", category = "C2")
# msigdb_go = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
# 
# find_common_gene_sets = function(gene1,gene2,gs_db) {
#   gs1 = gs_db %>% filter(gene_symbol==gene1)
#   gs2 = gs_db %>% filter(gene_symbol==gene2)
#   overlap_geneset = intersect(gs1$gs_name,gs2$gs_name)
#   gs_subset = gs_db %>% filter(gs_name %in% overlap_geneset)
#   gs_subset_sizes = gs_subset %>% group_by(gs_name) %>% summarise(n())
#   gs_subset_merged = merge(gs_subset,gs_subset_sizes,on="gs_name") %>% distinct(gs_name,.keep_all = TRUE) %>% arrange(`n()`)
#   gs_subset_merged_small = gs_subset_merged[,c("gs_name","n()")]
#   return(gs_subset_merged_small)
# }
# find_common_gene_sets("CD93","CD34",msigdb_go)
# 
# (msigdb_go %>% filter(gs_name=="GOBP_POSITIVE_REGULATION_OF_MESENCHYMAL_CELL_PROLIFERATION"))$gene_symbol
# 
# unique((msigdb_go %>% filter(grepl("STEM_CELL",gs_name)))$gs_name)
# test_pathsets = (msigdb_go %>% filter(grepl("STEM_CELL",gs_name)))
# test_pathsets %>% group_by(gs_name) %>% summarise(n()) %>% arrange(`n()`)# %>% distinct(gene_symbol,.keep_all = FALSE)
# 
# 
# (msigdb_go %>% filter(gs_name=="GOBP_POSITIVE_REGULATION_OF_MESENCHYMAL_CELL_PROLIFERATION"))$gene_symbol
# 
# ggsave(asc_new_marker_heatmap,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_top_cluster_marker_heatmap.png",dpi=300,width=24,height=12)
# 
# ## Prettify for figure
# 
# 
# 
# ## Perform ORA on genes that are positively enriched in each cluster (GO BP and KEGG Pathways)
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(msigdbr)
# 
# get_cluster_enrichment = function(so,all_markers,cluster) {
#   # Subset to cluster marker genes and take nominally significant genes
#   sub_df = all_markers[all_markers$cluster==cluster,]
#   #sig_genes = rownames(sub_df[sub_df$p_val < 0.05,])
#   sig_genes = rownames(sub_df[sub_df$p_val_adj < 0.1,])
#   bg_genes = rownames(so@assays$RNA$counts)
#   # GO Enrichment
#   print("GO Enrichment in progress...")
#   enrich_obj = enrichGO(
#     gene = sig_genes,
#     universe = bg_genes,
#     OrgDb = org.Hs.eg.db,
#     keyType = "SYMBOL",
#     ont="BP",
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.2,
#     minGSSize = 10,
#     maxGSSize = 500
#   )
#   # Kegg Enrichent
#   print("KEGG Enrichment in progress...")
#   sig_genes_uniprot = bitr(sig_genes,fromType="SYMBOL",toType="UNIPROT",OrgDb = org.Hs.eg.db)$UNIPROT 
#   sig_genes_uniprot = sig_genes_uniprot[!duplicated(sig_genes_uniprot)]
#   bg_genes_uniprot = bitr(bg_genes,fromType="SYMBOL",toType="UNIPROT",OrgDb = org.Hs.eg.db)$UNIPROT
#   bg_genes_uniprot = bg_genes_uniprot[!duplicated(bg_genes_uniprot)]
#   kegg_gson = gson_KEGG("hsa",keyType="uniprot")
#   #options(clusterProfiler.download.method = "wget")
#   kegg_obj = enricher(
#     gene=sig_genes_uniprot,
#     universe=bg_genes_uniprot,
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.2,
#     gson = kegg_gson,
#     minGSSize = 10,
#     maxGSSize = 500
#   )
#   # MSigDB oncogenic signature enrichment (upregulated gene sets only)
#   print("MSIGDB c6 Enrichment in progress...")
#   msigdb_c6_gene_sets = msigdbr(species = "Homo sapiens", category = "C6") %>% 
#     dplyr::filter(grepl("up-regulated",gs_description)) %>%
#     dplyr::select(gs_name, gene_symbol)
#   msig_c6_obj = enricher(
#     gene=sig_genes,
#     universe=bg_genes,
#     pvalueCutoff = 0.05,
#     pAdjustMethod = "BH",
#     qvalueCutoff = 0.2,
#     TERM2GENE=msigdb_c6_gene_sets,
#     minGSSize = 10,
#     maxGSSize = 500
#   )
#   
#   enrich_res = list(
#     "go"=enrich_obj,
#     "kegg"=kegg_obj,
#     "msigdb_c6"=msig_c6_obj
#   )
#   return(enrich_res)
# }
# 
# cluster_0_enrichement = get_cluster_enrichment(so,all_markers,0)
# cluster_1_enrichement = get_cluster_enrichment(so,all_markers,1)
# cluster_2_enrichement = get_cluster_enrichment(so,all_markers,2)
# cluster_3_enrichement = get_cluster_enrichment(so,all_markers,3)
# cluster_4_enrichement = get_cluster_enrichment(so,all_markers,4)
# 
# draw_enriched_bars = function(enrichment_res,enrich_type,title) {
#   plot = barplot(
#     enrichment_res[[enrich_type]], 
#     drop = TRUE, 
#     showCategory = 15, 
#     title = title,
#     font.size = 8
#   )
#   return(plot)
# }
# 
# ## Draw the GO term enrichments and Pathway enrichments side by side
# c0_go_bar = draw_enriched_bars(cluster_0_enrichement,"go","Parenchymal Breast\ncluster 0: GO-BP")
# c1_go_bar = draw_enriched_bars(cluster_1_enrichement,"go","Cutaneous Breast\ncluster 1: GO-BP")
# c2_go_bar = draw_enriched_bars(cluster_2_enrichement,"go","HNFS\ncluster 2: GO-BP")
# c3_go_bar = draw_enriched_bars(cluster_3_enrichement,"go","Mixed Types\ncluster 3: GO-BP")
# c4_go_bar = draw_enriched_bars(cluster_4_enrichement,"go","Immune Infiltrated\ncluster 4: GO-BP")
# 
# c0_kegg_bar = draw_enriched_bars(cluster_0_enrichement,"kegg", "Parenchymal Breast\ncluster 0: KEGG")
# c1_kegg_bar = draw_enriched_bars(cluster_1_enrichement,"kegg", "Cutaneous Breast\ncluster 1: KEGG")
# c2_kegg_bar = draw_enriched_bars(cluster_2_enrichement,"kegg", "HNFS\ncluster 2: KEGG")
# c3_kegg_bar = draw_enriched_bars(cluster_3_enrichement,"kegg", "Mixed Types\ncluster 3: KEGG")
# c4_kegg_bar = draw_enriched_bars(cluster_4_enrichement,"kegg", "Immune Infiltrated\ncluster 4: KEGG")
# 
# c0_msigdb_bar = draw_enriched_bars(cluster_0_enrichement,"msigdb_c6", "Parenchymal Breast\ncluster 0: Oncogenic Signatures")
# c1_msigdb_bar = draw_enriched_bars(cluster_1_enrichement,"msigdb_c6", "Cutaneous Breast\ncluster 1: Oncogenic Signatures")
# c2_msigdb_bar = draw_enriched_bars(cluster_2_enrichement,"msigdb_c6", "HNFS\ncluster 2: Oncogenic Signatures")
# c3_msigdb_bar = draw_enriched_bars(cluster_3_enrichement,"msigdb_c6", "Mixed Types\ncluster 3: Oncogenic Signatures")
# c4_msigdb_bar = draw_enriched_bars(cluster_4_enrichement,"msigdb_c6", "Immune Infiltrated\ncluster 4: Oncogenic Signatures")
# 
# go_kegg_msigdb_bars_combined = cowplot::plot_grid(
#   c0_go_bar,c1_go_bar,c2_go_bar,c3_go_bar,c4_go_bar,
#   c0_kegg_bar,c1_kegg_bar,c2_kegg_bar,c3_kegg_bar,c4_kegg_bar,
#   c0_msigdb_bar,c1_msigdb_bar,c2_msigdb_bar,c4_msigdb_bar,c3_msigdb_bar,
#   labels = "AUTO",
#   ncol = 5
# )
# go_kegg_msigdb_bars_combined
# ggsave(go_kegg_msigdb_bars_combined,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_cluster_top_markers_go_kegg_msigdb_enrichments.png",dpi=300,width=24,height=16)
# 
# 
# ## Check endothelial migration genes
# top_markers_dimplot = FeaturePlot(so, features = c(
#   "FLT1","APLNR","EBF3","NRP1",
#   "MAP4K2","SMYD2","VAV3","NRP2"
# ),pt.size = 2,ncol = 4,reduction = "pca")
# top_markers_dimplot
# 
# # Helper function to format expression subsets
# get_gene_subset_long = function(so,genes) {
#   gene_subset = FetchData(object=so,vars=genes,layer="vst_scaled")
#   gene_subset$og_id = rownames(gene_subset)
#   gene_subset_long = gene_subset %>% pivot_longer(cols = -c("og_id"),names_to = "gene_name",values_to = "count_value")
#   gene_subset_merged = merge(gene_subset_long,so@meta.data,on="og_id",all.x=TRUE)
#   gene_subset_merged$gene_name = factor(gene_subset_merged$gene_name,levels=genes)
#   return(gene_subset_merged)
# }
# 
# immuno_gene_subset_merged = get_gene_subset_long(so,genes=c("CTLA4","CD274","PDCD1"))
# immuno_labeller = labeller(gene_name=c("CTLA4"="CTLA4","CD274"="CD274 (PD-L1)","PDCD1"="PDCD1 (PD-1)"))
# 
# ## Make a boxplot of CTLA4 expressions
# immuno_boxplot = ggplot(immuno_gene_subset_merged, 
#                         aes(x = `seurat_clusters_renamed_str`, y = count_value)) +
#   geom_violin(width = 1.2) +
#   geom_boxplot(width = 0.1, outlier.shape = NA) +
#   facet_grid(cols = vars(gene_name), labeller = immuno_labeller) +
#   stat_compare_means(comparisons = list(c("Cluster 1", "Cluster 2"),
#                                         c("Cluster 1", "Cluster 3"),
#                                         c("Cluster 2", "Cluster 3"),
#                                         c("Cluster 2", "Cluster 4")
#                                         ),
#                      na.rm = TRUE, label = "p.format") +
#   stat_compare_means(method = "anova", label.y = 5) +
#   geom_quasirandom(aes(color = `Primary Site (Recombined)`)) +
#   theme_minimal() +
#   theme(strip.text.y = element_text(angle = 0)) + # Adjust facet label orientation
#   labs(x = "Cluster", y = "Normalized Expression Level", color = "Cluster") +
#   guides(color=guide_legend(title="Primary Site"))
# ggsave(immuno_boxplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_immuno_expr_boxplot.png",dpi=300,width=12,height=8)
# 
# ## FLT3, MYCL, and KRT
# clust3_4_gene_subsets_merged = get_gene_subset_long(so,genes=c("FLT3","MYCL","FGFR2","FGFR3","KRT1","IL18","IL37"))
# #clust3_4_gene_labeller = labeller(gene_name=c("CTLA4"="CTLA4","CD274"="CD274 (PD-L1)","PDCD1"="PDCD1 (PD-1)"))
# clust_gene_boxplot = ggplot(clust3_4_gene_subsets_merged, 
#                         aes(x = `seurat_clusters_renamed_str`, y = count_value)) +
#   geom_violin(width = 1.2) +
#   geom_boxplot(width = 0.1, outlier.shape = NA) +
#   facet_grid(cols = vars(gene_name)) +
#   stat_compare_means(method = "anova", label.y = 5) +
#   geom_quasirandom(aes(color = `Primary Site (Recombined)`)) +
#   theme_minimal() +
#   theme(strip.text.y = element_text(angle = 0)) + # Adjust facet label orientation
#   labs(x = "Cluster", y = "Normalized Expression Level", color = "Cluster") +
#   guides(color=guide_legend(title="Primary Site"))
# ggsave(clust_gene_boxplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_clust_3_4_expr_boxplot.png",dpi=300,width=24,height=6)
# 
# ## Best marker gene for each cluster
# clust_markers_long = get_gene_subset_long(so,genes=c("FLT1","FLT4","FGFR2","IGFL2","ST6GALNAC5"))
# clust_marker_boxplot = ggplot(clust_markers_long, 
#                             aes(x = `seurat_clusters_renamed_str`, y = count_value)) +
#   geom_violin(width = 1.2) +
#   geom_boxplot(width = 0.1, outlier.shape = NA) +
#   facet_grid(cols = vars(gene_name)) +
#   stat_compare_means(method = "anova", label.y = 5) +
#   geom_quasirandom(aes(color = `Primary Site (Recombined)`)) +
#   theme_minimal() +
#   theme(strip.text.y = element_text(angle = 0)) + # Adjust facet label orientation
#   labs(x = "Cluster", y = "Normalized Expression Level", color = "Cluster") +
#   guides(color=guide_legend(title="Primary Site"))
# ggsave(clust_marker_boxplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_all_clust_best_marker_boxplot.png",dpi=300,width=16,height=6)
# 
# # 
# # 
# # # Keratin expression by cutaneous status
# # # "KRT1","FLT3","MYCL","MYC",
# # # EGFR, ERBB3
# # kertain_gene_subset_merged = get_gene_subset_long(so,genes=c("FGFR1","FGFR2","FGFR3","FGFR4","FGFRL1","FGFR6"))
# # kertain_gene_subset_merged = kertain_gene_subset_merged %>% 
# #   mutate(cutaneous_viz = case_when(
# #     `CUTANEOUS AS (EHR_EXTRACTED)` == "0" ~ "Non-cutaneous AS",
# #     `CUTANEOUS AS (EHR_EXTRACTED)` == "1" ~ "Cutaneous AS",
# #     TRUE ~ "NA"
# #   ))
# # #kertain_gene_subset_merged$BX_EPITHELIOID
# # kertain_gene_subset_merged = kertain_gene_subset_merged[!is.na(kertain_gene_subset_merged$BX_EPITHELIOID),]
# # keratin_boxplot = ggplot(kertain_gene_subset_merged, 
# #                         aes(x = `Primary Site (Recombined)`, y = count_value)) +
# #   geom_violin(width = 1.2) +
# #   geom_boxplot(width = 0.1, outlier.shape = NA) +
# #   facet_grid(rows = vars(gene_name)) +
# #   stat_compare_means(comparisons = list(c("Breast (Parenchymal)", "HNFS"),
# #                                         c("Breast (Parenchymal)", "Breast (Cutaneous)")
# #                                         #c("Cluster 1", "Cluster 3"),
# #                                         #c("Cluster 2", "Cluster 3"),
# #                                         #c("Cluster 2", "Cluster 4")
# #   ),
# #   na.rm = TRUE, label = "p.format") +
# #   stat_compare_means(method = "anova", label.y = 5) +
# #   geom_quasirandom(aes(color = `cutaneous_viz`)) +
# #   theme_minimal() +
# #   theme(strip.text.y = element_text(angle = 0)) + # Adjust facet label orientation
# #   labs(x = "Cluster", y = "Normalized Expression Level", color = "Cluster") +
# #   guides(color=guide_legend(title="Cutaneous Status"))
# # ggsave(keratin_boxplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_keratin_expr_boxplot.png",dpi=300,width=12,height=8)
# # 
# # 
# #  ## Make biplot between marker genes
# plot_biplot = function(so,gene1,gene2,color_attr) {
#   gene_subset = so@assays$RNA$vst_scaled[c(gene1,gene2),]
#   gene_subset = as.data.frame(t(gene_subset),check.names=FALSE)
#   gene_subset$og_id = rownames(gene_subset)
#   gene_subset = merge(gene_subset,so@meta.data,on="og_id",how="left")
#   regression_result =  lm(as.formula(paste0(gene2, " ~ ", gene1)), data = gene_subset)
#   coef = coef(regression_result)[2]
#   p_value = summary(regression_result)$coefficients[2, 4]
#   genes_biplot = ggplot(gene_subset,aes_string(x=gene1,y=gene2)) +
#     geom_point(size=4,aes_string(color=color_attr)) +
#     geom_smooth(method = "lm", se = FALSE, color = "blue") +
#     #geom_point(size=4,aes(color=gene_subset[[color_attr]])) +
#     theme_minimal() +
#     geom_vline(xintercept = 0) +
#     geom_hline(yintercept = 0) +
#     annotate("text", x = Inf, y = Inf, label = paste0("Coef: ", signif(coef,3)), hjust = 1.1, vjust = 2, size = 4) +
#     annotate("text", x = Inf, y = Inf, label = paste0("P-val: ", signif(p_value,3)), hjust = 1.1, vjust = 3.5, size = 4)
#     #guides(fill=guide_legend(title="Cluster"))
#     #geom_smooth(method='lm')
#   return(genes_biplot)
# }
# 
# so_hnfs
# test_plot = FeatureScatter(object = so_hnfs, feature1 = "ERVMER61-1", feature2 = "METTL4", pt.size = 5)
# ggsave(test_plot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_erv_ifn_biplot.png",dpi=300,width=24,height=6)



# # Select a gene of interest
# gene_of_interest <- "ERVMER61-1"
# # Extract the expression matrix
# expr_matrix <- GetAssayData(so_hnfs, slot = "data") # Using log-normalized data
# expr_matrix <- as.matrix(expr_matrix)  # Ensure it’s a matrix
# expr_matrix_t <- t(expr_matrix) 
# # Compute correlation of the gene of interest with all other genes
# cor_values <- cor(expr_matrix_t, method = "pearson")[, gene_of_interest]
# 
# # Sort by correlation and find the most correlated gene (excluding itself)
# sorted_cor <- sort(cor_values, decreasing = TRUE)
# most_correlated_gene <- names(sorted_cor)[2]  # Skipping the first one, which is itself
# 
# # Print result
# cat("The most correlated gene with", gene_of_interest, "is", most_correlated_gene, "with correlation", sorted_cor[2])


#plot_biplot(so,"ERVMER61-1","IFNG","`Primary Site (Recombined)`")

#so@meta.data$seurat_clusters_renamed_str = paste0("Cluster ",so@meta.data$seurat_clusters_renamed)
# 
# myc_ctla4_biplot = plot_biplot(so,"MYC","CTLA4","`Primary Site (Recombined)`") + guides(color=guide_legend(title="Primary Site"))
# myc_pdl1_biplot = plot_biplot(so,"MYC","CD274","`Primary Site (Recombined)`") + guides(color=guide_legend(title="Primary Site")) + labs(y="CD274 (PD-L1)")
# myc_pd1_biplot = plot_biplot(so,"MYC","PDCD1","`Primary Site (Recombined)`") + guides(color=guide_legend(title="Primary Site")) + labs(y="PDCD1 (PD-1)")
# combined_myc_plot =  plot_grid(
#   myc_ctla4_biplot,myc_pdl1_biplot,myc_pd1_biplot,align = "hv",ncol = 3,labels = c('A', 'B',"C")#bp5,bp6,bp7,bp8,align = "hv",ncol = 2
# )
# ggsave(combined_myc_plot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_myc_immuno_biplot.png",dpi=300,width=24,height=6)
# 
# #bp1 = plot_biplot(so,"FLT4","MYC","seurat_clusters_renamed_str") + guides(color=guide_legend(title="Cluster"))
# bp1 = plot_biplot(so,"MYC","FLT4","`Primary Site (Recombined)`") + guides(color=guide_legend(title="Primary Site"))
# bp2 = plot_biplot(so,"MYC","FLT4","seurat_clusters_renamed_str")+ guides(color=guide_legend(title="Cluster"))
# bp3 = plot_biplot(so,"MYC","MYCL","`Primary Site (Recombined)`")+ guides(color=guide_legend(title="Primary Site"))
# bp4 = plot_biplot(so,"MYC","MYCL","`seurat_clusters_renamed_str`")+ guides(color=guide_legend(title="Cluster"))
# #bp5 = plot_biplot(so,"IL18","IL37","`Primary Site (Recombined)`")+ guides(color=guide_legend(title="Primary Site"))
# #bp6 = plot_biplot(so,"IL18","IL37","`seurat_clusters_renamed_str`")+ guides(color=guide_legend(title="Cluster"))
# #bp5 = plot_biplot(so,"IL18","IL37","seurat_clusters_renamed_str")+ guides(color=guide_legend(title="Cluster"))
# #bp5 = plot_biplot(so,"FLT1","FLT4","seurat_clusters_renamed_str")+ guides(color=guide_legend(title="Cluster"))
# #bp6 = plot_biplot(so,"FLT1","FLT3","seurat_clusters_renamed_str")+ guides(color=guide_legend(title="Cluster"))
# #bp7 = plot_biplot(so,"KRT1","KRT17","seurat_clusters_renamed_str")+ guides(color=guide_legend(title="Cluster"))
# #bp8 = plot_biplot(so,"CXCL1","FLT1","seurat_clusters_renamed_str")+ guides(color=guide_legend(title="Cluster"))
# 
# combined_biplots = plot_grid(
#   bp1,bp2,bp3,bp4,align = "hv",ncol = 2,labels = c('A', 'B',"C","D")#bp5,bp6,bp7,bp8,align = "hv",ncol = 2
# )
# ggsave(combined_biplots,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_myc_mycl_combined_biplots.png",dpi=300,width=12,height=10)
# 
# ## Make a biplot of all the vegf ligands 
# ## Along with their regression coefficient and p-values
# vegf_receptors = c("FLT1","KDR","FLT3","FLT4","NRP1","NRP2","PDGFRA","PDGFRB")
# vegf_ligands = c("VEGFA","VEGFB","VEGFC","FIGF","PIGF","PDGFA","PDGFB","PDGFC","PDGFD")
# vegf_subset_expr = data.frame(t(so@assays$RNA$vst_scaled[c(vegf_receptors,vegf_ligands),]))
# vegf_subset_expr$og_id = rownames(vegf_subset_expr)
# vegf_subset_expr_merged = merge(vegf_subset_expr,so@meta.data,on="og_id",how="left")
# vegf_subset_expr_merged
# 
# combinations = expand.grid(
#   receptor = vegf_receptors,
#   ligand = vegf_ligands,
#   stringsAsFactors = FALSE
# )
# 
# scatter_data = combinations %>%
#   rowwise() %>%
#   mutate(
#     receptor_value = list(vegf_subset_expr_merged[[receptor]]),
#     ligand_value = list(vegf_subset_expr_merged[[ligand]]),
#     receptor_label = receptor,
#     ligand_label = ligand,
#     cluster_value = list(vegf_subset_expr_merged[["seurat_clusters_renamed_str"]])
#   ) %>%
#   unnest(cols = c(receptor_value, ligand_value,cluster_value))
# 
# regression_results <- scatter_data %>%
#   group_by(receptor_label, ligand_label) %>%
#   summarise(
#     coef = coef(lm(ligand_value ~ receptor_value))[2],
#     p_value = summary(lm(ligand_value ~ receptor_value))$coefficients[2, 4],
#     highlight = ifelse(summary(lm(ligand_value ~ receptor_value))$coefficients[2, 4] < 0.05/72, TRUE, FALSE),
#     .groups = "drop"
#   )
# 
# scatter_data <- scatter_data %>%
#   left_join(regression_results, by = c("receptor_label", "ligand_label"))
# 
# vegf_biplots = ggplot(scatter_data, aes(x = receptor_value, y = ligand_value, color=cluster_value)) +
#   geom_point(alpha = 0.6) +
#   geom_smooth(method = "lm", se = FALSE, color = "blue") +
#   facet_grid(rows = vars(receptor_label), cols = vars(ligand_label), switch = "both") +
#   labs(
#     x = "Ligand Expression",
#     y = "Receptor Expression"
#   ) +
#   theme_minimal() +
#   geom_vline(xintercept = 0) + 
#   geom_hline(yintercept = 0) +
#   theme(
#     strip.text = element_text(size = 8),
#     axis.text = element_text(size = 6),
#     axis.title = element_text(size = 10),
#     panel.spacing = unit(1, "lines")
#   ) +
#   geom_text(
#     data = regression_results,
#     aes(
#       x = Inf, y = Inf,
#       label = paste0("Coef: ", signif(coef, 3), "\nP: ", signif(p_value, 3))
#     ),
#     inherit.aes = FALSE,
#     hjust = 1.1, vjust = 1.1,
#     size = 3
#   ) +
#   geom_rect(
#     data = regression_results %>% filter(highlight),  # Highlighted panels
#     aes(
#       xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
#     ),
#     inherit.aes = FALSE,
#     color = "red",
#     fill = NA,
#     size = 0.8
#   )
# 
# ggsave(vegf_biplots,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_VEGF_biplots.png",dpi=300,width=24,height=16)
# 
# 
# 
# 
# plot_biplot(so,"FLT1","NRP1","seurat_clusters")
# plot_biplot(so,"KRT1","KRT2","seurat_clusters")
# plot_biplot(so,"FLT4","CTSV","seurat_clusters")
# 
# 
# plot_biplot(so,"MYC","MYCL","PRIMARY SITE (Combined)")#"seurat_clusters")
# plot_biplot(so,"NRP1","NRP2","seurat_clusters")
# plot_biplot(so,"FLT1","FLT4","seurat_clusters")
# plot_biplot(so,"FLT1","FLG","PRIMARY SITE (Combined)")
# plot_biplot(so,"FLT3","MYCL","seurat_clusters")
# plot_biplot(so,"IL18","MYC","seurat_clusters")
# plot_biplot(so,"CXCL1","CXCL12","seurat_clusters")
# 
# 
# ## NRP1-NRP2
# nrp1_nrp2_biplot_by_site = plot_biplot(so,"NRP1","NRP2","PRIMARY SITE (Combined)")
# nrp1_nrp2_biplot_by_site
# 
# ## NRP1-MYC axis separates major clusters
# nrp_myc_biplot_by_site = plot_biplot(so,"NRP1","MYC","PRIMARY SITE (Combined)")
# ggsave(nrp_myc_biplot_by_site,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_NRP1_MYC_biplot.png",dpi=300,width=24,height=16)
# 
# ## TERT-FLT1 gets parenchymal breast samples
# tert_flt1_biplot_by_site = plot_biplot(so,"TERT","FLT1","PRIMARY SITE (Combined)")
# ggsave(tert_flt1_biplot_by_site,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_TERT_FLT1_biplot.png",dpi=300,width=24,height=16)
# 
# ## MYC-CTLA4 gets cutaneous breast samples
# myc_ctla4_biplot = plot_biplot(so,"CTLA4","MYC","PRIMARY SITE (Combined)")
# ggsave(myc_ctla4_biplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_MYC_CTLA4_biplot.png",dpi=300,width=24,height=16)
# 
# ## FGFR3?
# plot_biplot(so,"FGFR3","APOBEC3G","PRIMARY SITE (Combined)")
# plot_biplot(so,"FGFR3","PRKCA","PRIMARY SITE (Combined)")
# 
# 
# plot_biplot(so,"SEMA4D","MAPK1","PRIMARY SITE (Combined)")
# plot_biplot(so,"MAPK1","TBX1","PRIMARY SITE (Combined)")
# flt1_foxp3_biplot = plot_biplot(so,"FLT1","FOXP3","PRIMARY SITE (Combined)")
# ctla4_foxp3_biplot = plot_biplot(so,"CTLA4","FOXP3","PRIMARY SITE (Combined)")
# fgfr3_foxp3_biplot = plot_biplot(so,"FGFR3","FOXP3","PRIMARY SITE (Combined)")
# 
# ggsave(flt1_foxp3_biplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_FLT1_FOXP3_biplot.png",dpi=300,width=24,height=16)
# ggsave(ctla4_foxp3_biplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_CTLA4_FOXP3_biplot.png",dpi=300,width=24,height=16)
# ggsave(fgfr3_foxp3_biplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_FGFR3_FOXP3_biplot.png",dpi=300,width=24,height=16)
# 
# 
# VlnPlot(so,features = c(
#   "VEGFA","VEGFB","VEGFC","VEGFD",
#   "PDGFA","PDGFB","PDGFC",
#   "FLT1","KDR","FLT3","FLT4"
#   #"TIE1","MYC",
#   #"DAXX","ATRX","ARID1A","POT1","TERT"
#   ),slot="vst_scaled"
# )
# 
# so
# 
# curated_heatmap = DoHeatmap(so,features=c(
#   "FLT1","NRP1","FGFR1","APLNR","PDGFA","PDGFB","TERT","POT1","ATRX","DAXX",
#   "KDR","NRP2","TBX1",
#   "TIE1","FLT4","MYC",
#   "RAC1","CTLA4","BCR","TGFB1",#"IL2RA","ATR","CHEK1", 
#   "CHEK2", "MAD2L1", "BUB1B", "TP53", "BRCA1", "BRCA2", "RAD51", "NBN","MDM2","SMYD3","FAT1",
#   "KRAS","HRAS","NRAS",
#   "ERBB2","ERBB3","EGFR","MET","ALK","RET","PDGFRA","KIT","HLA-A","HLA-B","HLA-C",
#   "FGFR2","FGFR3","FGFR4","LRP5","ANGPT2","TEK",
#   "MAPK3","MAPK1",
#   #"ERK1","ERK2","MAPK","LRP6","ANG2","TIE2",
#   "FLT3","PDGFC",
#   "PIK3CA","PTPRB","PTPRD","PLCG1",
#   "VEGFA","VEGFB","VEGFC",
#   "PDCD1","FOXP3","CXCL10","CCR5",
#   "CD19","MS4A1","IL10","LAG3","CCL2","IL17A","TNF"
# ),slot = "vst_scaled",disp.min=-2.5,disp.max=2.5)
# curated_heatmap
# ggsave(curated_heatmap,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_curated_genes_heatmap.png",dpi=300,width=24,height=16)
# 
# 
# ## MYC-CTLA4 axis
# plot_biplot(so,"KRAS","HRAS","PRIMARY SITE (Combined)")
# plot_biplot(so,"MYC","MAPK1","PRIMARY SITE (Combined)")
# 
# plot_biplot(so,"KRAS","PERP","PRIMARY SITE (Combined)")
# 
# plot_biplot(so,"DLL4","KDR","PRIMARY SITE (Combined)")
# plot_biplot(so,"DLL4","FLT1","PRIMARY SITE (Combined)")
# plot_biplot(so,"DLL4","NOTCH1","PRIMARY SITE (Combined)")
# 
# plot_biplot(so,"VEGFA","FLT1","PRIMARY SITE (Combined)")
# 
# 
# ## MYC-TERT axis
# myc_tert_biplot = plot_biplot(so,"TERT","KDR","PRIMARY SITE (Combined)")
# myc_tert_biplot
# 
# ## MYC-CTLA4-KRAS axis
# myc_kras_biplot_by_site = plot_biplot(so,"MYC","KRAS","PRIMARY SITE (Combined)")
# ctla4_myc_biplot_by_site = plot_biplot(so,"MYC","CTLA4","PRIMARY SITE (Combined)")
# plot_grid(myc_kras_biplot_by_site,ctla4_myc_biplot_by_site)
# 
# tp53_smyd2_biplot_by_site = plot_biplot(so,"UGT2A2","PNPLA5","PRIMARY SITE (Combined)")
# tp53_smyd2_biplot_by_site
# 
# ## CTLA4-CLEC2A (cutaneous vs, rest)
# tp53_smyd2_biplot_by_site = plot_biplot(so,"CTLA4","CD274","PRIMARY SITE (Combined)")
# tp53_smyd2_biplot_by_site 
# 
# plot_biplot(so,"CDHR1","CHP2","RAAS_LAAS_Class")
# plot_biplot(so,"NRP1","NRP2","RAAS_LAAS_Class")
# 
# ## Confirm relation between VEGF receptors (FLT/NRP)
# nrp_biplot_by_site = plot_biplot(so,"NRP1","NRP2","PRIMARY SITE (Combined)")
# nrp_biplot_by_raas = plot_biplot(so,"NRP1","NRP2","RAAS (EHR_EXTRACTED)")
# plot_grid(nrp_biplot_by_site,nrp_biplot_by_raas)
# 
# flt_biplot_by_site = plot_biplot(so,"TERT","MYC","PRIMARY SITE (Combined)")
# flt_biplot_by_site
# 
# nrp_subset = so@assays$RNA$vst_scaled[c("NRP1","NRP2","VEGFA","KDR"),]
# nrp_subset = as.data.frame(t(nrp_subset),check.names=FALSE)
# nrp_subset$og_id = rownames(nrp_subset)
# nrp_subset = merge(nrp_subset,so@meta.data,on="og_id",how="left")
# color_attr = "PRIMARY SITE (Combined)"
# 
# nrp_biplot = ggplot(nrp_subset,aes(x=NRP1,y=NRP2)) +
#   geom_point(size=3,aes(color=nrp_subset[[color_attr]])) +
#   theme_minimal() +
#   geom_vline(xintercept = 0) + 
#   geom_hline(yintercept = 0) +
#   geom_smooth(method='lm')
# nrp_biplot
# 
# ## Write the seurat derived metadata
# write.table(so@meta.data,"data/processed/rna/02_seurat_processed_metadata.tsv",sep="\t")
# 
# 
# ## Show mutations overlaid on top of PCA
# mutation_dimplot = DimPlot(so, reduction = "pca", group.by=c(
#   "POT1_Mut_Germline","POT1_Mut_Onehot","TERT_Mut_Onehot","ATRX_Mut_Onehot",
#   "FLT1_Mut_Onehot","KDR_Mut_Onehot","FLT3_Mut_Onehot","FLT4_Mut_Onehot",
#   "POLE_Mut_Onehot","TP53_Mut_Onehot","PIK3CA_Mut_Onehot",
#   "PLCG1_Mut_Onehot","PTPRB_Mut_Onehot",
#   "KRAS_Mut_Onehot","NRAS_Mut_Onehot","HRAS_Mut_Onehot"
#   ),pt.size=3)
# mutation_dimplot
# ggsave(mutation_dimplot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_mutation_status_on_pcs.png",dpi=300,width=24,height=12)
# 
# 
# ## Visualize the top enriched genes in each cluster
# top_markers_plot = VlnPlot(so, features = c(
#   "FLT1","APLNR","EBF3","NRP1",
#   "MAP4K2","SMYD2","VAV3","NRP2",
#   "KRTAP24-1","CCDC172","UGT2A2","PNPLA5",
#   "ST6GALNAC5","C7","NCR3LG1","PTGS2",
#   "KLRF2","LCE5A","CLEC2A","SP8"
#   #"CARD18","SERPINB12","PGLYRP3","SEC14L1","IQCJ-SCHIP1","TSPAN18",
#   #"ZNF462","BMPR1A","GLI3","TRIB3","CDC25B","DUSP5"
#   ),pt.size = 2,ncol = 4)
# top_markers_plot
# ggsave(top_markers_plot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_top_markers_violin_plot.png",dpi=300,width=24,height=12)
# 
# ## The same content but on PCs
# top_markers_dimplot = FeaturePlot(so, features = c(
#   "FLT1","APLNR","EBF3","NRP1",
#   "MAP4K2","SMYD2","VAV3","NRP2",
#   "KRTAP24-1","CCDC172","UGT2A2","PNPLA5",
#   "ST6GALNAC5","C7","NCR3LG1","PTGS2",
#   "KLRF2","LCE5A","CLEC2A","SP8"
#   #"CARD18","SERPINB12","PGLYRP3","SEC14L1","IQCJ-SCHIP1","TSPAN18",
#   #"ZNF462","BMPR1A","GLI3","TRIB3","CDC25B","DUSP5"
# ),pt.size = 2,ncol = 4,reduction = "pca")
# top_markers_dimplot
# ggsave(top_markers_dimplot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_top_markers_pca_plot.png",dpi=300,width=24,height=12)
# 
# ## Show Continuous Feature overlaid on top of PCA
# continuous_clinical_dimplot = FeaturePlot(so, reduction = "pca", features=c("Age (Combined)","tmb"),pt.size=3)
# continuous_clinical_dimplot
# ggsave(continuous_clinical_dimplot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_continuous_clinical_feature_pca_plot.png",dpi=300,width=24,height=12)
# 
# ## Show gene expression overlaid on top of PCA
# genes_to_extract = c("MYC","TP53","PLCG1","PTPRB","POT1","HRAS","NRAS","KRAS","TERT","PIK3CA")#,"PIK3CA","TERT")  # Replace with your genes of interest
# expr_dim_plot = FeaturePlot(so,reduction = "pca", features=c(
#   "FLT1","KDR","FLT4","FLT3",
#   "PGF","PIGF","YAP1","MIB2",
#   "MYC","TP53","PLCG1","PTPRB","POT1","RAD51",
#   "HRAS","NRAS","KRAS","TERT","PIK3CA","GLI3","FAM83G","NFKB1"
# ),slot="vst_scaled",pt.size=2)#,cols = c("blue","red"))
# expr_dim_plot
# ggsave(expr_dim_plot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_gene_expression_pca_plot.png",dpi=300,width=24,height=12)
# 
# ## Combine all previous information and combine them into heatmap
# top_cluster_markers = all_markers %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 0.3) %>%
#   #dplyr::filter(!grepl("KRT",gene)) %>%
#   slice_head(n = 20) %>%
#   ungroup()
# 
# top_cluster_markers_heatmap = DoHeatmap(so, features = top_cluster_markers$gene,slot = "vst.scaled",disp.min=-3,disp.max=3)
# top_cluster_markers_heatmap
# ggsave(top_cluster_markers_heatmap,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_top_cluster_marker_heatmap.png",dpi=300,width=24,height=12)
# 
# ## Save the seurat derived clusters as metadata
# #write.csv("seurat")
# write.csv(so@meta.data,"reference_data/RNA_Seq/outputs/seurat/seurat_processed_metadata.csv")





