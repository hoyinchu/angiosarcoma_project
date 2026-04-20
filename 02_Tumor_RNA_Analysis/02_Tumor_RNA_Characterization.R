library(Seurat)
library(cowplot)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsci)
library(dendextend)


## Load count data and reformat IDs for downstream compatibility
count_data = read.csv("../data/processed/rna/00_filtered_gene_counts.csv",row.names = 1,check.names = FALSE)
vst_data = read.csv("../data/processed/rna/00_filtered_gene_counts.vst.csv",row.names = 1,check.names = FALSE)
tpm_data = read.csv("../data/processed/rna/00_filtered_gene_tpm.csv",row.names = 1,check.names = FALSE)
  
colnames(count_data) = gsub("-", ".", colnames(count_data))
colnames(count_data) = gsub("_", ".", colnames(count_data))
colnames(vst_data) = gsub("-", ".", colnames(vst_data))
colnames(vst_data) = gsub("_", ".", colnames(vst_data))
colnames(tpm_data) = gsub("-", ".", colnames(tpm_data))
colnames(tpm_data) = gsub("_", ".", colnames(tpm_data))

# Store a scaled version of the data row-wise (gene) for visualization
vst_data_scaled = as.data.frame(t(scale(t(vst_data))),check.names=FALSE)
vst_data_unscaled = as.data.frame((vst_data),check.names=FALSE)

#vst_data_scaled_pretty_names = gsub("\\.","_",colnames(vst_data_scaled))
#colnames(vst_data_scaled) =  vst_data_scaled_pretty_names

## Apply same formatting to metadata and reorder it
meta_data = read.csv("../data/processed/rna/00_filtered_sample_metadata.csv",check.names = FALSE)
meta_data$seurat_id = meta_data$`entity:sample_id`
meta_data$seurat_id = gsub("_", ".", meta_data$seurat_id)
meta_data$seurat_id = gsub("-", ".", meta_data$seurat_id)
meta_data[meta_data == ""] = NA
rownames(meta_data) = meta_data$seurat_id

# Create a seurat object and add meta data
#so = CreateSeuratObject(counts = vst_data,min.cells = 0, min.features = 0, meta.data = meta_data)
#so = CreateSeuratObject(counts = tpm_data,min.cells = 0, min.features = 0, meta.data = meta_data)
so = CreateSeuratObject(counts = count_data,min.cells = 0, min.features = 0, meta.data = meta_data)


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
hclust_dist_matrix = dist(t(genes_to_use_for_hclust), method = "euclidean")
hclust_res = hclust(hclust_dist_matrix, method = "ward.D2")

## Asign clusters via cut height
sample_clusters = cutree(hclust_res, k = NULL, h = 200)

# Get sample labels and primary site information
sample_labels = colnames(genes_to_use_for_hclust)
primary_sites = so@meta.data[sample_labels, ]$seurat_clusters_by_site_str#`Primary Site (Recombined)`
unique_clusters = unique(sample_clusters)

# Convert hclust object to dendrogram
dend = as.dendrogram(hclust_res)
# Set branch/group label size to a legible value (e.g., 0.8)
dend = color_branches(dend, k = length(unique_clusters), groupLabels = as.character(unique_clusters), labels_cex = 0.1)

# Color the leaves by primary site
dend_order = order.dendrogram(dend)
ordered_labels = primary_sites[dend_order]
leaf_colors = cluster_by_primary_site_palette[ordered_labels]
names(leaf_colors) = sample_labels

# --- UPDATES TO SHOW SQUARES INSTEAD OF TEXT ---
# 1. Hide the leaf text labels (replace with empty strings)
dend = set(dend, "labels", rep("", length(ordered_labels))) 
# 2. Set the leaf shape to a filled square (pch=22)
dend = set(dend, "leaves_pch", 22) 
# 3. Set the leaf fill color to the primary site color
dend = set(dend, "leaves_bg", leaf_colors)
# 4. Set the size of the square markers 
dend = set(dend, "leaves_cex", 0.5)
# 5. Set the marker border color
dend = set(dend, "leaves_col", "black") # Using black for border contrast
dend = set(dend, "hang_leaves", -2)
dend
# Plot colored dendrogram
#pdf("./outputs/plots/02_seurat_plots/02_gene_expression_cluster_dendrogram.pdf",width=6,height=3)
pdf("./outputs/plots/02_seurat_plots/02_gene_expression_cluster_dendrogram_with_legend.pdf",width=6,height=3)

par(mar=c(0, 3, 0, 0.5))
plot(dend,ylab = "Height")
# Add legend for primary sites
legend("topright",
       legend = unique(primary_sites),
       fill = cluster_by_primary_site_palette,
       pch = 22, # Ensure the legend key is a square
       title = "Primary Site",
       cex = 0.8)
dev.off()

# Add cluster assignments to metadata (rest of your code is unchanged)
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
so@meta.data$seurat_clusters_renamed = as.integer(so@meta.data$seurat_clusters)
so@meta.data$seurat_clusters_renamed = as.integer(so@meta.data$seurat_clusters_by_site)
so@meta.data$seurat_clusters_renamed_str = paste0("Cluster ",so@meta.data$seurat_clusters_renamed)

## Save Seurat Object
if (FALSE) {
  SaveSeuratRds(so,file = "../data/processed/rna/ASCSeuratObj2025.rds")
}

## Load Seurat Object for basic sumary statistics
so = readRDS("../data/processed/rna/ASCSeuratObj2025.rds")

# ## Plot estimate purity by site
estimate_purity_boxplot = ggplot(
  data=so@meta.data,aes(x=seurat_clusters_by_site_str,y=ESTIMATE_purity,color=seurat_clusters_by_site_str)
) + geom_boxplot() + pretty_plot() + L_border() + scale_color_npg() +
  geom_jitter(size=0.5)

estimate_purity_boxplot_clean = estimate_purity_boxplot + 
  theme(legend.position = "none",axis.title = element_blank(),
        axis.text.x = element_text(size=5), axis.ticks.x = element_line())
estimate_purity_boxplot_clean

cowplot::ggsave2("./outputs/plots/07_specificity_plots/ESTIMATE_purity_barplot_by_sites_with_legend.pdf",estimate_purity_boxplot,dpi=300,width=3.3,height=1.5)
cowplot::ggsave2("./outputs/plots/07_specificity_plots/ESTIMATE_purity_barplot_by_sites.pdf",estimate_purity_boxplot_clean,dpi=300,width=3.3,height=1.5)

## Plot cluster membership
clust_mem_df = table(so@meta.data$hcluster_by_expr_str,so@meta.data$`Primary Site (Recombined)`) %>% as.data.frame()
clust_mem_df_plot = ggplot(clust_mem_df,aes(x=Var1,y=Freq,fill=Var2)) + geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values=primary_site_palette) + pretty_plot() + L_border() +
  scale_y_continuous(expand=c(0,Inf)) +
  labs(x="Cluster",y="# of Samples")
clust_mem_df_plot_clean = clust_mem_df_plot + theme(legend.position = "none",axis.title = element_blank(),
        axis.text.x = element_text(size=5), axis.ticks.x = element_line())

clust_mem_df_plot_clean
cowplot::ggsave2("./outputs/plots/02_seurat_plots/cluster_membership_clean.pdf",clust_mem_df_plot_clean,dpi=300,width=1.85,height=1.5)
cowplot::ggsave2("./outputs/plots/02_seurat_plots/cluster_membership.pdf",clust_mem_df_plot,dpi=300,width=1.8,height=1.5)

## Do cluster membership but flip it
clust_mem_df_plot_site_x = ggplot(clust_mem_df,aes(x=Var2,y=Freq,fill=Var1)) + 
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values=seurart_cluster_by_expression_palette) + pretty_plot() + L_border() +
  scale_y_continuous(expand=c(0,Inf)) + 
  theme(axis.title = element_blank(),axis.text.x = element_text(size=5,angle = 45,vjust=1,hjust=1),
        axis.ticks.x = element_line()) +
  theme(legend.position = "right")
  #theme(legend.position = "none",axis.text.x = element_blank())
clust_mem_df_plot_site_x

clust_mem_df_plot_site_x
cowplot::ggsave2("./outputs/plots/02_seurat_plots/clust_mem_df_plot_site_x.pdf",clust_mem_df_plot_site_x,dpi=300,width=2.8,height=1.35)
cowplot::ggsave2("./outputs/plots/02_seurat_plots/clust_mem_df_plot_site_x_with_legend.pdf",clust_mem_df_plot_site_x,dpi=300,width=2.4,height=1.5)


## Number of parencymal breast samples by site
site_by_cluster_table = table(so@meta.data$seurat_clusters_by_site_str,so@meta.data$hcluster_by_expr_str)
site_totals = rowSums(site_by_cluster_table)
cluster_totals = colSums(site_by_cluster_table)

## Save the Seurat Object
if (FALSE) {
  saveRDS(so, "../data/processed/rna/ASCSeuratObj.rds")
}






