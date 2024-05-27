library(Seurat)
library(cowplot)
library(ggplot2)
library(tidyr)
library(dplyr)
library(fgsea)
library(msigdbr)

setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts")

## Load count data and reformat IDs for downstream compatibility
count_data = read.csv("data/processed/rna/00_filtered_gene_counts.csv",row.names = 1,check.names = FALSE)
vst_data = read.csv("data/processed/rna/00_filtered_gene_counts.vst.csv",row.names = 1,check.names = FALSE)
colnames(vst_data) = gsub("-", ".", colnames(vst_data))
colnames(vst_data) = gsub("_", ".", colnames(vst_data))

# Store a scaled version of the data row-wise (gene) for visualization
vst_data_scaled = as.data.frame(t(scale(t(vst_data))),check.names=FALSE)

## Apply same formatting to metadata and reorder it
meta_data = read.csv("data/processed/rna/00_filtered_sample_metadata.csv",check.names = FALSE)
meta_data[["seurat_id"]] = meta_data$`entity:sample_id`
meta_data[["seurat_id"]] = gsub("_", ".", meta_data[["seurat_id"]])
meta_data[["seurat_id"]] = gsub("-", ".", meta_data[["seurat_id"]])
meta_data[meta_data == ""] = NA
rownames(meta_data) = meta_data$seurat_id

# Create a seurat object and add meta data
so = CreateSeuratObject(counts = vst_data,min.cells = 0, min.features = 0, meta.data = meta_data)
so@meta.data[["og_id"]] = rownames(so@meta.data)
so@meta.data = merge(so@meta.data,meta_data,by.x="og_id",by.y="seurat_id",all.x=TRUE)
rownames(so@meta.data) = so@meta.data[["og_id"]]

so[["RNA"]]$vst_scaled = vst_data_scaled
so = SetIdent(so, value = rownames(so@meta.data))

# Adapted based on rocanja's recommendation in handling bulk data in seurat
# reference: https://github.com/satijalab/seurat/issues/7496

so = SetAssayData(object=so, layer="data", new.data = vst_data)
so = FindVariableFeatures(so,selection.method = "vst", nfeatures = 2000)
## Dimension reduction
so = ScaleData(so)
so = RunPCA(so, features = VariableFeatures(object = so))
## Run unsuperivsed clustering
so = FindNeighbors(so, dims = 1:10)
so = FindClusters(so, resolution = 1.5)
so = RunUMAP(so, dims = 1:10)

# Get the total variance:
get_var_explained = function(so) {
  mat = GetAssayData(so, assay = "RNA", slot = "scale.data")
  pca = so[["pca"]]
  total_variance = sum(matrixStats::rowVars(mat))
  eigen_values = (pca@stdev)^2  
  var_explained = eigen_values / total_variance
  return(var_explained)
}

# Plot PCA
pca_plot_by_batch = DimPlot(so, reduction = "pca",group.by="LC_BATCH",pt.size=4)
pca_plot_by_site = DimPlot(so, reduction = "pca",group.by="PRIMARY SITE (Combined)",pt.size=4)

pc_var_explained = get_var_explained(so)
pca_plot_by_batch = pca_plot_by_batch + 
  xlab(paste0("PC1 (",signif(pc_var_explained[1]*100,4),"%)")) +
  ylab(paste0("PC2 (",signif(pc_var_explained[2]*100,4),"%)"))
pca_plot_by_site = pca_plot_by_site + 
  xlab(paste0("PC1 (",signif(pc_var_explained[1]*100,4),"%)")) +
  ylab(paste0("PC2 (",signif(pc_var_explained[2]*100,4),"%)"))
pca_combined_plot = plot_grid(pca_plot_by_batch, pca_plot_by_site, labels = "AUTO")
pca_combined_plot
ggsave(pca_combined_plot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_pca_plot_by_batch_and_site.png",dpi=300,width = 16,height=8)

## Plot PCA Elbow
elbow_plot = ElbowPlot(so,reduction="pca")
elbow_plot
ggsave(elbow_plot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_pca_elbow_plot.png",dpi=300)

## Check the mean-variance relationships
top10_variable_genes = head(VariableFeatures(so),10)
mean_var_plot = VariableFeaturePlot(so)
mean_var_plot = LabelPoints(plot=mean_var_plot,points=top10_variable_genes,repel=TRUE)
mean_var_plot
ggsave(mean_var_plot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_mean_var_plot.png",dpi=300)

## Visualize the top loadings
gene_loading_plot = VizDimLoadings(so, dims = 1:6, reduction = "pca",nfeatures=15)
gene_loading_plot
ggsave(gene_loading_plot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_gene_loadings_plot.png",dpi=300,width=8,height=12)

## Visualize the genes on the two extreme ends of each PCs
dim_score_plot = DimHeatmap(so, dims = 1:15,balanced = TRUE,combine = TRUE,fast=FALSE)
dim_score_plot
ggsave(dim_score_plot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_dim_heatmap_plot.png",dpi=300,width=8,height=12)

dim_score_first_3_pcs_plot = DimHeatmap(so, dims = 1:3,balanced = TRUE,combine = TRUE,fast=FALSE,nfeatures=50,reduction="pca")
dim_score_first_3_pcs_plot
ggsave(dim_score_first_3_pcs_plot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_dim_heatmap_first_3_pcs_plot.png",dpi=300,width=8,height=12)

## Show clinical metadata overlaid on top of PCA
clinical_dimplot = DimPlot(so,reduction="pca",group.by = c(
  "seurat_clusters","PRIMARY SITE (Combined)","CUTANEOUS AS (EHR_EXTRACTED)","SEX (EHR_EXTRACTED)",
  "BX_SPINDLE_CELL","BX_NUCLEAR_GRADE","BX_VASOFORMATIVE","BX_EPITHELIOID",
  "RAAS_LAAS_Class","CUTANEOUS AS (EHR_EXTRACTED)"),pt.size=3)
clinical_dimplot
#ggsave(clinical_dimplot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_clinical_attributes_on_pcs.png",dpi=300,width=24,height=12)
ggsave(clinical_dimplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_clinical_attributes_on_pcs_res15.png",dpi=300,width=24,height=12)

## clinical markers
# General markers for angisoarcoma
general_asc_markers = c(
  "PECAM1", # PECAM1, most sensitive marker for endpthelial cells, almost universial positive
  "CD34", # Frequently positive in angiosarcoma, common endothelial cell marker
  "VWF", # Less used, endothelial differentiation marker
  "ERG", # Nuclear transcription factor sensitive for endothelial lineage
  "FLI1" #transcription factor expressed in most endothelial cells
)
cutaneous_asc_markers = c(
  "LYVE1", # (Lymphatic Vessel Endothelial Hyaluronan Receptor 1): Useful in distinguishing lymphatic from blood vessel origin, particularly relevant in cutaneous angiosarcoma.
  "PDPN" # (Podoplanin): Marks lymphatic endothelial cells, supporting a lymphatic origin which can sometimes be seen in angiosarcomas, especially those of the skin.
)
breast_asc_markers = c(
  "MYC" # typically associated with secondary angiosarcoma of the breast
)
epithelioid_asc_markers = c(
  "KRT5", # any keratin,
  "EMA" # epithelial membrane antigen
)
other_markers = c(
  "TERT","NRP1","NRP2","FLT1","KDR","FLT3","FLT4",
  "PROX1","FLT4", "ANGPT2",
  "TBX1","TIE1",
  "PLCG1","PTPRB","PTPRC","PTPRD",
  "VEGFA","DUSP1","DUSP10","DUSP26","DUSP12","MK2","BCL2L1","IL6","IL8","CPLA2","COX2","TNF","HBEGF","C5","CXCL12" #p38 MAPK
)

asc_clinical_markers = c(
  general_asc_markers,
  cutaneous_asc_markers,
  breast_asc_markers,
  epithelioid_asc_markers,
  other_markers
)


#Idents(so) = so@meta.data$`PRIMARY SITE (EHR_EXTRACTED)`
Idents(so) = so@meta.data$seurat_clusters
asc_clinical_marker_heatmap = DoHeatmap(so, features = asc_clinical_markers, 
          disp.max = 3.5,disp.min = -3.5,slot="vst_scaled",label=TRUE)
asc_clinical_marker_heatmap
ggsave(asc_clinical_marker_heatmap,file="02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_asc_clinical_markers_heatmap.png",dpi=300,width=24,height=12)

## Find markers associated with each cluster
Idents(so) = so@meta.data$seurat_clusters
all_markers = FindAllMarkers(
  so,
  slot="counts",
  test.use = "wilcox",
  only.pos = TRUE
)
## Save the output
write.csv(all_markers,"02_Tumor_RNA_Analysis/outputs/DEGs/02_Seurat_all_markers.csv")

## Load gene set libraries
fgsea_kegg_set = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c5_set = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c6_set = msigdbr(species = "Homo sapiens", category = "C6") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c8_set = msigdbr(species = "Homo sapiens", category = "C8") %>% split(x = .$gene_symbol, f = .$gs_name)

run_fgsea = function(degs,fgsea_sets) {
  deg_genes = degs %>%
    arrange(desc(avg_log2FC)) %>% 
    dplyr::select(gene, avg_log2FC)
  vec = deg_genes$avg_log2FC; names(vec) = deg_genes$gene
  fgseaRes = fgseaMultilevel(fgsea_sets, stats = vec)
  return(fgseaRes)
}

run_fora = function(degs,fgsea_sets) {
  genes = degs[degs$p_val_adj < 0.1,]$gene
  universe = rownames(so@assays$RNA$counts)
  fora_res = fora(fgsea_sets, genes, universe, minSize = 5, maxSize = 500)
  return(fora_res)
}

plot_fgsea = function(fgsea_res,title) {
  fgsea_res_sorted = fgsea_res %>% arrange((pval))
  p1 = ggplot(fgsea_res_sorted, aes(x = ES, y = -log10(padj), label=pathway)) + 
    geom_text_repel(data = fgsea_res_sorted[fgsea_res_sorted$pval < 0.05,]) +
    geom_hline(yintercept = -log10(0.05),linetype = "dashed") +
    geom_vline(xintercept = 0,linetype = "dashed") +
    geom_point() + labs(x = "Enrichment score", y = "-log10 padj", title=title)# +
    #pretty_plot(fontsize = 8) + L_border()
  return(p1)
}


clust0_degs_fora_res = run_fora(all_markers[all_markers$cluster==0,],fgsea_c8_set)
clust1_degs_fora_res = run_fora(all_markers[all_markers$cluster==1,],fgsea_c8_set)
clust2_degs_fora_res = run_fora(all_markers[all_markers$cluster==2,],fgsea_c8_set)
clust3_degs_fora_res = run_fora(all_markers[all_markers$cluster==3,],fgsea_c8_set)
clust4_degs_fora_res = run_fora(all_markers[all_markers$cluster==4,],fgsea_c8_set)

## Plot to see what genes are enriched in cluster 2 vs. 3
clust2_ora = all_markers[all_markers$cluster==2,]
clust3_ora = all_markers[all_markers$cluster==3,]
clust2_3_ora_merged =  merge(clust2_ora,clust3_ora,by="gene",suffixes = c("_cluster_2","_cluster_3"))
clust2_3_ora_merged$neg_log_pval_clust2 = -log(clust2_3_ora_merged$p_val_adj_cluster_2)
clust2_3_ora_merged$neg_log_pval_clust3 = -log(clust2_3_ora_merged$p_val_adj_cluster_3)
subset_to_viz = clust2_3_ora_merged[
  (clust2_3_ora_merged$neg_log_pval_clust3 > 10) | (clust2_3_ora_merged$neg_log_pval_clust2 > 7.5),
  ]
clust2_vs_clust3_ora_plot = ggplot(clust2_3_ora_merged,aes(x=neg_log_pval_clust2,y=neg_log_pval_clust3,label=gene)) + 
  geom_point() +
  geom_abline() +
  geom_text_repel(data=subset_to_viz)
ggsave("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_clust2_vs_clust3_ora_plot.png",clust2_vs_clust3_ora_plot)

# Cluster 0: Vascular endothelial cells
# Cluster 1: Lymphatic endothelial cells
## https://www.proteinatlas.org/humanproteome/single+cell+type/squamous+epithelial+cells
# Cluster 2: General squamous epithelial cells? (NK-cell rich? KLRF2, CLEC2A)
# Cluster 3: More specialized Squamous Epithelial cells (T-cell rich? CD82, CST6 / Secretory?)
# Cluster 4: Mixed Type

clust0_degs_fgsea_res = run_fgsea(all_markers[all_markers$cluster==0,],fgsea_c8_set)
clust1_degs_fgsea_res = run_fgsea(all_markers[all_markers$cluster==1,],fgsea_c8_set)
clust2_degs_fgsea_res = run_fgsea(all_markers[all_markers$cluster==2,],fgsea_c8_set)
clust3_degs_fgsea_res = run_fgsea(all_markers[all_markers$cluster==3,],fgsea_c8_set)
clust4_degs_fgsea_res = run_fgsea(all_markers[all_markers$cluster==4,],fgsea_c8_set)

clust0_fgsea_plot = plot_fgsea(clust0_degs_fgsea_res,"Cell Type GSEA, Cluster 0")
clust1_fgsea_plot = plot_fgsea(clust1_degs_fgsea_res,"Cell Type GSEA, Cluster 1")
clust2_fgsea_plot = plot_fgsea(clust2_degs_fgsea_res,"Cell Type GSEA, Cluster 2")
clust3_fgsea_plot = plot_fgsea(clust3_degs_fgsea_res,"Cell Type GSEA, Cluster 3")
clust4_fgsea_plot = plot_fgsea(clust4_degs_fgsea_res,"Cell Type GSEA, Cluster 4")

combined_celltype_gsea_plot = plot_grid(clust0_fgsea_plot,clust1_fgsea_plot,clust2_fgsea_plot,clust3_fgsea_plot,clust4_fgsea_plot)


## Check the top genes in each cluster
top_cluster_markers = all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.1) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  #dplyr::filter(!grepl("KRT",gene)) %>%
  slice_head(n = 20) %>%
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

## Plot the ASC markers combining clinical and new features
cluster_0_markers = c(
  "CD34","CD36","CD93","TERT",#"VWF","ERG","FLI1",
  "NRP1","FLT1","PDGFB","PIGF",
  "VEGFC","FIGF","ITGA5"
)
cluster_1_markers = c(
  "CD320",
  #"CD4","CD8A","CD8B",
  "NRP2","FLT4","LYVE1","PROX1","PECAM1","PDPN","MYC",
  "TBX1","TIE1"#"ANGPT2",
)
cluster_2_markers = c(
  "CD44",
  "CDH1",#"VCAM1",
  "FLT3",
  #"KRT1","FLG",
  "KRT5","SERPINB5","GATA3", #"CASP14","TP73","CDH3",
  "CLEC2A","KLRF2"
  #,"KRT10"
  #"KRT1","KRT2","KRT5","KRT10","KRT14","KRT15","KRT24","KRT77","KRT78","KRT80"
)
cluster_3_markers = c(
  #"KRT5",#"KRT14",
  "KRT6A","KRT16",
  "CD82",
  "CST6","MUC1","VIL1"
  #"KRT9","KRT16","KRT17"
)
cluster_4_markers = c(
  "EPCAM","CD163",
  "KRT8","KRT18",
  "VEGFA","VEGFB","VEGFC","FGF2"
  #"CD302"
  #,"KRT19","KRT7"
)
#keratin_markers = paste0("KRT",1:99)
#keratin_markers = c(keratin_markers,"KRT6A","KRT6B")

# cluster_4_markers = c(
#   "EGF","TGFB1",
#   #"VEGFA","VEGFB","VEGFC","VEGFD","PTPRD","CD144"#,"TNF","KDR","TP53"
# )

cd_markers = paste0("CD",1:400)
il_markers = paste0("IL",1:400)

asc_new_markers = c(
  cluster_0_markers,
  cluster_1_markers,
  cluster_2_markers,
  cluster_3_markers,
  cluster_4_markers
  #il_markers
  #cd_markers
  #keratin_markers
  #other_markers
)
Idents(so) = so@meta.data$seurat_clusters
asc_new_marker_heatmap = DoHeatmap(so, features = asc_new_markers, 
                                        disp.max = 3.5,disp.min = -3.5,slot="vst_scaled",label=TRUE)
asc_new_marker_heatmap
ggsave(asc_new_marker_heatmap,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_top_cluster_marker_heatmap.png",dpi=300,width=24,height=12)

## Perform ORA on genes that are positively enriched in each cluster (GO BP and KEGG Pathways)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

get_cluster_enrichment = function(so,all_markers,cluster) {
  # Subset to cluster marker genes and take nominally significant genes
  sub_df = all_markers[all_markers$cluster==cluster,]
  #sig_genes = rownames(sub_df[sub_df$p_val < 0.05,])
  sig_genes = rownames(sub_df[sub_df$p_val_adj < 0.1,])
  bg_genes = rownames(so@assays$RNA$counts)
  # GO Enrichment
  print("GO Enrichment in progress...")
  enrich_obj = enrichGO(
    gene = sig_genes,
    universe = bg_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont="BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 500
  )
  # Kegg Enrichent
  print("KEGG Enrichment in progress...")
  sig_genes_uniprot = bitr(sig_genes,fromType="SYMBOL",toType="UNIPROT",OrgDb = org.Hs.eg.db)$UNIPROT 
  sig_genes_uniprot = sig_genes_uniprot[!duplicated(sig_genes_uniprot)]
  bg_genes_uniprot = bitr(bg_genes,fromType="SYMBOL",toType="UNIPROT",OrgDb = org.Hs.eg.db)$UNIPROT
  bg_genes_uniprot = bg_genes_uniprot[!duplicated(bg_genes_uniprot)]
  kegg_gson = gson_KEGG("hsa",keyType="uniprot")
  #options(clusterProfiler.download.method = "wget")
  kegg_obj = enricher(
    gene=sig_genes_uniprot,
    universe=bg_genes_uniprot,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    gson = kegg_gson,
    minGSSize = 10,
    maxGSSize = 500
  )
  # MSigDB oncogenic signature enrichment (upregulated gene sets only)
  print("MSIGDB c6 Enrichment in progress...")
  msigdb_c6_gene_sets = msigdbr(species = "Homo sapiens", category = "C6") %>% 
    dplyr::filter(grepl("up-regulated",gs_description)) %>%
    dplyr::select(gs_name, gene_symbol)
  msig_c6_obj = enricher(
    gene=sig_genes,
    universe=bg_genes,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    TERM2GENE=msigdb_c6_gene_sets,
    minGSSize = 10,
    maxGSSize = 500
  )
  
  enrich_res = list(
    "go"=enrich_obj,
    "kegg"=kegg_obj,
    "msigdb_c6"=msig_c6_obj
  )
  return(enrich_res)
}

cluster_0_enrichement = get_cluster_enrichment(so,all_markers,0)
cluster_1_enrichement = get_cluster_enrichment(so,all_markers,1)
cluster_2_enrichement = get_cluster_enrichment(so,all_markers,2)
cluster_3_enrichement = get_cluster_enrichment(so,all_markers,3)
cluster_4_enrichement = get_cluster_enrichment(so,all_markers,4)

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

## Draw the GO term enrichments and Pathway enrichments side by side
c0_go_bar = draw_enriched_bars(cluster_0_enrichement,"go","Parenchymal Breast\ncluster 0: GO-BP")
c1_go_bar = draw_enriched_bars(cluster_1_enrichement,"go","Cutaneous Breast\ncluster 1: GO-BP")
c2_go_bar = draw_enriched_bars(cluster_2_enrichement,"go","HNFS\ncluster 2: GO-BP")
c3_go_bar = draw_enriched_bars(cluster_3_enrichement,"go","Mixed Types\ncluster 3: GO-BP")
c4_go_bar = draw_enriched_bars(cluster_4_enrichement,"go","Immune Infiltrated\ncluster 4: GO-BP")

c0_kegg_bar = draw_enriched_bars(cluster_0_enrichement,"kegg", "Parenchymal Breast\ncluster 0: KEGG")
c1_kegg_bar = draw_enriched_bars(cluster_1_enrichement,"kegg", "Cutaneous Breast\ncluster 1: KEGG")
c2_kegg_bar = draw_enriched_bars(cluster_2_enrichement,"kegg", "HNFS\ncluster 2: KEGG")
c3_kegg_bar = draw_enriched_bars(cluster_3_enrichement,"kegg", "Mixed Types\ncluster 3: KEGG")
c4_kegg_bar = draw_enriched_bars(cluster_4_enrichement,"kegg", "Immune Infiltrated\ncluster 4: KEGG")

c0_msigdb_bar = draw_enriched_bars(cluster_0_enrichement,"msigdb_c6", "Parenchymal Breast\ncluster 0: Oncogenic Signatures")
c1_msigdb_bar = draw_enriched_bars(cluster_1_enrichement,"msigdb_c6", "Cutaneous Breast\ncluster 1: Oncogenic Signatures")
c2_msigdb_bar = draw_enriched_bars(cluster_2_enrichement,"msigdb_c6", "HNFS\ncluster 2: Oncogenic Signatures")
c3_msigdb_bar = draw_enriched_bars(cluster_3_enrichement,"msigdb_c6", "Mixed Types\ncluster 3: Oncogenic Signatures")
c4_msigdb_bar = draw_enriched_bars(cluster_4_enrichement,"msigdb_c6", "Immune Infiltrated\ncluster 4: Oncogenic Signatures")

go_kegg_msigdb_bars_combined = cowplot::plot_grid(
  c0_go_bar,c1_go_bar,c2_go_bar,c3_go_bar,c4_go_bar,
  c0_kegg_bar,c1_kegg_bar,c2_kegg_bar,c3_kegg_bar,c4_kegg_bar,
  c0_msigdb_bar,c1_msigdb_bar,c2_msigdb_bar,c4_msigdb_bar,c3_msigdb_bar,
  labels = "AUTO",
  ncol = 5
)
go_kegg_msigdb_bars_combined
ggsave(go_kegg_msigdb_bars_combined,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_cluster_top_markers_go_kegg_msigdb_enrichments.png",dpi=300,width=24,height=16)



## Check endothelial migration genes
top_markers_dimplot = FeaturePlot(so, features = c(
  "FLT1","APLNR","EBF3","NRP1",
  "MAP4K2","SMYD2","VAV3","NRP2"
),pt.size = 2,ncol = 4,reduction = "pca")
top_markers_dimplot

## Make biplot between marker genes
plot_biplot = function(so,gene1,gene2,color_attr) {
  gene_subset = so@assays$RNA$vst_scaled[c(gene1,gene2),]
  gene_subset = as.data.frame(t(gene_subset),check.names=FALSE)
  gene_subset$og_id = rownames(gene_subset)
  gene_subset = merge(gene_subset,so@meta.data,on="og_id",how="left")
  genes_biplot = ggplot(gene_subset,aes_string(x=gene1,y=gene2)) +
    geom_point(size=4,aes(color=gene_subset[[color_attr]])) +
    theme_minimal() +
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0)# +
    #geom_smooth(method='lm')
  return(genes_biplot)
}

## NRP1-NRP2
nrp1_nrp2_biplot_by_site = plot_biplot(so,"NRP1","NRP2","PRIMARY SITE (Combined)")
nrp1_nrp2_biplot_by_site

## NRP1-MYC axis separates major clusters
nrp_myc_biplot_by_site = plot_biplot(so,"NRP1","MYC","PRIMARY SITE (Combined)")
ggsave(nrp_myc_biplot_by_site,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_NRP1_MYC_biplot.png",dpi=300,width=24,height=16)

## TERT-FLT1 gets parenchymal breast samples
tert_flt1_biplot_by_site = plot_biplot(so,"TERT","FLT1","PRIMARY SITE (Combined)")
ggsave(tert_flt1_biplot_by_site,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_TERT_FLT1_biplot.png",dpi=300,width=24,height=16)

## MYC-CTLA4 gets cutaneous breast samples
myc_ctla4_biplot = plot_biplot(so,"CTLA4","MYC","PRIMARY SITE (Combined)")
ggsave(myc_ctla4_biplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_MYC_CTLA4_biplot.png",dpi=300,width=24,height=16)

## FGFR3?
plot_biplot(so,"FGFR3","APOBEC3G","PRIMARY SITE (Combined)")
plot_biplot(so,"FGFR3","PRKCA","PRIMARY SITE (Combined)")


plot_biplot(so,"SEMA4D","MAPK1","PRIMARY SITE (Combined)")
plot_biplot(so,"MAPK1","TBX1","PRIMARY SITE (Combined)")
flt1_foxp3_biplot = plot_biplot(so,"FLT1","FOXP3","PRIMARY SITE (Combined)")
ctla4_foxp3_biplot = plot_biplot(so,"CTLA4","FOXP3","PRIMARY SITE (Combined)")
fgfr3_foxp3_biplot = plot_biplot(so,"FGFR3","FOXP3","PRIMARY SITE (Combined)")

ggsave(flt1_foxp3_biplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_FLT1_FOXP3_biplot.png",dpi=300,width=24,height=16)
ggsave(ctla4_foxp3_biplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_CTLA4_FOXP3_biplot.png",dpi=300,width=24,height=16)
ggsave(fgfr3_foxp3_biplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_FGFR3_FOXP3_biplot.png",dpi=300,width=24,height=16)


VlnPlot(so,features = c(
  "VEGFA","VEGFB","VEGFC","VEGFD",
  "PDGFA","PDGFB","PDGFC",
  "FLT1","KDR","FLT3","FLT4",
  #"TIE1","MYC",
  #"DAXX","ATRX","ARID1A","POT1","TERT"
  ),slot="vst_scaled"
)

so

curated_heatmap = DoHeatmap(so,features=c(
  "FLT1","NRP1","FGFR1","APLNR","PDGFA","PDGFB","TERT","POT1","ATRX","DAXX",
  "KDR","NRP2","TBX1",
  "TIE1","FLT4","MYC",
  "RAC1","CTLA4","BCR","TGFB1",#"IL2RA","ATR",
  "CHEK1", "CHEK2", "MAD2L1", "BUB1B", "TP53", "BRCA1", "BRCA2", "RAD51", "NBN","MDM2","SMYD3","FAT1",
  "KRAS","HRAS","NRAS",
  "ERBB2","ERBB3","EGFR","MET","ALK","RET","PDGFRA","KIT","HLA-A","HLA-B","HLA-C",
  "FGFR2","FGFR3","FGFR4","LRP5","ANGPT2","TEK",
  "MAPK3","MAPK1",
  #"ERK1","ERK2","MAPK","LRP6","ANG2","TIE2",
  "FLT3","PDGFC",
  "PIK3CA","PTPRB","PTPRD","PLCG1",
  "VEGFA","VEGFB","VEGFC",
  "PDCD1","FOXP3","CXCL10","CCR5",
  "CD19","MS4A1","IL10","LAG3","CCL2","IL17A","TNF"
),slot = "vst_scaled",disp.min=-2.5,disp.max=2.5)
curated_heatmap
ggsave(curated_heatmap,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_curated_genes_heatmap.png",dpi=300,width=24,height=16)


## MYC-CTLA4 axis
plot_biplot(so,"KRAS","HRAS","PRIMARY SITE (Combined)")
plot_biplot(so,"MYC","MAPK1","PRIMARY SITE (Combined)")

plot_biplot(so,"KRAS","PERP","PRIMARY SITE (Combined)")

plot_biplot(so,"DLL4","KDR","PRIMARY SITE (Combined)")
plot_biplot(so,"DLL4","FLT1","PRIMARY SITE (Combined)")
plot_biplot(so,"DLL4","NOTCH1","PRIMARY SITE (Combined)")

plot_biplot(so,"VEGFA","FLT1","PRIMARY SITE (Combined)")



## MYC-TERT axis
myc_tert_biplot = plot_biplot(so,"TERT","KDR","PRIMARY SITE (Combined)")
myc_tert_biplot

## MYC-CTLA4-KRAS axis
myc_kras_biplot_by_site = plot_biplot(so,"MYC","KRAS","PRIMARY SITE (Combined)")
ctla4_myc_biplot_by_site = plot_biplot(so,"MYC","CTLA4","PRIMARY SITE (Combined)")
plot_grid(myc_kras_biplot_by_site,ctla4_myc_biplot_by_site)

tp53_smyd2_biplot_by_site = plot_biplot(so,"UGT2A2","PNPLA5","PRIMARY SITE (Combined)")
tp53_smyd2_biplot_by_site

## CTLA4-CLEC2A (cutaneous vs, rest)
tp53_smyd2_biplot_by_site = plot_biplot(so,"CTLA4","CD274","PRIMARY SITE (Combined)")
tp53_smyd2_biplot_by_site 

plot_biplot(so,"CDHR1","CHP2","RAAS_LAAS_Class")
plot_biplot(so,"NRP1","NRP2","RAAS_LAAS_Class")

## Confirm relation between VEGF receptors (FLT/NRP)
nrp_biplot_by_site = plot_biplot(so,"NRP1","NRP2","PRIMARY SITE (Combined)")
nrp_biplot_by_raas = plot_biplot(so,"NRP1","NRP2","RAAS (EHR_EXTRACTED)")
plot_grid(nrp_biplot_by_site,nrp_biplot_by_raas)


flt_biplot_by_site = plot_biplot(so,"TERT","MYC","PRIMARY SITE (Combined)")
flt_biplot_by_site

nrp_subset = so@assays$RNA$vst_scaled[c("NRP1","NRP2","VEGFA","KDR"),]
nrp_subset = as.data.frame(t(nrp_subset),check.names=FALSE)
nrp_subset$og_id = rownames(nrp_subset)
nrp_subset = merge(nrp_subset,so@meta.data,on="og_id",how="left")
color_attr = "PRIMARY SITE (Combined)"

nrp_biplot = ggplot(nrp_subset,aes(x=NRP1,y=NRP2)) +
  geom_point(size=3,aes(color=nrp_subset[[color_attr]])) +
  theme_minimal() +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) +
  geom_smooth(method='lm')
nrp_biplot

## Write the seurat derived metadata
write.table(so@meta.data,"data/processed/rna/02_seurat_processed_metadata.tsv",sep="\t")


## Show mutations overlaid on top of PCA
mutation_dimplot = DimPlot(so, reduction = "pca", group.by=c(
  "POT1_Mut_Germline","POT1_Mut_Onehot","TERT_Mut_Onehot","ATRX_Mut_Onehot",
  "FLT1_Mut_Onehot","KDR_Mut_Onehot","FLT3_Mut_Onehot","FLT4_Mut_Onehot",
  "POLE_Mut_Onehot","TP53_Mut_Onehot","PIK3CA_Mut_Onehot",
  "PLCG1_Mut_Onehot","PTPRB_Mut_Onehot",
  "KRAS_Mut_Onehot","NRAS_Mut_Onehot","HRAS_Mut_Onehot"
  ),pt.size=3)
mutation_dimplot
ggsave(mutation_dimplot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_mutation_status_on_pcs.png",dpi=300,width=24,height=12)


## Visualize the top enriched genes in each cluster
top_markers_plot = VlnPlot(so, features = c(
  "FLT1","APLNR","EBF3","NRP1",
  "MAP4K2","SMYD2","VAV3","NRP2",
  "KRTAP24-1","CCDC172","UGT2A2","PNPLA5",
  "ST6GALNAC5","C7","NCR3LG1","PTGS2",
  "KLRF2","LCE5A","CLEC2A","SP8"
  #"CARD18","SERPINB12","PGLYRP3","SEC14L1","IQCJ-SCHIP1","TSPAN18",
  #"ZNF462","BMPR1A","GLI3","TRIB3","CDC25B","DUSP5"
  ),pt.size = 2,ncol = 4)
top_markers_plot
ggsave(top_markers_plot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_top_markers_violin_plot.png",dpi=300,width=24,height=12)

## The same content but on PCs
top_markers_dimplot = FeaturePlot(so, features = c(
  "FLT1","APLNR","EBF3","NRP1",
  "MAP4K2","SMYD2","VAV3","NRP2",
  "KRTAP24-1","CCDC172","UGT2A2","PNPLA5",
  "ST6GALNAC5","C7","NCR3LG1","PTGS2",
  "KLRF2","LCE5A","CLEC2A","SP8"
  #"CARD18","SERPINB12","PGLYRP3","SEC14L1","IQCJ-SCHIP1","TSPAN18",
  #"ZNF462","BMPR1A","GLI3","TRIB3","CDC25B","DUSP5"
),pt.size = 2,ncol = 4,reduction = "pca")
top_markers_dimplot
ggsave(top_markers_dimplot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_top_markers_pca_plot.png",dpi=300,width=24,height=12)

## Show Continuous Feature overlaid on top of PCA
continuous_clinical_dimplot = FeaturePlot(so, reduction = "pca", features=c("Age (Combined)","tmb"),pt.size=3)
continuous_clinical_dimplot
ggsave(continuous_clinical_dimplot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_continuous_clinical_feature_pca_plot.png",dpi=300,width=24,height=12)

## Show gene expression overlaid on top of PCA
genes_to_extract = c("MYC","TP53","PLCG1","PTPRB","POT1","HRAS","NRAS","KRAS","TERT","PIK3CA")#,"PIK3CA","TERT")  # Replace with your genes of interest
expr_dim_plot = FeaturePlot(so,reduction = "pca", features=c(
  "FLT1","KDR","FLT4","FLT3",
  "PGF","PIGF","YAP1","MIB2",
  "MYC","TP53","PLCG1","PTPRB","POT1","RAD51",
  "HRAS","NRAS","KRAS","TERT","PIK3CA","GLI3","FAM83G","NFKB1"
),slot="vst_scaled",pt.size=2)#,cols = c("blue","red"))
expr_dim_plot
ggsave(expr_dim_plot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_gene_expression_pca_plot.png",dpi=300,width=24,height=12)

## Combine all previous information and combine them into heatmap
top_cluster_markers = all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.3) %>%
  #dplyr::filter(!grepl("KRT",gene)) %>%
  slice_head(n = 20) %>%
  ungroup()

top_cluster_markers_heatmap = DoHeatmap(so, features = top_cluster_markers$gene,slot = "vst.scaled",disp.min=-3,disp.max=3)
top_cluster_markers_heatmap
ggsave(top_cluster_markers_heatmap,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_top_cluster_marker_heatmap.png",dpi=300,width=24,height=12)

## Save the seurat derived clusters as metadata
#write.csv("seurat")
write.csv(so@meta.data,"reference_data/RNA_Seq/outputs/seurat/seurat_processed_metadata.csv")





