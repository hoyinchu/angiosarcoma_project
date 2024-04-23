library(Seurat)
library(cowplot)
library(ggplot2)
library(tidyr)
library(dplyr)

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

## If available, add tumor WES derived metadata
meta_data_add_on = read.csv("data/processed/tumor_WES/00b_ASC_mutations_metadata_extended.csv",check.names = FALSE)
add_cols = c("Total_Mutations","TMB","SBS1","SBS7a","SBS7b","SBS10b","SBS15","SBS87","SBS1_rel","SBS7a_rel","SBS7b_rel","SBS10b_rel","SBS15_rel","SBS87_rel")
meta_data_add_on_subset = meta_data_add_on[,c("sample_alias",add_cols)] %>% drop_na()
meta_data = merge(meta_data,meta_data_add_on_subset,on="sample_alias",all.x=TRUE)

# Indicate whether the sample is blood based
#meta_data$is_blood = grepl("BLOOD",meta_data$sample_alias)
# 
# ## Add in GSVA scores (cancer)
# gsva_cancer_data = read.csv("reference_data/RNA_Seq/outputs/gsva/gsva_ssGSEA_score_vst_cancer_set.csv",row.names = 1,check.names = FALSE)
# gsva_cancer_data = as.data.frame(t(gsva_cancer_data),check.names=FALSE)
# gsva_cancer_data$id = rownames(gsva_cancer_data)
# meta_data = merge(meta_data,gsva_cancer_data,on="id")
# 
# ## Add in GSVA scores (GOBP)
# gsva_gobp_data = read.csv("reference_data/RNA_Seq/outputs/gsva/gsva_ssGSEA_score_vst_GO_BP_set.csv",row.names = 1,check.names = FALSE)
# gsva_gobp_data = as.data.frame(t(gsva_gobp_data),check.names=FALSE)
# gsva_gobp_data$id = rownames(gsva_gobp_data)
# meta_data = merge(meta_data,gsva_gobp_data,on="id")
# 
# ## Add in GSVA scores (Hallmark)
# gsva_hallmark_data = read.csv("reference_data/RNA_Seq/outputs/gsva/gsva_ssGSEA_score_vst_hallmark_set.csv",row.names = 1,check.names = FALSE)
# gsva_hallmark_data = as.data.frame(t(gsva_hallmark_data),check.names=FALSE)
# gsva_hallmark_data$id = rownames(gsva_hallmark_data)
# meta_data = merge(meta_data,gsva_hallmark_data,on="id")

meta_data[["seurat_id"]] = meta_data$`entity:sample_id`
meta_data[["seurat_id"]] = gsub("_", ".", meta_data[["seurat_id"]])
meta_data[["seurat_id"]] = gsub("-", ".", meta_data[["seurat_id"]])
meta_data[meta_data == ""] = NA
rownames(meta_data) = meta_data$seurat_id

# Create a seurat object and add meta data
#o = CreateSeuratObject(counts = tpm_data,min.cells = 0, min.features = 0, meta.data = meta_data)
so = CreateSeuratObject(counts = vst_data,min.cells = 0, min.features = 0, meta.data = meta_data)

so@meta.data[["og_id"]] = rownames(so@meta.data)
#so@meta.data = merge(so@meta.data,meta_data,by.x="og_id",by.y="seurat_id",all.x=TRUE)
#rownames(so@meta.data) = so@meta.data[["og_id"]]

so[["RNA"]]$vst_scaled = vst_data_scaled
so = SetIdent(so, value = rownames(so@meta.data))

# Using rocanja's recommendation in handling TPM data in seurat
# reference: https://github.com/satijalab/seurat/issues/7496
#so = SetAssayData(object=so, layer="data", new.data = log1p(tpm_data))
so = SetAssayData(object=so, layer="data", new.data = vst_data)

so = FindVariableFeatures(so,selection.method = "vst", nfeatures = 2000)
so = ScaleData(so)
so = RunPCA(so, features = VariableFeatures(object = so))
## Run unsuperivsed clustering
so = FindNeighbors(so, dims = 1:10)
so = FindClusters(so, resolution = 1.5)
#so = FindClusters(so, resolution = 3)
so = RunUMAP(so, dims = 1:10)

## TODO: subset to include only the single most representative sample from each patient?

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
ggsave(pca_combined_plot,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_pca_plot_by_batch_and_site.png",dpi=300,width = 16,height=8)

## Plot PCA Elbow
elbow_plot = ElbowPlot(so,reduction="pca")
elbow_plot
ggsave(elbow_plot,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_pca_elbow_plot.png",dpi=300)

## Check the mean-variance relationships
top10_variable_genes = head(VariableFeatures(so),10)
mean_var_plot = VariableFeaturePlot(so)
mean_var_plot = LabelPoints(plot=mean_var_plot,points=top10_variable_genes,repel=TRUE)
mean_var_plot
ggsave(mean_var_plot,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_mean_var_plot.png",dpi=300)

## Visualize the top loadings
gene_loading_plot = VizDimLoadings(so, dims = 1:6, reduction = "pca",nfeatures=15)
gene_loading_plot
ggsave(gene_loading_plot,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_gene_loadings_plot.png",dpi=300,width=8,height=12)

## Visualize the genes on the two extreme ends of each PCs
dim_score_plot = DimHeatmap(so, dims = 1:15,balanced = TRUE,combine = TRUE,fast=FALSE)
dim_score_plot
ggsave(dim_score_plot,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_dim_heatmap_plot.png",dpi=300,width=8,height=12)

dim_score_first_3_pcs_plot = DimHeatmap(so, dims = 1:3,balanced = TRUE,combine = TRUE,fast=FALSE,nfeatures=50,reduction="pca")
dim_score_first_3_pcs_plot
ggsave(dim_score_first_3_pcs_plot,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_dim_heatmap_first_3_pcs_plot.png",dpi=300,width=8,height=12)

## Show clinical metadata overlaid on top of PCA
clinical_dimplot = DimPlot(so,reduction="pca",group.by = c(
  "seurat_clusters","PRIMARY SITE (Combined)","CUTANEOUS AS (EHR_EXTRACTED)","sex_combined",
  "BX_SPINDLE_CELL","BX_NUCLEAR_GRADE","BX_VASOFORMATIVE","BX_EPITHELIOID",
  "RAAS_LAAS_Class","CUTANEOUS AS (EHR_EXTRACTED)"),pt.size=3)
clinical_dimplot
#ggsave(clinical_dimplot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_clinical_attributes_on_pcs.png",dpi=300,width=24,height=12)
ggsave(clinical_dimplot,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_clinical_attributes_on_pcs_res15.png",dpi=300,width=24,height=12)

## If available, also plot continuous features on top of PCA
clinical_feat_plot = FeaturePlot(so, features = c(
  "TMB","SBS1_rel","SBS7a_rel","SBS7b_rel","SBS10b_rel","SBS15_rel","SBS87_rel"
),pt.size = 3,ncol = 4,reduction = "pca")
ggsave(clinical_feat_plot,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_tmb_and_sbs_on_pcs_res15.png",dpi=300,width=24,height=12)


## Find markers associated with each cluster
Idents(so) = so@meta.data$seurat_clusters
all_markers = FindAllMarkers(
  so,
  slot="counts",
  test.use = "wilcox",
  only.pos = TRUE
)
## Save the output
write.csv(all_markers,"RNA_Analysis/outputs/DEGs/02_Seurat_all_markers.csv")

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
ggsave(top_cluster_markers_heatmap,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_top_cluster_marker_heatmap.png",dpi=300,width=24,height=12)

## overlap the top 5 in the PCs
top_5_cluster_markers = all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(p_val_adj < 0.1) %>%
  dplyr::filter(avg_log2FC > 2) %>%
  #dplyr::filter(!grepl("KRT",gene)) %>%
  slice_head(n = 5) %>%
  ungroup()

## The same content but on PCs
top_markers_dimplot = FeaturePlot(so, features = top_5_cluster_markers$gene,pt.size = 2,ncol = 5,reduction = "pca")
top_markers_dimplot
ggsave(top_markers_dimplot,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_top_5_markers_per_cluster_on_pca_plot.png",dpi=300,width=24,height=12)


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
ggsave(go_kegg_msigdb_bars_combined,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_cluster_top_markers_go_kegg_msigdb_enrichments.png",dpi=300,width=24,height=16)


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

## NRP1-MYC axis separates major clusters
nrp_myc_biplot_by_site = plot_biplot(so,"NRP1","MYC","PRIMARY SITE (Combined)")
ggsave(nrp_myc_biplot_by_site,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_NRP1_MYC_biplot.png",dpi=300,width=24,height=16)

## TERT-FLT1 gets parenchymal breast samples
tert_flt1_biplot_by_site = plot_biplot(so,"TERT","FLT1","PRIMARY SITE (Combined)")
ggsave(tert_flt1_biplot_by_site,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_TERT_FLT1_biplot.png",dpi=300,width=24,height=16)

## MYC-CTLA4 gets cutaneous breast samples
myc_ctla4_biplot = plot_biplot(so,"CTLA4","MYC","PRIMARY SITE (Combined)")
ggsave(myc_ctla4_biplot,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_MYC_CTLA4_biplot.png",dpi=300,width=24,height=16)

## FGFR3?
plot_biplot(so,"FGFR3","APOBEC3G","PRIMARY SITE (Combined)")
plot_biplot(so,"FGFR3","PRKCA","PRIMARY SITE (Combined)")


plot_biplot(so,"SEMA4D","MAPK1","PRIMARY SITE (Combined)")
plot_biplot(so,"MAPK1","TBX1","PRIMARY SITE (Combined)")
flt1_foxp3_biplot = plot_biplot(so,"FLT1","FOXP3","PRIMARY SITE (Combined)")
ctla4_foxp3_biplot = plot_biplot(so,"CTLA4","FOXP3","PRIMARY SITE (Combined)")
fgfr3_foxp3_biplot = plot_biplot(so,"FGFR3","FOXP3","PRIMARY SITE (Combined)")

ggsave(flt1_foxp3_biplot,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_FLT1_FOXP3_biplot.png",dpi=300,width=24,height=16)
ggsave(ctla4_foxp3_biplot,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_CTLA4_FOXP3_biplot.png",dpi=300,width=24,height=16)
ggsave(fgfr3_foxp3_biplot,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_FGFR3_FOXP3_biplot.png",dpi=300,width=24,height=16)


VlnPlot(so,features = c(
  "VEGFA","VEGFB","VEGFC","VEGFD",
  "PDGFA","PDGFB","PDGFC",
  "FLT1","KDR","FLT3","FLT4"
  )
)

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
ggsave(curated_heatmap,filename = "RNA_Analysis/outputs/plots/02_seurat_plots/02_curated_genes_heatmap.png",dpi=300,width=24,height=16)


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



# ## Given the high display of keratin related markers, we filter out markers that are significant when compared  between cluster 1 and 4 
# hnfs_clust_vs_cut_breast = FindMarkers(
#   so,
#   #cells.1 = 2,#WhichCells(so,idents=2),
#   #cells.2 = 4,#WhichCells(so,idents=4),
#   ident.1 = 2,
#   ident.2 = 4,
#   slot="counts",
#   test.use = "wilcox",
#   only.pos=TRUE,
#   )
# hnfs_clust_vs_cut_breast_sig = hnfs_clust_vs_cut_breast %>% dplyr::filter(p_val_adj < 0.05)
# #remove_genes = 

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



## Overlay enriched pathways in path view
# Back-translate to gene names
#enriched_uniprot = cluster_0_enrichement$kegg@
#bitr(sig_genes,fromType="SYMBOL",toType="UNIPROT",OrgDb = org.Hs.eg.db)$UNIPROT 
# 
# library(pathview)
# sub_df = all_markers[all_markers$cluster==0,]
# sig_genes = rownames(sub_df[sub_df$p_val < 0.05,])
# enriched_path_viz = pathview(gene.data=sig_genes, pathway.id="hsa04370", species = "hsa",gene.idtype="SYMBOL",out.suffix = "cluster0")
# enriched_path_viz = pathview(gene.data=sig_genes, pathway.id="hsa04080", species = "hsa",gene.idtype="SYMBOL",out.suffix = "cluster0")
# enriched_path_viz = pathview(gene.data=sig_genes, pathway.id="hsa04024", species = "hsa",gene.idtype="SYMBOL",out.suffix = "cluster0")
# enriched_path_viz = pathview(gene.data=sig_genes, pathway.id="hsa04015", species = "hsa",gene.idtype="SYMBOL",out.suffix = "cluster0.rap_signaling")
# enriched_path_viz = pathview(gene.data=sig_genes, pathway.id="hsa04218", species = "hsa",gene.idtype="SYMBOL",out.suffix = "cluster0.senescence")
# enriched_path_viz = pathview(gene.data=sig_genes, pathway.id="hsa04514", species = "hsa",gene.idtype="SYMBOL",out.suffix = "cluster0.adhesion")

#browseKEGG(cluster_0_enrichement$kegg,pathID="hsa04370")


## Overlay expression of genes from difference endothelial / angiogenesis pathways?

# Show some gsva scores that are correlated with PCs
gsva_bp_clinical_dimplot = FeaturePlot(
  so,reduction="pca",features=c(
    "GOBP_SKIN_EPIDERMIS_DEVELOPMENT","GOBP_RAP_PROTEIN_SIGNAL_TRANSDUCTION","GOBP_MORPHOGENESIS_OF_AN_ENDOTHELIUM",
    "GOBP_NEGATIVE_REGULATION_OF_CELLULAR_SENESCENCE","GOBP_NEGATIVE_REGULATION_OF_CELL_MIGRATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS","GOBP_ENDOTHELIAL_CELL_MATRIX_ADHESION",
    "GOBP_POSITIVE_REGULATION_OF_FIBROBLAST_MIGRATION","GOBP_POSITIVE_REGULATION_OF_CHROMATIN_BINDING","GOBP_REGULATION_OF_PLATELET_DERIVED_GROWTH_FACTOR_RECEPTOR_ALPHA_SIGNALING_PATHWAY",
    "GOBP_RENAL_SYSTEM_VASCULATURE_MORPHOGENESIS","GOBP_OVULATION_CYCLE","GOBP_SMOOTH_MUSCLE_TISSUE_DEVELOPMENT",
    "GOBP_MEIOTIC_CELL_CYCLE_PHASE_TRANSITION","GOBP_POSITIVE_REGULATION_OF_MITOTIC_CYTOKINESIS","GOBP_G2_MI_TRANSITION_OF_MEIOTIC_CELL_CYCLE",
    "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_BREAK_INDUCED_REPLICATION"
  ),pt.size=3
)
gsva_bp_clinical_dimplot
ggsave(gsva_bp_clinical_dimplot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_gsva_bp_corr_on_pcs.png",dpi=300,width=36,height=18)

## Show some cancer related gsva scores that correlate with PCs
gsva_cancer_clinical_dimplot = FeaturePlot(
  so,reduction="pca",features=c(
    "MYC","MYC_UP.V1_DN","MYC_UP.V1_UP",
    "KRAS","KRAS.600_UP.V1_DN","KRAS.600_UP.V1_UP",
    "KDR", "VEGF_A_UP.V1_DN", "VEGF_A_UP.V1_UP",
    "HRAS","TERT","POT1"
    #"SINGH_KRAS_DEPENDENCY_SIGNATURE","KRAS.50_UP.V1_DN","ESC_J1_UP_EARLY.V1_UP","PKCA_DN.V1_UP",
    #"CSR_EARLY_UP.V1_UP","E2F1_UP.V1_UP","BCAT_BILD_ET_AL_DN","CSR_LATE_UP.V1_UP",
    #"VEGF_A_UP.V1_UP","VEGF_A_UP.V1_DN","MYC_UP.V1_UP","PKCA_DN.V1_UP","YAP1_UP"
    #"JNK_DN.V1_UP","CORDENONSI_YAP_CONSERVED_SIGNATURE",
    #"RAF_UP.V1_UP",
  ),pt.size = 3,
  ncol=3
  #cols=c("blue","red")#,keep.scale="all"
)
gsva_cancer_clinical_dimplot

## Check POT1 and TERT expression by mutation status
Idents(so) = so@meta.data$POT1_Mut_Onehot
so_sub = subset(x = so, subset = primary_site_combined == "BREAST (PARENCHYMAL)")
VlnPlot(so,features=c("TERT","POT1","TRF1","TRF2",
                      "ATRX","DAXX","BRCA1","BRCA2",
                      "MRE11","RAD50","NBS1",
                      "BLM","RAD51","RAD52","WRN","FEN1",
                      "SMARCB1"
                      ),layer="data",pt.size = 3)

## Show some hallmark pathway set related gsva scores that correlate with PCs
gsva_hallmark_clinical_dimplot = FeaturePlot(
  so,reduction="pca",features=c(
    "HALLMARK_KRAS_SIGNALING_UP","HALLMARK_ESTROGEN_RESPONSE_LATE","HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_MITOTIC_SPINDLE","HALLMARK_G2M_CHECKPOINT","HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_UV_RESPONSE_UP","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_COAGULATION",
    "HALLMARK_MYC_TARGETS_V2","HALLMARK_E2F_TARGETS","HALLMARK_ANGIOGENESIS"
  ),pt.size = 3,ncol=3
  #cols=c("blue","red")#,kemep.scale="all"
)
gsva_hallmark_clinical_dimplot

## Show the relative log TPM expression of known angiosarcoma related genes
angiosarcoma_genes_dimplot= FeaturePlot(
  so,reduction="pca",features=c(
    "MYC","KDR","PLCG1","PTPRB","POT1","HRAS","KRAS","NRAS"
  ),pt.size = 3
  #cols=c("blue","red")#,keep.scale="all"
)
angiosarcoma_genes_dimplot


# gsva_up_clinical_dimplot = FeaturePlot(
#   so, reduction = "pca", features=c(
#     "KRAS.300_UP.V1_UP","KRAS.600.LUNG.BREAST_UP.V1_UP","KRAS.LUNG.BREAST_UP.V1_UP",
#     "MYC_UP.V1_UP","P53_DN.V1_UP","P53_DN.V2_UP","PDGF_UP.V1_UP","PGF_UP.V1_UP","PKCA_DN.V1_UP",
#     "SINGH_KRAS_DEPENDENCY_SIGNATURE","SRC_UP.V1_UP","VEGF_A_UP.V1_UP","YAP1_UP"
#     ),pt.size=3)
# gsva_clinical_dimplot
# 
# gsva_dn_clinical_dimplot = FeaturePlot(
#   so, reduction = "pca", features=c(
#     "KRAS.300_UP.V1_DN","KRAS.600.LUNG.BREAST_UP.V1_DN","KRAS.LUNG.BREAST_UP.V1_DN",
#     "MYC_UP.V1_DN","P53_DN.V1_DN","P53_DN.V2_DN","PDGF_UP.V1_DN","PGF_UP.V1_DN","PKCA_DN.V1_DN",
#     "SINGH_KRAS_DEPENDENCY_SIGNATURE","SRC_UP.V1_DN","VEGF_A_UP.V1_DN","YAP1_DN"
#   ),pt.size=3)
# gsva_dn_clinical_dimplot

## Check which features correlate with PCs
gobp_columns <- grep("^GOBP", names(so@meta.data), value = TRUE)
cancer_sig_columns <- colnames(gsva_cancer_data)
cancer_sig_columns <- cancer_sig_columns[cancer_sig_columns!="id"]
hallmark_columns <- colnames(gsva_hallmark_data)
hallmark_columns <- hallmark_columns[hallmark_columns!="id"]


# Define the function
correlatePCs <- function(seuratObj,columns, pc_num) {
  # Ensure the object has PCA results
    # Extract PC1 scores
    pc1_scores <- so@reductions$pca@cell.embeddings[, pc_num]
    
    # Filter metadata columns that start with "GOBP"
    #go_columns <- grep("^GOBP", names(seuratObj@meta.data), value = TRUE)
    
    # Initialize a list to store correlation results
    correlations <- list()
    # Loop through the GOBP columns and calculate correlation with PC1
    for(col in columns) {
      # Ensure that the column is numeric
      if(is.numeric(seuratObj@meta.data[[col]])) {
        # Calculate correlation
        cor_result <- cor.test(pc1_scores, seuratObj@meta.data[[col]], method = "pearson")
        
        # Store the correlation coefficient and p-value
        correlations[[col]] <- list(correlation = cor_result$estimate,
                                    p_value = cor_result$p.value)
      } else {
        warning(paste("Column", col, "is not numeric and was skipped."))
      }
    }
    # Return the list of correlations
    return(correlations)
}

convertCorListToDF <- function(cor_list) {
  # Extracting correlations, p-values, and names
  correlations <- sapply(cor_list, function(x) x$correlation)
  p_values <- sapply(cor_list, function(x) x$p_value)
  names_vector <- names(cor_list)
  # Create a dataframe from the extracted values
  cor_df <- data.frame(
    Names = names_vector,
    Correlation = correlations,
    P_Value = p_values,
    stringsAsFactors = FALSE  # To keep strings as character type
  )
  return(cor_df)
}

## Create the correlation dfs
gobp_corr_pc1 = correlatePCs(so,gobp_columns,1)
gobp_corr_pc2 = correlatePCs(so,gobp_columns,2)
cancer_corr_pc1 = correlatePCs(so,cancer_sig_columns,1)
cancer_corr_pc2 = correlatePCs(so,cancer_sig_columns,2)
hallmark_corr_pc1 = correlatePCs(so,hallmark_columns,1)
hallmark_corr_pc2 = correlatePCs(so,hallmark_columns,2)

## Check GO_BP correlations
pc1_corr_df = convertCorListToDF(gobp_corr_pc1)
pc1_corr_df_nom_sig = pc1_corr_df[pc1_corr_df$P_Value<0.05,]
pc1_corr_df_nom_sig = pc1_corr_df_nom_sig[order(pc1_corr_df_nom_sig$P_Value),]
pc1_corr_df_nom_sig_pos = pc1_corr_df_nom_sig[pc1_corr_df_nom_sig$Correlation > 0,]
pc1_corr_df_nom_sig_neg = pc1_corr_df_nom_sig[pc1_corr_df_nom_sig$Correlation < 0,]

pc2_corr_df = convertCorListToDF(gobp_corr_pc2)
pc2_corr_df_nom_sig = pc2_corr_df[pc2_corr_df$P_Value<0.05,]
pc2_corr_df_nom_sig = pc2_corr_df_nom_sig[order(pc2_corr_df_nom_sig$P_Value),]
pc2_corr_df_nom_sig_pos = pc2_corr_df_nom_sig[pc2_corr_df_nom_sig$Correlation > 0,]
pc2_corr_df_nom_sig_neg = pc2_corr_df_nom_sig[pc2_corr_df_nom_sig$Correlation < 0,]

## Check Cancer correlations
pc1_cancer_corr_df = convertCorListToDF(cancer_corr_pc1)
pc1_cancer_corr_df_nom_sig = pc1_cancer_corr_df[pc1_cancer_corr_df$P_Value<0.05,]
pc1_cancer_corr_df_nom_sig = pc1_cancer_corr_df_nom_sig[order(pc1_cancer_corr_df_nom_sig$P_Value),]
pc1_cancer_corr_df_nom_sig_pos = pc1_cancer_corr_df_nom_sig[pc1_cancer_corr_df_nom_sig$Correlation > 0,]
pc1_cancer_corr_df_nom_sig_neg = pc1_cancer_corr_df_nom_sig[pc1_cancer_corr_df_nom_sig$Correlation < 0,]

pc2_cancer_corr_df = convertCorListToDF(cancer_corr_pc2)
pc2_cancer_corr_df_nom_sig = pc2_cancer_corr_df[pc2_cancer_corr_df$P_Value<0.05,]
pc2_cancer_corr_df_nom_sig = pc2_cancer_corr_df_nom_sig[order(pc2_cancer_corr_df_nom_sig$P_Value),]
pc2_cancer_corr_df_nom_sig_pos = pc2_cancer_corr_df_nom_sig[pc2_cancer_corr_df_nom_sig$Correlation > 0,]
pc2_cancer_corr_df_nom_sig_neg = pc2_cancer_corr_df_nom_sig[pc2_cancer_corr_df_nom_sig$Correlation < 0,]

## Check Hallmark correlations
pc1_hallmark_corr_df = convertCorListToDF(hallmark_corr_pc1)
pc1_hallmark_corr_df_nom_sig = pc1_hallmark_corr_df[pc1_hallmark_corr_df$P_Value<0.05,]
pc1_hallmark_corr_df_nom_sig = pc1_hallmark_corr_df_nom_sig[order(pc1_hallmark_corr_df_nom_sig$P_Value),]
pc1_hallmark_corr_df_nom_sig_pos = pc1_hallmark_corr_df_nom_sig[pc1_hallmark_corr_df_nom_sig$Correlation > 0,]
pc1_hallmark_corr_df_nom_sig_neg = pc1_hallmark_corr_df_nom_sig[pc1_hallmark_corr_df_nom_sig$Correlation < 0,]

pc2_hallmark_corr_df = convertCorListToDF(hallmark_corr_pc2)
pc2_hallmark_corr_df_nom_sig = pc2_hallmark_corr_df[pc2_hallmark_corr_df$P_Value<0.05,]
pc2_hallmark_corr_df_nom_sig = pc2_hallmark_corr_df_nom_sig[order(pc2_hallmark_corr_df_nom_sig$P_Value),]
pc2_hallmark_corr_df_nom_sig_pos = pc2_hallmark_corr_df_nom_sig[pc2_hallmark_corr_df_nom_sig$Correlation > 0,]
pc2_hallmark_corr_df_nom_sig_neg = pc2_hallmark_corr_df_nom_sig[pc2_hallmark_corr_df_nom_sig$Correlation < 0,]

