library(Seurat)
library(cowplot)
library(ggplot2)
library(tidyr)
library(dplyr)

setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts")

## Load the previously created seurat data
so = readRDS("data/processed/rna/ASCSeuratObj2025.rds")
#so = readRDS("data/processed/rna/ASCSeuratObj2025_Purity.rds")

## Load the palettes
source("util_scripts/project_palettes.R")

# Plot PCA
pca_plot_by_batch = DimPlot(so, reduction = "pca",group.by="LC_BATCH",pt.size=4)
pca_plot_by_site = DimPlot(so, reduction = "pca",group.by="Primary Site (Recombined)",pt.size=4)

## Subset by site clusters
so_cut_breast = subset(so,subset = seurat_clusters_by_site_str == "Breast (Cutaneous)")
so_noncut_breast = subset(so,subset = seurat_clusters_by_site_str == "Breast (Parenchymal)")
so_hnfs = subset(so,subset = seurat_clusters_by_site_str == "HNFS")
so_heart = subset(so,subset = seurat_clusters_by_site_str == "Heart")
so_extremities = subset(so,subset = seurat_clusters_by_site_str == "Extremities")
so_others = subset(so,subset = seurat_clusters_by_site_str == "Others")

## Perform subcluster PCAs
so_cut_breast = FindVariableFeatures(so_cut_breast)
so_cut_breast = RunPCA(so_cut_breast,features=VariableFeatures(so_cut_breast,nfeatures = 2000),npcs = 2)

so_noncut_breast = FindVariableFeatures(so_noncut_breast)
so_noncut_breast = RunPCA(so_noncut_breast,features=VariableFeatures(so_noncut_breast,nfeatures = 2000),npcs = 2)

so_hnfs = FindVariableFeatures(so_hnfs)
so_hnfs = RunPCA(so_hnfs,features=VariableFeatures(so_hnfs,nfeatures = 2000),npcs = 2)

so_heart = FindVariableFeatures(so_heart)
so_heart = RunPCA(so_heart,features=VariableFeatures(so_heart,nfeatures = 2000),npcs = 2)

so_extremities = FindVariableFeatures(so_extremities)
so_extremities = RunPCA(so_extremities,features=VariableFeatures(so_extremities,nfeatures = 2000),npcs = 2)

so_others = FindVariableFeatures(so_others)
so_others = RunPCA(so_others,features=VariableFeatures(so_others,nfeatures = 2000),npcs = 2)

# so = CellCycleScoring(so, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
# ggsave("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_overall_cell_cycle_pca_plots.png",DimPlot(so,group.by = "Phase",pt.size = 5))

# Get the total variance:
get_var_explained = function(so) {
  mat = GetAssayData(so, assay = "RNA", slot = "scale.data")
  pca = so[["pca"]]
  total_variance = sum(matrixStats::rowVars(mat))
  eigen_values = (pca@stdev)^2  
  var_explained = eigen_values / total_variance
  return(var_explained)
}

## Helper function to make embedding plots
make_embedding_plot = function(df, color_col,title_text, palette=NULL) {
  embed_plot = ggplot(df,aes_string(x="PC_1",y="PC_2",color=color_col)) +
    geom_point(size=2) +
    theme_minimal() +
    #scale_color_npg() +
    labs(x="PC 1",y="PC 2",color=title_text,title=title_text) +
    theme(plot.title = element_text(hjust=0.5),legend.position="top",legend.title=element_blank())
  if (!is.null(palette)) {
    embed_plot = embed_plot + scale_color_manual(values = palette)
  }
  return(embed_plot)
}

## Plot PCA Dims 
plot_pca_by_clin = function(so,filename,partial=FALSE) {
  ## Define clinical information that are to be included
  clin_groups = c(
    #"seurat_clusters_by_expression_str",
    "hcluster_by_expr_str",
    "Primary Site (Recombined)","CUTANEOUS AS (EHR_EXTRACTED)","SEX (EHR_EXTRACTED)",
    "BX_SPINDLE_CELL","BX_NUCLEAR_GRADE","BX_VASOFORMATIVE","BX_EPITHELIOID",
    "RAAS_LAAS_Class","CUTANEOUS AS (EHR_EXTRACTED)","BATCH NUMBER (EHR_EXTRACTED)"
  ) 
  
  ## Grab embeddings
  pc_embeddings = Embeddings(so,reduction="pca")[,c("PC_1","PC_2")]
  clin_attr_subset = so@meta.data[,clin_groups]
  clin_attr_subset_merged = merge(clin_attr_subset,pc_embeddings,by = 'row.names', all.x = TRUE)
  #clin_attr_subset_merged$seurat_clusters_renamed = paste0("",so@meta.data$seurat_clusters_renamed)
  
  ## Rename some features
  clin_attr_subset_merged = clin_attr_subset_merged %>% 
    mutate(cutaneous_viz = case_when(
      `CUTANEOUS AS (EHR_EXTRACTED)` == "0" ~ "Non-cutaneous AS",
      `CUTANEOUS AS (EHR_EXTRACTED)` == "1" ~ "Cutaneous AS",
      TRUE ~ "NA"
    )) %>%
    mutate(sex_viz = case_when(
      `SEX (EHR_EXTRACTED)` == "F" ~ "Female",
      `SEX (EHR_EXTRACTED)` == "M" ~ "Male",
      TRUE~ "NA"
    )) %>%
    mutate(batch_num = as.factor(`BATCH NUMBER (EHR_EXTRACTED)`))
  
  ## Make a plot for each clin group
  ## Make sure the palette script is sourced
  #cluster_plot = make_embedding_plot(clin_attr_subset_merged,"seurat_clusters_by_expression_str","Unsupervised Cluster",seurart_cluster_by_expression_palette)
  cluster_plot = make_embedding_plot(clin_attr_subset_merged,"hcluster_by_expr_str","Unsupervised Cluster",seurart_cluster_by_expression_palette)
  primary_site_plot = make_embedding_plot(clin_attr_subset_merged,"`Primary Site (Recombined)`","Primary Sites",primary_site_palette)
  cutaneous_clin_plot = make_embedding_plot(clin_attr_subset_merged,"cutaneous_viz","Is Cutaneous",cutaneous_palette)
  sex_clin_plot = make_embedding_plot(clin_attr_subset_merged,"sex_viz","Sex",sex_clin_palette)
  bx_spindle_cell_plot = make_embedding_plot(clin_attr_subset_merged,"BX_SPINDLE_CELL","Spindle Cell Status",spindle_cell_palette)
  bx_nuclear_grade_plot = make_embedding_plot(clin_attr_subset_merged,"BX_NUCLEAR_GRADE","Nuclear Grade",nuclear_grade_palette)
  bx_vasoformative_grade_plot = make_embedding_plot(clin_attr_subset_merged,"BX_VASOFORMATIVE","Vasoformative",vasoformative_cell_palette)
  bx_epitheliod_grade_plot = make_embedding_plot(clin_attr_subset_merged,"BX_EPITHELIOID","Epithelioid",epitheliod_palette)
  raas_laas_clin_plot = make_embedding_plot(clin_attr_subset_merged,"RAAS_LAAS_Class","RAAS/LAAS",RAAS_class_palette)
  batch_plot = make_embedding_plot(clin_attr_subset_merged,"batch_num","Batch")
  
  ## Also make for continuous variables
  continuous_attributes = c("Age (Combined)","ESTIMATE_purity")#,"endMT_score1","angiogenesis_score1")#,"fibroblast_score1","lymphangiogenesis_score1","radioresistance_score1","radiation_exposure_score1",)
  
  ## Grab embeddings
  cont_attr_subset = so@meta.data[,continuous_attributes]
  clin_attr_cont_subset_merged = merge(cont_attr_subset,pc_embeddings,by = 'row.names', all.x = TRUE)
  
  ## Make a plot for each continous valued group
  age_plot = make_embedding_plot(clin_attr_cont_subset_merged,"`Age (Combined)`","Age",palette = NULL)
  purity_plot = make_embedding_plot(clin_attr_cont_subset_merged,"ESTIMATE_purity","ESTIMATE Purity",palette = NULL)
  
  if (partial) {
    combined_clin_plot = plot_grid(
      primary_site_plot,cluster_plot,
      ncol = 2
    )
    ggsave(combined_clin_plot,filename=filename,dpi=300,width=16,height=4)
  }
  else {
    ## Combine them
    combined_clin_plot = plot_grid(
      primary_site_plot,cluster_plot,
      cutaneous_clin_plot,sex_clin_plot,
      bx_spindle_cell_plot,bx_nuclear_grade_plot,
      bx_vasoformative_grade_plot,bx_epitheliod_grade_plot,
      raas_laas_clin_plot,batch_plot,
      age_plot,purity_plot,
      ncol = 2
    )
    ggsave(combined_clin_plot,filename=filename,dpi=300,width=14,height=18)
  }
}

## Plot PCA continuous attributes 
plot_pca_by_cont = function(so,filename) {
  ## Define clinical information that are to be included
  continuous_attributes = c("Age (Combined)","ESTIMATE_purity")#,"endMT_score1","angiogenesis_score1")#,"fibroblast_score1","lymphangiogenesis_score1","radioresistance_score1","radiation_exposure_score1",)
  
  ## Grab embeddings
  pc_embeddings = Embeddings(so,reduction="pca")[,c("PC_1","PC_2")]
  cont_attr_subset = so@meta.data[,continuous_attributes]
  clin_attr_subset_merged = merge(cont_attr_subset,pc_embeddings,by = 'row.names', all.x = TRUE)
  
  ## Make a plot for each continous valued group
  age_plot = make_embedding_plot(clin_attr_subset_merged,"`Age (Combined)`","Age",palette = NULL)
  purity_plot = make_embedding_plot(clin_attr_subset_merged,"ESTIMATE_purity","ESTIMATE Purity",palette = NULL)
  #endMT_plot = make_embedding_plot(clin_attr_subset_merged,"endMT_score1","EndMT Score",palette = NULL)
  #angiogenesis_plot = make_embedding_plot(clin_attr_subset_merged,"angiogenesis_score1","Angiogenesis Score",palette = NULL)
  
  ## Combine them
  combined_cont_plot = plot_grid(
    age_plot,purity_plot,
    #endMT_plot,angiogenesis_plot,
    ncol = 2
  )
  ggsave(combined_cont_plot,filename=filename,dpi=300,width=14,height=18)
}

# clinical_dimplot = DimPlot(so,reduction="pca",group.by = clin_groups,pt.size=3)
# clinical_dimplot_2cols = DimPlot(so,reduction="pca",group.by = clin_groups,pt.size=3, ncol = 2)
# clinical_dimplot_2cols

# primary_site_plot = make_embedding_plot(clin_attr_subset_merged,"`Primary Site (Recombined)`") +
#   scale_color_manual(values=primary_site_palette) +
#   labs(x="PC 1",y="PC 2",color="Primary Site")
# cluster_plot = make_embedding_plot(clin_attr_subset_merged,"seurat_clusters_renamed") +
#   scale_color_manual(values=seurart_cluster_palette) +
#   labs(x="PC 1",y="PC 2",color="Unsupervised Clusters")
# 
# combined_primary_clusters = plot_grid(cluster_plot,primary_site_plot)
# combined_primary_clusters
# ggsave("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_Fig2_cluster_primay_site_pc_plot.png",combined_primary_clusters,dpi=300,width=10,height=3)

## Plot PCA Dim
plot_pca_dims = function(so,filename) {
  so_var_exp = get_var_explained(so)
  pca_dim_plot = DimPlot(so, reduction = "pca",group.by=c(
    "LC_BATCH","Primary Site (Recombined)","CUTANEOUS (EHR_EXTRACTED)",
    "BX_VASOFORMATIVE","BX_EPITHELIOID","BX_SPINDLE_CELL","BX_NUCLEAR_GRADE",
    "RAAS_LAAS_Class","RACE (PRD)","SEX (EHR_EXTRACTED)"),pt.size=4)
  pca_dim_plot = pca_dim_plot + 
    xlab(paste0("PC1 (",signif(so_var_exp[1]*100,4),"%)")) +
    ylab(paste0("PC2 (",signif(so_var_exp[2]*100,4),"%)"))
  ggsave(pca_dim_plot,filename = filename,dpi=300,width = 16,height=8)
}

## Plot loadings
plot_pca_loadings = function(so,filename) {
  so_loading_plot = VizDimLoadings(so, dims = 1:2, reduction = "pca",nfeatures=15)
  ggsave(so_loading_plot,filename = filename,dpi=300,width=8,height=12)
}

## Plot loading dims
plot_pca_loading_dims = function(so,filename) {
  ## PC Loadings dims
  so_dim_loadings = DimHeatmap(so, dims = 1:2,balanced = TRUE,combine = TRUE,fast=FALSE,nfeatures=50,reduction="pca")
  ggsave(so_dim_loadings,filename = filename,dpi=300,width=8,height=12)
}

## Viasulize top loading features
plot_top_loading_features = function(so,filename,n=5) {
  loadings = Loadings(so)
  top_n_pc1 = names(loadings[order(loadings[, "PC_1"], decreasing = TRUE)[1:n], "PC_1"])
  bottom_n_pc1 = names(loadings[order(loadings[, "PC_1"], decreasing = FALSE)[1:n], "PC_1"])
  top_n_pc2 = names(loadings[order(loadings[, "PC_2"], decreasing = TRUE)[1:n], "PC_2"])
  bottom_n_pc2 = names(loadings[order(loadings[, "PC_2"], decreasing = FALSE)[1:n], "PC_2"])
  all_features = c(top_n_pc1,bottom_n_pc1,top_n_pc2,bottom_n_pc2)
  so_feature_plot = FeaturePlot(so,all_features,pt.size = 4)
  ggsave(so_feature_plot,filename = filename,dpi=300,width=16,height=12)
}

## Check the mean-variance relationships
plot_mean_variance = function(so,filename) {
  top10_variable_genes = head(VariableFeatures(so),10)
  mean_var_plot = VariableFeaturePlot(so)
  mean_var_plot = LabelPoints(plot=mean_var_plot,points=top10_variable_genes,repel=TRUE)
  ggsave(mean_var_plot,filename = filename,dpi=300)
}

## Check PCA ELbow plot
plot_pca_elbow_plot = function(so,filename) {
  elbow_plot = ElbowPlot(so,reduction="pca")
  ggsave(elbow_plot,filename=filename,dpi=300)
}

# Plot PCA related figures:
plot_pca_figs = function(so,dirname,prefix) {
  dir.create(file.path(dirname))
  # #plot_pca_dims(so,paste0(dirname,"/",prefix,"_pca_dims.png"))
  # plot_pca_by_clin(so,paste0(dirname,"/",prefix,"_pca_dims.png"))
  # plot_pca_by_cont(so,paste0(dirname,"/",prefix,"_pca_cont_meta.png"))
  # plot_pca_loadings(so,paste0(dirname,"/",prefix,"_pca_loadings.png"))
  # plot_pca_loading_dims(so,paste0(dirname,"/",prefix,"_pca_loading_dims.png"))
  # plot_top_loading_features(so,paste0(dirname,"/",prefix,"_pca_top_loading_features.png"))
  # plot_pca_elbow_plot(so,paste0(dirname,"/",prefix,"_pca_elbow_plot.png"))
  # plot_mean_variance(so,paste0(dirname,"/",prefix,"_mean_variance_plot.png"))
  # 
  plot_pca_by_clin(so,paste0(dirname,"/",prefix,"_pca_dims.pdf"))
  plot_pca_by_cont(so,paste0(dirname,"/",prefix,"_pca_cont_meta.pdf"))
  plot_pca_loadings(so,paste0(dirname,"/",prefix,"_pca_loadings.pdf"))
  plot_pca_loading_dims(so,paste0(dirname,"/",prefix,"_pca_loading_dims.pdf"))
  plot_top_loading_features(so,paste0(dirname,"/",prefix,"_pca_top_loading_features.pdf"))
  plot_pca_elbow_plot(so,paste0(dirname,"/",prefix,"_pca_elbow_plot.pdf"))
  plot_mean_variance(so,paste0(dirname,"/",prefix,"_mean_variance_plot.pdf"))
}


plot_pca_figs(so,"02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_overall_pca_plots","02_overall")
plot_pca_figs(so_noncut_breast,"02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_noncutaneous_breast_pca_plots","02_noncutaneous_breast")
plot_pca_figs(so_cut_breast,"02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_cutaneous_breast_pca_plots","02_cutaneous_breast")
plot_pca_figs(so_hnfs,"02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_hnfs_pca_plots","02_hnfs")
plot_pca_figs(so_heart,"02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_heart_pca_plots","02_heart")
plot_pca_figs(so_extremities,"02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_extremities_pca_plots","02_extremities")
plot_pca_figs(so_others,"02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_others_pca_plots","02_others")

## Make a special PCA plot with just the first two clinical groups
plot_pca_by_clin(so,"02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_PCA_by_site_and_hclust.png",partial=TRUE)
plot_pca_by_clin(so,"02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_PCA_by_site_and_hclust.pdf",partial=TRUE)

## FeaturePlot
subset_feature_plot_cut_breast = FeaturePlot(so_cut_breast,c("FLT4","COL22A1","MYC","COL14A1","DSC3","FBN1","KRT1","CLEC2A","FOXN1","IL8","JUN","SCQ"),pt.size = 4)
ggsave(subset_feature_plot_cut_breast,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_subset_feature_plot_cutaneous_breast.png",dpi=300,width=12,height=6)
ggsave(subset_feature_plot_cut_breast,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_subset_feature_plot_cutaneous_breast.pdf",dpi=300,width=12,height=6)





