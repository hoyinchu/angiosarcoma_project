library(Seurat)
library(cowplot)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(stringr)
library(ggpubr)

setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts")

## Load the palettes
source("util_scripts/project_palettes.R")

## Load enriched genes and pathways
all_deseq_res_pos = read.csv("02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_sites_combined_results_pos.tsv", sep = "\t")
deseq2_combined_fora = read.csv("02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_sites_combined_FORA_results.tsv", sep = "\t")
## Load the hierarchical clusters
so = readRDS("data/processed/rna/ASCSeuratObj2025.rds")
hclust_res_pos = read.csv("02_Tumor_RNA_Analysis/outputs/DEGs/02_DESEQ2_hclust_expr_sites_combined_results_pos.tsv", sep = "\t")
hclust_combined_fora = read.csv("02_Tumor_RNA_Analysis/outputs/DEGs/02_expr_hclust_combined_FORA_results.tsv", sep = "\t")


## Select Gene sets to Highlight per site
## Given marker df, fora df, cluster name, and pathway
## Return the overlapped genes sorted by significance
## Maybe also require that the top 5 genes are also significant?

pick_genes = function(marker_df,fora_df,cluster,pathway) {
  marker_subset = marker_df[marker_df$cluster==cluster,]
  fora_df_subset = fora_df[
    (fora_df$cluster==cluster)&(fora_df$pathway==pathway),"overlapGenes"
  ][1]
  gene_list = strsplit(fora_df_subset,",")[[1]]
  marker_df_subset = marker_subset[marker_subset$gene %in% gene_list,]
  marker_df_subset = marker_df_subset[order(marker_df_subset$padj),]
  #marker_df_subset = marker_df_subset[order(-marker_df_subset$log2FoldChange),]
  return(marker_df_subset)
}

## Given marker df and fora df, for a given cluster, pick the top 10 most overexpressed genes
## if the overexpressed gene is in an enriched pathway
pick_genes_any_enriched = function(marker_df,fora_df,cluster,term_set="Hallmarks",top_n=10) {
  marker_subset = marker_df[marker_df$cluster==cluster,]
  marker_subset = marker_subset[order(marker_subset$pad),]
  #marker_subset = head(marker_subset,100)
  ## Reduce computation burden by keeping only top 100
  enriched_pathways = fora_df[fora_df$padj < 0.1,]
  qualifying_genes = enriched_pathways %>% pull(overlapGenes) %>%
    str_split(pattern = ",") %>%
    unlist() %>%
    unique()
  marker_subset$is_qualified = marker_subset$gene %in% qualifying_genes
  top_subset = marker_subset[marker_subset$is_qualified,]
  top_subset = head(top_subset,top_n)
  return(top_subset)
}
pick_genes_any_enriched(all_deseq_res_pos,deseq2_combined_fora,"ParenchymalBreast")

make_gene_set_df = function(marker_df,fora_df,cluster,pathway,name,top_n=10) {
  #picked_genes = pick_genes(all_deseq_res_pos,deseq2_combined_fora,cluster,pathway)[1:top_n,"gene"]
  picked_genes = pick_genes_any_enriched(all_deseq_res_pos,deseq2_combined_fora,cluster)[1:top_n,"gene"]
  picked_df = data.frame(row.names = picked_genes)
  picked_df$`Gene Set` = name
  return(picked_df)
}

## Given marker df and fora df, and cluster name, and a mapping between cluster name to pathway gene sets
## Select overexpressed genes (from high to low by pval in marker df) such that
## the best n genes in the top n pathways are selected
pick_genes_enriched_by_set = function(marker_df,fora_df,cluster,term_set="Hallmarks",top_n_genes=5,top_n_pathway=2) {
  # Step 1: Filtering 
  marker_subset = marker_df[marker_df$cluster==cluster,]
  marker_subset = marker_subset[order(marker_subset$padj),]
  #fora_subset = fora_df[(fora_df$term_set==term_set),]
  #fora_subset = fora_subset[order(fora_subset$padj),]
  top_fora_pathways = fora_df
  #top_fora_pathways = head(fora_df[order(fora_df$padj),],top_n_pathway)
  
  # Step 4: Initialize a list to store the top genes
  all_top_genes = list()
  used_genes = c()
  for (i in 1:nrow(top_fora_pathways)) {
    overlap_genes = top_fora_pathways$overlapGenes[i] %>% 
      str_split(pattern = ",") %>%
      unlist()
    overlap_genes = setdiff(overlap_genes, used_genes) # Exclude genes already included in previous sets
    marker_subset$is_qualified = marker_subset$gene %in% overlap_genes
    top_subset = marker_subset[marker_subset$is_qualified,]
    top_subset_genes = head(top_subset,top_n_genes)$gene
    used_genes = c(used_genes, top_subset_genes)
    # Store the result
    all_top_genes[[i]] = data.frame(
      pathway = top_fora_pathways$pathway[i], 
      gene = top_subset_genes
    )
  }
  # Step 6: Combine all results into one dataframe
  final_df = bind_rows(all_top_genes)
  # Return the combined dataframe
  return(final_df)
}

get_top_fora = function(fora_df,clust,term_set_to_use,top_n) {
  hclust_fora_head = fora_df %>% 
    filter(cluster == clust, term_set == term_set_to_use) %>% 
    arrange(padj) %>% head(top_n)
  return(hclust_fora_head)
}
term_set_to_use = "Hallmarks"
hclust1_fora_head = hclust_combined_fora %>% 
  filter(cluster == "Cluster 1", term_set == term_set_to_use) %>% 
  arrange(padj)# %>% head(3)
hclust2_fora_head = hclust_combined_fora %>% 
  filter(cluster == "Cluster 2", term_set == term_set_to_use) %>% 
  arrange(padj)# %>% head(1)
hclust3_fora_head = hclust_combined_fora %>% 
  filter(cluster == "Cluster 3", term_set == term_set_to_use) %>% 
  arrange(padj)# %>% head(1)
hclust4_fora_head = hclust_combined_fora %>% 
  filter(cluster == "Cluster 4", term_set == term_set_to_use) %>% 
  arrange(padj)# %>% head(3)
hclust5_fora_head = hclust_combined_fora %>% 
  filter(cluster == "Cluster 5", term_set == term_set_to_use) %>% 
  arrange(padj)# %>% head(3)

## Pick 1 from hallmark/C6 depending on which cluster has the most of its members
check_cluster_member_counts = table(so@meta.data$hcluster_by_expr_str,so@meta.data$seurat_clusters_by_site_str)
gene_set_df = rbind(
  pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
                                                            "Cluster 1","Hallmarks",1),"ParenchymalBreast"),
  pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
                                                            "Cluster 1","C6",1),"ParenchymalBreast"),
  pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
                                                            "Cluster 2","Hallmarks",1),"HNFS"),
  pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
                                                            "Cluster 2","C6",1),"HNFS"),
  pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
                                                            "Cluster 3","Hallmarks",1),"HNFS"),
  pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
                                                            "Cluster 3","C6",1),"HNFS"),
  pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
                                                            "Cluster 4","Hallmarks",1),"CutaneousBreast"),
  pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
                                                            "Cluster 4","C6",1),"CutaneousBreast"),
  pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
                                                            "Cluster 5","Hallmarks",1),"CutaneousBreast"),
  pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
                                                            "Cluster 5","C6",1),"CutaneousBreast")#,
  # pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
  #                                                           "Cluster 2","Hallmarks",1),"Heart"),
  # pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
  #                                                           "Cluster 2","C6",1),"Heart"),
  # pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
  #                                                           "Cluster 5","Hallmarks",1),"Extremities"),
  # pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
  #                                                           "Cluster 5","C6",1),"Extremities"),
  # pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
  #                                                           "Cluster 2","Hallmarks",1),"Others"),
  # pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
  #                                                           "Cluster 2","C6",1),"Others")
)

# ## This is the pick 1 from hallmark and 1 from C6 gene set
# gene_set_df = rbind(
#   pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
#                                                             "Cluster 1","Hallmarks",1),"ParenchymalBreast"),
#   pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
#                                                             "Cluster 1","C6",1),"ParenchymalBreast"),
#   pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
#                                                             "Cluster 2","Hallmarks",1),"CutaneousBreast"),
#   pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
#                                                             "Cluster 2","C6",1),"CutaneousBreast"),
#   pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
#                                                             "Cluster 3","Hallmarks",1),"CutaneousBreast"),
#   pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
#                                                             "Cluster 3","C6",1),"CutaneousBreast"),
#   pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
#                                                             "Cluster 4","Hallmarks",1),"HNFS"),
#   # pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
#   #                                                           "Cluster 4","Hallmarks",2),"HNFS")[6:10,], # The most enriched pathway (estrogen) is already taken by cluster 2
#   # so we go for the next highest enriched one
#   pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
#                                                             "Cluster 4","C6",1),"HNFS"),
#   pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
#                                                             "Cluster 5","Hallmarks",1),"HNFS"),
#   pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
#                                                             "Cluster 5","C6",1),"HNFS")
#   # pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
#   #                                                           "Cluster 1","Hallmarks",1),"HNFS"),
#   # pick_genes_enriched_by_set(all_deseq_res_pos,get_top_fora(hclust_combined_fora,
#   #                                                           "Cluster 1","C6",1),"HNFS")
#   #pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust_combined_fora,2),"CutaneousBreast"),
#   #pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust3_fora_head,2),"CutaneousBreast"),
#   #pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust1_fora_head,2),"HNFS"),
#   #pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust4_fora_head,2),"HNFS"),
#   #pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust5_fora_head,2),"HNFS")
#   #pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust1_fora_head,2),"Heart"),
#   #pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust1_fora_head,2),"Extremities"),
#   #pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust2_fora_head,2),"Extremities"),
#   #pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust1_fora_head,2),"Others")
#   #pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust3_fora_head,2),"Others")
# )

# ## This is the pick two from hallmark geneset
# gene_set_df = rbind(
#   pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust1_fora_head,1),"ParenchymalBreast"),
#   pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust2_fora_head,2),"CutaneousBreast"),
#   pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust3_fora_head,2),"CutaneousBreast"),
#   pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust4_fora_head,2),"HNFS"),
#   pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust5_fora_head,2),"HNFS")
#   #pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust1_fora_head,2),"Heart"),
#   #pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust1_fora_head,2),"Extremities"),
#   #pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust2_fora_head,2),"Extremities"),
#   #pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust1_fora_head,2),"Others")
#   #pick_genes_enriched_by_set(all_deseq_res_pos,head(hclust3_fora_head,2),"Others")
# )

gene_set_df = gene_set_df[!duplicated(gene_set_df$gene),]
rownames(gene_set_df) = gene_set_df$gene
gene_set_df$`Gene Set` = gene_set_df$pathway
gene_set_df$`Gene Set` = gsub("\\.V1_UP", "", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("SINGH_KRAS_DEPENDENCY_SIGNATURE", "KRAS_DEP", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("\\HALLMARK_", "", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("EPITHELIAL_MESENCHYMAL_TRANSITION", "EMT", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("ESTROGEN_RESPONSE_LATE", "ESTROGEN", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("KRAS_SIGNALING_DN", "KRAS", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("MYC_TARGETS_V1", "MYC_TARGETS", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("HEME_METABOLISM", "HEME META", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("KRAS.600_UP", "KRAS UP", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("CHOLESTEROL_HOMEOSTASIS", "CHOLESTEROL", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("NFE2L2.V2", "NFE2L2", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("FATTY_ACID_METABOLISM", "FATTY_ACID", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("\\_", " ", gene_set_df$`Gene Set`)

gene_set_df

gene_set_df$pathway = NULL
gene_set_df$gene = NULL

## If using the hallmark / C6 mixed gene set
gene_set_df$`Gene Set` = factor(
  gene_set_df$`Gene Set`,
  levels = c(
    "EMT","TGFB UP",
    "MYC TARGETS","CSR LATE UP",
    "ESTROGEN","KRAS DEP",
    "CHOLESTEROL","NFE2L2",
    "HEME META","KRAS UP"
    #"MYC_TARGETS_V1","E2F_TARGETS",
    #"ESTROGEN_RESPONSE_LATE","KRAS_SIGNALING_DN",
    #"P53_PATHWAY","FATTY_ACID_METABOLISM"
    #"AKT_UP","MEK_UP","SIRNA_EIF4GI_UP"
  )
)

# ## If using C6 gene set
# gene_set_df$`Gene Set` = factor(
#   gene_set_df$`Gene Set`,
#   levels = c(
#     "TGFB_UP","PTEN_DN",
#     "KRAS_DEP","NRL_DN","P53_DN",
#     "CSR_LATE_UP","RB_P107_DN",
#     "AKT_UP","MEK_UP","SIRNA_EIF4GI_UP"
#   )
# )
## If using Hallmark gene set
# gene_set_df$`Gene Set` = factor(
#   gene_set_df$`Gene Set`,
#   levels = c(
#     "EMT","MYOGENESIS",
#     "MYC","E2F_TARGETS",
#     "ESTROGEN_RESPONSE","KRAS",
#     "P53_PATHWAY","METABOLISM"
#     #"MYC_TARGETS_V1","E2F_TARGETS",
#     #"ESTROGEN_RESPONSE_LATE","KRAS_SIGNALING_DN",
#     #"P53_PATHWAY","FATTY_ACID_METABOLISM"
#     #"AKT_UP","MEK_UP","SIRNA_EIF4GI_UP"
#   )
# )


# tissue_pathway_map = list(
#   "ParenchymalBreast" = "TGFB_UP.V1_UP",
#   "CutaneousBreast" = "MYC_UP.V1_UP",
#   "HNFS" = "PGF_UP.V1_UP",
#   "Heart" = "ESC_J1_UP_LATE.V1_UP",
#   "Extremities" = "KRAS.600.LUNG.BREAST_UP.V1_UP",
#   "Others" = "PTEN_DN.V1_UP"
# )
# pathway_name_map = list(
#   "TGFB_UP.V1_UP" = "TGF-Beta Signature",
#   "MYC_UP.V1_UP" = "MYC Signature",
#   "PGF_UP.V1_UP" = "PGF Signature",
#   "ESC_J1_UP_LATE.V1_UP" = "Embryonic Stem Cell Signature",
#   "KRAS.600.LUNG.BREAST_UP.V1_UP" = "KRAS Signature",
#   "PTEN_DN.V1_UP" = "PTEN Mutant Signature"
# )

# Combine all generated dataframes
# gene_set_dfs = lapply(names(tissue_hallmark_map), function(tissue) {
#   make_gene_set_df(all_deseq_res_pos, deseq2_combined_fora, tissue, tissue_hallmark_map[[tissue]], pathway_name_map[[tissue_hallmark_map[[tissue]]]])
# })

# gene_set_dfs = lapply(names(tissue_hallmark_map), function(tissue) {
#   make_gene_set_df(all_deseq_res_pos, deseq2_combined_fora, tissue, tissue_hallmark_map[[tissue]], pathway_name_map[[tissue_hallmark_map[[tissue]]]])
# })

#gene_set_df_old = do.call(rbind, gene_set_dfs)
#gene_set_df_old$`Gene Set` = factor(gene_set_df_old$`Gene Set`,levels = unname(pathway_name_map))

## Define color palette for gene sets
total_sets = length(unique(gene_set_df$`Gene Set`))
dynamic_pal_npg = pal_npg("nrc")(total_sets)
#gene_set_palette = setNames(dynamic_pal_npg, levels(gene_set_df$`Gene Set`))

## Define color palette for the experiemtnal hybrid set
# Get unique pathways
unique_pathways = unique(gene_set_df$`Gene Set`)
npg_palette = pal_npg("nrc")(length(unique_pathways))
gene_set_palette = setNames(npg_palette, unique_pathways)

# 
# for (i in 1:total_sets) {
#   
# }
# gene_set_palette = c(
#   "Epithelial Mesenchymal Transition"=pal_npg("nrc")(6)[1],
#   "MYC Targets"=pal_npg("nrc")(6)[2],
#   "P53 Pathway"=pal_npg("nrc")(6)[3],
#   "TGF-Beta Signaling"=pal_npg("nrc")(6)[4],
#   "Myogenesis"=pal_npg("nrc")(6)[5],
#   "Heme Metabolism"=pal_npg("nrc")(6)[6]
# )

#fora_c6_subset = deseq2_combined_fora[deseq2_combined_fora$term_set=="C6",]
#fora_c6_subset = fora_c6_subset[order(fora_c6_subset$padj),]
#fora_c6_subset_dedupped = fora_c6_subset %>% arrange(padj) %>% distinct(cluster, .keep_all = TRUE)


# ## Pick the top 5 genes by log2FC in a pathway that the site overexpressed in
# noncut_breast_genes = pick_genes(all_deseq_res_pos,deseq2_combined_fora,"ParenchymalBreast","TGFB_UP.V1_UP")[1:5,"gene"]
# 
# #noncut_breast_genes = pick_genes(all_deseq_res_pos,deseq2_combined_fora,"ParenchymalBreast","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")[1:5,"gene"]
# #pick_genes(all_deseq_res_pos,deseq2_combined_fora,"ParenchymalBreast","HALLMARK_MYOGENESIS")
# #pick_genes(all_deseq_res_pos,deseq2_combined_fora,"ParenchymalBreast","HALLMARK_TGF_BETA_SIGNALING")
# #pick_genes(all_deseq_res_pos,deseq2_combined_fora,"ParenchymalBreast","HALLMARK_ANGIOGENESIS")
# noncut_breast_df = data.frame(row.names = noncut_breast_genes)
# #noncut_breast_df$`Gene Set` = "Epithelial Mesenchymal Transition"
# noncut_breast_df$`Gene Set` = "TGF-Beta Upregulated"
# 
# 
# #cut_breast_genes = pick_genes(all_deseq_res_pos,deseq2_combined_fora,"CutaneousBreast","HALLMARK_MYC_TARGETS_V2")[1:5,"gene"]
# cut_breast_genes = pick_genes(all_deseq_res_pos,deseq2_combined_fora,"CutaneousBreast","HALLMARK_MYC_TARGETS_V2")[1:5,"gene"]
# 
# #pick_genes(all_deseq_res_pos,deseq2_combined_fora,"CutaneousBreast","HALLMARK_E2F_TARGETS")
# #pick_genes(all_deseq_res_pos,deseq2_combined_fora,"CutaneousBreast","HALLMARK_UV_RESPONSE_UP")
# cut_breast_df = data.frame(row.names = cut_breast_genes)
# cut_breast_df$`Gene Set` = "MYC Targets"
# 
# hnfs_genes = pick_genes(all_deseq_res_pos,deseq2_combined_fora,"HNFS","HALLMARK_P53_PATHWAY")[1:5,"gene"]
# #pick_genes(all_deseq_res_pos,deseq2_combined_fora,"Heart","HALLMARK_COAGULATION")
# hnfs_df = data.frame(row.names = hnfs_genes)
# hnfs_df$`Gene Set` = "P53 Pathway"
# 
# heart_genes = pick_genes(all_deseq_res_pos,deseq2_combined_fora,"Heart","HALLMARK_TGF_BETA_SIGNALING")[1:5,"gene"]
# heart_df = data.frame(row.names = heart_genes)
# heart_df$`Gene Set` = "TGF-Beta Signaling"
# 
# extremities_genes = pick_genes(all_deseq_res_pos,deseq2_combined_fora,"Extremities","HALLMARK_MYOGENESIS")[1:5,"gene"]
# extremities_df = data.frame(row.names = extremities_genes)
# extremities_df$`Gene Set` = "Myogenesis"
# 
# others_genes = pick_genes(all_deseq_res_pos,deseq2_combined_fora,"Others","HALLMARK_HEME_METABOLISM")[1:5,"gene"]
# others_df = data.frame(row.names = others_genes)
# others_df$`Gene Set` = "Heme Metabolism"
# 
# gene_set_df = rbind(
#   #noncut_breast_df,
#   #cut_breast_df,
#   #hnfs_df,
#   heart_df,
#   extremities_df,
#   others_df
# )
# 
# gene_set_df$`Gene Set` = factor(
#   gene_set_df$`Gene Set`,
#   levels = c("Epithelial Mesenchymal Transition",
#              "MYC Targets",
#              "P53 Pathway",
#              "TGF-Beta Signaling",
#              "Myogenesis",
#              "Heme Metabolism"
#   ))
# 
# gene_set_palette = c(
#   "Epithelial Mesenchymal Transition"=pal_npg("nrc")(6)[1],
#   "MYC Targets"=pal_npg("nrc")(6)[2],
#   "P53 Pathway"=pal_npg("nrc")(6)[3],
#   "TGF-Beta Signaling"=pal_npg("nrc")(6)[4],
#   "Myogenesis"=pal_npg("nrc")(6)[5],
#   "Heme Metabolism"=pal_npg("nrc")(6)[6]
# )
# 
# gene_set_df

## Sort the idendities,sites,and heatmap data
asc_cluster_idents = so@meta.data %>% arrange(seurat_clusters_by_site,`entity:sample_id`) %>% dplyr::select(seurat_clusters_by_site_str)
asc_cluster_idents$seurat_clusters_by_site_str = factor(
  asc_cluster_idents$seurat_clusters_by_site_str,
  levels = c(
    "Breast (Parenchymal)",
    "Breast (Cutaneous)",
    "HNFS",
    "Heart",
    "Extremities",
    "Others"
  )
)
asc_primary_sites =  so@meta.data %>% arrange(seurat_clusters_by_site,`entity:sample_id`) %>% dplyr::select(`Primary Site (Recombined)`)

asc_new_marker_heatmap_data = GetAssayData(so,slot="vst_scaled")[rownames(gene_set_df),rownames(asc_cluster_idents)]

## Get Heatmap Data
asc_gene_set_annotations = rowAnnotation(
  ` `=gene_set_df$`Gene Set`,
  annotation_legend_param = list(` ` = list(title = "Gene Set")),
  col=list(` `=gene_set_palette)
)

## Annotate Fold Change
FC_col_annot = colorRamp2(c(0, 5), c("white", "red"))
#logp_col_annot = colorRamp2(c(0, 5), c("white", "purple"))
purples_cols = brewer.pal(4,"Purples")
significance_pal = c(
  "adj. p < 0.01" = purples_cols[4],
  "adj. p < 0.1" = purples_cols[3],
  "nom. p < 0.05" = purples_cols[2],
  "nom. p >= 0.05" = purples_cols[1]
)
all_markers_dedupped = all_deseq_res_pos %>% arrange(-log2FoldChange) %>% distinct(gene, .keep_all = TRUE)
all_new_markers_FC = all_markers_dedupped %>% filter(gene %in% rownames(gene_set_df))
all_new_markers_FC = all_new_markers_FC %>% mutate(
  significance = case_when(
    padj < 0.01 ~ "adj. p < 0.01",
    padj < 0.1 ~ "adj. p < 0.1",
    pvalue < 0.05 ~ "nom. p < 0.05",
    TRUE ~ "nom. p >= 0.05"
  )
)
all_new_markers_FC

# Annotate the gene's FC in its highest expressed cluster relative to other clusters
asc_FC_annotations = rowAnnotation(
  `Gene log2FC`=all_new_markers_FC[rownames(gene_set_df),"log2FoldChange"],
  `Gene Significance`=all_new_markers_FC[rownames(gene_set_df),"significance"],
  annotation_name_rot=90,
  col=list(`Gene log2FC`=FC_col_annot,`Gene Significance`=significance_pal)
)

# Annotate each sample by their most representative somatic, germline, or both mutations
# This annotation require inputs from running all other scripts so if it is your first time 
# Running this block you should come back after running all other scripts first
rep_mut_df = read.csv("03_Tumor_WES_Analysis/outputs/tables/representative_mutation_table.csv",check.names = FALSE)
rep_mut_df = rep_mut_df[,c("sample_alias_cleaned","tumor_mutation_id","tumor_mutation_id_short","Hugo_Symbol")] %>%
  dplyr::rename(
    rep_tumor_mutation_id = "tumor_mutation_id",
    rep_tumor_mutation_id_short = "tumor_mutation_id_short",
    rep_Hugo_Symbol = "Hugo_Symbol"
  )
representative_somatic_muts_merged = merge(so@meta.data,rep_mut_df,by="sample_alias_cleaned",all.x=TRUE)

## Also add in sample TMB
tmb_df = read.csv("03_Tumor_WES_Analysis/outputs/tables/total_mutations_per_sample.csv")
tmb_df_smol = tmb_df %>%
  dplyr::select(sample_alias,TMB_all_mutations)
representative_somatic_muts_merged_with_tmb = merge(
  representative_somatic_muts_merged,tmb_df_smol,by.x="sample_alias_cleaned",
  by.y="sample_alias",all.x=TRUE
)
# assign tmb categories
representative_somatic_muts_merged_with_tmb = representative_somatic_muts_merged_with_tmb %>%
  mutate(
    tmb_cat = case_when(
      TMB_all_mutations <= 1 ~ "<= 1 mut/MB",
      TMB_all_mutations <= 3 ~ "<= 3 mut/MB",
      TMB_all_mutations <= 10 ~ "<= 10 mut/MB",
      TMB_all_mutations > 10 ~ "> 10 mut/MB",
      TRUE ~ "NA"
    )
  )
# TMB palette
orange_cols = brewer.pal(4,"Oranges")
tmb_palette = c(
  "<= 1 mut/MB" = orange_cols[1],
  "<= 3 mut/MB" = orange_cols[2],
  "<= 10 mut/MB" = orange_cols[3],
  "> 10 mut/MB" = orange_cols[4],
  "NA" = "gray"
)

germline_rep_mut_df = read.csv("05_Germline_WES_Analysis/outputs/tables/representative_germline_pv_table.csv",check.names=FALSE)
germline_cols_to_keep = c("Sample","inferred_ancestry_PCA","inferred_sex","Hugo_Symbol","STUDY ID","germline_gene_priority_score")
germline_rep_mut_df_smol = germline_rep_mut_df[,germline_cols_to_keep]
germline_rep_mut_df_smol = germline_rep_mut_df_smol %>%
  dplyr::rename(
    "germline_inferred_sex" = "inferred_sex",
    "germline_inferred_ancestry_PCA" = "inferred_ancestry_PCA",
    "germline_Hugo_Symbol" = "Hugo_Symbol"
  )
# Making sure only one sample per study id
length(unique(germline_rep_mut_df_smol$`STUDY ID`)) == length(unique(germline_rep_mut_df_smol$`Sample`))
representative_somatic_muts_merged_with_germline = merge(
  representative_somatic_muts_merged_with_tmb,
  germline_rep_mut_df_smol,
  by="STUDY ID",
  all.x=TRUE
)

representative_somatic_muts_merged = representative_somatic_muts_merged_with_germline %>%
  arrange(seurat_clusters_by_site,`entity:sample_id`) %>% 
  replace_na(list(rep_Hugo_Symbol="Not available"))

representative_germline_vars_merged = representative_somatic_muts_merged_with_germline %>%
  arrange(seurat_clusters_by_site,`entity:sample_id`) %>% 
  #select(`germline_Hugo_Symbol`) %>% 
  replace_na(list(germline_Hugo_Symbol="No PV Detected"))

## Only highlight genes that are interesting (somatic)
rep_som_mut_palette = c(
  "POT1"=pal_npg("nrc")(10)[1],
  "TP53"=pal_npg("nrc")(10)[2],
  "CFTR"=pal_npg("nrc")(10)[3],
  "KDR"=pal_npg("nrc")(10)[4],
  "PLCG1"=pal_npg("nrc")(10)[5],
  #"FLG"=pal_npg("nrc")(10)[3],
  #"FLT1"=pal_npg("nrc")(10)[1],
  #"FLT3"=pal_npg("nrc")(10)[3],
  #"FLT4"=pal_npg("nrc")(10)[4],
  #"FLG"=pal_npg("nrc")(10)[8],
  #"BRAF"=pal_npg("nrc")(10)[9],
  "Others"="black",
  "Not available"="gray"
)

rep_germ_var_palette = c(
  "POT1"=pal_npg("nrc")(10)[1],
  "TP53"=pal_npg("nrc")(10)[2],
  "CFTR"=pal_npg("nrc")(10)[3],
  "BRCA2"=pal_npg("nrc")(10)[6],
  "CHEK2"=pal_npg("nrc")(10)[7],
  #"FLG"=pal_npg("nrc")(10)[3],
  #"BRCA1"=pal_npg("nrc")(10)[7],
  #"MUTYH"=pal_npg("nrc")(10)[10],
  #"PKHD1"="red",
  #"USH2A"="pink",
  #"PAH"=pal_npg("nrc")(10)[11],
  #"GJB2"=pal_npg("nrc")(10)[12],
  "Others"="black",
  "No PV Detected"="lightblue"
)

rep_som_genes_to_highlight = names(rep_som_mut_palette)
rep_germ_genes_to_highlight = names(rep_germ_var_palette)

representative_somatic_muts_merged = representative_somatic_muts_merged %>%
  mutate(highlight_som_rep=case_when(
    rep_Hugo_Symbol %in% rep_som_genes_to_highlight ~ rep_Hugo_Symbol,
    TRUE ~ "Others"
  )) %>%
  mutate(`SEX (EHR_EXTRACTED)`=case_when(
    `SEX (EHR_EXTRACTED)`=="F" ~ "Female",
    `SEX (EHR_EXTRACTED)`=="M" ~ "Male",
    TRUE ~ NA
  ))
representative_germline_vars_merged = representative_germline_vars_merged %>%
  mutate(highlight_germ_rep=case_when(
    germline_Hugo_Symbol %in% rep_germ_genes_to_highlight ~ germline_Hugo_Symbol,
    TRUE ~ "Others"
  ))
representative_germline_vars_merged

bx_palette = c(
  "YES" = "red",
  "NO" = "lightpink"
)
representative_somatic_muts_merged = representative_somatic_muts_merged %>%
  mutate(is_epithelioid=ifelse(as.character(BX_EPITHELIOID)=="YES (FOCAL)","YES",as.character(BX_EPITHELIOID)))
representative_somatic_muts_merged = representative_somatic_muts_merged %>%
  mutate(is_spindle=ifelse(as.character(BX_SPINDLE_CELL)=="YES (FOCAL)","YES",as.character(BX_SPINDLE_CELL)))

nuclear_grade_palette = c(
  "HIGH" = "orange",
  "LOW" = "lightpink"
)

representative_somatic_muts_merged = representative_somatic_muts_merged %>%
  mutate(has_mets = case_when(
    `HAS METS AT DX (EHR_EXTRACTED)`==1 ~ "YES",
    `HAS METS AT DX (EHR_EXTRACTED)`==0 ~ "NO",
    TRUE ~ NA
  ))

mets_dx_palette = c(
  "YES" = "darkblue",
  "NO" = "lightblue"
)

representative_somatic_muts_merged = representative_somatic_muts_merged %>% mutate(
  cutaneous_viz = case_when(
    `CUTANEOUS AS (EHR_EXTRACTED)` == 0 ~ "Non-cutaneous AS",
    `CUTANEOUS AS (EHR_EXTRACTED)` == 1 ~ "Cutaneous AS",
    TRUE ~ NA
  )
)

representative_somatic_muts_merged

## Check if there is a purity bias by primary site
ggplot(representative_somatic_muts_merged, aes(x = hcluster_by_expr_str, y = ESTIMATE_purity)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_compare_means(method = "kruskal.test") 

## Top Annotation
age_col_annot = colorRamp2(c(20, 80), c("white", "purple"))
purity_col_annot = colorRamp2(c(0.4,1), c("white", "pink"))

asc_top_annotations = HeatmapAnnotation(
  `Primary Site`= asc_primary_sites$`Primary Site (Recombined)`,
  `Cutaneous` = representative_somatic_muts_merged$`cutaneous_viz`,
  "RAAS/LAAS" = representative_somatic_muts_merged$RAAS_LAAS_Class,
  `Age`= representative_somatic_muts_merged$`Age (Combined)`,
  "Sex" = representative_somatic_muts_merged$`SEX (EHR_EXTRACTED)`,
  "Sample Purity" = representative_somatic_muts_merged$ESTIMATE_purity,
  "Epithelioid" = representative_somatic_muts_merged$is_epithelioid,
  "Spindle Cell" = representative_somatic_muts_merged$is_spindle,
  "Nuclear Grade" = representative_somatic_muts_merged$BX_NUCLEAR_GRADE,
  "Expression Cluster" = representative_somatic_muts_merged$hcluster_by_expr_str,
  #"Mets at Dx" = representative_somatic_muts_merged$has_mets,
  #"Vasoformative" = representative_somatic_muts_merged$BX_VASOFORMATIVE,
  col=list(
    `Primary Site`=primary_site_palette,
    `Cutaneous` = cutaneous_palette,
    `Expression Cluster`=seurart_cluster_by_expression_palette,
    "RAAS/LAAS" = RAAS_class_palette,
    "Sex" = sex_clin_palette,
    "Age" = age_col_annot,
    "Sample Purity" = purity_col_annot,
    "Epithelioid" = bx_palette,
    "Spindle Cell"=bx_palette,
    "Nuclear Grade" = nuclear_grade_palette
    #"Mets at Dx" = mets_dx_palette
  )
)

## Bottom Annotation
asc_rep_mut_annnotation = HeatmapAnnotation(
  `Repr. Somatic Mut.` = representative_somatic_muts_merged$highlight_som_rep,
  `Repr. Germline Var.` = representative_germline_vars_merged$highlight_germ_rep,
  `TMB` = representative_germline_vars_merged$tmb_cat,
  col = list(
    `Repr. Somatic Mut.`=rep_som_mut_palette,
    `Repr. Germline Var.`=rep_germ_var_palette,
    `TMB` = tmb_palette
  )
)

## Create a column dendrogram based on highlighted genes
so_small_embeddings= Embeddings(so,reduction = "pca")[,1:2]
#so_small_embedding = FetchData(so,vars=rownames(gene_set_df),layer = "vst_scaled")
so_small_hclust = hclust(dist(so_small_embeddings))
so_small_dendro = as.dendrogram(so_small_hclust)

# Make the actual heatmap plot
complex_heatmap_fig = ComplexHeatmap::Heatmap(
  asc_new_marker_heatmap_data,
  #row_order = rownames(gene_set_df),
  cluster_row_slices = FALSE,
  cluster_rows=FALSE,
  cluster_columns=TRUE,
  #cluster_column_slices = T,
  #cluster_columns = so_small_dendro,
  #show_row_names=FALSE,
  show_column_names=FALSE,
  top_annotation=asc_top_annotations,
  bottom_annotation = asc_rep_mut_annnotation,
  right_annotation=asc_gene_set_annotations,
  left_annotation=asc_FC_annotations,
  column_split = asc_cluster_idents,
  cluster_column_slices = FALSE,
  show_column_dend = FALSE,
  #column_dend_reorder = FALSE,
  column_dend_reorder = TRUE,
  #column_split = FALSE,
  row_split = gene_set_df$`Gene Set`,
  heatmap_legend_param = list(
    title = "Scaled Expression", at = c(-3, 0, 3)
    #labels = c("neg_two", "zero", "pos_two")
  )
)
#pdf("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_cluster_gene_set_heatmap_dendro.pdf",width=24,height=18)
pdf("02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_cluster_gene_set_heatmap_dendro_cluster_by_site_top_5_pathway.pdf",width=24,height=18)
draw(complex_heatmap_fig, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()

## After drawing the heatmap calculate some basic summary statistics about each cluster / primary sites

#so = readRDS("data/processed/rna/ASCSeuratObj2025.rds")

## Number of parenchymal breast samples by site
site_by_cluster_table = table(so@meta.data$seurat_clusters_by_site_str,so@meta.data$hcluster_by_expr_str)
site_totals = rowSums(site_by_cluster_table)
cluster_totals = colSums(site_by_cluster_table)

## Number of KDR mutations per cluster
## Symbol level
table(representative_germline_vars_merged$rep_Hugo_Symbol,representative_germline_vars_merged$hcluster_by_expr_str)
## Mutation level
table(representative_germline_vars_merged$rep_tumor_mutation_id_short,representative_germline_vars_merged$hcluster_by_expr_str)

## Check the cluster assignments for cutaneous breast angiosarcomas
site_by_cluster_table
## Check MYC's enrichment in cluster 4
hclust_res_pos[hclust_res_pos$gene=="MYC",]
## Check MYC pathway enrichment in cluster 4
head(hclust_combined_fora[(hclust_combined_fora$cluster=="Cluster 4")&(hclust_combined_fora$term_set=="Hallmarks"),])
## Check CSR LAte up p-values
head(hclust_combined_fora[(hclust_combined_fora$cluster=="Cluster 4")&(hclust_combined_fora$term_set=="C6"),])
## Check expression of NRP2 and FLT4 in cluster 4
hclust_res_pos[hclust_res_pos$gene=="NRP2",]
hclust_res_pos[hclust_res_pos$gene=="FLT4",]

## Check FORA results for cluster 5
head(hclust_combined_fora[(hclust_combined_fora$cluster=="Cluster 5")&(hclust_combined_fora$term_set=="Hallmarks"),])
head(hclust_combined_fora[(hclust_combined_fora$cluster=="Cluster 5")&(hclust_combined_fora$term_set=="C6"),])
## Check expression of IL37 and CTLA4
deseq2_combined_fora[deseq2_combined_fora$gene=="IL37",]
all_deseq_res_pos[all_deseq_res_pos$gene=="CTLA4",]

## Check the cluster assignments for HNFS angiosarcomas
site_by_cluster_table
## Check FORA results for cluster 4
head(hclust_combined_fora[(hclust_combined_fora$cluster=="Cluster 2")&(hclust_combined_fora$term_set=="Hallmarks"),])
head(hclust_combined_fora[(hclust_combined_fora$cluster=="Cluster 3")&(hclust_combined_fora$term_set=="Hallmarks"),])
## Check expression of SRBF2
all_deseq_res_pos[all_deseq_res_pos$gene=="SREBF2",]
hclust_res_pos[hclust_res_pos$gene=="SREBF2",]


