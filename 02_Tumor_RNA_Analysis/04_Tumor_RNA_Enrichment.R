library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(EnhancedVolcano)
library(msigdbr)
library(fgsea)
library(forcats)
library(DESeq2)
library(cowplot)


## Load the previously created seurat data
so = readRDS("../data/processed/rna/ASCSeuratObj2025.rds")

## Load the palettes
source("../util_scripts/project_palettes.R")

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
raw_counts = as.matrix(GetAssayData(so, layer = "counts"))
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
  write.csv(deseq2_results,file = outpath,row.names = FALSE)
  return(deseq2_results)
}

## Make the tables by site
noncut_breast_deseq_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                                 "~ ESTIMATE_purity + rna_inferred_female + site.breast.noncut",c("site.breast.noncut","ParenchymalBreast","Rest"),
                                                 "./outputs/DEGs/02_DESEQ2_noncut_breast_vs_rest.csv")
cut_breast_deseq_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                              "~ ESTIMATE_purity + rna_inferred_female + site.breast.cut",c("site.breast.cut","CutaneousBreast","Rest"),
                                              "./outputs/DEGs/02_DESEQ2_cut_breast_vs_rest.csv")
hnfs_deseq_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                        "~ ESTIMATE_purity + rna_inferred_female + site.hnfs",c("site.hnfs","HNFS","Rest"),
                                        "./outputs/DEGs/02_DESEQ2_hnfs_vs_rest.csv")
heart_deseq_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                         "~ ESTIMATE_purity + rna_inferred_female + site.heart",c("site.heart","Heart","Rest"),
                                         "./outputs/DEGs/02_DESEQ2_heart_vs_rest.csv")
extremities_deseq_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                               "~ ESTIMATE_purity + rna_inferred_female + site.extremities",c("site.extremities","Extremities","Rest"),
                                               "./outputs/DEGs/02_DESEQ2_extremities_vs_rest.csv")
others_deseq_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                          "~ ESTIMATE_purity + rna_inferred_female + site.others",c("site.others","Others","Rest"),
                                          "./outputs/DEGs/02_DESEQ2_others_vs_rest.csv")

## All results regardless of whether they are positive
noncut_breast_deseq_res_all = noncut_breast_deseq_res %>% as.data.frame() %>% mutate(cluster="ParenchymalBreast")
cut_breast_deseq_res_all = cut_breast_deseq_res %>% as.data.frame() %>% mutate(cluster="CutaneousBreast")
hnfs_deseq_res_all = hnfs_deseq_res %>% as.data.frame() %>% mutate(cluster="HNFS")
heart_deseq_res_all = heart_deseq_res %>% as.data.frame() %>% mutate(cluster="Heart")
extremities_deseq_res_all = extremities_deseq_res %>% as.data.frame() %>% mutate(cluster="Extremities")
others_deseq_res_all = others_deseq_res %>% as.data.frame() %>% mutate(cluster="Others")
all_deseq_res_all = rbind(noncut_breast_deseq_res_all,cut_breast_deseq_res_all,hnfs_deseq_res_all,
                          heart_deseq_res_all,extremities_deseq_res_all,others_deseq_res_all)
write.table(all_deseq_res_all,file = "./outputs/DEGs/02_DESEQ2_sites_combined_results_all.csv", sep = "\t", quote = FALSE)


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
write.csv(all_deseq_res_pos,file = "./outputs/DEGs/02_DESEQ2_sites_combined_results_pos.csv", row.names = FALSE)

## Non-cutaneous breast
all_deseq_res_pos_filtered = all_deseq_res_pos %>% filter(baseMean > 100)
all_deseq_res_pos_filtered = all_deseq_res_pos_filtered %>% mutate(padj_sig = padj < 0.05)
noncut_breast_deseq_pos_filtered = all_deseq_res_pos_filtered %>% filter(cluster=="ParenchymalBreast")

#noncut_breast_deseq_pos_filtered_bonferroni_thresh = 0.05/(dim(noncut_breast_deseq_pos_filtered)[1])
noncut_breast_deseq_pos_filtered_bonferroni_thresh = 0.05/(dim(all_deseq_res_pos %>% filter(cluster=="ParenchymalBreast"))[1])
noncut_breast_deseq_pos_filtered$is_bonferroni_sig = noncut_breast_deseq_pos_filtered$pvalue < noncut_breast_deseq_pos_filtered_bonferroni_thresh

noncut_breast_volcano = ggplot(noncut_breast_deseq_pos_filtered,aes(x=log2FoldChange,y=-log10(padj),color=padj_sig)) +
  geom_point() +
  geom_text_repel(
    data = subset(all_deseq_res_pos_filtered_breast, (log2FoldChange > 3.5)|(-log10(padj)>12)),
    aes(label = gene),size=2,color="black",box.padding = 0.2
  ) + pretty_plot() + L_border() + geom_hline(yintercept = -log10(0.05),linetype="dashed") +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) + 
  theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
noncut_breast_volcano

cowplot::ggsave2("./outputs/plots/04_enrichment_plots/noncut_breast_volcano.png",noncut_breast_volcano,dpi=600,height=3,width=3.8)

# cutaneous breast
cut_breast_deseq_pos_filtered = all_deseq_res_pos_filtered %>% filter(cluster=="CutaneousBreast")

#cut_breast_deseq_pos_filtered_bonferroni_thresh = 0.05/(dim(cut_breast_deseq_pos_filtered)[1])
cut_breast_deseq_pos_filtered_bonferroni_thresh = 0.05/(dim(all_deseq_res_pos %>% filter(cluster=="CutaneousBreast"))[1])
cut_breast_deseq_pos_filtered$is_bonferroni_sig = cut_breast_deseq_pos_filtered$pvalue < cut_breast_deseq_pos_filtered_bonferroni_thresh

cut_breast_volcano = ggplot(cut_breast_deseq_pos_filtered,aes(x=log2FoldChange,y=-log10(padj),color=padj_sig)) +
  geom_point() +
  geom_text_repel(
    data = subset(cut_breast_deseq_pos_filtered, (log2FoldChange > 4)|(-log10(padj)>5)|(gene=="FGFR4")),
    aes(label = gene),size=2,max.overlaps = Inf,color="black",box.padding = 0.2
  )  + pretty_plot() + L_border() + geom_hline(yintercept = -log10(0.05),linetype="dashed") +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) + 
  theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
cut_breast_volcano
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/cut_breast_volcano.png",cut_breast_volcano,dpi=600,height=3,width=3.8)


# HNFS
hnfs_deseq_pos_filtered = all_deseq_res_pos_filtered %>% filter(cluster=="HNFS")

#hnfs_deseq_pos_filtered_bonferroni_thresh = 0.05/(dim(hnfs_deseq_pos_filtered)[1])
hnfs_deseq_pos_filtered_bonferroni_thresh = 0.05/(dim(all_deseq_res_pos %>% filter(cluster=="HNFS"))[1])
hnfs_deseq_pos_filtered$is_bonferroni_sig = hnfs_deseq_pos_filtered$pvalue < hnfs_deseq_pos_filtered_bonferroni_thresh

hnfs_volcano = ggplot(hnfs_deseq_pos_filtered,aes(x=log2FoldChange,y=-log10(padj),color=padj_sig)) +
  geom_point() +
  geom_text_repel(
    data = subset(hnfs_deseq_pos_filtered, (-log10(padj)>9)|(log2FoldChange>7.5)),
    aes(label = gene),size=2,max.overlaps = Inf,color="black",box.padding = 0.2
  )  + pretty_plot() + L_border() + geom_hline(yintercept = -log10(0.05),linetype="dashed") +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) + 
  theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
hnfs_volcano
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hnfs_volcano.png",hnfs_volcano,dpi=600,height=3,width=3.8)

# Heart
heart_deseq_pos_filtered = all_deseq_res_pos_filtered %>% filter(cluster=="Heart")

#heart_deseq_pos_filtered_bonferroni_thresh = 0.05/(dim(heart_deseq_pos_filtered)[1])
heart_deseq_pos_filtered_bonferroni_thresh = 0.05/(dim(all_deseq_res_pos %>% filter(cluster=="Heart"))[1])
heart_deseq_pos_filtered$is_bonferroni_sig = heart_deseq_pos_filtered$pvalue < heart_deseq_pos_filtered_bonferroni_thresh

heart_volcano = ggplot(heart_deseq_pos_filtered,aes(x=log2FoldChange,y=-log10(padj),color=padj_sig)) +
  geom_point() +
  geom_text_repel(
    data = subset(heart_deseq_pos_filtered, ((log2FoldChange > 5)&(padj_sig))|(-log10(padj)>6.4)),
    aes(label = gene),size=2,max.overlaps = Inf,color="black",box.padding = 0.2
  )+ pretty_plot() + L_border() + geom_hline(yintercept = -log10(0.05),linetype="dashed") +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) + 
  theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
heart_volcano
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/heart_volcano.png",heart_volcano,dpi=600,height=3,width=3.8)

# Extremities
extremities_deseq_pos_filtered = all_deseq_res_pos_filtered %>% filter(cluster=="Extremities")

extremities_deseq_pos_filtered_bonferroni_thresh = 0.05/(dim(extremities_deseq_pos_filtered)[1])
extremities_deseq_pos_filtered$is_bonferroni_sig = extremities_deseq_pos_filtered$pvalue < extremities_deseq_pos_filtered_bonferroni_thresh

extremities_volcano = ggplot(extremities_deseq_pos_filtered,aes(x=log2FoldChange,y=-log10(padj),color=padj_sig)) +
  geom_point() +
  geom_text_repel(
    data = subset(extremities_deseq_pos_filtered, padj_sig & (log2FoldChange > 4)|(-log10(padj)>5)),
    aes(label = gene),size=2,max.overlaps = Inf,color="black",box.padding = 0.2
  ) + pretty_plot() + L_border() + geom_hline(yintercept = -log10(0.05),linetype="dashed") +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) + 
  theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
extremities_volcano
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/extremities_volcano.png",extremities_volcano,dpi=600,height=3,width=3.8)

# Others
others_deseq_pos_filtered = all_deseq_res_pos_filtered %>% filter(cluster=="Others")

others_deseq_pos_filtered_bonferroni_thresh = 0.05/(dim(others_deseq_pos_filtered)[1])
others_deseq_pos_filtered$is_bonferroni_sig = others_deseq_pos_filtered$pvalue < others_deseq_pos_filtered_bonferroni_thresh

others_volcano = ggplot(others_deseq_pos_filtered,aes(x=log2FoldChange,y=-log10(padj),color=padj_sig)) +
  geom_point() +
  geom_hline(-log10(others_deseq_pos_filtered_bonferroni_thresh)) +
  geom_text_repel(
    data = subset(others_deseq_pos_filtered, (log2FoldChange>3.6)&padj_sig),
    aes(label = gene),size=2,max.overlaps = Inf,color="black",box.padding = 0.2
  ) + pretty_plot() + L_border() + geom_hline(yintercept = -log10(0.05),linetype="dashed") +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) + 
  theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
others_volcano
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/others_volcano.png",others_volcano,dpi=600,height=3,width=3.8)

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
noncut_breast_pos_degs = noncut_breast_deseq_pos_filtered %>% filter(padj < 0.05) %>% pull(gene)

noncut_breast_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,noncut_breast_pos_degs,fora_gene_universe) %>% as.data.frame()
noncut_breast_pos_degs_fora_hallmark$padj_sig = noncut_breast_pos_degs_fora_hallmark$padj < 0.05
noncut_breast_pos_degs_fora_hallmark = add_enriched_or(noncut_breast_pos_degs_fora_hallmark, n_degs = length(noncut_breast_pos_degs), universe_size = length(fora_gene_universe))
noncut_breast_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",noncut_breast_pos_degs_fora_hallmark$pathway)
noncut_breast_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",noncut_breast_pos_degs_fora_hallmark$pretty_pathway)

noncut_breast_pos_degs_fora_c6_up = fora(fgsea_c6_set_up_only,noncut_breast_pos_degs,fora_gene_universe) %>% as.data.frame()
noncut_breast_pos_degs_fora_c6_up$padj_sig = noncut_breast_pos_degs_fora_c6_up$padj < 0.05
noncut_breast_pos_degs_fora_c6_up = add_enriched_or(noncut_breast_pos_degs_fora_c6_up, n_degs = length(noncut_breast_pos_degs), universe_size = length(fora_gene_universe))
noncut_breast_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_UP","",noncut_breast_pos_degs_fora_c6_up$pathway)
noncut_breast_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_"," ",noncut_breast_pos_degs_fora_c6_up$pretty_pathway)

noncut_breast_pos_degs_fora_hallmark_plot = ggplot(noncut_breast_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(noncut_breast_pos_degs_fora_hallmark, padj_sig),
    aes(label = pretty_pathway), size=2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
noncut_breast_pos_degs_fora_hallmark_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/noncut_breast_pos_degs_fora_hallmark_plot.pdf",noncut_breast_pos_degs_fora_hallmark_plot,dpi=300,height=1.4,width=2)
#save_deg_csv(noncut_breast_pos_degs_fora_hallmark,"./outputs/DEGs/noncut_breast_pos_degs_fora_hallmark.csv")

noncut_breast_pos_degs_fora_c6_up_plot = ggplot(noncut_breast_pos_degs_fora_c6_up,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(noncut_breast_pos_degs_fora_c6_up, padj_sig | odds_ratio > 2.5),
    aes(label = pretty_pathway),max.overlaps = Inf, size=2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
noncut_breast_pos_degs_fora_c6_up_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/noncut_breast_pos_degs_fora_c6_up_plot.pdf",noncut_breast_pos_degs_fora_c6_up_plot,dpi=300,height=1.4,width=2)
#save_deg_csv(noncut_breast_pos_degs_fora_c6_up,"./outputs/DEGs/noncut_breast_pos_degs_fora_c6_up.csv")


## Cutaneous breast
cut_breast_pos_degs = cut_breast_deseq_pos_filtered %>% filter(padj < 0.05) %>% pull(gene)

cut_breast_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,cut_breast_pos_degs,fora_gene_universe) %>% as.data.frame()
cut_breast_pos_degs_fora_hallmark$padj_sig = cut_breast_pos_degs_fora_hallmark$padj < 0.05
cut_breast_pos_degs_fora_hallmark = add_enriched_or(cut_breast_pos_degs_fora_hallmark, n_degs = length(cut_breast_pos_degs), universe_size = length(fora_gene_universe))
cut_breast_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",cut_breast_pos_degs_fora_hallmark$pathway)
cut_breast_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",cut_breast_pos_degs_fora_hallmark$pretty_pathway)

cut_breast_pos_degs_fora_c6_up = fora(fgsea_c6_set_up_only,cut_breast_pos_degs,fora_gene_universe) %>% as.data.frame()
cut_breast_pos_degs_fora_c6_up$padj_sig = cut_breast_pos_degs_fora_c6_up$padj < 0.05
cut_breast_pos_degs_fora_c6_up = add_enriched_or(cut_breast_pos_degs_fora_c6_up, n_degs = length(cut_breast_pos_degs), universe_size = length(fora_gene_universe))
cut_breast_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_UP","",cut_breast_pos_degs_fora_c6_up$pathway)
cut_breast_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_"," ",cut_breast_pos_degs_fora_c6_up$pretty_pathway)


cut_breast_pos_degs_fora_hallmark_plot = ggplot(cut_breast_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(cut_breast_pos_degs_fora_hallmark, padj_sig),
    aes(label = pretty_pathway), size=2, max.overlaps = Inf, box.padding = 0.2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
cut_breast_pos_degs_fora_hallmark_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/cut_breast_pos_degs_fora_hallmark_plot.pdf",cut_breast_pos_degs_fora_hallmark_plot,dpi=300,height=1.4,width=2)


cut_breast_pos_degs_fora_c6_up_plot = ggplot(cut_breast_pos_degs_fora_c6_up,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(cut_breast_pos_degs_fora_c6_up, padj_sig),
    aes(label = pretty_pathway), size=2, max.overlaps = Inf,box.padding = 0.2
  )+ theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
cut_breast_pos_degs_fora_c6_up_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/cut_breast_pos_degs_fora_c6_up_plot.pdf",cut_breast_pos_degs_fora_c6_up_plot,dpi=300,height=1.4,width=2)

## HNFS
hnfs_pos_degs = hnfs_deseq_pos_filtered %>% filter(padj < 0.05) %>% pull(gene)

hnfs_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,hnfs_pos_degs,fora_gene_universe) %>% as.data.frame()
hnfs_pos_degs_fora_hallmark$padj_sig = hnfs_pos_degs_fora_hallmark$padj < 0.05
hnfs_pos_degs_fora_hallmark = add_enriched_or(hnfs_pos_degs_fora_hallmark, n_degs = length(hnfs_pos_degs), universe_size = length(fora_gene_universe))
hnfs_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",hnfs_pos_degs_fora_hallmark$pathway)
hnfs_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",hnfs_pos_degs_fora_hallmark$pretty_pathway)

hnfs_pos_degs_fora_c6_up = fora(fgsea_c6_set_up_only,hnfs_pos_degs,fora_gene_universe) %>% as.data.frame()
hnfs_pos_degs_fora_c6_up$padj_sig = hnfs_pos_degs_fora_c6_up$padj < 0.05
hnfs_pos_degs_fora_c6_up = add_enriched_or(hnfs_pos_degs_fora_c6_up, n_degs = length(hnfs_pos_degs), universe_size = length(fora_gene_universe))
hnfs_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_UP","",hnfs_pos_degs_fora_c6_up$pathway)
hnfs_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_"," ",hnfs_pos_degs_fora_c6_up$pretty_pathway)

hnfs_pos_degs_fora_hallmark_plot = ggplot(hnfs_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(hnfs_pos_degs_fora_hallmark, padj_sig & odds_ratio>3),
    aes(label = pretty_pathway),size=2, max.overlaps = Inf,box.padding = 0.2
  )+ theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
hnfs_pos_degs_fora_hallmark_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hnfs_pos_degs_fora_hallmark_plot.pdf",hnfs_pos_degs_fora_hallmark_plot,dpi=300,height=1.4,width=2)

hnfs_pos_degs_fora_c6_up_plot = ggplot(hnfs_pos_degs_fora_c6_up,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(hnfs_pos_degs_fora_c6_up, -log10(padj)>2.75 | (odds_ratio>3)&(-log10(padj)>2)),
    aes(label = pretty_pathway),size=2, max.overlaps = Inf,box.padding = 0.2
  )+ theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
hnfs_pos_degs_fora_c6_up_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hnfs_pos_degs_fora_c6_up_plot.pdf",hnfs_pos_degs_fora_c6_up_plot,dpi=300,height=1.4,width=2)


## Heart
heart_pos_degs = heart_deseq_pos_filtered %>% filter(padj < 0.05) %>% pull(gene)

heart_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,heart_pos_degs,fora_gene_universe) %>% as.data.frame()
heart_pos_degs_fora_hallmark$padj_sig = heart_pos_degs_fora_hallmark$padj < 0.05
heart_pos_degs_fora_hallmark = add_enriched_or(heart_pos_degs_fora_hallmark, n_degs = length(heart_pos_degs), universe_size = length(fora_gene_universe))
heart_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",heart_pos_degs_fora_hallmark$pathway)
heart_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",heart_pos_degs_fora_hallmark$pretty_pathway)

heart_pos_degs_fora_c6_up = fora(fgsea_c6_set_up_only,heart_pos_degs,fora_gene_universe) %>% as.data.frame()
heart_pos_degs_fora_c6_up$padj_sig = heart_pos_degs_fora_c6_up$padj < 0.05
heart_pos_degs_fora_c6_up = add_enriched_or(heart_pos_degs_fora_c6_up, n_degs = length(heart_pos_degs), universe_size = length(fora_gene_universe))
heart_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_UP","",heart_pos_degs_fora_c6_up$pathway)
heart_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_"," ",heart_pos_degs_fora_c6_up$pretty_pathway)

heart_pos_degs_fora_hallmark_plot = ggplot(heart_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(heart_pos_degs_fora_hallmark, padj_sig | (odds_ratio > 2)&(pval<0.05)),
    aes(label = pretty_pathway),size=2, max.overlaps = Inf,box.padding = 0.2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
heart_pos_degs_fora_hallmark_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/heart_pos_degs_fora_hallmark_plot.pdf",heart_pos_degs_fora_hallmark_plot,dpi=300,height=1.4,width=2)


heart_pos_degs_fora_c6_up_plot = ggplot(heart_pos_degs_fora_c6_up,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(heart_pos_degs_fora_c6_up, padj_sig | odds_ratio > 2.5),
    aes(label = pretty_pathway),size=2, max.overlaps = Inf,box.padding = 0.2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
heart_pos_degs_fora_c6_up_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/heart_pos_degs_fora_c6_up_plot.pdf",heart_pos_degs_fora_c6_up_plot,dpi=300,height=1.4,width=2)

## Extremities
extremities_pos_degs = extremities_deseq_pos_filtered %>% filter(padj < 0.05) %>% pull(gene)

extremities_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,extremities_pos_degs,fora_gene_universe) %>% as.data.frame()
extremities_pos_degs_fora_hallmark$padj_sig = extremities_pos_degs_fora_hallmark$padj < 0.05
extremities_pos_degs_fora_hallmark = add_enriched_or(extremities_pos_degs_fora_hallmark, n_degs = length(extremities_pos_degs), universe_size = length(fora_gene_universe))
extremities_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",extremities_pos_degs_fora_hallmark$pathway)
extremities_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",extremities_pos_degs_fora_hallmark$pretty_pathway)

extremities_pos_degs_fora_c6_up = fora(fgsea_c6_set_up_only,extremities_pos_degs,fora_gene_universe) %>% as.data.frame()
extremities_pos_degs_fora_c6_up$padj_sig = extremities_pos_degs_fora_c6_up$padj < 0.05
extremities_pos_degs_fora_c6_up = add_enriched_or(extremities_pos_degs_fora_c6_up, n_degs = length(extremities_pos_degs), universe_size = length(fora_gene_universe))
extremities_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_UP","",extremities_pos_degs_fora_c6_up$pathway)
extremities_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_"," ",extremities_pos_degs_fora_c6_up$pretty_pathway)

extremities_pos_degs_fora_hallmark_plot = ggplot(extremities_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(extremities_pos_degs_fora_hallmark, padj_sig),
    aes(label = pretty_pathway),size=2, max.overlaps = Inf,box.padding = 0.2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
extremities_pos_degs_fora_hallmark_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/extremities_pos_degs_fora_hallmark_plot.pdf",extremities_pos_degs_fora_hallmark_plot,dpi=300,height=1.4,width=2)

extremities_pos_degs_fora_c6_up_plot = ggplot(extremities_pos_degs_fora_c6_up,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(extremities_pos_degs_fora_c6_up, padj_sig | odds_ratio > 3.35),
    aes(label = pretty_pathway),size=2, max.overlaps = Inf,box.padding = 0.2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
extremities_pos_degs_fora_c6_up_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/extremities_pos_degs_fora_c6_up_plot.pdf",extremities_pos_degs_fora_c6_up_plot,dpi=300,height=1.4,width=2)

## Others
others_pos_degs = others_deseq_pos_filtered %>% filter(padj < 0.05) %>% pull(gene)

others_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,others_pos_degs,fora_gene_universe) %>% as.data.frame()
others_pos_degs_fora_hallmark$padj_sig = others_pos_degs_fora_hallmark$padj < 0.05
others_pos_degs_fora_hallmark = add_enriched_or(others_pos_degs_fora_hallmark, n_degs = length(others_pos_degs), universe_size = length(fora_gene_universe))
others_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",others_pos_degs_fora_hallmark$pathway)
others_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",others_pos_degs_fora_hallmark$pretty_pathway)

others_pos_degs_fora_c6_up = fora(fgsea_c6_set_up_only,others_pos_degs,fora_gene_universe) %>% as.data.frame()
others_pos_degs_fora_c6_up$padj_sig = others_pos_degs_fora_c6_up$padj < 0.05
others_pos_degs_fora_c6_up = add_enriched_or(others_pos_degs_fora_c6_up, n_degs = length(others_pos_degs), universe_size = length(fora_gene_universe))
others_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_UP","",others_pos_degs_fora_c6_up$pathway)
others_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_"," ",others_pos_degs_fora_c6_up$pretty_pathway)

others_pos_degs_fora_hallmark_plot = ggplot(others_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(others_pos_degs_fora_hallmark, padj_sig),
    aes(label = pretty_pathway),size=2, max.overlaps = Inf,box.padding = 0.2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
others_pos_degs_fora_hallmark_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/others_pos_degs_fora_hallmark_plot.pdf",others_pos_degs_fora_hallmark_plot,dpi=300,height=1.4,width=2)

others_pos_degs_fora_c6_up_plot = ggplot(others_pos_degs_fora_c6_up,aes(x=odds_ratio,y=-log10(pval),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(others_pos_degs_fora_c6_up, padj_sig | odds_ratio > 5),
    aes(label = pretty_pathway),size=2, max.overlaps = Inf,box.padding = 0.2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
others_pos_degs_fora_c6_up_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/others_pos_degs_fora_c6_up_plot.pdf",others_pos_degs_fora_c6_up_plot,dpi=300,height=1.4,width=2)

## Combine all DEGs and write it out
noncut_breast_pos_degs_fora_hallmark$cluster = "ParenchymalBreast"
noncut_breast_pos_degs_fora_hallmark$term_set = "Hallmarks"
noncut_breast_pos_degs_fora_c6_up$cluster = "ParenchymalBreast"
noncut_breast_pos_degs_fora_c6_up$term_set = "C6"

cut_breast_pos_degs_fora_hallmark$cluster = "CutaneousBreast"
cut_breast_pos_degs_fora_hallmark$term_set = "Hallmarks"
cut_breast_pos_degs_fora_c6_up$cluster = "CutaneousBreast"
cut_breast_pos_degs_fora_c6_up$term_set = "C6"

hnfs_pos_degs_fora_hallmark$cluster = "HNFS"
hnfs_pos_degs_fora_hallmark$term_set = "Hallmarks"
hnfs_pos_degs_fora_c6_up$cluster = "HNFS"
hnfs_pos_degs_fora_c6_up$term_set = "C6"

heart_pos_degs_fora_hallmark$cluster = "Heart"
heart_pos_degs_fora_hallmark$term_set = "Hallmarks"
heart_pos_degs_fora_c6_up$cluster = "Heart"
heart_pos_degs_fora_c6_up$term_set = "C6"

extremities_pos_degs_fora_hallmark$cluster = "Extremities"
extremities_pos_degs_fora_hallmark$term_set = "Hallmarks"
extremities_pos_degs_fora_c6_up$cluster = "Extremities"
extremities_pos_degs_fora_c6_up$term_set = "C6"

others_pos_degs_fora_hallmark$cluster = "Others"
others_pos_degs_fora_hallmark$term_set = "Hallmarks"
others_pos_degs_fora_c6_up$cluster = "Others"
others_pos_degs_fora_c6_up$term_set = "C6"

deg_hallmark_c6_combined = rbind(
  noncut_breast_pos_degs_fora_hallmark,cut_breast_pos_degs_fora_hallmark,
  hnfs_pos_degs_fora_hallmark,heart_pos_degs_fora_hallmark,
  extremities_pos_degs_fora_hallmark,others_pos_degs_fora_hallmark,
  noncut_breast_pos_degs_fora_c6_up,cut_breast_pos_degs_fora_c6_up,
  hnfs_pos_degs_fora_c6_up,heart_pos_degs_fora_c6_up,
  extremities_pos_degs_fora_c6_up,others_pos_degs_fora_c6_up
)

save_deg_csv(deg_hallmark_c6_combined,"./outputs/DEGs/deg_hallmark_c6_combined.csv")

## Repeat the analysis above but with HClust derived clusters instead
## Make the tables by site
hclust1_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                     "~ ESTIMATE_purity + rna_inferred_female + expr.hclust1",c("expr.hclust1","Cluster1","Rest"),
                                     "./outputs/DEGs/02_hclust_Cluster1_vs_rest.tsv")
hclust2_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                     "~ ESTIMATE_purity + rna_inferred_female + expr.hclust2",c("expr.hclust2","Cluster2","Rest"),
                                     "./outputs/DEGs/02_hclust_Cluster2_vs_rest.tsv")
hclust3_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                     "~ ESTIMATE_purity + rna_inferred_female + expr.hclust3",c("expr.hclust3","Cluster3","Rest"),
                                     "./outputs/DEGs/02_hclust_Cluster3_vs_rest.tsv")
hclust4_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                     "~ ESTIMATE_purity + rna_inferred_female + expr.hclust4",c("expr.hclust4","Cluster4","Rest"),
                                     "./outputs/DEGs/02_hclust_Cluster4_vs_rest.tsv")
hclust5_res = one_versus_rest_deseq2(raw_counts,deseq2_metadata,
                                     "~ ESTIMATE_purity + rna_inferred_female + expr.hclust5",c("expr.hclust5","Cluster5","Rest"),
                                     "./outputs/DEGs/02_hclust_Cluster5_vs_rest.tsv")

hclust1_res_all = hclust1_res %>% as.data.frame() %>% mutate(cluster="Cluster 1")
hclust2_res_all = hclust2_res %>% as.data.frame() %>% mutate(cluster="Cluster 2")
hclust3_res_all = hclust3_res %>% as.data.frame() %>% mutate(cluster="Cluster 3")
hclust4_res_all = hclust4_res %>% as.data.frame() %>% mutate(cluster="Cluster 4")
hclust5_res_all = hclust5_res %>% as.data.frame() %>% mutate(cluster="Cluster 5")
all_hclust_deseq_res_all = rbind(hclust1_res_all,hclust2_res_all,hclust3_res_all,hclust4_res_all,hclust5_res_all)
write.table(all_hclust_deseq_res_all,file = "./outputs/DEGs/02_DESEQ2_hclust_expr_sites_combined_results_all.tsv", sep = "\t", quote = FALSE)

## Subset DESEQ2 markers to positives only
hclust1_res_pos = hclust1_res[hclust1_res$log2FoldChange >= 0,] %>% as.data.frame() %>% mutate(cluster="Cluster 1")
hclust2_res_pos = hclust2_res[hclust2_res$log2FoldChange >= 0,]  %>% as.data.frame() %>% mutate(cluster="Cluster 2")
hclust3_res_pos = hclust3_res[hclust3_res$log2FoldChange >= 0,] %>% as.data.frame() %>% mutate(cluster="Cluster 3")
hclust4_res_pos = hclust4_res[hclust4_res$log2FoldChange >= 0,] %>% as.data.frame() %>% mutate(cluster="Cluster 4")
hclust5_res_pos = hclust5_res[hclust5_res$log2FoldChange >= 0,] %>% as.data.frame() %>% mutate(cluster="Cluster 5")
## Combine into one
all_hclust_deseq_res_pos = rbind(hclust1_res_pos,hclust2_res_pos,hclust3_res_pos,hclust4_res_pos,hclust5_res_pos)
## Write the combined results
write.table(all_hclust_deseq_res_pos,file = "./outputs/DEGs/02_DESEQ2_hclust_expr_sites_combined_results_pos.tsv", sep = "\t", quote = FALSE)


all_hclust_deseq_res_pos

## Plot hclust volcanoes
all_hclust_deseq_res_pos_filtered = all_hclust_deseq_res_pos %>% filter(baseMean > 100)
all_hclust_deseq_res_pos_filtered = all_hclust_deseq_res_pos_filtered %>% mutate(padj_sig = padj < 0.05)

hclust_1_deseq_pos_filtered = all_hclust_deseq_res_pos_filtered %>% filter(cluster == "Cluster 1")
hclust1_volcano = ggplot(hclust_1_deseq_pos_filtered,aes(x=log2FoldChange,y=-log10(padj),color=padj_sig)) +
  geom_point() +
  geom_text_repel(
    data = subset(hclust_1_deseq_pos_filtered, (log2FoldChange > 3)|(-log10(padj)>18)),
    aes(label = gene),size=2,color="black",box.padding = 0.2
  ) +
  pretty_plot() + L_border() + geom_hline(yintercept = -log10(0.05),linetype="dashed") +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) + 
  theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
hclust1_volcano
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hclust1_volcano.png",hclust1_volcano,dpi=600,height=3,width=3.8)

hclust_2_deseq_pos_filtered = all_hclust_deseq_res_pos_filtered %>% filter(cluster == "Cluster 2")
hclust2_volcano = ggplot(hclust_2_deseq_pos_filtered,aes(x=log2FoldChange,y=-log10(padj),color=padj_sig)) +
  geom_point() +
  geom_text_repel(
    data = subset(hclust_2_deseq_pos_filtered, (log2FoldChange > 5)|(-log10(padj)>7)),
    aes(label = gene),size=2,color="black",box.padding = 0.2,max.overlaps = Inf
  ) +
  pretty_plot() + L_border() + geom_hline(yintercept = -log10(0.05),linetype="dashed") +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) + 
  theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
hclust2_volcano
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hclust2_volcano.png",hclust2_volcano,dpi=600,height=3,width=3.8)

hclust_3_deseq_pos_filtered = all_hclust_deseq_res_pos_filtered %>% filter(cluster == "Cluster 3")
hclust3_volcano = ggplot(hclust_3_deseq_pos_filtered,aes(x=log2FoldChange,y=-log10(padj),color=padj_sig)) +
  geom_point() +
  geom_text_repel(
    data = subset(hclust_3_deseq_pos_filtered, (log2FoldChange > 10)|(-log10(padj)>15)),
    aes(label = gene),size=2,color="black",box.padding = 0.2,max.overlaps = Inf
  ) +
  pretty_plot() + L_border() + geom_hline(yintercept = -log10(0.05),linetype="dashed") +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) + 
  theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
hclust3_volcano
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hclust3_volcano.png",hclust3_volcano,dpi=600,height=3,width=3.8)

hclust_4_deseq_pos_filtered = all_hclust_deseq_res_pos_filtered %>% filter(cluster == "Cluster 4")
hclust4_volcano = ggplot(hclust_4_deseq_pos_filtered,aes(x=log2FoldChange,y=-log10(padj),color=padj_sig)) +
  geom_point() +
  geom_text_repel(
    data = subset(hclust_4_deseq_pos_filtered, (log2FoldChange > 3)|(-log10(padj)>8)),
    aes(label = gene),size=2,color="black",box.padding = 0.2,max.overlaps = Inf
  ) +
  pretty_plot() + L_border() + geom_hline(yintercept = -log10(0.05),linetype="dashed") +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) + 
  theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
hclust4_volcano
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hclust4_volcano.png",hclust4_volcano,dpi=600,height=3,width=3.8)

hclust_5_deseq_pos_filtered = all_hclust_deseq_res_pos_filtered %>% filter(cluster == "Cluster 5")
hclust5_volcano = ggplot(hclust_5_deseq_pos_filtered,aes(x=log2FoldChange,y=-log10(padj),color=padj_sig)) +
  geom_point() +
  geom_text_repel(
    data = subset(hclust_5_deseq_pos_filtered, (log2FoldChange > 3)|(-log10(padj)>5)),
    aes(label = gene),size=2,color="black",box.padding = 0.2,max.overlaps = Inf
  ) +
  pretty_plot() + L_border() + geom_hline(yintercept = -log10(0.05),linetype="dashed") +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray")) + 
  theme(legend.position = "none",axis.title.x = element_blank(),axis.title.y = element_blank())
hclust5_volcano
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hclust5_volcano.png",hclust5_volcano,dpi=600,height=3,width=3.8)

## Now also do pathway enrichment per hclust
## Declare the gene universe
fora_gene_universe = rownames(so@assays$RNA$counts)
hclust1_pos_degs = hclust_1_deseq_pos_filtered %>% filter(padj < 0.05) %>% pull(gene)

hclust1_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,hclust1_pos_degs,fora_gene_universe) %>% as.data.frame()
hclust1_pos_degs_fora_hallmark$padj_sig = hclust1_pos_degs_fora_hallmark$padj < 0.05
hclust1_pos_degs_fora_hallmark = add_enriched_or(hclust1_pos_degs_fora_hallmark, n_degs = length(hclust1_pos_degs), universe_size = length(fora_gene_universe))
hclust1_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",hclust1_pos_degs_fora_hallmark$pathway)
hclust1_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",hclust1_pos_degs_fora_hallmark$pretty_pathway)

hclust1_pos_degs_fora_c6_up = fora(fgsea_c6_set_up_only,hclust1_pos_degs,fora_gene_universe) %>% as.data.frame()
hclust1_pos_degs_fora_c6_up$padj_sig = hclust1_pos_degs_fora_c6_up$padj < 0.05
hclust1_pos_degs_fora_c6_up = add_enriched_or(hclust1_pos_degs_fora_c6_up, n_degs = length(hclust1_pos_degs), universe_size = length(fora_gene_universe))
hclust1_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_UP","",hclust1_pos_degs_fora_c6_up$pathway)
hclust1_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_"," ",hclust1_pos_degs_fora_c6_up$pretty_pathway)

hclust1_pos_degs_fora_hallmark_plot = ggplot(hclust1_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(hclust1_pos_degs_fora_hallmark, padj_sig),
    aes(label = pretty_pathway), size=2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
hclust1_pos_degs_fora_hallmark_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hclust1_pos_degs_fora_hallmark_plot.pdf",hclust1_pos_degs_fora_hallmark_plot,dpi=300,height=1.4,width=2)
#save_deg_csv(noncut_breast_pos_degs_fora_hallmark,"./outputs/DEGs/noncut_breast_pos_degs_fora_hallmark.csv")

hclust1_pos_degs_fora_c6_up_plot = ggplot(hclust1_pos_degs_fora_c6_up,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(hclust1_pos_degs_fora_c6_up, -log10(padj)>2),
    aes(label = pretty_pathway),max.overlaps = Inf, size=2
  ) +
  theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
hclust1_pos_degs_fora_c6_up_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hclust1_pos_degs_fora_c6_up_plot.pdf",hclust1_pos_degs_fora_c6_up_plot,dpi=300,height=1.4,width=2)

## hclust2
hclust2_pos_degs = hclust_2_deseq_pos_filtered %>% filter(padj < 0.05) %>% pull(gene)

hclust2_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,hclust2_pos_degs,fora_gene_universe) %>% as.data.frame()
hclust2_pos_degs_fora_hallmark$padj_sig = hclust2_pos_degs_fora_hallmark$padj < 0.05
hclust2_pos_degs_fora_hallmark = add_enriched_or(hclust2_pos_degs_fora_hallmark, n_degs = length(hclust2_pos_degs), universe_size = length(fora_gene_universe))
hclust2_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",hclust2_pos_degs_fora_hallmark$pathway)
hclust2_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",hclust2_pos_degs_fora_hallmark$pretty_pathway)

hclust2_pos_degs_fora_c6_up = fora(fgsea_c6_set_up_only,hclust2_pos_degs,fora_gene_universe) %>% as.data.frame()
hclust2_pos_degs_fora_c6_up$padj_sig = hclust2_pos_degs_fora_c6_up$padj < 0.05
hclust2_pos_degs_fora_c6_up = add_enriched_or(hclust2_pos_degs_fora_c6_up, n_degs = length(hclust2_pos_degs), universe_size = length(fora_gene_universe))
hclust2_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_UP","",hclust2_pos_degs_fora_c6_up$pathway)
hclust2_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_"," ",hclust2_pos_degs_fora_c6_up$pretty_pathway)

hclust2_pos_degs_fora_hallmark_plot = ggplot(hclust2_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(hclust2_pos_degs_fora_hallmark, padj_sig & (odds_ratio > 3)),
    aes(label = pretty_pathway), size=2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
hclust2_pos_degs_fora_hallmark_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hclust2_pos_degs_fora_hallmark_plot.pdf",hclust2_pos_degs_fora_hallmark_plot,dpi=300,height=1.4,width=2)
#save_deg_csv(noncut_breast_pos_degs_fora_hallmark,"./outputs/DEGs/noncut_breast_pos_degs_fora_hallmark.csv")

hclust2_pos_degs_fora_c6_up_plot = ggplot(hclust2_pos_degs_fora_c6_up,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(hclust2_pos_degs_fora_c6_up, (-log10(padj)>2) & (odds_ratio > 3) | (-log10(padj)>3)),
    aes(label = pretty_pathway),max.overlaps = Inf, size=2
  ) +
  theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
hclust2_pos_degs_fora_c6_up_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hclust2_pos_degs_fora_c6_up_plot.pdf",hclust2_pos_degs_fora_c6_up_plot,dpi=300,height=1.4,width=2)

## hclust3
hclust3_pos_degs = hclust_3_deseq_pos_filtered %>% filter(padj < 0.05) %>% pull(gene)

hclust3_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,hclust3_pos_degs,fora_gene_universe) %>% as.data.frame()
hclust3_pos_degs_fora_hallmark$padj_sig = hclust3_pos_degs_fora_hallmark$padj < 0.05
hclust3_pos_degs_fora_hallmark = add_enriched_or(hclust3_pos_degs_fora_hallmark, n_degs = length(hclust3_pos_degs), universe_size = length(fora_gene_universe))
hclust3_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",hclust3_pos_degs_fora_hallmark$pathway)
hclust3_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",hclust3_pos_degs_fora_hallmark$pretty_pathway)

hclust3_pos_degs_fora_c6_up = fora(fgsea_c6_set_up_only,hclust3_pos_degs,fora_gene_universe) %>% as.data.frame()
hclust3_pos_degs_fora_c6_up$padj_sig = hclust3_pos_degs_fora_c6_up$padj < 0.05
hclust3_pos_degs_fora_c6_up = add_enriched_or(hclust3_pos_degs_fora_c6_up, n_degs = length(hclust3_pos_degs), universe_size = length(fora_gene_universe))
hclust3_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_UP","",hclust3_pos_degs_fora_c6_up$pathway)
hclust3_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_"," ",hclust3_pos_degs_fora_c6_up$pretty_pathway)

hclust3_pos_degs_fora_hallmark_plot = ggplot(hclust3_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(hclust3_pos_degs_fora_hallmark, padj_sig & (odds_ratio > 3.5)),
    aes(label = pretty_pathway), size=2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
hclust3_pos_degs_fora_hallmark_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hclust3_pos_degs_fora_hallmark_plot.pdf",hclust3_pos_degs_fora_hallmark_plot,dpi=300,height=1.4,width=2)

hclust3_pos_degs_fora_c6_up_plot = ggplot(hclust3_pos_degs_fora_c6_up,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(hclust3_pos_degs_fora_c6_up, (-log10(padj)>3.5) | (odds_ratio > 4)),
    aes(label = pretty_pathway),max.overlaps = Inf, size=2
  ) +
  theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
hclust3_pos_degs_fora_c6_up_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hclust3_pos_degs_fora_c6_up_plot.pdf",hclust3_pos_degs_fora_c6_up_plot,dpi=300,height=1.4,width=2)

## hclust4
hclust4_pos_degs = hclust_4_deseq_pos_filtered %>% filter(padj < 0.05) %>% pull(gene)

hclust4_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,hclust4_pos_degs,fora_gene_universe) %>% as.data.frame()
hclust4_pos_degs_fora_hallmark$padj_sig = hclust4_pos_degs_fora_hallmark$padj < 0.05
hclust4_pos_degs_fora_hallmark = add_enriched_or(hclust4_pos_degs_fora_hallmark, n_degs = length(hclust4_pos_degs), universe_size = length(fora_gene_universe))
hclust4_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",hclust4_pos_degs_fora_hallmark$pathway)
hclust4_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",hclust4_pos_degs_fora_hallmark$pretty_pathway)

hclust4_pos_degs_fora_c6_up = fora(fgsea_c6_set_up_only,hclust4_pos_degs,fora_gene_universe) %>% as.data.frame()
hclust4_pos_degs_fora_c6_up$padj_sig = hclust4_pos_degs_fora_c6_up$padj < 0.05
hclust4_pos_degs_fora_c6_up = add_enriched_or(hclust4_pos_degs_fora_c6_up, n_degs = length(hclust4_pos_degs), universe_size = length(fora_gene_universe))
hclust4_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_UP","",hclust4_pos_degs_fora_c6_up$pathway)
hclust4_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_"," ",hclust4_pos_degs_fora_c6_up$pretty_pathway)

hclust4_pos_degs_fora_hallmark_plot = ggplot(hclust4_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(hclust4_pos_degs_fora_hallmark, (odds_ratio > 5)&(-log10(padj)>20)),
    aes(label = pretty_pathway), size=2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
hclust4_pos_degs_fora_hallmark_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hclust4_pos_degs_fora_hallmark_plot.pdf",hclust4_pos_degs_fora_hallmark_plot,dpi=300,height=1.4,width=2)

hclust4_pos_degs_fora_c6_up_plot = ggplot(hclust4_pos_degs_fora_c6_up,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(hclust4_pos_degs_fora_c6_up, (-log10(padj)>7) | (odds_ratio > 3.7)),
    aes(label = pretty_pathway),max.overlaps = Inf, size=2
  ) +
  theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
hclust4_pos_degs_fora_c6_up_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hclust4_pos_degs_fora_c6_up_plot.pdf",hclust4_pos_degs_fora_c6_up_plot,dpi=300,height=1.4,width=2)

## hclust5
hclust5_pos_degs = hclust_5_deseq_pos_filtered %>% filter(padj < 0.05) %>% pull(gene)

hclust5_pos_degs_fora_hallmark = fora(fgsea_hallmark_set,hclust5_pos_degs,fora_gene_universe) %>% as.data.frame()
hclust5_pos_degs_fora_hallmark$padj_sig = hclust5_pos_degs_fora_hallmark$padj < 0.05
hclust5_pos_degs_fora_hallmark = add_enriched_or(hclust5_pos_degs_fora_hallmark, n_degs = length(hclust5_pos_degs), universe_size = length(fora_gene_universe))
hclust5_pos_degs_fora_hallmark$pretty_pathway = gsub("\\HALLMARK_","",hclust5_pos_degs_fora_hallmark$pathway)
hclust5_pos_degs_fora_hallmark$pretty_pathway = gsub("\\_"," ",hclust5_pos_degs_fora_hallmark$pretty_pathway)

hclust5_pos_degs_fora_c6_up = fora(fgsea_c6_set_up_only,hclust5_pos_degs,fora_gene_universe) %>% as.data.frame()
hclust5_pos_degs_fora_c6_up$padj_sig = hclust5_pos_degs_fora_c6_up$padj < 0.05
hclust5_pos_degs_fora_c6_up = add_enriched_or(hclust5_pos_degs_fora_c6_up, n_degs = length(hclust5_pos_degs), universe_size = length(fora_gene_universe))
hclust5_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_UP","",hclust5_pos_degs_fora_c6_up$pathway)
hclust5_pos_degs_fora_c6_up$pretty_pathway = gsub("\\_"," ",hclust5_pos_degs_fora_c6_up$pretty_pathway)

hclust5_pos_degs_fora_hallmark_plot = ggplot(hclust5_pos_degs_fora_hallmark,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(hclust5_pos_degs_fora_hallmark, (odds_ratio > 3.7)),
    aes(label = pretty_pathway), size=2
  ) + theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
hclust5_pos_degs_fora_hallmark_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hclust5_pos_degs_fora_hallmark_plot.pdf",hclust5_pos_degs_fora_hallmark_plot,dpi=300,height=1.4,width=2)

hclust5_pos_degs_fora_c6_up_plot = ggplot(hclust5_pos_degs_fora_c6_up,aes(x=odds_ratio,y=-log10(padj),color=padj_sig)) + 
  geom_point() + pretty_plot() + L_border() +
  geom_text_repel(
    data = subset(hclust5_pos_degs_fora_c6_up, (-log10(padj)>1) | (odds_ratio > 10)),
    aes(label = pretty_pathway),max.overlaps = Inf, size=2
  ) +
  theme(legend.position = "none", axis.title.x=element_blank(),axis.title.y=element_blank()) +
  scale_color_manual(values=c("TRUE"="dodgerblue3","FALSE"="gray"))
hclust5_pos_degs_fora_c6_up_plot
cowplot::ggsave2("./outputs/plots/04_enrichment_plots/hclust5_pos_degs_fora_c6_up_plot.pdf",hclust5_pos_degs_fora_c6_up_plot,dpi=300,height=1.4,width=2)

## Combine all DEGs and write it out
hclust1_pos_degs_fora_hallmark$cluster = "Cluster 1"
hclust1_pos_degs_fora_hallmark$term_set = "Hallmarks"
hclust1_pos_degs_fora_c6_up$cluster = "Cluster 1"
hclust1_pos_degs_fora_c6_up$term_set = "C6"

hclust2_pos_degs_fora_hallmark$cluster = "Cluster 2"
hclust2_pos_degs_fora_hallmark$term_set = "Hallmarks"
hclust2_pos_degs_fora_c6_up$cluster = "Cluster 2"
hclust2_pos_degs_fora_c6_up$term_set = "C6"

hclust3_pos_degs_fora_hallmark$cluster = "Cluster 3"
hclust3_pos_degs_fora_hallmark$term_set = "Hallmarks"
hclust3_pos_degs_fora_c6_up$cluster = "Cluster 3"
hclust3_pos_degs_fora_c6_up$term_set = "C6"

hclust4_pos_degs_fora_hallmark$cluster = "Cluster 4"
hclust4_pos_degs_fora_hallmark$term_set = "Hallmarks"
hclust4_pos_degs_fora_c6_up$cluster = "Cluster 4"
hclust4_pos_degs_fora_c6_up$term_set = "C6"

hclust5_pos_degs_fora_hallmark$cluster = "Cluster 5"
hclust5_pos_degs_fora_hallmark$term_set = "Hallmarks"
hclust5_pos_degs_fora_c6_up$cluster = "Cluster 5"
hclust5_pos_degs_fora_c6_up$term_set = "C6"

hclust_deg_hallmark_c6_combined = rbind(
  hclust1_pos_degs_fora_hallmark,hclust2_pos_degs_fora_hallmark,
  hclust3_pos_degs_fora_hallmark,hclust4_pos_degs_fora_hallmark,
  hclust5_pos_degs_fora_hallmark,
  hclust1_pos_degs_fora_c6_up,hclust2_pos_degs_fora_c6_up,
  hclust3_pos_degs_fora_c6_up,hclust4_pos_degs_fora_c6_up,
  hclust5_pos_degs_fora_c6_up
)

save_deg_csv(hclust_deg_hallmark_c6_combined,"./outputs/DEGs/hclust_deg_hallmark_c6_combined.csv")
