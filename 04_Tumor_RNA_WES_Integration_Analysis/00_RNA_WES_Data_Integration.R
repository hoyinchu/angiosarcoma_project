library(Seurat)
library(cowplot)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)
library(EnhancedVolcano)

## Load MAF 
asc_maf_path = "../data/processed/tumor_WES/ASC_mutations.maf"
asc_maf_metadata_path = "../data/processed/tumor_WES/ASC_mutations_metadata.tsv"
asc_cn_call_path = "../data/processed/tumor_WES/combined_cnvkit_calls.csv"

asc_cn_df = read.csv(asc_cn_call_path)
asc_cn_df_filtered = asc_cn_df %>% filter(CN == "Amp" | CN == "DeepDel")
asc_maf = read.maf(maf=asc_maf_path,clinicalData=asc_maf_metadata_path,cnTable = asc_cn_df_filtered)

## Load SeuratObject
so = LoadSeuratRds("../data/processed/rna/ASCSeuratObj.rds")

## Look for intersecting samples
intersect_samples = intersect(so@meta.data$sample_alias,unique(asc_maf@clinical.data$sample_alias))

# 2. Create a metadata slice from MAF containing the Tumor_Sample_Barcode and the sample alias
so@meta.data$seurat_id = row.names(so@meta.data)
so_meta_short = so@meta.data[,c("seurat_id","sample_alias_cleaned")] %>% rename("sample_alias_cleaned" = "sample_alias")
maf_metadata = asc_maf@clinical.data[, .(Tumor_Sample_Barcode, sample_alias)]
mapping_df = merge(so_meta_short, maf_metadata, by = "sample_alias", all.x = FALSE)

## Subset to intersecting samples
so_subset = subset(so, subset = sample_alias %in% intersect_samples)

asc_maf_subset_tsb = unique(asc_maf@clinical.data[asc_maf@clinical.data$sample_alias %in% intersect_samples,] %>% pull(Tumor_Sample_Barcode))
asc_maf_subset = subsetMaf(asc_maf,tsb=asc_maf_subset_tsb)

## Get recurrently mutated genes
asc_subset_gene_summary = getGeneSummary(asc_maf_subset)
asc_subset_gene_summary

high_freq_genes = asc_subset_gene_summary %>% filter(AlteredSamples >= 5) %>% pull(Hugo_Symbol)

## Quick inspection of top mutated genes
oncoplot(maf = asc_maf_subset,
         genes = high_freq_genes,
         clinicalFeatures = c("Primary_Site_(Recombined)", "CUTANEOUS_AS_(EHR_EXTRACTED)"),
         #topBarData = "TMB",
         draw_titv = TRUE,
         sortByAnnotation = TRUE,
         fontSize = 0.8)


genes_to_check = c("POT1","TP53","KDR","MYC","PLCG1")
samples_with_muts = genesToBarcodes(maf = asc_maf_subset, genes = genes_to_check, justNames = TRUE)
mut_sample_barcodes = unique(unname(unlist(samples_with_muts)))
asc_maf_non_driver_subset = subsetMaf(
  maf = asc_maf_subset, 
  tsb = setdiff(unique(asc_maf_subset@clinical.data$Tumor_Sample_Barcode), mut_sample_barcodes)
)

oncoplot(maf = asc_maf_non_driver_subset,
         genes = high_freq_genes,
         clinicalFeatures = c("Primary_Site_(Recombined)", "CUTANEOUS_AS_(EHR_EXTRACTED)"),
         #topBarData = "TMB",
         draw_titv = TRUE,
         sortByAnnotation = TRUE,
         fontSize = 0.8)

## Add this as metadata to the seurat object
gene_mut_matrix = genesToBarcodes(maf = asc_maf_subset, genes = high_freq_genes, justNames = TRUE)
for(gene in high_freq_genes){
  mut_samples <- gene_mut_matrix[[gene]]
  mut_seurat_ids <- mapping_df$seurat_id[mapping_df$Tumor_Sample_Barcode %in% mut_samples]
  so[[gene]] <- ifelse(colnames(so) %in% mut_seurat_ids, "Mutant", "WT")
}

table(so@meta.data$PLCG1)

VlnPlot(so,features = c("rna_MYC","rna_MYCN","rna_POT1","rna_KDR","rna_ERBB2"),group.by = "MYC")

# Identify signature
Idents(so) = "POT1"
pot1_signature <- FindMarkers(
  object = so,
  ident.1 = "Mutant",
  ident.2 = "WT",
  group.by = "POT1",
  test.use = "LR",           # Logistic Regression
  latent.vars = "ESTIMATE_purity", # The adjustment variable
  logfc.threshold = 0
)

## Make volcano plot
pot1_signature$padj_significant = pot1_signature$p_val_adj < 0.05
pot1_signature$gene = rownames(pot1_signature)
pot1_sig_volcano = ggplot(pot1_signature,aes(x=avg_log2FC,y=-log2(p_val),color=padj_significant)) +
  geom_point() + pretty_plot() + L_border() + scale_color_manual(values=c("TRUE"="firebrick","FALSE"="black")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  geom_vline(xintercept = 0,linetype="dashed") +
  #geom_text_repel(data=pot1_signature %>% filter(padj_significant),aes(label=gene)) +
  theme(legend.position = "none",axis.title = element_blank())

pot1_sig_volcano
cowplot::ggsave2("./outputs/plots/POT1_volcano.png",pot1_sig_volcano,dpi=600,width=1.6,height=1.3)
write.csv(pot1_signature,"./outputs/tables/POT1_volcano.csv")

pot1_signature

myc_signature <- FindMarkers(
  object = so,
  ident.1 = "Mutant",
  ident.2 = "WT",
  group.by = "MYC",
  test.use = "LR",           # Logistic Regression
  latent.vars = "ESTIMATE_purity", # The adjustment variable
  logfc.threshold = 0
)

myc_signature$gene = rownames(myc_signature)
myc_signature$to_highlight = myc_signature$gene %in% c("MYC")
myc_sig_volcano = ggplot(myc_signature %>% arrange(to_highlight),aes(x=avg_log2FC,y=-log2(p_val),color=to_highlight)) +
  geom_point() + pretty_plot() + L_border() +
  geom_vline(xintercept=0,linetype="dashed") +
  geom_hline(yintercept=-log2(0.05),linetype="dashed") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
myc_sig_volcano

write.csv(myc_signature,"./outputs/tables/MYC_volcano.csv")

myc_signature_filtered = myc_signature %>% filter(p_val < 0.05,avg_log2FC > 0.5) %>% arrange(-avg_log2FC)


tp53_signature <- FindMarkers(
  object = so,
  ident.1 = "Mutant",
  ident.2 = "WT",
  group.by = "TP53",
  test.use = "LR",           # Logistic Regression
  latent.vars = "ESTIMATE_purity", # The adjustment variable
  logfc.threshold = 0.25
)

tp53_signature_filtered = tp53_signature %>% filter(p_val_adj < 0.05)
write.csv(tp53_signature,"../unused/tp53_padj_sig.csv")

kdr_signature <- FindMarkers(
  object = so,
  ident.1 = "Mutant",
  ident.2 = "WT",
  group.by = "KDR",
  test.use = "LR",           # Logistic Regression
  latent.vars = "ESTIMATE_purity", # The adjustment variable
  logfc.threshold = 0.25
)

plcg1_signature <- FindMarkers(
  object = so,
  ident.1 = "Mutant",
  ident.2 = "WT",
  group.by = "PLCG1",
  test.use = "LR",           # Logistic Regression
  latent.vars = "ESTIMATE_purity", # The adjustment variable
  logfc.threshold = 0.1
)

ptprb_signature <- FindMarkers(
  object = so,
  ident.1 = "Mutant",
  ident.2 = "WT",
  group.by = "PTPRB",
  test.use = "LR",           # Logistic Regression
  latent.vars = "ESTIMATE_purity", # The adjustment variable
  logfc.threshold = 0.1
)

# ## Adding metadata
# mutation_cols = c("Total_Mutations","TMB_all_mutations","Total_Nonsyn_Mutations","TMB_nonsyn")
# sbs_cols = names(tumor_WES_meta)[grepl("SBS",names(tumor_WES_meta))]
# one_hot_cols = names(tumor_WES_onehot_data)
# mutation_subset_cols_to_keep = c("Tumor_Sample_Barcode","sample_alias",mutation_cols,sbs_cols)
# tumor_WES_meta_subset = tumor_WES_meta[,mutation_subset_cols_to_keep]

## Add one-hot somatic mutations
tumor_WES_onehot_data = read.csv("../data/processed/tumor_WES/ASC_mutations_MutSig_Recurrent_Gene_One_Hot.tsv",sep="\t",check.names = FALSE)
## Change column names to disambiguate
colnames(tumor_WES_onehot_data) = paste(colnames(tumor_WES_onehot_data),"mutation", sep = "_")
names(tumor_WES_onehot_data)[names(tumor_WES_onehot_data) == 'Tumor_Sample_Barcode_mutation'] = 'Tumor_Sample_Barcode'
tumor_WES_meta_subset_merged = merge(tumor_WES_meta_subset,tumor_WES_onehot_data,on="Tumor_Sample_Barcode",all.x=TRUE)

## Add to Seurat Metadata
seurat_meta_expanded = merge(so_subset@meta.data,tumor_WES_meta_subset_merged,by="sample_alias",all.x=TRUE)
seurat_meta_expanded = seurat_meta_expanded %>% select("og_id", everything())
so_subset@meta.data = seurat_meta_expanded
rownames(so_subset@meta.data) = so_subset@meta.data$og_id

## Re-scale the data w.r.t. to the rest of the subset samples
vst_rescaled = so_subset@assays$RNA$counts
vst_rescaled = as.data.frame(t(scale(t(vst_rescaled))),check.names=FALSE)
so_subset@assays$RNA$vst_rescaled = vst_rescaled

FeaturePlot(so,features = c("rna_PLCG1"))

## Check SBS by cluster
VlnPlot(so_subset,features = c("SBS6","SBS7","SBS15","SBS87"))

## Check Shelterin-related gene expression in POT1-mutant samples
telomere_related_geneset = c(
  "POT1","TERT","ATRX","DAXX",
  "TERF1","TERF2","TINF2","RPA",
  "TPP1","WRN","ATM","ATR",
  "BRCA1","BRCA2","RAD51","BLM","BARD1",
)
so_subset = SetIdent(so_subset, value = so_subset@meta.data$POT1_mutation)
pot1_mutant_heatmap = DoHeatmap(so_subset, features = telomere_related_geneset, 
                                disp.max = 3.5,disp.min = -3.5,slot="vst_rescaled")
pot1_mutant_heatmap
ggsave(pot1_mutant_heatmap,filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_POT1_mutant_shelterin_expression_plot.png",dpi=300)

## Check telomere-related gene expression in POT1-mutant samples in parenchymal breast samples only
so_site_subset = subset(so_subset, idents = "BREAST (PARENCHYMAL)")
so_site_subset = SetIdent(so_site_subset, value = so_site_subset@meta.data$POT1_mutation)
pot1_mutant_parenchymal_breast_heatmap = DoHeatmap(so_site_subset, features = telomere_related_geneset,disp.max = 3.5,disp.min = -3.5,slot="vst_rescaled")
pot1_mutant_parenchymal_breast_heatmap
ggsave(pot1_mutant_parenchymal_breast_heatmap,filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_POT1_mutant_telomere_related_expression_in_parenchymal_breast_plot.png",dpi=300)


## Check Shelterin-related gene expression in KDR-mutant samples
so_subset = SetIdent(so_subset, value = so_subset@meta.data$KDR_mutation)
kdr_mutant_heatmap = DoHeatmap(so_subset, features = telomere_related_geneset,
                               disp.max = 3.5,disp.min = -3.5,slot="vst_rescaled")
kdr_mutant_heatmap
ggsave(kdr_mutant_heatmap,filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_KDR_mutant_shelterin_expression_plot.png",dpi=300)

## Check general gene markers per primary sites
so_subset = SetIdent(so_subset, value = so_subset@meta.data$`PRIMARY SITE (Combined)`)
primary_site_telomere_genes_heatmap = DoHeatmap(so_subset, features = telomere_related_geneset,
                               disp.max = 3.5,disp.min = -3.5,slot="vst_rescaled"
)
primary_site_telomere_genes_heatmap
ggsave(primary_site_telomere_genes_heatmap,
       filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_telomere_related_genes_expression_plot_by_primary_sites.png",dpi=300,height=10,width=20)

so_subset@meta.data[so_subset@meta.data$FLT4_mutation==1,]

# so_subset = SetIdent(so_subset, value = so_subset@meta.data$`PRIMARY SITE (Combined)`)
# breast_kdr_diff_expr = FindMarkers(so_subset, ident.1 = 1, group.by = 'POT1_mutation', subset.ident = "BREAST (PARENCHYMAL)")
# so_site_subset = subset(so_subset, idents = "BREAST (PARENCHYMAL)")
# so_site_subset = SetIdent(so_site_subset, value = so_site_subset@meta.data$POT1_mutation)
# pos_enrich_kdr_genes = breast_kdr_diff_expr#[(breast_kdr_diff_expr$p_val<0.05) & (breast_kdr_diff_expr$avg_log2FC>0),]
#DoHeatmap(so_site_subset, features = rownames(pos_enrich_kdr_genes)[1:50],disp.max = 3.5,disp.min = -3.5,slot="vst_rescaled")
#DoHeatmap(so_site_subset, features = c("TERT","ATRX","DAXX"),disp.max = 3.5,disp.min = -3.5,slot="vst_rescaled")

pot1_subset = subset(so_subset, subset = POT1 < 1,slot="vst_rescaled")

rna_mutation_overlay = DimPlot(so_subset, reduction = "pca",group.by=c(
  "seurat_clusters","PRIMARY SITE (Combined)",
  "TP53_mutation","KDR_mutation","POT1_mutation",
  "FLT4_mutation","ASXL1_mutation","ATRX_mutation",
  "PTPRO_mutation","PLCG1_mutation","PTPRB_mutation",
  "SETBP1_mutation","ARID1A_mutation","LRP2_mutation"
  ),ncol=5,pt.size=4)

rna_mutation_overlay

## Check the correlation between gene expression and TMB
rna_vst_expression = so_subset@assays$RNA$counts
tmb_nonsyn = so_subset@meta.data$TMB_nonsyn
names(tmb_nonsyn) = rownames(so_subset@meta.data)

calc_correlation = function(expr_mat,feature_vec) {
  gene_feature_correlations = as.data.frame(t(apply(expr_mat, 1, function(gene_data) {
    test_res = cor.test(gene_data, feature_vec, use = "complete.obs", method = "pearson")
    parital_res = c(
      gene = rownames(gene_data)[[1]],
      correlation = test_res$estimate,
      pvalue = test_res$p.value,
      conf_low = test_res$conf.int[1],
      conf_high = test_res$conf.int[2]
    )
    return(parital_res)
  })),check.names=FALSE)
  gene_feature_correlations_ordered = gene_feature_correlations[order(gene_feature_correlations$pvalue),]
  gene_feature_correlations$neg_logp = -log(gene_feature_correlations$pvalue)
  gene_feature_correlations$neg_logp_capped = pmin(gene_feature_correlations$neg_logp,20)
  gene_feature_correlations$gene = rownames(gene_feature_correlations)
  return(gene_feature_correlations)
}

plot_correlation_volcano = function(corr_df) {
  correlation_volcano = ggplot(corr_df,aes(x=correlation.cor,y=neg_logp_capped,label=gene)) +
    geom_text_repel(data = subset(corr_df, pvalue < 0.05), max.overlaps=50) +
    geom_point() +
    geom_hline(yintercept = -log(0.05)) +
    geom_vline(xintercept = 0) +
    theme_minimal() +
    labs(x="TMB vs. Gene Expression Pearson Correlation",y="-log(p-value)")
  return(correlation_volcano)
}

tmb_gene_correlations = calc_correlation(rna_vst_expression,tmb_nonsyn)
tmb_correlation_volcano = plot_correlation_volcano(tmb_gene_correlations)
tmb_correlation_volcano
ggsave(correlation_volcano,filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_tmb_rna_expression_correlation_volcano.png",dpi=300)

## Remove those with TMB > 10 and redo
tmb_nonsyn_no_hyper = tmb_nonsyn[tmb_nonsyn<10]
rna_vst_expression_no_hyper = rna_vst_expression[,colnames(rna_vst_expression) %in% names(tmb_nonsyn_no_hyper)]
tmb_gene_correlations_no_hyper = calc_correlation(rna_vst_expression_no_hyper,tmb_nonsyn_no_hyper)
correlation_no_hyper_volcano = plot_correlation_volcano(tmb_gene_correlations_no_hyper)

## Restrict to same site only
hnfs_only = rownames(so_subset@meta.data[so_subset@meta.data$`PRIMARY SITE (Combined)` == "HNFS",])
tmb_nonsyn_hfns_only = tmb_nonsyn[hnfs_only]
rna_vst_expression_hfns_only = rna_vst_expression[,colnames(rna_vst_expression) %in% names(tmb_nonsyn_hfns_only)]
tmb_gene_correlations_hfns_only = calc_correlation(rna_vst_expression_hfns_only,tmb_nonsyn_hfns_only)
correlation_hfns_volcano = plot_correlation_volcano(tmb_gene_correlations_hfns_only)
correlation_hfns_volcano
ggsave(correlation_hfns_volcano,filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_tmb_rna_expression_correlation_hfns_only_volcano.png",dpi=300)

## Make POT1-TPP1 plots
so_subset = SetIdent(so_subset, value = so_subset@meta.data$POT1_mutation)
pot1_tpp1_scatter_by_mut = FeatureScatter(so_subset,feature1 = "POT1",feature2 = "TPP1")
ggsave(pot1_tpp1_scatter_by_mut,filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_pot1_vs_tpp1_scatter_by_pot1_mutation.png",dpi=300)

pot1_tpp1_violin = VlnPlot(so_subset,features = c("TPP1"))
ggsave(pot1_tpp1_violin,filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_pot1_vs_tpp1_violin.png",dpi=300)

so_subset = SetIdent(so_subset, value = so_subset@meta.data$seurat_clusters)
pot1_tpp1_scatter_by_cluster = FeatureScatter(so_subset,feature1 = "POT1",feature2 = "TPP1")
ggsave(pot1_tpp1_scatter_by_cluster,filename = "04_Tumor_RNA_WES_Integration_Analysis/outputs/plots/00_pot1_vs_tpp1_scatter_by_cluster.png",dpi=300)

## Check correlation between SBS7 and gene expression?
sbs7_expr = so_subset@meta.data$SBS1
sbs7_gene_correlations = calc_correlation(rna_vst_expression,sbs7_expr)
sbs7_gene_correlations_volcano = plot_correlation_volcano(sbs7_gene_correlations)
sbs7_gene_correlations_volcano

##
# ir_genes = read.csv("data/public/ionizing_radiation_up_genes.csv")
# so_subset = AddModuleScore(so_subset,features = list(c(ir_genes$gene)),name="ir_signature")
# FeaturePlot(so_subset,reduction = "pca",features = "ir_signature1",pt.size=4)
# so_subset = SetIdent(so_subset, value = so_subset@meta.data$seurat_clusters)
# VlnPlot(so_subset,features = "ir_signature1")
# FeatureScatter(so_subset,feature1 = "FLT4",feature2 = "FLT1")


# ## Check if Cluster 2&3 samples are separated by keratin expression
# so_subset = UpdateSlots(so_subset)
# colnames(so_subset) = Cells(so_subset[["RNA_snn"]])
# so_subset = SetIdent(so_subset, value = so_subset@meta.data$seurat_clusters)
# so_subset_cluster_2_and_3 = subset(so_subset, subset = seurat_clusters == 2)
# 
# 
# 
# VlnPlot(so_subset_cluster_2_and_3,features = c("KRT1"))
# FeaturePlot(so_subset_cluster_2_and_3,reduction = "pca",features = c("TMB_nonsyn"))
# FeatureScatter(so_subset_cluster_2_and_3,feature1 = "TMB_nonsyn",feature2 = "KRT4")

