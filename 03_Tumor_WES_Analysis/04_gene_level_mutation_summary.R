library(maftools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts")

# Load processed MAF file and clinical metadata
asc_maf_path = "data/processed/tumor_WES/ASC_mutations.maf"
asc_maf_metadata_path = "data/processed/tumor_WES/ASC_mutations_metadata.tsv"
asc_maf = read.maf(maf=asc_maf_path,clinicalData=asc_maf_metadata_path)

make_oncoplot = function(maf,save_path="") {
  if (save_path!="") {
    pdf(file=save_path,height=6)
    #pdf(file = "Tumor_WES_Analysis/outputs/01_oncoplot_all_by_pathways.pdf",height=6)
  }
  oncoplot(maf = maf,
           top = 20,
           clinicalFeatures=c("PRIMARY_SITE_(Combined)","CUTANEOUS_AS_(EHR_EXTRACTED)"),#,"LOCAL_RECURRENCE_(EHR_EXTRACTED)"),
           topBarData="TMB_nonsyn",
           draw_titv = TRUE,
           sortByAnnotation = TRUE,
           fontSize = 0.6,
  )
  if (save_path!="") {
    dev.off()
  }
}

# Overall oncoplot
make_oncoplot(asc_maf,save_path = "03_Tumor_WES_Analysis/outputs/plots/04_oncoplot_all.pdf")

## Subset to specific samples
## Non-HNFS
asc_maf_non_hnfs_tsb = asc_maf@clinical.data[asc_maf@clinical.data$`PRIMARY_SITE_(Combined)` != "HNFS"]$Tumor_Sample_Barcode
asc_maf_non_hnfs = subsetMaf(asc_maf,tsb = asc_maf_non_hnfs_tsb)
make_oncoplot(asc_maf_non_hnfs,save_path = "03_Tumor_WES_Analysis/outputs/plots/04_oncoplot_non_HNFS.pdf")
somaticInteractions(maf = asc_maf_non_hnfs, top = 25, pvalue = c(0.05, 0.1))

## HNFS
asc_maf_hnfs_tsb = asc_maf@clinical.data[asc_maf@clinical.data$`PRIMARY_SITE_(Combined)` == "HNFS"]$Tumor_Sample_Barcode
asc_maf_hnfs = subsetMaf(asc_maf,tsb = asc_maf_hnfs_tsb)
make_oncoplot(asc_maf_hnfs,save_path = "03_Tumor_WES_Analysis/outputs/plots/04_oncoplot_HNFS.pdf")
somaticInteractions(maf = asc_maf_hnfs, top = 25, pvalue = c(0.05, 0.1))

## Highlight Mutsig Significant genes
mutsig_output_path = "data/processed/tumor_WES/mutsig/Apr16_2024_sig_genes.txt"
mutsig_gene = read.csv(mutsig_output_path,sep="\t")

## Add number of mutations
gene_mutsig_merged = merge(asc_maf@gene.summary,mutsig_gene,by.x="Hugo_Symbol",by.y="gene",all.x=TRUE)
qval_threshold = 0.1
gene_mutsig_merged[["neg_log_p"]]=-log(gene_mutsig_merged$p)
gene_mutsig_merged[["q_pass_threshold"]]=gene_mutsig_merged$q < qval_threshold
mutsig_gene_sub = gene_mutsig_merged[gene_mutsig_merged$MutatedSamples>=2,] %>%
  drop_na(q_pass_threshold) %>%
  mutate(q_pass_threshold = ifelse(q_pass_threshold, "q < 0.1", "q >= 0.1"))

# Plot mutsig significant genes
mutsig_plot = ggplot(mutsig_gene_sub,aes(x=MutatedSamples,y=neg_log_p,label=Hugo_Symbol,color=q_pass_threshold)) +
  geom_text_repel(data = subset(mutsig_gene_sub, (p < 0.05 & MutatedSamples>=5) | (MutatedSamples>=10) | (q < 0.5))) +
  geom_point(size=3) +
  theme_minimal() +
  geom_hline(yintercept=-log(0.05),linetype="dashed",color="gray") +
  geom_text(aes(0,-log(0.05),label = "p = 0.05", vjust = -1)) +
  labs(color='Mutsig q-value', y = "-log(Mutsig p-value)", x="Number of Samples with Mutation")
mutsig_plot
ggsave(mutsig_plot,file="03_Tumor_WES_Analysis/outputs/plots/04_MutSig_num_samples_by_significance_plot.png",width=10,height=12)

## Identify recurrently mutated genes with nominal MutSig significance
asc_maf_onehot = table(asc_maf@data$Tumor_Sample_Barcode,asc_maf@data$Hugo_Symbol)
asc_maf_onehot = (asc_maf_onehot>0)*1
recurrent_genes = gene_mutsig_merged %>% filter(MutatedSamples >= 5, p < 0.05)

# Plot genes that are enriched by a clinical attributes
make_clin_enrichment_table = function(maf_df,clin_feat,clin_feat_label) {
  clin_enrich = clinicalEnrichment(
    maf = asc_maf,
    clinicalFeature = c(clin_feat),
    minMut = 5,
  )
  clin_enrich_group_wise = clin_enrich$groupwise_comparision
  clin_enrich_group_wise_renamed = clin_enrich_group_wise
  clin_enrich_group_wise_renamed$p_bin = cut(
    clin_enrich_group_wise_renamed$p_value,
    breaks = c(0,0.01,0.05,0.1,1.1),
    labels = c("p<0.01","p<0.05","p<0.1","p>=0.1")
  )
  clin_enrich_group_wise_sig_genes = clin_enrich_group_wise_renamed %>% filter(p_value < 0.1 & OR > 1)
  return(clin_enrich_group_wise_sig_genes)
}

make_clin_enrichment_plot = function(enrichment_table,clin_feat_label) {
  enrich_dotplot = ggplot(enrichment_table,aes(x=Group1,y=Hugo_Symbol,size=OR,color=p_bin)) +
    geom_point() +
    #geom_text(size=2) +
    theme_minimal() + 
    labs(color="Fisher p-value", x=clin_feat_label, y = "Genes", title="Gene Mutation Frequency Enrichment (one vs. rest)")
  return(enrich_dotplot)
}

# Genes enriched by primary sites
primary_site_enrichment_table = make_clin_enrichment_table(asc_maf,"PRIMARY_SITE_(Combined)")
primary_site_enrichment_table_subset = primary_site_enrichment_table[primary_site_enrichment_table$Hugo_Symbol %in% recurrent_genes$Hugo_Symbol,]
primary_site_enrichment_plot = make_clin_enrichment_plot(primary_site_enrichment_table_subset,"Primary Site")
primary_site_enrichment_plot
ggsave(primary_site_enrichment_plot,file="03_Tumor_WES_Analysis/outputs/plots/04_gene_mutaton_frequency_enrichment_by_primary_sites.png",width=9)


## Perform ORA on these recurrent mutated genes (GO BP and KEGG Pathways)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

mutsig_bg_genes = mutsig_gene$gene
mutsig_recurrent_sig_genes = recurrent_genes$Hugo_Symbol
mutsig_enrich_obj = enrichGO(
  gene = mutsig_recurrent_sig_genes,
  universe = mutsig_bg_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont="BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500
)
mutsig_sig_genes_uniprot = bitr(mutsig_recurrent_sig_genes,fromType="SYMBOL",toType="UNIPROT",OrgDb = org.Hs.eg.db)$UNIPROT
mutsig_sig_genes_uniprot = mutsig_sig_genes_uniprot[!duplicated(mutsig_sig_genes_uniprot)]
mutsig_bg_genes_uniprot = bitr(mutsig_bg_genes,fromType="SYMBOL",toType="UNIPROT",OrgDb = org.Hs.eg.db)$UNIPROT
mutsig_bg_genes_uniprot = mutsig_bg_genes_uniprot[!duplicated(mutsig_bg_genes_uniprot)]
kegg_gson = gson_KEGG("hsa",keyType="uniprot")
mutsig_kegg_obj = enricher(
  gene=mutsig_sig_genes_uniprot,
  universe=mutsig_bg_genes_uniprot,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  gson = kegg_gson,
  minGSSize = 10,
  maxGSSize = 500
)
# MSigDB oncogenic signature enrichment (upregulated gene sets only)
msigdb_c6_gene_sets = msigdbr(species = "Homo sapiens", category = "C6") %>% 
  dplyr::filter(grepl("up-regulated",gs_description)) %>%
  dplyr::select(gs_name, gene_symbol)
mutsig_msig_c6_obj = enricher(
  gene=mutsig_recurrent_sig_genes,
  universe=mutsig_bg_genes,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2,
  TERM2GENE=msigdb_c6_gene_sets,
  minGSSize = 10,
  maxGSSize = 500
)

mutsig_enrich_res = list(
  "go"=mutsig_enrich_obj,
  "kegg"=mutsig_kegg_obj,
  "msigdb_c6"=mutsig_msig_c6_obj
)

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

mutsig_go_bar = draw_enriched_bars(mutsig_enrich_res,"go","MutSig Recurrently Mutated Genes:\nGO-BP")
mutsig_kegg_bar = draw_enriched_bars(mutsig_enrich_res,"kegg", "MutSig Recurrently Mutated Genes:\nKEGG")
mutsig_msigdb_bar = draw_enriched_bars(mutsig_enrich_res,"msigdb_c6", "MutSig Recurrently Mutated Genes: Oncogenic Signatures")

#mutsig_bars_combined =  cowplot::plot_grid(mutsig_go_bar,mutsig_kegg_bar,mutsig_msigdb_bar,labels="AUTO",ncol=3)
# Plot the enrichment results
mutsig_bars_combined =  cowplot::plot_grid(mutsig_go_bar,mutsig_kegg_bar,labels="AUTO",ncol=2)
ggsave(mutsig_bars_combined,file="03_Tumor_WES_Analysis/outputs/plots/04_MutSig_recurrent_gene_enrichment_analysis.png",width=12)

## Make a onehot matrix for samples carrying recurrent mutations
asc_maf_onehot = asc_maf_onehot[,recurrent_genes$Hugo_Symbol]
asc_maf_onehot = as.data.frame.matrix(asc_maf_onehot)
asc_maf_onehot = cbind(Tumor_Sample_Barcode=rownames(asc_maf_onehot),asc_maf_onehot)
write.table(asc_maf_onehot,file="data/processed/tumor_WES/ASC_mutations_MutSig_Recurrent_Gene_One_Hot.tsv",sep="\t",row.names=FALSE,quote = FALSE)

## Jaccard similarity?
jaccard_dist_matrix = as.matrix(proxy::dist(asc_maf_onehot, method = "Jaccard",convert_similarities = FALSE))
ComplexHeatmap::Heatmap(jaccard_dist_matrix)


# asc_maf_onehot$Tumor_Sample_Barcode = rownames(asc_maf_onehot)
# #asc_maf_onehot = cbind(Tumor_Sample_Barcode=rownames(asc_maf_onehot),asc_maf_onehot)
# asc_maf_onehot_merged = merge(asc_maf_onehot,asc_maf@clinical.data,by="Tumor_Sample_Barcode",all.x=TRUE)
# 
# 
# 
# asc_maf_onehot_merged$TP53 = factor(asc_maf_onehot_merged$TP53)
# asc_maf_onehot_merged$TMB_nonsyn = as.numeric(asc_maf_onehot_merged$TMB_nonsyn)
# 
# asc_maf_onehot_merged
# 
# wilcox.test(TMB_nonsyn ~ TP53, data = asc_maf_onehot_merged)
# 
# tp53_tmb_plot = ggplot(asc_maf_onehot_merged,aes(x=TP53,y=TMB_nonsyn)) +
#   geom_boxplot()
# ggsave(tp53_tmb_plot,file="03_Tumor_WES_Analysis/outputs/plots/04_tp53_tmb_plot.png",height=12,width=12)
# 
# tmb_sbs7_plot = ggplot(asc_maf_onehot_merged,aes(x=TMB_nonsyn,y=SBS7)) +
#   geom_point()
# ggsave(tmb_sbs7_plot,file="03_Tumor_WES_Analysis/outputs/plots/04_tmb_sbs7_plot.png",height=12,width=12)
# 
# asc_maf_onehot_merged$`Age_(Combined)` = as.numeric(asc_maf_onehot_merged$`Age_(Combined)`)
# age_tmb_plot = ggplot(asc_maf_onehot_merged,aes(x=`Age_(Combined)`,y=TMB_nonsyn,color=`PRIMARY_SITE_(Combined)`)) +
#   geom_point()
# ggsave(age_tmb_plot,file="03_Tumor_WES_Analysis/outputs/plots/04_age_tmb_plot.png",height=12,width=12)



# # See if certian mutations are associated with clinical attributes
# clin_enrich = clinicalEnrichment(
#   maf = asc_maf,
#   clinicalFeature = c("PRIMARY_SITE_(Combined)"),
#   minMut = 5,
# )
# clin_enrich_group_wise = clin_enrich$groupwise_comparision
# clin_enrich_group_wise_renamed = clin_enrich_group_wise
# clin_enrich_group_wise_renamed = clin_enrich_group_wise %>%
#   mutate(Group1 = case_when(
#     Group1==0 ~ "cluster 0\n(Parenchymal Breast)",
#     Group1==1 ~ "cluster 1\n(Cutaneous Breast)",
#     Group1==2 ~ "cluster 2",
#     Group1==3 ~ "cluster 3\n(HFNS)",
#     Group1==4 ~ "cluster 4",
#   ))
# clin_enrich_group_wise_renamed$p_bin = cut(
#   clin_enrich_group_wise_renamed$p_value,
#   breaks = c(0,0.01,0.05,0.1,1.1),
#   labels = c("p<0.01","p<0.05","p<0.1","p>=0.1")
#   )
# clin_enrich_group_wise_sig_genes = clin_enrich_group_wise_renamed %>% filter(p_value < 0.1 & OR > 1)
# enrich_dotplot = ggplot(clin_enrich_group_wise_sig_genes,aes(x=Group1,y=Hugo_Symbol,size=OR,color=p_bin)) +
#   geom_point() +
#   #geom_text(size=2) +
#   theme_minimal() + 
#   labs(color="Fisher p-value", x="Primary Site", y = "Genes", title="Gene Mutation Frequency Enrichment (one vs. rest)")
#   #scale_size(range=c(2,10))
# enrich_dotplot
# ggsave(enrich_dotplot,file="Tumor_WES_Analysis/outputs/plots/01_gene_mutaton_frequency_enrichment_by_seurat_clusters.png",width=12)


# # Add genes that are cancer implicated
# cosmic_table = read.csv("data/public/cosmic_cancer_gene_census.csv",sep=",",header=TRUE)
# oncokb_table = read.csv("data/public/oncokb_cancer_gene_list.tsv",sep="\t",header=TRUE,check.names = FALSE)
# oncokb_table_filtered = oncokb_table %>% filter(`OncoKB Annotated` == "Yes" & `COSMIC CGC (v99)`=="Yes")
# oncokb_filtered_gene_list = oncokb_table_filtered$`Hugo Symbol`
# reported_gene_list <- c("PTPRB","CIC","TRF1","TRF2","RAP1","TIN2","TPP1") # Add shelterin genes
# angiosarcoma_gene_list <- c(oncokb_filtered_gene_list,reported_gene_list)
# asc_oncokb_maf = subsetMaf(maf = asc_maf, genes = angiosarcoma_gene_list)
# make_oncoplot(asc_oncokb_maf,save_path = "Tumor_WES_Analysis/outputs/plots/01_oncoplot_oncoKB_only.pdf")

# Plot tmb by gene alteration frequency
# gene_mutation_by_tmb =  asc_maf@clinical.data
