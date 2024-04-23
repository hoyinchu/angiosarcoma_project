library(DESeq2)
library(ggplot2)
library(ggrepel)
library(tidyr)

## Load filtered data and GSVA data
setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts")
count_data = read.csv("data/processed/rna/00_filtered_gene_counts.csv",row.names = 1,check.names = FALSE)
vst_data = read.csv("data/processed/rna/00_filtered_gene_counts.vst.csv",row.names = 1,check.names = FALSE)
tpm_data = read.csv("data/processed/rna/00_filtered_gene_tpm.csv",row.names = 1,check.names = FALSE)
#meta_data = read.csv("reference_data/RNA_Seq/outputs/filtered_metadata.csv",check.names = FALSE)
meta_data = read.csv("data/processed/rna/00_filtered_sample_metadata.csv",check.names = FALSE)

# 
# gsva_vst_go_bp_data = read.csv("reference_data/RNA_Seq/outputs/gsva/gsva_ssGSEA_score_vst_GO_BP_set.csv",row.names = 1,check.names = FALSE)
# gsva_vst_cancer_data = read.csv("reference_data/RNA_Seq/outputs/gsva/gsva_ssGSEA_score_vst_cancer_set.csv",row.names = 1,check.names = FALSE)
# gsva_log_tpm_go_bp_data = read.csv("reference_data/RNA_Seq/outputs/gsva/gsva_ssGSEA_score_log_tpm_GO_BP_set.csv",row.names = 1,check.names = FALSE)
# gsva_log_tpm_cancer_data = read.csv("reference_data/RNA_Seq/outputs/gsva/gsva_ssGSEA_score_log_tpm_cancer_set.csv",row.names = 1,check.names = FALSE)
# 
# singscore_tpm_hallmark_data = read.csv("reference_data/RNA_Seq/outputs/singscore/tpm_hallmark_set_singscores.csv",row.names = 1,check.names = FALSE)
# singscore_tpm_go_bp_data = read.csv("reference_data/RNA_Seq/outputs/singscore/tpm_singscore_go_bp_set.csv",row.names = 1,check.names = FALSE)
# singscore_tpm_paired_cancer_data = read.csv("reference_data/RNA_Seq/outputs/singscore/tpm_paired_cancer_set_singscores.csv",row.names = 1,check.names = FALSE)
# singscore_tpm_bidirection_cancer_data = read.csv("reference_data/RNA_Seq/outputs/singscore/tpm_bidirection_cancer_set_singscores.csv",row.names = 1,check.names = FALSE)
# singscore_tpm_paired_cgp_data = read.csv("reference_data/RNA_Seq/outputs/singscore/tpm_cgp_paired_set_singscores.csv",row.names = 1,check.names = FALSE)
# 
# loh_singscore_paired_cancer_data = read.csv("reference_data/RNA_Seq/outputs/singscore/loh_raw_counts_bidirection_cancer_set_singscores.csv",row.names = 1,check.names = FALSE)

## Plot the library sizes by samples
lib_size_plot = ggplot(meta_data, aes(x=library_size,fill=LC_BATCH)) + 
  geom_histogram() +
  scale_x_log10() +
  xlab("Library Size") +
  ylab("# of Samples")
lib_size_plot
ggsave("RNA_Analysis/outputs/plots/00_qc_plots/01_library_size_histogram.png",lib_size_plot)

## Plot logged gene count distribution
total_gene_counts = rowSums(count_data)
total_gene_counts = data.frame("log_count"=log(total_gene_counts[order(-total_gene_counts)]))
gene_counts_plot = ggplot(total_gene_counts, aes(x=log_count)) + 
  geom_histogram() +
  xlab("Log Gene Read Counts") +
  ylab("# of Genes")
gene_counts_plot
ggsave("RNA_Analysis/outputs/plots/00_qc_plots/01_gene_read_counts_histogram.png",gene_counts_plot)

## Plot VST transformed library size and gene size
vst_library_sizes = as.data.frame(colSums(vst_data))
colnames(vst_library_sizes) = "library_size"
vst_library_sizes_plot = ggplot(vst_library_sizes, aes(x=library_size)) + 
  geom_histogram() +
  xlab("VST-transformed Library Sizes") +
  ylab("# of Samples")
vst_library_sizes_plot
ggsave("RNA_Analysis/outputs/plots/00_qc_plots/01_vst_transformed_library_size_histogram.png",vst_library_sizes_plot)

vst_gene_counts = as.data.frame(rowSums(vst_data))
colnames(vst_gene_counts) = "gene_counts"
vst_gene_counts_plot = ggplot(vst_gene_counts, aes(x=gene_counts)) + 
  geom_histogram() +
  xlab("VST-transformed Gene Counts") +
  ylab("# of Genes with Read Count")
vst_gene_counts_plot
ggsave("RNA_Analysis/outputs/plots/00_qc_plots/01_vst_transformed_gene_counts_histogram.png",vst_gene_counts_plot)

## Plot mean-variance relationship between genes
mean_var_plot = function(gene_data) {
  gene_data_mean = apply(gene_data,1,mean)
  gene_data_var = apply(gene_data,1,var)
  q99_mean= quantile(gene_data_mean,probs=c(0.99))
  q99_var= quantile(gene_data_var,probs=c(0.99))
  mean_var_df = as.data.frame(list("gene_mean"=gene_data_mean,"gene_var"=gene_data_var))
  mean_var_df$gene = rownames(mean_var_df)
  #return(mean_var_df)
  mean_var_plot = ggplot(mean_var_df,aes(x=gene_mean,y=gene_var)) +
    geom_point() +
    geom_text_repel(
      data = subset(mean_var_df, gene_var > q99_var | gene_mean > q99_mean),  # Only label the top q99 genes
      aes(label = gene),
      size=4
    )
  mean_var_plot
}

log_tom_mean_var_plot = mean_var_plot(log(tpm_data+1))
log_tom_mean_var_plot
ggsave(log_tom_mean_var_plot,file="RNA_Analysis/outputs/plots/00_qc_plots/01_filtered_log_tpm_mean_var.png")

vst_mean_var_plot = mean_var_plot(vst_data)
vst_mean_var_plot
ggsave(log_tom_mean_var_plot,file="RNA_Analysis/outputs/plots/00_qc_plots/01_filtered_vst_mean_var.png")

# singscore_mean_var_plot = mean_var_plot(singscore_tpm_paired_cancer_data)
# ggsave(singscore_mean_var_plot,file="reference_data/RNA_Seq/outputs/plots/01_singscore_paired_cancer_mean_var.png")
# 
# gapdh_tpm = as.numeric(tpm_data["GAPDH",])
# tpm_data_gapdh = sweep(tpm_data, 2, gapdh_tpm, FUN="/")
# log_tpm_data_gapdh = log(tpm_data_gapdh+1)
# 
# gapdh_mean_var_plot = mean_var_plot(log_tpm_data_gapdh)

# 
# # Function to perform min-max rank normalization
# # from https://medium.com/@josef.waples/min-max-rank-normalization-in-r-using-mtcars-7d3327dc7937sa
# min_max_rank_normalize = function(x) {
#   ranked = rank(x)
#   normalized = (ranked - min(ranked)) / (max(ranked) - min(ranked))
#   return(normalized)
# }
# 
# # Compute for each gene the proportion of reads per sample library size
# get_gene_ranks = function(gene_data,normalize_by="colSums",rank_by="mean_prop",descending=TRUE) {
#   normalized_df = NA
#   if (normalize_by=="colSums") {
#     ## Normalize counts by library size
#     normalized_df = sweep(gene_data,2,colSums(count_data),`/`)
#   } else if (normalize_by=="rank") {
#     normalized_df = as.data.frame(lapply(gene_data, min_max_rank_normalize))
#     rownames(normalized_df) = rownames(gene_data)
#   } else if (normalize_by=="none") {
#     normalized_df = gene_data
#   }
#   mean_prop_across_samples = apply(normalized_df,1,mean)
#   median_across_samples = apply(normalized_df,1,median)
#   prop_var_across_samples = apply(normalized_df,1,var)
#   rank_df = data.frame(list(
#     "gene"=names(mean_prop_across_samples),
#     "mean_prop"=mean_prop_across_samples,
#     "median"=median_across_samples,
#     "prop_var"=prop_var_across_samples
#   ))
#   if(descending){
#     rank_df = rank_df[order(-rank_df[[rank_by]]),]
#   }else{
#     rank_df = rank_df[order(rank_df[[rank_by]]),]
#   }
#   rank_df$rank = c(1:dim(rank_df)[[1]])
#   return(rank_df)
# }
# 
# make_rank_plot = function(gene_data,normalize_by="colSums") {
#   ## Plot genes on average that takes up the highest proportion of reads per sample
#   gene_ranks_df = get_gene_ranks(gene_data,normalize_by=normalize_by)
#   gene_ranks_df$label = rownames(gene_ranks_df)
#   rank_plot = ggplot(gene_ranks_df, aes(x=rank, y = mean_prop)) +
#     geom_point(aes(color = rank <= 10), size = 1.2) +  # Color points differently if they are in the top 10
#     geom_text_repel(
#       data = subset(gene_ranks_df, rank <= 10),  # Only label the top 10 genes
#       aes(label = label),
#       size=4
#     ) +
#     #ylim(c(0,0.08)) +
#     theme_minimal() +
#     theme(plot.background = element_rect(fill = "white")) +
#     xlab("Rank order") +
#     ylab("Mean read proportion across samples")
#   rank_plot
# }
# 
# ## Plot genes on average that takes up the highest proportion of reads per sample (count version)
# count_rank_plot = make_rank_plot(count_data)
# count_rank_plot
# ggsave("reference_data/RNA_Seq/outputs/plots/01_gene_prop_rank_plot.png",count_rank_plot)
# 
# ## Plot genes on average that takes up the highest proportion of reads per sample (tpm version)
# tpm_rank_plot = make_rank_plot(tpm_data)
# tpm_rank_plot
# ggsave("reference_data/RNA_Seq/outputs/plots/01_gene_prop_rank_plot_tpm.png",tpm_rank_plot)
# 
# ## Plot genes on average that takes up the highest proportion of reads per sample (vst version)
# vst_rank_plot = make_rank_plot(vst_data)
# vst_rank_plot
# ggsave("reference_data/RNA_Seq/outputs/plots/01_gene_prop_rank_plot_vst.png",vst_rank_plot)
# 
# ## Plot TPM Box by Sample
# tpm_data_long = tpm_data
# tpm_data_long$gene = rownames(tpm_data_long)
# tpm_data_long = pivot_longer(tpm_data_long,cols = -gene,names_to = "Sample",values_to = "Value") 
# tpm_data_long_merge = merge(tpm_data_long,meta_data[,c("id","BATCH")],all.x=TRUE,by.x="Sample",by.y="id")
# tpm_data_long_merge$Value = log(tpm_data_long_merge$Value + 1)
# 
# ## Sort by Batch
# tpm_data_long_merge = tpm_data_long_merge %>% arrange(BATCH, Sample)
# tpm_data_long_merge$Sample = factor(tpm_data_long_merge$Sample, levels = unique(tpm_data_long_merge$Sample))
#                       
# log_tpm_boxplot = ggplot(tpm_data_long_merge,aes(x=Sample,y=Value,color=BATCH)) + 
#   geom_boxplot(outlier.shape = NA) +
#   theme(axis.text.x=element_blank()) +
#   ylab("Log2 (TPM + 1)") +
#   xlab("Sample")
# 
# log_tpm_boxplot
# ggsave("reference_data/RNA_Seq/outputs/plots/01_sample_by_log_tpm_boxplot.png",log_tpm_boxplot,width=14)
# 
# 
# make_gsva_rank_plot = function(gsva_data,return_data=FALSE,rank_by="median",top_viz=10,normalize_by="rank") {
#   # Make rank normalized plot for GSVA-related data
#   if(rank_by=="mean"){
#     rank_by="mean_prop"
#   }
#   gsva_rank_df = get_gene_ranks(gsva_data,normalize_by=normalize_by,rank_by=rank_by,descending = FALSE)
#   gsva_rank_df$label = rownames(gsva_rank_df)
#   gsva_rank_df[["rank_by_inverse"]] = 1 - gsva_rank_df[[rank_by]]
#   gsva_rank_df_plot = ggplot(gsva_rank_df, aes_string(x="rank",y="rank_by_inverse")) +
#     geom_point(aes(color = (rank <= top_viz)|(rank >=  dim(gsva_rank_df)[[1]]-top_viz)), size = 1.2) +
#     geom_text_repel(
#       data = subset(gsva_rank_df, rank <= top_viz),  # Only label the top 10 genes
#       aes(label = label),
#       size=3
#     ) + 
#     geom_text_repel(
#       data = subset(gsva_rank_df, rank >=  dim(gsva_rank_df)[[1]]-top_viz),  # Only label the top 10 genes
#       aes(label = label),
#       size=3
#     ) + 
#     ylab(paste0("1 - (",rank_by," rank)")) +
#     #labs(fill='Top/Bottom 20 Ranks') 
#     guides(color=guide_legend(title="Top/Bottom 20 Ranks"))
#   if(return_data){
#     return(gsva_rank_df)
#   } else{
#     return(gsva_rank_df_plot)
#   }
# }
# 
# make_mean_var_plot = function(rank_df) {
#   ## Todo: get mean ratio between cutaneous vs. non-cutaneous
#   # Given a ranking of the pathway/genes, plot the mean vs. variance relationship
#   rank_df_mean_var_plot = ggplot(rank_df,aes(y=mean_prop,x=prop_var,color=mean_prop)) +
#     geom_point() +
#     geom_text_repel(
#       data = rank_df,  # Only label the top 10 genes
#       aes(label = gene),
#       size=3
#     ) +
#     scale_colour_gradient2(mid = "purple", midpoint = 0.5)
#   rank_df_mean_var_plot
# }
# 
# ## Plot GSVA Ranks
# log_tpm_go_bp_rank_plot = make_gsva_rank_plot(gsva_log_tpm_go_bp_data)
# ggsave("reference_data/RNA_Seq/outputs/plots/01_log_tpm_go_bp_rank_plot_median.png",log_tpm_go_bp_rank_plot,width=14)
# 
# log_tpm_cancer_rank_plot = make_gsva_rank_plot(gsva_log_tpm_cancer_data)
# ggsave("reference_data/RNA_Seq/outputs/plots/01_log_tpm_cancer_rank_plot_median.png",log_tpm_cancer_rank_plot,width=14)
# 
# gsva_vst_go_bp_rank_plot = make_gsva_rank_plot(gsva_vst_go_bp_data)
# gsva_vst_go_bp_rank_plot
# ggsave("reference_data/RNA_Seq/outputs/plots/01_vst_go_bp_rank_plot_median.png",gsva_vst_go_bp_rank_plot,width=14)
# 
# gsva_vst_cancer_rank_plot = make_gsva_rank_plot(gsva_vst_cancer_data)
# gsva_vst_cancer_rank_plot
# ggsave("reference_data/RNA_Seq/outputs/plots/01_vst_cancer_rank_plot_median.png",gsva_vst_cancer_rank_plot,width=14)
# 
# ## Plot singscore Ranks
# tpm_go_bp_singscore_rank_plot = make_gsva_rank_plot(singscore_tpm_go_bp_data,return_data=TRUE,normalize_by="none")
# tpm_go_bp_singscore_rank_plot
# ggsave("reference_data/RNA_Seq/outputs/plots/01_singscore_go_bp_rank_plot_median.png",tpm_go_bp_singscore_rank_plot,width=14)
# 
# tpm_paired_cancer_singscore_rank_df = make_gsva_rank_plot(singscore_tpm_paired_cancer_data,return_data=TRUE,normalize_by="none")
# tpm_paired_cancer_singscore_rank_plot = ggplot(head(tpm_paired_cancer_singscore_rank_df,20),
#                                                aes(y=reorder(label,-rank),x=order(mean_prop))) + geom_col()
# tpm_paired_cancer_singscore_rank_plot
# ggsave("reference_data/RNA_Seq/outputs/plots/01_singscore_paired_cancer_rank_plot_median.png",tpm_paired_cancer_singscore_rank_plot,width=14)
# 
# ## Plot the mean-variance relationships
# singscore_mean_var_tpm_paired_cancer_plot = make_mean_var_plot(tpm_paired_cancer_singscore_rank_df,normalize_by="c")
# singscore_mean_var_tpm_paired_cancer_plot
# ggsave("reference_data/RNA_Seq/outputs/plots/01_singscore_mean_var_tpm_paired_cancer_plot.png",singscore_mean_var_tpm_paired_cancer_plot,width=14)
# 
# 
# tpm_bidirection_cancer_singscore_rank_plot = make_gsva_rank_plot(singscore_tpm_bidirection_cancer_data)
# tpm_bidirection_cancer_singscore_rank_plot
# ggsave("reference_data/RNA_Seq/outputs/plots/01_singscore_paired_cancer_rank_plot_median.png",tpm_go_bp_singscore_rank_plot,width=14)
# 
# tpm_bidirection_hallmark_singscore_rank_plot = make_gsva_rank_plot(singscore_tpm_hallmark_data)
# tpm_bidirection_hallmark_singscore_rank_plot
# ggsave("reference_data/RNA_Seq/outputs/plots/01_singscore_paired_cancer_rank_plot_median.png",tpm_go_bp_singscore_rank_plot,width=14)
# 
# tpm_paired_cgp_singscore_rank_df = make_gsva_rank_plot(singscore_tpm_paired_cgp_data,return_data = TRUE,rank_by="median")
# tpm_paired_cgp_singscore_rank_plot = ggplot(head(tpm_paired_cgp_singscore_rank_df,20),
#                                                aes(y=reorder(label,-rank),x=order(median))) + geom_col()
# tpm_paired_cgp_singscore_rank_plot
# ggsave("reference_data/RNA_Seq/outputs/plots/01_singscore_paired_cgp_rank_plot_median.png",tpm_paired_cgp_singscore_rank_plot,width=14)
# 
# 
# ## Focus on the ranks by POT1-biallelic samples
# pot1_biallelic_samples = meta_data[(meta_data$POT1_Mut_Germline=="Detected")&(meta_data$POT1_Mut!="Undetected")&(meta_data$POT1_Mut!="Not Profiled"),]
# pot1_biallelic_samples = pot1_biallelic_samples[!duplicated(pot1_biallelic_samples$individual_alias),"id"]
# gsva_vst_go_bp_pot1 = gsva_vst_go_bp_data[,pot1_biallelic_samples]
# pot1_gsva_plot = make_gsva_rank_plot(gsva_vst_go_bp_pot1,rank_by="mean")
# pot1_gsva_rank_df = make_gsva_rank_plot(gsva_vst_go_bp_pot1,rank_by="mean",return_data = TRUE)
# pot1_gsva_plot
# ggsave("reference_data/RNA_Seq/outputs/plots/01_pot1_biallelic_samples_gsva_rank_plot.png",pot1_gsva_plot,width=14)
# 
# 
# ## Do a site-specific rank plot
# breast_cutaneous_samples = meta_data[meta_data$primary_site_combined=="BREAST (CUTANEOUS)","id"]
# gsva_vst_go_bp_breast_cut = gsva_vst_go_bp_data[,breast_cutaneous_samples]
# gsva_vst_cancer_breast_cut = gsva_vst_cancer_data[,breast_cutaneous_samples]
# breast_cut_gsva_go_bp_plot = make_gsva_rank_plot(gsva_vst_go_bp_breast_cut,rank_by="mean_prop")
# breast_cut_gsva_cancer_plot = make_gsva_rank_plot(gsva_vst_cancer_breast_cut,rank_by="mean_prop")
# 
# breast_cut_gsva_go_bp_plot
# breast_cut_gsva_cancer_plot
# 
# breast_non_cutaneous_samples = meta_data[meta_data$primary_site_combined=="BREAST (PARENCHYMAL)","id"]
# gsva_vst_go_bp_breast_noncut = gsva_vst_go_bp_data[,breast_non_cutaneous_samples]
# gsva_vst_cancer_breast_noncut = gsva_vst_cancer_data[,breast_non_cutaneous_samples]
# breast_noncut_gsva_go_bp_plot = make_gsva_rank_plot(gsva_vst_go_bp_breast_noncut,rank_by="mean_prop")
# breast_noncut_gsva_cancer_plot = make_gsva_rank_plot(gsva_vst_cancer_breast_noncut,rank_by="mean_prop")
# 
# breast_noncut_gsva_go_bp_plot
# breast_cut_gsva_cancer_plot

