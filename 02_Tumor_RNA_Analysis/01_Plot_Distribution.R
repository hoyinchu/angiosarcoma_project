library(DESeq2)
library(ggplot2)
library(ggrepel)
library(tidyr)

## Script for visualizing qc metrics

## Load filtered data 
count_data = read.csv("../data/processed/rna/00_filtered_gene_counts.csv",row.names = 1,check.names = FALSE)
vst_data = read.csv("../data/processed/rna/00_filtered_gene_counts.vst.csv",row.names = 1,check.names = FALSE)
tpm_data = read.csv("../data/processed/rna/00_filtered_gene_tpm.csv",row.names = 1,check.names = FALSE)
meta_data = read.csv("../data/processed/rna/00_filtered_sample_metadata.csv",check.names = FALSE)

## Plot the library sizes by samples
lib_size_plot = ggplot(meta_data, aes(x=library_size,fill=LC_BATCH)) + 
  geom_histogram() +
  scale_x_log10() +
  xlab("Library Size") +
  ylab("# of Samples")
lib_size_plot
#ggsave("RNA_Analysis/outputs/plots/00_qc_plots/01_library_size_histogram.png",lib_size_plot)

## Plot logged gene count distribution
total_gene_counts = rowSums(count_data)
total_gene_counts = data.frame("log_count"=log(total_gene_counts[order(-total_gene_counts)]))
gene_counts_plot = ggplot(total_gene_counts, aes(x=log_count)) + 
  geom_histogram() +
  xlab("Log Gene Read Counts") +
  ylab("# of Genes")
gene_counts_plot
#ggsave("RNA_Analysis/outputs/plots/00_qc_plots/01_gene_read_counts_histogram.png",gene_counts_plot)

## Plot VST transformed library size and gene size
vst_library_sizes = as.data.frame(colSums(vst_data))
colnames(vst_library_sizes) = "library_size"
vst_library_sizes_plot = ggplot(vst_library_sizes, aes(x=library_size)) + 
  geom_histogram() +
  xlab("VST-transformed Library Sizes") +
  ylab("# of Samples")
vst_library_sizes_plot
#ggsave("RNA_Analysis/outputs/plots/00_qc_plots/01_vst_transformed_library_size_histogram.png",vst_library_sizes_plot)

vst_gene_counts = as.data.frame(rowSums(vst_data))
colnames(vst_gene_counts) = "gene_counts"
vst_gene_counts_plot = ggplot(vst_gene_counts, aes(x=gene_counts)) + 
  geom_histogram() +
  xlab("VST-transformed Gene Counts") +
  ylab("# of Genes with Read Count")
vst_gene_counts_plot
#ggsave("RNA_Analysis/outputs/plots/00_qc_plots/01_vst_transformed_gene_counts_histogram.png",vst_gene_counts_plot)

## Plot mean-variance relationship between genes
mean_var_plot = function(gene_data) {
  gene_data_mean = apply(gene_data,1,mean)
  gene_data_var = apply(gene_data,1,var)
  q99_mean= quantile(gene_data_mean,probs=c(0.99))
  q99_var= quantile(gene_data_var,probs=c(0.99))
  mean_var_df = as.data.frame(list("gene_mean"=gene_data_mean,"gene_var"=gene_data_var))
  mean_var_df$gene = rownames(mean_var_df)
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
#ggsave(log_tom_mean_var_plot,file="RNA_Analysis/outputs/plots/00_qc_plots/01_filtered_log_tpm_mean_var.png")

vst_mean_var_plot = mean_var_plot(vst_data)
vst_mean_var_plot
#ggsave(log_tom_mean_var_plot,file="RNA_Analysis/outputs/plots/00_qc_plots/01_filtered_vst_mean_var.png")

