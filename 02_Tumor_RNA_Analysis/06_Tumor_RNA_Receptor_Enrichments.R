library(ggplot2)
library(dplyr)
library(tidyr)
library(BuenColors)


## Show expression of Receptor Tyrosine Kinase families Per Cluster
receptor_family_table = read.csv("../data/curated/receptor_family_tables_modified.csv")
genes_with_targeted_therapy = read.csv("../data/curated/genes_with_targeted_therapy.csv")

all_hclust_deseq_res_all = read.csv("./outputs/DEGs/02_DESEQ2_hclust_expr_sites_combined_results_all.tsv",sep="\t")
all_deseq_res_all = read.csv("./outputs/DEGs/02_DESEQ2_sites_combined_results_all.tsv",sep="\t")

## Subset to genes of interest
all_hclust_deseq_res_all_subset = all_hclust_deseq_res_all %>% filter(gene %in% genes_with_targeted_therapy$Gene.Symbol)
all_deseq_res_all_subset = all_deseq_res_all %>% filter(gene %in% genes_with_targeted_therapy$Gene.Symbol)

## Dotplot showing 
all_deseq_res_all_subset_plot_data = all_deseq_res_all_subset %>%
  mutate(stars = cut(padj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                     labels = c("***", "**", "*", ""))) %>%
  mutate(gene = factor(gene, levels = sort(unique(gene), decreasing = TRUE)))

# Generate the plot
all_deseq_res_all_subset_plot = ggplot(all_deseq_res_all_subset_plot_data , aes(x = cluster, y = gene)) +
  geom_point(aes(size = log10(baseMean+1), color = log2FoldChange)) +
  geom_text(aes(label = stars), vjust = 0.8, color = "black") + # Adds the stars
  scale_color_gradient2(low = "dodgerblue3", mid = "white", high = "firebrick") +
  labs(size = "BaseMean", color = "Log2FC") +
  pretty_plot() + L_border() +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
  theme(legend.position = "none") #+
  #theme(legend.position = "bottom")
all_deseq_res_all_subset_plot
ggsave("./outputs/plots/06_receptor_plots/all_deseq_res_all_subset_drug_target_plot_new.pdf",all_deseq_res_all_subset_plot,dpi=300,height=8.2,width=3.3)
ggsave("./outputs/plots/06_receptor_plots/all_deseq_res_all_subset_drug_target_plot_no_legend.pdf",all_deseq_res_all_subset_plot,dpi=300,height=8.2,width=3.3)

# Generate it for clusters
all_hclust_deseq_res_all_subset_plot_data = all_hclust_deseq_res_all_subset %>%
  mutate(stars = cut(padj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                     labels = c("***", "**", "*", ""))) %>%
  mutate(gene = factor(gene, levels = sort(unique(gene), decreasing = TRUE)))

all_hclust_deseq_res_all_subset_plot = ggplot(all_hclust_deseq_res_all_subset_plot_data , aes(x = cluster, y = gene)) +
  geom_point(aes(size = log10(baseMean+1), color = log2FoldChange)) +
  geom_text(aes(label = stars), vjust = 0.8, color = "black") + # Adds the stars
  scale_color_gradient2(low = "dodgerblue3", mid = "white", high = "firebrick") +
  labs(size = "BaseMean", color = "Log2FC") +
  pretty_plot() + L_border() +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) + 
  theme(legend.position = "none") #+
  #theme(legend.position = "bottom")
all_hclust_deseq_res_all_subset_plot
ggsave("./outputs/plots/06_receptor_plots/all_hclust_deseq_res_all_subset_plot_new.pdf",all_hclust_deseq_res_all_subset_plot,dpi=300,height=8.2,width=3.3)
ggsave("./outputs/plots/06_receptor_plots/all_hclust_deseq_res_all_subset_plot_no_text.pdf",all_hclust_deseq_res_all_subset_plot,dpi=300,height=8.2,width=3.3)


