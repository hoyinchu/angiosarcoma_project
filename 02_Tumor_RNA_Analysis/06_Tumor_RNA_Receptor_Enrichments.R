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





# 
# plot_receptor_families = function(expr_df_subset_merged) {
#   # Discretize significance levels and assign colors by fold change direction
#   # Categorize significance and create expr_df_subset_merged symbols
#   print(hclust_expr_df_subset_merged)
#   expr_df_subset_merged = expr_df_subset_merged %>%
#     mutate(
#       significance = case_when(
#         padj < 0.001 ~ "Adj p < 0.001",
#         padj < 0.01  ~ "Adj p < 0.01",
#         padj < 0.05  ~ "Adj p < 0.05",
#         TRUE         ~ "Non-significant"
#       ),
#       annotation = case_when(
#         padj < 0.001 ~ "***",
#         padj < 0.01  ~ "**",
#         padj < 0.05  ~ "*",
#         TRUE         ~ ""
#       ),
#       color_intensity = abs(log2FoldChange),  # Scale color intensity by absolute log2FoldChange
#       color_category = case_when(
#         log2FoldChange > 0 ~ "red",
#         log2FoldChange < 0 ~ "blue",
#         TRUE ~ "gray70"
#       )
#     )# %>%
#     #filter(significance != "Non-significant")  # Remove non-significant genes
#   
#   # Define color mapping for intensity scaling
#   color_scale = scale_color_gradient2(
#     low = "blue", high = "red", mid = "gray70", midpoint = 0,
#     limits = c(-4, 4),
#     oob = scales::squish,
#     name = "Log2 Fold Change"
#   )
#   # Create dot plot with family-based faceting and significance annotation
#   dotplot_by_family = ggplot(expr_df_subset_merged, 
#                               aes(y = cluster, x = gene, 
#                                   size = abs(log2FoldChange),  # Size reflects fold change
#                                   color = pmax(pmin(log2FoldChange, 4), -4))) +  # Color intensity based on fold change
#     geom_point() +
#     geom_text(aes(label = annotation), vjust = 0.9, size = 5, color = "black", fontface = "bold") +  # Add significance annotations
#     scale_size(range = c(2, 10)) +  # Adjust dot size range
#     color_scale +  # Apply color gradient
#     facet_grid(. ~ Pathway, scales = "free_x", space = "free_x") +  # Group by family
#     theme_minimal() +
#     theme(strip.text.y = element_text(angle = 0, face = "bold"),  # Make family labels readable
#           panel.spacing = unit(1, "lines"),
#           axis.text.x = element_text(angle = 45, hjust = 1)) +
#     labs(title = "",
#          x = "Gene", 
#          y = "Cluster", 
#          size = "Fold Change (abs)")
#   
#   return(dotplot_by_family)
# }
# 
# make_merged_expr_df = function(expr_df,annot_df) {
#   ## Merge the dataframe with expression information
#   receptor_genes_to_subset = genes_with_targeted_therapy$Gene.Symbol
#   expr_df_subset = expr_df[expr_df$gene %in% receptor_genes_to_subset,]
#   expr_df_subset_merged = merge(expr_df_subset,annot_df,by.x="gene",by.y="Gene.Symbol")
#   return(expr_df_subset_merged)
# }
# 
# hclust_expr_df_subset_merged = make_merged_expr_df(all_hclust_deseq_res_all,genes_with_targeted_therapy)
# hclust_expr_df_subset_merged$cluster = factor(hclust_expr_df_subset_merged$cluster, levels = rev(c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5")))
# # hclust_expr_df_subset_merged$Pathway = factor(hclust_expr_df_subset_merged$Pathway, levels = c(
# #   "VEGF Signaling","PDGF Signaling","FGF Signaling","Angiopoietin-TIE Signaling","HGF Signaling","Neuropilin","Immunotherapy"
# # ))
# 
# hclust_receptor_plot = plot_receptor_families(hclust_expr_df_subset_merged)
# hclust_receptor_plot
# #ggsave("02_Tumor_RNA_Analysis/outputs/plots/06_receptor_plots/receptor_enrichments_by_cluster.png",hclust_receptor_plot,dpi=300,height=10,width=8)
# ggsave("02_Tumor_RNA_Analysis/outputs/plots/06_receptor_plots/receptor_enrichments_by_cluster.png",hclust_receptor_plot,dpi=300,height=4,width=10)
# ggsave("02_Tumor_RNA_Analysis/outputs/plots/06_receptor_plots/receptor_enrichments_by_cluster.pdf",hclust_receptor_plot,dpi=300,height=4,width=10)
# 
# 
# ## Do the same for cluster by sites
# site_expr_df_subset_merged = make_merged_expr_df(all_deseq_res_all,genes_with_targeted_therapy)
# site_expr_df_subset_merged$cluster = factor(site_expr_df_subset_merged$cluster, levels = rev(c(
#   "ParenchymalBreast","CutaneousBreast","HNFS","Heart","Extremities","Others"
#   )))
# # site_expr_df_subset_merged$Pathway = factor(site_expr_df_subset_merged$Pathway, levels = c(
# #   "VEGF Signaling","PDGF Signaling","FGF Signaling","Angiopoietin-TIE Signaling","HGF Signaling","Neuropilin","Immunotherapy"
# # ))
# site_receptor_plot = plot_receptor_families(site_expr_df_subset_merged)
# site_receptor_plot
# #ggsave("02_Tumor_RNA_Analysis/outputs/plots/06_receptor_plots/receptor_enrichments_by_site.png",site_receptor_plot,dpi=300,height=10,width=8)
# ggsave("02_Tumor_RNA_Analysis/outputs/plots/06_receptor_plots/receptor_enrichments_by_site.png",site_receptor_plot,dpi=300,height=4,width=10)
# ggsave("02_Tumor_RNA_Analysis/outputs/plots/06_receptor_plots/receptor_enrichments_by_site.pdf",site_receptor_plot,dpi=300,height=4,width=10)


