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


## Load the palettes
source("../util_scripts/project_palettes.R")

## Load the hierarchical clusters
so = readRDS("../data/processed/rna/ASCSeuratObj2025.rds")
hclust_res_pos = read.csv("./outputs/DEGs/02_DESEQ2_hclust_expr_sites_combined_results_pos.tsv", sep = "\t")
hclust_combined_fora = read.csv("./outputs/DEGs/hclust_deg_hallmark_c6_combined.csv")

## Decide on what genes to show for the heatmap by
## getting the top gene set and show the top overexpressed genes
## If genes is already shown in another set skip it

## hclust1
hclust1_hallmark = hclust_combined_fora %>% 
  filter(padj < 0.01, term_set == "Hallmarks",cluster == "Cluster 1") %>%
  #arrange(-odds_ratio) %>%
  head(1)
hclust1_set_hallmark_genes = unlist(strsplit(hclust1_hallmark %>% pull(overlapGenes), ", "))
hclust1_top_hallmark_genes = hclust_res_pos %>% 
  filter(cluster=="Cluster 1", gene %in% hclust1_set_hallmark_genes) %>%
  arrange(padj) %>% head(5) %>% arrange(-log2FoldChange)
hclust1_top_hallmark_genes

hclust1_set_c6 = hclust_combined_fora %>% 
  filter(padj < 0.01, term_set == "C6",cluster == "Cluster 1") %>%
  #arrange(-odds_ratio) %>%
  head(1)
hclust1_set_c6_genes = unlist(strsplit(hclust1_set_c6 %>% pull(overlapGenes), ", "))
hclust1_top_c6_genes = hclust_res_pos %>% 
  filter(cluster=="Cluster 1", gene %in% hclust1_set_c6_genes) %>%
  arrange(padj) %>% head(5) %>% arrange(-log2FoldChange)
hclust1_top_c6_genes


## hclust2
hclust2_hallmark = hclust_combined_fora %>% 
  filter(padj < 0.01, term_set == "Hallmarks",cluster == "Cluster 2") %>%
  #arrange(-odds_ratio) %>%
  head(1)
hclust2_set_hallmark_genes = unlist(strsplit(hclust2_hallmark %>% pull(overlapGenes), ", "))
hclust2_top_hallmark_genes = hclust_res_pos %>% 
  filter(cluster=="Cluster 2", gene %in% hclust2_set_hallmark_genes) %>%
  arrange(padj) %>% head(5) %>% arrange(-log2FoldChange)
hclust2_top_hallmark_genes

hclust2_set_c6 = hclust_combined_fora %>% 
  filter(padj < 0.01, term_set == "C6",cluster == "Cluster 2") %>%
  #arrange(-odds_ratio) %>%
  head(1)
hclust2_set_c6_genes = unlist(strsplit(hclust2_set_c6 %>% pull(overlapGenes), ", "))
hclust2_top_c6_genes = hclust_res_pos %>% 
  filter(cluster=="Cluster 2", gene %in% hclust2_set_c6_genes) %>%
  arrange(padj) %>% head(5) %>% arrange(-log2FoldChange)
hclust2_top_c6_genes

## hclust3
hclust3_hallmark = hclust_combined_fora %>% 
  filter(padj < 0.01, term_set == "Hallmarks",cluster == "Cluster 3") %>%
  #arrange(-odds_ratio) %>%
  head(1)
hclust3_set_hallmark_genes = unlist(strsplit(hclust3_hallmark %>% pull(overlapGenes), ", "))
hclust3_top_hallmark_genes = hclust_res_pos %>% 
  filter(cluster=="Cluster 3", gene %in% hclust3_set_hallmark_genes) %>%
  arrange(padj) %>% head(5) %>% arrange(-log2FoldChange)
hclust3_top_hallmark_genes

hclust3_set_c6 = hclust_combined_fora %>% 
  filter(padj < 0.01, term_set == "C6",cluster == "Cluster 3") %>%
  #arrange(-odds_ratio) %>%
  head(1)
hclust3_set_c6_genes = unlist(strsplit(hclust3_set_c6 %>% pull(overlapGenes), ", "))
hclust3_top_c6_genes = hclust_res_pos %>% 
  filter(cluster=="Cluster 3", gene %in% hclust3_set_c6_genes) %>%
  arrange(padj) %>% head(5) %>% arrange(-log2FoldChange)
hclust3_top_c6_genes

## hclust4
hclust4_hallmark = hclust_combined_fora %>% 
  filter(padj < 0.01, term_set == "Hallmarks",cluster == "Cluster 4") %>%
  head(1)
hclust4_set_hallmark_genes = unlist(strsplit(hclust4_hallmark %>% pull(overlapGenes), ", "))
hclust4_top_hallmark_genes = hclust_res_pos %>% 
  filter(cluster=="Cluster 4", gene %in% hclust4_set_hallmark_genes) %>%
  arrange(padj) %>% head(5) %>% arrange(-log2FoldChange)
hclust4_top_hallmark_genes

hclust4_set_c6 = hclust_combined_fora %>% 
  filter(padj < 0.01, term_set == "C6",cluster == "Cluster 4") %>%
  head(1)
hclust4_set_c6_genes = unlist(strsplit(hclust4_set_c6 %>% pull(overlapGenes), ", "))
hclust4_top_c6_genes = hclust_res_pos %>% 
  filter(cluster=="Cluster 4", gene %in% hclust4_set_c6_genes) %>%
  arrange(padj) %>% head(5) %>% arrange(-log2FoldChange)
hclust4_top_c6_genes


## hclust5
hclust5_hallmark = hclust_combined_fora %>% 
  filter(padj < 0.01, term_set == "Hallmarks",cluster == "Cluster 5") %>%
  head(2) %>% ## ESR Late has already been picked by an earlier cluster so we move to the next best
  tail(1)
hclust5_set_hallmark_genes = unlist(strsplit(hclust5_hallmark %>% pull(overlapGenes), ", "))
hclust5_top_hallmark_genes = hclust_res_pos %>% 
  filter(cluster=="Cluster 5", gene %in% hclust5_set_hallmark_genes) %>%
  arrange(padj) %>% head(5) %>% arrange(-log2FoldChange)
hclust5_top_hallmark_genes

hclust5_set_c6 = hclust_combined_fora %>% 
  filter(padj < 0.05, term_set == "C6",cluster == "Cluster 5") %>%
  head(1)
hclust5_set_c6_genes = unlist(strsplit(hclust5_set_c6 %>% pull(overlapGenes), ", "))
hclust5_top_c6_genes = hclust_res_pos %>% 
  filter(cluster=="Cluster 5", gene %in% hclust5_set_c6_genes) %>%
  arrange(padj) %>% head(5) %>% arrange(-log2FoldChange)
hclust5_top_c6_genes


## Combine all into the gene_set_df
gene_set_df = rbind(
  # hclust1
  data.frame(gene = hclust1_top_hallmark_genes$gene, pathway = hclust1_hallmark$pretty_pathway[[1]]),
  data.frame(gene = hclust1_top_c6_genes$gene, pathway = hclust1_set_c6$pretty_pathway[[1]]),
  # Cut
  data.frame(gene = hclust2_top_hallmark_genes$gene, pathway = hclust2_hallmark$pretty_pathway[[1]]),
  data.frame(gene = hclust2_top_c6_genes$gene, pathway = hclust2_set_c6$pretty_pathway[[1]]),
  # HNFS
  data.frame(gene = hclust3_top_hallmark_genes$gene, pathway = hclust3_hallmark$pretty_pathway[[1]]),
  data.frame(gene = hclust3_top_c6_genes$gene, pathway = hclust3_set_c6$pretty_pathway[[1]]),
  # Heart
  data.frame(gene = hclust4_top_hallmark_genes$gene, pathway = hclust4_hallmark$pretty_pathway[[1]]),
  data.frame(gene = hclust4_top_c6_genes$gene, pathway = hclust4_set_c6$pretty_pathway[[1]]),
  # Extremities
  data.frame(gene = hclust5_top_hallmark_genes$gene, pathway = hclust5_hallmark$pretty_pathway[[1]]),
  data.frame(gene = hclust5_top_c6_genes$gene, pathway = hclust5_set_c6$pretty_pathway[[1]])
)
gene_set_df


rownames(gene_set_df) = gene_set_df$gene
gene_set_df$`Gene Set` = gene_set_df$pathway
gene_set_df$pathway = NULL
gene_set_df$gene = NULL

gene_set_df$`Gene Set` = gsub("\\.V1", "", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("EPITHELIAL MESENCHYMAL TRANSITION", "EMT", gene_set_df$`Gene Set`)
#gene_set_df$`Gene Set` = gsub("TGFB", "TGFB UP", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("MYC TARGETS V1", "MYC_TARGETS", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("MYC_TARGETS", "MYC", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("HEME METABOLISM", "HEME META", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("ESTROGEN RESPONSE LATE", "ESTROGEN", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("KRAS SIGNALING DN", "KRAS DN", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("ESC J1 EARLY", "ESC EARLY", gene_set_df$`Gene Set`)
#gene_set_df$`Gene Set` = gsub("ATF2 S", "ATF2 UP", gene_set_df$`Gene Set`)
gene_set_df$`Gene Set` = gsub("\\_", " ", gene_set_df$`Gene Set`)

gene_set_df

## If using the hallmark / C6 mixed gene set
gene_set_df$`Gene Set` = factor(
  gene_set_df$`Gene Set`,
  levels = c(
    "EMT","TGFB",
    "HEME META","EGFR",
    "ESTROGEN","ESC EARLY",
    "MYC","CSR LATE",
    "KRAS DN","P53 DN"
  )
)

## Define color palette for gene sets
total_sets = length(unique(gene_set_df$`Gene Set`))
#dynamic_pal_npg = pal_npg("nrc")(total_sets)
npg_colors = pal_npg("nrc")(10) 
dynamic_pal_npg = colorRampPalette(npg_colors)(total_sets)
gene_set_palette = setNames(dynamic_pal_npg, levels(gene_set_df$`Gene Set`))

## Define color palette for the experiemtnal hybrid set
# Get unique pathways
#unique_pathways = unique(gene_set_df$`Gene Set`)
#npg_palette = pal_npg("nrc")(length(unique_pathways))
#gene_set_palette = setNames(npg_palette, unique_pathways)


## Sort the idendities,sites,and heatmap data
asc_cluster_idents = so@meta.data %>% arrange(hcluster_by_expr_str,`entity:sample_id`) %>% dplyr::select(hcluster_by_expr_str)
asc_cluster_idents$hcluster_by_expr_str = factor(
  asc_cluster_idents$hcluster_by_expr_str,
  levels = c(
    "Cluster 1",
    "Cluster 2",
    "Cluster 3",
    "Cluster 4",
    "Cluster 5"
  )
)

# This time sort by hclust
asc_primary_sites =  so@meta.data %>% arrange(hcluster_by_expr_str,`entity:sample_id`) %>% dplyr::select(`Primary Site (Recombined)`)

asc_new_marker_heatmap_data = GetAssayData(so,layer="vst_scaled")[rownames(gene_set_df),rownames(asc_cluster_idents)]

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
all_markers_dedupped = hclust_res_pos %>% arrange(-log2FoldChange) %>% distinct(gene, .keep_all = TRUE)
all_new_markers_FC = all_markers_dedupped %>% filter(gene %in% rownames(gene_set_df))
rownames(all_new_markers_FC) = all_new_markers_FC$gene
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
# This annotation require inputs from running other scripts so if it is your first time 
# running this block you should come back after running all other scripts first
# Specifically, skip the rest of the RNA scripts and run things in the 
# 03_Tumor_WES_Analsysis folder then come back

rep_mut_df = read.csv("../03_Tumor_WES_Analysis/outputs/tables/representative_mutation_table.csv",check.names = FALSE)
rep_mut_df = rep_mut_df[,c("sample_alias_cleaned","tumor_mutation_id","tumor_mutation_id_short","Hugo_Symbol")] %>%
  dplyr::rename(
    rep_tumor_mutation_id = "tumor_mutation_id",
    rep_tumor_mutation_id_short = "tumor_mutation_id_short",
    rep_Hugo_Symbol = "Hugo_Symbol"
  )
representative_somatic_muts_merged = merge(so@meta.data,rep_mut_df,by="sample_alias_cleaned",all.x=TRUE)

## Also add in sample TMB
tmb_df = read.csv("../03_Tumor_WES_Analysis/outputs/tables/total_mutations_per_sample.csv")
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

germline_rep_mut_df = read.csv("../05_Germline_WES_Analysis/outputs/tables/representative_germline_pv_table.csv",check.names=FALSE)
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
  arrange(hcluster_by_expr_str,`entity:sample_id`) %>% 
  replace_na(list(rep_Hugo_Symbol="Not available"))

representative_germline_vars_merged = representative_somatic_muts_merged_with_germline %>%
  arrange(hcluster_by_expr_str,`entity:sample_id`) %>% 
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
  "Expression Cluster" = representative_somatic_muts_merged$hcluster_by_expr_str,
  `Primary Site`= asc_primary_sites$`Primary Site (Recombined)`,
  `Cutaneous` = representative_somatic_muts_merged$`cutaneous_viz`,
  "RAAS/LAAS" = representative_somatic_muts_merged$RAAS_LAAS_Class,
  `Age`= representative_somatic_muts_merged$`Age (Combined)`,
  "Sex" = representative_somatic_muts_merged$`SEX (EHR_EXTRACTED)`,
  "Sample Purity" = representative_somatic_muts_merged$ESTIMATE_purity,
  "Epithelioid" = representative_somatic_muts_merged$is_epithelioid,
  "Spindle Cell" = representative_somatic_muts_merged$is_spindle,
  "Nuclear Grade" = representative_somatic_muts_merged$BX_NUCLEAR_GRADE,
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
pdf("./outputs/plots/02_seurat_plots/02_cluster_gene_set_heatmap_dendro_cluster_by_hclust_top_5_pathway_test.pdf",width=24,height=20)
#png("./outputs/plots/02_seurat_plots/02_cluster_gene_set_heatmap_dendro_cluster_by_site_top_5_pathway_test.png",width=24,height=18)

draw(complex_heatmap_fig, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()

