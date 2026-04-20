library(dplyr)
library(tidyverse)

setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts")

germline_df = read.csv("../reference_data/ClinicalTables/Jun2023_ASC_Case_Control_Sample_Merged_Germline.tsv",sep="\t",check.names=FALSE)

## Convert disease status to int and factors
germline_df$is_ASC = recode(as.factor(germline_df$is_ASC),"True"=1,"False"=0)
germline_df$inferred_sex = recode(as.factor(germline_df$inferred_sex),"F"=0,"M"=1)

## Subset to gene columns
case_germline_df = germline_df %>% filter(is_ASC==1)
case_germline_cols = colnames(case_germline_df)[grepl( "Germline_" , colnames(case_germline_df))]

# Convert to long format
case_germline_df_subset = case_germline_df[,c("Sample",case_germline_cols)]
case_germline_df_subset_long = case_germline_df_subset %>%
  pivot_longer(cols=-c("Sample"),names_to = "Hugo_Symbol", values_to = "PV_count") %>%
  filter(PV_count > 0)
case_germline_df_subset_long$Hugo_Symbol = sapply(str_split(case_germline_df_subset_long$Hugo_Symbol, "_"), function(x) x[2])
case_germline_df_subset_long

# Merge information about each gene
# Check p-value from burden test
burden_result = read.csv("05_Germline_WES_Analysis/outputs/tables/germline_PV_gene_enrichment_min_2_in_case.tsv",sep="\t")
burden_result = burden_result %>% arrange(p.value)
burden_result$Hugo_Symbol = sapply(str_split(burden_result$term, "_"), function(x) x[2])
burden_significant_genes = burden_result[burden_result$adjusted_pval < 0.1,]$Hugo_Symbol
burden_nominal_genes = burden_result[burden_result$p.value < 0.05,]$Hugo_Symbol

# Check if it is a COSMIC tier 1 gene
cosmic_path = "data/public/cosmic_cancer_gene_census.csv"
cosmic_df = read.csv(cosmic_path)
cosmic_tier1_genes = cosmic_df[cosmic_df$Tier==1,]$Gene.Symbol
# Or if somatic pathogenic variant was also detected in the gene
somatic_germline_mut_table = read.csv("06_Germline_WES_Tumor_WES_Analysis/outputs/tables/germline_vs_tumor_gene_counts.csv")
somatic_gerline_mut_genes = somatic_germline_mut_table %>%
  filter(germline_pv_carrier_count > 0, somatic_nonsyn_carrier_count>0) %>%
  pull(Hugo_Symbol)
# Or if it is in a gene of interest
goi_list = c("CFTR")

case_germline_df_subset_long = case_germline_df_subset_long %>%
  mutate(
    is_burden_significant =  Hugo_Symbol %in% burden_significant_genes,
    is_burden_nominal = Hugo_Symbol %in% burden_nominal_genes,
    is_cosmic_tier1 = Hugo_Symbol %in% cosmic_tier1_genes,
    is_somatic_mutated = Hugo_Symbol %in% somatic_gerline_mut_genes
  )

# Add priority score
case_germline_df_subset_long = case_germline_df_subset_long %>%
  mutate(germline_gene_priority_score = (is_burden_significant*1) +
           (is_burden_nominal*1) +
           (is_cosmic_tier1*1) +
           (is_somatic_mutated+1)
  )
case_germline_df_subset_long

# Pick most prioritized gene per sample
case_germline_df_high_priority = case_germline_df_subset_long %>%
  arrange(-germline_gene_priority_score) %>%
  #filter(germline_gene_priority_score > 0) %>%
  distinct(Sample,.keep_all = TRUE)
dim(case_germline_df_high_priority)

## Add metadata back in
keep_meta_cols = c(
  "Sample",
  "PCA1","PCA2","PCA3","PCA4","PCA5","PCA6","PCA7","PCA8","PCA9","PCA10",
  "inferred_ancestry_PCA","inferred_sex","sample_alias","individual_alias",
  "STUDY ID"
)
case_germline_df_non_germline = case_germline_df[,keep_meta_cols]

case_germline_df_high_priority_merged = merge(
  case_germline_df_high_priority,
  case_germline_df_non_germline,
  by="Sample",all.x=TRUE
  )

write.csv(case_germline_df_high_priority_merged,
          "05_Germline_WES_Analysis/outputs/tables/representative_germline_pv_table.csv",row.names = FALSE)



