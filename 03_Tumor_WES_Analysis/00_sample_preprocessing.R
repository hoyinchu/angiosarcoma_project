library(dplyr)
library(tidyr)
## This scripts identify all the samples that passed/failed quality control

# Load all tumor normal pairs available
pair_df = read.csv("../data/raw/sample_tables/pair_Jan13_2024.tsv",sep="\t",check.names = FALSE)
pair_set_df = read.csv("../data/raw/sample_tables/pair_set_membership.tsv",sep="\t",check.names = FALSE)

# Format names
pair_df["case_sample"] = sub('.*entityName:([^,}]*)[},].*', '\\1', pair_df$case_sample)
pair_df["control_sample"] = sub('.*entityName:([^,}]*)[},].*', '\\1', pair_df$control_sample)

# A single sample's name formatted incorrectly
pair_df["case_sample"] = gsub("__", "_", pair_df$case_sample)

# Identify which bait was used for each tumor-normal pair
illumina_pairs = pair_set_df %>% filter(`membership:pair_set_id` == "whole_exome_illumina_coding_v1_bait_set")
twist_pairs = pair_set_df %>% filter(`membership:pair_set_id` == "broad_custom_exome_v1_bait_set")
pair_df[["is_illumina"]] = pair_df[["entity:pair_id"]] %in% illumina_pairs$pair
pair_df[["is_twist"]] = pair_df[["entity:pair_id"]] %in% twist_pairs$pair


# Apply FFPE filter
illumina_ffpe_qval_threshold = 30
twist_ffpe_qval_threshold = 35
pair_df[["pass_ffpe_filter"]] = (pair_df[["is_illumina"]] & (pair_df[["ffpe_OBF_q_val"]] > illumina_ffpe_qval_threshold)) | (pair_df[["is_twist"]] & (pair_df[["ffpe_OBF_q_val"]] > twist_ffpe_qval_threshold)) 

## Apply Other quality filters
oxog_qval_threshold = 30
frac_contam_threshold = 0.04
tumor_in_normal_threshold = 0.3
pair_df[["pass_oxog_filter"]] = pair_df[["oxoG_OBF_q_val"]] > oxog_qval_threshold
pair_df[["pass_frac_contam_filter"]] = pair_df[["fracContam"]] < frac_contam_threshold
pair_df[["pass_TiN_filter"]] = pair_df[["TiN"]] <= tumor_in_normal_threshold

# Samples with missing information failed pipeline for various technical factors
pair_df[is.na(pair_df)] = FALSE

# If a patient has multiple normal samples, select the one used in germline analysis
pair_df[["is_in_germline_analysis"]] = pair_df[["in_germline_analysis"]]=="true"

# Check the number of tumor-normal pairs that passed all filters
passed_all_data_filters = pair_df[["pass_ffpe_filter"]] & pair_df[["pass_oxog_filter"]] & pair_df[["pass_frac_contam_filter"]] & pair_df[["pass_TiN_filter"]]
pair_df[["passed_all_data_filters"]] = passed_all_data_filters
pair_df[["pass_final_filter"]] = pair_df[["passed_all_data_filters"]] & pair_df[["is_in_germline_analysis"]]

## Write the output so we only grab these pairs from terra
pair_df_subset = pair_df[,c(
  "entity:pair_id","case_sample","control_sample",
  "is_illumina","is_twist",
  "pass_ffpe_filter","pass_oxog_filter","pass_frac_contam_filter","pass_TiN_filter","passed_all_data_filters",
  "is_in_germline_analysis","pass_final_filter"
)]

pair_df_subset["data_type"] = "tumor_WES"
write.table(pair_df_subset,"../data/processed/qc_tables/tumor_WES_pair_status.tsv",sep="\t",row.names = FALSE)

pair_df_subset_filtered = pair_df_subset %>% filter(pass_final_filter==TRUE)
write.csv(pair_df_subset_filtered,"../data/processed/qc_tables/all_filter_passed_pairs.csv")


