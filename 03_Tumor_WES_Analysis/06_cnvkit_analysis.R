library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(maftools)
library(ggsci)
library(ggbeeswarm)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

# Define the directory path
cnvkit_data_dir = "../data/processed/tumor_WES/cnvkit"

# Create the combined dataframe with cleaned sample IDs
combined_cnvkit <- list.files(path = cnvkit_data_dir, 
                              pattern = "\\.gainloss\\.txt$", 
                              full.names = TRUE) %>%
  set_names(basename(.) %>% str_remove("\\.gainloss\\.txt$")) %>% 
  map_dfr(~read_tsv(.x, col_types = cols(chromosome = col_character())), .id = "sample_id")

combined_cnvkit = combined_cnvkit %>%
  mutate(cnv_call = case_when(
    log2 < -1.5 ~ "DeepDel",
    log2 < -0.5 ~ "ShallowDel",
    log2 > 1.5  ~ "Amp",
    log2 > 0.5  ~ "ShallowAmp",
    TRUE        ~ "Neutral"
  ))


## Rename the pair names by tumor sample names
pair_df_subset_filtered = read.csv("../data/processed/qc_tables/all_filter_passed_pairs.csv")
combined_cnvkit_mapped = combined_cnvkit %>%
  left_join(
    pair_df_subset_filtered %>% 
      dplyr::select(entity.pair_id, case_sample), 
    by = c("sample_id" = "entity.pair_id")
  ) %>%
  relocate(case_sample, .after = sample_id)

## Remove sex chromosome genes
## Remove Neutral calls
## Remove calls for pairs not in the filtered list
## Remove low confidence calls
combined_cnvkit_filtered = combined_cnvkit_mapped %>%
  filter(chromosome != "Y", chromosome != "X", cnv_call != "Neutral") %>%
  drop_na() %>% 
  filter(p_bintest < 0.001)

noise_check <- combined_cnvkit_filtered %>%
  filter(cn != 2) %>%
  group_by(sample_id) %>%
  summarise(
    low_probe_calls = sum(probes < 3),
    total_calls = n(),
    noise_ratio = low_probe_calls / total_calls
  )

ggplot(noise_check,aes(total_calls)) + geom_histogram()

sample_depth_stats <- combined_cnvkit_filtered %>%
  group_by(sample_id) %>%
  summarise(
    mean_depth = mean(depth),
    sd_depth = sd(depth),
    coefficient_of_variation = sd_depth / mean_depth # Lower is better/more uniform
  ) %>%
  arrange(mean_depth)

# Identify genes with dangerously low coverage
low_coverage_genes <- combined_cnvkit_filtered %>%
  filter(depth < 10)

## Rename columns and subset
combined_cnvkit_filtered_subset = combined_cnvkit_filtered %>%
  dplyr::select(gene,case_sample,cnv_call) %>%
  rename(c("gene"="Gene","case_sample"="Sample_name","cnv_call"="CN"))

ggplot(combined_cnvkit_filtered %>% filter(log2 > 1),aes(x=cn)) + geom_histogram()

write.csv(combined_cnvkit_filtered_subset,"../data/processed/tumor_WES/combined_cnvkit_calls.csv",row.names = FALSE,quote = FALSE)
