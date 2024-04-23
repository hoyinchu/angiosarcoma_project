library(dplyr)

## Load the preprocessed one-hot encoded sample by pathogenic variant carrier status dataframe. Takes a while
germline_df = read.csv("reference_data/ClinicalTables/Jun2023_ASC_Case_Control_Sample_Merged_Germline.tsv",sep="\t")

