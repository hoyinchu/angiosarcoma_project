library(ggplot2)
library(cowplot)
library(tidyr)
library(dplyr)

setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts")

# Import processed clinical Data
clin_df = read.csv("data/processed/clinical_data.tsv",sep="\t",check.names = FALSE)

# Number of total patients
num_total_patients = length(unique(clin_df$`STUDY ID`))

# Average age of onset

