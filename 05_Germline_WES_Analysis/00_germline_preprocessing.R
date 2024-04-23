library(dplyr)

setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional")
sample_df = read.csv("workspace_tables/sample_Mar03_2023.tsv",sep="\t",check.names = FALSE)
germline_table = read.csv("reference_data/IntegratedData/germline_df.tsv",sep="\t",check.names=FALSE)

# Check if sample is included after filtering 
normal_samples = sample_df %>% filter(sample_type=="Normal")
normal_samples_subset = normal_samples[,c("entity:sample_id","sample_alias")]
normal_samples_subset["passed_germline_filter"] = normal_samples_subset$sample_alias %in% germline_table$Sample
normal_samples_subset["data_type"] = "normal_WES"

write.table(normal_samples_subset,file="./scripts/Germline_WES_Analysis/germline_sample_status.tsv",sep="\t",row.names=FALSE)
