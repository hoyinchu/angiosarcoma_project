library(dplyr)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(stringr)

## Load the preprocessed one-hot encoded sample by pathogenic variant carrier status dataframe. 
germline_df = read.csv("../data/processed/germline/Jun2023_ASC_Case_Control_Sample_Merged_Germline.tsv",sep="\t",check.names=FALSE)
germline_case_subset = germline_df[germline_df$is_ASC=="True",]

## Load the preprocessed somatic variant table 
# Load processed MAF file and clinical metadata
asc_maf_path = "../data/processed/tumor_WES/ASC_mutations.maf"
asc_maf_metadata_path = "../data/processed/tumor_WES/ASC_mutations_metadata.tsv"
asc_maf = read.maf(maf=asc_maf_path,clinicalData=asc_maf_metadata_path)

## Make a long format dataframe that lists all the patients 
germline_cols =  names(germline_case_subset)[grep("^Germline", names(germline_case_subset))]
germline_cols_to_subset = c("Sample",germline_cols)
germline_subset = germline_case_subset[,germline_cols_to_subset]
germline_subset_long = germline_subset %>%
  gather(gene, germline_pv_carrier, -Sample) %>%
  filter(germline_pv_carrier == 1)
germline_subset_long$Hugo_Symbol <- sub("Germline_(.*)_Patho", "\\1", germline_subset_long$gene)
germline_subset_long$individual_alias = sapply(strsplit(germline_subset_long$Sample, "_"), function(x) paste(x[1:2], collapse = "_"))

# For every putative germline PV, look in the correponding tumor sample to see if a second hit occured
maf_symbol_by_sample = asc_maf@data[,c("Hugo_Symbol","Tumor_Sample_Barcode")]
maf_sample_by_patient = asc_maf@clinical.data[,c("Tumor_Sample_Barcode","individual_alias")]
maf_symbol_by_patient = left_join(maf_symbol_by_sample,maf_sample_by_patient,by=c("Tumor_Sample_Barcode"))
maf_symbol_by_patient$somatic_nonsyn_carrier = 1

# Merge the two
germline_somatic_merged = merge(germline_subset_long,maf_symbol_by_patient,by=c("individual_alias","Hugo_Symbol"),all.x=TRUE,all.y=TRUE)
germline_somatic_merged = germline_somatic_merged %>% dplyr::select(individual_alias,Hugo_Symbol,germline_pv_carrier,somatic_nonsyn_carrier,Tumor_Sample_Barcode)
germline_somatic_merged[is.na(germline_somatic_merged$germline_pv_carrier),"germline_pv_carrier"] = 0
germline_somatic_merged[is.na(germline_somatic_merged$somatic_nonsyn_carrier),"somatic_nonsyn_carrier"] = 0

## Per gene, how many patient has mutation in both tumor and germline in the same gene?
germline_somatic_merged_no_dup = germline_somatic_merged[!duplicated(germline_somatic_merged[c("individual_alias","Hugo_Symbol")]),]
germline_somatic_merged_no_dup$germline_and_somatic_carrier = germline_somatic_merged_no_dup$germline_pv_carrier * germline_somatic_merged_no_dup$somatic_nonsyn_carrier

# Count the number of patients with germline/somatic mutations in each gene
germline_versus_tumor_mut_counts = germline_somatic_merged_no_dup %>% group_by(Hugo_Symbol) %>% summarise(
  germline_pv_carrier_count=sum(germline_pv_carrier,na.rm=TRUE),
  somatic_nonsyn_carrier_count=sum(somatic_nonsyn_carrier,na.rm=TRUE),
  germline_and_somatic_count=sum(germline_and_somatic_carrier,na.rm=TRUE)
)

# For visualization purposes make a discretized cersion
germline_versus_tumor_mut_counts$germline_and_somatic_count_chr = as.factor(germline_versus_tumor_mut_counts$germline_and_somatic_count)
germline_versus_tumor_mut_counts = germline_versus_tumor_mut_counts %>% mutate(
  bicarrier_category = case_when(
    germline_versus_tumor_mut_counts$germline_and_somatic_count == 5 ~ "5",
    germline_versus_tumor_mut_counts$germline_and_somatic_count == 4 ~ "4",
    germline_versus_tumor_mut_counts$germline_and_somatic_count == 1 ~ "1",
    germline_and_somatic_count == 0 ~ "0"
  )
)

always_show_genes = germline_versus_tumor_mut_counts[germline_versus_tumor_mut_counts$germline_and_somatic_count > 0,]
germline_vs_tumor_label_subset = subset(germline_versus_tumor_mut_counts,
                                        ((germline_pv_carrier_count > 2 & somatic_nonsyn_carrier_count>1) |
                                          (germline_pv_carrier_count > 3) | 
                                          (somatic_nonsyn_carrier_count > 3) |
                                          (germline_and_somatic_count > 0)) & (!(Hugo_Symbol %in% always_show_genes$Hugo_Symbol)) 
)

germline_versus_tumor_gene_plot = ggplot(
  germline_versus_tumor_mut_counts,
  aes(x=germline_pv_carrier_count,y=somatic_nonsyn_carrier_count,label=Hugo_Symbol,color=bicarrier_category)) +
  geom_point(aes(size=germline_and_somatic_count)) +
  geom_label_repel(data=germline_vs_tumor_label_subset,max.overlaps = 10) +
  geom_label_repel(data=always_show_genes, box.padding = 0.5) +
  theme_minimal() +
  theme(legend.position = "top") +
  guides(size="none") +
  labs(
    x="# of Patients with Germline Pathogenic Variants in Gene",
    y="# of Patients with Non-synonomous Mutations in Gene",
    color="# of Patients with double hits"
    )

germline_versus_tumor_gene_plot
  
#ggsave("./outputs/plots/00_germline_versus_tumor_gene_plot.png",germline_versus_tumor_gene_plot,dpi=300,width=10,height=8)

germline_versus_tumor_mut_counts
## Make a rank plot of total number of patients with both germline and somatic variants in the same gene
germline_versus_tumor_mut_counts_filtered = germline_versus_tumor_mut_counts# %>% filter(germline_pv_carrier_count > 0)
germline_versus_tumor_mut_counts_filtered$bicarrier_category_rank = rank(-as.numeric(germline_versus_tumor_mut_counts_filtered$bicarrier_category),ties.method = "first")
rank_range = range(germline_versus_tumor_mut_counts_filtered$bicarrier_category_rank, na.rm = TRUE)
biallelic_rank_plot = ggplot(data=germline_versus_tumor_mut_counts_filtered,aes(x=bicarrier_category_rank,y=bicarrier_category,color=bicarrier_category>1)) +
  geom_point() + pretty_plot() + L_border() + scale_color_manual(values = c("TRUE"="firebrick","FALSE"="black")) +
  scale_x_continuous(breaks = rank_range) +
  #scale_y_discrete(expand = c(0,0.075)) +
  theme(legend.position = "none", axis.title = element_blank())
biallelic_rank_plot
cowplot::ggsave2("./outputs/plots/biallelic_rank_plot.pdf",biallelic_rank_plot,dpi=300,width=1.6,height=1.3)

## Alternatively plot it this way
highlight_genes = c("MUC16","TP53","KDR","PKHD1","USH2A","CFTR","FLG","PAH","POT1","TTN","CYP21A2","GJB2")
germline_versus_tumor_mut_counts$to_highlight = germline_versus_tumor_mut_counts$Hugo_Symbol %in% highlight_genes
germline_versus_tumor_gene_plot = ggplot(
  germline_versus_tumor_mut_counts %>% arrange(to_highlight),aes(x=germline_pv_carrier_count,y=somatic_nonsyn_carrier_count,color=to_highlight)
) + geom_point() + pretty_plot() + L_border() + scale_color_manual(values=c("TRUE"="firebrick","FALSE"="black")) +
  geom_text_repel(data=germline_versus_tumor_mut_counts %>% filter(to_highlight),aes(label=Hugo_Symbol),size=2) +
  theme(legend.position = "none", axis.title = element_blank())
germline_versus_tumor_gene_plot
cowplot::ggsave2("./outputs/plots/germline_versus_tumor_gene_plot.pdf",germline_versus_tumor_gene_plot,dpi=300,width=1.6,height=1.3)

## Write the count table 
write.csv(germline_versus_tumor_mut_counts,"./outputs/tables/germline_vs_tumor_gene_counts.csv",row.names = FALSE)

## Load the count table
germline_versus_tumor_mut_counts = read.csv("./outputs/tables/germline_vs_tumor_gene_counts.csv")

## Also add other clinically relevant information
clin_data = read.csv("../data/processed/clinical_data.tsv",sep="\t",check.names = FALSE)
clin_data = read.csv("../data/processed/clinical_data_manual_deident.tsv",sep="\t",check.names = FALSE)
clin_data$individual_alias = clin_data$`STUDY ID`
clin_data$germline_avail = as.numeric(clin_data$individual_alias %in% germline_case_subset$individual_alias)
#clin_data$germline_avail = as.numeric(clin_data$individual_alias %in% germline_case_alias_ids)

clin_data$somatic_avail = as.numeric(clin_data$individual_alias %in% unique(maf_symbol_by_patient$individual_alias))
clin_data$germline_somatic_avail = clin_data$germline_avail * clin_data$somatic_avail
clin_data = clin_data %>% mutate(
  germline_somatic_avail_status = case_when(
    germline_somatic_avail == 1 ~ "Somatic & Germline WES",
    (germline_avail == 1) & (somatic_avail == 0) ~ "Germline WES Only",
    (germline_avail == 0) & (somatic_avail == 1) ~ "Somatic WES Only",
    (germline_avail == 0) & (somatic_avail == 0) ~ "No WES available",
  )
) %>% mutate(germline_somatic_avail_status = factor(
  germline_somatic_avail_status,
  levels = c("Somatic & Germline WES","Somatic WES Only","Germline WES Only","No WES available")
)) 

pot1_biallelic_carriers = germline_somatic_merged_no_dup[
  (germline_somatic_merged_no_dup$Hugo_Symbol=="POT1")&
  (germline_somatic_merged_no_dup$germline_and_somatic_carrier==1),
  ]$individual_alias
pot1_somatic_only_carriers = germline_somatic_merged_no_dup[
  (germline_somatic_merged_no_dup$Hugo_Symbol=="POT1")&
    (germline_somatic_merged_no_dup$somatic_nonsyn_carrier==1)&
    (germline_somatic_merged_no_dup$germline_pv_carrier==0),
]$individual_alias
pot1_germline_only_carriers = germline_somatic_merged_no_dup[
  (germline_somatic_merged_no_dup$Hugo_Symbol=="POT1")&
    (germline_somatic_merged_no_dup$somatic_nonsyn_carrier==0)&
    (germline_somatic_merged_no_dup$germline_pv_carrier==1),
]$individual_alias

clin_data = clin_data %>% mutate(
  pot1_mutation_status = case_when(
    individual_alias %in% pot1_biallelic_carriers ~ "Somatic + Germline",
    individual_alias %in% pot1_somatic_only_carriers ~ "Somatic Only",
    individual_alias %in% pot1_germline_only_carriers ~ "Germline Only",
    TRUE ~ "None Detected"
  )
)  %>% mutate(pot1_mutation_status = factor(
  pot1_mutation_status,
  #levels = c("Somatic + Germline","Somatic Only","Germline Only","None Detected"),
  levels = c("None Detected","Somatic Only","Germline Only","Somatic + Germline")
  ))

clin_data 


## Make the age plot
pot1_age_plot = ggplot(clin_data %>% filter(`Primary Site (Recombined)` != "Missing or Unknown"), aes(x = pot1_mutation_status, y = `Age (Combined)`)) +
  pretty_plot() + 
  L_border() +
  geom_jitter(aes(fill = `Primary Site (Recombined)`), 
              shape = 21, width = 0.1, color = "black", size = 1, stroke = 0.3) +
  # 3. Keep boxplot colors static (black lines, white fill)
  geom_boxplot(width = 0.4, fill = "white", color = "black", outlier.shape = NA) +
  stat_compare_means(comparisons = list(
    c("Somatic Only", "None Detected"),
    c("Germline Only", "None Detected"),
    c("Somatic + Germline", "None Detected")
  )) +
  scale_fill_manual(values=primary_site_palette) +
  theme(axis.title = element_blank(), axis.text.x = element_blank()) +
  theme(legend.position = "none")

pot1_age_plot
cowplot::ggsave2("./outputs/plots/pot1_age_of_onset_woth_color.pdf",pot1_age_plot,dpi=300,width=1.8,height=1.2)
cowplot::ggsave2("./outputs/plots/pot1_age_of_onset_with_legend.pdf",pot1_age_plot,dpi=300,width=1.8,height=1.2)


clin_data %>%
  group_by(pot1_mutation_status) %>%
  summarize(mean_value = mean(`Age (Combined)`, na.rm = TRUE))
