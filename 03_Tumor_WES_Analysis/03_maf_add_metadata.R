library(dplyr)
library(ggplot2)
library(forcats)
library(ggpubr)

## This script prepares the metadata in a format that is compatible with downstream analysis

# Read MAF
maf_path = "../data/processed/tumor_WES/ASC_mutations.maf"
#maf_path = "../data/processed/tumor_WES/ASC_mutations_new.maf"

maf_df = read.csv(maf_path,sep="\t")

## Also read in cnTable
asc_cn_call_path = "../data/processed/tumor_WES/combined_cnvkit_calls.csv"
asc_cn_df = read.csv(asc_cn_call_path)

asc_maf_samples = unique(maf_df$Tumor_Sample_Barcode)
asc_cn_samples = unique(asc_cn_df$Sample_name)
meta_keep_samples = union(asc_maf_samples,asc_cn_samples)
#meta_keep_samples = intersect(asc_maf_samples,asc_cn_samples)

# clinical meta information
clin_sample_meta_path = "../data/processed/sample_clin_data.tsv"
sample_meta_df = read.csv(clin_sample_meta_path,check.names = FALSE,sep="\t")


# subset to samples in MAF
#sample_meta_df_subset = sample_meta_df[sample_meta_df$`entity:sample_id` %in% maf_df$Tumor_Sample_Barcode,]
sample_meta_df_subset = sample_meta_df[sample_meta_df$`entity:sample_id` %in% meta_keep_samples,]

sample_meta_df_subset

# Seurat metadata 
# seurat_meta_path = "data/processed/rna/02_seurat_processed_metadata.tsv"
# seurat_meta_data = read.csv(seurat_meta_path,sep="\t",check.names = FALSE)
# seurat_meta_short = seurat_meta_data[,c(setdiff(colnames(seurat_meta_data), colnames(sample_meta_df_subset)),"sample_alias","entity:sample_id")]
# seurat_meta_short[["paired_RNA_sample"]] = seurat_meta_short[["entity:sample_id"]]
# seurat_meta_short[["entity:sample_id"]] = NULL

# Add seurat information to sample meta df subset
# clin_seurat_merged = merge(sample_meta_df_subset,seurat_meta_short,by="sample_alias",all.x=TRUE)
# dim(clin_seurat_merged)

# Count the number of total mutations and non-syn mutations per sample
# Calculate the TMB as well assuming the size of the exome is 30Mb
exome_size_mb = 30
total_mutations_per_sample = maf_df %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(Total_Mutations = n()) %>%
  mutate(TMB_all_mutations = Total_Mutations / exome_size_mb)

maftools_nonsilent_classes = c(
  "Missense_Mutation","Splice_Site","Nonsense_Mutation",
  "Nonstop_Mutation", "Frame_Shift_Ins","Frame_Shift_Del",
  "In_Frame_Ins","In_Frame_Del"
  #"Start_Codon_SNP", "De_novo_Start_InFrame","De_novo_Start_OutOfFrame",
)
total_nonsyn_mutations_per_sample = maf_df %>%
  filter(Variant_Classification %in% maftools_nonsilent_classes) %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(Total_Nonsyn_Mutations = n()) %>%
  mutate(TMB_nonsyn = Total_Nonsyn_Mutations / exome_size_mb)

# Merge the two
total_mutations_per_sample_combined = merge(total_mutations_per_sample,total_nonsyn_mutations_per_sample,by="Tumor_Sample_Barcode",all.x=TRUE)

# Inspect the IQR of TMB
quantile(total_mutations_per_sample_combined$TMB_all_mutations,prob=c(.25,.5,.75))

## Add meta info back
sample_meta_df_subset_smol = sample_meta_df_subset %>%
  dplyr::select(c("entity:sample_id","sample_alias","individual_alias"))
total_mutations_per_sample_combined_with_meta = merge(
  total_mutations_per_sample_combined,sample_meta_df_subset_smol,
  by.x="Tumor_Sample_Barcode",by.y="entity:sample_id",all.x=TRUE
)

## Also add fraction genome altered
cnvkit_cns_file = read_tsv("../data/processed/tumor_WES/ASC_cnvkit_pair_cns.tsv")
cnvkit_cns_file$length = cnvkit_cns_file$end - cnvkit_cns_file$start
cnvkit_cns_log2_threshold = 1
cnvkit_fga_summary = cnvkit_cns_file %>%
  group_by(pair_id) %>%
  summarise(
    total_length = sum(length, na.rm = TRUE),
    altered_length = sum(length[abs(log2) > cnvkit_cns_log2_threshold], na.rm = TRUE),
    fga = altered_length / total_length
  ) %>%
  arrange(desc(fga))
cnvkit_fga_summary

pair_df = read.csv("../data/raw/sample_tables/pair_Jan13_2024.tsv",sep="\t",check.names = FALSE)
pair_df["case_sample"] = sub('.*entityName:([^,}]*)[},].*', '\\1', pair_df$case_sample)
pair_df["control_sample"] = sub('.*entityName:([^,}]*)[},].*', '\\1', pair_df$control_sample)
pair_df["case_sample"] = gsub("__", "_", pair_df$case_sample)
cnvkit_fga_summary_merged = merge(cnvkit_fga_summary,pair_df[,c("entity:pair_id","case_sample")],by.x="pair_id",by.y="entity:pair_id")


ggplot(cnvkit_fga_summary,aes(x=fga)) + geom_histogram()

total_mutations_per_sample_combined_with_meta_merged = merge(
  total_mutations_per_sample_combined_with_meta,
  cnvkit_fga_summary[,c("pair_id","fga")],
  by="pair_id",all.x=TRUE
)



if (FALSE) {
  write.csv(total_mutations_per_sample_combined_with_meta,file="./outputs/tables/total_mutations_per_sample.csv",row.names=FALSE)
}



# Read in the mutational signature decomposition outputs
sig_matrix = read.csv("../data/processed/tumor_WES/ASC_mutational_signature_relative_new.tsv",sep="\t",check.names = FALSE,row.names = 1)
sig_matrix_t = as.data.frame(t(sig_matrix),check.names=FALSE)
sig_matrix_t = cbind(Tumor_Sample_Barcode=rownames(sig_matrix_t),sig_matrix_t)
## Add the sum for SBS7
sig_matrix_t$SBS7 = sig_matrix_t$SBS7a + sig_matrix_t$SBS7b + sig_matrix_t$SBS7c + sig_matrix_t$SBS7d

# Add tmb and signature decomposition to original metadata
sample_meta_df_added = merge(sample_meta_df_subset,total_mutations_per_sample_combined,by.x="entity:sample_id",by.y="Tumor_Sample_Barcode",all.x=TRUE)
sample_meta_df_added = merge(sample_meta_df_added,sig_matrix_t,by.x="entity:sample_id",by.y="Tumor_Sample_Barcode",all.x=TRUE)

## Plot TMB by primary sites
tmb_by_site_plot = ggplot(sample_meta_df_added,aes(x=fct_reorder(`Primary Site (Recombined)`,-TMB_all_mutations),y=TMB_all_mutations,color =`Primary Site (Recombined)`)) +
  #geom_violin() +
  geom_boxplot(color="black",outlier.shape = NA)+
  #geom_point(size=0.5) +
  scale_y_continuous(trans='log10') +
  #labs(x="Primary Site",y="TMB") +
  #annotation_logticks(sides = 'l') +
  stat_compare_means(comparisons = list(
    c("Breast (Cutaneous)","Breast (Parenchymal)"),
    c("HNFS","Breast (Cutaneous)"),
    c("HNFS","Breast (Parenchymal)")
    ),label="p.format") +
  geom_jitter(width=0.2,size=0.5) + pretty_plot() +
  scale_color_manual(values=primary_site_palette) +
  L_border() +
  theme(axis.text.x = element_text(size=6))
  
tmb_by_site_plot = tmb_by_site_plot + theme(legend.position = "none",axis.title = element_blank())
tmb_by_site_plot
ggsave(filename="./outputs/plots/03_TMB_by_primary_sites.pdf",tmb_by_site_plot,dpi=300,height = 1.5,width=4.2)
#ggsave(filename="./outputs/plots/03_TMB_by_primary_sites_new.pdf",tmb_by_site_plot,dpi=300,height = 1.5,width=2.6)

## Get TMB mean of each site
sample_meta_df_added %>% group_by(`Primary Site (Recombined)`) %>% summarise(median(TMB_all_mutations,na.rm=TRUE))

# Change column name and order for downstream compatibility
processed_meta_data = sample_meta_df_added %>% rename("entity:sample_id"="Tumor_Sample_Barcode")
#write.csv(processed_meta_data, "../data/processed/tumor_WES/ASC_mutations_metadata.csv",row.names = FALSE)
write_tsv(processed_meta_data, "../data/processed/tumor_WES/ASC_mutations_metadata.tsv")

write.csv(processed_meta_data, "../data/processed/tumor_WES/ASC_mutations_metadata_new.csv",row.names = FALSE)
#write.csv(processed_meta_data, "../data/processed/tumor_WES/ASC_mutations_metadata_new_small.csv",row.names = FALSE)


#write.table(processed_meta_data, "../data/processed/tumor_WES/ASC_mutations_metadata.tsv", sep = "\t",quote = TRUE,row.names = FALSE)






