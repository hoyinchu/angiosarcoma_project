library(dplyr)
library(maftools)

## This script filters the MAF to only retain high confidence tumor variant calls
maf_path = "../data/raw/tumor_WES/CombinedMAF_20260202.maf"
maf_df = read.csv(maf_path,sep="\t")

# A single sample's name formatted incorrectly
maf_df["Tumor_Sample_Barcode"] = gsub("ASCProject_0194__", "ASCProject_0194_", maf_df$Tumor_Sample_Barcode)

## Sample-level filter
## Discard variants in samples that did not pass QC
filter_status_df = read.csv("../data/processed/qc_tables/tumor_WES_pair_status.tsv",check.names = FALSE,sep="\t")
maf_df_merged = merge(
  maf_df,
  filter_status_df,
  by.x=c("Tumor_Sample_Barcode","Matched_Norm_Sample_Barcode"),
  by.y=c("case_sample","control_sample"),
  all.x=TRUE
  )
filtered_maf = maf_df_merged[maf_df_merged$pass_final_filter,]

final_filter_passed_samples = filter_status_df %>% filter(pass_final_filter) %>% pull(case_sample)
maf_og_samples = unique(maf_df$Tumor_Sample_Barcode)
maf_filtered_samples = unique(filtered_maf$Tumor_Sample_Barcode)

## Variant level filters
# Filter out variants falling into the blacklist regions as provided by encode
# The ENCODE Blacklist: Identification of Problematic Regions of the Genome
# Haley M. Amemiya, Anshul Kundaje & Alan P. Boyle 
blacklist_regions = read.csv("../data/public/encode-hg19-blacklist.v2.bed",sep="\t",header = FALSE)
colnames(blacklist_regions) = c("Chromosome","Start_position","End_position","Blacklist_reason")
blacklist_regions["Chromosome"] = gsub("chr","",blacklist_regions$Chromosome)

check_overlap = function(chrom, start, end, blacklist) {
  overlap = blacklist$Chromosome == chrom & blacklist$Start_position <= end & blacklist$End_position >= start
  any(overlap)
}
filtered_maf = filtered_maf %>% 
  rowwise() %>% 
  filter(!check_overlap(Chromosome, Start_position, End_position, blacklist_regions))

## Filter out variants with a failure reason
filtered_maf = filtered_maf %>% filter(failure_reasons == "" | is.na(failure_reasons))

## Filter out variants suscepted to be an oxoG or ffpe artifact
filtered_maf = filtered_maf %>% filter(i_oxog_p_value < i_oxog_p_value_cutoff | is.na(i_oxog_p_value))
filtered_maf = filtered_maf %>% filter(i_ffpe_p_value < i_ffpe_p_value_cutoff | is.na(i_oxog_p_value))

## Filter out variants that are detected only on one strand
filtered_maf = filtered_maf %>% filter(detected_in_both_strands == "True")

filtered_maf["tumor_mutation_id"] = with(filtered_maf, paste(Hugo_Symbol,Chromosome, Start_position, End_position, Tumor_Seq_Allele1,Tumor_Seq_Allele2,Protein_Change,sep="_"))
filtered_maf["tumor_mutation_id_short"] = with(filtered_maf, paste(Hugo_Symbol,Protein_Change,sep="_"))

## Variants that occur in more than two individuals are subject to manual IGV review
sample_meta_df = read.csv("../data/processed/sample_clin_data.tsv",check.names = FALSE,sep="\t")
filtered_maf_merged = merge(filtered_maf,sample_meta_df,by.x="Tumor_Sample_Barcode",by.y="entity:sample_id",all.x=TRUE)
mut_unique_counts = filtered_maf_merged %>%
  group_by(tumor_mutation_id) %>% summarize(
    count=n(),
    unique_individuals = length(unique(individual_alias))
  )
mut_unique_counts = mut_unique_counts[order(-mut_unique_counts$unique_individuals),]
mut_unique_counts[mut_unique_counts$unique_individuals > 1,]

# igv_blacklist = c(
#   "TLN2_15_63063187_63063187_G_T_p.E134*",
#   "KCTD3_1_215792484_215792485_-_A_",
#   "ASPH_8_62588573_62588574_-_G_",
#   "ZNF654_3_88189337_88189338_-_TTT_p.293_293R>IW",
#   "TMEM232_5_109974038_109974039_-_CCCTC_p.Y121fs"
# )
# 
# filtered_maf_merged = filtered_maf_merged %>% filter(!tumor_mutation_id %in% igv_blacklist)

write.table(filtered_maf,file="../data/processed/tumor_WES/ASC_mutations.maf",sep="\t",row.names=FALSE,quote = FALSE)

## Sanity maftool check
maf_to_check = filtered_maf_merged
asc_maf = read.maf(maf=maf_to_check)
oncoplot(maf = asc_maf,
         top = 30,
         draw_titv = TRUE,
         fontSize = 0.6,
)

