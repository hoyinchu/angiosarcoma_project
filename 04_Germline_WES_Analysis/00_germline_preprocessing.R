library(dplyr)
sample_df = read.csv("../data/processed/germline/sample_Mar03_2023.tsv",sep="\t",check.names = FALSE)
germline_table = read.csv("../data/processed/germline/germline_df.tsv",sep="\t",check.names=FALSE)

# Check if sample is included after filtering 
normal_samples = sample_df %>% filter(sample_type=="Normal")
normal_samples_subset = normal_samples[,c("entity:sample_id","sample_alias","individual_alias","sampleMeanCoverage")]
normal_samples_subset["passed_germline_filter"] = normal_samples_subset$sample_alias %in% germline_table$Sample
normal_samples_subset["data_type"] = "normal_WES"
write.table(normal_samples_subset,file="../data/processed/germline/germline_sample_status.tsv",sep="\t",row.names=FALSE)

# Sample coverage
summary(normal_samples_subset$sampleMeanCoverage)

## After running BAMs through DeepVariant pipeline, the resulting VCF is filtered using notebooks
## Come back to this section after the filtering has been done 
## Get the filtered tables from notebooks
germline_filtered_variants_cases = read.csv("./outputs/tables/filtered_germline_variants_cases.tsv",check.names = FALSE,sep="\t")
germline_filtered_variants_cases

## Identify number of variants that were qualified 
dim(germline_filtered_variants_cases)
preselected_genes_cases = germline_filtered_variants_cases[germline_filtered_variants_cases$pass_germline_filter != "unselected_clinvar_pathogenic_high_conf",]

## Number of unique samples
length(unique(preselected_genes_cases$Sample))

## Convert into a MAF for visualization
maf_columns = c("Hugo_Symbol","Chromosome","Start_Position","End_Position","Reference_Allele",
                "Tumor_Seq_Allele2","Variant_Classification","Variant_Type","Tumor_Sample_Barcode","Protein_Change")

# Calculate the start / end position
germline_filtered_variants_cases_annotated = germline_filtered_variants_cases %>%
  mutate(Start_Position = as.numeric(sub(".*:(\\d+).*", "\\1", Location)),
         End_Position = as.numeric(ifelse(grepl("-", Location), sub(".*-(\\d+)", "\\1", Location), Start_Position)))
colnames(germline_filtered_variants_cases_annotated)

# Prettify protein change
germline_filtered_variants_cases_annotated = germline_filtered_variants_cases_annotated %>%
  mutate(Protein_Change = ifelse(Protein_position == "-" | Amino_acids == "-", "-",
                             paste0("p.",substr(Amino_acids, 1, 1), 
                                    as.numeric(sub("^(\\d+).*", "\\1", Protein_position)), 
                                    substr(Amino_acids, nchar(Amino_acids), nchar(Amino_acids)))))

# Unify variant class annotaiton
germline_filtered_variants_cases_annotated = germline_filtered_variants_cases_annotated %>% mutate(
  VARIANT_CLASS=case_when(
    VARIANT_CLASS=="SNV" ~ "SNP",
    VARIANT_CLASS=="deletion" ~ "DEL",
    VARIANT_CLASS=="insertion" ~ "INS"
  )
)

germline_filtered_variants_cases_annotated$Consequence_consolidated

# Change the way how some of the consequences are annotated
germline_filtered_variants_cases_annotated = germline_filtered_variants_cases_annotated %>% mutate(
  MAF_Consequence = case_when(
    Consequence_consolidated == "stop_gained" ~ "Nonsense_Mutation",
    (Consequence_consolidated == "frameshift_variant") & (VARIANT_CLASS=="INS") ~ "Frame_Shift_Ins",
    (Consequence_consolidated == "frameshift_variant") & (VARIANT_CLASS=="DEL") ~ "Frame_Shift_Del",
    Consequence_consolidated == "missense_variant" ~ "Missense_Mutation",
    Consequence_consolidated=="splice_donor_variant" ~ "Splice_Site",
    Consequence_consolidated=="splice_acceptor_variant" ~ "Splice_Site",
    Consequence_consolidated=="splice_region_variant" ~ "Splice_Site",
    Consequence_consolidated=="inframe_deletion" ~ "In_Frame_Del",
    Consequence_consolidated=="inframe_insertion" ~ "In_Frame_Ins",
    Consequence_consolidated=="intron_variant" ~ "Silent",
    Consequence_consolidated=="synonymous_variant" ~ "Silent",
    TRUE ~ Consequence_consolidated
  )
)

germline_filtered_variants_cases_annotated_subset = germline_filtered_variants_cases_annotated[
  ,c("SYMBOL","CHROM","Start_Position","End_Position","REF","ALT","MAF_Consequence","VARIANT_CLASS","Sample","Protein_Change")
]

colnames(germline_filtered_variants_cases_annotated_subset) = maf_columns
germline_filtered_variants_cases_annotated_subset
write.table(germline_filtered_variants_cases_annotated_subset,
            file = "./outputs/tables/filtered_germline_variants_cases_subset.maf", row.names=FALSE, sep="\t", quote = FALSE)

## Also change the clinical table so that it can be read by downstream library
# clinical meta information
clin_sample_meta_path = "../data/processed/sample_clin_data.tsv"
sample_meta_df = read.csv(clin_sample_meta_path,check.names = FALSE,sep="\t")
#germline_sample_clin_df = sample_meta_df[sample_meta_df$sample_alias %in% germline_filtered_variants_cases_annotated_subset$Tumor_Sample_Barcode,]
germline_sample_clin_df = sample_meta_df %>% mutate(Tumor_Sample_Barcode = sample_alias_cleaned)
## Sanitize the column names for those to be visualized
germline_sample_clin_df$Age = germline_sample_clin_df$`Age (Combined)`
germline_sample_clin_df$PrimarySite = germline_sample_clin_df$`Primary Site (Recombined)`

write.table(germline_sample_clin_df,"./outputs/tables/germline_sample_clin_df.tsv",sep="\t",quote=TRUE,row.names = FALSE)

