library(dplyr)
library(tidyr)
library(ggplot2)
library(maftools)
library(ComplexHeatmap)
library(scales)

maf_path = "../data/processed/tumor_WES/ASC_mutations.maf"
maf_meta_path = "../data/processed/tumor_WES/ASC_mutations_metadata.tsv"
maf_df = read.csv(maf_path,sep="\t",check.names=FALSE)
maf_meta = read.csv(maf_meta_path,sep="\t",check.names=FALSE)
maf_merged = merge(maf_df,maf_meta,by="Tumor_Sample_Barcode",all.x=TRUE)

## Get a set of variants that are non-silent (using maftools' classification)
maftools_nonsilent_classes = c(
  "Missense_Mutation","Splice_Site","Nonsense_Mutation",
  "Nonstop_Mutation", "Frame_Shift_Ins","Frame_Shift_Del",
  "In_Frame_Ins","In_Frame_Del"
  #"Start_Codon_SNP", "De_novo_Start_InFrame","De_novo_Start_OutOfFrame",
)
nonsilent_maf = maf_merged %>% filter(
  Variant_Classification %in% maftools_nonsilent_classes
)

# Load processed MAF file and clinical metadata
asc_maf_path = "../data/processed/tumor_WES/ASC_mutations.maf"
asc_maf_metadata_path = "../data/processed/tumor_WES/ASC_mutations_metadata.tsv"
asc_maf = read.maf(maf=asc_maf_path,clinicalData=asc_maf_metadata_path)

## Verify that the total number of variants are the same between the two
dim(asc_maf@data)[[1]] == dim(nonsilent_maf)[[1]]

## TODO:
## 3. Make a rank plot where x is ordered in total frequency, y is frequency, show see KDR p.771R and POT1 R117C pop up on top
make_count_df = function(maf,col,colnames) {
  count_table = table(maf[[col]])
  count_df = as.data.frame(count_table)
  colnames(count_df) = colnames
  return(count_df)
}

## For each variant, get its total frequency, short name, frequency by unique individual,
## For each gene, get its total frequency and frequency by unique individual
calc_mutation_frequency = function(maf) {
  mut_count = make_count_df(maf,"tumor_mutation_id",c("tumor_mutation_id","mut_count"))
  maf_one_mut_per_patient = maf[!duplicated(maf[c("individual_alias","tumor_mutation_id")]),]
  mut_count_unique = make_count_df(maf_one_mut_per_patient,"tumor_mutation_id",c("tumor_mutation_id","mut_count_unique"))

  gene_mut_count = make_count_df(maf,"Hugo_Symbol",c("Hugo_Symbol","gene_mut_count"))
  maf_one_gene_mut_per_patient =  maf[!duplicated(maf[c("individual_alias","Hugo_Symbol")]),]
  gene_mut_count_unique = make_count_df(maf_one_gene_mut_per_patient,"Hugo_Symbol",c("Hugo_Symbol","gene_mut_count_unique"))

  mapping_df = maf[,c("tumor_mutation_id","tumor_mutation_id_short","Hugo_Symbol")]
  mapping_df = mapping_df[!duplicated(mapping_df$tumor_mutation_id),]
  mapping_df_merged = merge(mapping_df,mut_count,by="tumor_mutation_id",all.x=TRUE)
  mapping_df_merged = merge(mapping_df_merged,mut_count_unique,by="tumor_mutation_id",all.x=TRUE)
  mapping_df_merged = merge(mapping_df_merged,gene_mut_count,by="Hugo_Symbol",all.x=TRUE)
  mapping_df_merged = merge(mapping_df_merged,gene_mut_count_unique,by="Hugo_Symbol",all.x=TRUE)
  return(mapping_df_merged)
}

mut_freq_df = calc_mutation_frequency(nonsilent_maf)
mut_freq_df_mut_recurrent_only = mut_freq_df[mut_freq_df$mut_count_unique > 1,]
## Exclude non missense change
mut_freq_df_mut_recurrent_only = mut_freq_df_mut_recurrent_only %>% filter(grepl("p\\.", tumor_mutation_id_short))

## Subset MAF to only these entries
maf_subset = maf_merged[maf_merged$tumor_mutation_id %in% mut_freq_df_mut_recurrent_only$tumor_mutation_id,]
maf_subset_unique = maf_subset[!duplicated(maf_subset[c("STUDY ID","tumor_mutation_id_short")]),]
mut_by_site_table = as.data.frame(table(maf_subset_unique$tumor_mutation_id_short,maf_subset_unique$`Primary Site (Recombined)`),check.names=FALSE)
colnames(mut_by_site_table) = c("tumor_mutation_id_short","primary_site","count")

# Set level by total count
tumor_mutation_id_short_order = mut_freq_df_mut_recurrent_only[order(-mut_freq_df_mut_recurrent_only$mut_count_unique),]$tumor_mutation_id_short
mut_by_site_table$tumor_mutation_id_short = factor(mut_by_site_table$tumor_mutation_id_short,levels=tumor_mutation_id_short_order)

## Also calculate the AlphaMissense scores for each of the variant
## Uniprot ID file = uniprot_to_gene_names.txt
## AlphaMissenese file is grep to only include entries with uniprot ids of interest
uniprot_map = read.csv("../data/public/uniprot_to_gene_names.txt",sep="\t",header=FALSE)
colnames(uniprot_map) = c("uniprot_id","id_type","gene_name")
mut_freq_df_mut_recurrent_only_merged = merge(mut_freq_df_mut_recurrent_only,
                                              uniprot_map,by.x="Hugo_Symbol",by.y="gene_name",how="left")

## This requires downloading the AlphaMissense csv
## You can get it from here:
##
alpha_missense_df = read.csv("../data/public/AlphaMissense_hg19.tsv",sep="\t",skip = 3)
alpha_missense_df_subset = alpha_missense_df %>% filter(
  uniprot_id %in% mut_freq_df_mut_recurrent_only_merged$uniprot_id,
)
write.csv(alpha_missense_df_subset,"../data/public/AlphaMissense_hg19_recurrent_subset.csv",row.names = FALSE)
# Rename mutation ID by protein name
mut_freq_df_mut_recurrent_only_merged$aa_change = strsplit(mut_freq_df_mut_recurrent_only_merged$tumor_mutation_id_short,"\\.") 
mut_freq_df_mut_recurrent_only_merged$aa_change = sapply(mut_freq_df_mut_recurrent_only_merged$aa_change, function(x) if (length(x) >= 2) x[2] else NA)
mut_freq_df_mut_recurrent_only_merged$uniprot_aa_change = paste0(mut_freq_df_mut_recurrent_only_merged$uniprot_id,
                                                                 "_",
                                                                 mut_freq_df_mut_recurrent_only_merged$aa_change
                                                                 )
alpha_missense_df_subset$uniprot_aa_change = paste0(
  alpha_missense_df_subset$uniprot_id,
  "_",
  alpha_missense_df_subset$protein_variant
)
# Make # of samples mutated vs. Alpha Missense prediction
mut_freq_df_mut_recurrent_only_with_score = merge(mut_freq_df_mut_recurrent_only_merged,
                                                  alpha_missense_df_subset,by="uniprot_aa_change",
                                                  all.x=TRUE) %>%
  drop_na() %>%
  arrange(-am_pathogenicity) %>%
  distinct(tumor_mutation_id,.keep_all = TRUE)
mut_freq_df_mut_recurrent_only_with_score


write.csv(uniprot_genes,"./outputs/uniprot_id_for_gene_with_recurrent_mutated_variants.csv",row.names = FALSE)

mut_count_alphamissense_plot = ggplot(mut_freq_df_mut_recurrent_only_with_score,aes(x=mut_count_unique,y=am_pathogenicity,label=tumor_mutation_id_short)) +
  geom_point() +
  ggrepel::geom_label_repel() +
  theme_minimal() +
  labs(x="# of Patients with Mutation",y="AlphaMissense Pathogenicity Score")
mut_count_alphamissense_plot
ggsave("./outputs/plots/05_variant_patients_by_am_pathogenicity.png",mut_count_alphamissense_plot,dpi=300,height=6,width=12)
ggsave("./outputs/plots/05_variant_patients_by_am_pathogenicity.pdf",mut_count_alphamissense_plot,dpi=300,height=6,width=12)

## Label the mutations marked as pathogenic by alphamissense
am_pathogenic_mutations = mut_freq_df_mut_recurrent_only_with_score[,c("tumor_mutation_id_short","am_pathogenicity","am_class")]
mut_by_site_table_am_merged = merge(mut_by_site_table,am_pathogenic_mutations,by="tumor_mutation_id_short",all.x=TRUE)
table(mut_by_site_table_am_merged$am_class)

am_label_data = mut_by_site_table_am_merged %>%
  group_by(tumor_mutation_id_short) %>%
  summarize(
    total_count = sum(count),
    am_class = first(am_class) # Assuming one am_class per mutation ID
  ) %>%
  filter(am_class == "pathogenic")

recurrent_mutation_plot = ggplot(mut_by_site_table_am_merged,aes(x=tumor_mutation_id_short,y=count,fill=primary_site)) +
  geom_bar(stat = "identity",position = "stack") +
  geom_text(data = am_label_data, 
            aes(x = tumor_mutation_id_short, y = total_count, label = "*"),
            inherit.aes = FALSE,    # Prevents searching for 'fill' in label_data
            vjust = -0.2,           # Adjust slightly above the bar
            size = 5,               # Size of the asterisk
            color = "black") +
  scale_y_continuous(breaks= pretty_breaks(),expand = expansion(mult = c(0, 0.15))) + 
  pretty_plot() + L_border() +
  labs(x="Gene",y="Number of Patients with Mutation",fill="Primary Site") +
  scale_fill_manual(values=primary_site_palette) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=5)) + 
  theme(legend.position="top") +
  guides(fill = guide_legend(nrow = 2))


recurrent_mutation_plot
ggsave(recurrent_mutation_plot,file="./outputs/plots/05_recurrent_mutation_plot.png",height=4,width=6)
ggsave(recurrent_mutation_plot,file="./outputs/plots/05_recurrent_mutation_plot.pdf",height=2,width=6)
ggsave(recurrent_mutation_plot,file="./outputs/plots/05_recurrent_mutation_plot_no_legend.pdf",height=1.3,width=5.2)
ggsave(recurrent_mutation_plot,file="./outputs/plots/05_recurrent_mutation_plot_with_legend.pdf",height=1.3,width=5.2)

## Also check the VAF for each variant
mut_by_site_table_am_merged

mut_by_site_table

temp_data = asc_maf@data
colnames(temp_data) = make.unique(colnames(temp_data))
asc_maf_vaf_subset = temp_data %>% 
  filter(tumor_mutation_id_short %in% mut_by_site_table_am_merged$tumor_mutation_id_short)

asc_maf_vaf_subset$t_VAF = asc_maf_vaf_subset$t_alt_count/asc_maf_vaf_subset$t_depth
asc_maf_vaf_subset$tumor_mutation_id_short = factor(asc_maf_vaf_subset$tumor_mutation_id_short,levels=unique(mut_by_site_table_am_merged$tumor_mutation_id_short))
asc_maf_vaf_subset_merged = merge(asc_maf_vaf_subset,asc_maf@clinical.data[,c("Tumor_Sample_Barcode","Primary_Site_(Recombined)")],by="Tumor_Sample_Barcode",all.x=TRUE)


recurrent_t_vaf_plot = ggplot(asc_maf_vaf_subset_merged,aes(x=tumor_mutation_id_short,y=t_VAF)) + 
  geom_boxplot() +
  geom_point(aes(color=`Primary_Site_(Recombined)`,alpha=0.9)) +
  scale_color_manual(values=primary_site_palette) +
  pretty_plot() + L_border() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,size=5),axis.title.x = element_blank(),axis.title.y = element_blank()) +
  theme(legend.position = "bottom") #theme(legend.position = "none")

recurrent_t_vaf_plot
ggsave(recurrent_t_vaf_plot,file="./outputs/plots/recurrent_t_vaf_plot_with_legend.pdf",height=3,width=5.2)
ggsave(recurrent_t_vaf_plot,file="./outputs/plots/recurrent_t_vaf_plot_no_legend.pdf",height=3,width=5.2)

asc_maf_vaf_subset_merged %>%
  group_by(tumor_mutation_id_short) %>%
  summarize(mean_value = median(t_VAF, na.rm = TRUE))

## Make a "representative somatic / germline mutation" per sample table that has somatic vs. germline mutation?
## Somatic Logic: Hotspot Mutation -> Mutation in Mutsig Significant Gene -> Mutation in COSMIC Cancer Tier 1 gene
## Germline Logic: POT1 -> Mutation in COSMIC Cancer Tier 1 gene -> Other recurrent mutations?

# Indicate whether a mutation is recurrent across patients (hotspots)
recurrent_mutations_short_ids = tumor_mutation_id_short_order
# or is in a MutSig significant gene
mutsig_output_path = "../data/processed/tumor_WES/Apr16_2024_sig_genes.txt"
mutsig_gene = read.csv(mutsig_output_path,sep="\t")
mutsig_q10_genes = mutsig_gene[mutsig_gene$q < 0.1,]$gene
# Or if it is a COSMIC tier 1 gene
cosmic_path = "../data/public/cosmic_cancer_gene_census.csv"
cosmic_df = read.csv(cosmic_path)
cosmic_tier1_genes = cosmic_df[cosmic_df$Tier==1,]$Gene.Symbol

# Or if germline pathogenic variant was also detected in the gene
# Getting this table requires running a script from a later step
## First run through the germline scripts (04)
## Then run the germline somatic integration scripts (05/00_two_hit_scan.R)
# So skip this step if it's the first time running this

somatic_germline_mut_table = read.csv("../06_Germline_WES_Tumor_WES_Analysis/outputs/tables/germline_vs_tumor_gene_counts.csv")
somatic_gerline_mut_genes = somatic_germline_mut_table %>%
  filter(germline_pv_carrier_count > 0, somatic_nonsyn_carrier_count>0) %>%
  pull(Hugo_Symbol)
# Or if it is in a gene of interest
goi_list = c("CFTR")
nonsilent_maf = nonsilent_maf %>% mutate(
  is_recurrent_mutation = tumor_mutation_id_short %in% recurrent_mutations_short_ids,
  is_mutsig_significant = Hugo_Symbol %in% mutsig_q10_genes,
  is_cosmic_tier1_gene = Hugo_Symbol %in% cosmic_tier1_genes,
  has_germline_somatic_mut = Hugo_Symbol %in% somatic_gerline_mut_genes,
  is_in_gene_of_interest = Hugo_Symbol %in% goi_list
)
## Assign priority scores to mutation based on these criteria
nonsilent_maf = nonsilent_maf %>% mutate(
  mut_priority_score = (is_recurrent_mutation*2) + 
    (is_mutsig_significant*1) + 
    (is_cosmic_tier1_gene*1) +
    (has_germline_somatic_mut*1) +
    (is_in_gene_of_interest * 3)
)
maf_high_score_muts_only = nonsilent_maf %>%
  arrange(-mut_priority_score) %>%
  #filter(mut_priority_score > 1) %>%
  distinct(Tumor_Sample_Barcode,.keep_all = TRUE)
dim(maf_high_score_muts_only)
maf_high_score_muts_only$Hugo_Symbol
maf_high_score_muts_only[,c("Tumor_Sample_Barcode","tumor_mutation_id","Hugo_Symbol","mut_priority_score")]

## Write down the "representative mutations" table somewhere
write.csv(maf_high_score_muts_only,"03_Tumor_WES_Analysis/outputs/tables/representative_mutation_table.csv",row.names = FALSE)

## Make lollipop plot for KDR
pdf(file="03_Tumor_WES_Analysis/outputs/plots/05_KDR_lollipop_plot.pdf",height=6)
kdr_lollipop = lollipopPlot(
  maf = asc_maf,
  gene = 'KDR',
  AACol = 'Protein_Change',
  #cBioPortal=TRUE,
  #repel=TRUE,
  showMutationRate = FALSE,
  printCount=TRUE,
  labelPos="all",
  showDomainLabel=FALSE
  #refSeqID="NM_015450",
  #ref.build = "hg19"
  #refSeqID="NM_001042594"
  #labelPos = 882
)
dev.off()

## Make lollipop plot for POT1
pdf(file="03_Tumor_WES_Analysis/outputs/plots/05_POT1_lollipop_plot.pdf",height=6)
pot1_lollipop = lollipopPlot(
  maf = asc_maf,
  gene = 'POT1',
  AACol = 'Protein_Change',
  #cBioPortal=TRUE,
  #repel=TRUE,
  showMutationRate = FALSE,
  printCount=TRUE,
  labelPos="all",
  showDomainLabel=FALSE
  #refSeqID="NM_015450",
  #ref.build = "hg19"
  #refSeqID="NM_001042594"
  #labelPos = 882
)
dev.off()

mut_freq_df_mut_recurrent_only_merged



