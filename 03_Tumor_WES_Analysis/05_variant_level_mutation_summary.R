library(dplyr)
library(tidyr)
library(ggplot2)
library(maftools)
library(ComplexHeatmap)
library(ggtree)
library(scales)

setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts")
maf_path = "data/processed/tumor_WES/ASC_mutations.maf"
maf_meta_path = "data/processed/tumor_WES/ASC_mutations_metadata.tsv"
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
asc_maf_path = "data/processed/tumor_WES/ASC_mutations.maf"
asc_maf_metadata_path = "data/processed/tumor_WES/ASC_mutations_metadata.tsv"
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

## Subset MAF to only these entries
maf_subset = maf_merged[maf_merged$tumor_mutation_id %in% mut_freq_df_mut_recurrent_only$tumor_mutation_id,]
maf_subset_unique = maf_subset[!duplicated(maf_subset[c("STUDY ID","tumor_mutation_id_short")]),]
mut_by_site_table = as.data.frame(table(maf_subset_unique$tumor_mutation_id_short,maf_subset_unique$`PRIMARY SITE (Combined)`),check.names=FALSE)
colnames(mut_by_site_table) = c("tumor_mutation_id_short","primary_site","count")

# Set level by total count
tumor_mutation_id_short_order = mut_freq_df_mut_recurrent_only[order(-mut_freq_df_mut_recurrent_only$mut_count_unique),]$tumor_mutation_id_short
mut_by_site_table$tumor_mutation_id_short = factor(mut_by_site_table$tumor_mutation_id_short,levels=tumor_mutation_id_short_order)
recurrent_mutation_plot = ggplot(mut_by_site_table,aes(x=tumor_mutation_id_short,y=count,fill=primary_site)) +
  geom_bar(stat = "identity",position = "stack") +
  theme_minimal() +
  scale_y_continuous(breaks= pretty_breaks()) + 
  labs(x="Gene",y="Number of Patients with Mutation",fill="Primary Site") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position="top")
recurrent_mutation_plot
ggsave(recurrent_mutation_plot,file="03_Tumor_WES_Analysis/outputs/plots/05_recurrent_mutation_plot.png",height=4,width=6)
ggsave(recurrent_mutation_plot,file="03_Tumor_WES_Analysis/outputs/plots/05_recurrent_mutation_plot.pdf",height=4,width=6)

## Also calculate the AlphaMissense scores for each of the variant
## Uniprot ID file = uniprot_to_gene_names.txt
## AlphaMissenese file is grep to only include entries with uniprot ids of interest
uniprot_map = read.csv("data/public/uniprot_to_gene_names.txt",sep="\t",header=FALSE)
colnames(uniprot_map) = c("uniprot_id","id_type","gene_name")
mut_freq_df_mut_recurrent_only_merged = merge(mut_freq_df_mut_recurrent_only,
                                              uniprot_map,by.x="Hugo_Symbol",by.y="gene_name",how="left")
alpha_missense_df = read.csv("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/reference_data/AlphaMissense/AlphaMissense_hg19_subset.tsv",sep="\t")
alpha_missense_df_subset = alpha_missense_df %>% filter(
  uniprot_id %in% uniprot_genes$uniprot_id,
)
# Rename mutation ID by protein name
mut_freq_df_mut_recurrent_only$aa_change = strsplit(mut_freq_df_mut_recurrent_only$tumor_mutation_id_short,"\\.") 
mut_freq_df_mut_recurrent_only$aa_change = sapply(mut_freq_df_mut_recurrent_only$aa_change, function(x) if (length(x) >= 2) x[2] else NA)
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

#write.csv(alpha_missense_df_subset,"03_Tumor_WES_Analysis/outputs/alpha_missense_scores_all_recurrent_mut.csv",row.names = FALSE)

write.csv(uniprot_genes,"03_Tumor_WES_Analysis/outputs/uniprot_id_for_gene_with_recurrent_mutated_variants.csv",row.names = FALSE)

mut_count_alphamissense_plot = ggplot(mut_freq_df_mut_recurrent_only_with_score,aes(x=mut_count_unique,y=am_pathogenicity,label=tumor_mutation_id_short)) +
  geom_point() +
  ggrepel::geom_label_repel() +
  theme_minimal() +
  labs(x="# of Patients with Mutation",y="AlphaMissense Pathogenicity Score")
ggsave("03_Tumor_WES_Analysis/outputs/plots/05_variant_patients_by_am_pathogenicity.png",mut_count_alphamissense_plot,dpi=300,height=6,width=12)
ggsave("03_Tumor_WES_Analysis/outputs/plots/05_variant_patients_by_am_pathogenicity.pdf",mut_count_alphamissense_plot,dpi=300,height=6,width=12)

## Make a "representative somatic / germline mutation" per sample table that has somatic vs. germline mutation?
## Somatic Logic: Hotspot Mutation -> Mutation in Mutsig Significant Gene -> Mutation in COSMIC Cancer Tier 1 gene
## Germline Logic: POT1 -> Mutation in COSMIC Cancer Tier 1 gene -> Other recurrent mutations?

# Indicate whether a mutation is recurrent across patients (hotspots)
recurrent_mutations_short_ids = tumor_mutation_id_short_order
# or is in a MutSig significant gene
mutsig_output_path = "data/processed/tumor_WES/mutsig/Apr16_2024_sig_genes.txt"
mutsig_gene = read.csv(mutsig_output_path,sep="\t")
mutsig_q10_genes = mutsig_gene[mutsig_gene$q < 0.1,]$gene
# Or if it is a COSMIC tier 1 gene
cosmic_path = "data/public/cosmic_cancer_gene_census.csv"
cosmic_df = read.csv(cosmic_path)
cosmic_tier1_genes = cosmic_df[cosmic_df$Tier==1,]$Gene.Symbol
# Or if germline pathogenic variant was also detected in the gene
# Getting this table requires running a script from a later step (00_two_hit_scan.R)
# So skip this step if it's the first time running this
somatic_germline_mut_table = read.csv("06_Germline_WES_Tumor_WES_Analysis/outputs/tables/germline_vs_tumor_gene_counts.csv")
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



## Check if recurrent mutated genes are enriched in certain pathways
library(msigdbr)
fgsea_kegg_set = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c5_set = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c6_set = msigdbr(species = "Homo sapiens", category = "C6") %>% split(x = .$gene_symbol, f = .$gs_name)
fgsea_c8_set = msigdbr(species = "Homo sapiens", category = "C8") %>% split(x = .$gene_symbol, f = .$gs_name)

test_genes = c("KDR","POT1","DEAF1","CGREF1","TPO","SYPL2","SPERT","PTPRO","PHF21B","OR10AG1","NYAP2","NCKAP5","MGRN1","MALT1","FAM194B","ERN2","DSCAM","DNAI1","COL19A1","CNR1","CHD3","ACTN2")
test_pathways = fgsea::fora(fgsea_c8_set,genes=test_genes,universe=unique(maf_merged$Hugo_Symbol))
test_pathways

# ## Add MutSig information for each gene
# mutsig_output_path = "data/processed/tumor_WES/mutsig/Apr16_2024_sig_genes.txt"
# mutsig_gene = read.csv(mutsig_output_path,sep="\t")
# mut_freq_df_mut_recurrent_only_merged = merge(mut_freq_df_mut_recurrent_only,mutsig_gene,by.x="Hugo_Symbol",by.y="gene",all.x=TRUE)
# mut_freq_df_mut_recurrent_only_merged[["mutsig_pval_sig"]] = mut_freq_df_mut_recurrent_only_merged$q < 0.1
# mut_freq_df_mut_recurrent_only_merged[is.na(mut_freq_df_mut_recurrent_only_merged)] = FALSE
# mut_freq_df_mut_recurrent_only_merged = mut_freq_df_mut_recurrent_only_merged[order(-mut_freq_df_mut_recurrent_only_merged$mut_count_unique),]
# 
# recurrent_mutation_plot = ggplot(mut_freq_df_mut_recurrent_only_merged,aes(y=reorder(tumor_mutation_id_short,mut_count_unique),x=mut_count_unique,fill=mutsig_pval_sig)) +
#   geom_bar(stat="identity") +
#   labs(y="Gene",x="Number of Patients with Mutation",fill="Gene MutSig q-val < 0.1") +
#   theme_minimal()
# recurrent_mutation_plot
# ggsave(recurrent_mutation_plot,file="03_Tumor_WES_Analysis/outputs/plots/05_recurrent_mutation_plot.png")

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



