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
tumor_mutation_id_short_order = mut_freq_df_mut_recurrent_only[order(mut_freq_df_mut_recurrent_only$mut_count_unique),]$tumor_mutation_id_short
mut_by_site_table$tumor_mutation_id_short = factor(mut_by_site_table$tumor_mutation_id_short,levels=tumor_mutation_id_short_order)
recurrent_mutation_plot = ggplot(mut_by_site_table,aes(y=tumor_mutation_id_short,x=count,fill=primary_site)) +
  geom_bar(stat = "identity",position = "stack") +
  theme_minimal() +
  scale_x_continuous(breaks= pretty_breaks()) + 
  labs(y="Gene",x="Number of Patients with Mutation",fill="Primary Site")
ggsave(recurrent_mutation_plot,file="03_Tumor_WES_Analysis/outputs/plots/05_recurrent_mutation_plot.png")


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



