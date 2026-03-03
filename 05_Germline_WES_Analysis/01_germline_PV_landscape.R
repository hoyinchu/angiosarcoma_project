library(maftools)

germline_maf_file_path = "./outputs/tables/filtered_germline_variants_cases_subset.maf"
germline_clin_file_path = "./outputs/tables/germline_sample_clin_df.tsv"

germline_case_maf = read.maf(germline_maf_file_path,clinicalData = germline_clin_file_path)

## Full table annotation
germline_case_table = read.csv("./outputs/tables/filtered_germline_variants_cases.tsv",sep="\t")
subset_genes = unique(germline_case_table[germline_case_table$pass_germline_filter!="unselected_clinvar_pathogenic_high_conf",]$SYMBOL)

germline_case_maf_subset = subsetMaf(germline_case_maf,genes = subset_genes)

## Rename column for visalization
germline_case_maf_subset@clinical.data$PrimarySite = as.character(
  germline_case_maf_subset@clinical.data$`Primary_Site_(Recombined)`
)
annotation_color_simple = list(
  PrimarySite = clean_palette
)

annotation_color = list(`Primary_Site_(Recombined)`=primary_site_palette_maftools)
landscape_save_path = "./outputs/plots/02_germline_PV_landscape.pdf"
pdf(file=landscape_save_path,height=4,width=6)
oncoplot(
  maf = germline_case_maf_subset,
  top = 100,
  topBarData = "Age",
  clinicalFeatures = "PrimarySite",  # Use the new simple name
  annotationColor = annotation_color_simple,
  color = variant_col,
  draw_titv = FALSE,
  bgCol = "#FFFFFF",
  borderCol = "#000000",
  cohortSize = 229,
  fontSize=0.7,
  titleFontSize=1,
  legendFontSize = 0.8,
  annotationFontSize=0.8
)
dev.off()

## Binomial exact test of # of patients carrying germline variants over 
binom.test(39,229)
germline_case_maf_subset@data$Tumor_Sample_Barcode

## Plot a dual lollipop focused on POT1
somatic_maf_path = "../data/processed/tumor_WES/ASC_mutations.maf"
somatic_maf_metadata_path = "../data/processed/tumor_WES/ASC_mutations_metadata.tsv"
somatic_maf = read.maf(maf=somatic_maf_path,clinicalData=somatic_maf_metadata_path)

pot1_dual_lollipop_save_path = "./outputs/plots/02_POT1_dual_lollipop.pdf"
pdf(file=pot1_dual_lollipop_save_path,height=6,width=12)
lollipopPlot2(
  m1 = somatic_maf,
  m2 = germline_case_maf_subset, gene = "POT1",
  #AACol1 = "amino_acid_change",
  #AACol2 = "amino_acid_change",
  m1_name = "Somatic",
  m2_name = "Germline",
  m1_label="all",
  m2_label="all"
  )
dev.off()


