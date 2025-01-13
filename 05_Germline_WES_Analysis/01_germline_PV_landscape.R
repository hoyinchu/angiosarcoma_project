library(maftools)
setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional/scripts")

germline_maf_file_path = "05_Germline_WES_Analysis/outputs/tables/filtered_germline_variants_cases_subset.maf"
germline_clin_file_path = "05_Germline_WES_Analysis/outputs/tables/germline_sample_clin_df.tsv"

germline_case_maf = read.maf(germline_maf_file_path,clinicalData = germline_clin_file_path)

## Full table annotation
germline_case_table = read.csv("05_Germline_WES_Analysis/outputs/tables/filtered_germline_variants_cases.tsv",sep="\t")
subset_genes = unique(germline_case_table[germline_case_table$pass_germline_filter!="unselected_clinvar_pathogenic_high_conf",]$SYMBOL)

germline_case_maf_subset = subsetMaf(germline_case_maf,genes = subset_genes)

landscape_save_path = "05_Germline_WES_Analysis/outputs/plots/02_germline_PV_landscape.pdf"
pdf(file=landscape_save_path,height=6,width=12)
oncoplot(
  germline_case_maf_subset,
  top = 100,
  topBarData="Age",
  clinicalFeatures="PrimarySite",
  draw_titv=FALSE,
  cohortSize=229,
  bgCol="#FFFFFF",
  borderCol="#000000"
)
dev.off()
dim(germline_case_maf@clinical.data)

## Binomial exact test of # of patients carrying germline variants over 
binom.test(39,229)
germline_case_maf_subset@data$Tumor_Sample_Barcode

## Plot a dual lollipop focused on POT1
somatic_maf_path = "data/processed/tumor_WES/ASC_mutations.maf"
somatic_maf_metadata_path = "data/processed/tumor_WES/ASC_mutations_metadata.tsv"
somatic_maf = read.maf(maf=somatic_maf_path,clinicalData=somatic_maf_metadata_path)

pot1_dual_lollipop_save_path = "05_Germline_WES_Analysis/outputs/plots/02_POT1_dual_lollipop.pdf"
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


