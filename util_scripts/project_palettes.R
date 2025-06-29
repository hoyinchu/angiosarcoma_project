library(ggsci)
library(RColorBrewer)
library(circlize)


primary_site_palette = c(
  "Other Visceral Organs"=pal_npg("nrc")(9)[1],
  "Hepatobiliary"=pal_npg("nrc")(9)[2],
  "Breast (Cutaneous)"=pal_npg("nrc")(9)[3],
  "HNFS"=pal_npg("nrc")(9)[4],
  "Extremities"=pal_npg("nrc")(9)[5],
  "Heart"=pal_npg("nrc")(9)[6],
  "Musculoskeletal"=pal_npg("nrc")(9)[7],
  "Breast (Parenchymal)"=pal_npg("nrc")(9)[8],
  "NA"=pal_npg("nrc")(9)[9],
  "Other Rare Sites"=pal_npg("nrc")(9)[9]
)

cluster_by_primary_site_palette = c(
  "Others"=pal_npg("nrc")(9)[1],
  "Breast (Cutaneous)"=pal_npg("nrc")(9)[3],
  "HNFS"=pal_npg("nrc")(9)[4],
  "Extremities"=pal_npg("nrc")(9)[5],
  "Heart"=pal_npg("nrc")(9)[6],
  "Breast (Parenchymal)"=pal_npg("nrc")(9)[8],
  "NA"=pal_npg("nrc")(9)[9]
)

seurart_cluster_palette =c(
  "1"=pal_npg("nrc")(5)[1],
  "2"=pal_npg("nrc")(5)[2],
  "3"=pal_npg("nrc")(5)[3],
  "4"=pal_npg("nrc")(5)[4],
  "5"=pal_npg("nrc")(5)[5]
)

seurart_cluster_by_expression_palette =c(
  "Cluster 1"=pal_npg("nrc")(5)[1],
  "Cluster 2"=pal_npg("nrc")(5)[2],
  "Cluster 3"=pal_npg("nrc")(5)[3],
  "Cluster 4"=pal_npg("nrc")(5)[4],
  "Cluster 5"=pal_npg("nrc")(5)[5]
)

cutaneous_palette = c(
  "Non-cutaneous AS"=pal_npg("nrc")(5)[1],
  "Cutaneous AS"=pal_npg("nrc")(5)[2],
  "Unknown"="gray80"
)

sex_clin_palette = c(
  "Female"=pal_npg("nrc")(5)[1],
  "Male"=pal_npg("nrc")(5)[2]
)
bx_palette = c(
  "YES" = "red",
  "NO" = "lightpink"
)
spindle_cell_palette = c(
  "YES" = "red",
  "NO" = "lightpink"
)
vasoformative_cell_palette = c(
  "YES" = "red",
  "NO" = "lightpink"
)
epitheliod_palette = c(
  "YES" = "red",
  "NO" = "lightpink"
)

RAAS_class_palette = c(
  "RAAS" = "#DC0000FF",
  "LAAS" = "#4DBBD5FF",
  "RAAS & LAAS" = "#00A087FF",
  "Non-RAAS/LAAS" = "#3C5488FF",
  "Unknown" = "#B09C85FF"
)

# gene_set_palette = c(
#   "Angiogenesis"=pal_npg("nrc")(5)[1],
#   "Lymphangiogenesis"=pal_npg("nrc")(5)[2],
#   "FGFR Enriched"=pal_npg("nrc")(5)[3],
#   "17q Amplified"=pal_npg("nrc")(5)[4],
#   "Radioresistant"=pal_npg("nrc")(5)[5]
# )

# TMB palette
orange_cols = brewer.pal(4,"Oranges")
tmb_palette = c(
  "<= 1 mut/MB" = orange_cols[1],
  "<= 3 mut/MB" = orange_cols[2],
  "<= 10 mut/MB" = orange_cols[3],
  "> 10 mut/MB" = orange_cols[4],
  "NA" = "gray"
)

## Annotate Fold Change
FC_col_annot = colorRamp2(c(0, 5), c("white", "red"))
#logp_col_annot = colorRamp2(c(0, 5), c("white", "purple"))
purples_cols = brewer.pal(4,"Purples")
significance_pal = c(
  "adj. p < 0.01" = purples_cols[4],
  "adj. p < 0.1" = purples_cols[3],
  "nom. p < 0.05" = purples_cols[2],
  "nom. p >= 0.05" = purples_cols[1]
)

## Only highlight genes that are interesting (somatic)
rep_som_mut_palette = c(
  "POT1"=pal_npg("nrc")(10)[1],
  "TP53"=pal_npg("nrc")(10)[2],
  "CFTR"=pal_npg("nrc")(10)[3],
  "KDR"=pal_npg("nrc")(10)[4],
  "PLCG1"=pal_npg("nrc")(10)[5],
  #"FLG"=pal_npg("nrc")(10)[3],
  #"FLT1"=pal_npg("nrc")(10)[1],
  #"FLT3"=pal_npg("nrc")(10)[3],
  #"FLT4"=pal_npg("nrc")(10)[4],
  #"FLG"=pal_npg("nrc")(10)[8],
  #"BRAF"=pal_npg("nrc")(10)[9],
  "Others"="black",
  "Not available"="gray"
)

rep_germ_var_palette = c(
  "POT1"=pal_npg("nrc")(10)[1],
  "TP53"=pal_npg("nrc")(10)[2],
  "CFTR"=pal_npg("nrc")(10)[3],
  "BRCA2"=pal_npg("nrc")(10)[6],
  "CHEK2"=pal_npg("nrc")(10)[7],
  #"FLG"=pal_npg("nrc")(10)[3],
  #"BRCA1"=pal_npg("nrc")(10)[7],
  #"MUTYH"=pal_npg("nrc")(10)[10],
  #"PKHD1"="red",
  #"USH2A"="pink",
  #"PAH"=pal_npg("nrc")(10)[11],
  #"GJB2"=pal_npg("nrc")(10)[12],
  "Others"="black",
  "No PV Detected"="lightblue"
)

nuclear_grade_palette = c(
  "HIGH" = "orange",
  "LOW" = "lightpink"
)

mets_dx_palette = c(
  "YES" = "darkblue",
  "NO" = "lightblue"
)
