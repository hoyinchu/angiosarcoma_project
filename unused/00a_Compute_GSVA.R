library(GSVA)
library(GSEABase)

## Load filtered data
setwd("/Users/hoyin/Desktop/DanaFarber/workspaces/CMI_Painter_Angiosarcoma_WES_analysis_mh_regional")
count_data = read.csv("reference_data/RNA_Seq/outputs/filtered_gene_counts.csv",row.names = 1,check.names = FALSE)
vst_data = read.csv("reference_data/RNA_Seq/outputs/filtered_gene_counts.vst.csv",row.names = 1,check.names = FALSE)
tpm_data = read.csv("reference_data/RNA_Seq/outputs/filtered_gene_tpm.csv",row.names = 1, check.names = FALSE)
meta_data = read.csv("reference_data/RNA_Seq/outputs/filtered_metadata.csv",check.names = FALSE)

## Load GO annotations for Biological Processes 
go_gmt_file_path = "reference_data/RNA_Seq/gene_sets/c5.go.bp.v2023.2.Hs.symbols.gmt"
go_bp_set = getGmt(go_gmt_file_path,geneIdType=SymbolIdentifier())

## Load Cancer related pathway sets
cancer_gmt_file_path = "reference_data/RNA_Seq/gene_sets/c6.all.v2023.2.Hs.symbols.gmt"
cancer_pathset = getGmt(cancer_gmt_file_path,geneIdType=SymbolIdentifier())

## Load Hallmark Signature Set
hallmark_gmt_file_path = "reference_data/RNA_Seq/gene_sets/h.all.v2023.2.Hs.symbols.gmt"
hallmark_pathset = getGmt(hallmark_gmt_file_path,geneIdType=SymbolIdentifier())
# 
# 
# ## Compute GSVA per-sample gene set enrichment score (broad Biological processes). Using the Poisson kcdf for count data
# gsvaPar_bp_set = gsvaParam(
#   exprData=as.matrix(count_data),
#   geneSet=go_bp_set,
#   minSize = 5,
#   maxSize = 500,
#   kcdf="Poisson"
# )
# gsva_go_bp_es = gsva(gsvaPar_bp_set, verbose=TRUE)
# write.csv(gsva.go.es,file="reference_data/RNA_Seq/outputs/gsva/gsva_ssGSEA_score_GO_BP_set.csv")
# 
# ## Compute GSVA per-sample gene set enrichment score (cancer pathways). Using the Poisson kcdf for count data
# gsvaPar_cancer_set = gsvaParam(
#   exprData=as.matrix(count_data),
#   geneSet=cancer_pathset,
#   minSize = 5,
#   maxSize = 500,
#   kcdf="Poisson"
# )
# gsva_cancer_es = gsva(gsvaPar_cancer_set, verbose=TRUE)
# write.csv(gsva_cancer_es,file="reference_data/RNA_Seq/outputs/gsva/gsva_ssGSEA_score_cancer_set.csv")

## Compute GSVA per-sample gene set enrichment score (broad Biological processes). Using the Gaussian kcdf for vst-transformed data
gsvaPar_vst_bp_set = gsvaParam(
  exprData=as.matrix(vst_data),
  geneSet=go_bp_set,
  minSize = 10,
  maxSize = 500,
  kcdf="Gaussian"
)
gsva_go_vst_bp_es = gsva(gsvaPar_vst_bp_set, verbose=TRUE)
write.csv(gsva_go_vst_bp_es,file="reference_data/RNA_Seq/outputs/gsva/gsva_ssGSEA_score_vst_GO_BP_set.csv")

## Compute GSVA per-sample gene set enrichment score (cancer pathways). Using the Gaussian kcdf for vst-transformed data
gsvaPar_vst_cancer_set = gsvaParam(
  exprData=as.matrix(vst_data),
  geneSet=cancer_pathset,
  minSize = 10,
  maxSize = 500,
  kcdf="Gaussian"
)
gsva_cancer_vst_es = gsva(gsvaPar_vst_cancer_set, verbose=TRUE)
write.csv(gsva_cancer_vst_es,file="reference_data/RNA_Seq/outputs/gsva/gsva_ssGSEA_score_vst_cancer_set.csv")

## Compute GSVA per-sample gene set enrichment score (hallmark pathways). Using the Gaussian kcdf for log transformed tpm data
gsvaPar_vst_hallmark_set = gsvaParam(
  exprData=as.matrix(vst_data),
  geneSet=hallmark_pathset,
  minSize = 10,
  maxSize = 500,
  kcdf="Gaussian"
)
gsva_hallmark_vst_es = gsva(gsvaPar_vst_hallmark_set, verbose=TRUE)
write.csv(gsva_hallmark_vst_es,file="reference_data/RNA_Seq/outputs/gsva/gsva_ssGSEA_score_vst_hallmark_set.csv")

# 
# ## Log transform TPM data to make it normal
# tpm_data_transformed = log1p(tpm_data)
# 
# ## Compute GSVA per-sample gene set enrichment score (broad Biological processes). Using the Gaussian kcdf for log transformed tpm data
# gsvaPar_tpm_bp_set = gsvaParam(
#   exprData=as.matrix(tpm_data_transformed),
#   geneSet=go_bp_set,
#   minSize = 5,
#   maxSize = 500,
#   kcdf="Gaussian"
# )
# gsva_go_tpm_bp_es = gsva(gsvaPar_tpm_bp_set, verbose=TRUE)
# write.csv(gsva_go_tpm_bp_es,file="reference_data/RNA_Seq/outputs/gsva/gsva_ssGSEA_score_log_tpm_GO_BP_set.csv")
# 
# ## Compute GSVA per-sample gene set enrichment score (cancer pathways). Using the Gaussian kcdf for log transformed tpm data
# gsvaPar_tpm_cancer_set = gsvaParam(
#   exprData=as.matrix(tpm_data_transformed),
#   geneSet=cancer_pathset,
#   minSize = 5,
#   maxSize = 500,
#   kcdf="Gaussian"
# )
# gsva_cancer_tpm_es = gsva(gsvaPar_tpm_cancer_set, verbose=TRUE)
# write.csv(gsva_cancer_tpm_es,file="reference_data/RNA_Seq/outputs/gsva/gsva_ssGSEA_score_log_tpm_cancer_set.csv")
# 
# ## Compute GSVA per-sample gene set enrichment score (hallmark pathways). Using the Gaussian kcdf for log transformed tpm data
# gsvaPar_tpm_hallmark_set = gsvaParam(
#   exprData=as.matrix(tpm_data_transformed),
#   geneSet=hallmark_pathset,
#   minSize = 5,
#   maxSize = 500,
#   kcdf="Gaussian"
# )
# gsva_hallmark_tpm_es = gsva(gsvaPar_tpm_hallmark_set, verbose=TRUE)
# write.csv(gsva_hallmark_tpm_es,file="reference_data/RNA_Seq/outputs/gsva/gsva_ssGSEA_score_log_tpm_hallmark_set.csv")


