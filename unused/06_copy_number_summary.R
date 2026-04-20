library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(maftools)
library(ggsci)
library(ggbeeswarm)
library(ggpubr)
library(RColorBrewer)

# Preprocessing
gistic_outdir = "../data/processed/tumor_WES/GISTIC/Feb2024/call-tumor_gistic/"

## TODO: clustering based on significance bands only
gistic_broad_values_by_arm = read.csv(paste0(gistic_outdir,"broad_values_by_arm.txt"),sep="\t",check.names = FALSE)
rownames(gistic_broad_values_by_arm) = gistic_broad_values_by_arm$`Chromosome Arm`
gistic_broad_values_by_arm$`Chromosome Arm` = NULL
gistic_broad_values_by_arm = gistic_broad_values_by_arm[!(rownames(gistic_broad_values_by_arm) %in% c("Xq","Xp","Yq","Yp")),]
gistic_broad_values_by_arm = t(gistic_broad_values_by_arm)

gistic_significance_by_arm = read.csv(paste0(gistic_outdir,"broad_significance_results.txt"),sep="\t",check.names = FALSE)
qval_cutoff = 0.1
sig_amp_arms = gistic_significance_by_arm[(gistic_significance_by_arm$`Amp q-value` < qval_cutoff),]$Arm
sig_del_arms = gistic_significance_by_arm[(gistic_significance_by_arm$`Del q-value` < qval_cutoff),]$Arm
#keep_arms = c(sig_amp_arms,sig_del_arms)
keep_arms = gistic_significance_by_arm[(gistic_significance_by_arm$`Amp q-value` < qval_cutoff) |
                                         (gistic_significance_by_arm$`Del q-value` < qval_cutoff),]$Arm
keep_arms = keep_arms[!keep_arms %in% c("Xp","Xq","Yq","Yp")]

gistic_broad_values_by_arm_subset = gistic_broad_values_by_arm[,keep_arms]
Heatmap(t(gistic_broad_values_by_arm_subset))

asc_gistic = readGistic(gisticDir = gistic_outdir)

# Load MAF file (with GISTIC)
all_lesions =paste0(gistic_outdir,"all_lesions.conf_99.txt")
amp_genes = paste0(gistic_outdir,"amp_genes.conf_99.txt")
del_genes = paste0(gistic_outdir,"del_genes.conf_99.txt")
scores_gis = paste0(gistic_outdir,"scores.gistic")

asc_maf_path = "data/processed/tumor_WES/ASC_mutations.maf"
asc_maf_metadata_path = "data/processed/tumor_WES/ASC_mutations_metadata.tsv"
asc_maf = read.maf(
  maf=asc_maf_path,
  clinicalData=asc_maf_metadata_path
  #gisticAllLesionsFile = all_lesions,
  #gisticAmpGenesFile = amp_genes,
  #gisticDelGenesFile = del_genes,
  #gisticScoresFile = scores_gis,
)


## G-Score plot
pdf("03_Tumor_WES_Analysis/outputs/plots/06_GISTIC_GScore_plot.pdf",width = 12,height=8)
gisticChromPlot(
  gistic = asc_gistic,
  markBands = "all",
  fdrCutOff = 0.1,
  ref.build = "hg19",
  cytobandOffset = 0.05,
  #mutGenesTxtSize =1,
  mutGenes = c("KDR","POT1","FAT1","HFE","MYC","TP53","TP53BP1","TP53TG3C","TP53TG3D","MUC16","STK11"),#VAV1
  maf = asc_maf
)
dev.off()


# Load meta data from RNA analysis
rna_analysis_metadata = read.csv("02_Tumor_RNA_Analysis/outputs/post_analysis_metadata.csv",check.names = FALSE)
# Load meta data from clinical data
clin_sample_meta_path = "data/processed/sample_clin_data.tsv"
sample_meta_df = read.csv(clin_sample_meta_path,check.names = FALSE,sep="\t")
# Since some patients have multiple samples we dedup here
sample_meta_df_cnv = sample_meta_df[!duplicated(sample_meta_df$sample_alias),]
# Construct CNV metadata
sample_meta_df_cnv_subset = sample_meta_df_cnv[sample_meta_df_cnv$sample_alias %in% rownames(gistic_broad_values_by_arm),]
# Add RNA cluster ID
rna_analysis_metadata_subset = rna_analysis_metadata[,c("sample_alias","seurat_clusters_renamed")]
sample_meta_df_cnv_subset = merge(sample_meta_df_cnv_subset,rna_analysis_metadata_subset,by="sample_alias",all.x=TRUE)
sample_meta_df_cnv_subset$seurat_clusters_renamed_disc = paste0("",sample_meta_df_cnv_subset$seurat_clusters_renamed)


## Define color palettes
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
  "Other Rare Sites"="gray"
)
cutaneous_palette = c(
  "Non-cutaneous AS"=pal_npg("nrc")(5)[1],
  "Cutaneous AS"=pal_npg("nrc")(5)[2]
)
sex_clin_palette = c(
  "F"=pal_npg("nrc")(5)[1],
  "M"=pal_npg("nrc")(5)[2]
)
bx_palette = c(
  "F"=pal_npg("nrc")(5)[1],
  "M"=pal_npg("nrc")(5)[2]
)
RAAS_class_palette = c(
  "RAAS" = "#DC0000FF",
  "LAAS" = "#4DBBD5FF",
  "RAAS & LAAS" = "#00A087FF",
  "Non-RAAS/LAAS" = "#3C5488FF",
  "Unknown" = "#B09C85FF"
)

# TMB palette
orange_cols = brewer.pal(4,"Oranges")
tmb_palette = c(
  "<= 1 mut/MB" = orange_cols[1],
  "<= 3 mut/MB" = orange_cols[2],
  "<= 10 mut/MB" = orange_cols[3],
  "> 10 mut/MB" = orange_cols[4],
  "NA" = "gray"
)
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
seurart_cluster_palette =c(
  "1"=pal_npg("nrc")(5)[1],
  "2"=pal_npg("nrc")(5)[2],
  "3"=pal_npg("nrc")(5)[3],
  "4"=pal_npg("nrc")(5)[4],
  "5"=pal_npg("nrc")(5)[5],
  "NA"="gray"
)
cluster_color_map = seurart_cluster_palette
names(cluster_color_map) = paste0("Cluster ",names(seurart_cluster_palette))
## Add sample rna cluster / clinical annotation
library(circlize)
age_col_annot = colorRamp2(c(20, 80), c("white", "purple"))
asc_cnv_top_annotations = HeatmapAnnotation(
  `Primary Site`= sample_meta_df_cnv_subset$`Primary Site (Recombined)`,
  "RAAS/LAAS" = sample_meta_df_cnv_subset$RAAS_LAAS_Class,
  `Age`= sample_meta_df_cnv_subset$`Age (Combined)`,
  "Sex" = sample_meta_df_cnv_subset$`SEX (EHR_EXTRACTED)`,
  "Nuclear Grade" = sample_meta_df_cnv_subset$BX_NUCLEAR_GRADE,
  "Vasoformative" = sample_meta_df_cnv_subset$BX_VASOFORMATIVE,
  "RNA Cluster ID" = sample_meta_df_cnv_subset$seurat_clusters_renamed_disc,
  col=list(
    `Primary Site`=primary_site_palette,
    "RAAS/LAAS" = RAAS_class_palette,
    "Sex" = sex_clin_palette,
    "Age" = age_col_annot,
    "Epithelioid" = bx_palette,
    "Spindle Cell"=bx_palette,
    "Nuclear Grade" = nuclear_grade_palette,
    "RNA Cluster ID" = seurart_cluster_palette
    #"Mets at Dx" = mets_dx_palette
  )
)

ComplexHeatmap::Heatmap(t(gistic_broad_values_by_arm_subset),top_annotation=asc_cnv_top_annotations)
dev.off()
## Do it more granularly on the sub arm level
# Load gistic by gene data
gistic_all_data_by_genes = read.csv(paste0(gistic_outdir,"all_data_by_genes.txt"),sep="\t")
gistic_all_thresholded_by_genes = read.csv(paste0(gistic_outdir,"all_thresholded.by_genes.txt"),sep="\t")
gistic_all_data_by_genes=gistic_all_thresholded_by_genes
gistic_all_data_by_genes$Cytoband_subarm = sub("\\..*", "", gistic_all_data_by_genes$Cytoband)

cytoband_max = gistic_all_data_by_genes %>%
  group_by(Cytoband_subarm) %>%
  summarize(across(starts_with("ASCProject"), median, na.rm = TRUE), .groups = "drop")
cytoband_max_df= as.data.frame(cytoband_max)
rownames(cytoband_max_df) = cytoband_max_df$Cytoband_subarm
cytoband_max_df = cytoband_max_df[!grepl("X|Y", rownames(cytoband_max_df)), ]
cytoband_max_df$Cytoband_subarm = NULL

# TODO: significant peaks only?
gistic_all_lesions = read.csv(paste0(gistic_outdir,"all_lesions.txt"),sep="\t")
gistic_all_lesions_filtered = gistic_all_lesions[gistic_all_lesions$q.values < 1e-3,]
gistic_all_lesions_filtered = gistic_all_lesions_filtered[grepl("CN",gistic_all_lesions_filtered$Unique.Name),]
gistic_all_lesions_filtered = gistic_all_lesions_filtered[!grepl("X",gistic_all_lesions_filtered$Descriptor),]
gistic_all_lesions_filtered = gistic_all_lesions_filtered[!grepl("Y",gistic_all_lesions_filtered$Descriptor),]

# For deletions encode as negatives
deletion_peak_rows = grepl("Deletion", gistic_all_lesions_filtered$Unique.Name)
asc_cnv_cols = grep("^ASCProject", colnames(gistic_all_lesions_filtered), value = TRUE)
#gistic_all_lesions_filtered[deletion_peak_rows, asc_cnv_cols] = gistic_all_lesions_filtered[deletion_peak_rows, asc_cnv_cols] * -1
gistic_all_lesions_filtered$name = paste0(gistic_all_lesions_filtered$Unique.Name,": ",gistic_all_lesions_filtered$Descriptor)
gistic_all_lesions_filtered_subset = gistic_all_lesions_filtered[,c("name",asc_cnv_cols)]

rownames(gistic_all_lesions_filtered_subset) = gistic_all_lesions_filtered_subset$name
gistic_all_lesions_filtered_subset$name = NULL

gistic_score = read.csv(paste0(gistic_outdir,"scores.gistic"),sep="\t")
gistic_amp_genes = read.csv(paste0(gistic_outdir,"amp_genes.txt"),sep="\t")

cytoband_max_df_subset = cytoband_max_df[c("7p12","7q31","21q22","1q21","17q21","4q13","12q13","15q15","6q23"),]

# TODO: look at thresholded intensity instead


# gistic_all_data_by_genes_subset = gistic_all_data_by_genes[grep("1q21",gistic_all_data_by_genes$Cytoband),]
# 
# gistic_all_data_by_genes_subset = gistic_all_data_by_genes %>%
#   filter(grepl(paste(c("1q21","17q21"), collapse = "|"), Cytoband))
## TODO: collapse into subarm level
# 
# rownames(gistic_all_data_by_genes_subset) = gistic_all_data_by_genes_subset$Gene.Symbol
# gistic_all_data_by_genes_subset$Gene.Symbol = NULL
# gistic_all_data_by_genes_subset$Gene.ID = NULL
# gistic_all_data_by_genes_subset$Cytoband = NULL


## Instead of clustering order by RNA cluster
column_order_by_clust = sample_meta_df_cnv_subset[order(sample_meta_df_cnv_subset$seurat_clusters_renamed_disc),]$sample_alias

#pdf("03_Tumor_WES_Analysis/outputs/plots/06_GISTIC_peaks_cluster_broad.pdf",width = 18,height=14)
pdf("03_Tumor_WES_Analysis/outputs/plots/06_GISTIC_peaks_cluster_all_lesions.pdf",width = 18,height=14)
cnv_complex_heatmap_fig = ComplexHeatmap::Heatmap(
  #t(gistic_broad_values_by_arm),
  gistic_all_lesions_filtered_subset,
  top_annotation=asc_cnv_top_annotations,
  show_column_names=FALSE,
  #column_order= column_order_by_clust,
  column_split = sample_meta_df_cnv_subset$seurat_clusters_renamed_disc,
  cluster_column_slices = FALSE,
  column_dend_reorder = F,
  heatmap_legend_param = list(
    title = "Copy Numbers"#, at = c(-10, 0, 10) 
    #labels = c("neg_two", "zero", "pos_two")
  )
  )
draw(cnv_complex_heatmap_fig, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",merge_legend = TRUE,padding=unit(c(5,5,5,24),"mm"))
dev.off()

gistic_all_lesions_filtered_subset_t = as.data.frame(t(gistic_all_lesions_filtered_subset))
gistic_all_lesions_filtered_subset_t$sample_alias = rownames(gistic_all_lesions_filtered_subset_t)

rna_analysis_metadata_subset_merged = merge(
  rna_analysis_metadata_subset,
  gistic_all_lesions_filtered_subset_t,
  by="sample_alias",all.x=TRUE
  )
rna_analysis_metadata_subset_merged$seurat_clusters_renamed_disc = paste0("",rna_analysis_metadata_subset_merged$seurat_clusters_renamed)

colnames(rna_analysis_metadata_subset_merged) = trimws(colnames(rna_analysis_metadata_subset_merged))


#Amplification Peak 16 - CN values: 6p22.2
cn_peak_compare_plot = ggplot(rna_analysis_metadata_subset_merged,aes(x=seurat_clusters_renamed_disc,y=`Amplification Peak 32 - CN values: 17q23.3`)) +
  geom_violin(width=1.2) +
  geom_boxplot(width=0.1,outlier.shape = NA) +
  geom_quasirandom() +
  theme_minimal() + 
  labs(x="RNA Cluster",y="Gistic CN Values (17q23.3)") +
  stat_compare_means(comparisons = list(c("1","2"),c("1","3"),c("1","4"),c("1","5")),na.rm=TRUE)
cn_peak_compare_plot
ggsave("03_Tumor_WES_Analysis/outputs/plots/06_17q23.3_Cluster_by_CN_Num.png",cn_peak_compare_plot,dpi=300,height=4,width=4)

ggplot(rna_analysis_metadata_subset_merged,aes(x=seurat_clusters_renamed_disc,y=`Amplification Peak 16 - CN values: 6p22.2`)) +
  geom_point()

ggplot(rna_analysis_metadata_subset_merged,aes(x=`Amplification Peak 32 - CN values: 17q23.3`,y=`Amplification Peak 16 - CN values: 6p22.2`, color=seurat_clusters_renamed_disc)) +
  geom_point()


# ComplexHeatmap::Heatmap(cytoband_max_df_subset,top_annotation=asc_cnv_top_annotations)
#gisticBubblePlot(gistic = asc_gistic)
## Load processed MAF file and clinical metadata
# asc_maf_path = "data/processed/tumor_WES/ASC_mutations.maf"
# asc_maf_metadata_path = "data/processed/tumor_WES/ASC_mutations_metadata.tsv"
# asc_maf = read.maf(maf=asc_maf_path,clinicalData=asc_maf_metadata_path,
#                    gisticScoresFile = paste0(gistic_outdir,"scores.gistic"),
#                    gisticAllLesionsFile = paste0(gistic_outdir,"all_lesions.txt"),
#                    gisticAmpGenesFile = paste0(gistic_outdir,"amp_genes.txt"),
#                    gisticDelGenesFile = paste0(gistic_outdir,"del_genes.txt")
#                    )
# oncoplot(asc_maf)





