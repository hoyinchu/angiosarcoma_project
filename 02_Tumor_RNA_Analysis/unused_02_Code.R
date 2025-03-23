
#so = RunUMAP(so, dims = 1:10)

## Perform hclust based on PCs / vst space?
# so_embeddings= Embeddings(so,reduction = "pca")[,1:10]
# so_embeddings#so_embeddings = FetchData(so,vars=VariableFeatures(so),layer = "vst_scaled")
# so_hclust = hclust(dist(so_embeddings))
# so_dendro = as.dendrogram(so_hclust)
# so_dendro
# so_dendro_assign = cutree(so_hclust, k =3)

## Plot the ASC markers combining clinical and new features
angiogenesis_markers = c(
  #"VWF","ERG","FLI1",
  #"CD36","CD93"
  #"ITGA5","FIGF"
  #"TGFBR1","PIGF",
  #"VEGFC"
  #"KDR","PIK3CA","PGF"
  #"STAT1"
  "FLT1","NRP1","KDR","APLN","APLNR",#"PDGFA","PDGFRA",
  "PDGFB","PDGFRB","FGF1","FGFR1","CD34"#"TBX2","TBX3","NOTCH4", #"HFE",
  #,"PECAM1","PECAM-1"
)

angiogenesis_markers = head(all_markers[all_markers$cluster==1,]$gene,10)
angiogenesis_df = data.frame(row.names = angiogenesis_markers)
angiogenesis_df$`Gene Set` = "Angiogenesis"

# stem_cell_proliferation_markers = c(
#   "CD34","FERMT2"
#   #"TGFB1","TGFBR1","CCNE1","SIX2"
#   #"TERT"
#   #"FGF2","PTPRC","NANOG","WNT1","WNT3"
# )
# stem_cell_df = data.frame(row.names = stem_cell_proliferation_markers)
# stem_cell_df$`Gene Set` = "Stem Cell Proliferation"

lymphoangiogensis_markers = c(
  #"CD320",
  #"LYVE1",
  "FLT4","NRP2","PROX1","PDPN","MYC", #"PECAM1",
  "TBX1","TPX2","TIE1","MAP4K2","CTLA4"
  #"CD274"
  #"CLEC14A",#"ANGPT2",
  #"CCBE1","EPHA2","VEGFC","VASH1"
)

lymphoangiogensis_markers = head(all_markers[all_markers$cluster==2,]$gene,10)

lymphoangiogensis_df = data.frame(row.names = lymphoangiogensis_markers)
lymphoangiogensis_df$`Gene Set` = "Lymphangiogenesis"

flt3_related_markers = c(
  "FGFR2","FGFR3","FGF22","FGFBP1",
  "FLT3","FLT3LG","MYCL",
  "IL37","IL18",#"CD44","CDH1",
  "LMO1","WNT7B","FBXW7","CXCL12","DKK2",
  #"KDR",
  #"VCAM1",
  #"FGF10"
  #"KRT1","FLG",
  #"KRT5","SERPINB5","GATA3", #"CASP14","TP73","CDH3",
  #"CLEC2A","KLRF2"
  #,"KRT10"
  #"KRT1"#"KRT5"#"KRT17","KRT24",
  "KRT1","KRT2","KRT5","KRT10","FLG"
  #"KRT14","KRT15","KRT77","KRT78","KRT80"
)

flt3_related_markers = head(all_markers[all_markers$cluster==3,]$gene,10)

flt3_related_df = data.frame(row.names = flt3_related_markers)
flt3_related_df$`Gene Set` = "FGFR Enriched"

skin_epithelial_markers = c(
  #"KRT5",#"KRT14","PEX","GLI2",
  "KRT17","KRT6A","KRT6B","KRT6C",
  "EGFR",#"BIRC5","TP63",
  "KRT16","GSDMA","CTSV","ERBB3","PTPRF",
  #"IL36RN","FLG","IL36B","SOX9","FOXE1","VDR","IL36G",
  #"RUNX2","CARM1","TGFB1","ANGIO","AAMP","ANGPT1","ANGPT2","CALR","CXCL9","CXCL10",
  #"EPO","FGF1","FGF2","HOXB4","PGF","SERPINB5","TIMP1",
  "ERBB2",#"IL37",
  "CD82",#"FOXC1",
  "CST6",
  "SOX9",
  "DPP4",
  "SMAD3",
  "IGFL1","IGFL2"
  #,"MUC1"
  #"SOX9"
  #,"VIL1"
  #"KRT9","KRT16","KRT17"
)

skin_epithelial_markers = head(all_markers[all_markers$cluster==4,]$gene,10)

skin_epithelial_df = data.frame(row.names = skin_epithelial_markers)
skin_epithelial_df$`Gene Set` = "17q Amplified"

other_markers = c(
  "SELP","IL6","IL8",#"IL6R","IL10",#"CTNNB1","CD163",#"TNFRSF1A","DLL4",
  "STAT3",#"SOCS3",
  "BCL3",#"IL21",
  "PTPRD","OLR1",
  "ST6GALNAC5","CXCL1",
  "CXCL2","CXCL3",
  #"EPCAM","CD163",
  "HIF1A",
  #"KRT8","KRT18","ICAM1",
  "THBS1","PTPRO","EGR1",
  "VEGFA"
  #,"KDR"
  
  #"MMP2","MMP9","CTNNB1","CDH1","ITGA","S100A4"
  #,"VEGFB"#,"VEGFC"
  #"CD302"
  #,"KRT19","KRT7"
)

# other_markers = all_markers[grepl("KRT",all_markers$gene),]$gene
# other_markers = other_markers[!duplicated(other_markers)]

other_markers = head(all_markers[all_markers$cluster==5,]$gene,10)

other_markers_df = data.frame(row.names = other_markers)
other_markers_df$`Gene Set` = "Radioresistant"

#gene_set_df = rbind(angiogenesis_df,stem_cell_df,lymphoangiogensis_df,epitheliod_df,diff_epitheliod_df,other_markers_df)
gene_set_df = rbind(
  angiogenesis_df,
  lymphoangiogensis_df,
  flt3_related_df,
  skin_epithelial_df,
  other_markers_df
)


# Keratin expression by cutaneous status
# "KRT1","FLT3","MYCL","MYC",
# EGFR, ERBB3
kertain_gene_subset_merged = get_gene_subset_long(so,genes=c("FGFR1","FGFR2","FGFR3","FGFR4","FGFRL1","FGFR6"))
kertain_gene_subset_merged = kertain_gene_subset_merged %>% 
  mutate(cutaneous_viz = case_when(
    `CUTANEOUS AS (EHR_EXTRACTED)` == "0" ~ "Non-cutaneous AS",
    `CUTANEOUS AS (EHR_EXTRACTED)` == "1" ~ "Cutaneous AS",
    TRUE ~ "NA"
  ))
#kertain_gene_subset_merged$BX_EPITHELIOID
kertain_gene_subset_merged = kertain_gene_subset_merged[!is.na(kertain_gene_subset_merged$BX_EPITHELIOID),]
keratin_boxplot = ggplot(kertain_gene_subset_merged, 
                         aes(x = `Primary Site (Recombined)`, y = count_value)) +
  geom_violin(width = 1.2) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  facet_grid(rows = vars(gene_name)) +
  stat_compare_means(comparisons = list(c("Breast (Parenchymal)", "HNFS"),
                                        c("Breast (Parenchymal)", "Breast (Cutaneous)")
                                        #c("Cluster 1", "Cluster 3"),
                                        #c("Cluster 2", "Cluster 3"),
                                        #c("Cluster 2", "Cluster 4")
  ),
  na.rm = TRUE, label = "p.format") +
  stat_compare_means(method = "anova", label.y = 5) +
  geom_quasirandom(aes(color = `cutaneous_viz`)) +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0)) + # Adjust facet label orientation
  labs(x = "Cluster", y = "Normalized Expression Level", color = "Cluster") +
  guides(color=guide_legend(title="Cutaneous Status"))
ggsave(keratin_boxplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_keratin_expr_boxplot.png",dpi=300,width=12,height=8)


## Make biplot between marker genes
plot_biplot = function(so,gene1,gene2,color_attr) {
  gene_subset = so@assays$RNA$vst_scaled[c(gene1,gene2),]
  gene_subset = as.data.frame(t(gene_subset),check.names=FALSE)
  gene_subset$og_id = rownames(gene_subset)
  gene_subset = merge(gene_subset,so@meta.data,on="og_id",how="left")
  regression_result =  lm(as.formula(paste0(gene2, " ~ ", gene1)), data = gene_subset)
  coef = coef(regression_result)[2]
  p_value = summary(regression_result)$coefficients[2, 4]
  genes_biplot = ggplot(gene_subset,aes_string(x=gene1,y=gene2)) +
    geom_point(size=4,aes_string(color=color_attr)) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    #geom_point(size=4,aes(color=gene_subset[[color_attr]])) +
    theme_minimal() +
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) +
    annotate("text", x = Inf, y = Inf, label = paste0("Coef: ", signif(coef,3)), hjust = 1.1, vjust = 2, size = 4) +
    annotate("text", x = Inf, y = Inf, label = paste0("P-val: ", signif(p_value,3)), hjust = 1.1, vjust = 3.5, size = 4)
  #guides(fill=guide_legend(title="Cluster"))
  #geom_smooth(method='lm')
  return(genes_biplot)
}
so@meta.data$seurat_clusters_renamed_str = paste0("Cluster ",so@meta.data$seurat_clusters_renamed)

myc_ctla4_biplot = plot_biplot(so,"MYC","CTLA4","`Primary Site (Recombined)`") + guides(color=guide_legend(title="Primary Site"))
myc_pdl1_biplot = plot_biplot(so,"MYC","CD274","`Primary Site (Recombined)`") + guides(color=guide_legend(title="Primary Site")) + labs(y="CD274 (PD-L1)")
myc_pd1_biplot = plot_biplot(so,"MYC","PDCD1","`Primary Site (Recombined)`") + guides(color=guide_legend(title="Primary Site")) + labs(y="PDCD1 (PD-1)")
combined_myc_plot =  plot_grid(
  myc_ctla4_biplot,myc_pdl1_biplot,myc_pd1_biplot,align = "hv",ncol = 3,labels = c('A', 'B',"C")#bp5,bp6,bp7,bp8,align = "hv",ncol = 2
)
ggsave(combined_myc_plot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_myc_immuno_biplot.png",dpi=300,width=24,height=6)

#bp1 = plot_biplot(so,"FLT4","MYC","seurat_clusters_renamed_str") + guides(color=guide_legend(title="Cluster"))
bp1 = plot_biplot(so,"MYC","FLT4","`Primary Site (Recombined)`") + guides(color=guide_legend(title="Primary Site"))
bp2 = plot_biplot(so,"MYC","FLT4","seurat_clusters_renamed_str")+ guides(color=guide_legend(title="Cluster"))
bp3 = plot_biplot(so,"MYC","MYCL","`Primary Site (Recombined)`")+ guides(color=guide_legend(title="Primary Site"))
bp4 = plot_biplot(so,"MYC","MYCL","`seurat_clusters_renamed_str`")+ guides(color=guide_legend(title="Cluster"))
#bp5 = plot_biplot(so,"IL18","IL37","`Primary Site (Recombined)`")+ guides(color=guide_legend(title="Primary Site"))
#bp6 = plot_biplot(so,"IL18","IL37","`seurat_clusters_renamed_str`")+ guides(color=guide_legend(title="Cluster"))
#bp5 = plot_biplot(so,"IL18","IL37","seurat_clusters_renamed_str")+ guides(color=guide_legend(title="Cluster"))
#bp5 = plot_biplot(so,"FLT1","FLT4","seurat_clusters_renamed_str")+ guides(color=guide_legend(title="Cluster"))
#bp6 = plot_biplot(so,"FLT1","FLT3","seurat_clusters_renamed_str")+ guides(color=guide_legend(title="Cluster"))
#bp7 = plot_biplot(so,"KRT1","KRT17","seurat_clusters_renamed_str")+ guides(color=guide_legend(title="Cluster"))
#bp8 = plot_biplot(so,"CXCL1","FLT1","seurat_clusters_renamed_str")+ guides(color=guide_legend(title="Cluster"))

combined_biplots = plot_grid(
  bp1,bp2,bp3,bp4,align = "hv",ncol = 2,labels = c('A', 'B',"C","D")#bp5,bp6,bp7,bp8,align = "hv",ncol = 2
)
ggsave(combined_biplots,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_myc_mycl_combined_biplots.png",dpi=300,width=12,height=10)

## Make a biplot of all the vegf ligands 
## Along with their regression coefficient and p-values
vegf_receptors = c("FLT1","KDR","FLT3","FLT4","NRP1","NRP2","PDGFRA","PDGFRB")
vegf_ligands = c("VEGFA","VEGFB","VEGFC","FIGF","PIGF","PDGFA","PDGFB","PDGFC","PDGFD")
vegf_subset_expr = data.frame(t(so@assays$RNA$vst_scaled[c(vegf_receptors,vegf_ligands),]))
vegf_subset_expr$og_id = rownames(vegf_subset_expr)
vegf_subset_expr_merged = merge(vegf_subset_expr,so@meta.data,on="og_id",how="left")
vegf_subset_expr_merged

combinations = expand.grid(
  receptor = vegf_receptors,
  ligand = vegf_ligands,
  stringsAsFactors = FALSE
)

scatter_data = combinations %>%
  rowwise() %>%
  mutate(
    receptor_value = list(vegf_subset_expr_merged[[receptor]]),
    ligand_value = list(vegf_subset_expr_merged[[ligand]]),
    receptor_label = receptor,
    ligand_label = ligand,
    cluster_value = list(vegf_subset_expr_merged[["seurat_clusters_renamed_str"]])
  ) %>%
  unnest(cols = c(receptor_value, ligand_value,cluster_value))

regression_results <- scatter_data %>%
  group_by(receptor_label, ligand_label) %>%
  summarise(
    coef = coef(lm(ligand_value ~ receptor_value))[2],
    p_value = summary(lm(ligand_value ~ receptor_value))$coefficients[2, 4],
    highlight = ifelse(summary(lm(ligand_value ~ receptor_value))$coefficients[2, 4] < 0.05/72, TRUE, FALSE),
    .groups = "drop"
  )

scatter_data <- scatter_data %>%
  left_join(regression_results, by = c("receptor_label", "ligand_label"))

vegf_biplots = ggplot(scatter_data, aes(x = receptor_value, y = ligand_value, color=cluster_value)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  facet_grid(rows = vars(receptor_label), cols = vars(ligand_label), switch = "both") +
  labs(
    x = "Ligand Expression",
    y = "Receptor Expression"
  ) +
  theme_minimal() +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) +
  theme(
    strip.text = element_text(size = 8),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 10),
    panel.spacing = unit(1, "lines")
  ) +
  geom_text(
    data = regression_results,
    aes(
      x = Inf, y = Inf,
      label = paste0("Coef: ", signif(coef, 3), "\nP: ", signif(p_value, 3))
    ),
    inherit.aes = FALSE,
    hjust = 1.1, vjust = 1.1,
    size = 3
  ) +
  geom_rect(
    data = regression_results %>% filter(highlight),  # Highlighted panels
    aes(
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ),
    inherit.aes = FALSE,
    color = "red",
    fill = NA,
    size = 0.8
  )

ggsave(vegf_biplots,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_VEGF_biplots.png",dpi=300,width=24,height=16)




plot_biplot(so,"FLT1","NRP1","seurat_clusters")
plot_biplot(so,"KRT1","KRT2","seurat_clusters")
plot_biplot(so,"FLT4","CTSV","seurat_clusters")


plot_biplot(so,"MYC","MYCL","PRIMARY SITE (Combined)")#"seurat_clusters")
plot_biplot(so,"NRP1","NRP2","seurat_clusters")
plot_biplot(so,"FLT1","FLT4","seurat_clusters")
plot_biplot(so,"FLT1","FLG","PRIMARY SITE (Combined)")
plot_biplot(so,"FLT3","MYCL","seurat_clusters")
plot_biplot(so,"IL18","MYC","seurat_clusters")
plot_biplot(so,"CXCL1","CXCL12","seurat_clusters")


## NRP1-NRP2
nrp1_nrp2_biplot_by_site = plot_biplot(so,"NRP1","NRP2","PRIMARY SITE (Combined)")
nrp1_nrp2_biplot_by_site

## NRP1-MYC axis separates major clusters
nrp_myc_biplot_by_site = plot_biplot(so,"NRP1","MYC","PRIMARY SITE (Combined)")
ggsave(nrp_myc_biplot_by_site,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_NRP1_MYC_biplot.png",dpi=300,width=24,height=16)

## TERT-FLT1 gets parenchymal breast samples
tert_flt1_biplot_by_site = plot_biplot(so,"TERT","FLT1","PRIMARY SITE (Combined)")
ggsave(tert_flt1_biplot_by_site,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_TERT_FLT1_biplot.png",dpi=300,width=24,height=16)

## MYC-CTLA4 gets cutaneous breast samples
myc_ctla4_biplot = plot_biplot(so,"CTLA4","MYC","PRIMARY SITE (Combined)")
ggsave(myc_ctla4_biplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_MYC_CTLA4_biplot.png",dpi=300,width=24,height=16)

## FGFR3?
plot_biplot(so,"FGFR3","APOBEC3G","PRIMARY SITE (Combined)")
plot_biplot(so,"FGFR3","PRKCA","PRIMARY SITE (Combined)")


plot_biplot(so,"SEMA4D","MAPK1","PRIMARY SITE (Combined)")
plot_biplot(so,"MAPK1","TBX1","PRIMARY SITE (Combined)")
flt1_foxp3_biplot = plot_biplot(so,"FLT1","FOXP3","PRIMARY SITE (Combined)")
ctla4_foxp3_biplot = plot_biplot(so,"CTLA4","FOXP3","PRIMARY SITE (Combined)")
fgfr3_foxp3_biplot = plot_biplot(so,"FGFR3","FOXP3","PRIMARY SITE (Combined)")

ggsave(flt1_foxp3_biplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_FLT1_FOXP3_biplot.png",dpi=300,width=24,height=16)
ggsave(ctla4_foxp3_biplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_CTLA4_FOXP3_biplot.png",dpi=300,width=24,height=16)
ggsave(fgfr3_foxp3_biplot,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_FGFR3_FOXP3_biplot.png",dpi=300,width=24,height=16)


VlnPlot(so,features = c(
  "VEGFA","VEGFB","VEGFC","VEGFD",
  "PDGFA","PDGFB","PDGFC",
  "FLT1","KDR","FLT3","FLT4"
  #"TIE1","MYC",
  #"DAXX","ATRX","ARID1A","POT1","TERT"
),slot="vst_scaled"
)

so

curated_heatmap = DoHeatmap(so,features=c(
  "FLT1","NRP1","FGFR1","APLNR","PDGFA","PDGFB","TERT","POT1","ATRX","DAXX",
  "KDR","NRP2","TBX1",
  "TIE1","FLT4","MYC",
  "RAC1","CTLA4","BCR","TGFB1",#"IL2RA","ATR","CHEK1", 
  "CHEK2", "MAD2L1", "BUB1B", "TP53", "BRCA1", "BRCA2", "RAD51", "NBN","MDM2","SMYD3","FAT1",
  "KRAS","HRAS","NRAS",
  "ERBB2","ERBB3","EGFR","MET","ALK","RET","PDGFRA","KIT","HLA-A","HLA-B","HLA-C",
  "FGFR2","FGFR3","FGFR4","LRP5","ANGPT2","TEK",
  "MAPK3","MAPK1",
  #"ERK1","ERK2","MAPK","LRP6","ANG2","TIE2",
  "FLT3","PDGFC",
  "PIK3CA","PTPRB","PTPRD","PLCG1",
  "VEGFA","VEGFB","VEGFC",
  "PDCD1","FOXP3","CXCL10","CCR5",
  "CD19","MS4A1","IL10","LAG3","CCL2","IL17A","TNF"
),slot = "vst_scaled",disp.min=-2.5,disp.max=2.5)
curated_heatmap
ggsave(curated_heatmap,filename = "02_Tumor_RNA_Analysis/outputs/plots/02_seurat_plots/02_curated_genes_heatmap.png",dpi=300,width=24,height=16)


## MYC-CTLA4 axis
plot_biplot(so,"KRAS","HRAS","PRIMARY SITE (Combined)")
plot_biplot(so,"MYC","MAPK1","PRIMARY SITE (Combined)")

plot_biplot(so,"KRAS","PERP","PRIMARY SITE (Combined)")

plot_biplot(so,"DLL4","KDR","PRIMARY SITE (Combined)")
plot_biplot(so,"DLL4","FLT1","PRIMARY SITE (Combined)")
plot_biplot(so,"DLL4","NOTCH1","PRIMARY SITE (Combined)")

plot_biplot(so,"VEGFA","FLT1","PRIMARY SITE (Combined)")


## MYC-TERT axis
myc_tert_biplot = plot_biplot(so,"TERT","KDR","PRIMARY SITE (Combined)")
myc_tert_biplot

## MYC-CTLA4-KRAS axis
myc_kras_biplot_by_site = plot_biplot(so,"MYC","KRAS","PRIMARY SITE (Combined)")
ctla4_myc_biplot_by_site = plot_biplot(so,"MYC","CTLA4","PRIMARY SITE (Combined)")
plot_grid(myc_kras_biplot_by_site,ctla4_myc_biplot_by_site)

tp53_smyd2_biplot_by_site = plot_biplot(so,"UGT2A2","PNPLA5","PRIMARY SITE (Combined)")
tp53_smyd2_biplot_by_site

## CTLA4-CLEC2A (cutaneous vs, rest)
tp53_smyd2_biplot_by_site = plot_biplot(so,"CTLA4","CD274","PRIMARY SITE (Combined)")
tp53_smyd2_biplot_by_site 

plot_biplot(so,"CDHR1","CHP2","RAAS_LAAS_Class")
plot_biplot(so,"NRP1","NRP2","RAAS_LAAS_Class")

## Confirm relation between VEGF receptors (FLT/NRP)
nrp_biplot_by_site = plot_biplot(so,"NRP1","NRP2","PRIMARY SITE (Combined)")
nrp_biplot_by_raas = plot_biplot(so,"NRP1","NRP2","RAAS (EHR_EXTRACTED)")
plot_grid(nrp_biplot_by_site,nrp_biplot_by_raas)

flt_biplot_by_site = plot_biplot(so,"TERT","MYC","PRIMARY SITE (Combined)")
flt_biplot_by_site

nrp_subset = so@assays$RNA$vst_scaled[c("NRP1","NRP2","VEGFA","KDR"),]
nrp_subset = as.data.frame(t(nrp_subset),check.names=FALSE)
nrp_subset$og_id = rownames(nrp_subset)
nrp_subset = merge(nrp_subset,so@meta.data,on="og_id",how="left")
color_attr = "PRIMARY SITE (Combined)"

nrp_biplot = ggplot(nrp_subset,aes(x=NRP1,y=NRP2)) +
  geom_point(size=3,aes(color=nrp_subset[[color_attr]])) +
  theme_minimal() +
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) +
  geom_smooth(method='lm')
nrp_biplot

## Write the seurat derived metadata
write.table(so@meta.data,"data/processed/rna/02_seurat_processed_metadata.tsv",sep="\t")


## Show mutations overlaid on top of PCA
mutation_dimplot = DimPlot(so, reduction = "pca", group.by=c(
  "POT1_Mut_Germline","POT1_Mut_Onehot","TERT_Mut_Onehot","ATRX_Mut_Onehot",
  "FLT1_Mut_Onehot","KDR_Mut_Onehot","FLT3_Mut_Onehot","FLT4_Mut_Onehot",
  "POLE_Mut_Onehot","TP53_Mut_Onehot","PIK3CA_Mut_Onehot",
  "PLCG1_Mut_Onehot","PTPRB_Mut_Onehot",
  "KRAS_Mut_Onehot","NRAS_Mut_Onehot","HRAS_Mut_Onehot"
),pt.size=3)
mutation_dimplot
ggsave(mutation_dimplot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_mutation_status_on_pcs.png",dpi=300,width=24,height=12)


## Visualize the top enriched genes in each cluster
top_markers_plot = VlnPlot(so, features = c(
  "FLT1","APLNR","EBF3","NRP1",
  "MAP4K2","SMYD2","VAV3","NRP2",
  "KRTAP24-1","CCDC172","UGT2A2","PNPLA5",
  "ST6GALNAC5","C7","NCR3LG1","PTGS2",
  "KLRF2","LCE5A","CLEC2A","SP8"
  #"CARD18","SERPINB12","PGLYRP3","SEC14L1","IQCJ-SCHIP1","TSPAN18",
  #"ZNF462","BMPR1A","GLI3","TRIB3","CDC25B","DUSP5"
),pt.size = 2,ncol = 4)
top_markers_plot
ggsave(top_markers_plot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_top_markers_violin_plot.png",dpi=300,width=24,height=12)

## The same content but on PCs
top_markers_dimplot = FeaturePlot(so, features = c(
  "FLT1","APLNR","EBF3","NRP1",
  "MAP4K2","SMYD2","VAV3","NRP2",
  "KRTAP24-1","CCDC172","UGT2A2","PNPLA5",
  "ST6GALNAC5","C7","NCR3LG1","PTGS2",
  "KLRF2","LCE5A","CLEC2A","SP8"
  #"CARD18","SERPINB12","PGLYRP3","SEC14L1","IQCJ-SCHIP1","TSPAN18",
  #"ZNF462","BMPR1A","GLI3","TRIB3","CDC25B","DUSP5"
),pt.size = 2,ncol = 4,reduction = "pca")
top_markers_dimplot
ggsave(top_markers_dimplot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_top_markers_pca_plot.png",dpi=300,width=24,height=12)

## Show Continuous Feature overlaid on top of PCA
continuous_clinical_dimplot = FeaturePlot(so, reduction = "pca", features=c("Age (Combined)","tmb"),pt.size=3)
continuous_clinical_dimplot
ggsave(continuous_clinical_dimplot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_continuous_clinical_feature_pca_plot.png",dpi=300,width=24,height=12)

## Show gene expression overlaid on top of PCA
genes_to_extract = c("MYC","TP53","PLCG1","PTPRB","POT1","HRAS","NRAS","KRAS","TERT","PIK3CA")#,"PIK3CA","TERT")  # Replace with your genes of interest
expr_dim_plot = FeaturePlot(so,reduction = "pca", features=c(
  "FLT1","KDR","FLT4","FLT3",
  "PGF","PIGF","YAP1","MIB2",
  "MYC","TP53","PLCG1","PTPRB","POT1","RAD51",
  "HRAS","NRAS","KRAS","TERT","PIK3CA","GLI3","FAM83G","NFKB1"
),slot="vst_scaled",pt.size=2)#,cols = c("blue","red"))
expr_dim_plot
ggsave(expr_dim_plot,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_gene_expression_pca_plot.png",dpi=300,width=24,height=12)

## Combine all previous information and combine them into heatmap
top_cluster_markers = all_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.3) %>%
  #dplyr::filter(!grepl("KRT",gene)) %>%
  slice_head(n = 20) %>%
  ungroup()

top_cluster_markers_heatmap = DoHeatmap(so, features = top_cluster_markers$gene,slot = "vst.scaled",disp.min=-3,disp.max=3)
top_cluster_markers_heatmap
ggsave(top_cluster_markers_heatmap,filename = "reference_data/RNA_Seq/outputs/plots/02a_Seurat_top_cluster_marker_heatmap.png",dpi=300,width=24,height=12)

## Save the seurat derived clusters as metadata
#write.csv("seurat")
write.csv(so@meta.data,"reference_data/RNA_Seq/outputs/seurat/seurat_processed_metadata.csv")