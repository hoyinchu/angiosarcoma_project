# 02_Tumor_RNA_Analysis

Tumor RNA-seq analysis of 122 QC-passing samples.

| Script | What it does | Figures |
| --- | --- | --- |
| `00_RNA_Data_Preprocessing.R` | Sample QC, gene filtering, VST normalization | – |
| `01_Plot_Distribution.R` | Exploratory distribution plots; unused | – |
| `02_Tumor_RNA_Characterization.R` | Seurat object, hierarchical expression clusters, purity by site | Supp. Figs. 5, 7 |
| `03_Tumor_RNA_PCA.R` | PCA colored by site, cutaneous status, cluster and other variables | Fig. 2a–c; Supp. Fig. 6 |
| `04_Tumor_RNA_Enrichment.R` | DESeq2 one-vs-rest by site and by cluster (adjusted for purity and inferred sex); Hallmark/C6 enrichment | Supp. Figs. 10–15 (sites), 17–21 (clusters) |
| `05_Tumor_RNA_Heatmap.R` | Heatmap of top pathway genes by site, with mutation and germline annotations | Fig. 2d |
| `05b_Tumor_RNA_Heatmap_by_hclust.R` | Same, grouped by expression cluster | Supp. Fig. 16 |
| `06_Tumor_RNA_Receptor_Enrichments.R` | Expression of genes targeted by therapies | Supp. Fig. 8 |
| `07_Tumor_RNA_Specificity.R` | Checks DEGs against GTEx tissue-specific genes | Supp. Fig. 9 |
| `convert_gtex_scanpy_obj.ipynb` | Converts the GTEx snRNA-seq atlas (download separately) for `07` | – |

`05*` need outputs from `03_Tumor_WES_Analysis` and `05_Germline_WES_Tumor_WES_Analysis`.

Outputs:
- `outputs/DEGs/`: DESeq2 results per site (`02_DESEQ2_<site>_vs_rest.*`) and cluster (`02_hclust_Cluster<N>_vs_rest.tsv`), combined tables, pathway enrichment, Seurat markers (also `SD6`), GTEx DEGs.
- `outputs/post_analysis_metadata.csv`: per-sample metadata with cluster assignments.
- `outputs/plots/`: subfolders numbered by producing script.

Sample names in the Seurat object use dots (`RP.1447.ASCProject.P2HVYU.T1.v1.RNA.OnPrem`); match on its `sample_alias` column (`repo_sample_alias` in `data/id_mapping/`).
