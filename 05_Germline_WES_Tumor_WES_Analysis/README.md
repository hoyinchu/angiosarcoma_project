# 05_Germline_WES_Tumor_WES_Analysis

Genes hit by both a germline PV and a somatic mutation. Run after `04_Germline_WES_Analysis`.

| Script | What it does | Figures |
| --- | --- | --- |
| `00_two_hit_scan.R` | Germline vs. somatic hits per gene; biallelic ranking; POT1 age of onset | Fig. 5a–c |
| `01_pick_representative_PV.R` | One representative PV per patient and gene, used in the Fig. 2d heatmap | – |

Outputs: `outputs/tables/germline_vs_tumor_gene_counts.csv` and the figure plots in `outputs/plots/`.
