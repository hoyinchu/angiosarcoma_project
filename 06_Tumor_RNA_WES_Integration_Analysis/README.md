# 06_Tumor_RNA_WES_Integration_Analysis

`00_RNA_WES_Data_Integration.R` compares expression between altered and non-altered tumors for each gene altered in ≥5 of the 57 samples with both RNA-seq and WES (Seurat `FindMarkers`, likelihood-ratio test), then runs Hallmark enrichment.

| Figure | Output |
| --- | --- |
| Fig. 5e | `outputs/plots/POT1_volcano.png`, `outputs/tables/POT1_volcano.csv` |
| Fig. 5f | `outputs/plots/POT1_pos_degs_hallmark.pdf` *(Zenodo only)* |
| Supp. Fig. 25 | `outputs/plots/<gene>_sig_volcano*.png`, `<gene>_mut_pos_degs_fora_hallmark_plot.pdf` *(Zenodo only)* |

Fig. 5d is a schematic. Other outputs: per-gene results (`outputs/tables/by_gene/`, Zenodo only), `key_enriched_pathways_fora.csv`, `MYC_volcano.csv`, and exploratory `outputs/plots/00_*` plots not in the paper.
