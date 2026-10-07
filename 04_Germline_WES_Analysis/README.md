# 04_Germline_WES_Analysis

Germline pathogenic variants (PVs) in 229 patients and burden testing against controls.

| Script | What it does | Figures |
| --- | --- | --- |
| `00_germline_preprocessing.R` | Selects germline samples; formats filtered case variants | – |
| `01_germline_PV_landscape.R` | PV landscape; POT1 germline/somatic lollipop | Fig. 4a, 4c |
| `02_germline_gene_enrichment.R` | Case-control gene burden test (≥2 case carriers) | Fig. 4b |

Inputs in `data/processed/germline/` are on Zenodo only.

Outputs (`outputs/tables/`): filtered variants in cases and controls, burden test results, representative PV per patient and gene. The released PV table is `data/supplementary_data/SD8_Germline_Pathogenic_Variants.tsv`.
