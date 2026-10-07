# 03_Tumor_WES_Analysis

Somatic analysis of 94 QC-passing tumor WES samples.

| Script | What it does | Figures |
| --- | --- | --- |
| `00_sample_preprocessing.R` | Tumor/normal pair QC (FFPE, OxoG, contamination, tumor-in-normal) | – |
| `01_maf_preprocessing.R` | Filters the MAF to QC-passing pairs → `data/processed/tumor_WES/ASC_mutations.maf` | – |
| `02_mutational_signatures.R` | Trinucleotide spectra; COSMIC signatures (SBS7 by site) | Supp. Fig. 22 |
| `03_maf_add_metadata.R` | Adds clinical, signature and CN metadata; TMB | Fig. 3a |
| `04_gene_level_mutation_summary.R` | MutSig significance, oncoplots with CN, cutaneous vs. non-cutaneous | Fig. 3b–d |
| `05_variant_level_mutation_summary.R` | Recurrent missense variants, AlphaMissense, VAF | Fig. 3e; Supp. Fig. 23 |
| `06_cnvkit_analysis.R` | Gene-level amplification/deletion calls → `combined_cnvkit_calls.csv` | – |

Run `06` before `03` and `04`; `05` needs `05_Germline_WES_Tumor_WES_Analysis` output.

Outputs: per-sample mutation counts and representative mutations (`outputs/tables/`), recurrent genes and AlphaMissense scores (`outputs/`), plots by script number (`outputs/plots/`; `unused/` is not in the paper). `onco_matrix.txt` *(Zenodo only)* is the oncoplot matrix.

`ASC_mutations.maf` keeps all variant classes; filter to nonsynonymous before comparing with other resources.
