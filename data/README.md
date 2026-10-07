# data

Items marked *Zenodo only* are in the full snapshot at https://zenodo.org/records/20416385, not on GitHub. Raw sequencing data are in dbGaP ([phs001931](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001931.v1.p1), controlled access).

## Supplementary data (start here)

De-identified data released with the paper, in `supplementary_data/`. The file numbers here do **not** match the paper's Supplementary Data numbers; use the last column.

| File | Contents | Paper |
| --- | --- | --- |
| `SD1_clinical_data_manual_deident.tsv` | Per-patient clinical data | Supp. Data 1 |
| `SD2_other_cancer_data.tsv` | Other (non-angiosarcoma) cancers | Supp. Data 2 |
| `SD3_metastatic_occurence_count_table.csv` | Metastatic site counts | – |
| `SD4_treatment_data.tsv` | Treatments | Supp. Data 3 |
| `SD5_ASCSeuratObj2025.rds` | Seurat object: tumor RNA-seq expression and sample metadata | – |
| `SD6_Seurat_all_markers_DEGs.csv` | RNA-seq marker genes | – (DESeq2 results for Supp. Data 4–6 are in `02_Tumor_RNA_Analysis/outputs/DEGs/`) |
| `SD7_Mutsig_Output.tsv` | MutSig significantly mutated genes | Supp. Data 7 |
| `SD8_Germline_Pathogenic_Variants.tsv` | Germline pathogenic variants | Supp. Data 9 |
| `SD9_Selection_of_Sarcoma_Genes.tsv` | 107 sarcoma-related genes | Supp. Data 8 |

`supplementary_table/` holds copies of `SD8` and `SD9`. Results for Supp. Data 10–11 (expression by mutation status) are in `06_Tumor_RNA_WES_Integration_Analysis/outputs/tables/by_gene/` (*Zenodo only*).

**Known issues:** in `SD8` (and its copy), the `VAF` column repeats `DP`; compute VAF from `AD` (alt / total). In `processed/tumor_WES/ASC_mutations.maf`, `NCBI_Build` reads `__UNKNOWN__`; all coordinates are GRCh37/hg19.

## Glossary

| Term | Meaning |
| --- | --- |
| HNFS | Head, neck, face and scalp |
| RAAS / LAAS | Radiation- / lymphedema-associated angiosarcoma (LAAS includes Stewart–Treves) |
| `(EHR_EXTRACTED)`, MedR | From abstracted medical records |
| `(PRD)` | Patient-reported data (survey) |
| `(Combined)` | Medical record value if available, otherwise patient-reported |
| `(Wagner2024)` | Site categories following Wagner et al., JAMA Netw. Open 2024 |
| `(Recombined)` | Site grouping used in the paper's figures |
| Blank values | Not available (e.g. blank `SEX (EHR_EXTRACTED)` for patients without abstracted medical records); shown as blank categories in `01_Clinical_Data_Analysis/outputs/tables/` |
| QC | Tumor WES: FFPE, OxoG, contamination and tumor-in-normal filters; RNA: contamination, rRNA/mRNA rate, chimeric reads, genes detected |

## Other folders

- `processed/`: analysis-ready inputs.
  - `tumor_WES/`: `ASC_mutations.maf`, `ASC_mutations_metadata.csv`, trinucleotide counts, MutSig output, `combined_cnvkit_calls.csv`. CNVkit segments and signature tables are *Zenodo only*.
  - `rna/` (counts, TPM, VST, Seurat objects), `germline/`, `qc_tables/` (QC pass/fail), `*.tsv` (clinical, sample, treatment tables): *Zenodo only*.
- `raw/`: sample/pair tables, unfiltered MAFs, CNA segments, source spreadsheets. *Zenodo only.*
- `public/`: reference resources (COSMIC, OncoKB, GTEx, ENCODE blacklist, UniProt). `AlphaMissense_hg19.tsv` and `gencode.v19.annotation.gtf` are *Zenodo only*.
- `curated/`: curated gene lists (targeted therapies, receptor families).
- `id_mapping/`: ID maps, below.

## ID mapping: this repository vs. cBioPortal

Most patients here have legacy numeric IDs (`ASCProject_0005`); cBioPortal uses de-identified IDs (`PRLB6M`). All maps share the same sequencing barcodes, so they can be joined.

| File | One row per | Key columns |
| --- | --- | --- |
| `patient_id_map.tsv` | Patient (303) | `repo_patient_id`, `cbioportal_patient_id`, `id_type` (`legacy_numeric` / `same_as_cbioportal`), `in_SD1_clinical_data`, `has_sequencing_in_repo`, `notes` |
| `sample_id_map.tsv` | Sequenced sample (529) | `repo_sample_barcode`, `repo_sample_alias`, cBioPortal sample/patient IDs, `in_cbioportal`, `data_type`, `sample_type`, `cbioportal_has_mutations`/`_rna`, `passed_qc`, `qc_criteria`, `notes` |
| `sample_sequencing_data_map.tsv` | Biopsy or normal (429) | IDs, WES and RNA barcodes side by side, tumor/normal pair and its QC, matched normal, RNA QC, `notes` |
| `patient_sequencing_data_map.tsv` | Patient (303) | `canonical_normal_wes_sample`, `analysis_tumor_normal_pairs`, `analysis_rna_samples`, and all tumor WES / normal WES / RNA samples with barcodes |

- **`passed_qc`** means the sample was used in the paper: tumor WES in a QC-passing pair (94, = `ASC_mutations.maf`), normal WES passing germline QC (229), RNA passing QC (122; 121 after removing a duplicate).
- **`in_cbioportal`**: 129 tumor WES and 157 RNA samples. Normals are never on cBioPortal. DNA and RNA from one biopsy share a cBioPortal sample ID, so `repo_sample_alias` and `cbioportal_sample_id` are not unique in `sample_id_map.tsv`; join on `repo_sample_barcode`.
- **Canonical normal**: the patient's normal that passed germline QC. At most one per patient, and it is the normal in every QC-passing tumor/normal pair. Blank for 3 patients whose only normal failed QC and 9 with no normal.
- 47 patients in `patient_id_map.tsv` are only on cBioPortal (no SD1 data or sequencing here). `ASCProject_0036_BLOOD` (sample type `Unknown`) belongs to a patient excluded from the study.
- In `patient_sequencing_data_map.tsv`, lists are `;`-separated and each `*_samples` column matches the order of its `*_barcodes` column. When a tumor has several normals and none passed QC, all are listed with `;` in `sample_sequencing_data_map.tsv`.

## FAQ: Differences from cBioPortal

cBioPortal study [`angs_painter_2025`](https://www.cbioportal.org/study/summary?id=angs_painter_2025) is a provisional release prepared separately from the paper's data freeze.

**IDs.** 226 of 254 patients in SD1 have legacy numeric IDs that differ on cBioPortal; use `id_mapping/`. Formats also differ: cBioPortal patient `P2HVYU`, sample `ASCProject_P2HVYU_T1`; here `ASCProject_P2HVYU` and barcodes like `RP-1447_ASCProject_P2HVYU_T1_v1_RNA_OnPrem` (dots instead of `-`/`_` in the Seurat object). Timepoint suffixes (`T1`, `BLOOD`, …) match. Some `*_BLOOD` samples are tumor cfDNA. A few patients have two IDs here (e.g. `ASCProject_0301` = `ASCProject_PY6FJH`); see `notes`.

**Counts.**

| | This repository | cBioPortal |
| --- | --- | --- |
| Patients with clinical data | 254 (228 on cBioPortal) | 274 |
| Tumor WES | 143 pairs, 94 pass QC | 129 samples, 81 with mutations (13 of the 94 QC-passing tumors have none there) |
| Tumor RNA-seq | 159 samples, 122 pass QC | 157 samples |
| Total samples | – | 328, incl. 140 `*_PLACEHOLDER` (clinical data only) |

Use the files here to reproduce the paper.

**Mutations.** `ASC_mutations.maf` has all variant classes (23,022 variants, incl. silent and non-coding); cBioPortal shows protein-altering and splice variants only (9,759). Filter to nonsynonymous classes before comparing; calls then largely agree for samples with mutations on cBioPortal. Both use hg19.

**Copy number.** Gene-level CNVkit calls here (`combined_cnvkit_calls.csv`: Amp, ShallowAmp, DeepDel, …; segments in `ASC_cnvkit_pair_cns.tsv`, *Zenodo only*) vs. GISTIC discrete calls on cBioPortal; not directly comparable.

**Clinical values.** cBioPortal bins age into 5-year groups and uses its own primary site labels (closest to `PRIMARY SITE (EHR_EXTRACTED)`); the paper uses `Age (Combined)` and `Primary Site (Recombined)`. Patients with only patient-reported data have blank `*_EHR_EXTRACTED` columns here; use the `(PRD)` or `(Combined)` columns.
