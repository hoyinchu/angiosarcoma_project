# 00_Preprocessing

Notebooks that turn raw clinical spreadsheets and sequencing sample tables (`data/raw/`) into the tables in `data/processed/`. Inputs and outputs are on Zenodo only. No figures.

| Notebook | Output |
| --- | --- |
| `00_Process_Clinical_Data.ipynb` | Per-patient clinical table merging patient-reported and medical-record data. Released only in de-identified form (`data/processed/clinical_data_manual_deident.tsv`, also `SD1`). |
| `01_Process_Sequencing_Sample_Data.ipynb` | `sample_data.tsv`, `sample_clin_data.tsv`: one row per sequenced sample, joined with clinical data |
| `02_Process_Treatment_Data.ipynb` | `treatment_data.tsv` (also `SD4`) |
| `03_Process_Other_Cancer_Data.ipynb` | `other_cancer_data.tsv` (also `SD2`) |
