# Angiosarcoma Project Code Directory

![project_logo](images/Fig1_Repo_ver.png)

Analysis code for Chu, Hollyer, Borden et al., *Patient-partnered multiomics reveals the molecular architecture of angiosarcoma*, Nat. Commun. 17, 9032 (2026). https://doi.org/10.1038/s41467-026-75810-2 (citation: [`CITATION.cff`](CITATION.cff))

## About the study

Angiosarcoma is a rare cancer of the cells lining blood and lymph vessels. Through the patient-partnered [Angiosarcoma Project](https://ascproject.org) ([Count Me In](https://joincountmein.org)), 254 patients shared survey answers, medical records and/or tumor and saliva/blood samples. Main findings:

- Subcutaneous angiosarcomas often show TGF-β and receptor tyrosine kinase signaling, with driver mutations in *KDR*, *PLCG1* and *POT1*.
- Cutaneous angiosarcomas are enriched for MYC-driven programs, UV mutational signatures, immune checkpoint gene expression, and *TP53*, *FLT4* and *BRAF* mutations.
- Inherited (germline) *POT1* pathogenic variants carry a 92.7-fold higher risk of angiosarcoma; patients with both inherited and tumor *POT1* variants develop it decades earlier.

Shared data are de-identified: patients are identified only by study IDs, and dates are given as days from diagnosis.

## Code

Run scripts from inside their own directory, in numeric order. Some scripts need outputs from another directory; each directory's README notes these, along with the figures each script produces and its outputs.

| Directory | Contents | Figures |
| --- | --- | --- |
| [`00_Preprocessing`](00_Preprocessing/) | Builds analysis-ready clinical, sample and treatment tables | – |
| [`01_Clinical_Data_Analysis`](01_Clinical_Data_Analysis/) | Cohort demographics, metastases, prior cancers, treatments | Fig. 1; Supp. Figs. 1–4 |
| [`02_Tumor_RNA_Analysis`](02_Tumor_RNA_Analysis/) | RNA-seq QC, clustering, differential expression, pathway enrichment | Fig. 2; Supp. Figs. 5–21 |
| [`03_Tumor_WES_Analysis`](03_Tumor_WES_Analysis/) | Somatic mutations, signatures, copy number | Fig. 3; Supp. Figs. 22–23 |
| [`04_Germline_WES_Analysis`](04_Germline_WES_Analysis/) | Germline pathogenic variants and burden test | Fig. 4 |
| [`05_Germline_WES_Tumor_WES_Analysis`](05_Germline_WES_Tumor_WES_Analysis/) | Germline + somatic (two-hit) scan | Fig. 5a–c |
| [`06_Tumor_RNA_WES_Integration_Analysis`](06_Tumor_RNA_WES_Integration_Analysis/) | Expression differences by mutation status | Fig. 5e–f; Supp. Fig. 25 |
| [`data`](data/) | Supplementary data, inputs, references, cBioPortal ID maps | – |
| [`util_scripts`](util_scripts/) | Shared color palettes | – |

## Data

- **Start with [`data/supplementary_data/`](data/supplementary_data/)**: the de-identified data released with the paper (`SD1`–`SD9`; their numbering differs from the paper's Supplementary Data, see [`data/README.md`](data/README.md)).
- **Missing a file?** Large files are not on GitHub. Download the full snapshot (~1.75 GB, same layout) from Zenodo: https://zenodo.org/records/20416385
- **Raw sequencing data** (WES, RNA-seq) are in dbGaP under controlled access: [phs001931](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001931.v1.p1).

## FAQ: Differences from cBioPortal ([`angs_painter_2025`](https://www.cbioportal.org/study/summary?id=angs_painter_2025))

- **IDs differ:** Legacy IDs (`ASCProject_0005`) were renamed for cBioPortal (`PRLB6M`). Convert with [`data/id_mapping/`](data/id_mapping/).
- **Counts differ:** cBioPortal is a later, provisional release; this repository applies extra QC.
- **Mutation counts differ:** `ASC_mutations.maf` includes silent and non-coding variants.
- **Copy number differs:** CNVkit-based gene calls here vs. GISTIC calls on cBioPortal.

Details: [full FAQ](data/README.md#faq-differences-from-cbioportal).

## Contact

- Code: hoyinchu2016@gmail.com
- Data access: data@ascproject.org or Saud AlDubayan (Saud_Aldubayan@DFCI.HARVARD.EDU)
