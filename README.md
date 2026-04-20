# Angiosarcoma Project Code Directory

Code repository for analysis performed in the angiosarcoma project manuscript. Each directory contains code related to a particular type of analysis and are numbered in order of execution unless otherwise indicated in the code files.

Unless otherwise specificied, the analysis code should be ran in ascending order by directory as well as by the prefix in each script in each directory.

- `00_Preprocessing`: Scripts used to preprocess extracted clinical data into data table for downstream analysis.
- `01_Clinical_Data_Analysis` Scripts used to compute various demographic statistics as well as further polishing of the data tables for downstream analysis
- `02_Tumor_RNA_Analysis` Scripts used for conducting Tumor RNA analysis. Includes loading of the gene count table, QC and normalization, differential gene expression and other analysis
- `03_Tumor_WES_Analysis` Scripts used for analyzing tumor WES data. Includes scripts for further QCing somatic variant, computing gene and variant level statistics 
- `04_Germline_WES_Analysis` Scripts used for performing germline WES analysis. Includes script for genomic landscape visualization and burden testing
- `05_Germline_WES_Tumor_WES_Analysis` Scripts used for performing joint Germline and Tumor Analysis. Includes scripts for performing biallelic mutation scan
- `06_Tumor_RNA_WES_Integration_Analysis` Scripts used for performing join Tumor WES and Tumor RNA analysis. Inlucdes scripts for performing differential expression conditioned on mutation status of individual genes. 

If there any question about the data or the scripts, please reach out to hoyinchu2016@gmail.com.
