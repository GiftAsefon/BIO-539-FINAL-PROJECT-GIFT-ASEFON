# Data Sources

This directory contains data for the lung adenocarcinoma analysis project.

## Data Source

The data used in this project was obtained from The Cancer Genome Atlas (TCGA) Lung Adenocarcinoma (LUAD) dataset, accessed through cBioPortal.

## Files

The original data files are located in the `luad_tcga_gdc` directory:

- `data_clinical_patient.txt`: Patient-level clinical information
- `data_clinical_sample.txt`: Sample-level clinical information
- `data_mutations.txt`: Somatic mutations
- `data_cna.txt`: Copy number alterations
- `data_mrna_seq_fpkm.txt`: Gene expression data (FPKM values)

## Preprocessed Data

For reproducibility, preprocessed data is saved as R objects:
- `preprocessed_data.RData`: Contains all loaded raw data
- `survival_data.RData`: Contains prepared survival data

## How to Download

If you need to download the original TCGA-LUAD dataset:

1. Visit [cBioPortal](https://www.cbioportal.org/)
2. Search for "TCGA Lung Adenocarcinoma"
3. Select the "TCGA PanCancer Atlas" or "TCGA Firehose Legacy" dataset
4. Navigate to the "Download" tab
5. Select the data types you need and download as TSV files
6. Place the downloaded files in the `luad_tcga_gdc` directory