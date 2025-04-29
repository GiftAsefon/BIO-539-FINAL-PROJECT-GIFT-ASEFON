# Lung Adenocarcinoma Genomic Analysis

## Research Question

What genomic alterations in lung adenocarcinoma are associated with clinical outcomes, and how do mutation profiles correlate with patient survival?

## Overview

This project analyzes genomic and clinical data from lung adenocarcinoma patients in The Cancer Genome Atlas (TCGA) dataset. The analysis includes:

1. **Clinical characteristics analysis**: Patient demographics, stage distribution, and risk factors
2. **Mutation profiling**: Identification of key driver mutations and their frequencies
3. **Survival analysis**: Correlating clinical factors and genetic alterations with survival outcomes
4. **Gene expression analysis**: Patterns of gene expression in lung adenocarcinoma
5. **Integrated analysis**: Multivariate models of prognostic factors

## Key Findings

- Demographic analysis revealed lung adenocarcinoma patients had a median age of 67 years (n=585), with slight female predominance (48% vs 41% male) and majority White racial composition, while smoking history showed most patients had 20-50 pack-years exposure.
- The most frequently mutated genes were TP53 (50%), KRAS (27%), and EGFR (13%).
- TP53 mutations were significantly associated with worse survival outcomes (p = 0.0393).
- Pathologic stage was strongly associated with survival (p < 0.001).
- Multivariate analysis identified Stage III showing the highest hazard ratio (HR = 3.02) as independent prognostic factors.

## Repository Structure

- `data/`: Contains data files and instructions for obtaining them
- `scripts/`: R scripts for each component of the analysis
  - `01_data_preprocessing.R`: Data loading and preprocessing
  - `02_clinical_analysis.R`: Demographic and clinical analysis
  - `03_mutation_analysis.R`: Mutation profiling
  - `04_survival_analysis.R`: Survival analysis by clinical factors
  - `05_expression_analysis.R`: Gene expression analysis
  - `06_integrated_analysis.R`: Integration of multiple data types
- `results/`: Generated figures and tables
  - `figures/`: Visualizations organized by analysis type
  - `tables/`: Data tables with summary statistics

## Requirements

- R version 4.0.0 or higher
- Required R packages:
  - tidyverse
  - maftools
  - survival
  - survminer
  - pheatmap
  - dplyr
  - ggplot2
  - gridExtra
  - org.Hs.eg.db
  - AnnotationDbi
  - ComplexHeatmap

## Reproduction Instructions

1. Clone this repository
2. Follow the instructions in `data/README.md` to obtain the TCGA-LUAD dataset
3. Run the scripts in the `scripts/` directory in numerical order
4. View the results in the `results/` directory

## Author

Gift Asefon

## Acknowledgments

- TCGA Research Network for providing the data
- cBioPortal for data access
- Dr Rachel Schwartz for intense Big Data Analysis Training