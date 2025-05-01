###############################################################
# SCRIPT: 01_data_preprocessing.R
# PURPOSE: Loading and preprocessing data for lung adenocarcinoma analysis
# AUTHOR: GIFT ASEFON
# DATE: 04-28-2025
###############################################################

# Install Bioconductor packages (uncomment if needed)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c("maftools", "ComplexHeatmap", "org.Hs.eg.db", "AnnotationDbi"))
# 
# # Install CRAN packages (uncomment if needed)
# install.packages(c("tidyverse", "survival", "survminer", "pheatmap", "dplyr", "ggplot2", "gridExtra"))

# Load necessary libraries
# Design choice: Using tidyverse for data manipulation and visualization
# because it provides consistent syntax and improved readability
library(tidyverse)

# Specialized bioinformatics packages for mutation analysis
library(maftools)      # For mutation analysis
library(survival)      # For survival analysis 
library(survminer)     # For survival visualization
library(ComplexHeatmap) # For complex heatmaps
library(org.Hs.eg.db)  # For gene annotation
library(AnnotationDbi) # For annotation database interface
library(pheatmap)      # For heatmap visualization
library(gridExtra)     # For arranging multiple plots

# Set working directory to project root (use relative paths from there)
# setwd("path/to/lung_adenocarcinoma_analysis")

# Read Clinical Data
# Design choice: Using read_delim from tidyverse instead of read.delim
# because it provides better defaults and returns a tibble which is more
# consistent with tidyverse operations
clinical_patient <- read_delim("data/luad_tcga_gdc/data_clinical_patient.txt", 
                               comment = "#", 
                               delim = "\t") %>%
  # Convert column names to uppercase for consistency with existing code
  rename_with(toupper)

clinical_sample <- read_delim("data/luad_tcga_gdc/data_clinical_sample.txt", 
                              comment = "#", 
                              delim = "\t") %>%
  rename_with(toupper)

# Merge sample and patient data
# Design choice: Using inner_join to ensure we only keep samples
# that have matching patient data, which ensures data integrity
clinical_merged <- clinical_sample %>%
  inner_join(clinical_patient, by = "PATIENT_ID")

# Read Mutation Data
# Design choice: For mutations, we're keeping read.delim to maintain
# compatibility with the maftools package which expects specific formats
mutations <- read.delim("data/luad_tcga_gdc/data_mutations.txt", 
                        comment.char = "#", 
                        stringsAsFactors = FALSE)

# Convert to MAF object for analysis with maftools
# This step is essential for using maftools visualization functions
maf <- read.maf(maf = mutations)

# Read Copy Number Data
cna <- read_delim("data/luad_tcga_gdc/data_cna.txt", delim = "\t")

# Read Gene Expression Data
expression <- read_delim("data/luad_tcga_gdc/data_mrna_seq_fpkm.txt", delim = "\t")

# Define key genes for analysis
# Design choice: Focusing on these specific genes because they are
# established drivers in lung adenocarcinoma based on literature
key_genes <- c("EGFR", "KRAS", "TP53", "STK11", "KEAP1")

# Create directories for results if they don't exist
# Design choice: Organizing results by analysis type for better navigation
dir.create("results/figures/clinical", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/mutations", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/survival", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/expression", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/integrated", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Save preprocessed data for use in other scripts
# Design choice: Saving as a single R object file to simplify data loading
# in subsequent analysis scripts and ensure consistency
save(clinical_patient, clinical_sample, clinical_merged, mutations, maf, 
     cna, expression, key_genes, file = "data/preprocessed_data.RData")

cat("Data preprocessing complete. Data saved to 'data/preprocessed_data.RData'\n")
