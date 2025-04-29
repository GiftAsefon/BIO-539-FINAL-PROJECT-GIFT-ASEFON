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
library(tidyverse)
library(maftools)
library(survival)
library(survminer)
library(ComplexHeatmap)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(pheatmap)
library(ggplot2)
library(gridExtra)

# Set working directory to project root (use relative paths from there)
# setwd("path/to/lung_adenocarcinoma_analysis")

# Read Clinical Data
clinical_patient <- read.delim("data/luad_tcga_gdc/data_clinical_patient.txt", comment.char="#")
clinical_sample <- read.delim("data/luad_tcga_gdc/data_clinical_sample.txt", comment.char="#")

# Merge sample and patient data
clinical_merged <- clinical_sample %>%
  inner_join(clinical_patient, by="PATIENT_ID")

# Read Mutation Data
mutations <- read.delim("data/luad_tcga_gdc/data_mutations.txt", 
                        comment.char = "#", 
                        stringsAsFactors = FALSE)

# Convert to MAF object for analysis with maftools
maf <- read.maf(maf = mutations)

# Read Copy Number Data
cna <- read.delim("data/luad_tcga_gdc/data_cna.txt")

# Read Gene Expression Data
expression <- read.delim("data/luad_tcga_gdc/data_mrna_seq_fpkm.txt")

# Define key genes for analysis
key_genes <- c("EGFR", "KRAS", "TP53", "STK11", "KEAP1")

# Create directories for results if they don't exist
dir.create("results/figures/clinical", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/mutations", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/survival", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/expression", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures/integrated", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Save preprocessed data for use in other scripts
save(clinical_patient, clinical_sample, clinical_merged, mutations, maf, 
     cna, expression, key_genes, file = "data/preprocessed_data.RData")

cat("Data preprocessing complete. Data saved to 'data/preprocessed_data.RData'\n")
