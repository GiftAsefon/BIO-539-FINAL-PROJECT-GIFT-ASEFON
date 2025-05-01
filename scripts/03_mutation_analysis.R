###############################################################
# SCRIPT: 03_mutation_analysis.R
# PURPOSE: Analysis of genetic mutations in lung adenocarcinoma
# AUTHOR: GIFT ASEFON
# DATE: 04-28-2025
###############################################################

# Load libraries
library(maftools)  # For mutation analysis and visualization
library(tidyverse) # For data manipulation and visualization

# Load preprocessed data
load("data/preprocessed_data.RData")

# Generate mutation summary plot
# Design choice: Using plotmafSummary to get a comprehensive overview of the
# mutation landscape, removing extreme outliers for better visualization
png("results/figures/mutations/mutation_summary.png", width = 10, height = 8, units = "in", res = 300)
plotmafSummary(maf, rmOutlier = TRUE, addStat = 'median', 
               textSize = 0.8, # Ensure text is readable
               showBarcodes = FALSE) # Hide individual sample IDs to avoid clutter
dev.off()

# Create oncoplot for top 20 mutated genes
# Design choice: Focusing on top 20 genes to identify the most relevant
# genetic alterations across the cohort while keeping the plot readable
png("results/figures/mutations/oncoplot_top20.png", width = 12, height = 10, units = "in", res = 300)
# Adding annotations for clinical variables would enhance this plot
oncoplot(maf, top = 20, 
         fontSize = 0.8, # Adjust font size for readability
         showTumorSampleBarcodes = FALSE) # Hide individual sample IDs to reduce clutter
dev.off()

# Lollipop plots for key lung cancer genes
# Design choice: Creating separate lollipop plots for each key gene to
# examine the distribution of mutations across protein domains
for (gene in key_genes) {
  # Check if the gene has mutations in our dataset
  if (gene %in% unique(mutations$Hugo_Symbol)) {
    png(paste0("results/figures/mutations/lollipop_", gene, ".png"), 
        width = 10, height = 6, units = "in", res = 300)
    
    # Add more detailed title and labels for better interpretation
    lollipopPlot(maf, gene = gene, 
                 showMutationRate = TRUE, # Show mutation rate in the title
                 labelPos = "all", # Label all domains
                 showDomainLabel = TRUE) # Show domain names
    
    dev.off()
    
    # Print progress to console
    cat(paste0("Created lollipop plot for ", gene, "\n"))
  } else {
    cat(paste0("No mutations found for ", gene, " in this dataset\n"))
  }
}

# Mutation co-occurrence/mutual exclusivity
# Design choice: Analyzing interactions between mutations to identify
# potential functional relationships between different gene alterations
png("results/figures/mutations/somatic_interactions.png", width = 12, height = 10, units = "in", res = 300)
# Looking at top 25 genes to capture more potential interactions
somaticInteractions(maf, top = 25, 
                    pvalue = 0.05, # Only show significant interactions
                    returnAll = FALSE) # Focus on significant interactions
dev.off()

# Create mutation frequency table using tidyverse
# Design choice: Using tidyverse for more readable data manipulation
# and to calculate percentages more efficiently
mutations_tibble <- as_tibble(mutations)
unique_samples <- n_distinct(mutations_tibble$Tumor_Sample_Barcode)

mutation_counts <- mutations_tibble %>%
  filter(Hugo_Symbol %in% key_genes) %>%
  count(Hugo_Symbol, name = "Frequency") %>%
  mutate(Percentage = round(100 * Frequency / unique_samples, 2)) %>%
  arrange(desc(Frequency)) %>%
  rename(Gene = Hugo_Symbol)

# Write the results to a CSV file
write_csv(mutation_counts, "results/tables/mutation_frequency.csv")

# Create a more detailed mutation summary for key genes
# Design choice: Adding more detailed analysis for the key genes
# to provide deeper insights into mutation patterns
key_gene_details <- tibble()

for (gene in key_genes) {
  # Count mutations by type for this gene
  gene_muts <- mutations_tibble %>%
    filter(Hugo_Symbol == gene) %>%
    count(Variant_Classification) %>%
    mutate(Gene = gene)
  
  # Add to our results table
  key_gene_details <- bind_rows(key_gene_details, gene_muts)
}

# Reshape data for a better table format
mutation_details_wide <- key_gene_details %>%
  pivot_wider(
    id_cols = Gene,
    names_from = Variant_Classification,
    values_from = n,
    values_fill = 0
  ) %>%
  # Add the total frequency from our earlier calculation
  left_join(mutation_counts, by = "Gene")

# Save the detailed mutation table
write_csv(mutation_details_wide, "results/tables/mutation_details.csv")

# Create a mutation type distribution plot
# Design choice: Visualizing mutation types to understand the
# mutational processes active in lung adenocarcinoma
png("results/figures/mutations/mutation_types.png", width = 10, height = 8, units = "in", res = 300)
# Plot the distribution of mutation types
plotmafSummary(maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()

cat("Mutation analysis complete. Results saved to 'results/figures/mutations/' and 'results/tables/'\n")