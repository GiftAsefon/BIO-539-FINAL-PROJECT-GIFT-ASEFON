###############################################################
# SCRIPT: 03_mutation_analysis.R
# PURPOSE: Analysis of genetic mutations in lung adenocarcinoma
# AUTHOR: GIFT ASEFON
# DATE: 04-28-2025
###############################################################

# Load libraries
library(maftools)
library(dplyr)

# Load preprocessed data
load("data/preprocessed_data.RData")

# Generate mutation summary plot
png("results/figures/mutations/mutation_summary.png", width = 10, height = 8, units = "in", res = 300)
plotmafSummary(maf, rmOutlier = TRUE, addStat = 'median')
dev.off()

# Create oncoplot for top 20 mutated genes
png("results/figures/mutations/oncoplot_top20.png", width = 12, height = 10, units = "in", res = 300)
oncoplot(maf, top = 20)
dev.off()

# Lollipop plots for key lung cancer genes
for (gene in key_genes) {
  if (gene %in% unique(mutations$Hugo_Symbol)) {
    png(paste0("results/figures/mutations/lollipop_", gene, ".png"), width = 10, height = 6, units = "in", res = 300)
    lollipopPlot(maf, gene = gene)
    dev.off()
  }
}

# Mutation co-occurrence/mutual exclusivity
png("results/figures/mutations/somatic_interactions.png", width = 12, height = 10, units = "in", res = 300)
somaticInteractions(maf, top = 25)
dev.off()

# Create mutation frequency table
mutation_counts <- data.frame(
  Gene = character(),
  Frequency = integer(),
  Percentage = numeric(),
  stringsAsFactors = FALSE
)

for (gene in key_genes) {
  count <- sum(mutations$Hugo_Symbol == gene, na.rm = TRUE)
  percentage <- round(100 * count / length(unique(mutations$Tumor_Sample_Barcode)), 1)
  
  mutation_counts <- rbind(mutation_counts, data.frame(
    Gene = gene,
    Frequency = count,
    Percentage = percentage,
    stringsAsFactors = FALSE
  ))
}

# Sort by frequency
mutation_counts <- mutation_counts[order(mutation_counts$Frequency, decreasing = TRUE), ]
write.csv(mutation_counts, "results/tables/mutation_frequency.csv", row.names = FALSE)

cat("Mutation analysis complete. Results saved to 'results/figures/mutations/' and 'results/tables/'\n")