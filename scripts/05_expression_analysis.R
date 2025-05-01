###############################################################
# SCRIPT: 05_expression_analysis.R
# PURPOSE: Analysis of gene expression patterns
# AUTHOR: GIFT ASEFON
# DATE: 04-28-2025
###############################################################

# Load libraries
library(tidyverse)     # For data manipulation and visualization
library(pheatmap)      # For heatmap visualization
library(org.Hs.eg.db)  # For gene annotation
library(AnnotationDbi) # For interfacing with annotation databases

# Load preprocessed data
load("data/preprocessed_data.RData")

# Process expression data - Using tidyverse approach to avoid dimnames error
# Design choice: Converting to a tidy format first, then reconstructing the matrix
# to ensure dimensions match correctly
cat("Processing expression data...\n")

# First convert to tibble and ensure gene IDs are the first column
expr_data <- as_tibble(expression)
gene_column_name <- colnames(expr_data)[1]

# Create a tidy version of the expression data
expr_tidy <- expr_data %>%
  # Convert to long format for easier manipulation
  pivot_longer(cols = -all_of(gene_column_name), 
               names_to = "sample_id", 
               values_to = "expression")

# Convert back to matrix format properly
expr_matrix <- expr_tidy %>%
  pivot_wider(id_cols = all_of(gene_column_name), 
              names_from = sample_id, 
              values_from = expression) %>%
  column_to_rownames(gene_column_name) %>%
  as.matrix()

cat("Expression matrix created with dimensions:", paste(dim(expr_matrix), collapse=" x "), "\n")

# Check for and remove zero-variance genes
# Design choice: Removing genes with no variability across samples as they
# provide no information for differential expression or clustering analyses
var_genes <- apply(expr_matrix, 1, var, na.rm = TRUE)
zero_var <- which(var_genes == 0)
cat("Number of zero-variance genes:", length(zero_var), "out of", nrow(expr_matrix), "\n")

# Remove zero-variance genes to focus analysis on variable genes
if (length(zero_var) > 0) {
  expr_matrix_filtered <- expr_matrix[-zero_var, ]
  cat("Removed zero-variance genes. Remaining genes:", nrow(expr_matrix_filtered), "\n")
} else {
  expr_matrix_filtered <- expr_matrix
  cat("No zero-variance genes found\n")
}

# Log-transform expression values for better visualization
# Design choice: Log transformation makes the distribution more normal
# and reduces the impact of extreme values, improving visualization
expr_matrix_log <- log2(expr_matrix_filtered + 1)
cat("Applied log2 transformation to expression values\n")

# Get top variable genes for analysis and visualization
# Design choice: Focusing on the most variable genes to identify
# patterns of biological interest and reduce computational load
var_genes_log <- apply(expr_matrix_log, 1, var, na.rm = TRUE)
top_var_genes <- names(sort(var_genes_log, decreasing = TRUE))[1:50]
cat("Top 10 most variable genes:", paste(head(top_var_genes, 10), collapse=", "), "\n")

# Save the list of top variable genes for reference
write_csv(
  tibble(
    Gene = top_var_genes,
    Variance = var_genes_log[top_var_genes]
  ),
  "results/tables/top_variable_genes.csv"
)

# Create heatmap of top variable genes
# Design choice: Using hierarchical clustering to identify patient
# subgroups based on expression patterns of the most variable genes
cat("Creating heatmap of top variable genes...\n")
png("results/figures/expression/top_variable_genes_heatmap.png", 
    width = 12, height = 10, units = "in", res = 300)
pheatmap(
  expr_matrix_log[top_var_genes, ],
  # Scale rows (genes) to focus on relative expression patterns
  scale = "row", 
  # Use correlation distance for clustering genes with similar patterns
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  # Hide column names (samples) to avoid crowding
  show_colnames = FALSE,
  # Add informative labels
  main = "Expression Patterns of Top 50 Variable Genes",
  # Custom color palette from blue (low) to red (high)
  color = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100)
)
dev.off()

# Run Principal Component Analysis (PCA) on the filtered matrix
# Design choice: PCA reduces dimensionality and helps visualize
# the major sources of variation in the expression data
cat("Performing PCA analysis...\n")
pca_result <- prcomp(t(expr_matrix_log), scale = TRUE)

# Calculate variance explained by each principal component
pca_var <- summary(pca_result)$importance[2,] * 100
cat("Variance explained by PC1:", round(pca_var[1], 1), "%\n")
cat("Variance explained by PC2:", round(pca_var[2], 1), "%\n")
cat("Cumulative variance explained by PC1-2:", round(sum(pca_var[1:2]), 1), "%\n")

# Create PCA plot data
pca_data <- tibble(
  PC1 = pca_result$x[,1], 
  PC2 = pca_result$x[,2],
  Sample = colnames(expr_matrix_filtered)
)

# Create PCA plot
# Design choice: Visualizing sample relationships in reduced dimensional space
# to identify potential batch effects or natural sample groupings
cat("Creating PCA plot...\n")
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.7, color = "#0072B2") +
  theme_minimal() +
  labs(
    title = "Principal Component Analysis of Gene Expression Data",
    subtitle = paste0("PC1: ", round(pca_var[1], 1), "% variance, PC2: ", round(pca_var[2], 1), "% variance"),
    x = paste0("PC1 (", round(pca_var[1], 1), "%)"),
    y = paste0("PC2 (", round(pca_var[2], 1), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

# Save PCA plot
ggsave("results/figures/expression/expression_pca.png", 
       plot = pca_plot, 
       width = 10, height = 8, dpi = 300)

# Create a data frame for the top variable genes boxplots
# Design choice: Using a tidy data format for easier visualization
# with ggplot2, focusing on the top 15 genes for clarity
cat("Creating expression boxplots for top variable genes...\n")
top_genes_for_plot <- top_var_genes[1:15]  # Use top 15 for clearer visualization

# Convert to tidy format for plotting
expr_data_long <- expr_tidy %>%
  filter(!!sym(gene_column_name) %in% top_genes_for_plot) %>%
  # Calculate log2 transformed values
  mutate(log_expression = log2(expression + 1)) %>%
  # Ensure genes are ordered by variance
  mutate(!!gene_column_name := factor(!!sym(gene_column_name), 
                                      levels = top_genes_for_plot))

# Create and save boxplot of top variable genes
boxplot <- ggplot(expr_data_long, aes(x = !!sym(gene_column_name), 
                                      y = log_expression, 
                                      fill = !!sym(gene_column_name))) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.size = 2) +
  theme_minimal() +
  labs(
    title = "Expression of Top 15 Variable Genes in Lung Adenocarcinoma",
    subtitle = paste("n =", length(unique(expr_data_long$sample_id)), "samples"),
    x = "Gene ID",
    y = "log2(FPKM + 1)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave("results/figures/expression/top_variable_genes_expression.png", 
       plot = boxplot, 
       width = 14, height = 8, dpi = 300)

# Get gene symbols for top variable genes (if they are Entrez IDs)
# Design choice: Converting Entrez IDs to gene symbols improves readability
# and makes it easier to interpret results in the context of the literature
cat("Mapping gene IDs to gene symbols (if applicable)...\n")
if (grepl("^\\d+$", top_var_genes[1])) {
  cat("Detected numeric IDs, attempting to map to gene symbols...\n")
  top_var_entrez <- as.character(top_var_genes[1:15])
  
  # Try to map Entrez IDs to gene symbols
  gene_symbols <- mapIds(
    org.Hs.eg.db, 
    keys = top_var_entrez,
    column = "SYMBOL", 
    keytype = "ENTREZID",
    multiVals = "first"
  )
  
  # Create a table of Gene IDs and symbols
  gene_table <- tibble(
    EntrezID = top_var_entrez,
    Symbol = gene_symbols
  )
  
  # Print the mapping
  print(gene_table)
  
  # Save this for reference
  write_csv(gene_table, "results/tables/top_variable_genes_mapping.csv")
  
  cat("Gene ID to symbol mapping complete\n")
} else {
  cat("Gene IDs do not appear to be Entrez IDs, skipping symbol mapping\n")
}

# Create a summary of expression analysis results
# Design choice: Providing a concise summary of key findings to facilitate
# interpretation and inclusion in the paper
cat("Summarizing expression analysis results...\n")
sink("results/tables/expression_analysis_summary.txt")
cat("EXPRESSION ANALYSIS SUMMARY\n\n")
cat("Total genes analyzed:", nrow(expr_matrix_filtered), "\n")
cat("Total samples:", ncol(expr_matrix_filtered), "\n")
cat("Top 10 most variable genes:", paste(head(top_var_genes, 10), collapse=", "), "\n\n")
cat("Principal Component Analysis:\n")
cat("Variance explained by PC1:", round(pca_var[1], 1), "%\n")
cat("Variance explained by PC2:", round(pca_var[2], 1), "%\n")
cat("Cumulative variance explained by PC1-2:", round(sum(pca_var[1:2]), 1), "%\n")
sink()

cat("Expression analysis complete. Results saved to 'results/figures/expression/' and 'results/tables/'\n")