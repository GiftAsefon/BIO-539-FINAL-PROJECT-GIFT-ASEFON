###############################################################
# SCRIPT: 05_expression_analysis.R
# PURPOSE: Analysis of gene expression patterns
# AUTHOR: GIFT ASEFON
# DATE: 04-28-2025
###############################################################

# Load libraries
library(pheatmap)
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Load preprocessed data
load("data/preprocessed_data.RData")

# Process expression data
gene_names <- expression[,1]  # First column contains gene names
expr_matrix <- as.matrix(expression[,-1])  # Remove gene names column
rownames(expr_matrix) <- gene_names  # Set row names to gene names

# Check for zero-variance genes
var_genes <- apply(expr_matrix, 1, var, na.rm = TRUE)
zero_var <- which(var_genes == 0)
cat("Number of zero-variance genes:", length(zero_var), "\n")

# Remove zero-variance genes
if (length(zero_var) > 0) {
  expr_matrix_filtered <- expr_matrix[-zero_var, ]
} else {
  expr_matrix_filtered <- expr_matrix
}

# Get top variable genes
var_genes <- apply(expr_matrix_filtered, 1, var, na.rm = TRUE)
top_var_genes <- names(sort(var_genes, decreasing = TRUE))[1:50]
cat("Top 10 most variable genes:", paste(head(top_var_genes, 10), collapse=", "), "\n")

# Create heatmap of top variable genes
png("results/figures/expression/top_variable_genes_heatmap.png", width = 12, height = 10, units = "in", res = 300)
pheatmap(expr_matrix[top_var_genes, ],
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         show_colnames = FALSE,
         main = "Top 50 Variable Genes Expression")
dev.off()

# Run PCA on the filtered matrix
pca_result <- prcomp(t(expr_matrix_filtered), scale = TRUE)

# Create PCA plot data
pca_data <- data.frame(PC1 = pca_result$x[,1], 
                       PC2 = pca_result$x[,2],
                       sample = colnames(expr_matrix_filtered))

# Create PCA plot
png("results/figures/expression/expression_pca.png", width = 10, height = 8, units = "in", res = 300)
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA of Gene Expression Data",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)"))
dev.off()

# Create a data frame for the top variable genes boxplots
expr_data <- data.frame()

for (gene in top_var_genes[1:15]) {  # Use top 15 for clearer visualization
  # Extract expression values
  gene_expr <- as.numeric(expr_matrix[gene, ])
  
  # Create temporary data frame
  temp_df <- data.frame(
    Gene = rep(gene, length(gene_expr)),
    Expression = gene_expr,
    Sample = colnames(expr_matrix)
  )
  
  # Add to main data frame
  expr_data <- rbind(expr_data, temp_df)
}

# Create and save boxplot of top variable genes
png("results/figures/expression/top_variable_genes_expression.png", width = 14, height = 8, units = "in", res = 300)
p <- ggplot(expr_data, aes(x = Gene, y = log2(Expression + 1), fill = Gene)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Expression of Top 15 Variable Genes in Lung Adenocarcinoma",
       x = "Gene ID",
       y = "log2(FPKM + 1)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
print(p)
dev.off()

# Get gene symbols for top variable genes (if they are Entrez IDs)
if (grepl("^\\d+$", top_var_genes[1])) {
  top_var_entrez <- as.character(top_var_genes[1:15])
  gene_symbols <- mapIds(org.Hs.eg.db, 
                         keys = top_var_entrez,
                         column = "SYMBOL", 
                         keytype = "ENTREZID",
                         multiVals = "first")
  
  # Create a table of Gene IDs and symbols
  gene_table <- data.frame(
    EntrezID = top_var_entrez,
    Symbol = gene_symbols
  )
  
  # Print the mapping
  print(gene_table)
  
  # Save this 
  write.csv(gene_table, "results/tables/top_variable_genes_mapping.csv", row.names = FALSE)
}

cat("Expression analysis complete. Results saved to 'results/figures/expression/' and 'results/tables/'\n")