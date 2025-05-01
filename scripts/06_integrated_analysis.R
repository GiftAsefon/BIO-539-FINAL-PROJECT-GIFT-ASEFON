###############################################################
# SCRIPT: 06_integrated_analysis.R
# PURPOSE: Integration of mutation, expression, and survival data
# AUTHOR: GIFT ASEFON
# DATE: 04-28-2025
###############################################################

# Load libraries
library(survival)      # For survival analysis
library(survminer)     # For survival visualization
library(tidyverse)     # For data manipulation and visualization
library(ggplot2)       # For advanced plotting
library(grid)          # For grid graphics
library(gridExtra)     # For arranging multiple plots

# Load preprocessed data
load("data/preprocessed_data.RData")
load("data/survival_data.RData")

# Function to convert sample ID to patient ID
# Design choice: Creating a helper function to standardize the conversion
# from sample IDs to patient IDs, ensuring consistency across analyses
sample_to_patient <- function(sample_id) {
  # Extract the patient ID portion (before the first hyphen)
  # This handles TCGA barcodes where patient IDs are the first part
  return(gsub("-01.*$", "", sample_id))
}

# Create mutation status data frame for key genes
# Design choice: Combining mutation data with survival data to analyze
# the impact of specific mutations on patient outcomes
cat("Creating mutation status data frame...\n")

# Initialize data frame with patient IDs
mutation_status <- tibble(PATIENT_ID = unique(clinical_patient$PATIENT_ID))

# For each key gene, determine which patients have mutations
for (gene in key_genes) {
  # Get samples with mutations in this gene
  mutated_samples <- mutations$Tumor_Sample_Barcode[mutations$Hugo_Symbol == gene]
  
  # Convert sample IDs to patient IDs
  mutated_patients <- unique(sapply(mutated_samples, sample_to_patient))
  
  # Count mutations for this gene
  cat(paste0("Found ", length(mutated_patients), " patients with ", gene, " mutations\n"))
  
  # Add mutation status to data frame (1 = mutated, 0 = wild-type)
  mutation_status <- mutation_status %>%
    mutate(!!gene := ifelse(PATIENT_ID %in% mutated_patients, 1, 0))
}

# Merge with survival data
# Design choice: Using inner_join to ensure we only analyze patients
# with both mutation and survival data available
surv_mutations <- inner_join(surv_data, mutation_status, by = "PATIENT_ID")
cat("Created survival dataset with mutation status for", nrow(surv_mutations), "patients\n")

# Create Kaplan-Meier plots for each key gene
# Design choice: Analyzing each mutation separately to understand
# its individual impact on patient survival
for (gene in key_genes) {
  # Check if we have enough mutations to analyze
  mutations_count <- sum(surv_mutations[[gene]], na.rm = TRUE)
  cat(paste0("Analyzing survival by ", gene, " mutation (", mutations_count, " mutations)\n"))
  
  if (mutations_count >= 5) {
    # Create a clean dataset for this analysis
    gene_data <- surv_mutations %>%
      # Create a factor for the mutation status with informative labels
      mutate(
        mutation_status = factor(
          !!sym(gene),
          levels = c(0, 1),
          labels = c("Wild-type", "Mutant")
        )
      )
    
    # Create survival formula and fit
    formula <- as.formula(paste0("Surv(time, status) ~ mutation_status"))
    
    tryCatch({
      fit_gene <- survfit(formula, data = gene_data)
      
      # Log-rank test for significance
      log_rank <- survdiff(formula, data = gene_data)
      p_value <- 1 - pchisq(log_rank$chisq, df = 1)
      cat(paste0("Log-rank test for ", gene, ": p-value = ", round(p_value, 4), "\n"))
      
      # Create plot using survminer for better visualization
      tryCatch({
        gene_plot <- ggsurvplot(
          fit_gene,
          data = gene_data,
          # Add risk table and p-value
          risk.table = TRUE,
          pval = TRUE,
          # Add confidence intervals
          conf.int = TRUE,
          # Use a professional color palette
          palette = c("#2166AC", "#B2182B"),
          # Improve axis labels
          xlab = "Time (months)",
          ylab = "Overall Survival Probability",
          # Add informative title
          title = paste0("Survival by ", gene, " Mutation Status"),
          # Use a clean theme
          ggtheme = theme_minimal(),
          # Format legend
          legend.labs = c("Wild-type", "Mutant"),
          legend.title = paste(gene, "Status"),
          # Format risk table
          risk.table.height = 0.25,
          tables.theme = theme_cleantable()
        )
        
        # Save the plot
        ggsave(
          paste0("results/figures/integrated/survival_by_",   gene, "_mutation.png"),
          width = 10, height = 8, units = "in", res = 300
        )
        plot(fit_gene, 
             col = c("#2166AC", "#B2182B"),
             xlab = "Time (months)",
             ylab = "Survival Probability",
             main = paste0("Survival by ", gene, " Mutation Status"),
             lwd = 2)
        
        # Add legend
        legend("topright", 
               legend = c("Wild-type", "Mutant"),
               col = c("#2166AC", "#B2182B"),
               lty = 1,
               lwd = 2)
        
        # Add p-value
        text(x = max(fit_gene$time) * 0.7, y = 0.2, 
             labels = paste("p =", round(p_value, 4)),
             cex = 1.2)
        
        dev.off()
      })
      
      # Save p-value to a text file for reference
      cat(paste0("Log-rank test for ", gene, ": p-value = ", round(p_value, 4), "\n"), 
          file = "results/tables/mutation_survival_pvalues.txt", 
          append = TRUE, sep = "")
      
    }, error = function(e) {
      cat(paste0("Error creating survival analysis for ", gene, ": ", e$message, "\n"))
    })
  } else {
    cat(paste0("Skipping ", gene, " - not enough mutations (", mutations_count, ")\n"))
  }
}

# Create more comprehensive dataset for Cox analysis
# Design choice: Incorporating clinical and molecular data in a multivariate model
# to identify independent prognostic factors while controlling for confounders
cat("\nPreparing data for Cox proportional hazards model...\n")

# Start with mutation and survival data
cox_data <- surv_mutations

# Add stage information
# Design choice: Pathologic stage is a critical prognostic factor to include in multivariate analysis
cox_data <- cox_data %>%
  left_join(
    select(clinical_patient, PATIENT_ID, PATH_STAGE),
    by = "PATIENT_ID"
  )

# Add gender information
cox_data <- cox_data %>%
  left_join(
    select(clinical_patient, PATIENT_ID, SEX),
    by = "PATIENT_ID"
  )

# Add smoking information
cox_data <- cox_data %>%
  left_join(
    select(clinical_patient, PATIENT_ID, SMOKING_PACK_YEARS),
    by = "PATIENT_ID"
  )

# Simplify stage for analysis
# Design choice: Consolidating stage information into major categories
# to ensure sufficient group sizes for statistical power
cox_data <- cox_data %>%
  mutate(
    STAGE_SIMPLE = case_when(
      grepl("I", PATH_STAGE) ~ "Stage I",
      grepl("II", PATH_STAGE) ~ "Stage II",
      grepl("III", PATH_STAGE) ~ "Stage III",
      grepl("IV", PATH_STAGE) ~ "Stage IV",
      TRUE ~ NA_character_
    )
  ) %>%
  # Ensure proper ordering of stages
  mutate(
    STAGE_SIMPLE = factor(
      STAGE_SIMPLE, 
      levels = c("Stage I", "Stage II", "Stage III", "Stage IV")
    )
  )

# Create a categorical variable for smoking
# Design choice: Categorizing smoking history to aid in interpretation
# and account for potential non-linear effects
cox_data <- cox_data %>%
  mutate(
    SMOKING_CAT = case_when(
      SMOKING_PACK_YEARS == 0 ~ "Non-smoker",
      SMOKING_PACK_YEARS > 0 & SMOKING_PACK_YEARS <= 20 ~ "Light",
      SMOKING_PACK_YEARS > 20 & SMOKING_PACK_YEARS <= 40 ~ "Moderate",
      SMOKING_PACK_YEARS > 40 ~ "Heavy",
      TRUE ~ NA_character_
    ),
    # Ensure proper factor ordering
    SMOKING_CAT = factor(
      SMOKING_CAT, 
      levels = c("Non-smoker", "Light", "Moderate", "Heavy")
    )
  )

# Create formula for Cox model
# Design choice: Including all key genes and clinical factors in the model
# to identify independent prognostic factors
cox_formula <- as.formula(paste0(
  "Surv(time, status) ~ ", 
  paste(key_genes, collapse = " + "), 
  " + STAGE_SIMPLE + SEX + SMOKING_CAT"
))

# Remove rows with missing data for the Cox model
# Design choice: Using complete cases to ensure valid model fitting,
# though this could introduce selection bias if missingness is informative
cox_data_complete <- cox_data %>%
  drop_na(time, status, all_of(c(key_genes, "STAGE_SIMPLE", "SEX", "SMOKING_CAT")))

cat("Cox model will use", nrow(cox_data_complete), "out of", nrow(cox_data), "patients with complete data\n")

# Run the Cox model
cat("Fitting Cox proportional hazards model...\n")
cox_result <- coxph(cox_formula, data = cox_data_complete)
cox_summary <- summary(cox_result)

# Save the detailed results
capture.output(cox_summary, file = "results/tables/cox_regression_results.txt")

# Extract HR and CI for forest plot
# Design choice: Creating a clean data frame with model results for visualization
cox_df <- tibble(
  variable = rownames(cox_summary$conf.int),
  HR = cox_summary$conf.int[, "exp(coef)"],
  lower = cox_summary$conf.int[, "lower .95"],
  upper = cox_summary$conf.int[, "upper .95"],
  p_value = cox_summary$coefficients[, "Pr(>|z|)"]
) %>%
  # Add significance markers for easy interpretation
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Create forest plot
# Design choice: Using a forest plot to visualize hazard ratios and confidence intervals,
# which is the standard visualization for Cox regression results
p <- ggplot(cox_df, aes(x = HR, y = variable)) +
  # Add points for hazard ratios
  geom_point(size = 3) +
  # Add error bars for confidence intervals
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  # Add reference line at HR = 1 (no effect)
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  # Use log scale for better visualization of ratios
  scale_x_log10(
    breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8),
    labels = c("0.1", "0.25", "0.5", "1", "2", "4", "8")
  ) +
  # Add informative labels
  labs(
    x = "Hazard Ratio (95% CI)",
    y = "",
    title = "Cox Proportional Hazards Model",
    subtitle = "Multivariate Analysis of Survival Factors in Lung Adenocarcinoma"
  ) +
  # Use a clean theme
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title.x = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  # Add significance markers
  annotate(
    "text", 
    x = max(cox_df$upper) * 1.1, 
    y = 1:nrow(cox_df), 
    label = cox_df$significance, 
    size = 5
  )

# Save the forest plot
ggsave(
  "results/figures/integrated/cox_model_forest_plot.png", 
  plot = p, 
  width = 10, height = 8, dpi = 300
)
# Create a nice table for the Cox results
cox_table <- data.frame(
  Variable = cox_df$variable,
  HR = round(cox_df$HR, 2),
  CI_95 = paste0(round(cox_df$lower, 2), "-", round(cox_df$upper, 2)),
  P_value = round(cox_df$p_value, 4),
  Significance = cox_df$significance
)
write.csv(cox_table, "results/tables/cox_regression_table.csv", row.names = FALSE)

cat("Integrated analysis complete. Results saved to 'results/figures/integrated/' and 'results/tables/'\n")

