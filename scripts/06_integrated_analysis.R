###############################################################
# SCRIPT: 06_integrated_analysis.R
# PURPOSE: Integration of mutation, expression, and survival data
# AUTHOR: GIFT ASEFON
# DATE: 04-28-2025
###############################################################

# Load libraries
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

# Load preprocessed data
load("data/preprocessed_data.RData")
load("data/survival_data.RData")

# Function to convert sample ID to patient ID
sample_to_patient <- function(sample_id) {
  return(gsub("-01.*$", "", sample_id))
}

# Create mutation status for key genes
mutation_status <- data.frame(PATIENT_ID = unique(clinical_patient$PATIENT_ID))

for (gene in key_genes) {
  # Get patients with mutations in this gene
  mutated_samples <- mutations$Tumor_Sample_Barcode[mutations$Hugo_Symbol == gene]
  mutated_patients <- unique(sapply(mutated_samples, sample_to_patient))
  
  # Add mutation status to data frame (1 = mutated, 0 = wild-type)
  mutation_status[[gene]] <- ifelse(mutation_status$PATIENT_ID %in% mutated_patients, 1, 0)
}

# Merge with survival data
surv_mutations <- merge(surv_data, mutation_status, by = "PATIENT_ID")

# Create Kaplan-Meier plots for each key gene
for (gene in key_genes) {
  # Check if we have enough mutations to analyze
  mutations_count <- sum(surv_mutations[[gene]], na.rm = TRUE)
  cat(paste0("Found ", mutations_count, " mutations for ", gene, "\n"))
  
  if (mutations_count >= 5) {
    # Create factor for visualization
    surv_mutations[[paste0(gene, "_status")]] <- factor(
      surv_mutations[[gene]],
      levels = c(0, 1),
      labels = c("Wild-type", "Mutant")
    )
    
    # Create survival formula and fit
    formula <- as.formula(paste0("Surv(time, status) ~ ", gene, "_status"))
    
    tryCatch({
      fit_gene <- survfit(formula, data = surv_mutations)
      
      # Create plot using base R
      png(paste0("results/figures/integrated/survival_by_", gene, "_mutation_base.png"), 
          width = 10, height = 8, units = "in", res = 300)
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
      
      # Add p-value from log-rank test
      log_rank <- survdiff(formula, data = surv_mutations)
      p_value <- 1 - pchisq(log_rank$chisq, df = 1)
      text(x = max(fit_gene$time) * 0.7, y = 0.2, 
           labels = paste("p =", round(p_value, 4)),
           cex = 1.2)
      
      dev.off()
      
      # Try with ggsurvplot if possible
      tryCatch({
        gene_plot <- ggsurvplot(fit_gene,
                                data = surv_mutations,
                                risk.table = TRUE,
                                pval = TRUE,
                                title = paste0("Survival by ", gene, " Mutation Status"),
                                palette = c("#2166AC", "#B2182B"),
                                xlab = "Time (months)",
                                ggtheme = theme_minimal())
        
        # Save the ggsurvplot
        png(paste0("results/figures/integrated/survival_by_", gene, "_mutation.png"), 
            width = 10, height = 8, units = "in", res = 300)
        print(gene_plot)
        dev.off()
        
        cat(paste0("Successfully created ggsurvplot for ", gene, "\n"))
      }, error = function(e) {
        cat(paste0("Error with ggsurvplot for ", gene, ": ", e$message, "\n"))
        cat("Base R plot has been created as an alternative\n")
      })
      
      # Save p-value to a text file for reference
      cat(paste0("Log-rank test for ", gene, ": p-value = ", round(p_value, 4), "\n"), 
          file = "results/tables/mutation_survival_pvalues.txt", 
          append = TRUE, sep = "")
      
      cat(paste0("Successfully created survival plots for ", gene, " - p-value = ", round(p_value, 4), "\n"))
    }, error = function(e) {
      cat(paste0("Error creating survival analysis for ", gene, ": ", e$message, "\n"))
    })
  } else {
    cat(paste0("Skipping ", gene, " - not enough mutations (", mutations_count, ")\n"))
  }
}

# Create more comprehensive dataset for Cox analysis
cox_data <- surv_mutations

# Add stage information
stage_data <- data.frame(
  PATIENT_ID = clinical_patient$PATIENT_ID,
  PATH_STAGE = clinical_patient$PATH_STAGE
)
cox_data <- merge(cox_data, stage_data, by = "PATIENT_ID", all.x = TRUE)

# Add gender information
gender_data <- data.frame(
  PATIENT_ID = clinical_patient$PATIENT_ID,
  SEX = clinical_patient$SEX
)
cox_data <- merge(cox_data, gender_data, by = "PATIENT_ID", all.x = TRUE)

# Add smoking information
smoking_data <- data.frame(
  PATIENT_ID = clinical_patient$PATIENT_ID,
  SMOKING_PACK_YEARS = clinical_patient$SMOKING_PACK_YEARS
)
cox_data <- merge(cox_data, smoking_data, by = "PATIENT_ID", all.x = TRUE)

# Simplify stage for analysis
cox_data$STAGE_SIMPLE <- NA
cox_data$STAGE_SIMPLE[grep("I", cox_data$PATH_STAGE)] <- "Stage I"
cox_data$STAGE_SIMPLE[grep("II", cox_data$PATH_STAGE)] <- "Stage II"
cox_data$STAGE_SIMPLE[grep("III", cox_data$PATH_STAGE)] <- "Stage III"
cox_data$STAGE_SIMPLE[grep("IV", cox_data$PATH_STAGE)] <- "Stage IV"
cox_data$STAGE_SIMPLE <- factor(cox_data$STAGE_SIMPLE, 
                                levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))

# Create a categorical variable for smoking
cox_data$SMOKING_CAT <- cut(cox_data$SMOKING_PACK_YEARS,
                            breaks = c(-Inf, 0, 20, 40, Inf),
                            labels = c("Non-smoker", "Light", "Moderate", "Heavy"))

# Create formula for Cox model
cox_formula <- as.formula(paste0("Surv(time, status) ~ ", 
                                 paste(key_genes, collapse = " + "), 
                                 " + STAGE_SIMPLE + SEX + SMOKING_CAT"))

# Remove rows with missing data for the Cox model
cox_data_complete <- cox_data[complete.cases(cox_data[, 
                                                      c("time", "status", key_genes, 
                                                        "STAGE_SIMPLE", "SEX", "SMOKING_CAT")]), ]

# Run the Cox model
cox_result <- coxph(cox_formula, data = cox_data_complete)
cox_summary <- summary(cox_result)

# Save the results
capture.output(cox_summary, file = "results/tables/cox_regression_results.txt")

# Extract HR and CI for forest plot
cox_df <- data.frame(
  variable = rownames(cox_summary$conf.int),
  HR = cox_summary$conf.int[, "exp(coef)"],
  lower = cox_summary$conf.int[, "lower .95"],
  upper = cox_summary$conf.int[, "upper .95"],
  p_value = cox_summary$coefficients[, "Pr(>|z|)"]
)

# Add significance markers
cox_df$significance <- ifelse(cox_df$p_value < 0.001, "***",
                              ifelse(cox_df$p_value < 0.01, "**",
                                     ifelse(cox_df$p_value < 0.05, "*", "ns")))

# Create forest plot
p <- ggplot(cox_df, aes(x = HR, y = variable)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  scale_x_log10() +
  labs(x = "Hazard Ratio (95% CI)", y = "",
       title = "Cox Proportional Hazards Model") +
  theme_minimal() +
  annotate("text", x = max(cox_df$upper) * 1.1, y = 1:nrow(cox_df), 
           label = cox_df$significance, size = 5)

# Save the forest plot
png("results/figures/integrated/cox_model_forest_plot.png", width = 10, height = 8, units = "in", res = 300)
print(p)
dev.off()

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

