###############################################################
# SCRIPT: 04_survival_analysis.R
# PURPOSE: Survival analysis by clinical factors
# AUTHOR: GIFT ASEFON
# DATE: 04-28-2025
###############################################################

###############################################################
# OVERVIEW:
# This script performs survival analysis on lung adenocarcinoma patients using
# the TCGA dataset. The primary goals are:
# 1. Evaluate the impact of pathologic stage on patient survival
# 2. Examine survival differences between genders
# 3. Assess the relationship between smoking history and survival outcomes
# 4. Generate publication-quality survival curves
# 5. Perform statistical tests (log-rank) to assess significance of findings
###############################################################

# Load libraries
library(survival)    # For survival analysis functions (survfit, Surv, survdiff)
library(survminer)   # For enhanced survival curve visualization with ggsurvplot
library(tidyverse)   # For data manipulation and visualization

# Load preprocessed data
# This data was prepared in 01_data_preprocessing.R and contains
# clinical information for lung adenocarcinoma patients
load("data/preprocessed_data.RData")

###############################################################
# STEP 1: PREPARE SURVIVAL DATA
###############################################################

# Create survival data from clinical patient data
# Design choice: Using tidyverse for data preparation to improve readability
# and make the code more maintainable
surv_data <- clinical_patient %>%
  # Select only the columns needed for survival analysis
  select(PATIENT_ID, OS_MONTHS, OS_STATUS) %>%
  # Create a proper status indicator (1=event/death, 0=censored/alive)
  # This is the format expected by the survival package
  mutate(status = ifelse(OS_STATUS == "1:DECEASED", 1, 0)) %>%
  # Rename for clarity - 'time' is a more standard name for follow-up time
  rename(time = OS_MONTHS) %>%
  # Remove missing values or zero/negative survival times
  # This is critical for valid survival analysis
  filter(!is.na(time), time > 0)

# Save survival data for other analyses
# This ensures consistent survival data across all analyses
save(surv_data, file = "data/survival_data.RData")

###############################################################
# STEP 2: SURVIVAL ANALYSIS BY PATHOLOGIC STAGE
###############################################################

# Create stage data and merge with survival data
# Design choice: Using tidyverse join operations instead of base R merge
# This improves code readability and maintains consistency with the tidyverse approach
surv_stage <- clinical_patient %>%
  # Select only the columns needed for this analysis
  select(PATIENT_ID, PATH_STAGE) %>%
  # Remove patients with missing or empty stage information
  filter(!is.na(PATH_STAGE), PATH_STAGE != "") %>%
  # Join with survival data
  inner_join(surv_data, by = "PATIENT_ID")

# Create survival object
# The survfit function creates a Kaplan-Meier survival curve
# stratified by pathologic stage
fit_stage <- survfit(Surv(time, status) ~ PATH_STAGE, data = surv_stage)

# Create survival plot using base R
# Design choice: Creating a base R plot first ensures we have a reliable
# visualization even if the ggsurvplot version fails
png("results/figures/survival/survival_by_stage_base.png", width = 10, height = 8, units = "in", res = 300)
plot(fit_stage, 
     # Use different colors for each stage to improve readability
     col = c("blue", "green", "orange", "red"),
     xlab = "Time (months)",
     ylab = "Survival Probability",
     main = "Survival by Pathologic Stage",
     # Increase line width for better visibility
     lwd = 2)
# Add a legend to identify each stage
legend("topright", 
       legend = levels(factor(surv_stage$PATH_STAGE)),
       col = c("blue", "green", "orange", "red"),
       lty = 1,
       lwd = 2)
dev.off()

# Perform log-rank test for stage differences
# Design choice: The log-rank test is the standard method for comparing
# survival curves in different groups
log_rank_stage <- survdiff(Surv(time, status) ~ PATH_STAGE, data = surv_stage)
# Calculate p-value from chi-square distribution
# The degrees of freedom is the number of groups minus 1
stage_p_value <- 1 - pchisq(log_rank_stage$chisq, 
                            df = n_distinct(surv_stage$PATH_STAGE) - 1)
cat("Log-rank test for stage differences: p-value =", stage_p_value, "\n")

# Try creating a more advanced plot with ggsurvplot
# Design choice: ggsurvplot creates publication-quality plots with risk tables
# and p-values incorporated, but may fail in some environments
tryCatch({
  # Create ggsurvplot with enhanced features
  stage_plot <- ggsurvplot(fit_stage, 
                           data = surv_stage,
                           # Add a risk table showing number of patients at risk over time
                           risk.table = TRUE,
                           # Add p-value from log-rank test
                           pval = TRUE,
                           title = "Survival by Pathologic Stage",
                           # Use journal-quality color palette
                           palette = "jco",
                           xlab = "Time (months)",
                           # Use minimal theme for cleaner appearance
                           ggtheme = theme_minimal(),
                           # Adjust risk table height for better proportions
                           risk.table.height = 0.25)
  
  # Save the plot
  png("results/figures/survival/survival_by_stage.png", width = 10, height = 8, units = "in", res = 300)
  print(stage_plot)
  dev.off()
}, error = function(e) {
  # Catch and report any errors that occur
  # This ensures the script continues running even if ggsurvplot fails
  cat("Error with ggsurvplot for stage:", e$message, "\n")
})

###############################################################
# STEP 3: SURVIVAL ANALYSIS BY GENDER
###############################################################

# Create gender data and merge with survival data
# Design choice: Using tidyverse approach for consistency
surv_gender <- clinical_patient %>%
  # Select only the relevant columns
  select(PATIENT_ID, SEX) %>%
  # Remove patients with missing gender information
  filter(!is.na(SEX)) %>%
  # Join with survival data
  inner_join(surv_data, by = "PATIENT_ID")

# Create survival object stratified by gender
fit_gender <- survfit(Surv(time, status) ~ SEX, data = surv_gender)

# Create survival plot using base R
# Design choice: Using distinctive, color-blind friendly colors for gender
png("results/figures/survival/survival_by_gender_base.png", width = 10, height = 8, units = "in", res = 300)
plot(fit_gender, 
     col = c("#0072B2", "#D55E00"),  # Blue and orange - color-blind friendly
     xlab = "Time (months)",
     ylab = "Survival Probability",
     main = "Survival by Gender",
     lwd = 2)
legend("topright", 
       legend = levels(factor(surv_gender$SEX)),
       col = c("#0072B2", "#D55E00"),
       lty = 1,
       lwd = 2)
dev.off()

# Perform log-rank test for gender differences
# Design choice: For binary variables like gender, df=1
log_rank_gender <- survdiff(Surv(time, status) ~ SEX, data = surv_gender)
gender_p_value <- 1 - pchisq(log_rank_gender$chisq, df = 1)
cat("Log-rank test for gender differences: p-value =", gender_p_value, "\n")

# Try ggsurvplot version with the same error handling approach
tryCatch({
  gender_plot <- ggsurvplot(fit_gender, 
                            data = surv_gender,
                            risk.table = TRUE,
                            pval = TRUE,
                            title = "Survival by Gender",
                            # Use the same color scheme as the base R plot
                            palette = c("#0072B2", "#D55E00"),
                            xlab = "Time (months)",
                            ggtheme = theme_minimal())
  
  # Save the plot
  png("results/figures/survival/survival_by_gender.png", width = 10, height = 8, units = "in", res = 300)
  print(gender_plot)
  dev.off()
}, error = function(e) {
  cat("Error with ggsurvplot for gender:", e$message, "\n")
})

###############################################################
# STEP 4: SURVIVAL ANALYSIS BY SMOKING STATUS
###############################################################

# Create smoking data with categories and merge with survival data
# Design choice: Using tidyverse for data transformation and joining
surv_smoking <- clinical_patient %>%
  # Select only the relevant columns
  select(PATIENT_ID, SMOKING_PACK_YEARS) %>%
  # Create categories for smoking pack years
  # Design choice: Using case_when instead of cut for more readable code
  mutate(SMOKING_GROUP = case_when(
    is.na(SMOKING_PACK_YEARS) ~ NA_character_,
    SMOKING_PACK_YEARS == 0 ~ "Non-smoker",
    SMOKING_PACK_YEARS > 0 & SMOKING_PACK_YEARS <= 20 ~ "Light",
    SMOKING_PACK_YEARS > 20 & SMOKING_PACK_YEARS <= 40 ~ "Moderate",
    SMOKING_PACK_YEARS > 40 ~ "Heavy"
  )) %>%
  # Convert to factor with proper ordering
  mutate(SMOKING_GROUP = factor(SMOKING_GROUP, 
                                levels = c("Non-smoker", "Light", "Moderate", "Heavy"))) %>%
  # Remove patients with missing smoking information
  filter(!is.na(SMOKING_GROUP)) %>%
  # Join with survival data
  inner_join(surv_data, by = "PATIENT_ID")

# Create survival object stratified by smoking group
fit_smoking <- survfit(Surv(time, status) ~ SMOKING_GROUP, data = surv_smoking)

# Create survival plot using base R
# Design choice: Using a different color palette to distinguish from other plots
png("results/figures/survival/survival_by_smoking_base.png", width = 10, height = 8, units = "in", res = 300)
plot(fit_smoking, 
     col = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),  # ColorBrewer palette
     xlab = "Time (months)",
     ylab = "Survival Probability",
     main = "Survival by Smoking History",
     lwd = 2)
legend("topright", 
       legend = levels(surv_smoking$SMOKING_GROUP),
       col = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),
       lty = 1,
       lwd = 2)
dev.off()

# Perform log-rank test for smoking differences
# Design choice: The degrees of freedom is the number of smoking groups minus 1
log_rank_smoking <- survdiff(Surv(time, status) ~ SMOKING_GROUP, data = surv_smoking)
smoking_p_value <- 1 - pchisq(log_rank_smoking$chisq, 
                              df = n_distinct(surv_smoking$SMOKING_GROUP) - 1)
cat("Log-rank test for smoking differences: p-value =", smoking_p_value, "\n")

# Try ggsurvplot version
tryCatch({
  smoking_plot <- ggsurvplot(fit_smoking, 
                             data = surv_smoking,
                             risk.table = TRUE,
                             pval = TRUE,
                             title = "Survival by Smoking History",
                             # Using a predefined palette for consistency
                             palette = "jco",
                             xlab = "Time (months)",
                             ggtheme = theme_minimal())
  
  # Save the plot
  png("results/figures/survival/survival_by_smoking.png", width = 10, height = 8, units = "in", res = 300)
  print(smoking_plot)
  dev.off()
}, error = function(e) {
  cat("Error with ggsurvplot for smoking:", e$message, "\n")
})

###############################################################
# STEP 5: SUMMARIZE SURVIVAL ANALYSIS RESULTS
###############################################################

# Create and save table with survival p-values
# Design choice: Using tidyverse to create and save the summary table
# This maintains consistency with the rest of the code
tibble(
  Factor = c("Pathologic Stage", "Gender", "Smoking Status"),
  P_value = c(stage_p_value, gender_p_value, smoking_p_value)
) %>%
  # Add a column indicating statistical significance
  mutate(Significant = ifelse(P_value < 0.05, "Yes", "No")) %>%
  # Add formatted p-values for better readability in reports
  mutate(P_value_formatted = case_when(
    P_value < 0.001 ~ "<0.001",
    P_value < 0.01 ~ sprintf("%.3f", P_value),
    TRUE ~ sprintf("%.2f", P_value)
  )) %>%
  # Save to CSV for easy import into other software
  write_csv("results/tables/survival_pvalues.csv")

# Print completion message with information about saved results
cat("Survival analysis by clinical factors complete. Results saved to 'results/figures/survival/' and 'results/tables/'\n")
