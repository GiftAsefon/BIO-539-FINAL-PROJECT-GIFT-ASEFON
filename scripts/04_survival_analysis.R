###############################################################
# SCRIPT: 04_survival_analysis.R
# PURPOSE: Survival analysis by clinical factors
# AUTHOR: GIFT ASEFON
# DATE: 04-28-2024
###############################################################

# Load libraries
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

# Load preprocessed data
load("data/preprocessed_data.RData")

# Create survival data from clinical patient data
surv_data <- data.frame(
  PATIENT_ID = clinical_patient$PATIENT_ID,
  time = clinical_patient$OS_MONTHS,
  status = ifelse(clinical_patient$OS_STATUS == "1:DECEASED", 1, 0)
)

# Remove missing values
surv_data <- surv_data[!is.na(surv_data$time) & surv_data$time > 0, ]

# Save survival data for other analyses
save(surv_data, file = "data/survival_data.RData")

# --- Survival by Stage ---
# Create stage data
path_stage_data <- data.frame(
  PATIENT_ID = clinical_patient$PATIENT_ID,
  PATH_STAGE = clinical_patient$PATH_STAGE
)

# Merge with survival data
surv_stage <- merge(surv_data, path_stage_data, by = "PATIENT_ID")
surv_stage <- surv_stage[!is.na(surv_stage$PATH_STAGE) & surv_stage$PATH_STAGE != "", ]

# Create survival object
fit_stage <- survfit(Surv(time, status) ~ PATH_STAGE, data = surv_stage)

# Base R plot
png("results/figures/survival/survival_by_stage_base.png", width = 10, height = 8, units = "in", res = 300)
plot(fit_stage, 
     col = c("blue", "green", "orange", "red"),
     xlab = "Time (months)",
     ylab = "Survival Probability",
     main = "Survival by Pathologic Stage",
     lwd = 2)
legend("topright", 
       legend = levels(factor(surv_stage$PATH_STAGE)),
       col = c("blue", "green", "orange", "red"),
       lty = 1,
       lwd = 2)
dev.off()

# Log-rank test for stage differences
log_rank_stage <- survdiff(Surv(time, status) ~ PATH_STAGE, data = surv_stage)
stage_p_value <- 1 - pchisq(log_rank_stage$chisq, df = length(levels(factor(surv_stage$PATH_STAGE))) - 1)
cat("Log-rank test for stage differences: p-value =", stage_p_value, "\n")

# Try ggsurvplot version
tryCatch({
  stage_plot <- ggsurvplot(fit_stage, 
                           data = surv_stage,
                           risk.table = TRUE,
                           pval = TRUE,
                           title = "Survival by Pathologic Stage",
                           palette = "jco",
                           xlab = "Time (months)",
                           ggtheme = theme_minimal(),
                           risk.table.height = 0.25)
  
  # Save the plot
  png("results/figures/survival/survival_by_stage.png", width = 10, height = 8, units = "in", res = 300)
  print(stage_plot)
  dev.off()
}, error = function(e) {
  cat("Error with ggsurvplot for stage:", e$message, "\n")
})

# --- Survival by Gender ---
# Create gender data
gender_data <- data.frame(
  PATIENT_ID = clinical_patient$PATIENT_ID,
  SEX = clinical_patient$SEX
)

# Merge with survival data
surv_gender <- merge(surv_data, gender_data, by = "PATIENT_ID")
surv_gender <- surv_gender[!is.na(surv_gender$SEX), ]

# Create survival object
fit_gender <- survfit(Surv(time, status) ~ SEX, data = surv_gender)

# Create survival plot using base R
png("results/figures/survival/survival_by_gender_base.png", width = 10, height = 8, units = "in", res = 300)
plot(fit_gender, 
     col = c("#0072B2", "#D55E00"),
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

# Log-rank test for gender differences
log_rank_gender <- survdiff(Surv(time, status) ~ SEX, data = surv_gender)
gender_p_value <- 1 - pchisq(log_rank_gender$chisq, df = 1)
cat("Log-rank test for gender differences: p-value =", gender_p_value, "\n")

# Try ggsurvplot version
tryCatch({
  gender_plot <- ggsurvplot(fit_gender, 
                            data = surv_gender,
                            risk.table = TRUE,
                            pval = TRUE,
                            title = "Survival by Gender",
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

# --- Survival by Smoking Status ---
# Create smoking data with categories
smoking_data <- data.frame(
  PATIENT_ID = clinical_patient$PATIENT_ID,
  SMOKING_PACK_YEARS = clinical_patient$SMOKING_PACK_YEARS
)

# Create categories for smoking pack years
smoking_data$SMOKING_GROUP <- cut(smoking_data$SMOKING_PACK_YEARS, 
                                  breaks = c(-1, 0, 20, 40, Inf),
                                  labels = c("Non-smoker", "Light", "Moderate", "Heavy"))

# Merge with survival data
surv_smoking <- merge(surv_data, smoking_data, by = "PATIENT_ID")
surv_smoking <- surv_smoking[!is.na(surv_smoking$SMOKING_GROUP), ]

# Create survival object
fit_smoking <- survfit(Surv(time, status) ~ SMOKING_GROUP, data = surv_smoking)

# Create survival plot using base R
png("results/figures/survival/survival_by_smoking_base.png", width = 10, height = 8, units = "in", res = 300)
plot(fit_smoking, 
     col = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"),
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

# Log-rank test for smoking differences
log_rank_smoking <- survdiff(Surv(time, status) ~ SMOKING_GROUP, data = surv_smoking)
smoking_p_value <- 1 - pchisq(log_rank_smoking$chisq, df = length(levels(surv_smoking$SMOKING_GROUP)) - 1)
cat("Log-rank test for smoking differences: p-value =", smoking_p_value, "\n")

# Try ggsurvplot version
tryCatch({
  smoking_plot <- ggsurvplot(fit_smoking, 
                             data = surv_smoking,
                             risk.table = TRUE,
                             pval = TRUE,
                             title = "Survival by Smoking History",
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

# Create and save table with survival p-values
survival_pvalues <- data.frame(
  Factor = c("Pathologic Stage", "Gender", "Smoking Status"),
  P_value = c(stage_p_value, gender_p_value, smoking_p_value)
)
write.csv(survival_pvalues, "results/tables/survival_pvalues.csv", row.names = FALSE)

cat("Survival analysis by clinical factors complete. Results saved to 'results/figures/survival/' and 'results/tables/'\n")