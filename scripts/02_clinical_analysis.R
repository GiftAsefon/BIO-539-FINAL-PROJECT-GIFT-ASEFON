###############################################################
# SCRIPT: 02_clinical_analysis.R
# PURPOSE: Demographic and clinical characteristics analysis
# AUTHOR: GIFT ASEFON
# DATE: 04-28-2025
###############################################################

# Load libraries
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(dplyr)

# Load preprocessed data
load("data/preprocessed_data.RData")

# Basic summary of patient demographics
demographic_summary <- clinical_patient %>%
  summarize(
    total_patients = n(),
    median_age = median(AGE, na.rm = TRUE),
    male_count = sum(SEX == "Male", na.rm = TRUE),
    female_count = sum(SEX == "Female", na.rm = TRUE),
    smoker_years_median = median(SMOKER_YEARS, na.rm = TRUE),
    smoking_pack_years_median = median(SMOKING_PACK_YEARS, na.rm = TRUE),
    alive_count = sum(VITAL_STATUS == "Alive", na.rm = TRUE),
    deceased_count = sum(VITAL_STATUS == "Dead", na.rm = TRUE)
  )

print(demographic_summary)
write.csv(demographic_summary, "results/tables/demographic_summary.csv", row.names = FALSE)

# Find stage-related columns
stage_columns <- grep("stage|tumor|grade", colnames(clinical_patient), 
                      value = TRUE, ignore.case = TRUE)
print(paste("Stage-related columns:", paste(stage_columns, collapse=", ")))

# Create a list to store plots
plot_list <- list()

# 1. Age Distribution plot
age_data <- clinical_patient$AGE
age_data <- age_data[!is.na(age_data)]

age_plot <- ggplot(data.frame(AGE = age_data), aes(x = AGE)) +
  geom_histogram(binwidth = 5, fill = "#56B4E9", color = "white", alpha = 0.8) +
  geom_density(aes(y = after_stat(count) * 5), color = "#0072B2", linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Age Distribution of Lung Adenocarcinoma Patients",
    x = "Age at Diagnosis (Years)",
    y = "Number of Patients"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

plot_list[["age"]] <- age_plot

# 2. Gender Distribution plot
sex_table <- table(clinical_patient$SEX)
gender_data <- data.frame(
  SEX = names(sex_table),
  n = as.numeric(sex_table),
  stringsAsFactors = FALSE
) %>%
  mutate(percentage = n / sum(n) * 100)

gender_plot <- ggplot(gender_data, aes(x = "", y = n, fill = SEX)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = paste0(round(percentage), "%")), 
            position = position_stack(vjust = 0.5), 
            size = 4, color = "white", fontface = "bold") +
  coord_polar("y", start = 0) +
  theme_minimal() +
  labs(
    title = "Gender Distribution",
    fill = "Gender"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = c("Male" = "#0072B2", "Female" = "#D55E00"))

plot_list[["gender"]] <- gender_plot

# 3. Smoking Pack Years Distribution
smoking_data <- clinical_patient$SMOKING_PACK_YEARS
smoking_data <- smoking_data[!is.na(smoking_data)]

smoking_plot <- ggplot(data.frame(SMOKING_PACK_YEARS = smoking_data), aes(x = SMOKING_PACK_YEARS)) +
  geom_histogram(binwidth = 10, fill = "#56B4E9", color = "white", alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Distribution of Smoking History (Pack Years)",
    x = "Smoking Pack Years",
    y = "Number of Patients"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

plot_list[["smoking"]] <- smoking_plot

# 4. Vital Status Distribution
vital_table <- table(clinical_patient$VITAL_STATUS)
vital_data <- data.frame(
  VITAL_STATUS = names(vital_table),
  n = as.numeric(vital_table),
  stringsAsFactors = FALSE
)
vital_data$percentage <- vital_data$n / sum(vital_data$n) * 100

vital_plot <- ggplot(vital_data, aes(x = "", y = n, fill = VITAL_STATUS)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = paste0(round(percentage), "%")), 
            position = position_stack(vjust = 0.5), 
            size = 4, color = "white", fontface = "bold") +
  coord_polar("y", start = 0) +
  theme_minimal() +
  labs(
    title = "Vital Status Distribution",
    fill = "Status"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values = c("Alive" = "#009E73", "Dead" = "#CC79A7"))

plot_list[["vital"]] <- vital_plot

# 5. Race Distribution
race_table <- table(clinical_patient$RACE)
race_data <- data.frame(
  RACE = names(race_table),
  n = as.numeric(race_table),
  stringsAsFactors = FALSE
)
# Remove empty race entries (if any)
race_data <- race_data[race_data$RACE != "", ]
# Sort by frequency
race_data <- race_data[order(race_data$n, decreasing = TRUE), ]

race_plot <- ggplot(race_data, aes(x = reorder(RACE, -n), y = n, fill = n)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  scale_fill_gradient(low = "#56B4E9", high = "#0072B2") +
  theme_minimal() +
  labs(
    title = "Distribution by Race",
    x = "Race",
    y = "Number of Patients"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

plot_list[["race"]] <- race_plot

# 6. Ethnicity Distribution
ethnicity_table <- table(clinical_patient$ETHNICITY)
ethnicity_data <- data.frame(
  ETHNICITY = names(ethnicity_table),
  n = as.numeric(ethnicity_table),
  stringsAsFactors = FALSE
)
# Remove empty ethnicity entries (if any)
ethnicity_data <- ethnicity_data[ethnicity_data$ETHNICITY != "", ]
# Sort by frequency
ethnicity_data <- ethnicity_data[order(ethnicity_data$n, decreasing = TRUE), ]

ethnicity_plot <- ggplot(ethnicity_data, aes(x = reorder(ETHNICITY, -n), y = n, fill = n)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = n), vjust = -0.3, size = 3.5) +
  scale_fill_gradient(low = "#56B4E9", high = "#0072B2") +
  theme_minimal() +
  labs(
    title = "Distribution by Ethnicity",
    x = "Ethnicity",
    y = "Number of Patients"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

plot_list[["ethnicity"]] <- ethnicity_plot

# Arrange all plots in a grid and save combined plot
combined_plot <- arrangeGrob(
  age_plot, gender_plot,
  smoking_plot, vital_plot,
  race_plot, ethnicity_plot,
  ncol = 2
)

# Save individual plots
ggsave("results/figures/clinical/age_distribution.png", plot = age_plot, width = 8, height = 6, dpi = 300)
ggsave("results/figures/clinical/gender_distribution.png", plot = gender_plot, width = 8, height = 6, dpi = 300)
ggsave("results/figures/clinical/smoking_distribution.png", plot = smoking_plot, width = 8, height = 6, dpi = 300)
ggsave("results/figures/clinical/vital_status_distribution.png", plot = vital_plot, width = 8, height = 6, dpi = 300)
ggsave("results/figures/clinical/race_distribution.png", plot = race_plot, width = 8, height = 6, dpi = 300)
ggsave("results/figures/clinical/ethnicity_distribution.png", plot = ethnicity_plot, width = 8, height = 6, dpi = 300)

# Save combined plot
ggsave(
  filename = "results/figures/clinical/demographics_combined.png",
  plot = combined_plot,
  width = 12,
  height = 15,
  dpi = 300
)

cat("Clinical analysis complete. Results saved to 'results/figures/clinical/' and 'results/tables/'\n")
