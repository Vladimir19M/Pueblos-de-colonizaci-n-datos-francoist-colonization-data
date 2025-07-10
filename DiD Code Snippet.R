
# TFG Code Matching 
# Vladimir Miculaiciuc Kudriavtsev
# The code used for the actual results is not identical, However the functional form remains the same.

#------------------------------------------
# 1. Libraries
#------------------------------------------
library(readr)
library(dplyr)
library(fixest)
library(ggplot2)
library(broom)
library(modelsummary)
library(car)        
library(kableExtra)
library(did)    

#------------------------------------------
# 2. Load Data
#------------------------------------------
data <- read_delim(
  "data/Match Pairs1.csv",
  delim = ";",
  locale = locale(decimal_mark = ",", grouping_mark = "."),
  show_col_types = FALSE
)

#------------------------------------------
# 3. Event Study Function
#------------------------------------------
run_event_study_with_controls <- function(data, cohort_var, treatment_year, ref_period = -10) {
  data %>%
    mutate(
      treat = !!sym(cohort_var),
      time_to_treat = year - treatment_year,
      log_pop = log(population + 1)
    ) %>%
    filter(treat == 1 | is_treated == 0) %>%
    feols(
      log_pop ~ i(time_to_treat, treat, ref = ref_period) +
        log_area + sqrt_urban_dist + river_distance_m + Latitude +
        sqrt_urban_dist:year,
      weights = ~ weight,
      cluster = ~ municipality_ID
    )
}

#------------------------------------------
# 4. Define Cohorts
#------------------------------------------
cohorts <- tibble::tibble(
  var = c("treatment_age_50", "treatment_age_60", "treatment_age_70"),
  year = c(1940, 1950, 1960),
  name = c("1940s", "1950s", "1960s")
)

#------------------------------------------
# 5. Estimate Models
#------------------------------------------
event_models <- cohorts %>%
  rowwise() %>%
  mutate(
    model = list(run_event_study_with_controls(data, var, year)),
    tidy = list(tidy(model, conf.int = TRUE) %>% mutate(cohort = name))
  )

#------------------------------------------
# 6. Clean Results
#------------------------------------------
combined_event <- bind_rows(event_models$tidy) %>%
  filter(grepl("time_to_treat", term)) %>%
  mutate(
    time = as.numeric(sub("time_to_treat::(-?\\d+):treat", "\\1", term)),
    significant = p.value < 0.05
  ) %>%
  filter(!is.na(time), time >= -20, time <= 30)

#------------------------------------------
# 7. Pre-Trend Hypothesis Tests
#------------------------------------------
pretrend_tests <- lapply(event_models$model, function(m) {
  tryCatch(wald(m, keep = "time_to_treat::-1|time_to_treat::-2"), error = function(e) NULL)
})

#------------------------------------------
# 8. Plot Event Study
#------------------------------------------
if (nrow(combined_event) > 0) {
  ggplot(combined_event, aes(x = time, y = estimate, color = cohort)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, color = "black") +
    geom_line() +
    geom_point(aes(shape = significant)) +
    scale_x_continuous(breaks = seq(-20, 30, 10)) +
    labs(
      title = "Event Study with Controls (Log Population)",
      subtitle = "Years Relative to Treatment (Ref = -10)",
      x = "Years Relative to Treatment",
      y = "Estimated Effect"
    ) +
    theme_minimal() +
    facet_wrap(~cohort, ncol = 1)
} else {
  warning("No event study results to plot")
}

#------------------------------------------
# 9. Placebo Test (Fake 1930)
#------------------------------------------
placebo_model <- run_event_study_with_controls(data, "treatment_age_50", 1930)

placebo_tidy <- tidy(placebo_model, conf.int = TRUE) %>%
  filter(grepl("time_to_treat", term)) %>%
  mutate(
    time = as.numeric(sub("time_to_treat::(-?\\d+):treat", "\\1", term)),
    significant = p.value < 0.05,
    cohort = "Placebo 1930"
  ) %>%
  filter(!is.na(time), time >= -20, time <= 30)

#------------------------------------------
# 10. Plot Placebo vs 1960s
#------------------------------------------
placebo_compare <- combined_event %>%
  filter(cohort == "1960s") %>%
  bind_rows(placebo_tidy)

ggplot(placebo_compare, aes(x = time, y = estimate, color = cohort)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "black") +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = cohort), alpha = 0.05, color = NA) +
  geom_point(aes(shape = significant)) +
  scale_x_continuous(breaks = seq(-20, 30, 10)) +
  labs(
    title = "Event Study vs Placebo (Log Population)",
    subtitle = "Placebo (1930) vs Actual Treatment (1960s)",
    x = "Years Relative to Treatment",
    y = "Estimated Effect"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

#------------------------------------------
# 11. Print Pre-Trend & Model Summaries
#------------------------------------------
cat("\n=== Pre-Trend Tests ===\n")
for (i in seq_along(pretrend_tests)) {
  if (!is.null(pretrend_tests[[i]])) {
    cat("\nCohort:", cohorts$name[i], "\n")
    print(pretrend_tests[[i]])
  }
}

cat("\n=== Model Summaries ===\n")
for (i in seq_along(event_models$model)) {
  cat("\n-----------------------\n")
  cat("Cohort:", cohorts$name[i], "\n")
  print(summary(event_models$model[[i]]))
}

cat("\n=== Placebo Model Summary ===\n")
print(summary(placebo_model))

#------------------------------------------
# 12. Callaway & Sant’Anna Estimation
#------------------------------------------
data <- data %>%
  mutate(
    treatment_year = case_when(
      treatment_age_50 == 1 ~ 1940,
      treatment_age_60 == 1 ~ 1950,
      treatment_age_70 == 1 ~ 1960,
      is_treated == 0 ~ 0,
      TRUE ~ NA_real_
    )
  )

cs_data <- data %>%
  filter(!is.na(treatment_year)) %>%
  mutate(
    id = municipality_ID,
    G = ifelse(is_treated == 1, treatment_year, 0),
    Y = log(population + 1),
    T = year
  )

att_gt_obj <- att_gt(
  yname = "Y",
  tname = "T",
  idname = "id",
  gname = "G",
  xformla = ~ log_area + sqrt_urban_dist + river_distance_m + Latitude,
  data = cs_data,
  panel = TRUE,
  est_method = "dr"
)

# Dynamic ATT
es_did <- aggte(att_gt_obj, type = "dynamic")
ggdid(es_did) +
  labs(
    title = "Callaway & Sant’Anna Event Study",
    subtitle = "Dynamic ATT by Relative Time",
    x = "Years Relative to Treatment",
    y = "ATT on Log Population"
  )
print(es_did)

#------------------------------------------
# 13. Overall ATT
#------------------------------------------
overall_att <- aggte(att_gt_obj, type = "simple")
print(overall_att)
