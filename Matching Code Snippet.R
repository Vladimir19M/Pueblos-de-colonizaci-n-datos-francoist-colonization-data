
# TFG Code Matching 
# Vladimir Miculaiciuc Kudriavtsev
# The code used for the actual results is not identical, However the functional form remains the same.

#------------------------------------------
# 1. Libraries
#------------------------------------------
library(tidyr)
library(dplyr)
library(readr)
library(MatchIt)
library(cobalt)
library(ggplot2)
library(gridExtra)
library(knitr)
library(kableExtra)

#------------------------------------------
# 2. Data / Load database classify columns
#------------------------------------------
df <- read_csv2("DataBase1.csv")

match_data <- df %>%
  filter(!is.na(Pop_1930)) %>%
  mutate(
    river_distance_m = as.numeric(river_distance_m),
    Latitude = as.numeric(Latitude),
    area_km2 = as.numeric(`area_km^2`),
    urban_dist = as.numeric(urban_dist_km),
    longitude = as.numeric(longitude),
    log_Pop_1930 = log(Pop_1930 + 1),
    log_area = log(area_km2 + 1),
    sqrt_urban_dist = sqrt(urban_dist)
  ) %>%
  drop_na(river_distance_m, Latitude, urban_dist, area_km2, longitude, Pop_1930)

#------------------------------------------
# 3. Propensity score Estimation
#------------------------------------------
ps_model <- glm(
  is_treated ~ log_Pop_1930 + river_distance_m + Latitude + 
    sqrt_urban_dist + log_area + longitude,
  data = match_data,
  family = binomial
)

match_data$pscore <- predict(ps_model, type = "response")

#------------------------------------------
# 4. Trim Data
#------------------------------------------
ps_trim <- quantile(match_data$pscore, probs = c(0.05, 0.95), na.rm = TRUE)

match_data_trimmed <- match_data %>%
  filter(pscore >= ps_trim[1], pscore <= ps_trim[2])

#------------------------------------------
# 5. matching specifications
#------------------------------------------
# Improved Matching: 1:3, caliper = 0.1 ; (Used in the paper)
match_improved <- matchit(
  is_treated ~ log_Pop_1930 + river_distance_m + Latitude + 
    sqrt_urban_dist + log_area + longitude,
  data = match_data_trimmed,
  method = "nearest",
  distance = "logit",
  ratio = 3,
  caliper = 0.1,
  std.caliper = TRUE
)

matched_df_improved <- match.data(match_improved)

# Original 1:1 Matching (for comparison)
match_original <- matchit(
  is_treated ~ log_Pop_1930 + river_distance_m + Latitude + 
    sqrt_urban_dist + log_area + longitude,
  data = match_data_trimmed,
  method = "nearest",
  distance = "logit"
)

#------------------------------------------
# 6. Covariates Balance Assessment
#------------------------------------------
# Improved matching
p_improved <- love.plot(
  match_improved,
  binary = "std",
  thresholds = c(m = 0.1),
  var.order = "unadjusted",
  title = "Balance After Improved Matching (1:3, Caliper = 0.1)",
  colors = c("red", "blue")
) + theme_minimal()

# Original matching
p_original <- love.plot(
  match_original,
  binary = "std",
  thresholds = c(m = 0.1),
  var.order = "unadjusted",
  title = "Balance After Original Matching (1:1)",
  colors = c("red", "blue")
) + theme_minimal()

# comparison
grid.arrange(p_original, p_improved, ncol = 2)

#------------------------------------------
# 7. Different Diagnostics 
#------------------------------------------
# Propensity Score Distribution
ggplot(match_data_trimmed, aes(x = pscore, fill = factor(is_treated))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9"),
                    labels = c("Control", "Treated"),
                    name = "Group") +
  labs(title = "Propensity Score Distribution (Trimmed Sample)",
       x = "Propensity Score", y = "Density") +
  theme_minimal()

# Sample Sizes
cat("Sample sizes after trimming:\n")
print(table(match_data_trimmed$is_treated))
cat("\nMatched sample sizes (improved matching):\n")
print(table(matched_df_improved$is_treated))

#------------------------------------------
#8. Covariate Balance Table
#------------------------------------------
balance_table <- summary(match_improved)$sum.all %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Covariate")

print(kable(balance_table, caption = "Covariate Balance After Improved Matching") %>%
        kable_styling(full_width = FALSE))

#------------------------------------------
#9. Save Outputs
#------------------------------------------
write.csv(matched_df_improved, "Match Pairs1.csv", row.names = FALSE)
