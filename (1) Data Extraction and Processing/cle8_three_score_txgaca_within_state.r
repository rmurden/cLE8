# =============================================================================
# County-Level LE8 (cLE8) Three-Score Analysis
# States: Texas (TX), Georgia (GA), California (CA)
# Data:   CDC PLACES 2024 release
#
# Purpose:
#   Compute three versions of the cLE8 score for each county using the
#   within-state distribution:
#
#     Score 1 (point):  ntile scoring using published point estimate prevalences
#     Score 2 (lower):  ntile scoring using lower bounds of published 95% CIs
#     Score 3 (upper):  ntile scoring using upper bounds of published 95% CIs
#
#   
# Components (all from CDC PLACES 2024, all inverse-scored):
#   1. Physical Activity  (LPA_AdjPrev)
#   2. Nicotine Exposure  (CSMOKING_AdjPrev)
#   3. Sleep Health       (SLEEP_AdjPrev)
#   4. Diet               (FOODINSECU_AdjPrev - food insecurity proxy)
#   5. Blood Pressure     (BPHIGH_AdjPrev)
#   6. Blood Lipids       (HIGHCHOL_AdjPrev)
#   7. Blood Glucose      (DIABETES_AdjPrev)
#   8. Body Weight        (OBESITY_AdjPrev)
#
# Output:
#   cle8_three_scores_txgaca_within_state.csv
#
#
# =============================================================================

library(tidyverse); library(readr)

# =============================================================================
# SECTION 1: CONFIGURATION (set path to cLE8 Github)
# =============================================================================

places_2024_path <- "Data/PLACES__County_Data_(GIS_Friendly_Format),_2024_release_20260622.csv"

target_states <- c("Texas", "Georgia", "California")

state_abbr_map <- c(
  "Texas"      = "TX",
  "Georgia"    = "GA",
  "California" = "CA"
)

stroke_belt_class <- c(
  "TX" = "partial",
  "GA" = "full",
  "CA" = "none"
)

# =============================================================================
# SECTION 2: COMPONENT METADATA
# =============================================================================

component_meta <- list(
  lpa = list(
    label     = "Physical Activity",
    point_col = "LPA_AdjPrev",
    ci_col    = "LPA_Adj95CI",
    invert    = TRUE,
    domain    = "behavior"
  ),
  smoking = list(
    label     = "Nicotine Exposure",
    point_col = "CSMOKING_AdjPrev",
    ci_col    = "CSMOKING_Adj95CI",
    invert    = TRUE,
    domain    = "behavior"
  ),
  sleep = list(
    label     = "Sleep Health",
    point_col = "SLEEP_AdjPrev",
    ci_col    = "SLEEP_Adj95CI",
    invert    = TRUE,
    domain    = "behavior"
  ),
  diet = list(
    label     = "Diet (Food Insecurity)",
    point_col = "FOODINSECU_AdjPrev",
    ci_col    = "FOODINSECU_Adj95CI",
    invert    = TRUE,
    domain    = "behavior"
  ),
  bp = list(
    label     = "Blood Pressure",
    point_col = "BPHIGH_AdjPrev",
    ci_col    = "BPHIGH_Adj95CI",
    invert    = TRUE,
    domain    = "factor"
  ),
  chol = list(
    label     = "Blood Lipids",
    point_col = "HIGHCHOL_AdjPrev",
    ci_col    = "HIGHCHOL_Adj95CI",
    invert    = TRUE,
    domain    = "factor"
  ),
  glucose = list(
    label     = "Blood Glucose",
    point_col = "DIABETES_AdjPrev",
    ci_col    = "DIABETES_Adj95CI",
    invert    = TRUE,
    domain    = "factor"
  ),
  obesity = list(
    label     = "Body Weight",
    point_col = "OBESITY_AdjPrev",
    ci_col    = "OBESITY_Adj95CI",
    invert    = TRUE,
    domain    = "factor"
  )
)

comp_names     <- names(component_meta)
behavior_comps <- comp_names[sapply(comp_names, function(x) component_meta[[x]]$domain == "behavior")]
factor_comps   <- comp_names[sapply(comp_names, function(x) component_meta[[x]]$domain == "factor")]

# =============================================================================
# SECTION 3: HELPER FUNCTIONS
# =============================================================================

parse_places_ci <- function(ci_vec) {
  clean <- str_remove_all(ci_vec, "[()]")
  mat   <- str_split_fixed(clean, ",", n = 2)
  list(
    lo = as.numeric(str_trim(mat[, 1])),
    hi = as.numeric(str_trim(mat[, 2]))
  )
}

score_component <- function(x, invert = TRUE) {
  r <- ntile(x, 10)
  if (invert) (11L - r) * 10L else r * 10L
}

composite_score <- function(score_mat) {
  rowMeans(score_mat, na.rm = FALSE)
}

subscore <- function(score_mat, cols) {
  rowMeans(score_mat[, cols, drop = FALSE], na.rm = FALSE)
}

cvh_category <- function(score) {
  case_when(
    score >= 80 ~ "High",
    score >= 50 ~ "Moderate",
    score >= 0  ~ "Low",
    TRUE        ~ NA_character_
  )
}

# =============================================================================
# SECTION 4: LOAD AND SUBSET DATA
# =============================================================================

cat("Loading 2024 PLACES data...\n")

places_raw <- read.csv(
  places_2024_path,
  show_col_types = FALSE,
  locale = locale(encoding = "UTF-8")
)

names(places_raw) <- str_remove(names(places_raw), "^\xef\xbb\xbf")

df <- places_raw |>
  filter(StateDesc %in% target_states) |>
  mutate(StateAbbr = state_abbr_map[StateDesc]) |>
  arrange(StateAbbr, CountyFIPS)

cat(sprintf("  %d counties loaded\n", nrow(df)))

# =============================================================================
# SECTION 5: EXTRACT POINT ESTIMATES AND CI BOUNDS
# =============================================================================

cat("\nExtracting prevalence values and CI bounds...\n")

n_county <- nrow(df)

pt_mat <- matrix(NA_real_, nrow = n_county, ncol = length(comp_names),
                 dimnames = list(df$CountyFIPS, comp_names))
lo_mat <- matrix(NA_real_, nrow = n_county, ncol = length(comp_names),
                 dimnames = list(df$CountyFIPS, comp_names))
hi_mat <- matrix(NA_real_, nrow = n_county, ncol = length(comp_names),
                 dimnames = list(df$CountyFIPS, comp_names))

for (comp in comp_names) {
  meta              <- component_meta[[comp]]
  pt                <- as.numeric(df[[meta$point_col]])
  parsed            <- parse_places_ci(df[[meta$ci_col]])
  pt_mat[, comp]    <- pt
  lo_mat[, comp]    <- parsed$lo
  hi_mat[, comp]    <- parsed$hi
}

# matrix of CI widths in % points

width_mat <- hi_mat - lo_mat 

# =============================================================================
# SECTION 6: WITHIN-STATE SCORING
# =============================================================================
# ntile() is applied within each state's subset of rows. 
# The three score matrices are built by scoring each state
# separately and binding results back into full-dataset order.

cat("\nComputing within-state scores for each scenario...\n")

score_point <- matrix(NA_real_, nrow = n_county, ncol = length(comp_names),
                      dimnames = list(df$CountyFIPS, comp_names))
score_lower <- matrix(NA_real_, nrow = n_county, ncol = length(comp_names),
                      dimnames = list(df$CountyFIPS, comp_names))
score_upper <- matrix(NA_real_, nrow = n_county, ncol = length(comp_names),
                      dimnames = list(df$CountyFIPS, comp_names))

for (state in c("TX", "GA", "CA")) {
  idx <- which(df$StateAbbr == state)
  n_state <- length(idx)
  cat(sprintf("  Scoring %s: %d counties\n", state, n_state))

  for (comp in comp_names) {
    inv <- component_meta[[comp]]$invert
    score_point[idx, comp] <- score_component(pt_mat[idx, comp], invert = inv)
    score_lower[idx, comp] <- score_component(lo_mat[idx, comp], invert = inv)
    score_upper[idx, comp] <- score_component(hi_mat[idx, comp], invert = inv)
  }
}

# =============================================================================
# SECTION 7: COMPOSITE AND SUBSCORES
# =============================================================================

cat("Computing composite and domain subscores...\n")

cle8_point <- composite_score(score_point)
cle8_lower <- composite_score(score_lower)
cle8_upper <- composite_score(score_upper)

behavior_point <- subscore(score_point, behavior_comps)
behavior_lower <- subscore(score_lower, behavior_comps)
behavior_upper <- subscore(score_upper, behavior_comps)

factor_point <- subscore(score_point, factor_comps)
factor_lower <- subscore(score_lower, factor_comps)
factor_upper <- subscore(score_upper, factor_comps)

# cle8_range measures how much the score moves numerically

cle8_range     <- abs(cle8_upper - cle8_lower)
behavior_range <- abs(behavior_upper - behavior_lower)
factor_range   <- abs(factor_upper - factor_lower)

cat_point  <- cvh_category(cle8_point)
cat_lower  <- cvh_category(cle8_lower)
cat_upper  <- cvh_category(cle8_upper)

# cat_stable is TRUE if a county's CVH category 
# (Low, Moderate, High) is the same across all
# 3 score scenarios

cat_stable <- (cat_point == cat_lower) & (cat_point == cat_upper)

# =============================================================================
# SECTION 8: ASSEMBLE OUTPUT
# =============================================================================

cat("Assembling output dataset...\n")

output <- df |>
  select(CountyFIPS, StateDesc, StateAbbr, CountyName) |>
  mutate(
    scoring_referent = "within_state",
    stroke_belt      = stroke_belt_class[StateAbbr],

    cle8_point    = cle8_point,
    cle8_lower    = cle8_lower,
    cle8_upper    = cle8_upper,
    cle8_range    = cle8_range,

    cvh_cat_point  = cat_point,
    cvh_cat_lower  = cat_lower,
    cvh_cat_upper  = cat_upper,
    cvh_cat_stable = cat_stable,

    behavior_point = behavior_point,
    behavior_lower = behavior_lower,
    behavior_upper = behavior_upper,
    behavior_range = behavior_range,

    factor_point = factor_point,
    factor_lower = factor_lower,
    factor_upper = factor_upper,
    factor_range = factor_range
  )

for (comp in comp_names) {
  output[[paste0(comp, "_point")]]      <- score_point[, comp]
  output[[paste0(comp, "_lower")]]      <- score_lower[, comp]
  output[[paste0(comp, "_upper")]]      <- score_upper[, comp]
  output[[paste0(comp, "_range")]]      <- abs(score_upper[, comp] - score_lower[, comp])
  output[[paste0(comp, "_ci_width_pp")]] <- width_mat[, comp]
}

cat(sprintf("  Output: %d rows x %d columns\n", nrow(output), ncol(output)))

# =============================================================================
# SECTION 9: DIAGNOSTICS
# =============================================================================

cat("\n=== WITHIN-STATE SUMMARY STATISTICS ===\n\n")

for (state in c("TX", "GA", "CA")) {
  sub <- output |> filter(StateAbbr == state)
  cat(sprintf("%s — %s stroke belt (%d counties):\n",
              state, stroke_belt_class[state], nrow(sub)))

  cat(sprintf("  %-12s  %6s  %6s  %6s  %8s\n",
              "", "Point", "Lower", "Upper", "Range"))
  cat("  ", strrep("-", 48), "\n", sep = "")

  for (score_name in c("cle8", "behavior", "factor")) {
    pt  <- mean(sub[[paste0(score_name, "_point")]], na.rm = TRUE)
    lo  <- mean(sub[[paste0(score_name, "_lower")]], na.rm = TRUE)
    hi  <- mean(sub[[paste0(score_name, "_upper")]], na.rm = TRUE)
    rng <- mean(sub[[paste0(score_name, "_range")]], na.rm = TRUE)
    cat(sprintf("  %-12s  %6.1f  %6.1f  %6.1f  %8.1f\n",
                score_name, pt, lo, hi, rng))
  }

  cat("  CVH category (point score): ")
  print(table(sub$cvh_cat_point))

  n_stable <- sum(sub$cvh_cat_stable, na.rm = TRUE)
  cat(sprintf("  Category stable across 3 scores: %d / %d (%.1f%%)\n\n",
              n_stable, nrow(sub), 100 * n_stable / nrow(sub)))
}


# =============================================================================
# SECTION 10: Load 2024 CHR&R data and extract % 65 and Older
# =============================================================================
chrr_2024_path <- here::here("analytic_data2024.csv")

chrr_2024<-read.csv(chrr_2024_path)


# Check column names
names(chrr_2024)
names(chrr_2024)[1:10]

###############################
# % Population Age 65 and Older
###############################

# Load required packages
library(dplyr)
library(stringr)

# Check the exact variable names related to age 65+
grep("65",  names(chrr_2024), value = TRUE, ignore.case = TRUE)

# Check the first 30 unique values of the 65+ raw variable
unique( chrr_2024[["% 65 and Older raw value"]])[1:30]

# Check how many missing values are already present
table(is.na(chrr_2024[["% 65 and Older raw value"]]))

# Look at the first 20 observations
head(chrr_2024 %>%select( `5-digit FIPS Code`,`% 65 and Older raw value`),20)


###############################
# Create 65+ dataset
###############################

pct_65plus <- chrr_2024 %>%
  select(County.Code = `5-digit FIPS Code`,pct_65_older = `% 65 and Older raw value`) %>%
  mutate(County.Code = str_pad(as.character(County.Code),width = 5,side = "left",pad = "0"),
         pct_65_older = as.numeric(str_replace(as.character(pct_65_older), "%", "")))

view(pct_65plus)

###############################
# Check resulting dataset
###############################

# Variable names
names(pct_65plus)

# Number of rows and columns
dim(pct_65plus)

# First 10 counties
head(pct_65plus, 10)

# Summary of 65+ variable
summary(pct_65plus$pct_65_older)

# Number of missing values after conversion
sum(is.na(pct_65plus$pct_65_older))

# Check for duplicated County.Code
sum(duplicated(pct_65plus$County.Code))

# Check structure
str(pct_65plus)

 

# =============================================================================
# SECTION 11: EXPORT
# =============================================================================


# Save diagnostic summaries to CSV

# State-level score summaries
state_summary <- output |>
  group_by(StateAbbr) |>
  summarise(
    n_counties = n(),
    cle8_mean  = round(mean(cle8_point, na.rm = TRUE), 2),
    cle8_median  = round(median(cle8_point, na.rm = TRUE), 2),
    cle8_min = round(min(cle8_point, na.rm = TRUE), 2),
    cle8_max = round(max(cle8_point, na.rm = TRUE), 2),
    range_mean = round(mean(cle8_range, na.rm = TRUE), 2),
    range_median = round(median(cle8_range, na.rm = TRUE), 2),
    n_stable = sum(cvh_cat_stable, na.rm = TRUE),
    pct_stable = round(100 * mean(cvh_cat_stable, na.rm = TRUE), 1),
    n_low  = sum(cvh_cat_point == "Low", na.rm = TRUE),
    n_moderate = sum(cvh_cat_point == "Moderate", na.rm = TRUE),
    n_high = sum(cvh_cat_point == "High", na.rm = TRUE),
    pct_low  = round(100 * mean(cvh_cat_point == "Low", na.rm = TRUE), 1),
    pct_moderate = round(100 * mean(cvh_cat_point == "Moderate", na.rm = TRUE), 1),
    pct_high = round(100 * mean(cvh_cat_point == "High", na.rm = TRUE), 1)
  )




write_csv(state_summary, "diagnostics_state_summary.csv")

write_csv(output, "cle8_three_scores_txgaca_within_state.csv")




