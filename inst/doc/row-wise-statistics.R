## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load-packages------------------------------------------------------------
library(tidymatrix)
library(dplyr, warn.conflicts = FALSE)

## ----data---------------------------------------------------------------------
tm <- tidymatrix(big5_responses, big5_respondents, big5_items) |>
  activate(rows) |>
  filter(completion_min > 3.5) |>
  transform_matrix(\(x, flip) ifelse(flip, 6L - x, x), flip = reversed)

## ----ttest--------------------------------------------------------------------
tm_fm <- tm |>
  activate(rows) |>
  filter(gender %in% c("Female", "Male"))

ttest <- tm_fm |>
  activate(columns) |>
  compute_ttest(group_col = "gender", control = "Male", treatment = "Female")

ttest |>
  select(item_id, trait, p.value, log2fc, p.adj) |>
  arrange(p.value)

## ----ttest-traits-------------------------------------------------------------
ttest |>
  filter(p.adj < 0.05) |>
  count(trait)

## ----wilcox-------------------------------------------------------------------
tm_fm |>
  activate(columns) |>
  compute_wilcox(group_col = "gender", control = "Male", treatment = "Female") |>
  filter(p.adj < 0.05) |>
  select(item_id, trait, median_diff, p.value, p.adj)

## ----correlation--------------------------------------------------------------
age_cor <- tm |>
  activate(columns) |>
  compute_correlation(var = "age")

age_cor |>
  filter(p.adj < 0.05) |>
  select(item_id, trait, correlation, p.adj) |>
  arrange(correlation)

## ----lm-simple----------------------------------------------------------------
tm |>
  activate(columns) |>
  compute_lm_simple(predictor = "life_satisfaction") |>
  select(item_id, trait, slope, r.squared, p.adj) |>
  arrange(slope)

## ----lm-multiple--------------------------------------------------------------
tm_fm <- tm_fm |>
  activate(rows) |>
  mutate(gender = factor(gender, levels = c("Male", "Female")))

tm_fm |>
  activate(columns) |>
  filter(trait == "Neuroticism") |>
  compute_lm(~ age + gender + education, coef = "genderFemale") |>
  select(item_id, item_text, estimate, se, p.adj)

## ----lm-overall---------------------------------------------------------------
tm_fm |>
  activate(columns) |>
  compute_lm(~ age + gender + education) |>
  select(item_id, trait, r.squared, p.value) |>
  arrange(desc(r.squared)) |>
  head()

## ----anova--------------------------------------------------------------------
tm |>
  activate(columns) |>
  compute_anova(group_col = "education") |>
  filter(p.adj < 0.05) |>
  select(item_id, trait, f.statistic, p.adj)

## ----kruskal------------------------------------------------------------------
tm |>
  activate(columns) |>
  compute_kruskal(group_col = "occupation") |>
  filter(p.adj < 0.05) |>
  select(item_id, trait, statistic, p.adj)

## ----custom-------------------------------------------------------------------
cohens_d <- function(values, respondents) {
  f <- values[respondents$gender == "Female"]
  m <- values[respondents$gender == "Male"]
  pooled_sd <- sqrt(
    ((length(f) - 1) * var(f) + (length(m) - 1) * var(m)) /
      (length(f) + length(m) - 2)
  )
  list(
    mean_female = mean(f),
    mean_male = mean(m),
    d = (mean(f) - mean(m)) / pooled_sd
  )
}

effect_sizes <- tm_fm |>
  activate(columns) |>
  compute_across(cohens_d)

effect_sizes |>
  group_by(trait) |>
  summarise(mean_d = mean(d), min_d = min(d), max_d = max(d))

## ----longstring---------------------------------------------------------------
longstring <- function(answers, items) {
  runs <- rle(answers[order(items$position)])
  list(longstring = max(runs$lengths))
}

tidymatrix(big5_responses, big5_respondents, big5_items) |>
  activate(rows) |>
  compute_across(longstring) |>
  arrange(desc(longstring)) |>
  select(respondent_id, longstring, completion_min) |>
  head(8)

## ----trait-scores-------------------------------------------------------------
tm |>
  activate(rows) |>
  compute_across(\(answers, items) tapply(answers, items$trait, mean)) |>
  select(respondent_id, Agreeableness:Openness)

## ----add-to-data--------------------------------------------------------------
tm_stats <- tm_fm |>
  activate(columns) |>
  compute_ttest(
    group_col = "gender", control = "Male", treatment = "Female",
    add_to_data = TRUE, prefix = "gender"
  ) |>
  compute_correlation(var = "age", add_to_data = TRUE, prefix = "age") |>
  compute_across(cohens_d, add_to_data = TRUE, prefix = "gender")

tm_stats |>
  activate(columns) |>
  select(item_id, trait, gender_p.adj, gender_d, age_correlation, age_p.adj)

## ----filter-by-stats----------------------------------------------------------
tm_stats |>
  activate(columns) |>
  filter(gender_p.adj < 0.05, abs(gender_d) > 0.2)

