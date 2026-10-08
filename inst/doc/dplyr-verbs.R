## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load-packages------------------------------------------------------------
library(tidymatrix)
library(dplyr, warn.conflicts = FALSE)

tm <- tidymatrix(big5_responses, big5_respondents, big5_items)

## ----activate-----------------------------------------------------------------
tm |> activate(rows) |> active()
tm |> activate(columns)

## ----filter-rows--------------------------------------------------------------
tm_clean <- tm |>
  activate(rows) |>
  filter(completion_min > 3.5)

nrow(tm$matrix)
nrow(tm_clean$matrix)

## ----filter-rows-2------------------------------------------------------------
tm_clean |>
  activate(rows) |>
  filter(age >= 30, age < 40, education %in% c("Master", "PhD"))

## ----filter-cols--------------------------------------------------------------
tm_clean |>
  activate(columns) |>
  filter(trait == "Openness", !reversed)

## ----mutate-rows--------------------------------------------------------------
tm_clean <- tm_clean |>
  activate(rows) |>
  mutate(
    age_group = cut(
      age,
      breaks = c(17, 29, 49, 64, Inf),
      labels = c("18-29", "30-49", "50-64", "65+")
    ),
    satisfied = life_satisfaction >= 7
  )

tm_clean

## ----mutate-cols--------------------------------------------------------------
tm_clean <- tm_clean |>
  activate(columns) |>
  mutate(label = paste0(item_id, if_else(reversed, " (R)", "")))

tm_clean |>
  activate(columns) |>
  pull(label)

## ----select-------------------------------------------------------------------
tm_clean |>
  activate(rows) |>
  select(respondent_id, age, gender, education)

tm_clean |>
  activate(columns) |>
  rename(text = item_text) |>
  relocate(label, .after = item_id)

## ----pull---------------------------------------------------------------------
tm_clean |>
  activate(rows) |>
  pull(occupation) |>
  table()

## ----arrange-cols-------------------------------------------------------------
tm_ordered <- tm_clean |>
  activate(columns) |>
  arrange(position)

colnames(tm_ordered$matrix)[1:10]

## ----arrange-rows-------------------------------------------------------------
tm_by_age <- tm_clean |>
  activate(rows) |>
  arrange(desc(age))

tm_by_age |> activate(rows) |> select(respondent_id, age, occupation)
rownames(tm_by_age$matrix)[1:5]

## ----slice--------------------------------------------------------------------
# first five respondents
tm_clean |> activate(rows) |> slice(1:5)

# first three items
tm_clean |> activate(columns) |> slice_head(n = 3)

# a random 10% sample of respondents
set.seed(1)
tm_clean |> activate(rows) |> slice_sample(prop = 0.1)

## ----score--------------------------------------------------------------------
tm_scored <- tm_clean |>
  activate(rows) |>
  transform_matrix(\(x, flip) ifelse(flip, 6L - x, x), flip = reversed)

trait_scores <- tm_scored |>
  activate(columns) |>
  group_by(trait) |>
  summarise(n_items = n())

trait_scores
dim(trait_scores$matrix)
round(head(trait_scores$matrix), 2)

## ----summarise-rows-----------------------------------------------------------
by_age <- tm_scored |>
  activate(rows) |>
  group_by(age_group) |>
  summarise(n = n(), mean_age = mean(age))

by_age

## ----summarise-both-----------------------------------------------------------
age_trait <- by_age |>
  activate(columns) |>
  group_by(trait) |>
  summarise()

m <- age_trait$matrix
dimnames(m) <- list(
  age_trait$row_data$age_group,
  age_trait$col_data$trait
)
round(m, 2)

## ----matrix-fn----------------------------------------------------------------
tm_scored |>
  activate(rows) |>
  group_by(gender) |>
  summarise(n = n(), .matrix_fn = median)

## ----count--------------------------------------------------------------------
tm_clean |>
  activate(rows) |>
  count(education)

tm_clean |>
  activate(columns) |>
  group_by(trait, reversed) |>
  tally()

## ----chain--------------------------------------------------------------------
tm_clean |>
  activate(rows) |>
  filter(
    !occupation %in% c("Student", "Retired"),
    education %in% c("Bachelor", "Master", "PhD")
  ) |>
  select(respondent_id, age, gender, education, occupation) |>
  activate(columns) |>
  filter(trait %in% c("Conscientiousness", "Neuroticism")) |>
  arrange(position) |>
  mutate(short = substr(item_text, 1, 25))

## ----invalidation-------------------------------------------------------------
tm_pca <- tm_scored |>
  activate(columns) |>
  compute_prcomp(n_components = 2)

list_analyses(tm_pca)

tm_pca_filtered <- tm_pca |>
  activate(columns) |>
  filter(trait != "Openness")

list_analyses(tm_pca_filtered)

