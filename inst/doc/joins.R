## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load-packages------------------------------------------------------------
library(tidymatrix)
library(dplyr, warn.conflicts = FALSE)

tm <- tidymatrix(big5_responses, big5_respondents, big5_items) |>
  activate(rows) |>
  filter(completion_min > 3.5)

## ----countries-in-data--------------------------------------------------------
tm |> activate(rows) |> pull(country) |> table()

## ----countries----------------------------------------------------------------
big5_countries

## ----left-join----------------------------------------------------------------
tm_left <- tm |>
  activate(rows) |>
  left_join(big5_countries, by = "country")

tm_left

tm_left |>
  activate(rows) |>
  filter(is.na(country_name)) |>
  select(respondent_id, country, country_name)

identical(dim(tm_left$matrix), dim(tm$matrix))

## ----inner-join---------------------------------------------------------------
tm_inner <- tm |>
  activate(rows) |>
  inner_join(big5_countries, by = "country")

nrow(tm$matrix)
nrow(tm_inner$matrix)

## ----inner-join-use-----------------------------------------------------------
tm_inner |>
  activate(rows) |>
  count(region)

## ----full-join----------------------------------------------------------------
tm_full <- tm |>
  activate(rows) |>
  full_join(big5_countries, by = "country")

nrow(tm_full$matrix)

tm_full |>
  activate(rows) |>
  filter(is.na(respondent_id)) |>
  select(respondent_id, country, country_name)

tail(tm_full$matrix[, 1:6], 3)

## ----semi-anti----------------------------------------------------------------
# respondents from countries in the table
tm_semi <- tm |>
  activate(rows) |>
  semi_join(big5_countries, by = "country")

nrow(tm_semi$matrix)
ncol(tm_semi$row_data)  # no new columns

# respondents from countries missing from the table
tm |>
  activate(rows) |>
  anti_join(big5_countries, by = "country") |>
  select(respondent_id, country, age, gender)

## ----trait-table--------------------------------------------------------------
trait_info <- data.frame(
  trait = c("Extraversion", "Agreeableness", "Conscientiousness",
            "Neuroticism", "Openness"),
  abbreviation = c("E", "A", "C", "N", "O"),
  high_pole = c("outgoing, energetic", "friendly, compassionate",
                "organised, dependable", "anxious, moody",
                "curious, imaginative"),
  low_pole = c("reserved, quiet", "critical, detached",
               "careless, spontaneous", "calm, stable",
               "conventional, practical")
)

tm_items <- tm |>
  activate(columns) |>
  left_join(trait_info, by = "trait")

tm_items |>
  activate(columns) |>
  select(item_id, trait, abbreviation, high_pole)

## ----short-form---------------------------------------------------------------
short_form <- data.frame(
  item_id = c("E1", "E5", "A1", "A5", "C1", "C5", "N1", "N5", "O1", "O5")
)

tm_short <- tm |>
  activate(columns) |>
  semi_join(short_form, by = "item_id")

colnames(tm_short$matrix)

## ----different-keys-----------------------------------------------------------
country_codes <- data.frame(
  iso2 = c("EE", "FI", "LV", "LT", "SE", "DE"),
  eu_member_since = c(2004, 1995, 2004, 2004, 1995, 1958)
)

tm |>
  activate(rows) |>
  left_join(country_codes, by = join_by(country == iso2)) |>
  select(respondent_id, country, eu_member_since)

## ----unique-keys--------------------------------------------------------------
anyDuplicated(big5_countries$country) == 0

## ----join-then-analyse--------------------------------------------------------
tm_analysed <- tm |>
  activate(rows) |>
  left_join(big5_countries, by = "country") |>
  activate(columns) |>
  compute_hclust(k = 5, method = "ward.D2")

list_analyses(tm_analysed)

