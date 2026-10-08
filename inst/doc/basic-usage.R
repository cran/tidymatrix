## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)
has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)
has_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)

## ----load-packages------------------------------------------------------------
library(tidymatrix)
library(dplyr, warn.conflicts = FALSE)

## ----data-matrix--------------------------------------------------------------
big5_responses[1:5, 1:8]

## ----data-rows----------------------------------------------------------------
head(big5_respondents)

## ----data-cols----------------------------------------------------------------
head(big5_items)

## ----create-------------------------------------------------------------------
tm <- tidymatrix(big5_responses, big5_respondents, big5_items)
tm

## ----filter-rows--------------------------------------------------------------
tm_clean <- tm |>
  activate(rows) |>
  filter(completion_min > 3.5)

tm_clean |> 
  activate(matrix)

## ----filter-cols--------------------------------------------------------------
tm_clean |>
  activate(columns) |>
  filter(trait == "Neuroticism") |>
  arrange(position)

## ----summarise----------------------------------------------------------------
by_country <- tm_clean |>
  activate(rows) |>
  group_by(country) |>
  summarise(n = n(), mean_age = mean(age))

by_country

by_country |>
  activate(matrix) |>
  transform_matrix(round, digits = 2)

## ----countries----------------------------------------------------------------
big5_countries

## ----join---------------------------------------------------------------------
tm_clean |>
  activate(rows) |>
  left_join(big5_countries, by = "country") |>
  select(respondent_id, country, country_name, region)

## ----anti-join----------------------------------------------------------------
tm_clean |>
  activate(rows) |>
  anti_join(big5_countries, by = "country") |>
  pull(country) |>
  table()

## ----reverse------------------------------------------------------------------
tm_scored <- tm_clean |>
  activate(rows) |>
  transform_matrix(\(x, flip) ifelse(flip, 6L - x, x), flip = reversed)

## ----matrix-ops---------------------------------------------------------------
tm_scored |>
  activate(columns) |>
  scale() |>                 # z-score each item
  activate(matrix) |>
  clip_values(min = -2, max = 2) |>
  activate(rows) |>
  add_stats(mean, sd)        # per-respondent mean and sd into row metadata

## ----pca----------------------------------------------------------------------
tm_items <- tm_scored |>
  activate(columns) |>
  compute_prcomp(n_components = 2) |>
  compute_hclust(k = 5, method = "ward.D2") |>
  compute_mds(k = 2)

tm_items |>
  activate(columns) |>
  select(item_id, trait, column_pca_PC1, column_pca_PC2, column_hclust_cluster)

list_analyses(tm_items)

## ----pca-check----------------------------------------------------------------
tm_items |>
  activate(columns) |>
  as_tibble() |>
  count(trait, column_hclust_cluster)

## ----ttest--------------------------------------------------------------------
gender_tests <- tm_scored |>
  activate(rows) |>
  filter(gender %in% c("Female", "Male")) |>
  activate(columns) |>
  compute_ttest(group_col = "gender", control = "Male", treatment = "Female")

gender_tests |>
  filter(p.adj < 0.05) |>
  select(item_id, trait, item_text, p.value, p.adj)

## ----plot-pca, eval = has_ggplot2---------------------------------------------
library(ggplot2)

tm_items |>
  activate(columns) |>
  as_tibble() |>
  ggplot(aes(column_pca_PC1, column_pca_PC2, colour = trait, label = item_id)) +
  geom_text() +
  labs(x = "PC1", y = "PC2", colour = NULL) +
  theme_minimal()

## ----to-long------------------------------------------------------------------
tm_scored |>
  to_long() |>
  select(respondent_id, gender, item_id, trait, value) |>
  head()

## ----heatmap, eval = has_pheatmap, fig.height = 6-----------------------------
tm_scored |>
  activate(columns) |>
  scale() |>
  compute_hclust(method = "ward.D2") |>
  activate(rows) |>
  compute_hclust(method = "ward.D2") |>
  plot_pheatmap(
    col_names = "item_id",
    row_annotation = c("gender", "age"),
    col_annotation = "trait",
    show_rownames = FALSE
  )

