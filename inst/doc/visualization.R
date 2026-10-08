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

## ----load-ggplot2, eval = has_ggplot2-----------------------------------------
library(ggplot2)
theme_set(theme_minimal())

## ----data---------------------------------------------------------------------
tm <- tidymatrix(big5_responses, big5_respondents, big5_items) |>
  activate(rows) |>
  filter(completion_min > 3.5) |>
  transform_matrix(\(x, flip) ifelse(flip, 6L - x, x), flip = reversed) |>
  compute_across(
    \(answers, items) tapply(answers, items$trait, mean),
    add_to_data = TRUE,
    prefix = "score"
  ) |>
  mutate(age_group = cut(
    age,
    breaks = c(17, 29, 49, 64, Inf),
    labels = c("18-29", "30-49", "50-64", "65+")
  )) |>
  compute_prcomp(n_components = 2, store = FALSE) |>
  activate(columns) |>
  compute_prcomp(n_components = 2, store = FALSE)

## ----pca-rows, eval = has_ggplot2---------------------------------------------
tm |>
  activate(rows) |>
  as_tibble() |>
  ggplot(aes(row_pca_PC1, row_pca_PC2, colour = score_Neuroticism)) +
  geom_point() +
  scale_colour_viridis_c() +
  labs(x = "PC1", y = "PC2", colour = "Neuroticism")

## ----pca-cols, eval = has_ggplot2---------------------------------------------
tm |>
  activate(columns) |>
  as_tibble() |>
  ggplot(aes(column_pca_PC1, column_pca_PC2, colour = trait)) +
  geom_text(aes(label = item_id)) +
  labs(x = "PC1", y = "PC2", colour = NULL)

## ----scores-age, eval = has_ggplot2-------------------------------------------
tm |>
  activate(rows) |>
  as_tibble() |>
  ggplot(aes(age, score_Conscientiousness)) +
  geom_jitter(height = 0.05, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x) +
  labs(x = "Age", y = "Conscientiousness score")

## ----volcano, eval = has_ggplot2----------------------------------------------
tm |>
  activate(rows) |>
  filter(gender %in% c("Female", "Male")) |>
  activate(columns) |>
  compute_ttest(group_col = "gender", control = "Male", treatment = "Female") |>
  ggplot(aes(log2fc, -log10(p.value), colour = trait)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_text(aes(label = item_id)) +
  labs(
    x = "log2(mean women / mean men)",
    y = "-log10(p-value)",
    colour = NULL
  )

## ----to-long------------------------------------------------------------------
long <- to_long(tm)
dim(long)
long |> select(respondent_id, gender, age_group, item_id, trait, value)

## ----long-bar, eval = has_ggplot2, fig.height = 5-----------------------------
long |>
  filter(trait == "Neuroticism", gender != "Non-binary") |>
  count(item_id, gender, value) |>
  group_by(item_id, gender) |>
  mutate(share = n / sum(n)) |>
  ggplot(aes(factor(value), share, fill = gender)) +
  geom_col(position = "dodge") +
  facet_wrap(~item_id) +
  labs(x = "Answer (re-scored)", y = "Share of respondents", fill = NULL)

## ----long-line, eval = has_ggplot2--------------------------------------------
long |>
  group_by(trait, item_id, age_group) |>
  summarise(mean = mean(value), .groups = "drop") |>
  ggplot(aes(age_group, mean, group = item_id, colour = trait)) +
  geom_line() +
  facet_wrap(~trait, nrow = 1) +
  labs(x = "Age group", y = "Mean answer") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45))

## ----heatmap, eval = has_pheatmap, fig.height = 7-----------------------------
tm |>
  activate(columns) |>
  scale() |>
  activate(matrix) |>
  clip_values(min = -2, max = 2) |>
  activate(columns) |>
  compute_hclust(method = "ward.D2") |>
  activate(rows) |>
  compute_hclust(method = "ward.D2") |>
  plot_pheatmap(
    col_names = "item_id",
    row_annotation = c("age_group", "gender"),
    col_annotation = "trait",
    show_rownames = FALSE,
    treeheight_row = 20
  )

## ----heatmap-agg, eval = has_pheatmap, fig.height = 3-------------------------
tm |>
  activate(rows) |>
  group_by(age_group) |>
  summarise(n = n()) |>
  activate(columns) |>
  arrange(trait, item_id) |>
  plot_pheatmap(
    row_names = "age_group",
    col_names = "item_id",
    row_annotation = FALSE,
    col_annotation = "trait",
    row_cluster = FALSE,
    col_cluster = FALSE,
    display_numbers = TRUE,
    number_format = "%.1f",
    fontsize_number = 6
  )

## ----export-matrix------------------------------------------------------------
m <- tm |>
  activate(matrix) |>
  pull_active()

m[1:3, 1:6]

## ----export-metadata----------------------------------------------------------
tm |>
  activate(rows) |>
  as_tibble() |>
  select(respondent_id, starts_with("score_"), starts_with("row_pca"))

## ----export-file--------------------------------------------------------------
out <- tempfile(fileext = ".csv")

tm |>
  activate(columns) |>
  select(item_id, trait, item_text) |>
  activate(rows) |>
  select(respondent_id, age, gender, country) |>
  to_long() |>
  write.csv(out, row.names = FALSE)

read.csv(out, nrows = 3)

