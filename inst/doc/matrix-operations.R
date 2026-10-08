## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load-packages------------------------------------------------------------
library(tidymatrix)
library(dplyr, warn.conflicts = FALSE)

tm <- tidymatrix(big5_responses, big5_respondents, big5_items)

## ----add-stats-rows-----------------------------------------------------------
tm <- tm |>
  activate(rows) |>
  add_stats(mean, sd, .names = c("resp_mean", "resp_sd"))

tm |>
  activate(rows) |>
  arrange(resp_sd) |>
  select(respondent_id, resp_mean, resp_sd, completion_min)

## ----add-stats-speed----------------------------------------------------------
tm |>
  activate(rows) |>
  arrange(completion_min) |>
  select(respondent_id, resp_sd, completion_min)

## ----clean--------------------------------------------------------------------
tm <- tm |>
  activate(rows) |>
  filter(completion_min > 3.5)

## ----add-stats-cols-----------------------------------------------------------
tm |>
  activate(columns) |>
  add_stats(.fns = list(
    item_mean = mean,
    item_sd = sd,
    pct_agree = \(x) mean(x >= 4)
  )) |>
  select(item_id, item_text, item_mean, item_sd, pct_agree)

## ----reverse------------------------------------------------------------------
tm_scored <- tm |>
  activate(rows) |>
  transform_matrix(\(x, flip) ifelse(flip, 6L - x, x), flip = reversed)

# E4 is positively keyed, E5 is reverse-keyed
tm$matrix[1:5, c("E4", "E5")]
tm_scored$matrix[1:5, c("E4", "E5")]

## ----reverse-check------------------------------------------------------------
round(cor(tm$matrix[, paste0("E", 1:6)])[5:6, 1:4], 2)
round(cor(tm_scored$matrix[, paste0("E", 1:6)])[5:6, 1:4], 2)

## ----transform-matrix---------------------------------------------------------
tm_scored |>
  activate(matrix) |>
  transform_matrix(\(m) (m - 1) / 4 * 100) |>
  pull_active() |>
  head(3)

## ----transform-cols-----------------------------------------------------------
tm_scored |>
  activate(columns) |>
  transform_matrix(\(x) rank(x) / length(x)) |>
  activate(matrix) |>
  pull_active() |>
  head(3) |>
  round(2)

## ----transform-args-----------------------------------------------------------
tm_scored |>
  activate(columns) |>
  transform_matrix(\(x) x / 4) |>
  activate(matrix) |>
  transform_matrix(round, digits = 1) |>
  pull_active() |>
  head(3)

## ----scale-cols---------------------------------------------------------------
tm_z <- tm_scored |>
  activate(columns) |>
  scale()

round(colMeans(tm_z$matrix)[1:6], 3)
round(apply(tm_z$matrix, 2, sd)[1:6], 3)

## ----center-rows--------------------------------------------------------------
tm_ipsative <- tm |>
  activate(rows) |>
  center()

round(rowMeans(tm_ipsative$matrix)[1:5], 3)

## ----clip---------------------------------------------------------------------
range(tm_z$matrix)

tm_clipped <- tm_z |>
  activate(matrix) |>
  clip_values(min = -2, max = 2)

range(tm_clipped$matrix)

## ----log----------------------------------------------------------------------
counts <- tidymatrix(
  matrix(c(0, 3, 12, 150, 7, 0, 1024, 45, 2), nrow = 3),
  data.frame(gene = c("g1", "g2", "g3")),
  data.frame(sample = c("s1", "s2", "s3"))
)

counts |>
  activate(matrix) |>
  log_transform(base = 2, offset = 1) |>
  pull_active() |>
  round(2)

## ----transpose----------------------------------------------------------------
tm_t <- t(tm_scored |> activate(rows))

dim(tm_scored$matrix)
dim(tm_t$matrix)

active(tm_t)
names(tm_t$row_data)

## ----extract------------------------------------------------------------------
m <- tm_scored |>
  activate(matrix) |>
  pull_active()

class(m)
dim(m)

## ----invalidation-------------------------------------------------------------
tm_pca <- tm_scored |>
  activate(columns) |>
  compute_prcomp(n_components = 2)

list_analyses(tm_pca)

tm_pca <- tm_pca |>
  activate(columns) |>
  scale()

list_analyses(tm_pca)

