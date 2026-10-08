## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)

## ----load-packages------------------------------------------------------------
library(tidymatrix)
library(dplyr, warn.conflicts = FALSE)

## ----data---------------------------------------------------------------------
tm <- tidymatrix(big5_responses, big5_respondents, big5_items) |>
  activate(rows) |>
  filter(completion_min > 3.5) |>
  transform_matrix(\(x, flip) ifelse(flip, 6L - x, x), flip = reversed)

tm

## ----pca-rows-----------------------------------------------------------------
tm <- tm |>
  activate(rows) |>
  compute_prcomp(n_components = 5)

tm |>
  activate(rows) |>
  select(respondent_id, row_pca_PC1:row_pca_PC5)

## ----pca-variance-------------------------------------------------------------
pca <- get_analysis(tm, "row_pca")
round(summary(pca)$importance[, 1:7], 3)

## ----pca-loadings-------------------------------------------------------------
tm <- tm |>
  activate(columns) |>
  mutate(
    loading_PC1 = pca$rotation[item_id, "PC1"],
    loading_PC2 = pca$rotation[item_id, "PC2"]
  )

tm |>
  activate(columns) |>
  as_tibble() |>
  group_by(trait) |>
  summarise(across(c(loading_PC1, loading_PC2), mean))

## ----pca-scores---------------------------------------------------------------
tm |>
  activate(rows) |>
  as_tibble() |>
  summarise(across(
    row_pca_PC1:row_pca_PC3,
    \(pc) cor(pc, life_satisfaction)
  ))

## ----pca-cols-----------------------------------------------------------------
tm <- tm |>
  activate(columns) |>
  compute_prcomp(n_components = 2)

tm |>
  activate(columns) |>
  as_tibble() |>
  group_by(trait) |>
  summarise(
    PC1 = mean(column_pca_PC1),
    PC2 = mean(column_pca_PC2)
  )

## ----hclust-items-------------------------------------------------------------
tm <- tm |>
  activate(columns) |>
  compute_hclust(k = 5, method = "ward.D2", name = "item_clusters")

tm |>
  activate(columns) |>
  as_tibble() |>
  count(trait, item_clusters_cluster)

## ----dendrogram, fig.height = 5-----------------------------------------------
hc <- get_analysis(tm, "item_clusters")
plot(hc, main = "Items", xlab = "", sub = "")

## ----hclust-raw---------------------------------------------------------------
tidymatrix(big5_responses, big5_respondents, big5_items) |>
  activate(rows) |>
  filter(completion_min > 3.5) |>
  activate(columns) |>
  compute_hclust(k = 5, method = "ward.D2") |>
  as_tibble() |>
  count(trait, column_hclust_cluster)

## ----hclust-rows--------------------------------------------------------------
tm <- tm |>
  activate(rows) |>
  compute_hclust(k = 4, method = "ward.D2", name = "resp_hclust")

tm |>
  activate(rows) |>
  as_tibble() |>
  count(resp_hclust_cluster)

## ----kmeans-------------------------------------------------------------------
set.seed(42)
tm <- tm |>
  activate(rows) |>
  compute_kmeans(centers = 3, nstart = 25, name = "resp_kmeans")

km <- get_analysis(tm, "resp_kmeans")
km$size

## ----profile------------------------------------------------------------------
profile <- tm |>
  activate(rows) |>
  group_by(resp_kmeans_cluster) |>
  summarise(
    n = n(),
    age = mean(age),
    life_satisfaction = mean(life_satisfaction)
  ) |>
  activate(columns) |>
  group_by(trait) |>
  summarise()

m <- profile$matrix
dimnames(m) <- list(
  paste("cluster", profile$row_data$resp_kmeans_cluster),
  profile$col_data$trait
)

profile$row_data
round(m, 2)

## ----compare------------------------------------------------------------------
tm |>
  activate(rows) |>
  as_tibble() |>
  with(table(hclust = resp_hclust_cluster, kmeans = resp_kmeans_cluster))

## ----manage-------------------------------------------------------------------
list_analyses(tm)

check_analyses(tm)

tm_small <- remove_analysis(tm, "resp_hclust")
list_analyses(tm_small)

## ----invalidate---------------------------------------------------------------
tm_filtered <- tm |>
  activate(rows) |>
  filter(age >= 30)

list_analyses(tm_filtered)

tm_filtered |>
  activate(rows) |>
  select(respondent_id, row_pca_PC1, resp_kmeans_cluster)

