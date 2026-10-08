library(dplyr)

# Regression tests for bugs that silently produced wrong results.

test_that("count(sort = TRUE) permutes the matrix along with the counts", {
  # Group means are deliberately distinct so a desynchronised matrix is visible:
  # group "a" has 2 columns of 1, "b" 3 columns of 10, "c" 4 columns of 100.
  mat <- matrix(c(rep(1, 2 * 2), rep(10, 2 * 3), rep(100, 2 * 4)), nrow = 2)
  tm <- tidymatrix(
    mat,
    data.frame(r = 1:2),
    data.frame(g = c(rep("a", 2), rep("b", 3), rep("c", 4)))
  )

  result <- tm |> activate(columns) |> count(g, sort = TRUE)

  expect_equal(result$col_data$g, c("c", "b", "a"))
  expect_equal(result$col_data$n, c(4L, 3L, 2L))
  expect_equal(as.vector(result$matrix[1, ]), c(100, 10, 1))
})

test_that("tally(sort = TRUE) permutes the matrix along with the counts", {
  mat <- matrix(c(rep(1, 2 * 2), rep(10, 2 * 3), rep(100, 2 * 4)), nrow = 2)
  tm <- tidymatrix(
    mat,
    data.frame(r = 1:2),
    data.frame(g = c(rep("a", 2), rep("b", 3), rep("c", 4)))
  )

  result <- tm |> activate(columns) |> group_by(g) |> tally(sort = TRUE)

  expect_equal(result$col_data$g, c("c", "b", "a"))
  expect_equal(as.vector(result$matrix[1, ]), c(100, 10, 1))
})

test_that("arrange() invalidates stored analyses", {
  mat <- matrix(rnorm(60), nrow = 6, ncol = 10)
  tm <- tidymatrix(mat, data.frame(gene = paste0("G", 1:6)))

  tm <- tm |> activate(rows) |> compute_hclust(k = 2)
  expect_equal(list_analyses(tm), "row_hclust")

  expect_warning(
    tm_sorted <- tm |> arrange(desc(gene)),
    "arrange"
  )
  expect_equal(list_analyses(tm_sorted), character(0))
})

test_that("summarize applies .matrix_fn to single-member groups", {
  tm <- tidymatrix(matrix(c(1, 2, 3, 4), 2, 2), data.frame(g = c("a", "b")))

  result <- tm |>
    activate(rows) |>
    group_by(g) |>
    summarise(.matrix_fn = function(v) 0)

  expect_equal(as.vector(result$matrix), c(0, 0, 0, 0))
})

test_that("summarize applies .matrix_fn to single-member column groups", {
  tm <- tidymatrix(
    matrix(c(1, 2, 3, 4), 2, 2),
    data.frame(r = 1:2),
    data.frame(g = c("a", "b"))
  )

  result <- tm |>
    activate(columns) |>
    group_by(g) |>
    summarise(.matrix_fn = function(v) 0)

  expect_equal(as.vector(result$matrix), c(0, 0, 0, 0))
})

test_that("compute_lm handles formulas too long for a single deparse line", {
  set.seed(2)
  mat <- matrix(rnorm(40), nrow = 2)
  col_data <- data.frame(
    a_very_long_variable_name_one = rnorm(20),
    a_very_long_variable_name_two = rnorm(20),
    a_very_long_variable_name_three = rnorm(20),
    a_very_long_variable_name_four = rnorm(20)
  )
  tm <- tidymatrix(mat, data.frame(id = 1:2), col_data)

  result <- tm |>
    activate(rows) |>
    compute_lm(
      ~ a_very_long_variable_name_one + a_very_long_variable_name_two +
        a_very_long_variable_name_three + a_very_long_variable_name_four
    )

  expect_true("r.squared" %in% names(result))
  expect_false("error" %in% names(result))
  expect_equal(nrow(result), 2)
})

test_that("compute_* with add_to_data errors clearly on colliding column names", {
  mat <- matrix(rnorm(600, mean = 10, sd = 2), nrow = 60, ncol = 10)
  tm <- tidymatrix(
    mat,
    data.frame(gene_id = paste0("Gene_", 1:60)),
    data.frame(
      condition = rep(c("Control", "Treatment"), each = 5),
      age = rnorm(10, 50, 10)
    )
  )

  once <- tm |>
    activate(rows) |>
    compute_ttest(group_col = "condition", add_to_data = TRUE)

  # Both tests produce a `p.value` column; the second must not silently
  # rename it out from under the p-value adjustment.
  expect_error(
    compute_correlation(once, var = "age", add_to_data = TRUE),
    "already present in row_data"
  )

  # Supplying a prefix resolves it.
  twice <- once |> compute_correlation(var = "age", add_to_data = TRUE, prefix = "age_cor")
  expect_true(all(c("p.value", "age_cor_p.value") %in% names(twice$row_data)))
})

test_that("plot_pheatmap lets the caller override derived annotations", {
  skip_if_not_installed("pheatmap")

  mat <- matrix(rnorm(50), nrow = 10, ncol = 5)
  tm <- tidymatrix(
    mat,
    data.frame(gene = paste0("Gene_", 1:10), type = rep(c("A", "B"), 5)),
    data.frame(sample = paste0("S", 1:5), cond = c("x", "x", "y", "y", "y"))
  )

  my_rows <- data.frame(custom = factor(rep(c("p", "q"), 5)))
  rownames(my_rows) <- paste0("Gene_", 1:10)

  # Previously this failed: "annotation_row matched by multiple actual arguments"
  p <- plot_pheatmap(
    tm,
    row_names = "gene",
    col_names = "sample",
    annotation_row = my_rows
  )
  expect_s3_class(p, "pheatmap")
})

# Joins keep matrix type and dimnames ----------------------------------------

test_that("row joins keep integer storage and row names", {
  mat <- matrix(1:6, nrow = 3, dimnames = list(c("r1", "r2", "r3"), c("a", "b")))
  tm <- tidymatrix(mat, data.frame(id = c("x", "y", "z")))
  ext <- data.frame(id = c("z", "x"), extra = 1:2)

  result <- tm |> activate(rows) |> inner_join(ext, by = "id")

  expect_type(result$matrix, "integer")
  expect_equal(rownames(result$matrix), c("r1", "r3"))
  expect_equal(colnames(result$matrix), c("a", "b"))
  expect_equal(unname(result$matrix[, "a"]), c(1L, 3L))
})

test_that("full joins add NA rows with NA row names", {
  mat <- matrix(1:4, nrow = 2, dimnames = list(c("r1", "r2"), c("a", "b")))
  tm <- tidymatrix(mat, data.frame(id = c("x", "y")))
  ext <- data.frame(id = c("y", "new"), extra = 1:2)

  result <- tm |> activate(rows) |> full_join(ext, by = "id")

  expect_type(result$matrix, "integer")
  expect_equal(nrow(result$matrix), 3)
  expect_equal(rownames(result$matrix)[1:2], c("r1", "r2"))
  expect_true(is.na(rownames(result$matrix)[3]))
  expect_true(all(is.na(result$matrix[3, ])))
})

test_that("column joins keep storage type and column names", {
  mat <- matrix(c(TRUE, FALSE, TRUE, TRUE), nrow = 2,
                dimnames = list(c("r1", "r2"), c("a", "b")))
  tm <- tidymatrix(mat, NULL, data.frame(sample = c("s1", "s2")))
  ext <- data.frame(sample = "s2", extra = 1)

  result <- tm |> activate(columns) |> semi_join(ext, by = "sample")

  expect_type(result$matrix, "logical")
  expect_equal(colnames(result$matrix), "b")
  expect_equal(rownames(result$matrix), c("r1", "r2"))
})

# summarize() invalidates stored analyses -----------------------------------

test_that("summarize removes stored analyses but keeps metadata columns", {
  set.seed(1)
  mat <- matrix(rnorm(40), nrow = 10)
  tm <- tidymatrix(mat, data.frame(g = rep(c("a", "b"), 5))) |>
    activate(rows) |>
    compute_prcomp(n_components = 1)

  expect_equal(list_analyses(tm), "row_pca")

  expect_warning(
    result <- tm |> activate(rows) |> group_by(g) |> summarise(pc = mean(row_pca_PC1)),
    "summarize"
  )
  expect_length(list_analyses(result), 0)
  expect_true("pc" %in% names(result$row_data))
})

test_that("count removes stored analyses", {
  mat <- matrix(rnorm(40), nrow = 10)
  tm <- tidymatrix(mat, data.frame(g = rep(c("a", "b"), 5))) |>
    activate(columns) |>
    compute_hclust()

  expect_warning(
    result <- tm |> activate(rows) |> count(g),
    "summarize"
  )
  expect_length(list_analyses(result), 0)
})

# Deterministic group order in two-group tests -------------------------------

make_group_tm <- function(groups) {
  set.seed(2)
  mat <- matrix(rnorm(3 * length(groups), mean = 10), nrow = 3)
  tidymatrix(mat, data.frame(gene = c("g1", "g2", "g3")),
             data.frame(group = groups))
}

test_that("compute_ttest infers groups independent of data order", {
  tm1 <- make_group_tm(rep(c("B", "A"), each = 4))
  tm2 <- make_group_tm(rep(c("A", "B"), each = 4))

  expect_message(
    r1 <- tm1 |> activate(rows) |> compute_ttest(group_col = "group"),
    "treatment = 'B' vs control = 'A'"
  )
  expect_message(
    r2 <- tm2 |> activate(rows) |> compute_ttest(group_col = "group"),
    "treatment = 'B' vs control = 'A'"
  )

  # Same direction: log2(mean B / mean A)
  m <- tm1$matrix
  expected <- log2(rowMeans(m[, 1:4]) / rowMeans(m[, 5:8]))
  expect_equal(r1$log2fc, unname(expected))
})

test_that("compute_ttest uses factor level order", {
  groups <- factor(rep(c("low", "high"), each = 4), levels = c("low", "high"))
  tm <- make_group_tm(groups)

  expect_message(
    tm |> activate(rows) |> compute_ttest(group_col = "group"),
    "treatment = 'high' vs control = 'low'"
  )
})

test_that("compute_ttest infers the other group when one is given", {
  tm <- make_group_tm(rep(c("A", "B"), each = 4))

  expect_message(
    tm |> activate(rows) |> compute_ttest(group_col = "group", treatment = "A"),
    "treatment = 'A' vs control = 'B'"
  )
})

test_that("explicit groups are silent and may pick two of several groups", {
  tm <- make_group_tm(rep(c("A", "B", "C"), each = 3))

  expect_silent(
    res <- tm |>
      activate(rows) |>
      compute_ttest(group_col = "group", control = "A", treatment = "C")
  )
  expect_equal(nrow(res), 3)
})

test_that("compute_ttest errors on unknown groups and columns", {
  tm <- make_group_tm(rep(c("A", "B"), each = 4))

  expect_error(
    tm |> activate(rows) |>
      compute_ttest(group_col = "group", control = "A", treatment = "X"),
    "not found"
  )
  expect_error(
    tm |> activate(rows) |> compute_ttest(group_col = "nope"),
    "must be the name of a column"
  )
})

test_that("missing group labels are ignored", {
  tm <- make_group_tm(c("A", "A", "A", NA, "B", "B", "B", NA))

  res <- tm |>
    activate(rows) |>
    compute_ttest(group_col = "group", control = "A", treatment = "B")
  expect_false(anyNA(res$p.value))

  res_w <- tm |>
    activate(rows) |>
    compute_wilcox(group_col = "group", control = "A", treatment = "B")
  expect_false(anyNA(res_w$median_diff))
})

# transform_matrix() arguments can use metadata ------------------------------

test_that("transform_matrix evaluates arguments in col_data when rows active", {
  mat <- matrix(1:6, nrow = 2)
  tm <- tidymatrix(mat, NULL, data.frame(flip = c(TRUE, FALSE, TRUE)))

  result <- tm |>
    activate(rows) |>
    transform_matrix(\(x, flip) ifelse(flip, 6L - x, x), flip = flip)

  expect_equal(result$matrix, matrix(c(5L, 4L, 3L, 4L, 1L, 0L), nrow = 2))
})

test_that("transform_matrix evaluates arguments in row_data when columns active", {
  mat <- matrix(c(1, 2, 3, 4), nrow = 2)
  tm <- tidymatrix(mat, data.frame(scale = c(10, 100)))

  result <- tm |>
    activate(columns) |>
    transform_matrix(\(x, s) x * s, s = scale)

  expect_equal(result$matrix, matrix(c(10, 200, 30, 400), nrow = 2))
})

test_that("transform_matrix can reach the calling environment with .env", {
  mat <- matrix(1:4, nrow = 2)
  tm <- tidymatrix(mat, NULL, data.frame(k = c(100, 200)))
  k <- 1

  result <- tm |>
    activate(rows) |>
    transform_matrix(\(x, add) x + add, add = .env$k)

  expect_equal(result$matrix, mat + 1)
})

test_that("transform_matrix evaluates arguments normally when matrix active", {
  mat <- matrix(1:4, nrow = 2)
  tm <- tidymatrix(mat)
  offset <- 10

  result <- tm |> activate(matrix) |> transform_matrix(`+`, offset)
  expect_equal(result$matrix, mat + 10)
})

test_that("summarize keeps the names of the dimension that is not aggregated", {
  mat <- matrix(1:12, nrow = 4, dimnames = list(paste0("r", 1:4), c("a", "b", "c")))
  tm <- tidymatrix(mat, data.frame(g = c("x", "x", "y", "y")),
                   data.frame(h = c("u", "u", "v")))

  by_rows <- tm |> activate(rows) |> group_by(g) |> summarise()
  expect_equal(colnames(by_rows$matrix), c("a", "b", "c"))

  by_cols <- tm |> activate(columns) |> group_by(h) |> summarise()
  expect_equal(rownames(by_cols$matrix), paste0("r", 1:4))
})
