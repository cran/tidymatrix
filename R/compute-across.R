#' Apply a function across matrix dimensions with metadata access
#'
#' Applies a user-defined function to each row (when rows are active) or each
#' column (when columns are active), providing both the values and the opposite
#' dimension's metadata. This enables complex statistical modeling where you
#' need access to all annotations.
#'
#' @param .data A tidymatrix object with rows or columns active (not matrix)
#' @param fn A function with signature `function(values, metadata, ...)` where:
#'   - `values`: Numeric vector of row/column values
#'   - `metadata`: Complete col_data (rows active) or row_data (columns active)
#'   - `...`: Additional arguments from `compute_across()`
#'   Must return a named list or named vector
#' @param add_to_data Logical. If TRUE, adds results to row_data/col_data and
#'   returns modified tidymatrix. If FALSE (default), returns data.frame
#' @param prefix Character. Optional prefix for result column names when
#'   `add_to_data = TRUE`
#' @param return_tibble Logical. If TRUE (default), returns tibble. If FALSE,
#'   returns data.frame. Only applies when `add_to_data = FALSE`
#' @param ... Additional arguments passed to `fn`
#'
#' @return If `add_to_data = FALSE`: data.frame/tibble with one row per
#'   matrix row/column, containing identifiers and computed statistics.
#'   If `add_to_data = TRUE`: modified tidymatrix with results added to metadata.
#'
#' @export
#'
#' @examples
#' # T-test example
#' mat <- matrix(rnorm(100, mean = 10), nrow = 10, ncol = 10)
#' col_data <- data.frame(
#'   sample = paste0("S", 1:10),
#'   condition = rep(c("Control", "Treatment"), each = 5)
#' )
#' row_data <- data.frame(gene = paste0("Gene", 1:10))
#' tm <- tidymatrix(mat, row_data, col_data)
#'
#' # Run t-test on each row
#' results <- tm |>
#'   activate(rows) |>
#'   compute_across(
#'     fn = function(vals, meta) {
#'       test <- t.test(vals ~ meta$condition)
#'       list(
#'         p.value = test$p.value,
#'         log2fc = log2(mean(vals[meta$condition == "Treatment"]) /
#'                       mean(vals[meta$condition == "Control"]))
#'       )
#'     }
#'   )
#'
#' # Linear model with multiple predictors
#' col_data2 <- data.frame(
#'   sample = paste0("S", 1:10),
#'   condition = rep(c("Control", "Treatment"), each = 5),
#'   batch = factor(rep(1:2, 5)),
#'   age = rnorm(10, 50, 10)
#' )
#' tm2 <- tidymatrix(mat, row_data, col_data2)
#'
#' lm_results <- tm2 |>
#'   activate(rows) |>
#'   compute_across(
#'     fn = function(vals, meta) {
#'       fit <- lm(vals ~ condition + batch + age, data = meta)
#'       summ <- summary(fit)
#'       coef_summ <- coef(summ)
#'
#'       list(
#'         condition_pval = coef_summ["conditionTreatment", "Pr(>|t|)"],
#'         condition_coef = coef_summ["conditionTreatment", "Estimate"],
#'         r.squared = summ$r.squared
#'       )
#'     }
#'   )
#'
#' # Add results to metadata
#' tm_with_stats <- tm |>
#'   activate(rows) |>
#'   compute_across(
#'     fn = function(vals, meta) {
#'       test <- t.test(vals ~ meta$condition)
#'       list(p.value = test$p.value)
#'     },
#'     add_to_data = TRUE,
#'     prefix = "ttest"
#'   )
compute_across <- function(.data,
                          fn,
                          add_to_data = FALSE,
                          prefix = NULL,
                          return_tibble = TRUE,
                          ...) {
  # Validation
  if (!is_tidymatrix(.data)) {
    stop("compute_across() can only be used on tidymatrix objects", call. = FALSE)
  }

  if (.data$active == "matrix") {
    stop(
      "Cannot use compute_across() when matrix is active.\n",
      "  Use activate(rows) or activate(columns) first.",
      call. = FALSE
    )
  }

  if (!is.function(fn)) {
    stop("fn must be a function", call. = FALSE)
  }

  # Extract components based on active dimension
  if (.data$active == "rows") {
    n_iter <- nrow(.data$matrix)
    metadata <- .data$col_data
    id_data <- .data$row_data
  } else {  # columns active
    n_iter <- ncol(.data$matrix)
    metadata <- .data$row_data
    id_data <- .data$col_data
  }

  # Apply function to each row/column
  results_list <- vector("list", n_iter)
  failed <- logical(n_iter)

  for (i in seq_len(n_iter)) {
    # Extract values
    values <- if (.data$active == "rows") {
      .data$matrix[i, ]
    } else {
      .data$matrix[, i]
    }

    # Apply user function with error handling
    result <- tryCatch(fn(values, metadata, ...), error = function(e) e)
    if (inherits(result, "error")) {
      failed[i] <- TRUE
      result <- list(error = conditionMessage(result))
    }

    # Ensure result is a list
    if (!is.list(result)) {
      result <- as.list(result)
    }

    results_list[[i]] <- result
  }

  # If every call failed there are no result columns to return, and the
  # wrappers (compute_ttest() etc.) would fail later with an obscure error
  dim_name <- if (.data$active == "rows") "row" else "column"
  if (n_iter > 0 && all(failed)) {
    stop(
      "fn failed for every ", dim_name, ". First error: ",
      results_list[[1]]$error,
      call. = FALSE
    )
  }
  if (any(failed)) {
    warning(
      "fn failed for ", sum(failed), " of ", n_iter, " ", dim_name, "s (",
      paste(utils::head(which(failed), 5), collapse = ", "),
      if (sum(failed) > 5) ", ..." else "",
      "); their results are NA. First error: ",
      results_list[[which(failed)[1]]]$error,
      call. = FALSE
    )
  }

  # Combine results into data.frame
  results_df <- dplyr::bind_rows(results_list)

  # Apply prefix if specified
  if (!is.null(prefix) && prefix != "") {
    names(results_df) <- paste0(prefix, "_", names(results_df))
  }

  # Return or add to metadata
  if (add_to_data) {
    target <- if (.data$active == "rows") .data$row_data else .data$col_data

    # bind_cols() would silently rename colliding columns (p.value...4 etc.),
    # which breaks any downstream lookup by name. Fail with something the
    # caller can act on instead.
    clashes <- intersect(names(target), names(results_df))
    if (length(clashes) > 0) {
      stop(
        "Result column(s) already present in ",
        if (.data$active == "rows") "row_data" else "col_data", ": ",
        paste(clashes, collapse = ", "), "\n",
        "  Use `prefix` to give this analysis its own column names, e.g.\n",
        "  prefix = \"cor\" produces cor_", clashes[1], ".",
        call. = FALSE
      )
    }

    # Add results to appropriate metadata
    if (.data$active == "rows") {
      .data$row_data <- dplyr::bind_cols(.data$row_data, results_df)
    } else {
      .data$col_data <- dplyr::bind_cols(.data$col_data, results_df)
    }
    return(.data)
  } else {
    # Combine with identifiers and return
    if (return_tibble) {
      output <- dplyr::bind_cols(tibble::as_tibble(id_data), results_df)
    } else {
      output <- dplyr::bind_cols(id_data, results_df)
    }
    return(output)
  }
}

#' Convenience wrapper for t-tests
#'
#' Performs t-tests comparing two groups for each row (or column) of a
#' tidymatrix. Automatically detects groups and applies multiple testing
#' correction.
#'
#' @param .data A tidymatrix object with rows or columns active
#' @param group_col Character. Name of column in metadata containing group labels
#' @param control Character. Label for the control group. If NULL, it is
#'   inferred: when `group_col` is a factor, the first level present in the
#'   data; otherwise the first value in sorted order. If only `treatment` is
#'   given, `control` is the other group.
#' @param treatment Character. Label for the treatment group. If NULL, the
#'   group that is not `control`. When either group is inferred, a message
#'   reports the comparison being made. Set both explicitly to silence it.
#' @param log2 Logical. If TRUE (default), computes log2 fold change. If FALSE,
#'   computes raw fold change
#' @param adjust Character. Method for p-value adjustment. Default "fdr".
#'   See `?p.adjust` for options. Use "none" for no adjustment
#' @param ... Additional arguments passed to `t.test()`
#'
#' @inheritParams compute_across
#' @return A data.frame/tibble with columns: identifiers, p.value, log2fc (or fc),
#'   and p.adj (if adjustment applied)
#'
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(100, mean = 10), nrow = 10, ncol = 10)
#' col_data <- data.frame(
#'   sample = paste0("S", 1:10),
#'   condition = rep(c("Control", "Treatment"), each = 5)
#' )
#' row_data <- data.frame(gene = paste0("Gene", 1:10))
#' tm <- tidymatrix(mat, row_data, col_data)
#'
#' # Simple t-test
#' results <- tm |>
#'   activate(rows) |>
#'   compute_ttest(group_col = "condition")
#'
#' # With specific group labels
#' results <- tm |>
#'   activate(rows) |>
#'   compute_ttest(
#'     group_col = "condition",
#'     control = "Control",
#'     treatment = "Treatment"
#'   )
compute_ttest <- function(.data,
                      group_col,
                      control = NULL,
                      treatment = NULL,
                      log2 = TRUE,
                      adjust = "fdr",
                      add_to_data = FALSE,
                      prefix = NULL,
                      return_tibble = TRUE,
                      ...) {
  # Get metadata for group detection
  if (.data$active == "rows") {
    meta <- .data$col_data
  } else {
    meta <- .data$row_data
  }

  groups <- resolve_two_groups(meta, group_col, control, treatment)
  control <- groups$control
  treatment <- groups$treatment

  # Define test function
  test_fn <- function(vals, metadata) {
    # which() drops rows with a missing group label instead of turning them
    # into NA values
    g1 <- vals[which(metadata[[group_col]] == control)]
    g2 <- vals[which(metadata[[group_col]] == treatment)]

    test <- stats::t.test(g1, g2, ...)

    fc_value <- if (log2) {
      log2(mean(g2) / mean(g1))
    } else {
      mean(g2) / mean(g1)
    }

    list(
      p.value = test$p.value,
      fc = fc_value
    )
  }

  # Apply test
  results <- compute_across(.data, test_fn,
                           add_to_data = add_to_data,
                           prefix = prefix,
                           return_tibble = return_tibble)

  # If adding to data, we need to handle the renaming and p-value adjustment differently
  if (add_to_data) {
    # Modify the tidymatrix metadata
    if (.data$active == "rows") {
      metadata <- results$row_data
    } else {
      metadata <- results$col_data
    }

    # Rename fold change column
    fc_col <- if (!is.null(prefix)) paste0(prefix, "_fc") else "fc"
    if (log2) {
      names(metadata)[names(metadata) == fc_col] <-
        if (!is.null(prefix)) paste0(prefix, "_log2fc") else "log2fc"
    }

    # Adjust p-values if requested
    if (!is.null(adjust) && adjust != "none") {
      pval_col <- if (!is.null(prefix)) paste0(prefix, "_p.value") else "p.value"
      padj_col <- if (!is.null(prefix)) paste0(prefix, "_p.adj") else "p.adj"
      metadata[[padj_col]] <- stats::p.adjust(metadata[[pval_col]], method = adjust)
    }

    # Update the metadata in the tidymatrix
    if (.data$active == "rows") {
      results$row_data <- metadata
    } else {
      results$col_data <- metadata
    }

    return(results)
  } else {
    # Rename fold change column
    if (log2) {
      names(results)[names(results) == "fc"] <- "log2fc"
    }

    # Adjust p-values if requested
    if (!is.null(adjust) && adjust != "none") {
      results$p.adj <- stats::p.adjust(results$p.value, method = adjust)
    }

    return(results)
  }
}

#' Convenience wrapper for linear models
#'
#' Fits linear models for each row (or column) of a tidymatrix. Can return
#' statistics for a specific coefficient or overall model fit statistics.
#'
#' @param .data A tidymatrix object with rows or columns active
#' @param formula_rhs Right-hand side of formula (e.g., `~ condition + batch + age`)
#' @param coef Character. Name of coefficient to extract. If NULL, returns overall
#'   model statistics
#' @param adjust Character. Method for p-value adjustment. Default "fdr".
#'   Use "none" for no adjustment
#' @param ... Additional arguments passed to `lm()`
#'
#' @inheritParams compute_across
#' @return A data.frame/tibble with model statistics. If `coef` specified:
#'   estimate, p.value, se, and p.adj. If `coef = NULL`: r.squared, f.statistic,
#'   p.value, and p.adj
#'
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' col_data <- data.frame(
#'   sample = paste0("S", 1:10),
#'   condition = factor(rep(c("Control", "Treatment"), each = 5)),
#'   batch = factor(rep(1:2, 5)),
#'   age = rnorm(10, 50, 10)
#' )
#' row_data <- data.frame(gene = paste0("Gene", 1:10))
#' tm <- tidymatrix(mat, row_data, col_data)
#'
#' # Extract specific coefficient
#' results <- tm |>
#'   activate(rows) |>
#'   compute_lm(
#'     formula_rhs = ~ condition + batch + age,
#'     coef = "conditionTreatment"
#'   )
#'
#' # Overall model statistics
#' results <- tm |>
#'   activate(rows) |>
#'   compute_lm(formula_rhs = ~ condition + batch + age)
compute_lm <- function(.data,
                   formula_rhs,
                   coef = NULL,
                   adjust = "fdr",
                   add_to_data = FALSE,
                   prefix = NULL,
                   return_tibble = TRUE,
                   ...) {
  # Define model function
  model_fn <- function(vals, metadata) {
    # deparse() splits long formulas across several strings, which would build
    # one formula per line; deparse1() always collapses to a single string.
    fit <- stats::lm(
      stats::as.formula(paste("vals", deparse1(formula_rhs))),
      data = metadata,
      ...
    )
    summ <- summary(fit)

    if (!is.null(coef)) {
      # Extract specific coefficient
      coef_summ <- stats::coef(summ)
      if (!coef %in% rownames(coef_summ)) {
        stop("Coefficient '", coef, "' not found in model", call. = FALSE)
      }

      list(
        estimate = stats::coef(fit)[coef],
        p.value = coef_summ[coef, "Pr(>|t|)"],
        se = coef_summ[coef, "Std. Error"]
      )
    } else {
      # Overall model stats
      fstat <- summ$fstatistic
      list(
        r.squared = summ$r.squared,
        f.statistic = fstat[1],
        p.value = stats::pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
      )
    }
  }

  # Apply model
  results <- compute_across(.data, model_fn,
                           add_to_data = add_to_data,
                           prefix = prefix,
                           return_tibble = return_tibble)

  # Handle p-value adjustment
  if (add_to_data) {
    # Get the metadata
    if (.data$active == "rows") {
      metadata <- results$row_data
    } else {
      metadata <- results$col_data
    }

    # Adjust p-values if present and requested
    pval_col <- if (!is.null(prefix)) paste0(prefix, "_p.value") else "p.value"
    if (!is.null(adjust) && adjust != "none" && pval_col %in% names(metadata)) {
      padj_col <- if (!is.null(prefix)) paste0(prefix, "_p.adj") else "p.adj"
      metadata[[padj_col]] <- stats::p.adjust(metadata[[pval_col]], method = adjust)
    }

    # Update the metadata in the tidymatrix
    if (.data$active == "rows") {
      results$row_data <- metadata
    } else {
      results$col_data <- metadata
    }

    return(results)
  } else {
    # Adjust p-values if present and requested
    if (!is.null(adjust) && adjust != "none" && "p.value" %in% names(results)) {
      results$p.adj <- stats::p.adjust(results$p.value, method = adjust)
    }

    return(results)
  }
}

#' Convenience wrapper for Wilcoxon rank-sum test
#'
#' Performs Wilcoxon rank-sum tests (Mann-Whitney U test) comparing two groups
#' for each row (or column) of a tidymatrix. This is a non-parametric alternative
#' to the t-test that doesn't assume normal distribution.
#'
#' @param .data A tidymatrix object with rows or columns active
#' @param group_col Character. Name of column in metadata containing group labels
#' @param control Character. Label for the control group. If NULL, it is
#'   inferred: when `group_col` is a factor, the first level present in the
#'   data; otherwise the first value in sorted order. If only `treatment` is
#'   given, `control` is the other group.
#' @param treatment Character. Label for the treatment group. If NULL, the
#'   group that is not `control`. When either group is inferred, a message
#'   reports the comparison being made. Set both explicitly to silence it.
#' @param log2 Logical. If TRUE (default), computes log2 fold change. If FALSE,
#'   computes raw fold change
#' @param adjust Character. Method for p-value adjustment. Default "fdr".
#'   See `?p.adjust` for options. Use "none" for no adjustment
#' @param ... Additional arguments passed to `wilcox.test()`
#'
#' @inheritParams compute_across
#' @return A data.frame/tibble with columns: identifiers, p.value, log2fc (or fc),
#'   median_diff, and p.adj (if adjustment applied)
#'
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(100, mean = 10), nrow = 10, ncol = 10)
#' col_data <- data.frame(
#'   sample = paste0("S", 1:10),
#'   condition = rep(c("Control", "Treatment"), each = 5)
#' )
#' row_data <- data.frame(gene = paste0("Gene", 1:10))
#' tm <- tidymatrix(mat, row_data, col_data)
#'
#' # Wilcoxon test
#' results <- tm |>
#'   activate(rows) |>
#'   compute_wilcox(group_col = "condition")
compute_wilcox <- function(.data,
                        group_col,
                        control = NULL,
                        treatment = NULL,
                        log2 = TRUE,
                        adjust = "fdr",
                        add_to_data = FALSE,
                        prefix = NULL,
                        return_tibble = TRUE,
                        ...) {
  # Get metadata for group detection
  if (.data$active == "rows") {
    meta <- .data$col_data
  } else {
    meta <- .data$row_data
  }

  groups <- resolve_two_groups(meta, group_col, control, treatment)
  control <- groups$control
  treatment <- groups$treatment

  # Define test function
  test_fn <- function(vals, metadata) {
    # which() drops rows with a missing group label instead of turning them
    # into NA values
    g1 <- vals[which(metadata[[group_col]] == control)]
    g2 <- vals[which(metadata[[group_col]] == treatment)]

    test <- suppressWarnings(stats::wilcox.test(g1, g2, ...))

    fc_value <- if (log2) {
      log2(mean(g2) / mean(g1))
    } else {
      mean(g2) / mean(g1)
    }

    list(
      p.value = test$p.value,
      fc = fc_value,
      median_diff = stats::median(g2) - stats::median(g1)
    )
  }

  # Apply test
  results <- compute_across(.data, test_fn,
                           add_to_data = add_to_data,
                           prefix = prefix,
                           return_tibble = return_tibble)

  # Handle renaming and p-value adjustment
  if (add_to_data) {
    # Get the metadata
    if (.data$active == "rows") {
      metadata <- results$row_data
    } else {
      metadata <- results$col_data
    }

    # Rename fold change column
    fc_col <- if (!is.null(prefix)) paste0(prefix, "_fc") else "fc"
    if (log2) {
      names(metadata)[names(metadata) == fc_col] <-
        if (!is.null(prefix)) paste0(prefix, "_log2fc") else "log2fc"
    }

    # Adjust p-values if requested
    if (!is.null(adjust) && adjust != "none") {
      pval_col <- if (!is.null(prefix)) paste0(prefix, "_p.value") else "p.value"
      padj_col <- if (!is.null(prefix)) paste0(prefix, "_p.adj") else "p.adj"
      metadata[[padj_col]] <- stats::p.adjust(metadata[[pval_col]], method = adjust)
    }

    # Update the metadata in the tidymatrix
    if (.data$active == "rows") {
      results$row_data <- metadata
    } else {
      results$col_data <- metadata
    }

    return(results)
  } else {
    # Rename fold change column
    if (log2) {
      names(results)[names(results) == "fc"] <- "log2fc"
    }

    # Adjust p-values if requested
    if (!is.null(adjust) && adjust != "none") {
      results$p.adj <- stats::p.adjust(results$p.value, method = adjust)
    }

    return(results)
  }
}

#' Convenience wrapper for ANOVA
#'
#' Performs one-way ANOVA for each row (or column) of a tidymatrix to test
#' for differences across multiple groups.
#'
#' @param .data A tidymatrix object with rows or columns active
#' @param group_col Character. Name of column in metadata containing group labels
#' @param adjust Character. Method for p-value adjustment. Default "fdr".
#'   Use "none" for no adjustment
#'
#' @inheritParams compute_across
#' @return A data.frame/tibble with columns: identifiers, f.statistic, p.value,
#'   df_between, df_within, and p.adj (if adjustment applied)
#'
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(150), nrow = 10, ncol = 15)
#' col_data <- data.frame(
#'   sample = paste0("S", 1:15),
#'   condition = rep(c("A", "B", "C"), each = 5)
#' )
#' row_data <- data.frame(gene = paste0("Gene", 1:10))
#' tm <- tidymatrix(mat, row_data, col_data)
#'
#' # ANOVA
#' results <- tm |>
#'   activate(rows) |>
#'   compute_anova(group_col = "condition")
compute_anova <- function(.data,
                       group_col,
                       adjust = "fdr",
                       add_to_data = FALSE,
                       prefix = NULL,
                       return_tibble = TRUE) {
  # Define ANOVA function
  anova_fn <- function(vals, metadata) {
    fit <- stats::aov(vals ~ metadata[[group_col]])
    summ <- summary(fit)[[1]]

    list(
      f.statistic = summ["metadata[[group_col]]", "F value"],
      p.value = summ["metadata[[group_col]]", "Pr(>F)"],
      df_between = summ["metadata[[group_col]]", "Df"],
      df_within = summ["Residuals", "Df"]
    )
  }

  # Apply ANOVA
  results <- compute_across(.data, anova_fn,
                           add_to_data = add_to_data,
                           prefix = prefix,
                           return_tibble = return_tibble)

  # Handle p-value adjustment
  if (add_to_data) {
    # Get the metadata
    if (.data$active == "rows") {
      metadata <- results$row_data
    } else {
      metadata <- results$col_data
    }

    # Adjust p-values if requested
    if (!is.null(adjust) && adjust != "none") {
      pval_col <- if (!is.null(prefix)) paste0(prefix, "_p.value") else "p.value"
      padj_col <- if (!is.null(prefix)) paste0(prefix, "_p.adj") else "p.adj"
      metadata[[padj_col]] <- stats::p.adjust(metadata[[pval_col]], method = adjust)
    }

    # Update the metadata in the tidymatrix
    if (.data$active == "rows") {
      results$row_data <- metadata
    } else {
      results$col_data <- metadata
    }

    return(results)
  } else {
    # Adjust p-values if requested
    if (!is.null(adjust) && adjust != "none") {
      results$p.adj <- stats::p.adjust(results$p.value, method = adjust)
    }

    return(results)
  }
}

#' Convenience wrapper for Kruskal-Wallis test
#'
#' Performs Kruskal-Wallis tests for each row (or column) of a tidymatrix.
#' This is a non-parametric alternative to one-way ANOVA.
#'
#' @param .data A tidymatrix object with rows or columns active
#' @param group_col Character. Name of column in metadata containing group labels
#' @param adjust Character. Method for p-value adjustment. Default "fdr".
#'   Use "none" for no adjustment
#'
#' @inheritParams compute_across
#' @return A data.frame/tibble with columns: identifiers, statistic, p.value,
#'   df, and p.adj (if adjustment applied)
#'
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(150), nrow = 10, ncol = 15)
#' col_data <- data.frame(
#'   sample = paste0("S", 1:15),
#'   condition = rep(c("A", "B", "C"), each = 5)
#' )
#' row_data <- data.frame(gene = paste0("Gene", 1:10))
#' tm <- tidymatrix(mat, row_data, col_data)
#'
#' # Kruskal-Wallis test
#' results <- tm |>
#'   activate(rows) |>
#'   compute_kruskal(group_col = "condition")
compute_kruskal <- function(.data,
                         group_col,
                         adjust = "fdr",
                         add_to_data = FALSE,
                         prefix = NULL,
                         return_tibble = TRUE) {
  # Define Kruskal-Wallis function
  kruskal_fn <- function(vals, metadata) {
    test <- stats::kruskal.test(vals ~ metadata[[group_col]])

    list(
      statistic = test$statistic,
      p.value = test$p.value,
      df = test$parameter
    )
  }

  # Apply test
  results <- compute_across(.data, kruskal_fn,
                           add_to_data = add_to_data,
                           prefix = prefix,
                           return_tibble = return_tibble)

  # Handle p-value adjustment
  if (add_to_data) {
    # Get the metadata
    if (.data$active == "rows") {
      metadata <- results$row_data
    } else {
      metadata <- results$col_data
    }

    # Adjust p-values if requested
    if (!is.null(adjust) && adjust != "none") {
      pval_col <- if (!is.null(prefix)) paste0(prefix, "_p.value") else "p.value"
      padj_col <- if (!is.null(prefix)) paste0(prefix, "_p.adj") else "p.adj"
      metadata[[padj_col]] <- stats::p.adjust(metadata[[pval_col]], method = adjust)
    }

    # Update the metadata in the tidymatrix
    if (.data$active == "rows") {
      results$row_data <- metadata
    } else {
      results$col_data <- metadata
    }

    return(results)
  } else {
    # Adjust p-values if requested
    if (!is.null(adjust) && adjust != "none") {
      results$p.adj <- stats::p.adjust(results$p.value, method = adjust)
    }

    return(results)
  }
}

#' Convenience wrapper for correlation tests
#'
#' Performs correlation tests between each row (or column) and a continuous
#' variable in the metadata.
#'
#' @param .data A tidymatrix object with rows or columns active
#' @param var Character. Name of continuous variable in metadata to correlate with
#' @param method Character. Correlation method: "pearson" (default), "spearman",
#'   or "kendall"
#' @param adjust Character. Method for p-value adjustment. Default "fdr".
#'   Use "none" for no adjustment
#'
#' @inheritParams compute_across
#' @return A data.frame/tibble with columns: identifiers, correlation, p.value,
#'   and p.adj (if adjustment applied)
#'
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' col_data <- data.frame(
#'   sample = paste0("S", 1:10),
#'   age = rnorm(10, 50, 10)
#' )
#' row_data <- data.frame(gene = paste0("Gene", 1:10))
#' tm <- tidymatrix(mat, row_data, col_data)
#'
#' # Pearson correlation
#' results <- tm |>
#'   activate(rows) |>
#'   compute_correlation(var = "age", method = "pearson")
#'
#' # Spearman correlation
#' results <- tm |>
#'   activate(rows) |>
#'   compute_correlation(var = "age", method = "spearman")
compute_correlation <- function(.data,
                             var,
                             method = "pearson",
                             adjust = "fdr",
                             add_to_data = FALSE,
                             prefix = NULL,
                             return_tibble = TRUE) {
  # Define correlation function
  cor_fn <- function(vals, metadata) {
    test <- stats::cor.test(vals, metadata[[var]], method = method)

    list(
      correlation = as.numeric(test$estimate),
      p.value = test$p.value
    )
  }

  # Apply correlation test
  results <- compute_across(.data, cor_fn,
                           add_to_data = add_to_data,
                           prefix = prefix,
                           return_tibble = return_tibble)

  # Handle p-value adjustment
  if (add_to_data) {
    # Get the metadata
    if (.data$active == "rows") {
      metadata <- results$row_data
    } else {
      metadata <- results$col_data
    }

    # Adjust p-values if requested
    if (!is.null(adjust) && adjust != "none") {
      pval_col <- if (!is.null(prefix)) paste0(prefix, "_p.value") else "p.value"
      padj_col <- if (!is.null(prefix)) paste0(prefix, "_p.adj") else "p.adj"
      metadata[[padj_col]] <- stats::p.adjust(metadata[[pval_col]], method = adjust)
    }

    # Update the metadata in the tidymatrix
    if (.data$active == "rows") {
      results$row_data <- metadata
    } else {
      results$col_data <- metadata
    }

    return(results)
  } else {
    # Adjust p-values if requested
    if (!is.null(adjust) && adjust != "none") {
      results$p.adj <- stats::p.adjust(results$p.value, method = adjust)
    }

    return(results)
  }
}

#' Convenience wrapper for simple linear regression
#'
#' Performs simple linear regression (one predictor) for each row (or column)
#' of a tidymatrix. For multiple predictors, use `compute_lm()` instead.
#'
#' @param .data A tidymatrix object with rows or columns active
#' @param predictor Character. Name of predictor variable in metadata
#' @param adjust Character. Method for p-value adjustment. Default "fdr".
#'   Use "none" for no adjustment
#'
#' @inheritParams compute_across
#' @return A data.frame/tibble with columns: identifiers, slope, intercept,
#'   r.squared, p.value, and p.adj (if adjustment applied)
#'
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' col_data <- data.frame(
#'   sample = paste0("S", 1:10),
#'   age = rnorm(10, 50, 10)
#' )
#' row_data <- data.frame(gene = paste0("Gene", 1:10))
#' tm <- tidymatrix(mat, row_data, col_data)
#'
#' # Simple linear regression
#' results <- tm |>
#'   activate(rows) |>
#'   compute_lm_simple(predictor = "age")
compute_lm_simple <- function(.data,
                           predictor,
                           adjust = "fdr",
                           add_to_data = FALSE,
                           prefix = NULL,
                           return_tibble = TRUE) {
  # Define regression function
  lm_fn <- function(vals, metadata) {
    fit <- stats::lm(vals ~ metadata[[predictor]])
    summ <- summary(fit)
    coef_summ <- stats::coef(summ)

    list(
      slope = stats::coef(fit)[2],
      intercept = stats::coef(fit)[1],
      r.squared = summ$r.squared,
      p.value = coef_summ[2, "Pr(>|t|)"],
      se = coef_summ[2, "Std. Error"]
    )
  }

  # Apply regression
  results <- compute_across(.data, lm_fn,
                           add_to_data = add_to_data,
                           prefix = prefix,
                           return_tibble = return_tibble)

  # Handle p-value adjustment
  if (add_to_data) {
    # Get the metadata
    if (.data$active == "rows") {
      metadata <- results$row_data
    } else {
      metadata <- results$col_data
    }

    # Adjust p-values if requested
    if (!is.null(adjust) && adjust != "none") {
      pval_col <- if (!is.null(prefix)) paste0(prefix, "_p.value") else "p.value"
      padj_col <- if (!is.null(prefix)) paste0(prefix, "_p.adj") else "p.adj"
      metadata[[padj_col]] <- stats::p.adjust(metadata[[pval_col]], method = adjust)
    }

    # Update the metadata in the tidymatrix
    if (.data$active == "rows") {
      results$row_data <- metadata
    } else {
      results$col_data <- metadata
    }

    return(results)
  } else {
    # Adjust p-values if requested
    if (!is.null(adjust) && adjust != "none") {
      results$p.adj <- stats::p.adjust(results$p.value, method = adjust)
    }

    return(results)
  }
}

#' Determine the control and treatment groups for two-group tests
#'
#' Groups that are not given explicitly are inferred in a deterministic way:
#' for a factor, in the order of its levels (unused levels dropped); otherwise
#' in sorted order. Missing values are ignored. This makes the direction of
#' fold changes independent of the row order of the metadata.
#'
#' @param meta Metadata data.frame containing `group_col`
#' @param group_col Name of the grouping column
#' @param control,treatment Group labels, or NULL to infer
#' @return A list with elements `control` and `treatment`
#' @noRd
resolve_two_groups <- function(meta, group_col, control = NULL,
                               treatment = NULL) {
  if (!is.character(group_col) || length(group_col) != 1 ||
      !group_col %in% names(meta)) {
    stop(
      "`group_col` must be the name of a column in the metadata.\n",
      "  Available columns: ", paste(names(meta), collapse = ", "),
      call. = FALSE
    )
  }

  values <- meta[[group_col]]
  if (is.factor(values)) {
    groups <- levels(droplevels(values[!is.na(values)]))
  } else {
    groups <- sort(unique(as.character(values[!is.na(values)])))
  }

  for (g in list(control, treatment)) {
    if (!is.null(g) && !as.character(g) %in% groups) {
      stop(
        "Group '", g, "' not found in ", group_col, ".\n",
        "  Available groups: ", paste(groups, collapse = ", "),
        call. = FALSE
      )
    }
  }

  if (!is.null(control) && !is.null(treatment)) {
    if (identical(as.character(control), as.character(treatment))) {
      stop("`control` and `treatment` must be different groups", call. = FALSE)
    }
    return(list(control = control, treatment = treatment))
  }

  if (length(groups) != 2) {
    stop(
      "Expected 2 groups in ", group_col, ", found ", length(groups),
      if (length(groups) > 0) paste0(" (", paste(groups, collapse = ", "), ")"),
      ".\n  Filter the data first, or set both `control` and `treatment`.",
      call. = FALSE
    )
  }

  if (is.null(control) && is.null(treatment)) {
    control <- groups[1]
    treatment <- groups[2]
  } else if (is.null(control)) {
    control <- setdiff(groups, as.character(treatment))
  } else {
    treatment <- setdiff(groups, as.character(control))
  }

  message(
    "Comparing ", group_col, ": treatment = '", treatment,
    "' vs control = '", control, "'"
  )

  list(control = control, treatment = treatment)
}
