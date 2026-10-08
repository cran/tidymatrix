#' Matrix Operations for tidymatrix
#'
#' Functions for transforming and manipulating the matrix component of
#' tidymatrix objects.
#'
#' @name matrix-operations
NULL

#' Center the active dimension of an object
#'
#' Generic for centering. See \code{\link{center.tidymatrix}} for the
#' tidymatrix method.
#'
#' @param x An object to center
#' @param ... Arguments passed to methods
#'
#' @return An object of the same class as \code{x}, centered
#' @examples
#' tm <- tidymatrix(big5_responses, big5_respondents, big5_items)
#'
#' # Center each item (column) on its mean
#' tm |>
#'   activate(columns) |>
#'   center()
#' @export
center <- function(x, ...) {
  UseMethod("center")
}

# Note: scale() is already a base R generic, so we don't need to create it
# Note: transform() is also a base R generic
# Note: t() is a base R generic, so we register t.tidymatrix directly

#' Transpose a tidymatrix
#'
#' Transpose the matrix and swap row and column metadata. This operation
#' works regardless of which component is active. Uses base R's \code{t()}
#' generic, so it works correctly even when tidyverse is loaded.
#'
#' @param x A tidymatrix object
#'
#' @return A tidymatrix object with transposed matrix and swapped metadata
#' @export
#'
#' @examples
#' mat <- matrix(1:12, nrow = 3, ncol = 4)
#' row_data <- data.frame(gene = paste0("G", 1:3))
#' col_data <- data.frame(sample = paste0("S", 1:4))
#' tm <- tidymatrix(mat, row_data, col_data)
#'
#' # Transpose: genes × samples → samples × genes
#' tm_t <- t(tm)
#' dim(tm_t$matrix)  # Now 4 × 3
t.tidymatrix <- function(x) {
  # Transpose the matrix
  x$matrix <- t(x$matrix)

  # Swap row and column metadata
  temp <- x$row_data
  x$row_data <- x$col_data
  x$col_data <- temp

  # Update active if it was rows or columns
  if (x$active == "rows") {
    x$active <- "columns"
  } else if (x$active == "columns") {
    x$active <- "rows"
  }

  # Remove all stored analyses since structure changed
  x <- remove_all_analyses(x, "transpose")

  x
}

#' Scale rows or columns
#'
#' Perform z-score scaling (center and scale to unit variance) on the active
#' dimension. Requires rows or columns to be active.
#'
#' @param x A tidymatrix object with rows or columns active
#' @param center If TRUE (default), center to mean = 0
#' @param scale If TRUE (default), scale to sd = 1
#'
#' @return A tidymatrix object with scaled matrix
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(20, mean = 10, sd = 5), nrow = 4, ncol = 5)
#' tm <- tidymatrix(mat)
#'
#' # Scale rows (z-score per row)
#' tm_scaled <- tm |>
#'   activate(rows) |>
#'   scale()
#'
#' # Scale columns
#' tm_scaled <- tm |>
#'   activate(columns) |>
#'   scale()
#'
#' # Only center, don't scale
#' tm_centered <- tm |>
#'   activate(rows) |>
#'   scale(center = TRUE, scale = FALSE)
scale.tidymatrix <- function(x, center = TRUE, scale = TRUE) {
  if (x$active == "matrix") {
    stop(
      "Cannot scale when matrix is active.\n",
      "  Use activate(rows) or activate(columns) to specify which dimension to scale.",
      call. = FALSE
    )
  }

  if (x$active == "rows") {
    # Scale each row
    x$matrix <- t(scale(t(x$matrix), center = center, scale = scale))
  } else if (x$active == "columns") {
    # Scale each column
    x$matrix <- scale(x$matrix, center = center, scale = scale)
  }

  # Remove stored analyses since matrix values changed
  x <- remove_all_analyses(x, "scale")

  x
}

#' Center rows or columns
#'
#' Center the active dimension by subtracting the mean. Requires rows or
#' columns to be active.
#'
#' @param x A tidymatrix object with rows or columns active
#' @param ... Not used
#'
#' @return A tidymatrix object with centered matrix
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(20, mean = 10), nrow = 4, ncol = 5)
#' tm <- tidymatrix(mat)
#'
#' # Center rows (mean of each row = 0)
#' tm_centered <- tm |>
#'   activate(rows) |>
#'   center()
#'
#' # Center columns (mean of each column = 0)
#' tm_centered <- tm |>
#'   activate(columns) |>
#'   center()
center.tidymatrix <- function(x, ...) {
  if (x$active == "matrix") {
    stop(
      "Cannot center when matrix is active.\n",
      "  Use activate(rows) or activate(columns) to specify which dimension to center.",
      call. = FALSE
    )
  }

  if (x$active == "rows") {
    # Center each row
    row_means <- rowMeans(x$matrix, na.rm = TRUE)
    x$matrix <- x$matrix - row_means
  } else if (x$active == "columns") {
    # Center each column
    col_means <- colMeans(x$matrix, na.rm = TRUE)
    x$matrix <- sweep(x$matrix, 2, col_means, "-")
  }

  # Remove stored analyses since matrix values changed
  x <- remove_all_analyses(x, "center")

  x
}

#' Apply a function to the matrix
#'
#' Applies a function to the matrix with behavior determined by the active
#' component:
#' \itemize{
#'   \item \code{activate(matrix)}: \code{fn} is applied to the entire matrix
#'     at once. Use functions like \code{log}, \code{exp}, \code{sqrt}, or
#'     anonymous functions like \code{\\(x) x^2}.
#'   \item \code{activate(rows)}: \code{fn} is applied independently to each
#'     row vector and must return a vector of the same length.
#'   \item \code{activate(columns)}: \code{fn} is applied independently to each
#'     column vector and must return a vector of the same length.
#' }
#'
#' Additional arguments in \code{...} are passed on to \code{fn}. With rows
#' or columns active they are evaluated with the metadata of the \emph{other}
#' dimension as a data mask, in the same way as in \code{dplyr::mutate()}.
#' A row vector has one element per matrix column, so with rows active a
#' column of \code{col_data} lines up element by element with the vector that
#' \code{fn} receives (and likewise for columns and \code{row_data}). This
#' makes it possible to transform values depending on their metadata, e.g.
#' to reverse-score some questionnaire items (see examples). Use
#' \code{.env$x} to refer to a variable \code{x} in the calling
#' environment when a metadata column has the same name. With the matrix
#' active, arguments are evaluated normally.
#'
#' @param .data A tidymatrix object
#' @param fn A function to apply. When matrix is active, receives the full
#'   matrix. When rows or columns are active, receives one row or column
#'   vector at a time. Must return values with the same dimensions.
#' @param ... Additional arguments passed to \code{fn}. With rows active they
#'   may refer to columns of \code{col_data}; with columns active, to columns
#'   of \code{row_data}. Avoid argument names that partially match \code{fn}
#'   (such as \code{f}), as R would match them to \code{fn}.
#'
#' @return A tidymatrix object with the transformed matrix
#' @export
#'
#' @examples
#' mat <- matrix(1:12, nrow = 3, ncol = 4)
#' tm <- tidymatrix(mat)
#'
#' # Element-wise: apply log to entire matrix
#' tm |>
#'   activate(matrix) |>
#'   transform_matrix(log)
#'
#' # Row-wise: rank values within each row
#' tm |>
#'   activate(rows) |>
#'   transform_matrix(rank)
#'
#' # Column-wise: min-max normalize each column to [0, 1]
#' tm |>
#'   activate(columns) |>
#'   transform_matrix(\(x) (x - min(x)) / (max(x) - min(x)))
#'
#' # Passing extra arguments: round to 2 decimal places
#' tm |>
#'   activate(matrix) |>
#'   transform_matrix(round, digits = 2)
#'
#' # Arguments can use the metadata of the other dimension. Reverse-score
#' # the questionnaire items flagged in the column metadata (1 <-> 5):
#' big5 <- tidymatrix(big5_responses, big5_respondents, big5_items)
#' big5 |>
#'   activate(rows) |>
#'   transform_matrix(\(x, flip) ifelse(flip, 6L - x, x), flip = reversed)
transform_matrix <- function(.data, fn, ...) {
  if (!is_tidymatrix(.data)) {
    stop("transform_matrix() can only be used on tidymatrix objects", call. = FALSE)
  }

  if (!is.function(fn)) {
    stop("fn must be a function", call. = FALSE)
  }

  original_dim <- dim(.data$matrix)
  original_dimnames <- dimnames(.data$matrix)

  # Evaluate the extra arguments once, with the metadata of the other
  # dimension as data mask: its rows line up with the elements of each
  # row/column vector that fn receives.
  mask <- switch(.data$active,
    rows = .data$col_data,
    columns = .data$row_data,
    NULL
  )
  if (!is.null(mask)) {
    mask <- as.data.frame(dplyr::ungroup(mask))
  }
  args <- lapply(rlang::enquos(...), rlang::eval_tidy, data = mask)

  call_fn <- function(x) do.call(fn, c(list(x), args))

  if (.data$active == "matrix") {
    .data$matrix <- call_fn(.data$matrix)
  } else if (.data$active == "rows") {
    rows <- lapply(seq_len(nrow(.data$matrix)), function(i) call_fn(.data$matrix[i, ]))
    .data$matrix <- do.call(rbind, rows)
  } else if (.data$active == "columns") {
    cols <- lapply(seq_len(ncol(.data$matrix)), function(j) call_fn(.data$matrix[, j]))
    .data$matrix <- do.call(cbind, cols)
  }

  # Coerce back to matrix if fn stripped the class (e.g. returned a plain vector)
  if (!is.matrix(.data$matrix)) {
    .data$matrix <- matrix(.data$matrix, nrow = original_dim[1], ncol = original_dim[2])
  }

  if (!identical(dim(.data$matrix), original_dim)) {
    stop(
      "fn must return values with the same dimensions as the input.\n",
      "  Expected: ", original_dim[1], " x ", original_dim[2], "\n",
      "  Got: ", dim(.data$matrix)[1], " x ", dim(.data$matrix)[2],
      call. = FALSE
    )
  }

  dimnames(.data$matrix) <- original_dimnames

  # Remove stored analyses since matrix values changed
  .data <- remove_all_analyses(.data, "transform_matrix")

  .data
}

#' Add statistics to metadata
#'
#' Compute statistics for the active dimension and add them as columns to
#' the corresponding metadata. Requires rows or columns to be active.
#'
#' @param .data A tidymatrix object with rows or columns active
#' @param ... Functions to compute statistics (unquoted names like mean, var, sd)
#' @param .fns Alternative way to specify functions as a list
#' @param .names Names for the new columns. If NULL, uses function names.
#'
#' @return A tidymatrix object with statistics added to metadata
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(20), nrow = 4, ncol = 5)
#' tm <- tidymatrix(mat)
#'
#' # Add row statistics
#' tm <- tm |>
#'   activate(rows) |>
#'   add_stats(mean, var, sd)
#'
#' # Now row_data has columns: mean, var, sd
#'
#' # Add column statistics
#' tm <- tm |>
#'   activate(columns) |>
#'   add_stats(median, min, max)
#'
#' # Custom names
#' tm <- tm |>
#'   activate(rows) |>
#'   add_stats(mean, var, .names = c("row_mean", "row_var"))
add_stats <- function(.data, ..., .fns = NULL, .names = NULL) {
  if (!inherits(.data, "tidymatrix")) {
    stop(".data must be a tidymatrix object", call. = FALSE)
  }

  if (.data$active == "matrix") {
    stop(
      "Cannot add_stats when matrix is active.\n",
      "  Use activate(rows) or activate(columns) to specify which dimension.",
      call. = FALSE
    )
  }

  # Capture function arguments
  fns <- rlang::enquos(...)

  # Use .fns if provided
  if (!is.null(.fns)) {
    if (!is.list(.fns)) {
      stop(".fns must be a list of functions", call. = FALSE)
    }
    fn_list <- .fns
    fn_names <- names(.fns)
    if (is.null(fn_names)) {
      fn_names <- paste0("stat", seq_along(.fns))
    }
  } else {
    # Extract function names from enquos
    fn_names <- names(fns)
    if (is.null(fn_names) || any(fn_names == "")) {
      # Use deparsed expressions as names
      fn_names <- vapply(fns, function(x) rlang::as_label(x), character(1))
    }
    fn_list <- lapply(fns, rlang::eval_tidy)
  }

  # Override with custom names if provided
  if (!is.null(.names)) {
    if (length(.names) != length(fn_list)) {
      stop(".names must have same length as number of functions", call. = FALSE)
    }
    fn_names <- .names
  }

  # Compute statistics
  if (.data$active == "rows") {
    # Compute for each row
    for (i in seq_along(fn_list)) {
      fn <- fn_list[[i]]
      stat_name <- fn_names[i]
      .data$row_data[[stat_name]] <- apply(.data$matrix, 1, fn)
    }
  } else if (.data$active == "columns") {
    # Compute for each column
    for (i in seq_along(fn_list)) {
      fn <- fn_list[[i]]
      stat_name <- fn_names[i]
      .data$col_data[[stat_name]] <- apply(.data$matrix, 2, fn)
    }
  }

  .data
}

#' Log transform matrix
#'
#' Apply log transformation to matrix values. Requires matrix to be active.
#' This is a convenience wrapper around \code{\link{transform_matrix}()}.
#'
#' @param .data A tidymatrix object with matrix active
#' @param base Logarithm base (2, 10, or "natural" for ln)
#' @param offset Value to add before log transform (default 1, for log(x + 1))
#'
#' @return A tidymatrix object with log-transformed matrix
#' @export
#'
#' @examples
#' mat <- matrix(abs(rnorm(20)), nrow = 4, ncol = 5)
#' tm <- tidymatrix(mat)
#'
#' # Log2 transform with pseudocount
#' tm_log <- tm |>
#'   activate(matrix) |>
#'   log_transform(base = 2, offset = 1)
#'
#' # Natural log
#' tm_ln <- tm |>
#'   activate(matrix) |>
#'   log_transform(base = "natural")
log_transform <- function(.data, base = 2, offset = 1) {
  if (!inherits(.data, "tidymatrix")) {
    stop(".data must be a tidymatrix object", call. = FALSE)
  }

  if (.data$active != "matrix") {
    stop(
      "Cannot log_transform when matrix is not active.\n",
      "  Use activate(matrix) first.",
      call. = FALSE
    )
  }

  # Apply log transformation
  if (base == 2) {
    .data$matrix <- log2(.data$matrix + offset)
  } else if (base == 10) {
    .data$matrix <- log10(.data$matrix + offset)
  } else if (base == "natural") {
    .data$matrix <- log(.data$matrix + offset)
  } else {
    .data$matrix <- log(.data$matrix + offset, base = base)
  }

  # Remove stored analyses since matrix values changed
  .data <- remove_all_analyses(.data, "log_transform")

  .data
}

#' Clip matrix values
#'
#' Cap matrix values at specified minimum and/or maximum. Requires matrix
#' to be active.
#'
#' @param .data A tidymatrix object with matrix active
#' @param min Minimum value (values below are set to this)
#' @param max Maximum value (values above are set to this)
#'
#' @return A tidymatrix object with clipped matrix values
#' @export
#'
#' @examples
#' mat <- matrix(rnorm(20), nrow = 4, ncol = 5)
#' tm <- tidymatrix(mat)
#'
#' # Clip to [-2, 2] range
#' tm_clipped <- tm |>
#'   activate(matrix) |>
#'   clip_values(min = -2, max = 2)
#'
#' # Only set floor
#' tm_floor <- tm |>
#'   activate(matrix) |>
#'   clip_values(min = 0)
clip_values <- function(.data, min = NULL, max = NULL) {
  if (!inherits(.data, "tidymatrix")) {
    stop(".data must be a tidymatrix object", call. = FALSE)
  }

  if (.data$active != "matrix") {
    stop(
      "Cannot clip_values when matrix is not active.\n",
      "  Use activate(matrix) first.",
      call. = FALSE
    )
  }

  if (is.null(min) && is.null(max)) {
    warning("No min or max specified, returning unchanged", call. = FALSE)
    return(.data)
  }

  # Clip values
  if (!is.null(min)) {
    .data$matrix[.data$matrix < min] <- min
  }
  if (!is.null(max)) {
    .data$matrix[.data$matrix > max] <- max
  }

  # Remove stored analyses since matrix values changed
  .data <- remove_all_analyses(.data, "clip_values")

  .data
}
