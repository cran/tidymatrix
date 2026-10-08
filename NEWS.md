# tidymatrix 0.1.0

## Major Features

### Analysis Integration
* Added comprehensive support for analytical methods with automatic metadata integration
* PCA: `compute_prcomp()` performs PCA, adds scores to metadata and stores
  the full `prcomp` object (use `store = FALSE` to only add the scores)
* Clustering:
  - `compute_hclust()`: Hierarchical clustering with dendrogram storage
  - `compute_kmeans()`: K-means clustering
* Other embeddings: `compute_mds()`, `compute_tsne()` (Rtsne) and
  `compute_umap()` (umap)
* Analysis results are added as columns to metadata with clear naming:
  - PCA: `{name}_PC1`, `{name}_PC2`, etc.
  - Clustering: `{name}_cluster`

### Analysis Management
* `get_analysis()`: Retrieve stored analysis objects (prcomp, hclust, kmeans)
* `list_analyses()`: List all stored analyses
* `remove_analysis()`: Remove specific or all stored analyses
* `check_analyses()`: Check validity of stored analyses

### Analysis Invalidation
* Stored analysis objects are automatically removed when data is modified:
  `filter()`, `slice()`, `arrange()`, joins, `summarize()`/`count()`/`tally()`,
  `t()` and all matrix transformations
* Metadata columns (PC scores, cluster assignments) are preserved even when analysis objects are removed
* Clear warnings indicate which analyses were removed during data modification

## Core Features (from initial development)

### Data Structure
* `tidymatrix()`: Create tidymatrix objects combining matrix data with row and column metadata
* `activate()`: Switch context between rows, columns, or matrix
* Automatic validation of matrix-metadata alignment

### dplyr Integration
* `filter()`: Filter rows or columns with corresponding matrix subsetting
* `select()`: Select metadata columns
* `mutate()`: Add or modify metadata columns
* `arrange()`: Reorder rows or columns
* `slice()`, `slice_head()`, `slice_tail()`, `slice_sample()`: Subset by position
* `rename()`: Rename metadata columns
* `pull()`: Extract metadata columns as vectors
* `relocate()`: Reorder metadata columns

### Grouping and Aggregation
* `group_by()`: Group by metadata variables (creates grouped_tidymatrix)
* `summarize()`: Aggregate grouped data with matrix aggregation
  - Default `mean()` for numeric matrices
  - Required `.matrix_fn` parameter for non-numeric matrices
  - Type-aware error messages with helpful suggestions
* `count()`: Count observations by group with matrix aggregation
* `tally()`: Count within existing groups
* `ungroup()`: Remove grouping

## Design Principles

1. **Tidy principle**: Same data type in, same data type out
2. **Explicit naming**: Clear, prefixed column names prevent conflicts
3. **Smart defaults**: Sensible defaults with explicit requirements where needed
4. **Helpful errors**: Type-aware error messages guide users
5. **Analysis objects**: Store full R objects (prcomp, hclust) for later use
6. **Metadata persistence**: Metadata columns persist even when analysis objects are invalidated

## Row- and column-wise statistics
* `compute_across()` applies any function to every row or column, with access
  to the metadata of the other dimension
* Wrappers: `compute_ttest()`, `compute_wilcox()`, `compute_anova()`,
  `compute_kruskal()`, `compute_correlation()`, `compute_lm()`,
  `compute_lm_simple()`
* In two-group tests, groups that are not given explicitly are inferred in a
  deterministic order (factor levels, otherwise sorted values) and a message
  reports the comparison. Missing group labels are ignored.

## Matrix operations
* `transform_matrix()`, `scale()`, `center()`, `clip_values()`,
  `log_transform()`, `t()` and `add_stats()`
* Extra arguments to `transform_matrix()` are evaluated with the metadata of
  the other dimension as data mask, so values can be transformed depending on
  their metadata (e.g. reverse-scoring questionnaire items)

## Joins
* All six dplyr joins on rows or columns. The matrix keeps its storage type and
  row/column names; rows or columns that only exist in the joined table are
  filled with `NA`

## Data and documentation
* New example datasets: a simulated Big Five personality survey
  (`big5_responses`, `big5_respondents`, `big5_items`) and a country table
  (`big5_countries`)
* Vignettes rewritten around the survey data: a getting-started tour plus
  vignettes on dplyr verbs, joins, matrix operations, PCA and clustering, other
  embeddings, row-wise statistics, and plotting and exporting
