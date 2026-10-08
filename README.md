
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tidymatrix

Many datasets are a matrix with two tables attached: survey answers with
information about the respondents and the questions, gene expression
with gene and sample annotations, sales per shop and product. Kept as
three separate objects, every filter or reordering has to be applied to
all three at once, and it is easy to end up with a matrix that no longer
lines up with its metadata.

`tidymatrix` keeps the three pieces together in one object and lets you
manipulate it with the dplyr verbs you already know. Much like
[tidygraph](https://tidygraph.data-imaginist.com/) does for networks,
you `activate()` the part you want to work on (rows, columns or the
matrix itself) and the rest follows along.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("raivokolde/tidymatrix")
```

## Example

The package ships with `big5`, a simulated personality survey: 400
people answered 30 statements on a 1–5 scale, and each statement
measures one of the Big Five personality traits.

``` r
library(tidymatrix)
library(dplyr, warn.conflicts = FALSE)

tm <- tidymatrix(big5_responses, big5_respondents, big5_items)
tm
#> # A tidymatrix: 400 x 30 matrix
#> # Active: matrix
#> #
#> # Row data: 400 rows x 8 columns
#> # Column data: 30 rows x 5 columns
#> #
#> # Matrix preview:
#>      E1 E2 E3 E4 E5 E6 A1 A2 A3 A4 A5 A6 C1 C2 C3 C4 C5 C6 N1 N2 N3 N4 N5 N6 O1
#> R001  4  3  5  3  3  4  3  5  3  4  4  3  1  2  1  2  4  5  3  3  3  2  3  2  2
#> R002  2  4  2  2  3  4  3  2  3  3  2  4  3  4  3  2  4  4  5  4  5  5  2  2  5
#> R003  2  2  3  4  4  3  3  3  4  5  4  3  2  3  1  3  4  4  4  4  3  2  4  2  1
#> R004  2  3  3  3  3  4  3  2  3  2  2  3  3  2  4  2  4  4  5  4  3  3  1  1  4
#> R005  1  2  2  1  5  5  5  5  5  5  2  1  2  5  4  3  2  2  4  2  3  4  3  2  2
#> R006  1  3  3  3  4  3  5  4  5  3  3  2  3  3  1  4  5  5  5  5  5  5  1  1  2
#>      O2 O3 O4 O5 O6
#> R001  4  1  3  2  3
#> R002  3  4  5  3  1
#> R003  4  1  2  5  4
#> R004  3  4  4  2  2
#> R005  3  2  3  4  3
#> R006  1  1  1  3  5
```

Filtering the metadata subsets the matrix to match. Here we drop
respondents who rushed through the survey and keep the neuroticism items
in the order they were asked:

``` r
tm_clean <- tm |>
  activate(rows) |>
  filter(completion_min > 3.5) |>
  activate(columns) |>
  filter(trait == "Neuroticism") |>
  arrange(position)

tm_clean
#> # A tidymatrix: 388 x 6 matrix
#> # Active: columns
#> #
#> # Row data: 388 rows x 8 columns
#> # Column data: 6 rows x 5 columns
#> #
#> # Active data (columns):
#>   item_id       trait reversed                   item_text position
#> 1      N1 Neuroticism    FALSE  I get stressed out easily.        4
#> 2      N2 Neuroticism    FALSE       I worry about things.        9
#> 3      N3 Neuroticism    FALSE      My mood changes often.       14
#> 4      N4 Neuroticism    FALSE     I get irritated easily.       19
#> 5      N5 Neuroticism     TRUE I stay calm under pressure.       24
#> 6      N6 Neuroticism     TRUE         I seldom feel blue.       29
```

Operations on the matrix can use the metadata too. With rows active,
`transform_matrix()` gets one respondent’s answers at a time, and extra
arguments can refer to item metadata, here to flip the reverse-keyed
items. `add_stats()` stores per-respondent summaries back into the row
metadata:

``` r
tm_clean |>
  activate(rows) |>
  transform_matrix(\(x, flip) ifelse(flip, 6L - x, x), flip = reversed) |>
  add_stats(mean) |>
  as_tibble() |>
  select(respondent_id, age, gender, neuroticism = mean)
#> # A tibble: 388 x 4
#>    respondent_id   age gender neuroticism
#>    <chr>         <int> <chr>        <dbl>
#>  1 R001             52 Female        3   
#>  2 R002             32 Female        4.5 
#>  3 R003             66 Male          3.17
#>  4 R004             42 Male          4.17
#>  5 R005             48 Female        3.33
#>  6 R006             62 Female        5   
#>  7 R007             26 Female        4.17
#>  8 R008             38 Female        5   
#>  9 R009             55 Male          3   
#> 10 R010             57 Female        2.33
#> # i 378 more rows
```

Analyses such as PCA and clustering store their results in the metadata
and keep the full result object for later:

``` r
tm |>
  activate(columns) |>
  compute_prcomp(n_components = 2) |>
  compute_hclust(k = 5, method = "ward.D2") |>
  list_analyses()
#> [1] "column_pca"    "column_hclust"
```

## What else is there

- **dplyr verbs**: `filter()`, `select()`, `mutate()`, `arrange()`,
  `slice()`, `group_by()` and `summarise()`, plus all six joins.
- **Matrix operations**: `transform_matrix()`, `scale()`, `center()`,
  `log_transform()`, `clip_values()`, `t()`, `add_stats()`.
- **Statistics for every row or column**: t-tests, Wilcoxon tests,
  ANOVA, Kruskal–Wallis, linear models, correlations, or any function of
  your own with `compute_across()`.
- **Dimensionality reduction and clustering**: PCA, MDS, t-SNE, UMAP,
  hierarchical clustering and k-means.
- **Getting data out**: `as_tibble()`, `to_long()` for ggplot2, and
  `plot_pheatmap()` for annotated heatmaps.

The [Get
started](https://raivokolde.github.io/tidymatrix/articles/basic-usage.html)
vignette is a quick tour of all of these, with links to a detailed
vignette for each topic.

## Related packages

- [tidygraph](https://tidygraph.data-imaginist.com/): tidy API for graph
  manipulation, the inspiration for this package
- [SummarizedExperiment](https://bioconductor.org/packages/SummarizedExperiment/):
  Bioconductor container for annotated matrix data
