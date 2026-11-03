# Create a [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html) containing pseudo-bulk samples

Create a
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
containing pseudo-bulk samples

## Usage

``` r
pseudobulk(
  data,
  group_by,
  ...,
  aggregation_functions = list(counts = "rowSums2", .default = "rowMeans2"),
  col_data = NULL,
  make_colnames = TRUE,
  verbose = TRUE
)
```

## Arguments

- data:

  a
  [`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
  or an object of a related class

- group_by:

  an unquoted expression that can refer to columns in the
  [`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).
  All observations with the same factor level are aggregated. The
  argument follows the same logic as
  [`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html).
  The argument must wrapped using [`vars()`](vars.md).

- ...:

  named expressions that summarize columns in
  [`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).
  Each expression must produce a value of length 1. The arguments follow
  the same logic as
  [`dplyr::summarize()`](https://dplyr.tidyverse.org/reference/summarise.html).

- aggregation_functions:

  a named list with functions that are used to aggregate the assays in
  the `data`.

- col_data:

  additional data with `ncol(data)` rows. The `group_by` and named
  arguments can refer to the columns of the `col_data` in addition to
  the columns in `colData(data)` (assuming `data` is a
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)).

- make_colnames:

  a boolean that decides if the column names are the concatenated values
  of `group_by`. Default: `TRUE`

- verbose:

  a boolean that indicates if information about the process are printed
  Default: `TRUE`.

## Value

a
[`SingleCellExperiment`](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
object

## Examples

``` r
 suppressPackageStartupMessages(library(SingleCellExperiment))

 data <- data.frame(sample = sample(c("samp1", "samp2", "samp3"), size = 50, replace = TRUE),
       celltype = sample(c("T cells", "B cells", "Macrophages"), size = 50, replace = TRUE),
       size = rnorm(n = 50, mean = 40, sd = 15))
 Y <- matrix(rnbinom(n = 100 * 50, mu = 3, size = 1/3.1), nrow = 100, ncol = 50)
 sce <- SingleCellExperiment(Y, colData = data)
 aggr_sce <- pseudobulk(sce, group_by = vars(sample, celltype), size = mean(size))
#> Aggregating assay '1' using 'rowMeans2'.

 aggr_sce
#> class: SingleCellExperiment 
#> dim: 100 9 
#> metadata(0):
#> assays(1): ''
#> rownames: NULL
#> rowData names(0):
#> colnames(9): samp3.T cells samp2.Macrophages ... samp3.Macrophages
#>   samp1.Macrophages
#> colData names(3): sample celltype size
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):

 colData(aggr_sce)
#> DataFrame with 9 rows and 3 columns
#>                        sample    celltype      size
#>                   <character> <character> <numeric>
#> samp3.T cells           samp3     T cells   35.6695
#> samp2.Macrophages       samp2 Macrophages   39.6905
#> samp3.B cells           samp3     B cells   32.4984
#> samp2.T cells           samp2     T cells   30.1540
#> samp1.T cells           samp1     T cells   33.9702
#> samp1.B cells           samp1     B cells   38.4399
#> samp2.B cells           samp2     B cells   42.2837
#> samp3.Macrophages       samp3 Macrophages   55.6694
#> samp1.Macrophages       samp1 Macrophages   20.7775
```
