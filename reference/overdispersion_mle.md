# Estimate the Overdispersion for a Vector of Counts

Estimate the Overdispersion for a Vector of Counts

## Usage

``` r
overdispersion_mle(
  y,
  mu,
  model_matrix = NULL,
  do_cox_reid_adjustment = !is.null(model_matrix),
  global_estimate = FALSE,
  subsample = FALSE,
  max_iter = 200,
  verbose = FALSE,
  use_nr_overdisp_impl = FALSE,
  do_parallel = 0
)
```

## Arguments

- y:

  a numeric or integer vector or matrix with the counts for which the
  overdispersion is estimated

- mu:

  a numeric vector of either length 1 or `length(y)` or if `y` is a
  matrix, a matrix with the same dimensions. Contains the predicted
  value for that sample. If missing: `mean(y)` / `rowMeans(y)`

- model_matrix:

  a numeric matrix that specifies the experimental design. It can be
  produced using
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html).
  Default: `NULL`

- do_cox_reid_adjustment:

  the classical maximum likelihood estimator of the `overdisperion` is
  biased towards small values. McCarthy *et al.* (2012) showed that it
  is preferable to optimize the Cox-Reid adjusted profile likelihood.\
  `do_cox_reid_adjustment` can be either be `TRUE` or `FALSE` to
  indicate if the adjustment is added during the optimization of the
  `overdispersion` parameter. Default: `TRUE` if a model matrix is
  provided, otherwise `FALSE`

- global_estimate:

  flag to decide if a single overdispersion for a whole matrix is
  calculated instead of one estimate per row. This parameter has no
  affect if `y` is a vector. Default: `FALSE`

- subsample:

  the estimation of the overdispersion is the slowest step when fitting
  a Gamma-Poisson GLM. For datasets with many samples, the estimation
  can be considerably sped up without loosing much precision by fitting
  the overdispersion only on a random subset of the samples. Default:
  `FALSE` which means that the data is not subsampled. If set to `TRUE`,
  at most 1,000 samples are considered. Otherwise the parameter just
  specifies the number of samples that are considered for each gene to
  estimate the overdispersion.

- max_iter:

  the maximum number of iterations for each gene

- verbose:

  a boolean that indicates if information about the individual steps are
  printed while fitting the GLM. Default: `FALSE`.

- use_nr_overdisp_impl:

  boolean flag to indicate if the alternate C++-based overdispersion
  estimation routine should be used

- do_parallel:

  integer to indicate degree of parallelization to be enabled, 0
  indicates no parallelism

## Value

The function returns a list with the following elements:

- `estimate`:

  the numerical estimate of the overdispersion.

- `iterations`:

  the number of iterations it took to calculate the result.

- `message`:

  additional information about the fitting process.

## Details

The function is optimized to be fast on many small counts. To achieve
this, the frequency table of the counts is calculated and used to avoid
repetitive calculations. If there are probably many unique counts the
optimization is skipped. Currently the heuristic is to skip if more than
half of the counts are expected to be unique. The estimation is based on
the largest observed count in `y`.

An earlier version of this package (\< 1.1.1) used a separate set of
functions for the case of many small counts based on a paper by Bandara
et al. (2019). However, this didn't bring a sufficient performance
increase and meant an additional maintenance burden.

## See also

[`glm_gp()`](glm_gp.md)

## Examples

``` r
 set.seed(1)
 # true overdispersion = 2.4
 y <- rnbinom(n = 10, mu = 3, size = 1/2.4)
 # estimate = 1.7
 overdispersion_mle(y)
#> $estimate
#> [1] 1.464695
#> 
#> $iterations
#> [1] 3
#> 
#> $message
#> [1] "both X-convergence and relative convergence (5)"
#> 


 # true overdispersion = 0
 y <- rpois(n = 10, lambda = 3)
 # estimate = 0
 overdispersion_mle(y)
#> $estimate
#> [1] 0
#> 
#> $iterations
#> [1] 0
#> 
#> $message
#> [1] "Even for very small theta, no maximum identified"
#> 
 # with different mu, overdispersion estimate changes
 overdispersion_mle(y, mu = 15)
#> $estimate
#> [1] 1.619791
#> 
#> $iterations
#> [1] 4
#> 
#> $message
#> [1] "relative convergence (4)"
#> 
 # Cox-Reid adjustment changes the result
 overdispersion_mle(y, mu = 15, do_cox_reid_adjustment = FALSE)
#> $estimate
#> [1] 1.619791
#> 
#> $iterations
#> [1] 4
#> 
#> $message
#> [1] "relative convergence (4)"
#> 


 # Many very small counts, true overdispersion = 50
 y <- rnbinom(n = 1000, mu = 0.01, size = 1/50)
 summary(y)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   0.000   0.000   0.000   0.013   0.000   3.000 
 # estimate = 31
 overdispersion_mle(y, do_cox_reid_adjustment = TRUE)
#> $estimate
#> [1] 31.37144
#> 
#> $iterations
#> [1] 3
#> 
#> $message
#> [1] "relative convergence (4)"
#> 

 # Function can also handle matrix input
 Y <- matrix(rnbinom(n = 10 * 3, mu = 4, size = 1/2.2), nrow = 10, ncol = 3)
 Y
#>       [,1] [,2] [,3]
#>  [1,]    0    5    2
#>  [2,]    0    0    2
#>  [3,]    0    3    0
#>  [4,]    0    7    0
#>  [5,]    4    0   12
#>  [6,]    7    2    1
#>  [7,]    0    0    6
#>  [8,]    0    3   28
#>  [9,]    0    6   15
#> [10,]    3    3    0
 as.data.frame(overdispersion_mle(Y))
#>      estimate iterations                                         message
#> 1   1.5655499          4                        relative convergence (4)
#> 2   7.9216048          6 both X-convergence and relative convergence (5)
#> 3  11.1896472          6 both X-convergence and relative convergence (5)
#> 4  16.8756932          6 both X-convergence and relative convergence (5)
#> 5   2.7529728          4 both X-convergence and relative convergence (5)
#> 6   0.5998486          3                        relative convergence (4)
#> 7  15.9185152          6 both X-convergence and relative convergence (5)
#> 8   4.4064051          4                        relative convergence (4)
#> 9   2.9537844          6 both X-convergence and relative convergence (5)
#> 10  0.8310498          5                        relative convergence (4)
```
