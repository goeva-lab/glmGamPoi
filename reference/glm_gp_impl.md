# Internal Function to Fit a Gamma-Poisson GLM

Internal Function to Fit a Gamma-Poisson GLM

## Usage

``` r
glm_gp_impl(
  Y,
  model_matrix,
  offset = 0,
  size_factors = c("normed_sum", "deconvolution", "poscounts", "ratio"),
  overdispersion = TRUE,
  overdispersion_shrinkage = TRUE,
  ridge_penalty = 0,
  do_cox_reid_adjustment = TRUE,
  subsample = FALSE,
  perf_optim,
  verbose = FALSE
)
```

## Arguments

- Y:

  any matrix-like object (e.g.
  [`matrix()`](https://rdrr.io/r/base/matrix.html), `DelayedArray()`,
  `HDF5Matrix()`) with one column per sample and row per gene.

- model_matrix:

  a numeric matrix that specifies the experimental design. It can be
  produced using
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html).
  Default: `NULL`

- offset:

  Constant offset in the model in addition to `log(size_factors)`. It
  can either be a single number, a vector of length `ncol(data)` or a
  matrix with the same dimensions as `dim(data)`. Default: `0`.

- size_factors:

  in large scale experiments, each sample is typically of different size
  (for example different sequencing depths). A size factor is an
  internal mechanism of GLMs to correct for this effect.\
  `size_factors` is either a numeric vector with positive entries that
  has the same lengths as columns in the data that specifies the size
  factors that are used. Or it can be a string that species the method
  that is used to estimate the size factors (one of
  `c("normed_sum", "deconvolution", "poscounts", "ratio")`). Note that
  `"normed_sum"` and `"poscounts"` are fairly simple methods and can
  lead to suboptimal results. For the best performance on data with many
  zeros, I recommend to use `size_factors = "deconvolution"` which calls
  [`scran::calculateSumFactors()`](https://rdrr.io/pkg/scran/man/computeSumFactors.html).
  However, you need to separately install the `scran` package from
  Bioconductor for this method to work. For small datasets common for
  bulk RNA-seq experiments, I recommend to use `size_factors = "ratio"`,
  which uses the same procedure as DESeq2 and edgeR. Also note that
  `size_factors = 1` and `size_factors = FALSE` are equivalent. If only
  a single gene is given, no size factor is estimated (ie.
  `size_factors = 1`). Default: `"normed_sum"`.

- overdispersion:

  the simplest count model is the Poisson model. However, the Poisson
  model assumes that \\variance = mean\\. For many applications this is
  too rigid and the Gamma-Poisson allows a more flexible mean-variance
  relation (\\variance = mean + mean^2 \* overdispersion\\).\
  `overdispersion` can either be

  - a single boolean that indicates if an overdispersion is estimated
    for each gene.

  - a numeric vector of length `nrow(data)` fixing the overdispersion to
    those values.

  - the string `"global"` to indicate that one dispersion is fit across
    all genes.

  Note that `overdispersion = 0` and `overdispersion = FALSE` are
  equivalent and both reduce the Gamma-Poisson to the classical Poisson
  model. Default: `TRUE`.

- overdispersion_shrinkage:

  the overdispersion can be difficult to estimate with few replicates.
  To improve the overdispersion estimates, we can share information
  across genes and shrink each individual overdispersion estimate
  towards a global overdispersion estimate. Empirical studies show
  however that the overdispersion varies based on the mean expression
  level (lower expression level =\> higher dispersion). If
  `overdispersion_shrinkage = TRUE`, a median trend of dispersion and
  expression level is fit and used to estimate the variances of a quasi
  Gamma Poisson model (Lund et al. 2012). Default: `TRUE`.

- ridge_penalty:

  to avoid overfitting, we can penalize fits with large coefficient
  estimates. Instead of directly minimizing the deviance per gene
  (\\\sum dev(y_i, X_i b)\\), we will minimize \\\sum dev(y_i, X_i b) +
  N \* \sum (penalty_p \* b_p)^2\\.\
  `ridge_penalty` can be

  - a scalar in which case all parameters except the intercept are
    penalized.

  - a vector which has to have the same length as columns in the model
    matrix

  - a matrix with the same number of columns as columns in the model
    matrix. This gives maximum flexibility for expert users and allows
    for full Tikhonov regularization.

  Default: `ridge_penalty = 0`, which is internally replaced with a
  small positive number for numerical stability.

- do_cox_reid_adjustment:

  the classical maximum likelihood estimator of the `overdisperion` is
  biased towards small values. McCarthy *et al.* (2012) showed that it
  is preferable to optimize the Cox-Reid adjusted profile likelihood.\
  `do_cox_reid_adjustment` can be either be `TRUE` or `FALSE` to
  indicate if the adjustment is added during the optimization of the
  `overdispersion` parameter. Default: `TRUE`.

- subsample:

  the estimation of the overdispersion is the slowest step when fitting
  a Gamma-Poisson GLM. For datasets with many samples, the estimation
  can be considerably sped up without loosing much precision by fitting
  the overdispersion only on a random subset of the samples. Default:
  `FALSE` which means that the data is not subsampled. If set to `TRUE`,
  at most 1,000 samples are considered. Otherwise the parameter just
  specifies the number of samples that are considered for each gene to
  estimate the overdispersion.

- perf_optim:

  a value that indicates if attempts should be made to optimize
  performance / memory use while fitting the GLM. Can be provided as a
  boolean (turns on/off all optimizations) or as a list to select only
  specific optimizations. For more details on specific optimization flag
  names and parameters, see
  [set_perf_optim_opts](set_perf_optim_opts.md).

  Note: various optimizations here can lead to BREAKING changes in
  output shape and/or values, and were designed/intended for application
  on single-cell transcriptomics data (e.g. sparse matrix, with l.t.
  ~10% density)

  Default: `FALSE`, meaning no optimizations are applied.

- verbose:

  a boolean that indicates if information about the individual steps are
  printed while fitting the GLM. Default: `FALSE`.

## Value

a list with four elements

- `Beta` the coefficient matrix

- `overdispersion` the vector with the estimated overdispersions

- `Mu` a matrix with the corresponding means for each gene and sample

- `size_factors` a vector with the size factor for each sample

- `ridge_penalty` a vector with the ridge penalty

## See also

[`glm_gp()`](glm_gp.md) and
[`overdispersion_mle()`](overdispersion_mle.md)
