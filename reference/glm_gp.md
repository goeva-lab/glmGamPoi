# Fit a Gamma-Poisson Generalized Linear Model

This function provides a simple to use interface to fit Gamma-Poisson
generalized linear models. It works equally well for small scale (a
single model) and large scale data (e.g. thousands of rows and columns,
potentially stored on disk). The function automatically determines the
appropriate size factors for each sample and efficiently finds the best
overdispersion parameter for each gene.

## Usage

``` r
glm_gp(
  data,
  design = ~1,
  col_data = NULL,
  reference_level = NULL,
  offset = 0,
  size_factors = c("normed_sum", "deconvolution", "poscounts", "ratio"),
  overdispersion = TRUE,
  overdispersion_shrinkage = TRUE,
  ridge_penalty = 0,
  do_cox_reid_adjustment = TRUE,
  subsample = FALSE,
  on_disk = NULL,
  use_assay = NULL,
  perf_optim = FALSE,
  verbose = FALSE
)
```

## Arguments

- data:

  any matrix-like object (e.g.
  [matrix](https://rdrr.io/r/base/matrix.html), sparse matrix
  ([dgCMatrix](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)),
  [DelayedArray](https://rdrr.io/pkg/DelayedArray/man/DelayedArray-class.html),
  [HDF5Matrix](https://rdrr.io/pkg/HDF5Array/man/HDF5Array-class.html))
  or anything that can be cast to a
  [SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  (e.g. `MSnSet`, `eSet` etc.) with one column per sample and row per
  gene.

- design:

  a specification of the experimental design used to fit the
  Gamma-Poisson GLM. It can be a
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) with one
  row for each sample and one column for each coefficient.\
  Alternatively, `design` can be a `formula`. The entries in the formula
  can refer to global objects, columns in the `col_data` parameter, or
  the `colData(data)` of `data` if it is a `SummarizedExperiment`.\
  The third option is that `design` is a vector where each element
  specifies to which condition a sample belongs.\
  Default: `design = ~ 1`, which means that all samples are treated as
  if they belong to the same condition. Note that this is the fasted
  option.

- col_data:

  a dataframe with one row for each sample in `data`. Default: `NULL`.

- reference_level:

  a single string that specifies which level is used as reference when
  the model matrix is created. The reference level becomes the intercept
  and all other coefficients are calculated with respect to the
  `reference_level`. Default: `NULL`.

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

- on_disk:

  a boolean that indicates if the dataset is loaded into memory or if it
  is kept on disk to reduce the memory usage. Processing in memory can
  be significantly faster than on disk. Default: `NULL` which means that
  the data is only processed in memory if `data` is an in-memory data
  structure.

- use_assay:

  Specify which assay to use. Default: `NULL`, which means that if
  available 'counts' are used. Otherwise an error is thrown except if
  there is only a single assay.

- perf_optim:

  a value that indicates if attempts should be made to optimize
  performance / memory use while fitting the GLM. Can be provided as a
  boolean (turns on/off all optimizations) or as a list to select only
  specific optimizations. Note: various optimizations here can lead to
  BREAKING changes in output shape and/or values, and were
  designed/intended for application on single-cell transcriptomics data
  (e.g. sparse matrix, with l.t. ~10% density) A list of currently
  specified optimizations is provided below:

  - `offset_as_vec`: if `offset` is provided as either a constant or a
    vector of per-sample offsets, `Offset` will be a vector instead of a
    matrix to avoid generating a dense n_sample \* n_genes matrix

  - `cast_dgC_Y_to_dgR`: data matrix will be internally converted to
    row-major format to improve memory access patterns as internal steps
    fetch data row-wise on a per-gene basis

  - `do_parallel`: if value greater than zero, then enables
    parallelization for model parameter fitting processes the number of
    threads used is equal to the provided integer

  - `delay_mu`: Mu matrix is never explicitly formed, instead Mu field
    of output is provided as a function which can be called to
    materialize full matrix. IMPORTANT NOTE:
    [predict.glmGamPoi2](predict.glmGamPoi2.md) and
    [residuals.glmGamPoi2](residuals.glmGamPoi2.md) methods will lead to
    a memory spike if this option is enabled, as they require the full
    Mu matrix to be present in memory.

  - `use_nr_overdisp_impl`: uses a C++-implemented newton-raphson
    procedure instead of nlmimb for estimating overdispersions, which
    allows for parallelizing that step

  Default: `FALSE`, meaning no optimizations are made

- verbose:

  a boolean that indicates if information about the individual steps are
  printed while fitting the GLM. Default: `FALSE`.

## Value

The method returns a list with the following elements:

- `Beta`:

  a matrix with dimensions `nrow(data) x n_coefficients` where
  `n_coefficients` is based on the `design` argument. It contains the
  estimated coefficients for each gene.

- `overdispersions`:

  a vector with length `nrow(data)`. The overdispersion parameter for
  each gene. It describes how much more the counts vary than one would
  expect according to the Poisson model.

- `overdispersion_shrinkage_list`:

  a list with additional information from the quasi-likelihood
  shrinkage. For details see
  [`overdispersion_shrinkage()`](overdispersion_shrinkage.md).

- `deviances`:

  a vector with the deviance of the fit for each row. The deviance is a
  measure how well the data is fit by the model. A smaller deviance
  means a better fit.

- `Mu`:

  a matrix with the same dimensions as `dim(data)`. If the calculation
  happened on disk, than `Mu` is a `HDF5Matrix`. It contains the
  estimated mean value for each gene and sample.

- `size_factors`:

  a vector with length `ncol(data)`. The size factors are the inferred
  correction factors for different sizes of each sample. They are also
  sometimes called the exposure factor.

- `Offset`:

  a matrix with the same dimensions as `dim(data)`. If the calculation
  happened on disk, than `Offset` is a `HDF5Matrix`. It contains the
  `log(size_factors) + offset` from the function call.

- `data`:

  a `SummarizedExperiment` that contains the input counts and the
  `col_data`

- `model_matrix`:

  a matrix with dimensions `ncol(data) x n_coefficients`. It is build
  based on the `design` argument.

- `design_formula`:

  the formula that used to fit the model, or `NULL` otherwise

- `ridge_penalty`:

  a vector with the specification of the ridge penalty.

## Details

The method follows the following steps:

1.  The size factors are estimated.\
    If `size_factors = "normed_sum"`, the column-sum for each cell is
    calculated and the resulting size factors are normalized so that
    their geometric mean is `1`. If `size_factors = "poscounts"`, a
    slightly adapted version of the procedure proposed by Anders and
    Huber (2010) in equation (5) is used. To handle the large number of
    zeros the geometric means are calculated for \\Y + 0.5\\ and ignored
    during the calculation of the median. Columns with all zeros get a
    default size factor of \\0.001\\. If
    `size_factors = "deconvolution"`, the method
    [`scran::calculateSumFactors()`](https://rdrr.io/pkg/scran/man/computeSumFactors.html)
    is called. If `size_factors = "ratio"`, the unmodified procedure
    from Anders and Huber (2010) in equation (5) is used.

2.  The dispersion estimates are initialized based on the moments of
    each row of \\Y\\.

3.  The coefficients of the model are estimated.\
    If all samples belong to the same condition (i.e. `design = ~ 1`),
    the betas are estimated using a quick Newton-Raphson algorithm. This
    is similar to the behavior of `edgeR`. For more complex designs, the
    general Fisher-scoring algorithm is used. Here, the code is based on
    a fork of the internal function `fitBeta()` from `DESeq2`. It does
    however contain some modification to make the fit more robust and
    faster.

4.  The mean for each gene and sample is calculate.\
    Note that this step can be very IO intensive if `data` is or
    contains a DelayedArray.

5.  The overdispersion is estimated.\
    The classical method for estimating the overdispersion for each gene
    is to maximize the Gamma-Poisson log-likelihood by iterating over
    each count and summing the the corresponding log-likelihood. It is
    however, much more efficient for genes with many small counts to
    work on the contingency table of the counts. Originally, this
    approach had already been used by Anscombe (1950). In this package,
    I have implemented an extension of their method that can handle
    general offsets.\
    See also [`overdispersion_mle()`](overdispersion_mle.md).

6.  The beta coefficients are estimated once more with the updated
    overdispersion estimates

7.  The mean for each gene and sample is calculated again.

This method can handle not just in memory data, but also data stored on
disk. This is essential for large scale datasets with thousands of
samples, as they sometimes encountered in modern single-cell RNA-seq
analysis. `glmGamPoi` relies on the `DelayedArray` and `beachmat`
package to efficiently implement the access to the on-disk data.

## References

- McCarthy, D. J., Chen, Y., & Smyth, G. K. (2012). Differential
  expression analysis of multifactor RNA-Seq experiments with respect to
  biological variation. Nucleic Acids Research, 40(10), 4288–4297.
  <https://doi.org/10.1093/nar/gks042>.

- Anders Simon, & Huber Wolfgang. (2010). Differential expression
  analysis for sequence count data. Genome Biology.
  <https://doi.org/10.1016/j.jcf.2018.05.006>.

- Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of
  fold change and dispersion for RNA-seq data with DESeq2. Genome
  Biology, 15(12), 550. <https://doi.org/10.1186/s13059-014-0550-8>.

- Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2009). edgeR: A
  Bioconductor package for differential expression analysis of digital
  gene expression data. Bioinformatics, 26(1), 139–140.
  <https://doi.org/10.1093/bioinformatics/btp616>.

- Lun ATL, Pagès H, Smith ML (2018). “beachmat: A Bioconductor C++ API
  for accessing high-throughput biological data from a variety of R
  matrix types.” PLoS Comput. Biol., 14(5), e1006135. doi:
  [10.1371/journal.pcbi.1006135.](https://doi.org/10.1371/journal.pcbi.1006135).

- Lund, S. P., Nettleton, D., McCarthy, D. J., & Smyth, G. K. (2012).
  Detecting differential expression in RNA-sequence data using
  quasi-likelihood with shrunken dispersion estimates. Statistical
  Applications in Genetics and Molecular Biology, 11(5).
  <https://doi.org/10.1515/1544-6115.1826>.

- Lun ATL, Bach K and Marioni JC (2016). Pooling across cells to
  normalize single-cell RNA sequencing data with many zero counts.
  Genome Biol. 17:75 <https://doi.org/10.1186/s13059-016-0947-7>

## See also

[`overdispersion_mle()`](overdispersion_mle.md) and
[`overdispersion_shrinkage()`](overdispersion_shrinkage.md) for the
internal functions that do the work. For differential expression
analysis, see [`test_de()`](test_de.md).

## Examples

``` r
 set.seed(1)
 # The simplest example
 y <- rnbinom(n = 10, mu = 3, size = 1/2.4)
 c(glm_gp(y, size_factors = FALSE))
#> $Beta
#>      Intercept
#> [1,]  1.360977
#> 
#> $overdispersions
#> [1] 1.709539
#> 
#> $overdispersion_shrinkage_list
#> $overdispersion_shrinkage_list$dispersion_trend
#> [1] 1.709539
#> 
#> $overdispersion_shrinkage_list$ql_disp_estimate
#> [1] 1
#> 
#> $overdispersion_shrinkage_list$ql_disp_trend
#> [1] 1
#> 
#> $overdispersion_shrinkage_list$ql_disp_shrunken
#> [1] 1
#> 
#> $overdispersion_shrinkage_list$ql_df0
#> [1] Inf
#> 
#> 
#> $deviances
#> [1] 10.0644
#> 
#> $Mu
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]  3.9  3.9  3.9  3.9  3.9  3.9  3.9  3.9  3.9   3.9
#> 
#> $size_factors
#>  [1] 1 1 1 1 1 1 1 1 1 1
#> 
#> $Offset
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]    0    0    0    0    0    0    0    0    0     0
#> 
#> $ridge_penalty
#> Intercept 
#>     1e-11 
#> attr(,"target")
#> [1] 0
#> 
#> $data
#> class: SummarizedExperiment 
#> dim: 1 10 
#> metadata(0):
#> assays(1): counts
#> rownames: NULL
#> rowData names(0):
#> colnames: NULL
#> colData names(0):
#> 
#> $model_matrix
#>       Intercept
#>  [1,]         1
#>  [2,]         1
#>  [3,]         1
#>  [4,]         1
#>  [5,]         1
#>  [6,]         1
#>  [7,]         1
#>  [8,]         1
#>  [9,]         1
#> [10,]         1
#> attr(,"assign")
#> [1] 0
#> attr(,"intercept_position")
#> [1] 1
#> 
#> $design_formula
#> ~1
#> attr(,"variables")
#> list()
#> attr(,"factors")
#> integer(0)
#> attr(,"term.labels")
#> character(0)
#> attr(,"order")
#> integer(0)
#> attr(,"intercept")
#> [1] 1
#> attr(,"response")
#> [1] 0
#> attr(,"predvars")
#> list()
#> attr(,"dataClasses")
#> named character(0)
#> attr(,"constructed_from")
#> [1] "formula"
#> 

 # Fitting a whole matrix
 model_matrix <- cbind(1, rnorm(5))
 true_Beta <- cbind(rnorm(n = 30), rnorm(n = 30, mean = 3))
 sf <- exp(rnorm(n = 5, mean = 0.7))
 model_matrix
#>      [,1]        [,2]
#> [1,]    1  0.37739565
#> [2,]    1  0.13333636
#> [3,]    1  0.80418951
#> [4,]    1 -0.05710677
#> [5,]    1  0.50360797
 Y <- matrix(rnbinom(n = 30 * 5, mu = sf * exp(true_Beta %*% t(model_matrix)), size = 1/2.4),
             nrow = 30, ncol = 5)

 fit <- glm_gp(Y, design = model_matrix, size_factors = sf, verbose = TRUE)
#> Make initial dispersion estimate
#> Make initial beta estimate
#> Estimate beta
#> Estimate dispersion
#> Fit dispersion trend
#> Shrink dispersion estimates
#> Estimate beta again
 summary(fit)
#> glmGamPoiFit object:
#> The data had 30 rows and 5 columns.
#> A model with 2 coefficient was fitted.
#> 
#> Beta:
#>           Min 1st Qu. Median 3rd Qu.   Max
#> Coef_1 -70.84 -0.3757  0.656    1.58  3.98
#> Coef_2  -4.38 -0.0224  0.800    2.63 89.26
#> 
#> deviance:
#>       Min 1st Qu. Median 3rd Qu.  Max
#>  7.99e-11     2.3   3.61    6.14 11.1
#> 
#> overdispersion:
#>  Min 1st Qu. Median 3rd Qu. Max
#>    0   0.987    1.9    4.21 8.3
#> 
#> Shrunken quasi-likelihood overdispersion:
#>   Min 1st Qu. Median 3rd Qu.  Max
#>  1.19    1.23   1.29    1.42 1.66
#> 
#> size_factors:
#>   Min 1st Qu. Median 3rd Qu.  Max
#>  1.53    3.84   5.46    7.07 7.38
#> 
#> Mu:
#>       Min 1st Qu. Median 3rd Qu. Max
#>  4.02e-33    2.61     10    40.4 615
#> 

 # Fitting a model with covariates
 data <- data.frame(fav_food = sample(c("apple", "banana", "cherry"), size = 50, replace = TRUE),
 city = sample(c("heidelberg", "paris", "new york"), size = 50, replace = TRUE),
 age = rnorm(n = 50, mean = 40, sd = 15))
 Y <- matrix(rnbinom(n = 100 * 50, mu = 3, size = 1/3.1), nrow = 100, ncol = 50)
 fit <- glm_gp(Y, design = ~ fav_food + city + age, col_data = data)
 summary(fit)
#> glmGamPoiFit object:
#> The data had 100 rows and 50 columns.
#> A model with 6 coefficient was fitted.
#> The design formula is: Y~fav_food + city + age
#> 
#> Beta:
#>                  Min 1st Qu.  Median 3rd Qu.  Max
#>      Intercept -1.52   0.376  0.8652   1.519 2.96
#> fav_foodbanana -1.64  -0.457  0.0465   0.538 1.82
#> fav_foodcherry -1.44  -0.484 -0.0158   0.475 2.16
#> ...
#> 
#> deviance:
#>   Min 1st Qu. Median 3rd Qu.  Max
#>  27.4    37.3   40.6    43.8 53.2
#> 
#> overdispersion:
#>   Min 1st Qu. Median 3rd Qu. Max
#>  1.33    2.54   3.11    3.62 5.3
#> 
#> Shrunken quasi-likelihood overdispersion:
#>    Min 1st Qu. Median 3rd Qu.  Max
#>  0.792   0.912  0.977    1.03 1.17
#> 
#> size_factors:
#>    Min 1st Qu. Median 3rd Qu.  Max
#>  0.713   0.851  0.988    1.09 1.58
#> 
#> Mu:
#>    Min 1st Qu. Median 3rd Qu.  Max
#>  0.141    1.58   2.54    3.96 17.8
#> 

 # Specify 'ridge_penalty' to penalize extreme Beta coefficients
 fit_reg <- glm_gp(Y, design = ~ fav_food + city + age, col_data = data, ridge_penalty = 1.5)
 summary(fit_reg)
#> glmGamPoiFit object:
#> The data had 100 rows and 50 columns.
#> A model with 6 coefficient was fitted.
#> The design formula is: Y~fav_food + city + age
#> 
#> Beta:
#>                    Min  1st Qu.    Median 3rd Qu.    Max
#>      Intercept -0.8814  0.57505  0.990486 1.49495 2.8742
#> fav_foodbanana -0.0226 -0.00858  0.000350 0.00894 0.0280
#> fav_foodcherry -0.0240 -0.00696 -0.000732 0.00783 0.0427
#> ...
#> 
#> deviance:
#>   Min 1st Qu. Median 3rd Qu.  Max
#>  28.3    36.8   40.4    43.2 51.2
#> 
#> overdispersion:
#>   Min 1st Qu. Median 3rd Qu.  Max
#>  1.61    2.96   3.52    4.07 6.42
#> 
#> Shrunken quasi-likelihood overdispersion:
#>    Min 1st Qu. Median 3rd Qu.  Max
#>  0.785   0.915  0.983    1.03 1.17
#> 
#> size_factors:
#>    Min 1st Qu. Median 3rd Qu.  Max
#>  0.713   0.851  0.988    1.09 1.58
#> 
#> Mu:
#>    Min 1st Qu. Median 3rd Qu.  Max
#>  0.425    2.19   2.83    3.58 13.4
#> 

```
