# Test for Differential Expression

Conduct a quasi-likelihood ratio test for a Gamma-Poisson fit.

## Usage

``` r
test_de(
  fit,
  contrast,
  reduced_design = NULL,
  full_design = fit$model_matrix,
  subset_to = NULL,
  pseudobulk_by = NULL,
  pval_adjust_method = "BH",
  sort_by = NULL,
  decreasing = FALSE,
  n_max = Inf,
  max_lfc = 10,
  compute_lfc_se = FALSE,
  perf_optim = attr(fit, "perf_optim"),
  verbose = FALSE
)
```

## Arguments

- fit:

  object of class `glmGamPoi`. Usually the result of calling
  `glm_gp(data, ...)`

- contrast:

  The contrast to test. Can be a single column name (quoted or as a
  string) that is removed from the full model matrix of `fit`. Or a
  complex contrast comparing two or more columns: e.g. `A - B`,
  `"A - 3 * B"`, `(A + B) / 2 - C` etc. For complicated experimental
  design that involved nested conditions, you specify the condition
  level to compare using the `cond()` helper function.\
  Only one of `contrast` or `reduced_design` must be specified.

- reduced_design:

  a specification of the reduced design used as a comparison to see what
  how much better `fit` describes the data. Analogous to the `design`
  parameter in [`glm_gp()`](glm_gp.md), it can be either a `formula`, a
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html), or a
  `vector`.\
  Only one of `contrast` or `reduced_design` must be specified.

- full_design:

  option to specify an alternative `full_design` that can differ from
  `fit$model_matrix`. Can be a `formula` or a `matrix`. Default:
  `fit$model_matrix`

- subset_to:

  a vector with the same length as `ncol(fit$data)` or an expression
  that evaluates to such a vector. The expression can reference columns
  from `colData(fit$data)`. A typical use case in single cell analysis
  would be to subset to a specific cell type (e.g.
  `subset_to = cell_type == "T-cells"`). Note that if this argument is
  set a new the model for the `full_design` is re-fit.\
  Default: `NULL` which means that the data is not subset.

- pseudobulk_by:

  *DEPRECTATED*, please use the [pseudobulk](pseudobulk.md) function
  instead.\
  A vector with the same length as `ncol(fit$data)` that is used to
  split the columns into different groups (calls
  [base::split](https://rdrr.io/r/base/split.html)). `pseudobulk_by` can
  also be an expression that evaluates to a vector. The expression can
  reference columns from `colData(fit$data)`.\
  The counts are summed across the groups to create "pseudobulk"
  samples. This is typically used in single cell analysis if the cells
  come from different samples to get a proper estimate of the
  differences. This is particularly powerful in combination with the
  `subset_to` parameter to analyze differences between samples for
  subgroups of cells. Note that this does a fresh fit for both the full
  and the reduced design. Default: `NULL` which means that the data is
  not aggregated.

- pval_adjust_method:

  one of the p-value adjustment method from
  [p.adjust.methods](https://rdrr.io/r/stats/p.adjust.html). Default:
  `"BH"`.

- sort_by:

  the name of the column or an expression used to sort the result. If
  `sort_by` is `NULL` the table is not sorted. Default: `NULL`

- decreasing:

  boolean to decide if the result is sorted increasing or decreasing
  order. Default: `FALSE`.

- n_max:

  the maximum number of rows to return. Default: `Inf` which means that
  all rows are returned

- max_lfc:

  set the maximum absolute log fold change that is returned. Large log
  fold changes occur for lowly expressed genes because the ratio of two
  small numbers can be impractically large. For example, limiting the
  range of log fold changes can clarify the patterns in a volcano plot.
  Default: `10` which corresponds to a thousand-fold (2^10) increase in
  expression.

- compute_lfc_se:

  Compute standard errors for the log fold changes, and add a column
  `lfc_se` to the returned dataframe. Only has an effect when using
  `contrast` instead of `reduced_design`. Default: `FALSE`

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

a `data.frame` with the following columns

- name:

  the rownames of the input data

- pval:

  the p-values of the quasi-likelihood ratio test

- adj_pval:

  the adjusted p-values returned from
  [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)

- f_statistic:

  the F-statistic: \\F = (Dev\_{full} - Dev\_{red}) / (df_1 \*
  disp\_{ql-shrunken})\\

- df1:

  the degrees of freedom of the test:
  `ncol(design) - ncol(reduced_design)`

- df2:

  the degrees of freedom of the fit: `ncol(data) - ncol(design) + df_0`

- lfc:

  the log2-fold change. If the alternative model is specified by
  `reduced_design` will be `NA`.

- lfc_se:

  the standard error of the log2-fold change. Only calculated if
  `compute_lfc_se == TRUE`.

## Details

The `cond()` helper function simplifies the specification of a contrast
for complex experimental designs. Instead of working out which
combination of coefficients corresponds to a research question, you can
simply specify the two conditions that you want to compare.

You can only call the `cond` function inside the `contrast` argument.
The arguments are the selected factor levels for each covariate. To
compare two conditions, simply subtract the two `cond` calls.
Internally, the package calls
[model.matrix](https://rdrr.io/r/stats/model.matrix.html) using the
provided values and the original formula from the fit to produce a
vector. Subtracting two of these vectors produces a contrast vector.
Missing covariates are filled with the first factor level or zero for
numerical covariates.

## References

- Lund, S. P., Nettleton, D., McCarthy, D. J., & Smyth, G. K. (2012).
  Detecting differential expression in RNA-sequence data using
  quasi-likelihood with shrunken dispersion estimates. Statistical
  Applications in Genetics and Molecular Biology, 11(5).
  <https://doi.org/10.1515/1544-6115.1826>.

## See also

[`glm_gp()`](glm_gp.md)

## Examples

``` r
 # Make Data
 Y <- matrix(rnbinom(n = 30 * 100, mu = 4, size = 0.3), nrow = 30, ncol = 100)
 annot <- data.frame(mouse = sample(LETTERS[1:6], size = 100, replace = TRUE),
        celltype = sample(c("Tcell", "Bcell", "Macrophages"), size = 100, replace = TRUE),
        cont1 = rnorm(100), cont2 = rnorm(100, mean = 30))
 annot$condition <- ifelse(annot$mouse %in% c("A", "B", "C"), "ctrl", "treated")
 head(annot)
#>   mouse    celltype       cont1    cont2 condition
#> 1     A       Bcell -0.09198166 29.68407      ctrl
#> 2     E       Tcell  0.71066364 30.24982   treated
#> 3     F Macrophages  0.11891539 30.75372   treated
#> 4     C       Bcell  0.39235734 30.86076      ctrl
#> 5     D Macrophages  0.56161326 30.12915   treated
#> 6     F Macrophages -0.43386063 31.45817   treated
 se <- SummarizedExperiment::SummarizedExperiment(Y, colData = annot)

 # Fit model
 fit <- glm_gp(se, design = ~ condition + celltype + cont1 + cont2)
 # Test with reduced design
 res <- test_de(fit, reduced_design = ~ celltype + cont1 + cont2)
 head(res)
#>    name       pval  adj_pval f_statistic df1      df2 lfc
#> 1 row_1 0.32262026 0.7069035  0.98313465   1 200.6624  NA
#> 2 row_2 0.01055799 0.3167396  6.66272989   1 200.6624  NA
#> 3 row_3 0.31275609 0.7069035  1.02414686   1 200.6624  NA
#> 4 row_4 0.35345175 0.7069035  0.86503247   1 200.6624  NA
#> 5 row_5 0.84719780 0.9689670  0.03722747   1 200.6624  NA
#> 6 row_6 0.43617202 0.8178225  0.60876644   1 200.6624  NA
 # Test with contrast argument, the results are identical
 res2 <- test_de(fit, contrast = conditiontreated)
 head(res2)
#>    name       pval  adj_pval f_statistic df1      df2        lfc
#> 1 row_1 0.32262026 0.7069035  0.98313465   1 200.6624 -0.5155925
#> 2 row_2 0.01055799 0.3167396  6.66272989   1 200.6624 -1.4145875
#> 3 row_3 0.31275609 0.7069035  1.02414686   1 200.6624  0.5810427
#> 4 row_4 0.35345175 0.7069035  0.86503247   1 200.6624  0.5542879
#> 5 row_5 0.84719780 0.9689670  0.03722747   1 200.6624  0.1033688
#> 6 row_6 0.43617202 0.8178225  0.60876644   1 200.6624  0.4541280
 # Test with explicit specification of the conditions
 # The results are still identical
 res3 <- test_de(fit, contrast = cond(condition = "treated", celltype = "Bcell") -
                                    cond(condition = "ctrl", celltype = "Bcell"))
 head(res3)
#>    name       pval  adj_pval f_statistic df1      df2        lfc
#> 1 row_1 0.32262026 0.7069035  0.98313465   1 200.6624 -0.5155925
#> 2 row_2 0.01055799 0.3167396  6.66272989   1 200.6624 -1.4145875
#> 3 row_3 0.31275609 0.7069035  1.02414686   1 200.6624  0.5810427
#> 4 row_4 0.35345175 0.7069035  0.86503247   1 200.6624  0.5542879
#> 5 row_5 0.84719780 0.9689670  0.03722747   1 200.6624  0.1033688
#> 6 row_6 0.43617202 0.8178225  0.60876644   1 200.6624  0.4541280


 # The column names of fit$Beta are valid variables in the contrast argument
 colnames(fit$Beta)
#> [1] "Intercept"           "conditiontreated"    "celltypeMacrophages"
#> [4] "celltypeTcell"       "cont1"               "cont2"              
 # You can also have more complex contrasts:
 # the following compares cont1 vs cont2:
 test_de(fit, cont1 - cont2, n_max = 4)
#>    name      pval  adj_pval f_statistic df1      df2         lfc
#> 1 row_1 0.4796013 0.9570646  0.50164237   1 200.6624 -0.36336041
#> 2 row_2 0.6326804 0.9570646  0.22914405   1 200.6624 -0.21585497
#> 3 row_3 0.8449570 0.9570646  0.03834159   1 200.6624 -0.08005402
#> 4 row_4 0.5521653 0.9570646  0.35464709   1 200.6624  0.27864863
 # You can also sort the output
 test_de(fit, cont1 - cont2, n_max = 4,
         sort_by = "pval")
#>      name       pval  adj_pval f_statistic df1      df2        lfc
#> 15 row_15 0.05822805 0.9570646    3.628501   1 200.6624 -0.8412637
#> 6   row_6 0.07718129 0.9570646    3.155617   1 200.6624  0.7416995
#> 16 row_16 0.09699879 0.9570646    2.780155   1 200.6624 -0.7272020
#> 30 row_30 0.21748931 0.9570646    1.530459   1 200.6624 -0.6968706
 test_de(fit, cont1 - cont2, n_max = 4,
         sort_by = - abs(f_statistic))
#>      name       pval  adj_pval f_statistic df1      df2        lfc
#> 15 row_15 0.05822805 0.9570646    3.628501   1 200.6624 -0.8412637
#> 6   row_6 0.07718129 0.9570646    3.155617   1 200.6624  0.7416995
#> 16 row_16 0.09699879 0.9570646    2.780155   1 200.6624 -0.7272020
#> 30 row_30 0.21748931 0.9570646    1.530459   1 200.6624 -0.6968706

 # If the data has multiple samples, it is a good
 # idea to aggregate the cell counts by samples.
 # This is called forming a "pseudobulk".
 se_reduced <- pseudobulk(se, group_by = vars(mouse, condition, celltype),
                          cont1 = mean(cont1), cont2 = min(cont2))
#> Aggregating assay '1' using 'rowMeans2'.
 fit_reduced <- glm_gp(se_reduced, design = ~ condition + celltype)
 test_de(fit_reduced, contrast = "conditiontreated", n_max = 4)
#>    name       pval  adj_pval f_statistic df1      df2        lfc
#> 1 row_1 0.81391213 0.9254585  0.05540476   1 79832644 -0.1409044
#> 2 row_2 0.02147207 0.2549152  5.28806491   1 79832644 -1.3418906
#> 3 row_3 0.37747047 0.8711526  0.77892384   1 79832644  0.5205163
#> 4 row_4 0.63116221 0.8916901  0.23048801   1 79832644 -0.2868613
 test_de(fit_reduced, contrast = cond(condition = "treated", celltype = "Macrophages") -
                                      cond(condition = "ctrl", celltype = "Macrophages"),
         n_max = 4)
#>    name       pval  adj_pval f_statistic df1      df2        lfc
#> 1 row_1 0.81391213 0.9254585  0.05540476   1 79832644 -0.1409044
#> 2 row_2 0.02147207 0.2549152  5.28806491   1 79832644 -1.3418906
#> 3 row_3 0.37747047 0.8711526  0.77892384   1 79832644  0.5205163
#> 4 row_4 0.63116221 0.8916901  0.23048801   1 79832644 -0.2868613

```
