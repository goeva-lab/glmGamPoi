# Estimate the Size Factors

Estimate the Size Factors

## Usage

``` r
estimate_size_factors(Y, method, verbose = FALSE)
```

## Arguments

- Y:

  any matrix-like object
  ([`base::matrix()`](https://rdrr.io/r/base/matrix.html),
  `DelayedArray`, `HDF5Matrix`,
  [`Matrix::Matrix()`](https://rdrr.io/pkg/Matrix/man/Matrix.html),
  etc.)

- method:

  one of `c("normed_sum", "deconvolution", "poscounts")`

## Value

a vector with one size factor per column of `Y`
