# Convert glmGamPoi object to a list

Convert glmGamPoi object to a list

## Usage

``` r
# S3 method for class 'glmGamPoi2'
as.list(x, ...)
```

## Arguments

- x:

  an object with class glmGamPoi

- ...:

  not used

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

- `Mu`:

  a matrix with the same dimensions as `dim(data)`. If the calculation
  happened on disk, than `Mu` is a `HDF5Matrix`. It contains the
  estimated mean value for each gene and sample.

- `size_factors`:

  a vector with length `ncol(data)`. The size factors are the inferred
  correction factors for different sizes of each sample. They are also
  sometimes called the exposure factor.

- `model_matrix`:

  a matrix with dimensions `ncol(data) x n_coefficients`. It is build
  based on the `design` argument.
