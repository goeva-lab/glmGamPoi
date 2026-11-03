# Make a quick first guess where reasonable beta would be

Make a quick first guess where reasonable beta would be

## Usage

``` r
estimate_betas_roughly(
  Y,
  model_matrix,
  offset_matrix,
  pseudo_count = 1,
  ridge_penalty = NULL
)
```

## Value

a matrix with one column for each coefficient
