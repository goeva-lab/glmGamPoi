# Estimate the Betas for Fixed Dispersions

Estimate the Betas for Fixed Dispersions

## Usage

``` r
estimate_betas_fisher_scoring(
  Y,
  model_matrix,
  offset_matrix,
  dispersions,
  beta_mat_init,
  ridge_penalty,
  try_recovering_convergence_problems = TRUE,
  max_iter = 1000,
  do_parallel
)
```

## Value

a list with two elements

- `Beta` a matrix with one column for each coefficient

- `iterations` the number of iterations
