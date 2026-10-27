# Estimate the Betas for Fixed Dispersions

Estimate the Betas for Fixed Dispersions

## Usage

``` r
estimate_betas_group_wise(
  Y,
  offset_matrix,
  dispersions,
  beta_group_init = NULL,
  beta_mat_init = NULL,
  groups,
  model_matrix,
  max_iter = 100,
  do_parallel
)
```

## Value

a list with three elements

- `Beta` a matrix with one column per group and a row for each gene

- `iterations` the number of iterations from the Newton-Raphson method

- `deviances` the deviance for each gene (sum of the deviance per group)
