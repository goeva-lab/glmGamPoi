# Extract Residuals of Gamma Poisson Model

Extract Residuals of Gamma Poisson Model

## Usage

``` r
# S3 method for class 'glmGamPoi2'
residuals(
  object,
  type = c("deviance", "pearson", "randomized_quantile", "working", "response"),
  feat.sub = NULL,
  obs.sub = NULL,
  ...
)
```

## Arguments

- object:

  a fit of type `glmGamPoi`. It is usually produced with a call to
  [`glm_gp()`](glm_gp.md).

- type:

  the type of residual that is calculated. See details for more
  information. Default: `"deviance"`.

- feat.sub:

  subset of features to compute residuals for, ignored if `NULL`.
  Default: `NULL`

- obs.sub:

  subset of observations to compute residuals for, ignored if `NULL`.
  Default: `NULL`

- ...:

  currently ignored.

## Value

a matrix with the same size as `fit$data`. If `fit$data` contains a
`DelayedArray` than the result will be a `DelayedArray` as well.

## Details

This method can calculate a range of different residuals:

- deviance:

  The deviance for the Gamma-Poisson model is \$\$dev = 2 \* (1/theta \*
  \log((1 + m \* theta) / (1 + y \* theta)) - y \log((m + y \* theta) /
  (y + y \* m \* theta)))\$\$ and the residual accordingly is \$\$res =
  sign(y - m) \sqrt{dev}.\$\$

- pearson:

  The Pearson residual is \\res = (y - m) / \sqrt{m + m^2 \* theta}.\\

- randomized_quantile:

  The randomized quantile residual was originally developed by Dunn &
  Smyth, 1995. Please see that publication or
  [`statmod::qresiduals()`](https://rdrr.io/pkg/statmod/man/qresiduals.html)
  for more information.

- working:

  The working residuals are \\res = (y - m) / m\\.

- response:

  The response residuals are \\res = y - m\\

## See also

[`glm_gp()`](glm_gp.md) and \`stats::residuals.glm()
