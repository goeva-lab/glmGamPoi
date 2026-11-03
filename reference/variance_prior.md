# Estimate the scale and df for a Inverse Chisquare distribution that generate the true gene variances

This function implements Smyth's 2004 variance shrinkage. It also
supports covariates that are fitted to log(s2) with natural splines.
This is based on the 2012 Lund et al. quasi-likelihood paper.

## Usage

``` r
variance_prior(s2, df, covariate = NULL, abundance_trend = NULL)
```

## Arguments

- s2:

  vector of observed variances. Must not contain `0`'s.

- df:

  vector or single number with the degrees of freedom

- covariate:

  a vector with the same length as s2. `covariate` is used to regress
  out the trend in `s2`. If `covariate = NULL`, it is ignored.

- abundance_trend:

  logical that decides if the additional abundance trend is fit to the
  data. If `NULL` the abundance trend is fitted if there are more than
  10 observations and the `covariate` is not `NULL`. Default: `NULL`

## Value

a list with three elements:

- variance0:

  estimate of the scale of the inverse Chisquared distribution. If
  covariate is `NULL` a single number, otherwise a vector of
  `length(covariate)`

- df0:

  estimate of the degrees of freedom of the inverse Chisquared
  distribution

- var_pos:

  the shrunken variance estimates: a combination of `s2` and `variance0`

## See also

[`limma::squeezeVar()`](https://rdrr.io/pkg/limma/man/squeezeVar.html)
