# Predict 'link' or 'response' values for Gamma-Poisson GLMs

Predict `mu` (i.e., `type = "response"`) or `log(mu)` (i.e.,
`type = "link"`) from a 'glmGamPoi' fit (created by `glm_gp(...)`) with
the corresponding estimate of the standard error. If `newdata` is
`NULL`, `mu` is returned for the original input data.

## Usage

``` r
# S3 method for class 'glmGamPoi2'
predict(
  object,
  newdata = NULL,
  type = c("link", "response"),
  se.fit = FALSE,
  offset = mean(object$Offset),
  perf_optim = attr(object, "perf_optim"),
  on_disk = NULL,
  verbose = FALSE,
  feat.sub = NULL,
  obs.sub = NULL,
  ...
)
```

## Arguments

- object:

  a `glmGamPoi` fit object (produced by [`glm_gp()`](glm_gp.md)).

- newdata:

  a specification of the new data for which the expression for each gene
  is predicted. `newdata` should be a

  data.frame

  : if the original fit was specified with a `formula`, provide a
    `data.frame` with one column for each variable in the formula. For
    example, if `glm_gp(se, design = ~ age + batch + treatment)`, then
    the data.frame needs a `age`, `batch`, and `treatment` column that
    contain the same data types as the original fit.

  vector

  : if the original fit was specified using a vector, you need to again
    provide a vector with the same format.

  matrix

  : if `newdata` is a matrix, it is applied directly as
    `Mu <- exp(object$Beta %*% t(newdata) + object$offset_matrix)`. So
    make sure, that it is constructed correctly.

  NULL

  : if `newdata` is `NULL`, the predicted values for the original input
    data are returned.

- type:

  either 'link' or 'response'. The default is 'link', which returns the
  predicted values before the link function
  ([`exp()`](https://rdrr.io/r/base/Log.html)) is applied. Thus, the
  values can be positive and negative numbers. However, often the
  predicted values are easier to interpret **after** the link function
  is applied (i.e., type = "response"), because then the values are on
  the same scale as the original counts.

- se.fit:

  boolean that indicates if in addition to the mean the standard error
  of the mean is returned.

- offset:

  count models (in particular for sequencing experiments) usually have a
  sample specific size factor (`offset = log(size factor)`). It defines
  how big we expect the predicted results are. If `newdata` is `NULL`,
  the `offset` is ignored, because the
  [`predict()`](https://rdrr.io/r/stats/predict.html) returns a result
  based on the pre-calculated `object$Mu`. If `newdata` is not `NULL`,
  by default the `offset` is `mean(object$Offset)`, which puts the in
  the same size as the average sample.

- perf_optim:

  see [`glm_gp()`](glm_gp.md) parameter of same name

- on_disk:

  a boolean that indicates if the results are `HDF5Matrix`'s from the
  `HDF5Array` package. If `newdata` is `NULL`, `on_disk` is ignored.
  Otherwise, if `on_disk = NULL`, the result is calculated on disk
  depending if `offset` is stored on disk.

- verbose:

  a boolean that indicates if information about the individual steps are
  printed while predicting. Default: `FALSE`.

- feat.sub:

  subset of features to compute residuals for, ignored if `NULL`.
  Default: `NULL`

- obs.sub:

  subset of observations to compute residuals for, ignored if `NULL`.
  Default: `NULL`

- ...:

  currently ignored.

## Value

If `se.fit == FALSE`, a matrix with dimensions
`nrow(object$data) x nrow(newdata)`.\
If `se.fit == TRUE`, a list with three entries

- fit:

  the predicted values as a matrix with dimensions
  `nrow(object$data) x nrow(newdata)`. This is what would be returned if
  `se.fit == FALSE`.

- se.fit:

  the associated standard errors for each `fit`. Also a matrix with
  dimensions `nrow(object$data) x nrow(newdata)`.

- residual.scale:

  Currently fixed to 1. In the future, this might become the values from
  `object$overdispersion_shrinkage_list$ql_disp_shrunken`.

## Details

For `se.fit = TRUE`, the function sticks very close to the behavior of
[`stats::predict.glm()`](https://rdrr.io/r/stats/predict.glm.html) for
fits from [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html).

**Note**: If `type = "link"`, the results are computed using the natural
logarithm as the link function. This differs from the `lfc` estimate
provided by [`test_de`](test_de.md), which are on the log2 scale.

## See also

[`stats::predict.lm()`](https://rdrr.io/r/stats/predict.lm.html) and
[`stats::predict.glm()`](https://rdrr.io/r/stats/predict.glm.html)

## Examples

``` r

 set.seed(1)
 # The simplest example
 y <- rnbinom(n = 10, mu = 3, size = 1/2.4)
 fit <- glm_gp(y, size_factors = FALSE)
 predict(fit, type = "response")
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,]  3.9  3.9  3.9  3.9  3.9  3.9  3.9  3.9  3.9   3.9
 predict(fit, type = "link", se.fit = TRUE)
#> $fit
#>          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]
#> [1,] 1.360977 1.360977 1.360977 1.360977 1.360977 1.360977 1.360977 1.360977
#>          [,9]    [,10]
#> [1,] 1.360977 1.360977
#> 
#> $se.fit
#>           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#> [1,] 0.4433902 0.4433902 0.4433902 0.4433902 0.4433902 0.4433902 0.4433902
#>           [,8]      [,9]     [,10]
#> [1,] 0.4433902 0.4433902 0.4433902
#> 
#> $residual.scale
#> [1] 1
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

 head(predict(fit, type = "response"))
#>           [,1]      [,2]       [,3]      [,4]       [,5]
#> [1,] 33.475583 2.0574718 615.186930  1.584906 99.2279246
#> [2,]  0.455028 0.3188184   0.118600  1.640317  0.3831607
#> [3,] 15.165606 6.6531144   8.963827 23.754177 16.2688751
#> [4,] 13.633900 2.4374501  38.725366  4.319617 23.2661147
#> [5,] 19.468804 2.1342826 130.061943  2.582380 42.7843259
#> [6,] 14.617643 3.8217793  21.359215  9.111371 20.4934837
 pred <- predict(fit, type = "link", se.fit = TRUE, verbose = TRUE)
#> 'newdata' is NULL, use 'Mu = object$Mu'
#> 'type = "link"', return 'fit = log(Mu)'
#> 'se.fit' is TRUE, calculate the standard error for each Mu estimate
 head(pred$fit)
#>            [,1]       [,2]      [,3]      [,4]       [,5]
#> [1,]  3.5108163  0.7214780  6.421926 0.4605253  4.5974195
#> [2,] -0.7873964 -1.1431337 -2.131998 0.4948892 -0.9593008
#> [3,]  2.7190301  1.8950851  2.193197 3.1677584  2.7892538
#> [4,]  2.6125593  0.8909525  3.656495 1.4631667  3.1469980
#> [5,]  2.9688134  0.7581306  4.868011 0.9487115  3.7561718
#> [6,]  2.6822292  1.3407161  3.061483 2.2095232  3.0201070
 head(pred$se.fit)
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 0.6781300 0.8891160 1.1899997 1.2214878 0.7343845
#> [2,] 1.1828846 1.2359141 2.3815916 1.7109261 1.4376720
#> [3,] 0.4991307 0.6150965 0.9063817 0.8425566 0.5591235
#> [4,] 0.6013115 0.7677328 1.0623098 1.0507128 0.6587581
#> [5,] 0.5333646 0.7024383 0.9260574 0.9631675 0.5741142
#> [6,] 0.2372328 0.3155340 0.4005399 0.4304783 0.2514965


 # Fitting a model with covariates
 data <- data.frame(fav_food = sample(c("apple", "banana", "cherry"), size = 50, replace = TRUE),
                    city = sample(c("heidelberg", "paris", "new york"), size = 50, replace = TRUE),
                    age = rnorm(n = 50, mean = 40, sd = 15))
 Y <- matrix(rnbinom(n = 4 * 50, mu = 3, size = 1/3.1), nrow = 4, ncol = 50)
 fit <- glm_gp(Y, design = ~ fav_food + city + age, col_data = data)
#> Warning: 1 columns contain too many zeros to calculate a size factor. The size factor will be fixed to 0.001
 predict(fit)[, 1:3]
#>           [,1]     [,2]     [,3]
#> [1,] 0.4308067 2.068513 3.078971
#> [2,] 0.4521625 2.369426 1.858950
#> [3,] 1.1158029 2.329241 1.596998
#> [4,] 0.7281115 1.762526 2.519676

 nd <- data.frame(fav_food = "banana", city = "paris", age = 29)
 predict(fit, newdata = nd)
#>               1
#> [1,] 0.08024906
#> [2,] 0.81635137
#> [3,] 0.86360403
#> [4,] 0.33975910

 nd <- data.frame(fav_food = "banana", city = "paris", age = 29:40)
 predict(fit, newdata = nd, se.fit = TRUE, type = "response")
#> $fit
#>             1        2        3        4        5        6        7        8
#> [1,] 1.083557 1.125058 1.168149 1.212890 1.259345 1.307579 1.357661 1.409661
#> [2,] 2.262231 2.256000 2.249786 2.243590 2.237410 2.231247 2.225102 2.218973
#> [3,] 2.371693 2.295073 2.220929 2.149180 2.079749 2.012561 1.947544 1.884627
#> [4,] 1.404609 1.418367 1.432259 1.446287 1.460453 1.474757 1.489202 1.503788
#>             9       10       11       12
#> [1,] 1.463652 1.519711 1.577918 1.638353
#> [2,] 2.212862 2.206767 2.200689 2.194627
#> [3,] 1.823742 1.764825 1.707811 1.652638
#> [4,] 1.518516 1.533389 1.548408 1.563574
#> 
#> $se.fit
#>              1         2         3         4         5         6         7
#> [1,] 0.5579656 0.5702174 0.5830179 0.5964286 0.6105174 0.6253588 0.6410345
#> [2,] 1.2372901 1.2151296 1.1940080 1.1739592 1.1550165 1.1372127 1.1205793
#> [3,] 0.8656589 0.8268159 0.7901980 0.7557049 0.7232390 0.6927055 0.6640123
#> [4,] 0.7748521 0.7708379 0.7672470 0.7641195 0.7614966 0.7594209 0.7579357
#>              8         9        10        11        12
#> [1,] 0.6576335 0.6752523 0.6939954 0.7139753 0.7353124
#> [2,] 1.1051465 1.0909424 1.0779929 1.0663207 1.0559456
#> [3,] 0.6370696 0.6117897 0.5880874 0.5658793 0.5450841
#> [4,] 0.7570852 0.7569137 0.7574656 0.7587853 0.7609159
#> 
#> $residual.scale
#> [1] 1
#> 

```
