# glmGamPoi2 performance optimizations

## `glmGamPoi2` fork

The `glmGamPoi2` package is a fork of the [existing glmGamPoi
package](https://github.com/const-ae/glmGamPoi) which aims to optimize
performance and memory use. Multiple “knobs” are provided which can be
used to enable the various implementations provided.

Primarily, the `glm_gp` entrypoint now has an added `perf_optim`
parameter, which can be used to enable/disable specific optimizations at
will.

## Currently supported parameters

We demonstrate below the various parameters currently introduced via the
`perf_optim` parameter. First, we load the library, initialize a
synthetic dataset, and generate a reference object:

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`glmGamPoi2`](https://github.com/goeva-lab/glmGamPoi)`)`\
\
`counts`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`[`rnbinom`](https://rdrr.io/r/stats/NegBinomial.html)`(``256``, mu ``=`` ``4``, size ``=`` ``0.25``)``, ``32``, ``8``)`\
`fit_ref`` ``<-`` `[`glm_gp`](../reference/glm_gp.md)`(``counts``, perf_optim ``=`` ``FALSE``)`

### `offset_as_vec`

The `offset_as_vec` parameter allows representing the
per-gene/per-observation offsets as a vector when they are not
configured to vary across any one observation:

\
[`head`](https://rdrr.io/r/utils/head.html)`(``fit_ref``[[``"Offset"``]``]``)`` ``# we can see this is a matrix`\
`#>            [,1]      [,2]     [,3]       [,4]      [,5]       [,6]        [,7]`\
`#> [1,] -0.1999854 0.3962198 0.161028 -0.3650651 0.1145079 -0.1784792 -0.09680116`\
`#> [2,] -0.1999854 0.3962198 0.161028 -0.3650651 0.1145079 -0.1784792 -0.09680116`\
`#> [3,] -0.1999854 0.3962198 0.161028 -0.3650651 0.1145079 -0.1784792 -0.09680116`\
`#> [4,] -0.1999854 0.3962198 0.161028 -0.3650651 0.1145079 -0.1784792 -0.09680116`\
`#> [5,] -0.1999854 0.3962198 0.161028 -0.3650651 0.1145079 -0.1784792 -0.09680116`\
`#> [6,] -0.1999854 0.3962198 0.161028 -0.3650651 0.1145079 -0.1784792 -0.09680116`\
`#>           [,8]`\
`#> [1,] 0.1685752`\
`#> [2,] 0.1685752`\
`#> [3,] 0.1685752`\
`#> [4,] 0.1685752`\
`#> [5,] 0.1685752`\
`#> [6,] 0.1685752`\
\
`fit_vec`` ``<-`` `[`glm_gp`](../reference/glm_gp.md)`(``counts``, perf_optim ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``offset_as_vec ``=`` ``TRUE``)``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``fit_vec``[[``"Offset"``]``]``)`` ``# we can see this is now a vector`\
`#> [1] -0.1999854  0.3962198  0.1610280 -0.3650651  0.1145079 -0.1784792`

### `delay_mu`

The `delay_mu` parameter allows representing the prediction matrix `Mu`
via a [thunk](https://en.wikipedia.org/wiki/Thunk) interface instead of
eagerly computing a dense representation:

\
[`head`](https://rdrr.io/r/utils/head.html)`(``fit_ref``[[``"Mu"``]``]``)`` ``# we can see this is a matrix`\
`#>          [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]     [,8]`\
`#> [1,] 1.556401 2.825206 2.233097 1.319557 2.131593 1.590236 1.725575 2.250015`\
`#> [2,] 3.140700 5.701054 4.506222 2.662768 4.301394 3.208977 3.482081 4.540360`\
`#> [3,] 2.719190 4.935922 3.901447 2.305401 3.724109 2.778303 3.014755 3.931003`\
`#> [4,] 3.871690 7.027959 5.555034 3.282520 5.302532 3.955857 4.292526 5.597117`\
`#> [5,] 4.504430 8.176520 6.462878 3.818973 6.169111 4.602352 4.994042 6.511839`\
`#> [6,] 3.097058 5.621833 4.443605 2.625767 4.241623 3.164385 3.433695 4.477269`\
\
`fit_delayed`` ``<-`` `[`glm_gp`](../reference/glm_gp.md)`(``counts``, perf_optim ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``delay_mu ``=`` ``TRUE``)``)`\
\
[`is.function`](https://rdrr.io/r/base/is.function.html)`(``fit_delayed``[[``"Mu"``]``]``)`` ``# we can see this is now a _function_`\
`#> [1] TRUE`\
\
`# demonstration of the thunk'd Mu API`\
[`all.equal`](https://rdrr.io/r/base/all.equal.html)`(``fit_delayed``[[``"Mu"``]``]``(``)``, ``fit_ref``[[``"Mu"``]``]``)`\
`#> [1] TRUE`\
[`all.equal`](https://rdrr.io/r/base/all.equal.html)`(``fit_delayed``[[``"Mu"``]``]``(``i ``=`` ``1L``)``, ``fit_ref``[[``"Mu"``]``]``[``1``, , drop ``=`` ``FALSE``]``)`\
`#> [1] TRUE`\
[`all.equal`](https://rdrr.io/r/base/all.equal.html)`(``fit_delayed``[[``"Mu"``]``]``(``j ``=`` ``1L``)``, ``fit_ref``[[``"Mu"``]``]``[``, ``1``, drop ``=`` ``FALSE``]``)`\
`#> [1] TRUE`\
[`all.equal`](https://rdrr.io/r/base/all.equal.html)`(``fit_delayed``[[``"Mu"``]``]``(``1L``, ``1L``)``, ``fit_ref``[[``"Mu"``]``]``[``1``, ``1``, drop ``=`` ``FALSE``]``)`\
`#> [1] TRUE`

### `use_nr_overdisp_impl` and `do_parallel`

The `use_nr_overdisp_impl` and `do_parallel` flags do not change output
shape, but modify behaviour of internal routines:

- `use_nr_overdisp_impl`: enables the use of an alternative
  overdispersion estimation routine written in `C++`
- `do_parallel`: enables parallelization of various steps across the
  package

Please note that while testing/benchmarking has found that these do not
*significantly* affect result outputs, there are edge cases (especially
for ill-conditioned fits with high deviance on the resulting fit) where
enabling `use_nr_overdisp_impl` can lead to slightly different results.

\
`fit_alt`` ``<-`` `[`glm_gp`](../reference/glm_gp.md)`(``counts``, perf_optim ``=`` `[`list`](https://rdrr.io/r/base/list.html)`(``do_parallel ``=`` ``2L``, use_nr_overdisp_impl ``=`` ``TRUE``)``)`\
\
`# we can see these have no effect on resulting model parameters`\
[`all.equal`](https://rdrr.io/r/base/all.equal.html)`(``fit_alt``[[``"overdispersions"``]``]``, ``fit_ref``[[``"overdispersions"``]``]``)`\
`#> [1] TRUE`\
[`all.equal`](https://rdrr.io/r/base/all.equal.html)`(``fit_alt``[[``"Beta"``]``]``, ``fit_ref``[[``"Beta"``]``]``)`\
`#> [1] TRUE`

## Changing package-wide defaults

Finally, to improve usability when one wishes to change parameters
package wide, the `set_perf_optim_opts` helper is provided. Notably, it
has a `v2` parameter, which allows for a coarse-grain toggle of all
recommended optimizations:

\
[`set_perf_optim_opts`](../reference/set_perf_optim_opts.md)`(``v2 ``=`` ``TRUE``)`\
\
[`get_perf_optim_opts`](../reference/set_perf_optim_opts.md)`(``)`\
`#> $offset_as_vec`\
`#> [1] TRUE`\
`#> `\
`#> $do_parallel`\
`#> [1] 4`\
`#> `\
`#> $delay_mu`\
`#> [1] TRUE`\
`#> `\
`#> $use_nr_overdisp_impl`\
`#> [1] TRUE`\
`#> `\
`#> $cast_dgC_Y_to_dgR`\
`#> [1] FALSE`\
\
`fit_v2`` ``<-`` `[`glm_gp`](../reference/glm_gp.md)`(``counts``)`\
\
`# we can see that all optimizations have now been enabled`\
[`head`](https://rdrr.io/r/utils/head.html)`(``fit_v2``[[``"Offset"``]``]``)`\
`#> [1] -0.1999854  0.3962198  0.1610280 -0.3650651  0.1145079 -0.1784792`\
[`is.function`](https://rdrr.io/r/base/is.function.html)`(``fit_v2``[[``"Mu"``]``]``)`\
`#> [1] TRUE`
