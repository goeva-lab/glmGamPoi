# configure optimization defaults

configure optimization defaults

## Usage

``` r
set_perf_optim_opts(v2 = TRUE)
set_perf_optim_opts(offset_as_vec = TRUE, do_parallel = 8L, delay_mu = TRUE, use_nr_overdisp_impl = TRUE)

get_perf_optim_opts()
```

## Arguments

- v2:

  boolean, coarse grained flag to enable all implemented optimizations.

  calling `set_perf_optim_opts(v2 = TRUE)` is equivalent to:

      set_perf_optim_opts(offset_as_vec = TRUE, do_parallel = parallel::detectCores(), delay_mu = TRUE, use_nr_overdisp_impl = TRUE)

- offset_as_vec:

  boolean, enables representing the model offsets as a per-gene vector
  when they do not vary across any one observation

- do_parallel:

  integer, sets level of parallelization (see details)

- delay_mu:

  boolean, enables representing `Mu` via a
  [thunk](https://en.wikipedia.org/wiki/Thunk) interface instead of
  eagerly computing a dense representation

- use_nr_overdisp_impl:

  boolean, enables using an alternate overdispersion estimation method
  implemented in C++

## Value

- `set_perf_optim_opts` returns NULL

- `get_perf_optim_opts` returns the currently configured list of
  optimization-relevant defaults

## Details

setting `do_parallel` to `0` disables parallelism entirely, while
setting it to `1` uses the parallelism code path but does not actually
use more than one thread (e.g. this is only useful for testing).

currently, the following steps support parallelization:

- beta estimation

- overdispersion estimation, ***but*** only when `use_nr_overdisp_impl`
  is `TRUE`
