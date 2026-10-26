.perf_optim.offset_as_vec <- FALSE
.perf_optim.do_parallel <- 0L
.perf_optim.delay_mu <- FALSE
.perf_optim.use_nr_overdisp_impl <- FALSE
.perf_optim.cast_dgC_Y_to_dgR <- FALSE

handle_perf_optim_parameter <- function(param) {
  default_opts <- list(
    "offset_as_vec" = list(.perf_optim.offset_as_vec, function(e) (is.logical(e) && (length(e) == 1)), "logical of length 1"),
    "do_parallel" = list(.perf_optim.do_parallel, function(e) (is.integer(e) && (length(e) == 1)), "integer of length 1"),
    "delay_mu" = list(.perf_optim.delay_mu, function(e) (is.logical(e) && (length(e) == 1)), "logical of length 1"),
    "use_nr_overdisp_impl" = list(.perf_optim.use_nr_overdisp_impl, function(e) (is.logical(e) && (length(e) == 1)), "logical of length 1"),
    "cast_dgC_Y_to_dgR" = list(.perf_optim.cast_dgC_Y_to_dgR, function(e) (is.logical(e) && (length(e) == 1)), "logical of length 1")
  )

  out <- if (is.logical(param) && (length(param) == 1)) {
    list(
      "offset_as_vec" = param,
      "do_parallel" = if (param) parallel::detectCores() else 0L,
      "delay_mu" = param,
      "use_nr_overdisp_impl" = param,
      "cast_dgC_Y_to_dgR" = param
    )
  } else if (is.list(param) && (length(union(names(default_opts), names(param)))) == length(default_opts)) {
    for (nm in names(param)) {
      if (!default_opts[[nm]][[2]](param[[nm]])) {
        stop(sprintf(
          "got perf_optim parameter (`%s`, at index - `%s`) of wrong type/shape, must be %s",
          param[[nm]],
          nm,
          default_opts[[nm]][[3]]
        ))
      }
    }
    utils::modifyList(lapply(default_opts, function(e) e[[1]]), param)
  } else {
    stop(sprintf(
      "got perf_optim parameter (`%s`) of wrong type/shape, must be list with names `%s` (or a subset thereof) or logical of length 1",
      paste0(sprintf("%s=%s", names(param), param), collapse = ", "),
      paste0(names(default_opts), collapse = ", ")
    ))
  }

  if (out[["do_parallel"]] == 1L) {
    warning(paste0(
      c(
        "got perf_optim$do_parallel=1, meaning no parallelization is enabled while still incurring cost of setting up parallelization-safe machinery.",
        "this is only useful for internal testing, if you wish to just disable parallelization, you should instead set perf_optim$do_parallel=0."
      ),
      collapse = "\n"
    ))
  }
  n_cores <- parallel::detectCores()
  if (out[["do_parallel"]] > n_cores) {
    warning(paste0(
      c(
        sprintf("got perf_optim$do_parallel=%s, which is greater than the number of cores detected by `parallel::detectCores() (%s)`.", out[["do_parallel"]], n_cores),
        "unless parallel::detectCores is returning an incorrect value, there is generally no reason to do this as spawning more threads than physical cores will usually harm performance."
      ),
      collapse = "\n"
    ))
  }

  out
}

#' configure optimization defaults
#'
#' @details
#' setting `do_parallel` to `0` disables parallelism entirely,
#' while setting it to `1` uses the parallelism code path but does
#' not actually use more than one thread (e.g. this is only useful for testing).
#'
#' currently, the following steps support parallelization:
#' - beta estimation
#' - overdispersion estimation, ***but*** only when `use_nr_overdisp_impl` is `TRUE`
#'
#' @usage
#' set_perf_optim_opts(v2 = TRUE)
#' set_perf_optim_opts(offset_as_vec = TRUE, do_parallel = 8L, delay_mu = TRUE, use_nr_overdisp_impl = TRUE)
#'
#' @param v2 boolean, coarse grained flag to enable all implemented optimizations.
#'
#'  calling `set_perf_optim_opts(v2 = TRUE)` is equivalent to:
#'
#'  ```
#'  set_perf_optim_opts(offset_as_vec = TRUE, do_parallel = parallel::detectCores(), delay_mu = TRUE, use_nr_overdisp_impl = TRUE)
#'  ```
#' @param offset_as_vec boolean, enables representing the model offsets as a per-gene vector when they do not vary across any one observation
#' @param do_parallel integer, sets level of parallelization (see details)
#' @param delay_mu boolean, enables representing `Mu` via a [thunk](https://en.wikipedia.org/wiki/Thunk) interface instead of eagerly computing a dense representation
#' @param use_nr_overdisp_impl boolean, enables using an alternate overdispersion estimation method implemented in C++
#' @returns
#'  - `set_perf_optim_opts` returns NULL
#'  - `get_perf_optim_opts` returns the currently configured list of optimization-relevant defaults
#' @export
set_perf_optim_opts <- function(v2, offset_as_vec, do_parallel, delay_mu, use_nr_overdisp_impl, cast_dgC_Y_to_dgR = FALSE) {
  opts <- handle_perf_optim_parameter(list(
    "offset_as_vec" = if (missing(offset_as_vec)) v2 else offset_as_vec,
    "do_parallel" = if (missing(do_parallel)) {
      if (v2) parallel::detectCores() else 0L
    } else {
      do_parallel
    },
    "delay_mu" = if (missing(delay_mu)) v2 else delay_mu,
    "use_nr_overdisp_impl" = if (missing(use_nr_overdisp_impl)) v2 else delay_mu,
    "cast_dgC_Y_to_dgR" = if (missing(cast_dgC_Y_to_dgR)) FALSE else cast_dgC_Y_to_dgR
  ))
  for (nm in names(opts)) {
    utils::assignInMyNamespace(sprintf(".perf_optim.%s", nm), opts[[nm]])
  }
}

#' @rdname set_perf_optim_opts
#' @export
get_perf_optim_opts <- function() {
  list(
    "offset_as_vec" = .perf_optim.offset_as_vec,
    "do_parallel" = .perf_optim.do_parallel,
    "delay_mu" = .perf_optim.delay_mu,
    "use_nr_overdisp_impl" = .perf_optim.use_nr_overdisp_impl,
    "cast_dgC_Y_to_dgR" = .perf_optim.cast_dgC_Y_to_dgR
  )
}
