test_that("set_perf_optim_opts/get_perf_optim_opts work", {
  tfo <- function(out, offset_as_vec, delay_mu) {
    if (offset_as_vec) {
      expect_vector(out[["Offset"]], numeric(), ncol(out[["data"]]))
      expect_named(out[["Offset"]], colnames(out[["data"]]))
    } else {
      expect_equal(class(out[["Offset"]])[[1]], "matrix")
      expect_equal(dimnames(out[["Offset"]]), dimnames(out[["data"]]))
    }
    if (delay_mu) {
      expect_true(is.function(out[["Mu"]]))
    } else {
      expect_equal(class(out[["Mu"]])[[1]], "matrix")
    }
  }

  prv_opts <- get_perf_optim_opts()
  new_opts <- list(offset_as_vec = !prv_opts[["offset_as_vec"]], do_parallel = 8L * (prv_opts[["do_parallel"]] + 1L), delay_mu = !prv_opts[["delay_mu"]], use_nr_overdisp_impl = !prv_opts[["use_nr_overdisp_impl"]], cast_dgC_Y_to_dgR = FALSE)

  tfo(glm_gp(matrix(rnbinom(n = 128, mu = 3, size = 0.4), 8, 16)), prv_opts[["offset_as_vec"]], prv_opts[["delay_mu"]])

  expect_no_error(do.call(set_perf_optim_opts, new_opts))
  expect_equal(get_perf_optim_opts(), new_opts)
  tfo(glm_gp(matrix(rnbinom(n = 128, mu = 3, size = 0.4), 8, 16)), new_opts[["offset_as_vec"]], new_opts[["delay_mu"]])

  expect_no_error(do.call(set_perf_optim_opts, prv_opts))
  expect_equal(get_perf_optim_opts(), prv_opts)
  tfo(glm_gp(matrix(rnbinom(n = 128, mu = 3, size = 0.4), 8, 16)), prv_opts[["offset_as_vec"]], prv_opts[["delay_mu"]])
})
