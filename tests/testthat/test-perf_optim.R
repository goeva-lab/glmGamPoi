test_that("perf_optim$offset_as_vec doesn't change anything but offset matrix", {
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol = 10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE), cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- ~ group + cont1 + cont2

  set.seed(1)
  res <- glm_gp(Y, design = design, col_data = annot, perf_optim = list(offset_as_vec = TRUE))
  expect_vector(res[["Offset"]], ptype = numeric(), size = ncol(Y))

  set.seed(1)
  res_orig <- glm_gp(Y, design = design, col_data = annot)

  res[["Offset"]] <- matrix(res[["Offset"]], nrow = nrow(Y), ncol = ncol(Y), byrow = TRUE)
  attr(res, "perf_optim") <- attr(res_orig, "perf_optim")

  expect_equal(res, res_orig, tolerance = 0)
})

test_that("perf_optim$cast_dgC_Y_to_dgR doesn't change anything", {
  Y <- as(matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol = 10), "CsparseMatrix")
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE), cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- ~ group + cont1 + cont2

  set.seed(1)
  res <- glm_gp(Y, design = design, col_data = annot, perf_optim = list(cast_dgC_Y_to_dgR = TRUE))

  set.seed(1)
  res_orig <- glm_gp(Y, design = design, col_data = annot)

  attr(res, "perf_optim") <- attr(res_orig, "perf_optim")

  expect_equal(res, res_orig, tolerance = 0)

  set.seed(1)
  res_dgr <- glm_gp(Y, design = design, col_data = annot, perf_optim = list(offset_as_vec = TRUE, mu_dropout_thresh = 1e-1, cast_dgC_Y_to_dgR = TRUE))
  expect_s4_class(res_dgr[["Mu"]], "RsparseMatrix")

  set.seed(1)
  res_dgc <- glm_gp(Y, design = design, col_data = annot, perf_optim = list(offset_as_vec = TRUE, mu_dropout_thresh = 1e-1))
  expect_s4_class(res_dgc[["Mu"]], "CsparseMatrix")

  res_dgr[["Mu"]] <- as(res_dgr[["Mu"]], "CsparseMatrix")
  attr(res_dgr, "perf_optim") <- attr(res_dgc, "perf_optim")
  expect_equal(res_dgr, res_dgc, tolerance = 1e-7)
})


test_that("perf_optim$mu_dropout_thresh doesn't change anything if set to value smaller than minimum of Mu", {
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol = 10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE), cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- ~ group + cont1 + cont2

  set.seed(1)
  res1 <- glm_gp(Y, design = design, col_data = annot, perf_optim = list(offset_as_vec = TRUE))
  expect_true(all(res1[["Mu"]] > 0))

  set.seed(1)
  res2 <- glm_gp(Y, design = design, col_data = annot, perf_optim = list(offset_as_vec = TRUE, mu_dropout_thresh = min(res1[["Mu"]]) - .Machine[["double.eps"]]))
  res2[["Mu"]] <- as.matrix(res2[["Mu"]])
  attr(res2, "perf_optim") <- attr(res1, "perf_optim")

  expect_equal(res2[["Mu"]] > 0, res1[["Mu"]] > 0, tolerance = 0)
  expect_equal(res1, res2, tolerance = 1e-8)
})

test_that("perf_optim$do_parallel doesn't change anything", {
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol = 10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE), cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- ~ group + cont1 + cont2

  set.seed(1)
  res1 <- glm_gp(Y, design = design, col_data = annot)
  lapply(
    c(1L, 2L, 4L, 16L, 32L),
    function(p) {
      set.seed(1)
      res2 <- glm_gp(Y, design = design, col_data = annot, perf_optim = list(do_parallel = p))
      attr(res2, "perf_optim") <- attr(res1, "perf_optim")
      expect_equal(res1, res2, tolerance = 0)
    }
  )
})

test_that("perf_optim$delay_mu doesn't change anything", {
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol = 10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE), cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- ~ group + cont1 + cont2

  set.seed(1)
  res1 <- glm_gp(Y, design = design, col_data = annot)

  set.seed(1)
  res2 <- glm_gp(Y, design = design, col_data = annot, perf_optim = list(delay_mu = TRUE))

  res2[["Mu"]] <- res2[["Mu"]]()
  attr(res2, "perf_optim") <- attr(res1, "perf_optim")

  expect_equal(res1, res2, tolerance = 1e-7)
})

test_that("perf_optim$delay_mu doesn't change anything (global overdispersion)", {
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol = 10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE), cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- ~ group + cont1 + cont2

  set.seed(1)
  res1 <- glm_gp(Y, design = design, col_data = annot, overdispersion = "global")

  set.seed(1)
  res2 <- glm_gp(Y, design = design, col_data = annot, overdispersion = "global", perf_optim = list(delay_mu = TRUE))

  res2[["Mu"]] <- res2[["Mu"]]()
  attr(res2, "perf_optim") <- attr(res1, "perf_optim")

  expect_equal(res1, res2, tolerance = 1e-8)
})

test_that("perf_optim$use_nr_overdisp_impl doesn't significantly change resulting Beta values (w/ cox-reid adjustment)", {
  skip("known issue, default overdispersion estimator often ends up falling back to version w/o cox-reid adjustment, which is not true for C++ NR-based estimator")

  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol = 10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE), cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- ~ group + cont1 + cont2

  set.seed(1)
  res1 <- glm_gp(Y, design = design, col_data = annot)

  set.seed(1)
  res2 <- glm_gp(Y, design = design, col_data = annot, perf_optim = list(use_nr_overdisp_impl = TRUE))

  expect_equal(res1[["Beta"]], res2[["Beta"]], tolerance = 1e-8)
})

test_that("perf_optim$use_nr_overdisp_impl doesn't significantly change resulting Beta values (w/o cox-reid adjustment)", {
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol = 10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE), cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- ~ group + cont1 + cont2

  set.seed(1)
  res1 <- glm_gp(Y, design = design, col_data = annot, do_cox_reid_adjustment = FALSE)

  set.seed(1)
  res2 <- glm_gp(Y, design = design, col_data = annot, do_cox_reid_adjustment = FALSE, perf_optim = list(use_nr_overdisp_impl = TRUE))

  expect_equal(res1[["Beta"]], res2[["Beta"]], tolerance = 1e-8)
})