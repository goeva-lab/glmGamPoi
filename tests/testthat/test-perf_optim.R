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
})

test_that("perf_optim$do_parallel doesn't change anything", {
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol = 10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE), cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- ~ group + cont1 + cont2

  set.seed(1)
  res1 <- glm_gp(Y, design = design, col_data = annot)
  for (p in c(1L, 2L, 4L, 16L, 32L)) {
    set.seed(1)
    res2 <- glm_gp(Y, design = design, col_data = annot, perf_optim = list(do_parallel = p))
    attr(res2, "perf_optim") <- attr(res1, "perf_optim")
    expect_equal(res1, res2, tolerance = 0)
  }
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

  expect_equal(res1, res2, tolerance = 1e-7)
})

local({
  Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol = 10)
  annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE), cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
  design <- ~ group + cont1 + cont2

  for (cr_adj in c(FALSE, TRUE)) {
    test_that(sprintf("perf_optim$use_nr_overdisp_impl doesn't significantly change resulting Beta values (%s cox-reid adjustment)", if (cr_adj) "w/" else "w/o"), {
      set.seed(1)
      res1 <- glm_gp(Y, design = design, col_data = annot, do_cox_reid_adjustment = FALSE)

      set.seed(1)
      res2 <- glm_gp(Y, design = design, col_data = annot, do_cox_reid_adjustment = FALSE, perf_optim = list(use_nr_overdisp_impl = TRUE))

      logliks <- function(res) vapply(seq_len(nrow(Y)), function(i) conventional_loglikelihood_fast(Y[i, ], res[["Mu"]][i, ], log(res[["overdispersions"]][[i]]), res[["model_matrix"]], cr_adj), numeric(1))

      # assert new estimation method yields about as good or better estimates of overdispersion (w/ a tolerance of 1e-9 towards worsening)
      expect_all_true(logliks(res1) <= (logliks(res2) + 1e-8))

      skip_if(cr_adj, "known issue, default overdispersion estimator often ends up falling back to version w/o cox-reid adjustment, which is not true for C++ NR-based estimator")
      expect_equal(res1[["Beta"]], res2[["Beta"]], tolerance = 1e-8)
    })
  }
})

local({
  y <- matrix(rnbinom(n = 1e4 * 4, mu = 5, size = 1 / 0.7), nrow = 4)
  annot <- data.frame(cont1 = runif(1e4))
  mm <- handle_design_parameter(~cont1, y, annot, NULL)$model_matrix
  off <- combine_size_factors_and_offset(0, "normed_sum", y)$offset_matrix
  beta <- estimate_betas_roughly(y, mm, off)

  test_w_mu <- function(mu) {
    r1 <- overdispersion_mle(y, mu, mm, use_nr_overdisp_impl = TRUE)

    set.seed(1)
    r2 <- overdispersion_mle(y, mu, mm, subsample = 1000, use_nr_overdisp_impl = TRUE)
    # check that shuffle respects seeding
    set.seed(1)
    r3 <- overdispersion_mle(y, mu, mm, subsample = 1000, use_nr_overdisp_impl = TRUE)

    expect_equal(r1$estimate, r2$estimate, tolerance = 0.1)
    expect_equal(r2$estimate, r3$estimate, tolerance = 0)
  }

  test_that("use_nr_overdisp_impl w/ n_subsamples parameter works", {
    test_w_mu(calculate_mu(beta, mm, off))
  })
  test_that("use_nr_overdisp_impl w/ n_subsamples parameter works (w/ delayed Mu)", {
    test_w_mu(list(offset_matrix = off, beta_mat = beta))
  })
})


for (fn in c("predict", "residuals")) {
  local({
    Y <- matrix(rnbinom(n = 30 * 10, mu = 4, size = 0.3), nrow = 30, ncol = 10)
    rownames(Y) <- sprintf("gene_%s", seq_len(nrow(Y)))
    colnames(Y) <- sprintf("cell_%s", seq_len(ncol(Y)))

    annot <- data.frame(group = sample(c("A", "B"), size = 10, replace = TRUE), cont1 = rnorm(10), cont2 = rnorm(10, mean = 30))
    design <- ~ group + cont1 + cont2

    res <- glm_gp(Y, design = design, col_data = annot)
    fn.type <- get(sprintf("%s.%s", fn, class(res)))

    for (pred.type in as.list(eval(formals(fn.type)[["type"]]))) {
      test_that(sprintf("%s w/ type=%s (feat|obs).sub arguments works", fn, pred.type), {
        skip_if(pred.type == "randomized_quantile", "subsetting w/ randomized_quantile residual type is unsupported")

        fn.w <- function(...) fn.type(..., type = pred.type)
        ref <- fn.w(res)

        # using expect_all_equal instead of expect_equal because otherwise this generates
        # an obnoxious "number" of checks that pollutes reporter metrics
        expect_all_equal(
          vapply(
            list(
              fn.w(res, feat.sub = seq_len(nrow(Y))),
              fn.w(res, obs.sub = seq_len(ncol(Y))),
              fn.w(res, feat.sub = rownames(Y)),
              Reduce(rbind, lapply(seq_len(nrow(Y)), function(i) fn.w(res, feat.sub = i))),
              Reduce(cbind, lapply(seq_len(ncol(Y)), function(j) fn.w(res, obs.sub = j))),
              Reduce(rbind, lapply(rownames(Y), function(i) fn.w(res, feat.sub = i))),
              Reduce(cbind, lapply(colnames(Y), function(j) fn.w(res, obs.sub = j)))
            ),
            function(e) {
              cmp <- waldo::compare(ref, e, tolerance = 0)
              if (length(cmp) == 0) {
                ""
              } else {
                cmp
              }
            },
            character(1)
          ),
          ""
        )
      })
    }
  })
}
