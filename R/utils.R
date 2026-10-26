delayed_matrix_apply_block <- function(Y, Mu, overdispersion, FUN) {
  stopifnot(nrow(Y) == nrow(Mu), ncol(Y) == ncol(Mu))
  res_sink <- HDF5Array::HDF5RealizationSink(dim(Y))
  on.exit(
    {
      DelayedArray::close(res_sink)
    },
    add = TRUE
  )

  res_grid <- DelayedArray::defaultAutoGrid(res_sink)

  for (coord1 in seq_len(ncol(res_grid))) {
    for (coord2 in seq_len(nrow(res_grid))) {
      sel <- res_grid[[coord2, coord1]]
      Y_block <- DelayedArray::read_block(Y, sel)
      Mu_block <- DelayedArray::read_block(Mu, sel)
      res_block <- FUN(Y_block, Mu_block, overdispersion[seq(sel@ranges@width[1]) - 1 + sel@ranges@start[1]])
      DelayedArray::write_block(res_sink, res_grid[[coord2, coord1]], res_block)
    }
  }

  as(res_sink, "DelayedArray")
}


delayed_matrix_multiply <- function(x, y) {
  res_sink <- HDF5Array::HDF5RealizationSink(c(nrow(x), ncol(y)))
  on.exit(
    {
      DelayedArray::close(res_sink)
    },
    add = TRUE
  )

  res_grid <- DelayedArray::defaultAutoGrid(res_sink)

  row_ticks <- cumsum(vapply(
    seq_len(dim(res_grid)[1]),
    function(idx) {
      dim(res_grid[[idx, 1L]])[1]
    },
    FUN.VALUE = 0L
  ))
  col_ticks <- cumsum(vapply(
    seq_len(dim(res_grid)[2]),
    function(idx) {
      dim(res_grid[[1L, idx]])[2]
    },
    FUN.VALUE = 0L
  ))

  x_grid <- DelayedArray::ArbitraryArrayGrid(tickmarks = list(row_ticks, ncol(x)))
  y_grid <- DelayedArray::ArbitraryArrayGrid(tickmarks = list(nrow(y), col_ticks))

  for (coord1 in seq_len(ncol(res_grid))) {
    for (coord2 in seq_len(nrow(res_grid))) {
      x_block <- DelayedArray::read_block(x, x_grid[[coord2]])
      y_block <- DelayedArray::read_block(y, y_grid[[coord1]])
      res_block <- x_block %*% y_block
      DelayedArray::write_block(res_sink, res_grid[[coord2, coord1]], res_block)
    }
  }

  as(res_sink, "DelayedArray")
}


add_vector_to_each_column <- function(matrix, vector) {
  stopifnot(length(vector) == 1 || length(vector) == nrow(matrix))
  matrix + vector
}


add_vector_to_each_row <- function(matrix, vector) {
  stopifnot(length(vector) == 1 || length(vector) == ncol(matrix))

  # This is a work-around for https://github.com/Bioconductor/DelayedArray/issues/123
  if (is(matrix, "ConstantMatrix") && all(matrix == 0)) {
    return(matrix(vector, nrow = nrow(matrix), ncol = ncol(matrix), byrow = TRUE))
  }

  if (length(vector) == 1) {
    matrix + vector
  } else {
    t(t(matrix) + vector)
  }
}


subtract_vector_from_each_column <- function(matrix, vector) {
  stopifnot(length(vector) == 1 || length(vector) == nrow(matrix))
  matrix - vector
}


subtract_vector_from_each_row <- function(matrix, vector) {
  stopifnot(length(vector) == 1 || length(vector) == ncol(matrix))
  if (length(vector) == 1) {
    matrix - vector
  } else {
    t(t(matrix) - vector)
  }
}


multiply_vector_to_each_column <- function(matrix, vector) {
  stopifnot(length(vector) == 1 || length(vector) == nrow(matrix))
  matrix * vector
}


multiply_vector_to_each_row <- function(matrix, vector) {
  stopifnot(length(vector) == 1 || length(vector) == ncol(matrix))
  if (length(vector) == 1) {
    matrix * vector
  } else {
    t(t(matrix) * vector)
  }
}


divide_each_column_by_vector <- function(matrix, vector) {
  stopifnot(length(vector) == 1 || length(vector) == nrow(matrix))
  matrix / vector
}


divide_each_row_by_vector <- function(matrix, vector) {
  stopifnot(length(vector) == 1 || length(vector) == ncol(matrix))
  if (length(vector) == 1) {
    matrix / vector
  } else {
    t(t(matrix) / vector)
  }
}


#' Solve the equation `Y = A B` for `A` or `B`
#'
#' @param Y the left side of the equation
#' @param A,B the known matrix on the right side of the equation
#' @param w a vector with weights. If `NULL` it is ignored,
#'   otherwise it must be of length 1 or have the same length as
#'   columns in `Y`. Default: `NULL`
#'
#' @keywords internals
solve_lm_for_A <- function(Y, B, w = NULL) {
  if (nrow(B) == 0) {
    matrix(numeric(0), nrow = nrow(Y), ncol = 0)
  } else if (nrow(Y) == 0) {
    matrix(numeric(0), nrow = 0, ncol = nrow(B))
  } else if (is.null(w)) {
    qrx <- qr(t(B))
    Q <- qr.Q(qrx)
    R <- qr.R(qrx)

    t(backsolve(R, t(Y %*% Q)))
  } else {
    stopifnot(length(w) == 1 || length(w) == ncol(Y))
    sqrt_w <- sqrt(w)
    qrx <- qr(t(B) * sqrt_w)
    Q <- qr.Q(qrx)
    R <- qr.R(qrx)

    t(backsolve(R, t(t(t(Y) * sqrt_w) %*% Q)))
  }
}


#' @rdname solve_lm_for_A
solve_lm_for_B <- function(Y, A, w = NULL) {
  if (ncol(A) == 0) {
    matrix(numeric(0), nrow = 0, ncol = ncol(Y))
  } else if (ncol(Y) == 0) {
    matrix(numeric(0), nrow = ncol(A), ncol = 0)
  } else if (is.null(w)) {
    qrx <- qr(A)
    Q <- qr.Q(qrx)
    R <- qr.R(qrx)

    backsolve(R, t(Q) %*% Y)
  } else {
    stopifnot(length(w) == 1 || length(w) == nrow(Y))
    sqrt_w <- sqrt(w)
    qrx <- qr(A * sqrt_w)
    Q <- qr.Q(qrx)
    R <- qr.R(qrx)

    backsolve(R, t(Q) %*% (Y * sqrt_w))
  }
}

# sparsity aware element-wise division (e.g. hadamard/schur operation)
# performs specific optimizations to conserve sparsity/layout of output iff:
# - `lhs` is `dgCMatrix` or `dgRMatrix`
# - `rhs` is entirely non-zero
div_mtx_elemwise <- function(lhs, rhs) {
  stopifnot(dim(lhs) == dim(rhs)) # elementwise operation

  # division of zero by zero does not preserve
  # sparsity so to ensure correctness we skip this
  # optimization if zero values are present in the rhs
  if (all(rhs != 0)) {
    if (is(lhs, "dsparseMatrix")) {
      t.lhs <- as(lhs, "TsparseMatrix")
      lhs@x <- lhs@x / c(rhs)[t.lhs@i + (t.lhs@j * nrow(lhs)) + 1]
      return(lhs)
    }
  }
  lhs / rhs
}

# sparsity aware column-wise division (e.g. elements of `rhs` represent a divisor for each column of `lhs`)
# performs specific optimizations to conserve sparsity/layout of output iff:
# - `lhs` is `dgCMatrix` or `dgRMatrix`
# - `rhs` is entirely non-zero
div_mtx_colwise <- function(lhs, rhs) {
  stopifnot(ncol(lhs) == length(rhs)) # column-wise operation

  # division of zero by zero does not preserve
  # sparsity so to ensure correctness we skip this
  # optimization if zero values are present in the rhs
  if (all(rhs != 0)) {
    if (is(lhs, "dgCMatrix")) {
      lhs@x <- lhs@x / rep(rhs, diff(lhs@p))
      return(lhs)
    } else if (is(lhs, "dgRMatrix")) {
      lhs@x <- lhs@x / rhs[lhs@j + 1]
      return(lhs)
    }
  }

  t(t(lhs) / rhs)
}


handle_Mu_rowmeans <- local({
  vec_off_fn <- compiler::cmpfun(function(gene_idx, beta_t, mm, offs) mean(calculate_mu(matrix(beta_t[, gene_idx], nrow = 1), mm, offs)), options = list(optimize = 3L))
  mtx_off_fn <- compiler::cmpfun(function(gene_idx, beta_t, mm, offs_t) mean(calculate_mu(matrix(beta_t[, gene_idx], nrow = 1), mm, offs_t[, gene_idx])), options = list(optimize = 3L))

  function(Mu, row.names) {
    if (is.function(Mu)) {
      beta_t <- t(attr(Mu, "betas"))
      mm <- attr(Mu, "model_matrix")
      offs <- attr(Mu, "offsets")

      if (is.vector(offs)) {
        vapply(setNames(seq_len(ncol(beta_t)), nm = row.names), vec_off_fn, numeric(1), beta_t = beta_t, mm = mm, offs = offs)
      } else {
        offs_t <- t(offs)
        vapply(setNames(seq_len(ncol(beta_t)), nm = row.names), mtx_off_fn, numeric(1), beta_t = beta_t, mm = mm, offs_t = offs_t)
      }
    } else {
      DelayedMatrixStats::rowMeans2(Mu)
    }
  }
})


handle_sub_param <- function(check.against, sub.param) {
  if (is.numeric(sub.param) && !is.integer(sub.param)) {
    stop("indexing parameter must be an integer value")
  }

  if (is.integer(check.against) && (length(check.against) == 1)) {
    if (is.integer(sub.param)) {
      stopifnot(all((sub.param > 0) & (sub.param <= check.against)))
      return(sub.param)
    }

    stopifnot(is.logical(sub.param) && (length(sub.param) == check.against))
    return(which(sub.param, useNames = FALSE))
  }

  stopifnot(is.vector(check.against, mode = "character"))

  if (is.character(sub.param)) {
    sub.idx <- pmatch(sub.param, check.against, duplicates.ok = TRUE)
    stopifnot(!anyNA(sub.idx))

    return(sub.idx)
  }

  if (is.logical(sub.param)) {
    stopifnot(length(sub.param) == length(check.against))
    return(which(sub.param, useNames = FALSE))
  }

  stopifnot(is.integer(sub.param) && all(sub.param > 0 & sub.param <= check.against))
  sub.param
}
