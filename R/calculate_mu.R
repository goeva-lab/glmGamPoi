#' @importFrom bit64 integer64
calculate_mu <- local({
  mu_col_w_dropout_ovec <- compiler::cmpfun(
    function(j, model_matrix.t, Beta, offset_vec, dropout_thresh) {
      col <- exp((Beta %*% model_matrix.t[, j]) + offset_vec[[j]])
      mask <- col >= dropout_thresh
      list(col[mask], which(mask), sum(mask))
    },
    options = list(optimize = 3)
  )
  mu_col_w_dropout_oconst <- compiler::cmpfun(
    function(j, model_matrix.t, Beta, offset_vec, dropout_thresh) {
      col <- exp((Beta %*% model_matrix.t[, j]) + offset_vec)
      mask <- col >= dropout_thresh
      list(col[mask], which(mask), sum(mask))
    },
    options = list(optimize = 3)
  )
  mu_row_w_dropout <- compiler::cmpfun(
    function(i, model_matrix, Beta.t, offset_vec, dropout_thresh) {
      row <- exp(Matrix::tcrossprod(Beta.t[, i], model_matrix) + offset_vec)
      mask <- row >= dropout_thresh
      list(row[mask], which(mask), sum(mask))
    },
    options = list(optimize = 3)
  )

  function(Beta, model_matrix, offset_matrix, dropout_thresh = 0, dropout_by_row = FALSE) {
    make_offset_hdf5_mat <- is(offset_matrix, "DelayedMatrix")
    if (make_offset_hdf5_mat) {
      mu <- exp(delayed_matrix_multiply(DelayedArray::DelayedArray(Beta), DelayedArray::DelayedArray(t(model_matrix))) + offset_matrix)
      mu <- HDF5Array::writeHDF5Array(mu)
      mu
    } else {
      if (is.matrix(offset_matrix)) {
        exp(Matrix::tcrossprod(Beta, model_matrix) + offset_matrix)
      } else {
        stopifnot(is.vector(offset_matrix, mode = "numeric"))

        if (dropout_thresh > 0) {
          if (dropout_by_row) {
            out <- lapply(
              seq_len(nrow(Beta)),
              if ((length(offset_matrix) == 1) || (length(offset_matrix) == nrow(model_matrix))) {
                mu_row_w_dropout
              } else {
                stop("wrong size offset vector")
              },
              model_matrix = model_matrix,
              Beta.t = t(Beta),
              offset_vec = offset_matrix,
              dropout_thresh = dropout_thresh
            )

            Matrix::sparseMatrix(
              x = unlist(lapply(out, function(e) e[[1]]), recursive = FALSE),
              j = unlist(lapply(out, function(e) e[[2]]), recursive = FALSE),
              p = bit64::c.integer64(0, cumsum(bit64::as.integer64(unlist(lapply(out, function(e) e[[3]]), recursive = FALSE)))),
              dims = c(nrow(Beta), length(offset_matrix)),
              repr = "R"
            )
          } else {
            out <- lapply(
              seq_len(nrow(model_matrix)),
              if (length(offset_matrix) == 1) {
                mu_col_w_dropout_oconst
              } else if (length(offset_matrix) == nrow(model_matrix)) {
                mu_col_w_dropout_ovec
              } else {
                stop("wrong size offset vector")
              },
              model_matrix.t = t(model_matrix),
              Beta = Beta,
              offset_vec = offset_matrix,
              dropout_thresh = dropout_thresh
            )

            Matrix::sparseMatrix(
              x = unlist(lapply(out, function(e) e[[1]]), recursive = FALSE),
              i = unlist(lapply(out, function(e) e[[2]]), recursive = FALSE),
              p = bit64::c.integer64(0, cumsum(bit64::as.integer64(unlist(lapply(out, function(e) e[[3]]), recursive = FALSE)))),
              dims = c(nrow(Beta), length(offset_matrix)),
              repr = "C"
            )
          }
        } else {
          exp(add_vector_to_each_row(Matrix::tcrossprod(Beta, model_matrix), offset_matrix))
        }
      }
    }
  }
})

mk_delayed_mu <- function(Beta, model_matrix, offset, dropout_thresh, dropout_by_row) {
  function(obs = NULL, feats = NULL) {
    if(!is.null(feats)){
      stopifnot(is.integer(feats) && all(feats > 0 & feats <= nrow(model_matrix)))
      Beta <- Beta[feats, ]
      if(!is.vector(offset)){
        offset <- offset[feats, ]
      }
    }
    if(!is.null(obs)){
      stopifnot(is.integer(obs) && all(obs > 0 & obs <= nrow(model_matrix)))
      model_matrix <- model_matrix[obs, ]
      if(is.vector(offset)){
        offset <- offset[obs]
      }else{
        offset <- offset[, obs]
      }
    }
    calculate_mu(Beta, model_matrix, offset, dropout_thresh, dropout_by_row)
  }
}