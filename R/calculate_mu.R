calculate_mu <- function(Beta, model_matrix, offset_matrix) {
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
      exp(add_vector_to_each_row(Matrix::tcrossprod(Beta, model_matrix), offset_matrix))
    }
  }
}

mk_delayed_mu <- function(Beta, model_matrix, offset) {
  function(obs = NULL, feat = NULL) {
    if (!is.null(feat)) {
      feat <- handle_sub_param(nrow(Beta), feat)
      Beta <- Beta[feat, , drop = FALSE]
      if (!is.vector(offset)) { offset <- offset[feat, , drop = FALSE] }
    }
    if (!is.null(obs)) {
      obs <- handle_sub_param(nrow(model_matrix), obs)
      model_matrix <- model_matrix[obs, , drop = FALSE]
      offset <- if (is.vector(offset)) { offset[obs] } else { offset[, obs, drop = FALSE] }
    }
    calculate_mu(Beta, model_matrix, offset)
  }
}
