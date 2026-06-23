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
  function(obs = NULL, feats = NULL) {
    if (!is.null(feats)) {
      stopifnot(is.integer(feats) && all(feats > 0 & feats <= nrow(Beta)))
      Beta <- Beta[feats, ]
      if (!is.vector(offset)) {
        offset <- offset[feats, ]
      }
    }
    if (!is.null(obs)) {
      stopifnot(is.integer(obs) && all(obs > 0 & obs <= nrow(model_matrix)))
      model_matrix <- model_matrix[obs, ]
      if (is.vector(offset)) {
        offset <- offset[obs]
      } else {
        offset <- offset[, obs]
      }
    }
    calculate_mu(Beta, model_matrix, offset)
  }
}
