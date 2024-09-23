#' Convert a sparse matrix to a SparseMatrix object
#'
#' @param y Sparse matrix (can be a base matrix, dgCMatrix, or SparseMatrix).
convert_to_sparse <- function(y) {
  if (is(y, "matrix") || is(y, "dgCMatrix")) {
    y <- as(y, "SparseMatrix")
  } else if (!is(y, "SparseMatrix")) {
    stop("y must be a base matrix, dgCMatrix, or SparseMatrix")
  }
  return(y)
}

#' Divide a by b, resulting in zero if b is zero
#'
#' @param a Numerator
#' @param b Denominator
safe_divide <- function(a, b) {
  ifelse(b == 0, 0, a / b)
}
