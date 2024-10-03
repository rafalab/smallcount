#' Convert a sparse matrix to a SparseMatrix object
#'
#' @param y Sparse matrix (can be a base matrix, dgCMatrix, or SparseMatrix).
#'
#' @import methods
#' @noRd
convert_to_sparse <- function(y) {
  if (is(y, "matrix") || is(y, "dgCMatrix")) {
    y <- as(y, "SparseMatrix")
  } else if (!is(y, "SparseMatrix")) {
    stop("y must be a base matrix, dgCMatrix, or SparseMatrix")
  }
  return(y)
}

#' Return a non-null default value or compute column sums
#'
#' @param y SparseMatrix object
#' @param default Default column sums
#'
#' @importFrom SparseArray colSums
#' @noRd
colsums_with_default <- function(y, default) {
  if (!is.null(default)) {
    return(default)
  }
  colSums(y)
}

#' Return a non-null default value or compute row-wise rates
#'
#' @param y SparseMatrix object
#' @param default Default row-wise rates
#'
#' @importFrom SparseArray colSums
#' @noRd
row_rates_with_default <- function(y, default) {
  if (!is.null(default)) {
    return(default)
  }
  rsums <- rowSums(y)
  rsums / sum(rsums)
}

#' Divide a by b, returning zero if b is zero
#'
#' @param a Numerator
#' @param b Denominator
#' @noRd
safe_divide <- function(a, b) {
  ifelse(b == 0, 0, a / b)
}
