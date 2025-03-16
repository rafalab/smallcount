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

#' Get the row indices of non-zero matrix values
#'
#' @param y SparseMatrix object
#' @param nz_ind Linear indices of non-zero values
#' @noRd
get_nz_rows <- function(y, nz_ind) {
  (nz_ind - 1) %% nrow(y) + 1
}

#' Get the column indices of non-zero matrix values
#'
#' @param y SparseMatrix object
#' @param nz_ind Linear indices of non-zero values
#' @noRd
get_nz_cols <- function(y, nz_ind) {
  (nz_ind - 1) %/% nrow(y) + 1
}

#' Calculate mu values from row rates and column sums
#'
#' @param y SparseMatrix object
#' @param rate Vector of row rates
#' @param n Vector of column sums
#' @return Vector of mu values for non-zero entries
#' @noRd
calculate_mu <- function(y, nz_ind, rate, n) {
  nz_rows <- get_nz_rows(y, nz_ind)
  nz_cols <- get_nz_cols(y, nz_ind)
  rate[nz_rows] * n[nz_cols]
}
