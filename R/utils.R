#' Convert a sparse matrix to a SparseMatrix object
#'
#' @param y Sparse matrix (can be a matrix, dgCMatrix, or SparseMatrix)
#' 
#' @return SparseMatrix object
#'
#' @importFrom methods as is
#' @keywords internal
.convertToSparse <- function(y) {
    if (is(y, "matrix") || is(y, "dgCMatrix")) {
        y <- as(y, "SparseMatrix")
    } else if (!is(y, "SparseMatrix")) {
        stop("y must be a matrix, dgCMatrix, or SparseMatrix")
    }
    return(y)
}

#' Return a non-null default value or compute column sums
#'
#' @param y SparseMatrix object
#' @param default Default column sums
#' 
#' @return Column sums, or default value if provided
#'
#' @importFrom SparseArray colSums
#' @keywords internal
.colsumsWithDefault <- function(y, default) {
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
#' @return Row-wise rates (row sums / total sum), or default value if provided
#'
#' @importFrom SparseArray colSums
#' @keywords internal
.rowRatesWithDefault <- function(y, default) {
    if (!is.null(default)) {
        return(default)
    }
    rsums <- rowSums(y)
    rsums / sum(rsums)
}

# Divide a by b, returning zero if b is zero
.safeDivide <- function(a, b) {
    ifelse(b == 0, 0, a / b)
}

#' Get the row indices of non-zero matrix values
#'
#' @param y SparseMatrix object
#' @param nz_ind Linear indices of non-zero values
#' 
#' @return Row indices of non-zero matrix values
#'
#' @keywords internal
.nzrows <- function(y, nz_ind) {
    (nz_ind - 1) %% nrow(y) + 1
}

#' Get the column indices of non-zero matrix values
#'
#' @param y SparseMatrix object
#' @param nz_ind Linear indices of non-zero values
#' 
#' @return Column indices of non-zero matrix values
#'
#' @keywords internal
.nzcols <- function(y, nz_ind) {
    (nz_ind - 1) %/% nrow(y) + 1
}

#' Calculate mu values from row rates and column sums
#'
#' @param y SparseMatrix object
#' @param rate Vector of row rates
#' @param n Vector of column sums
#' 
#' @return Expected values for non-zero matrix entries
#'
#' @keywords internal
.calculateMu <- function(y, nz_ind, rate, n) {
    nz_rows <- .nzrows(y, nz_ind)
    nz_cols <- .nzcols(y, nz_ind)
    rate[nz_rows] * n[nz_cols]
}
