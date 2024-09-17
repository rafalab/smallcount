#' Computes PCA for a residual matrix from its cross product
#'
#' @param rtr Cross product of residual matrix.
#' @param y Original sparse matrix, scaled by inverse standard deviation.
#' @param rate Row-wise rates.
#' @param k Number of principal components to return.
#'
#' @import RSpectra
#' @importFrom Matrix crossprod
compute_pca <- function(rtr, y, rate, k) {
  e <- RSpectra::eigs_sym(rtr, k = k)
  x <- t(sweep(
    as.matrix(Matrix::crossprod(e$vectors, y)),
    1,
    colSums(e$vectors * sqrt(rate)),
    "-"
  ))
  list(sdev = sqrt(e$values / (ncol(y) - 1)), rotation = e$vectors, x = x)
}

#' Principal Component Analysis on Pearson residuals
#'
#' @param y Sparse matrix (can be a base matrix, dgCMatrix, or SparseMatrix).
#' @param k Number of principal components to return.
#' @param rate Row-wise rates.
#' @param n Total counts in each column.
#'
#' @importFrom SparseArray colSums rowSums
#' 
#' @examples
#' data(tenx_subset)
#' dim(tenx_subset)
#' system.time({pc <- pca_poisson_residuals(tenx_subset, k = 10)})
#' plot(pc$x[,1:2])
#' barplot(pc$sdev)
#' @export
pca_poisson_residuals <- function(y, k = 50, rate = NULL, n = NULL) {
  if (is(y, "matrix") || is(y, "dgCMatrix")) {
    y <- as(y, "SparseMatrix")
  } else if (!is(y, "SparseMatrix")) {
    stop("y must be class matrix, dgCMatrix, or SparseMatrix")
  }

  n <- if (is.null(n)) colSums(y) else n
  rate <- if (is.null(rate)) rowSums(y) else rate
  total <- sum(n)
  rate <- rate / total
  sqrt_rate <- sqrt(rate)

  nz_ind <- nzwhich(y, arr.ind = TRUE)
  y[nz_ind] <- y[nz_ind] / (sqrt_rate[nz_ind[, 1]] * sqrt(n)[nz_ind[, 2]])
  scaled_y2 <- tcrossprod(y)

  # Cross product of residuals
  rtr <- scaled_y2 - total * outer(sqrt_rate, sqrt_rate)
  compute_pca(rtr, y, rate, k)
}
