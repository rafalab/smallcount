#' Principal component analysis on residuals
#'
#' @param y Sparse matrix (can be a base matrix, dgCMatrix, or SparseMatrix).
#' @param k Number of principal components to return. Default 50.
#' @param rate Row-wise rates.
#' @param n Total counts in each column.
#' @param residual Residual type ("raw" or "pearson"). Default raw.
#'
#' @examples
#' data(tenx_subset)
#' dim(tenx_subset)
#' system.time({
#'   pc <- pca_poisson_residuals(tenx_subset, k = 10, residual = "pearson")
#' })
#' plot(pc$x[,1:2])
#' barplot(pc$sdev)
#' @export
pca_poisson_residuals <- function(y, k = 50, rate = NULL, n = NULL,
                                  residual = "raw") {
  y <- convert_to_sparse(y)
  n <- colsums_with_default(y, n)
  rate <- row_rates_with_default(y, rate)

  if (residual == "raw") {
    pca_poisson_raw_residuals(y, k, rate, n)
  } else if (residual == "pearson") {
    pca_poisson_pearson_residuals(y, k, rate, n)
  } else {
    stop(paste0('Invalid residual type: "', residual,
                '". Only "raw" and "pearson" are supported.'))
  }
}

#' Principal component analysis on raw residuals
#'
#' @inheritParams pca_poisson_residuals
#' @param k Number of principal components to return.
pca_poisson_raw_residuals <- function(y, k, rate, n) {
  scaled_y2 <- tcrossprod(y)
  yu <- (y %*% as.matrix(n)) %*% rate

  # Cross product of residuals
  rtr <- scaled_y2 - yu - t(yu) + sum(n ^ 2) * outer(rate, rate)
  compute_pca(rtr, k, y, rate, n)
}

#' Principal component analysis on Pearson residuals
#'
#' @inheritParams pca_poisson_raw_residuals
#' 
#' @importFrom SparseArray nzwhich
pca_poisson_pearson_residuals <- function(y, k, rate, n) {
  sqrt_rate <- sqrt(rate)
  nz_ind <- nzwhich(y, arr.ind = TRUE)
  y[nz_ind] <- y[nz_ind] / (sqrt_rate[nz_ind[, 1]] * sqrt(n)[nz_ind[, 2]])
  scaled_y2 <- tcrossprod(y)

  # Cross product of residuals
  rtr <- scaled_y2 - sum(n) * outer(sqrt_rate, sqrt_rate)
  compute_pca(rtr, k, y, sqrt_rate, sqrt(n))
}

#' Compute the eigendecomposition of r \%*\% t(r) for a residual matrix r
#'
#' @param rtr Cross product of residual matrix r = y - offset1 \%*\% offset2.
#' @param k Number of principal components to return.
#' @param y Sparse count matrix, potentially scaled.
#' @param offset1,offset2 Vectors whose product is the difference between y and
#'   the residual matrix.
#'
#' @import RSpectra
compute_pca <- function(rtr, k, y, offset1 = NULL, offset2 = NULL) {
  e <- RSpectra::eigs_sym(rtr, k = k)
  x <- crossprod(y, e$vectors)
  if (!is.null(offset1) && !is.null(offset2)) {
    x <- x - (offset2 %*% crossprod(offset1, e$vectors))
  }
  list(sdev = sqrt(e$values / (ncol(y) - 1)), rotation = e$vectors, x = x)
}
