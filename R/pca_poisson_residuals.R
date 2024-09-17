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
#' @param y dgCMatrix sparse Matrix.
#' @param k Number of principal components to return.
#' @param rate Row-wise rates.
#' @param n Total counts in each column.
#'
#' @importFrom Matrix colSums rowSums tcrossprod
#'
#' @examples
#' data(tenx_subset)
#' dim(tenx_subset)
#' system.time({pc <- pca_poisson_residuals(tenx_subset, k = 10)})
#' plot(pc$x[,1:2])
#' barplot(pc$sdev)
#' @export
pca_poisson_residuals <- function(y, k = 50, rate = NULL, n = NULL){
  if (!is(y, "dgCMatrix")) stop("y must be class dgCMatrix")

  if (is.null(n)) n <- Matrix::colSums(y)
  total <- sum(n)
  if (is.null(rate)) rate <- Matrix::rowSums(y)
  ind <- which(rate>0L) ##before dividing to obtain rate, keep index of all 0 rows
  rate <- rate / total

  ## rep(n, diff(y@p)) * rate[y@i+1] is mu_hat
  y@x <- y@x / sqrt(rep(n, diff(y@p)) * rate[y@i+1])
  y2 <- Matrix::tcrossprod(y) ## this is slowest step, but no way around it i don't think

  uy <- matrix(0, length(rate), length(rate))
  uy[ind, ind] <- outer(1/sqrt(rate[ind]), sqrt(rate[ind])) * rate[ind]*total

  u2 <- total * outer(sqrt(rate), sqrt(rate))

  ## cross product of residuals
  rtr <- as.matrix(y2) - 2*uy + u2
  compute_pca(rtr, y, rate, k)
}

#' Principal Component Analysis on Pearson residuals
#'
#' @param y SparseMatrix object.
#' @param k Number of principal components to return.
#' @param rate Row-wise rates.
#' @param n Total counts in each column.
#'
#' @importFrom SparseArray colSums rowSums
#' 
#' @examples
#' data(tenx_subset_new)
#' dim(tenx_subset_new)
#' system.time({pc <- pca_poisson_residuals_new(tenx_subset_new, k = 10)})
#' plot(pc$x[,1:2])
#' barplot(pc$sdev)
#' @export
pca_poisson_residuals_new <- function(y, k = 50, rate = NULL, n = NULL) {
  if (!is(y, "SparseMatrix")) stop("y must be class SparseMatrix")

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
