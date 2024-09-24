#' Compute Poisson Dispersion
#'
#' @param y dgCMatrix sparse Matrix.
#' @param rate Row-wise rates.
#' @param n Total counts in each column.
#'
#' @import methods
#' @importFrom Matrix rowSums colSums
#' @export
poisson_dispersion <- function(y, rate = NULL, n = NULL){
  if(!is(y, "dgCMatrix")) stop("y must be class dgCMatrix")

  if(is.null(n)) n <- Matrix::colSums(y)
  if(is.null(rate)) rate <- Matrix::rowSums(y)/sum(n)
  ## change y to contain the quantities needed to compute overdispersion
  mu_hat <- rep(n, diff(y@p)) * rate[y@i+1]
  y@x <- (y@x - mu_hat)^2 / mu_hat

  return(Matrix::rowSums(y)/(ncol(y)-1))
}

#' Compute Poisson Dispersion
#'
#' @param y Sparse matrix (can be a base matrix, dgCMatrix, or SparseMatrix).
#' @param rate Row-wise rates.
#' @param n Total counts in each column.
#'
#' @importFrom SparseArray nzwhich rowSums
#' @export
poisson_dispersion_new <- function(y, rate = NULL, n = NULL) {
  y <- convert_to_sparse(y)
  n <- colsums_with_default(y, n)
  rate <- row_rates_with_default(y, rate)

  nz_ind <- nzwhich(y, arr.ind = TRUE)
  y[nz_ind] <- y[nz_ind]^2 / (rate[nz_ind[, 1]] * n[nz_ind[, 2]])
  (rowSums(y) - sum(n) * rate) / (ncol(y) - 1)
}
