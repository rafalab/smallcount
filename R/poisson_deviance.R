#' Poisson Deviance
#'
#' @param y Sparse matrix (can be a base matrix, dgCMatrix, or SparseMatrix).
#' @param rate Row-wise rates.
#' @param n Total counts in each column.
#'
#' @importFrom SparseArray rowSums nzwhich
#'
#' @examples
#' data("tenx_subset")
#' dev <- poissonDeviance(tenx_subset)
#' hist(dev, nclass = 50)
#' @export
poissonDeviance <- function(y, rate = NULL, n = NULL) {
  y <- convert_to_sparse(y)
  n <- colsums_with_default(y, n)
  rate <- row_rates_with_default(y, rate)

  nz_ind <- nzwhich(y)
  mu <- calculate_mu(y, nz_ind, rate, n)
  y[nz_ind] <- y[nz_ind] * log(y[nz_ind] / mu)
  2 * rowSums(y)
}
