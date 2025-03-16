#' Poisson Dispersion
#'
#' @param y Sparse matrix (can be a matrix, dgCMatrix, or SparseMatrix).
#' @param rate Row-wise rates.
#' @param n Total counts in each column.
#'
#' @importFrom SparseArray nzwhich rowSums
#'
#' @examples
#' data("tenx_subset")
#' disp <- poissonDispersion(tenx_subset)
#' hist(disp, nclass = 50)
#' @export
poissonDispersion <- function(y, rate = NULL, n = NULL) {
  y <- convert_to_sparse(y)
  n <- colsums_with_default(y, n)
  rate <- row_rates_with_default(y, rate)

  nz_ind <- nzwhich(y)
  mu <- calculate_mu(y, nz_ind, rate, n)
  y[nz_ind] <- y[nz_ind]^2 / mu
  (rowSums(y) - sum(n) * rate) / (ncol(y) - 1)
}
