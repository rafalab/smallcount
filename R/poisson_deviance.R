#' Poisson Deviance
#'
#' @param y dgCMatrix sparse Matrix.
#' @param rate Row-wise rates.
#' @param n Total counts in each column.
#'
#' @import methods
#' @importFrom Matrix colSums rowSums
#' @export
#'
#' @examples
#' data("tenx_subset")
#' dev <- poisson_deviance(tenx_subset)
#' hist(dev, nclass = 50)
poisson_deviance <- function(y, rate = NULL, n = NULL){
  if(!is(y, "dgCMatrix")) stop("y must be class dgCMatrix")

  if(is.null(n)) n <- Matrix::colSums(y)
  if(is.null(rate)) rate <- Matrix::rowSums(y)/sum(n)

  ## The number of non-zeros in each column
  nnz <- diff(y@p)
  ## change y to contain the quantities needed to compute deviance
  mu_hat <- rep(n, nnz) * rate[y@i+1]
  y@x <- y@x * log(y@x / mu_hat) - (y@x - mu_hat)

  s <- vector("numeric", nrow(y))
  for(i in 1:ncol(y)){
    if(nnz[i] > 0){
      ind <- y@i[(y@p[i]+1):y@p[i+1]] + 1
      s[-ind] <- s[-ind] + rate[-ind]*n[i]
    } else{
      s <- s + rate*n[i]
    }
  }

  return(2 * Matrix::rowSums(y) + 2 * s)
}

#' Poisson Deviance
#'
#' @param y Sparse matrix (can be a base matrix, dgCMatrix, or SparseMatrix).
#' @param rate Row-wise rates.
#' @param n Total counts in each column.
#'
#' @importFrom SparseArray rowSums nzwhich
#' @export
poisson_deviance_new <- function(y, rate = NULL, n = NULL) {
  y <- convert_to_sparse(y)
  n <- colsums_with_default(y, n)
  rate <- row_rates_with_default(y, rate)

  nz_ind <- nzwhich(y, arr.ind = TRUE)
  y[nz_ind] <- y[nz_ind] * log(y[nz_ind] / (rate[nz_ind[, 1]] * n[nz_ind[, 2]]))
  2 * rowSums(y)
}
