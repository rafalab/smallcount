#' Compute Poisson Dispersion
#'
#' @param y A tgCMatrix sparse Matrix.
#' @param rate The rowwise rates.
#' @param n The total counts in each column.
#'
#' @importFrom Matrix rowSums colSums
#' @export
#'
poisson_dispersion <- function(y, rate = NULL, n = NULL){

  if(!is(y, "dgCMatrix")) stop("y must be class dgCMatrix")

  if(is.null(n)) n <- Matrix::colSums(y)
  if(is.null(rate)) rate <- Matrix::rowSums(y)/sum(n)
  ## change y to contain the quantities needed to compute overdispersion
  mu_hat <- rep(n, diff(y@p)) * rate[y@i+1]
  y@x <- (y@x - mu_hat)^2 / mu_hat

  return(Matrix::rowSums(y)/(ncol(y)-1))
}

