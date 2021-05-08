poisson_deviance <- function(y, rate = NULL, n = NULL){

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
