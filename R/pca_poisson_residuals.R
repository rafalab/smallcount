#' Principal Component Analysis on Pearson residuals
#'
#' @param y A tgCMatrix sparse Matrix.
#' @param k Number of principal components to return.
#' @param rate The rowwise rates.
#' @param n The total counts in each column.
#'
#' @importFrom Matrix colSums rowSums crossprod
#' @export
#'
#' @example
#' data(tenx_subset)
#' dim(tenx_subset)
#' system.time({pc <- pca_poisson_residuals(tenx_subset, k = 10)})
#' plot(pc$x[,1:2])
#' barplot(pc$sdev)

pca_poisson_residuals <- function(y, k = 50, rate = NULL, n = NULL){

  if(!is(y, "dgCMatrix")) stop("y must be class dgCMatrix")

  if(is.null(n)) n <- Matrix::colSums(y)
  total <- sum(n)
  if(is.null(rate)) rate <- Matrix::rowSums(y)
  ind <- which(rate>0L) ##before dividing to obtain rate, keep index of all 0 rows
  rate <- rate / total

  ## rep(n, diff(y@p)) * rate[y@i+1] is mu_hat
  y@x <- y@x / sqrt(rep(n, diff(y@p)) * rate[y@i+1])
  y2 <-Matrix::tcrossprod(y) ## this is slowest step, but no way around it i don't think

  uy <- matrix(0, length(rate), length(rate))
  uy[ind, ind] <- outer(1/sqrt(rate[ind]), sqrt(rate[ind])) * rate[ind]*total

  u2 <- total * outer(sqrt(rate), sqrt(rate))

  ## cross product of residuals
  rtr <- as.matrix(y2) - 2*uy + u2
  ##
  e <- RSpectra::eigs_sym(rtr, k = k)


  return(list(sdev = sqrt(e$values/(ncol(y)-1)),
              rotation = e$vectors,
              x = t(sweep(as.matrix(Matrix::crossprod(e$vectors, y)),
                        1,
                        colSums(e$vectors* sqrt(rate)), "-"))))
}
