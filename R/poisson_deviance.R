#' Poisson Deviance
#'
#' @param y Sparse matrix (can be a matrix, dgCMatrix, or SparseMatrix)
#' @param rate Row-wise rates
#' @param n Total counts in each column
#' 
#' @return Row-wise Poisson deviance
#'
#' @importFrom SparseArray rowSums nzwhich
#'
#' @examples
#' data("tenx_subset")
#' dev <- poissonDeviance(tenx_subset)
#' hist(dev, nclass = 50)
#' @export
poissonDeviance <- function(y, rate = NULL, n = NULL) {
    y <- .convertToSparse(y)
    n <- .colsumsWithDefault(y, n)
    rate <- .rowRatesWithDefault(y, rate)

    nz_ind <- nzwhich(y)
    mu <- .calculateMu(y, nz_ind, rate, n)
    # y[nz_ind] <- y[nz_ind] * log(y[nz_ind] / mu)
    y@SVT <- cppPoissonDevianceTransformation(y@SVT, mu)
    y@type <- "double"
    2 * rowSums(y)
}
