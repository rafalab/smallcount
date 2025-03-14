#' Compute the eigendecomposition of `R %*% t(R)` for a residual matrix `R`
#'
#' @param rtr Cross product of residual matrix `R = y - offset1 %*% offset2`.
#' @param k Number of principal components to return.
#' @param y Sparse matrix.
#' @param offset1,offset2 Vectors whose product is the difference between `y`
#'   and the residual matrix.
#'
#' @import RSpectra
#' @noRd
compute_pca <- function(rtr, k, y, offset1 = NULL, offset2 = NULL) {
  e <- RSpectra::eigs_sym(rtr, k = k)
  x <- crossprod(y, e$vectors)
  if (!is.null(offset1) && !is.null(offset2)) {
    x <- x - (offset2 %*% crossprod(offset1, e$vectors))
  }
  list(sdev = sqrt(e$values / (ncol(y) - 1)), rotation = e$vectors, x = x)
}

#' Principal component analysis on raw residuals
#'
#' @inheritParams poissonPca
#' @param k Number of principal components to return.
#' @param row_offset,col_offset Vectors whose product subtracted from `y` gives
#'   the residual matrix.
#'
#' @noRd
raw_residuals_pca <- function(y, k, row_offset, col_offset) {
  y2 <- tcrossprod(y)
  yu <- (y %*% as.matrix(col_offset)) %*% row_offset
  u2 <- sum(col_offset ^ 2) * outer(row_offset, row_offset)
  
  # Cross product of residuals
  rtr <- y2 - yu - t(yu) + u2
  compute_pca(rtr, k, y, row_offset, col_offset)
}

#' Principal component analysis on Pearson residuals
#'
#' @inheritParams raw_residuals_pca
#'
#' @importFrom SparseArray nzwhich
#' @noRd
poisson_pearson_residuals_pca <- function(y, k) {
  n <- colSums(y)
  sqrt_rate <- sqrt(rowSums(y) / sum(n))
  nz_ind <- nzwhich(y, arr.ind = TRUE)
  y[nz_ind] <- y[nz_ind] / (sqrt_rate[nz_ind[, 1]] * sqrt(n)[nz_ind[, 2]])
  scaled_y2 <- tcrossprod(y)
  
  # Cross product of residuals
  rtr <- scaled_y2 - sum(n) * outer(sqrt_rate, sqrt_rate)
  compute_pca(rtr, k, y, sqrt_rate, sqrt(n))
}

#' Principal component analysis on deviance residuals
#'
#' @inheritParams raw_residuals_pca
#'
#' @importFrom SparseArray nzwhich
#' @noRd
poisson_deviance_residuals_pca <- function(y, k) {
  n <- colSums(y)
  rate <- rowSums(y) / sum(n)
  ys <- nzvals(y)
  nz_ind <- nzwhich(y, arr.ind = TRUE)
  mu <- rate[nz_ind[, 1]] * n[nz_ind[, 2]]
  
  deviance <- 2 * (ys * log(ys / mu) - ys + mu)
  y[nz_ind] <- sign(ys - mu) * sqrt(deviance) + sqrt(2 * mu)
  
  raw_residuals_pca(y, k, sqrt(2 * rate), sqrt(n))  
}

# Map of residual types to PCA functions
RESIDUAL_PCA <- list(pearson = poisson_pearson_residuals_pca,
                     deviance = poisson_deviance_residuals_pca)

# Converts a NULL value or a character string to a CountTransform object.
get_count_transform <- function(transform, center, scale, coef) {
  if (is.null(transform)) {
    return(identity_transform(center, scale))
  }
  switch(
    transform,
    id = identity_transform(center, scale),
    log1p = log1p_transform(center, scale),
    cpm_log1p = cpm_log1p_transform(center, scale),
    med_log1p = scaled_log1p_transform(coef, center, scale),
    NULL
  )
}

#' Principal Component Analysis on Poisson data
#'
#' @param y Sparse matrix (can be a base matrix, dgCMatrix, or SparseMatrix).
#' @param k Number of principal components to return. Default 50.
#' @param transform CountTransform object or character(1) specifying a
#'   transformation to apply to `y` before PCA. Arguments `center` and `scale`
#'   are only applied when `transform` is a character string for a non-residual
#'   transformation. Accepted strings are: `"pearson"` for Pearson residuals,
#'   `"deviance"` for deviance residuals, `"id"` (or `NULL`) for the identity
#'   function, `"log1p"` for `log(x + 1)`, `"cpm_log1p"` for `log(x/1e6 + 1)`,
#'   and `"med_log1p"` for `log(x/median(colSums(y)) + 1)`.
#' @inheritParams CountTransform
#'
#' @importFrom stats median
#' @importFrom SparseArray nzwhich
#'
#' @examples
#' data(tenx_subset)
#' dim(tenx_subset)
#' system.time({
#'   pc <- poissonPca(tenx_subset, k = 10, transform = "pearson")
#' })
#' plot(pc$x[,1:2])
#' barplot(pc$sdev)
#' @export
poissonPca <- function(y, k = 50, transform = NULL, center = FALSE,
                        scale = FALSE) {
  y <- convert_to_sparse(y)

  if (is.character(transform) && (transform %in% names(RESIDUAL_PCA))) {
    return(RESIDUAL_PCA[[transform]](y, k))
  } else if (is.character(transform) || is.null(transform)) {
    coef <- median(colSums(y))
    transform <- get_count_transform(transform, center, scale, coef)
  }

  if (!is(transform, "CountTransform")) {
    stop(paste0("Invalid transform: ", transform,
                ". See documentation for supported options."))
  }

  tmatrix <- TransformedMatrix(y, transform)
  if (is.null(tmatrix@row_offset)) {
    uncentered_pca <- compute_pca(tcrossprod(tmatrix@y), k, tmatrix@y)
    return(uncentered_pca)
  }
  raw_residuals_pca(tmatrix@y, k, tmatrix@row_offset, tmatrix@col_offset)
}
