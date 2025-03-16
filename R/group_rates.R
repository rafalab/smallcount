#' Row-wise Rates for Groups
#'
#' @param y Sparse matrix (can be a matrix, dgCMatrix, or SparseMatrix).
#' @param g Factor specifying the group for each column.
#'
#' @importFrom SparseArray nzvals nzwhich
#' @export
groupRates <- function(y, g) {
  y <- .convertToSparse(y)

  if(!is.factor(g)){
    warning("Coercing g into a factor")
    g <- as.factor(g)
  }

  # Compute row sums for each group
  nz_ind <- nzwhich(y)
  nz_rows <- .nzrows(y, nz_ind)
  nz_cols <- .nzcols(y, nz_ind)
  group_sums <- tapply(nzvals(y),
                       list(factor(nz_rows, levels=1:nrow(y)),
                            g[nz_cols]),
                       sum, default = 0)

  # Standardize column sums to 1
  rates <- sweep(group_sums, 2, colSums(group_sums), FUN = .safeDivide)
  colnames(rates) <- levels(g)
  rownames(rates) <- rownames(y)
  return(rates)
}
