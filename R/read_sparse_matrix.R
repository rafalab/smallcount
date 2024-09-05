#' Read a SparseMatrix object from a file
#'
#' @import Rcpp
#' @import SparseArray
#'
#' @param file character(1) path to a sparse matrix file (.h5, .mtx, and .csv
#'   extensions are supported).
#' @param representation character(1) internal SparseMatrix representation,
#'   either "svt" or "coo". (Default "svt")
#'
#' @return SparseMatrix object
#' @export
#' @useDynLib smallcount
read_sparse_matrix <- function(file, representation="svt") {
  cppReadSparseMatrix(file, representation)
}
