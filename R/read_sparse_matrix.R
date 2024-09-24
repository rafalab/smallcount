unzip_matrix_file <- function(file, unzip_func) {
  # Get the file extension preceding the .gz or .bz2 extension.
  unzipped_ext <- sub(".*\\.(.*)\\.(gz|bz2)", "\\1", file)
  # Special case for .tgz and .tbz2 files.
  if (grepl(".*\\.(tgz|tbz2)", file)) {
    unzipped_ext <- "tar"
  }
  # Unzip the compressed file.
  output_file <- tempfile(fileext = paste0(".", unzipped_ext))
  unzip_func(file, output_file, remove = FALSE)
  # Decompress tarballs in a temporary directory.
  if (unzipped_ext == "tar") {
    dir_suffix <- paste0(sample(letters, 10), collapse = "")
    temp_dir <- file.path(tempdir(), dir_suffix)
    dir.create(temp_dir, recursive = TRUE)
    utils::untar(output_file, exdir = temp_dir)
    output_file <- list.files(temp_dir, full.names = TRUE)[1]
  }
  output_file
}

#' Unzip (and untar) a .(t)gz or .(t)bz2 sparse matrix file
#' @note This function is a no-op for all other file extensions.
#'
#' @param file character(1) path to a potentially compressed sparse matrix file.
#'
#' @return character(1) path to unzipped file contents.
#' 
#' @import R.utils
#' @import tools
get_decompressed_matrix_file <- function(file) {
  output_file <- file
  compressed_file_ext <- tolower(tools::file_ext(file))
  if (compressed_file_ext == "gz" || compressed_file_ext == "tgz") {
    output_file <- unzip_matrix_file(file, R.utils::gunzip)
  } else if (compressed_file_ext == "bz2" || compressed_file_ext == "tbz2") {
    output_file <- unzip_matrix_file(file, R.utils::bunzip2)
  }
  output_file
}

#' Read a SparseMatrix object from a file
#'
#' @param file character(1) path to a sparse matrix file (.h5, .mtx, and .csv
#'   extensions are supported, compressed or decompressed).
#' @param representation character(1) internal SparseMatrix representation,
#'   either "svt" or "coo". (Default "svt").
#'
#' @import Rcpp
#' @import SparseArray
#'
#' @examples
#' data("tenx_subset")  # Original dataset
#' new <- read_sparse_matrix(system.file(
#'   "extdata/tenx_subset.csv.gz", package = "smallcount"))
#' identical(new, tenx_subset)
#' @export
#' @useDynLib smallcount
read_sparse_matrix <- function(file, representation = "svt") {
  decompressed_file <- get_decompressed_matrix_file(file)
  cppReadSparseMatrix(decompressed_file, representation)
}
