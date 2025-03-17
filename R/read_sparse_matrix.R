# End directory names with a trailing path separator.
standardize_directory_name <- function(dir) {
    if (substr(dir, nchar(dir), nchar(dir)) != .Platform$file.sep) {
        return(paste0(dir, .Platform$file.sep))
    }
    dir
}

# Returns a subdirectory of tempdir()
generate_temp_dir <- function() {
    dir_suffix <- paste0(sample(letters, 20), collapse = "")
    temp_dir <- file.path(tempdir(), dir_suffix)
    dir.create(temp_dir, recursive = TRUE)
    standardize_directory_name(temp_dir)
}

unzip_helper <- function(file, unzip_func, temp_dir) {
    if (!file.exists(file)) {
        stop("Invalid file. File \"", file, "\" does not exist.")
    }
    # Get the file extension preceding the .gz or .bz2 extension.
    unzipped_ext <- sub(".*\\.(.*)\\.(gz|bz2)", "\\1", file)
    # Special case for .tgz and .tbz2 files.
    if (grepl(".*\\.(tgz|tbz2)", file)) {
        unzipped_ext <- "tar"
    }

    # Unzip the compressed file.
    file_name <- tools::file_path_sans_ext(basename(file), compression = TRUE)
    file_ext <- paste0(".", unzipped_ext)
    output_file <- file.path(temp_dir, paste0(file_name, file_ext))
    file.create(output_file)
    unzip_func(file, output_file, overwrite = TRUE, remove = FALSE)

    # Decompress tarballs in a temporary directory.
    if (unzipped_ext == "tar") {
        temp_dir <- generate_temp_dir()
        utils::untar(output_file, exdir = temp_dir)
        output_file <- list.files(temp_dir,
            full.names = TRUE,
            include.dirs = TRUE
        )[1]
    }
    output_file
}

#' Decompress (and untar) .(t)gz or .(t)bz2 sparse matrix data at a file path
#' @note This function is a no-op for all other file extensions.
#'
#' @param filepath character(1) path to potentially compressed matrix data.
#' @param temp_dir character(1) temporary directory for unzipped files.
#'
#' @return character(1) path to unzipped file contents.
#'
#' @import R.utils
#' @import tools
#' @keywords internal
unzip_file <- function(filepath, temp_dir = tempdir()) {
    output_file <- filepath
    compressed_file_ext <- tolower(tools::file_ext(filepath))
    if (compressed_file_ext == "gz" || compressed_file_ext == "tgz") {
        output_file <- unzip_helper(filepath, R.utils::gunzip, temp_dir)
    } else if (compressed_file_ext == "bz2" || compressed_file_ext == "tbz2") {
        output_file <- unzip_helper(filepath, R.utils::bunzip2, temp_dir)
    }
    output_file
}

validate_tenx_directory <- function(directory, prefix) {
    tenx_file_list <- "(matrix.mtx)|(barcodes.tsv)|(features.tsv)|(genes.tsv)"
    tenx_file_pattern <- paste0("^", prefix, tenx_file_list)
    temp_dir <- generate_temp_dir()
    # Potentially unzip the 10x Genomics files in directory.
    for (file in list.files(directory,
        pattern = tenx_file_pattern,
        full.names = TRUE
    )) {
        unzipped_file <- unzip_file(file, temp_dir)
        unzipped_directory <- standardize_directory_name(dirname(unzipped_file))
        # Update directory if files were unzipped into temp_dir.
        if (directory != unzipped_directory) {
            directory <- temp_dir
        }
    }

    prefix_directory <- paste0(directory, prefix)
    if (!file.exists(paste0(prefix_directory, "matrix.mtx"))) {
        stop("Invalid file directory. Matrix could not be found.")
    } else if (!file.exists(paste0(prefix_directory, "barcodes.tsv"))) {
        stop("Invalid file directory. Barcodes could not be found.")
    } else if (!file.exists(paste0(prefix_directory, "features.tsv")) &&
        !file.exists(paste0(prefix_directory, "genes.tsv"))) {
        stop("Invalid file directory. Features/genes could not be found.")
    }
    return(prefix_directory)
}

#' Check that a path is a valid file, file directory, or file prefix
#'
#' @param filepath character(1) path to matrix data.
#'
#' @return character(1) potentially updated file path.
#'
#' @import tools
#' @keywords internal
validate_sample <- function(filepath) {
    file_ext <- tolower(tools::file_ext(filepath))
    if (file_ext == "h5" || file_ext == "csv") {
        # .h5 or .csv files.
        if (!file.exists(filepath)) {
            stop("Invalid file. File \"", filepath, "\" does not exist.")
        }
        return(filepath)
    }

    prefix <- ""
    directory <- standardize_directory_name(filepath)
    # Handle the case where path is a file prefix.
    if (!dir.exists(filepath)) {
        prefix <- basename(filepath)
        directory <- standardize_directory_name(dirname(filepath))
    }
    directory <- validate_tenx_directory(directory, prefix)
    return(directory)
}

#' Load data from a 10X Genomics experiment
#'
#' Creates a \code{\link[SparseArray]{SparseMatrix}} from the CellRanger output
#' directories for 10X Genomics data.
#'
#' @param sample character(1) directory name corresponding to a 10X sample. The
#'   directory should contain a matrix file, a gene/feature annotation file, and
#'   a barcode annotation file.
#'
#'   Alternatively, the string may contain a path to a HDF5 file in the sparse
#'   matrix format generated by 10X.
#'
#'   Alternatively, the string may contain a prefix of names for the three-file
#'   system described above, where the rest of the name of each file follows the
#'   standard 10X output.
#' @param col.names logical(1) indicating whether the columns of the matrix
#'   should be named with the cell barcodes.
#' @param row.names character(1) specifying whether to use Ensembl IDs ("id") or
#'   gene symbols ("symbol") as row names. If using symbols, the Ensembl ID will
#'   be appended to disambiguate in case the same symbol corresponds to multiple
#'   Ensembl IDs.
#' @param genome character(1) specifying the genome for HDF5 files output by
#'   CellRanger v2.
#'
#' @return A \code{\link[SparseArray]{SparseMatrix}} object containing count
#'   data for each gene (row) and cell (column) in \code{sample}.
#'
#' @details The signature of this function and its corresponding documentation
#' has largely been adapted from the \code{Read10xCounts} function in the
#' \pkg{DropletUtils} package.
#'
#' Note that user-level manipulation of sparse matrices requires loading of the
#' \pkg{SparseArray} package. Otherwise, calculation of \code{rowSums},
#' \code{colSums}, etc. will result in errors.
#'
#' @import Rcpp
#' @import Rhdf5lib
#' @import SparseArray
#'
#' @examples
#' data("tenx_subset") # Original dataset
#' new <- readSparseMatrix(system.file(
#'     "extdata/tenx_subset.csv.gz",
#'     package = "smallcount"
#' ))
#' identical(new, tenx_subset)
#'
#' @references Zheng GX, Terry JM, Belgrader P, and others (2017). Massively
#' parallel digital transcriptional profiling of single cells. \emph{Nat Commun}
#' 8:14049.
#'
#' 10X Genomics (2017). Gene-Barcode Matrices.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/2.2/output/matrices}
#'
#' 10X Genomics (2018). Feature-Barcode Matrices.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices}
#'
#' 10X Genomics (2018). HDF5 Gene-Barcode Matrix Format.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/2.2/advanced/h5_matrices}
#'
#' 10X Genomics (2018). HDF5 Feature Barcode Matrix Format.
#' \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices}
#'
#' @export
#' @useDynLib smallcount
readSparseMatrix <- function(sample,
                             col.names = FALSE,
                             row.names = c("id", "symbol"),
                             genome = NULL) {
    sample <- unzip_file(sample)
    sample <- validate_sample(sample)
    id_row_names <- match.arg(row.names) == "id"
    genome <- ifelse(is.null(genome), "", genome)
    features_tsv <- file.exists(paste0(sample, "features.tsv"))
    cppReadSparseMatrix(sample, col.names, id_row_names, genome, features_tsv)
}
