setClassUnion("numericOrNull", c("numeric", "NULL"))
setClassUnion("functionOrNull", c("function", "NULL"))

#' Preprocessing Transformations on Sparse Count Matrices
#'
#' Abstraction for the preprocessing transformations that are applied to count
#' data before statistical methods like PCA (e.g., log1p, scaling,
#' row-centering). Note that these transformations may not necessarily preserve
#' the sparsity of the matrix.
#'
#' @slot func Transformation to be applied to sparse count matrix. Should be a
#'   function with a single parameter.
#' @slot center_rows Whether transformed rows should be shifted to have mean
#'   zero.
#' @slot center_cols Whether transformed columns should be shifted to have mean
#'   zero.
#' @slot scale Whether transformed rows should be scaled to have unit variance.
#'
#' @export
setClass(
    "CountTransform",
    slots = c(
        func = "functionOrNull",
        center_rows = "logical",
        center_cols = "logical",
        scale = "logical"
    )
)

# Helper function so that users can specify a single logical value for whether
# the rows of the data should be centered.
.getDuplicatedArgument <- function(arg, name) {
    if (length(arg == 1)) c(arg, FALSE) else arg
}

#' CountTransform Constructor
#'
#' @param func Transformation to be applied to sparse count matrix. Should be a
#'   function with a single parameter.
#' @param center Whether transformed data should be shifted to have mean zero.
#'   Can be specified either as a vector of two logical values, specifying
#'   whether the rows/columns should be centered, respectively, or a single
#'   logical value specifying whether the columns should be centered (for
#'   consistency with [stats::prcomp()]).
#' @param scale Whether transformed rows should be scaled to have unit variance.
#'
#' @return CountTransform object
#' @export
#'
#' @examples
#' triple <- CountTransform(function(x) 3 * x, center = FALSE, scale = FALSE)
CountTransform <- function(func, center = FALSE, scale = FALSE) {
    center <- .getDuplicatedArgument(center, "center")
    new("CountTransform",
        func = func, center_rows = center[1],
        center_cols = center[2], scale = scale
    )
}

#' Transformed Count Matrix
#'
#' Representation of a sparse count matrix after a CountTransform is applied.
#'
#' @slot y SparseMatrix object.
#' @slot row_offset,col_offset Vectors whose product
#'   `outer(row_offset, col_offset)` represents the residual between `y` and a
#'   dense transformation of `y` (e.g., row-centered `y`).
#'
#' @export
setClass(
    "TransformedMatrix",
    slots = c(
        y = "SparseMatrix",
        row_offset = "numericOrNull",
        col_offset = "numericOrNull"
    )
)

#' TransformedMatrix Constructor.
#'
#' @param y SparseMatrix object.
#' @param transform Transformation to apply to `y`.
#'
#' @return TransformedMatrix object.
#' @export
#'
#' @examples
#' mat <- as(matrix(c(1:9), nrow = 3, ncol = 3), "SVT_SparseMatrix")
#' triple <- CountTransform(function(x) 3 * x, center = FALSE, scale = FALSE)
#' tripled_mat <- TransformedMatrix(mat, triple)
TransformedMatrix <- function(y, transform) {
    # Apply the transformation.
    if (!is.null(transform@func)) {
        y[nzwhich(y)] <- transform@func(y[nzwhich(y)])
    }

    # Scale the rows if requested.
    if (transform@scale) {
        # Calculate the standard deviations of the rows.
        sds <- sqrt((rowSums(y^2) - rowSums(y)^2 / ncol(y)) / (ncol(y) - 1))
        nz_ind <- nzwhich(y)
        nz_rows <- .nzrows(y, nz_ind)
        y[nz_ind] <- y[nz_ind] / sds[nz_rows]
    }

    # Store the row/column centers if requested.
    col_offset <- NULL
    row_offset <- NULL
    if (transform@center_rows && transform@center_cols) {
        col_offset <- colSums(y)
        row_offset <- rowSums(y) / sum(col_offset)
    } else if (transform@center_rows) {
        col_offset <- rep(1 / ncol(y), ncol(y))
        row_offset <- rowSums(y)
    } else if (transform@center_cols) {
        col_offset <- colSums(y)
        row_offset <- rep(1 / nrow(y), nrow(y))
    }
    new("TransformedMatrix",
        y = y, row_offset = row_offset,
        col_offset = col_offset
    )
}

setMethod("as.matrix", "TransformedMatrix", function(x, ...) {
    if (is.null(x@row_offset) || is.null(x@row_offset)) {
        as.matrix(x@y)
    } else {
        as.matrix(x@y) - outer(x@row_offset, x@col_offset)
    }
})

setMethod("as.array", "TransformedMatrix", function(x, ...) {
    as.array(as.matrix(x))
})
