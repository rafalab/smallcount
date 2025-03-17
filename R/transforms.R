#' Identity Transformation
#'
#' CountTransform mapping x to itself.
#'
#' @note Internally, the identity function is represented by a NULL value.
#'
#' @param center Whether transformed data should be shifted to have mean zero.
#'   Can be specified either as a vector of two logical values, specifying
#'   whether the rows/columns should be centered, respectively, or a single
#'   logical value specifying whether the columns should be centered (for
#'   consistency with [stats::prcomp()]).
#' @param scale Whether transformed rows should be scaled to have unit variance.
#'
#' @return CountTransform object.
#' @export
identity_transform <- function(center = FALSE, scale = FALSE) {
    CountTransform(NULL, center, scale)
}

#' Log1p Transformation
#'
#' CountTransform mapping x to log(x + 1).
#'
#' @inheritParams identity_transform
#'
#' @return CountTransform object.
#' @export
log1p_transform <- function(center = FALSE, scale = FALSE) {
    CountTransform(log1p, center, scale)
}

#' Scaled Log1p Transformation
#'
#' CountTransform mapping x to log(`coef` * x + 1).
#'
#' @param coef Scaling coefficient.
#' @inheritParams identity_transform
#'
#' @return CountTransform object.
#' @export
scaled_log1p_transform <- function(coef, center = FALSE, scale = FALSE) {
    CountTransform(function(y) log1p(coef * y), center, scale)
}

#' CPM Log1p Transformation
#'
#' CountTransform mapping x to log(x/1000000 + 1).
#'
#' @inheritParams identity_transform
#'
#' @return CountTransform object.
#' @export
cpm_log1p_transform <- function(center = FALSE, scale = FALSE) {
    CountTransform(function(y) log1p(y / 1e6), center, scale)
}
