NROW <- 5
NCOL <- 20
TOL <- 1e-5

compute_pearson_residuals <- function(y) {
    n <- colSums(y)
    rate <- rowSums(y) / sum(n)
    rate_n <- outer(rate, n)
    (y - rate_n) / sqrt(rate_n)
}

compute_deviance_residuals <- function(y) {
    n <- colSums(y)
    rate <- rowSums(y) / sum(n)
    rate_n <- outer(rate, n)
    safe_log <- function(a) {
        ifelse(a == 0, 0, log(a))
    }
    deviance <- 2 * (y * safe_log(y / rate_n) - (y - rate_n))
    sign(y - rate_n) * sqrt(deviance)
}

# Generates a count matrix with Poisson data.
generate_data <- function(nrow = NROW, ncol = NCOL, lambda = 1, seed = 12345) {
    set.seed(seed)
    data <- rpois(n = nrow * ncol, lambda = lambda)
    rm(.Random.seed, envir = globalenv())
    matrix(data, nrow = nrow, ncol = ncol)
}

# Verifies that two PCA results are equivalent, up to some tolerance.
# In particular, this function checks the variance, the rotation matrix, and the
# rotated data corresponding to the principal components.
validate_principal_components <- function(pc_old, pc_new, tol = TOL) {
    # Principal components with low variance have convergence issues.
    significant_pcs_old <- which(pc_old$sdev > tol)
    significant_pcs_new <- which(pc_new$sdev > tol)
    expect_gt(length(significant_pcs_old), 0)
    expect_equal(significant_pcs_new, significant_pcs_old)

    # Principal components are equal up to sign flips, so compare absolute values.
    max_abs_diff <- function(x1, x2) {
        max(abs(abs(x1) - abs(x2)))
    }
    expect_equal(dim(pc_old$rotation), c(NROW, NROW))
    expect_equal(dim(pc_old$rotation), dim(pc_new$rotation))
    expect_lt(max_abs_diff(pc_new$rotation, pc_old$rotation), tol)

    expect_equal(dim(pc_old$x), c(NCOL, NROW))
    expect_equal(dim(pc_old$x), dim(pc_new$x))
    x_new <- pc_new$x[, significant_pcs_new]
    x_old <- pc_old$x[, significant_pcs_old]
    expect_lt(max_abs_diff(x_new, x_old), tol)
}

test_that("Throws error for invalid residual type", {
    y <- generate_data()
    expect_error(
        pc <- poissonPca(y, k = NROW, transform = "foo"),
        "Invalid transform"
    )
})

test_that("Computes PCA on raw residuals of SparseMatrix", {
    y <- generate_data()

    n <- colSums(y)
    rate <- rowSums(y) / sum(n)
    residuals <- y - outer(rate, n)
    pc_old <- prcomp(t(residuals), center = FALSE)

    # Expect warning because all principal components are computed.
    expect_warning(
        pc_new <- poissonPca(y, k = NROW, center = c(TRUE, TRUE)),
        "all eigenvalues"
    )
    validate_principal_components(pc_old, pc_new)
})

test_that("Computes PCA on Pearson residuals of SparseMatrix", {
    y <- generate_data()

    residuals <- compute_pearson_residuals(y)
    pc_old <- prcomp(t(residuals), center = FALSE)

    # Expect warning because all principal components are computed.
    expect_warning(
        pc_new <- poissonPca(y, k = NROW, transform = "pearson"),
        "all eigenvalues"
    )
    validate_principal_components(pc_old, pc_new)
})

test_that("Computes PCA on deviance residuals of SparseMatrix", {
    y <- generate_data()

    residuals <- compute_deviance_residuals(y)
    pc_old <- prcomp(t(residuals), center = FALSE)

    # Expect warning because all principal components are computed.
    expect_warning(
        pc_new <- poissonPca(y, k = NROW, transform = "deviance"),
        "all eigenvalues"
    )
    validate_principal_components(pc_old, pc_new)
})
