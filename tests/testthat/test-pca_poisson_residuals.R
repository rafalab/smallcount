NROW <- 5
NCOL <- 20
TOL <- 1e-5

compute_pearson_residuals <- function(y) {
  n <- colSums(y)
  rate <- rowSums(y) / sum(n)
  rate_n <- outer(rate, n)
  (y - rate_n) / sqrt(rate_n)
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
  expect_equal(length(significant_pcs_old), NROW - 1)
  expect_equal(significant_pcs_new, significant_pcs_old)

  # Principal components are equal up to sign flips, so compare absolute values.
  max_abs_diff <- function(x1, x2) {
    max(abs(x1) - abs(x2))
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

test_that("Computes PCA on Pearson residuals of SparseMatrix", {
  y <- generate_data()

  residuals <- compute_pearson_residuals(y)
  pc_old <- prcomp(t(residuals), center = FALSE)

  sparse_y <- as(y, "SparseMatrix")
  # Expect warning because all principal components are computed.
  expect_warning(pc_new <- pca_poisson_residuals(sparse_y, k = NROW),
                 "all eigenvalues")

  validate_principal_components(pc_old, pc_new)
})
