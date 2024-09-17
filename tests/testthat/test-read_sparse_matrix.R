MATRIX_FILENAME <- "small_dense_square"

# Verifies that a SparseMatrix object read from a test file has shape 3x3 and
# contains the values 1-9 row-wise.
validate_test_matrix <- function(sparse_matrix) {
  expected_dim <- c(3, 3)
  expected_matrix <- t(matrix(c(1:9), nrow=3, ncol=3))

  dim <- sparse_matrix@dim
  expect_equal(dim, expected_dim)
  expect_equal(matrix(sparse_matrix, nrow=dim[1], ncol=dim[2]), expected_matrix)
}

test_that("Reads .h5 file", {
  test_file <- test_path("testdata", paste0(MATRIX_FILENAME, ".h5"))

  svt_matrix <- read_sparse_matrix(test_file)
  validate_test_matrix(svt_matrix)

  coo_matrix <- read_sparse_matrix(test_file, representation="coo")
  validate_test_matrix(coo_matrix)
})

test_that("Reads bzipped .mtx file", {
  test_file <- test_path("testdata", paste0(MATRIX_FILENAME, ".tbz2"))

  svt_matrix <- read_sparse_matrix(test_file)
  validate_test_matrix(svt_matrix)

  coo_matrix <- read_sparse_matrix(test_file, representation="coo")
  validate_test_matrix(coo_matrix)
})

test_that("Reads gzipped .csv file", {
  test_file <- test_path("testdata", paste0(MATRIX_FILENAME, ".tar.gz"))
  expected_dimnames <- list(c("r1", "r2", "r3"), c("c1", "c2", "c3"))

  svt_matrix <- read_sparse_matrix(test_file)
  validate_test_matrix(svt_matrix)
  expect_equal(svt_matrix@dimnames, expected_dimnames)

  coo_matrix <- read_sparse_matrix(test_file, representation="coo")
  validate_test_matrix(coo_matrix)
  expect_equal(coo_matrix@dimnames, expected_dimnames)
})
