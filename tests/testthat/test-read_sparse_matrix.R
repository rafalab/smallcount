MATRIX_FILENAME <- "small_dense_square"

# Verifies that a SparseMatrix object read from a test file has shape 3x3 and
# contains the values 1-9 row-wise.
validate_test_matrix <- function(sparse_matrix) {
  expected_dim <- c(3, 3)
  expected_matrix <- t(matrix(c(1:9), nrow=3, ncol=3))
  expected_dimnames <- list(c("r1", "r2", "r3"), c("c1", "c2", "c3"))

  dim <- sparse_matrix@dim
  expect_equal(dim, expected_dim)
  expect_equal(matrix(sparse_matrix, nrow=dim[1], ncol=dim[2]), expected_matrix)
  expect_equal(sparse_matrix@dimnames, expected_dimnames)
}

test_that("Reads .h5 file (Cell Ranger v3)", {
  matrix_file <- test_path("testdata", paste0(MATRIX_FILENAME, "_v3.h5"))

  svt_matrix <- read_sparse_matrix(matrix_file, col.names=TRUE)
  validate_test_matrix(svt_matrix)
})

test_that("Reads .h5 file (Cell Ranger v2)", {
  matrix_file <- test_path("testdata", paste0(MATRIX_FILENAME, "_v2.h5"))

  svt_matrix <- read_sparse_matrix(matrix_file, col.names=TRUE,
                                   row.names="symbol", genome="genome")
  validate_test_matrix(svt_matrix)
})

test_that("Reads bzipped .mtx directory (Cell Ranger v3)", {
  matrix_file <- test_path("testdata", paste0(MATRIX_FILENAME, "_v3.tbz2"))

  svt_matrix <- read_sparse_matrix(matrix_file, col.names=TRUE, row.names="id")
  validate_test_matrix(svt_matrix)
})

test_that("Reads .mtx directory with prefix (Cell Ranger v2)", {
  matrix_file <- test_path("testdata", paste0(MATRIX_FILENAME, "_v2/prefix_"))

  svt_matrix <- read_sparse_matrix(matrix_file, col.names=TRUE,
                                   row.names="symbol")
  validate_test_matrix(svt_matrix)
})

test_that("Reads gzipped .csv file", {
  matrix_file <- test_path("testdata", paste0(MATRIX_FILENAME, ".tar.gz"))

  svt_matrix <- read_sparse_matrix(matrix_file)
  validate_test_matrix(svt_matrix)
})
