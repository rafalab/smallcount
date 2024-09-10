test_that("Reads .h5 file", {
  test_file <- test_path("testdata", "small_dense_square.h5")
  expected_dim <- c(3, 3)
  expected_matrix <- t(matrix(c(1:9), nrow=3, ncol=3))

  svt_matrix <- read_sparse_matrix(test_file)
  dim <- svt_matrix@dim
  expect_equal(dim, expected_dim)
  expect_equal(matrix(svt_matrix, nrow=dim[1], ncol=dim[2]), expected_matrix)

  coo_matrix <- read_sparse_matrix(test_file, representation="coo")
  dim <- coo_matrix@dim
  expect_equal(dim, expected_dim)
  expect_equal(matrix(coo_matrix, nrow=dim[1], ncol=dim[2]), expected_matrix)
})

test_that("Reads .mtx file", {
  test_file <- test_path("testdata", "small_dense_square.tbz2")
  expected_dim <- c(3, 3)
  expected_matrix <- t(matrix(c(1:9), nrow=3, ncol=3))

  svt_matrix <- read_sparse_matrix(test_file)
  dim <- svt_matrix@dim
  expect_equal(dim, expected_dim)
  expect_equal(matrix(svt_matrix, nrow=dim[1], ncol=dim[2]), expected_matrix)

  coo_matrix <- read_sparse_matrix(test_file, representation="coo")
  dim <- coo_matrix@dim
  expect_equal(dim, expected_dim)
  expect_equal(matrix(coo_matrix, nrow=dim[1], ncol=dim[2]), expected_matrix)
})

test_that("Reads .csv file", {
  test_file <- test_path("testdata", "small_dense_square.tar.gz")
  expected_dim <- c(3, 3)
  expected_matrix <- t(matrix(c(1:9), nrow=3, ncol=3))

  svt_matrix <- read_sparse_matrix(test_file)
  dim <- svt_matrix@dim
  expect_equal(dim, expected_dim)
  expect_equal(matrix(svt_matrix, nrow=dim[1], ncol=dim[2]), expected_matrix)

  coo_matrix <- read_sparse_matrix(test_file, representation="coo")
  dim <- coo_matrix@dim
  expect_equal(dim, expected_dim)
  expect_equal(matrix(coo_matrix, nrow=dim[1], ncol=dim[2]), expected_matrix)
})
