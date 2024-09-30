test_that("Identity transformation succeeds", {
  data <- rbind(
    c(1, 2, 3),
    c(4, 5, 6),
    c(7, 8, 9)
  )
  sparse_data <- as(data, "SparseMatrix")

  tmatrix <- TransformedMatrix(sparse_data, identity_transform())
  expect_equal(as.matrix(tmatrix), data)
})

test_that("Row centering succeeds", {
  data <- rbind(
    c(1, 2, 3),
    c(2, 4, 6),
    c(3, 6, 9)
  )
  sparse_data <- as(data, "SparseMatrix")

  row_centered_data <- rbind(
    c(-1, 0, 1),
    c(-2, 0, 2),
    c(-3, 0, 3)
  )

  row_center_transform <- identity_transform(center=TRUE)
  tmatrix <- TransformedMatrix(sparse_data, row_center_transform)
  expect_equal(as.matrix(tmatrix), row_centered_data)
})

test_that("Row normalization succeeds", {
  data <- rbind(
    c(1, 0, -1),
    c(3, 1, -1),
    c(6, 1, -4)
  )
  sparse_data <- as(data, "SparseMatrix")

  row_normalized_data <- rbind(
    c(1, 0, -1),
    c(1, 0, -1),
    c(1, 0, -1)
  )

  row_normalize_transform <- identity_transform(center=TRUE, scale=TRUE)
  tmatrix <- TransformedMatrix(sparse_data, row_normalize_transform)
  expect_equal(as.matrix(tmatrix), row_normalized_data)
})

test_that("Column centering succeeds", {
  data <- rbind(
    c(1, 2, 3),
    c(2, 4, 6),
    c(3, 6, 9)
  )
  sparse_data <- as(data, "SparseMatrix")

  col_centered_data <- rbind(
    c(-1, -2, -3),
    c( 0,  0,  0),
    c( 1,  2,  3)
  )

  col_center_transform <- identity_transform(center=c(FALSE, TRUE))
  tmatrix <- TransformedMatrix(sparse_data, col_center_transform)
  expect_equal(as.matrix(tmatrix), col_centered_data)
})

test_that("Log1p transformations succeed", {
  data <- rbind(
    c(1, 2, 3),
    c(4, 5, 6),
    c(7, 8, 9)
  )
  sparse_data <- as(data, "SparseMatrix")

  tmatrix1 <- TransformedMatrix(sparse_data, log1p_transform())
  expect_equal(as.matrix(tmatrix1), log(data + 1))

  tmatrix2 <- TransformedMatrix(sparse_data, cpm_log1p_transform())
  expect_equal(as.matrix(tmatrix2), log(data/1000000 + 1))

  tmatrix3 <- TransformedMatrix(sparse_data, scaled_log1p_transform(3))
  expect_equal(as.matrix(tmatrix3), log(3 * data + 1))
})
