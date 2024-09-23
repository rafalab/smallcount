test_that("Computes group rates", {
  counts <- rbind(
    c(1, 0, 1, 0, 0),
    c(0, 1, 0, 1, 0),
    c(1, 0, 1, 0, 0),
    c(0, 0, 0, 0, 0)
  )
  groups <- as.factor(c("A", "B", "A", "B", "C"))
  rates <- group_rates(counts, groups)

  expected_rates <- rbind(
    c(0.5, 0, 0),
    c(0.0, 1, 0),
    c(0.5, 0, 0),
    c(0.0, 0, 0)
  )
  colnames(expected_rates) <- levels(groups)

  expect_equal(rates, expected_rates)
})
