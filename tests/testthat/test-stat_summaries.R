# Generates a count matrix with Poisson data.
generate_data <- function(nrow, ncol, lambda = 3, seed = 12345) {
  set.seed(seed)
  data <- rpois(n = nrow * ncol, lambda = lambda)
  rm(.Random.seed, envir = globalenv())
  matrix(data, nrow = nrow, ncol = ncol)
}

# Computes the Poisson deviance with a brute-force implementation.
brute_force_deviance <- function(y) {
  n <- colSums(y)
  rate <- rowSums(y) / sum(n)
  mu <- outer(rate, n)
  safe_log <- function(a) {
    ifelse(a == 0, 0, log(a))
  }
  deviance <- 2 * (y * safe_log(y / mu) - (y - mu))
  rowSums(deviance)
}

# Computes the Poisson dispersion with a brute-force implementation.
brute_force_dispersion <- function(y) {
  n <- colSums(y)
  rate <- rowSums(y) / sum(n)
  mu <- outer(rate, n)
  safe_divide <- function(a, b) {
    ifelse(b == 0, 0, a / b)
  }
  dispersion <- safe_divide((y - mu)^2, mu) / (ncol(y) - 1)
  rowSums(dispersion)
}

test_that("Computes group rates", {
  counts <- rbind(
    c(1, 0, 1, 0, 0),
    c(0, 1, 0, 1, 0),
    c(1, 0, 1, 0, 0),
    c(0, 0, 0, 0, 0)
  )
  groups <- as.factor(c("A", "B", "A", "B", "C"))
  rates <- groupRates(counts, groups)

  expected_rates <- rbind(
    c(0.5, 0, 0),
    c(0.0, 1, 0),
    c(0.5, 0, 0),
    c(0.0, 0, 0)
  )
  colnames(expected_rates) <- levels(groups)

  expect_equal(rates, expected_rates)
})

test_that("Computes deviance for 1D Poisson data", {
  ncol <- 100
  counts <- generate_data(nrow = 1, ncol = ncol)
  # Calculate deviance, giving equal weight to each column
  n <- rep(sum(counts) / ncol, ncol)
  deviance <- poissonDeviance(counts, n = n)

  poisson_glm <- glm(t(counts) ~ 1, family = poisson)
  expected_deviance <- poisson_glm$null.deviance
  expect_equal(deviance, expected_deviance)
})

test_that("Computes deviance for Poisson matrix", {
  counts <- generate_data(nrow = 50, ncol = 100)
  deviance <- poissonDeviance(counts)

  expected_deviance <- brute_force_deviance(counts)
  expect_equal(deviance, expected_deviance)
})

test_that("Computes dispersion for 1D Poisson data", {
  ncol <- 100
  counts <- generate_data(nrow = 1, ncol = ncol)
  # Calculate dispersion, giving equal weight to each column
  n <- rep(sum(counts) / ncol, ncol)
  dispersion <- poissonDispersion(counts, n = n)

  quasi_poisson_glm <- glm(t(counts) ~ 1, family = quasipoisson)
  expected_dispersion <- summary(quasi_poisson_glm)$dispersion
  expect_equal(dispersion, expected_dispersion)
})

test_that("Computes dispersion for Poisson matrix", {
  counts <- generate_data(nrow = 50, ncol = 100)
  dispersion <- poissonDispersion(counts)

  expected_dispersion <- brute_force_dispersion(counts)
  expect_equal(dispersion, expected_dispersion)
})
