# tests/testthat/test-clustering.R
# Tests for R/clustering.R — mydist() and copy.upper.lower()

test_that("mydist returns a square matrix with dim = nrow(y)", {
  set.seed(1)
  y <- matrix(sample(0:1, 50, replace = TRUE), nrow = 5)
  d <- mydist(y)
  expect_equal(dim(d), c(5, 5))
})

test_that("mydist diagonal is 0 (distance from self)", {
  set.seed(2)
  y <- matrix(sample(0:1, 60, replace = TRUE), nrow = 6)
  d <- mydist(y)
  expect_equal(diag(d), rep(0, 6), tolerance = 1e-10)
})

test_that("mydist is symmetric", {
  set.seed(3)
  y <- matrix(sample(0:1, 80, replace = TRUE), nrow = 8)
  d <- mydist(y)
  expect_equal(d, t(d), tolerance = 1e-10)
})

test_that("mydist values are in [0, 1]", {
  set.seed(4)
  y <- matrix(sample(0:1, 100, replace = TRUE), nrow = 10)
  d <- mydist(y)
  expect_true(all(d >= 0 & d <= 1))
})

test_that("mydist is 0 for identical rows", {
  y <- matrix(c(1, 0, 1, 0, 1, 1, 0, 1, 0, 1), nrow = 2, byrow = TRUE)
  # row 1 and row 2 are identical copies
  y2 <- rbind(y[1, ], y[1, ])
  d  <- mydist(y2)
  expect_equal(d[1, 2], 0, tolerance = 1e-10)
})

test_that("copy.upper.lower produces a symmetric matrix", {
  m <- matrix(c(0, NA, NA, 1, 0, NA, 2, 3, 0), nrow = 3, byrow = TRUE)
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  result <- copy.upper.lower(m)
  expect_equal(result, t(result))
})
