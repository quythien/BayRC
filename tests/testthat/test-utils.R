# tests/testthat/test-utils.R
# Tests for R/utils.R — adjust.to.2pi()

test_that("adjust.to.2pi keeps values already in [0, 2pi)", {
  expect_equal(adjust.to.2pi(0),   0)
  expect_equal(adjust.to.2pi(pi),  pi)
  expect_equal(adjust.to.2pi(2*pi - 0.001), 2*pi - 0.001)
})

test_that("adjust.to.2pi wraps negative angles into [0, 2pi)", {
  result <- adjust.to.2pi(-pi)
  expect_gte(result, 0)
  expect_lt(result, 2*pi)
  expect_equal(result, pi, tolerance = 1e-10)
})

test_that("adjust.to.2pi wraps angles > 2pi back into [0, 2pi)", {
  result <- adjust.to.2pi(3 * pi)
  expect_gte(result, 0)
  expect_lt(result, 2*pi)
  expect_equal(result, pi, tolerance = 1e-10)
})

test_that("adjust.to.2pi wraps large positive multiples correctly", {
  result <- adjust.to.2pi(5 * pi)   # 5pi = 2*2pi + pi
  expect_gte(result, 0)
  expect_lt(result, 2*pi)
  expect_equal(result, pi, tolerance = 1e-10)
})

test_that("adjust.to.2pi is identity on exact 0", {
  expect_equal(adjust.to.2pi(0), 0)
})
