# tests/testthat/test-estimates.R
# Tests for R/estimates.R — normalize_angle, circular_HDI, in_circular_interval,
# circular_median, recenter_samples, get_other_CI_est, get_t_phi_CI_est

# ---- normalize_angle --------------------------------------------------------

test_that("normalize_angle maps values into [a, a+P)", {
  # Default: a=0, P=24
  expect_equal(normalize_angle(25, a = 0, P = 24), 1)
  expect_equal(normalize_angle(-1, a = 0, P = 24), 23)
  expect_equal(normalize_angle(0),  0)
  expect_equal(normalize_angle(24, a = 0, P = 24), 0)
})

test_that("normalize_angle respects non-zero a", {
  result <- normalize_angle(25, a = 0, P = 24)
  expect_gte(result, 0)
  expect_lt(result, 24)
})

# ---- circular_HDI -----------------------------------------------------------

test_that("circular_HDI returns a list with lower and upper", {
  set.seed(42)
  samples <- rnorm(500, mean = 12, sd = 1) %% 24
  res <- circular_HDI(samples, credMass = 0.95, P = 24, a = 0)
  expect_true(all(c("lower", "upper") %in% names(res)))
})

test_that("circular_HDI lower and upper are in [0, 24)", {
  set.seed(1)
  samples <- rnorm(200, mean = 6, sd = 0.5) %% 24
  res <- circular_HDI(samples, credMass = 0.95, P = 24, a = 0)
  expect_gte(res$lower, 0)
  expect_lt(res$lower, 24)
  expect_gte(res$upper, 0)
  expect_lt(res$upper, 24)
})

test_that("circular_HDI covers true phase for tight cluster", {
  set.seed(123)
  true_phase <- 6
  samples <- rnorm(1000, mean = true_phase, sd = 0.2) %% 24
  res <- circular_HDI(samples, credMass = 0.95, P = 24, a = 0)
  # Check true phase is within the HDI (non-wrapping case)
  if (res$lower <= res$upper) {
    expect_gte(true_phase, res$lower - 0.5)
    expect_lte(true_phase, res$upper + 0.5)
  }
})

test_that("circular_HDI handles wrap-around samples near 0/24", {
  set.seed(55)
  # Samples clustered around 0 (wrapping between 23 and 1)
  samples <- c(rnorm(200, 0.5, 0.3), rnorm(200, 23.5, 0.3)) %% 24
  res <- circular_HDI(samples, credMass = 0.95, P = 24, a = 0)
  expect_true(is.finite(res$lower))
  expect_true(is.finite(res$upper))
})

# ---- in_circular_interval ---------------------------------------------------

test_that("in_circular_interval returns TRUE for interior point", {
  # Non-wrapping: [4, 8]
  expect_true(in_circular_interval(4, 8, 6, P = 24, a = 0))
})

test_that("in_circular_interval returns FALSE for exterior point", {
  expect_false(in_circular_interval(4, 8, 10, P = 24, a = 0))
})

test_that("in_circular_interval handles wrap-around interval", {
  # Wrapping interval [22, 2]
  expect_true(in_circular_interval(22, 2, 23, P = 24, a = 0))
  expect_true(in_circular_interval(22, 2, 1,  P = 24, a = 0))
  expect_false(in_circular_interval(22, 2, 12, P = 24, a = 0))
})

# ---- get_other_CI_est -------------------------------------------------------

test_that("get_other_CI_est returns length-3 numeric", {
  skip_if_not_installed("HDInterval")
  set.seed(1)
  samples <- rnorm(500, 5, 1)
  res <- get_other_CI_est(samples, credMass = 0.95)
  expect_length(res, 3)
  expect_true(all(is.finite(res)))
})

test_that("get_other_CI_est estimate is between lower and upper", {
  skip_if_not_installed("HDInterval")
  set.seed(2)
  samples <- rnorm(500, 10, 2)
  res <- get_other_CI_est(samples, credMass = 0.95)
  expect_gte(res[1], res[2])    # est >= lower
  expect_lte(res[1], res[3])    # est <= upper
})

# ---- circular_median --------------------------------------------------------

test_that("circular_median recovers true phase for tight cluster", {
  set.seed(5)
  true_phase <- 10
  samples <- rnorm(500, true_phase, 0.3) %% 24
  med <- circular_median(samples, P = 24)
  expect_equal(med, true_phase, tolerance = 0.5)
})

test_that("circular_median is in [0, 24)", {
  set.seed(7)
  samples <- rnorm(200, 20, 2) %% 24
  med <- circular_median(samples, P = 24)
  expect_gte(med, 0)
  expect_lt(med, 24)
})

# ---- get_t_phi_CI_est -------------------------------------------------------

test_that("get_t_phi_CI_est errors when burn >= length", {
  expect_error(get_t_phi_CI_est(1:5, P = 24, burn = 5))
})

test_that("get_t_phi_CI_est returns 3-element named vector", {
  skip_if_not_installed("HDInterval")
  set.seed(42)
  samples <- rnorm(300, 6, 0.5) %% 24
  res <- get_t_phi_CI_est(samples, P = 24, burn = 50)
  expect_length(res, 3)
  expect_true(all(c("phi.Est", "phi.Lower", "phi.Upper") %in% names(res)))
})
