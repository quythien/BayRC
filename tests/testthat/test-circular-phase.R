test_that("normalize_angle is periodic, so adding P cannot reorder bounds", {
  # adding P is the identity after normalisation
  expect_equal(BayRC:::normalize_angle(1.05 + 24), 1.05)
  expect_equal(BayRC:::normalize_angle(-10 + 24, a = -12), -10)
})

test_that("circular_width is positive and wraps correctly", {
  expect_equal(circular_width(2, 6), 4)
  # an arc from ZT 20.97 forward to ZT 1.05 crosses the seam
  expect_equal(circular_width(20.97, 1.05), 4.08, tolerance = 1e-8)
  expect_true(circular_width(20.97, 1.05) > 0)
  # plain subtraction is negative across the seam
  expect_lt(1.05 - 20.97, 0)
})

test_that("a phase posterior at ZT 0 gives an interval containing its estimate", {
  set.seed(15213)
  P <- 24
  # samples tight around ZT 0, so roughly half land just below 24 and half
  # just above 0
  samples <- (rnorm(4000, mean = 0, sd = 1)) %% P
  est <- BayRC:::get_t_phi_CI_est(samples, P = P, credMass = 0.95, burn = 1)

  expect_named(est, c("phi.Est", "phi.Lower", "phi.Upper"))
  expect_true(all(is.finite(est)))

  # the arc wraps, so the bounds come back out of order on purpose
  expect_lt(est[["phi.Upper"]], est[["phi.Lower"]])

  # the interval must still contain its own point estimate, and the circular
  # width must be positive and far short of a full cycle
  expect_true(BayRC:::in_circular_interval(est[["phi.Lower"]], est[["phi.Upper"]],
                                   point = est[["phi.Est"]], P = P))
  w <- circular_width(est[["phi.Lower"]], est[["phi.Upper"]], P = P)
  expect_gt(w, 0)
  expect_lt(w, P)

  # a plain range test does not apply across the seam
  expect_false(est[["phi.Est"]] >= est[["phi.Lower"]] &&
               est[["phi.Est"]] <= est[["phi.Upper"]])
})

test_that("a phase posterior away from the seam behaves like an ordinary interval", {
  set.seed(15213)
  P <- 24
  samples <- (rnorm(4000, mean = 12, sd = 1)) %% P
  est <- BayRC:::get_t_phi_CI_est(samples, P = P, credMass = 0.95, burn = 1)

  expect_lt(est[["phi.Lower"]], est[["phi.Upper"]])
  expect_true(BayRC:::in_circular_interval(est[["phi.Lower"]], est[["phi.Upper"]],
                                   point = est[["phi.Est"]], P = P))
  expect_equal(circular_width(est[["phi.Lower"]], est[["phi.Upper"]], P = P),
               est[["phi.Upper"]] - est[["phi.Lower"]])
})
