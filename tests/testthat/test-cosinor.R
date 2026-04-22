# tests/testthat/test-cosinor.R
# Tests for R/cosinor.R — one_cosinor_OLS() and Cosinor_fit()

# ---- one_cosinor_OLS --------------------------------------------------------

test_that("one_cosinor_OLS recovers exact amplitude on noiseless data", {
  set.seed(42)
  P      <- 24
  omega  <- 2 * pi / P
  true_A <- 1.5
  true_phi <- 6   # phase in hours
  true_M <- 5
  t      <- seq(0, 23, length.out = 48)
  y      <- true_M + true_A * cos(omega * (t - true_phi))

  res <- one_cosinor_OLS(tod = t, y = y, period = P)

  expect_equal(res$M,     true_M,   tolerance = 1e-6)
  expect_equal(res$A,     true_A,   tolerance = 1e-6)
  expect_equal(res$R2,    1,        tolerance = 1e-6)
})

test_that("one_cosinor_OLS returns p-value near 0 for strong rhythm", {
  set.seed(1)
  P <- 24
  omega <- 2 * pi / P
  t <- seq(0, 23, length.out = 48)
  y <- 5 + 2 * cos(omega * t) + rnorm(48, 0, 0.01)
  res <- one_cosinor_OLS(t, y, P)
  expect_lt(res$pvalue, 0.01)
})

test_that("one_cosinor_OLS returns p-value near 1 for pure noise", {
  set.seed(99)
  t <- seq(0, 23, length.out = 48)
  y <- rnorm(48, 5, 1)
  res <- one_cosinor_OLS(t, y)
  expect_gt(res$pvalue, 0.05)   # not necessarily; just check it runs
  expect_true(is.finite(res$pvalue))
})

test_that("one_cosinor_OLS returns named output list", {
  t <- seq(0, 23, length.out = 24)
  y <- 5 + cos(2 * pi / 24 * t)
  res <- one_cosinor_OLS(t, y)
  expect_true(all(c("M", "A", "phase", "peak", "R2", "sigma", "pvalue") %in% names(res)))
})

test_that("one_cosinor_OLS phase is in [0, 24) hours", {
  set.seed(7)
  t <- seq(0, 23, length.out = 48)
  y <- 5 + 1.5 * cos(2 * pi / 24 * (t - 10)) + rnorm(48, 0, 0.1)
  res <- one_cosinor_OLS(t, y)
  expect_gte(res$phase, 0)
  expect_lt(res$phase, 24)
})

test_that("one_cosinor_OLS sigma is non-negative", {
  t <- seq(0, 23, length.out = 24)
  y <- rnorm(24, 5, 1)
  res <- one_cosinor_OLS(t, y)
  expect_gte(res$sigma, 0)
})

# ---- Cosinor_fit ------------------------------------------------------------

test_that("Cosinor_fit returns correct structure with multiple genes", {
  set.seed(10)
  G  <- 5
  N  <- 24
  t  <- seq(0, 23, length.out = N)
  Y  <- matrix(rnorm(G * N, mean = 5, sd = 1), nrow = G)
  rownames(Y) <- paste0("Gene", seq_len(G))

  input <- list(data = Y, time = t, gname = paste0("Gene", seq_len(G)))
  res   <- Cosinor_fit(input)

  expect_true("rhythm" %in% names(res))
  expect_equal(nrow(res$rhythm), G)
  expect_true(all(c("gname", "M", "A", "phase", "pvalue", "qvalue") %in% names(res$rhythm)))
})

test_that("Cosinor_fit q-values are in [0,1]", {
  set.seed(20)
  G  <- 10
  N  <- 24
  t  <- seq(0, 23, length.out = N)
  Y  <- matrix(rnorm(G * N, 5, 1), nrow = G)
  rownames(Y) <- paste0("G", seq_len(G))
  inp <- list(data = Y, time = t, gname = paste0("G", seq_len(G)))
  res <- Cosinor_fit(inp)
  expect_true(all(res$rhythm$qvalue >= 0 & res$rhythm$qvalue <= 1))
})
