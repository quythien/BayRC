# tests/testthat/test-bfdr.R
# Tests for bfdr_from_posterior() and transition_classify() in R/internal.R

# ---- bfdr_from_posterior ----------------------------------------------------

test_that("bfdr_from_posterior returns required list elements", {
  probs <- c(0.95, 0.9, 0.8, 0.7, 0.2, 0.1)
  res   <- bfdr_from_posterior(probs, alpha = 0.05)
  expect_true(all(c("threshold", "rhythmic_genes", "n_rhythmic",
                    "bfdr_values", "sorted_probs") %in% names(res)))
})

test_that("bfdr_from_posterior threshold is in (0, 1]", {
  probs <- runif(100, 0, 1)
  res   <- bfdr_from_posterior(probs, alpha = 0.05)
  expect_gte(res$threshold, 0)
  expect_lte(res$threshold, 1)
})

test_that("bfdr_from_posterior n_rhythmic matches logical vector sum", {
  probs <- c(0.97, 0.96, 0.95, 0.3, 0.1, 0.05)
  res   <- bfdr_from_posterior(probs, alpha = 0.05)
  expect_equal(res$n_rhythmic, sum(res$rhythmic_genes))
})

test_that("bfdr_from_posterior returns threshold=1 when no gene passes", {
  probs <- rep(0.01, 10)
  res   <- bfdr_from_posterior(probs, alpha = 0.05)
  expect_equal(res$threshold, 1.0)
  expect_equal(res$n_rhythmic, 0)
})

test_that("bfdr_from_posterior includes all genes when all have prob=1", {
  probs <- rep(1.0, 20)
  res   <- bfdr_from_posterior(probs, alpha = 0.05)
  expect_equal(res$n_rhythmic, 20)
})

test_that("bfdr_from_posterior is monotone non-decreasing in alpha", {
  probs <- c(0.9, 0.85, 0.7, 0.6, 0.3, 0.1)
  n1    <- bfdr_from_posterior(probs, alpha = 0.05)$n_rhythmic
  n2    <- bfdr_from_posterior(probs, alpha = 0.10)$n_rhythmic
  n3    <- bfdr_from_posterior(probs, alpha = 0.20)$n_rhythmic
  expect_lte(n1, n2)
  expect_lte(n2, n3)
})

test_that("bfdr_from_posterior BFDR curve has length equal to input", {
  probs <- runif(50)
  res   <- bfdr_from_posterior(probs, alpha = 0.10)
  expect_length(res$bfdr_values, 50)
})

# ---- transition_classify ----------------------------------------------------

test_that("transition_classify returns required list elements", {
  set.seed(1)
  pA <- runif(50, 0.3, 1)
  pB <- runif(50, 0.3, 1)
  suppressMessages(
    res <- transition_classify(pA, pB, bfdr_alpha = 0.20)
  )
  expect_true(all(c("gain_genes", "loss_genes", "cons_genes",
                    "gain_loss_status", "p_gain", "p_loss") %in% names(res)))
})

test_that("transition_classify p_gain + p_loss + (pA*pB) uses independence form", {
  pA <- c(0.9, 0.1)
  pB <- c(0.1, 0.9)
  # p_gain[1] = (1-0.9)*0.1 = 0.01; p_loss[1] = 0.9*(1-0.1)=0.81
  suppressMessages(
    res <- transition_classify(pA, pB, bfdr_alpha = 0.50)
  )
  expect_equal(res$p_gain[1], (1 - pA[1]) * pB[1], tolerance = 1e-10)
  expect_equal(res$p_loss[1], pA[1] * (1 - pB[1]), tolerance = 1e-10)
})

test_that("transition_classify gains dominant when pB >> pA", {
  set.seed(3)
  n  <- 50
  pA <- rep(0.05, n)  # condition A is mostly arrhythmic
  pB <- rep(0.95, n)  # condition B is mostly rhythmic
  suppressMessages(
    res <- transition_classify(pA, pB, bfdr_alpha = 0.20)
  )
  expect_gt(length(res$gain_genes), length(res$loss_genes))
})

test_that("transition_classify loss dominant when pA >> pB", {
  set.seed(4)
  n  <- 50
  pA <- rep(0.95, n)
  pB <- rep(0.05, n)
  suppressMessages(
    res <- transition_classify(pA, pB, bfdr_alpha = 0.20)
  )
  expect_gt(length(res$loss_genes), length(res$gain_genes))
})
