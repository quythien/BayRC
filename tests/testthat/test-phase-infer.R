# tests/testthat/test-phase-infer.R
# Tests for phase_infer() in R/internal.R
# Tests use known-answer inputs so expected outputs are analytically derived.

# Helper: build a G x K phase matrix where all K samples equal a fixed value
make_fixed_phi <- function(G, K, value) {
  matrix(rep(value, G * K), nrow = G, ncol = K)
}

# Helper: run phase_infer with all genes marked Maintained and suppress prints
run_phase_infer <- function(phi1, phi2, shift = 2, bfdr_alpha = 0.25, P = 24) {
  G <- nrow(phi1)
  status <- rep("Maintained", G)
  suppressMessages(
    phase_infer(phi1, phi2, status,
                P = P, shift = shift, bfdr_alpha = bfdr_alpha,
                compute_hdi = FALSE)
  )
}

# ---- Known-answer: zero phase difference ------------------------------------

test_that("phase_infer: all phase diffs = 0 → all genes phase-conserved", {
  G <- 20; K <- 200; P <- 24
  phi <- make_fixed_phi(G, K, 6)   # both conditions peak at 6h
  res <- run_phase_infer(phi, phi, shift = 2, bfdr_alpha = 0.25)

  # prob_shift = mean(|0| >= 2) = 0 for every gene
  expect_true(all(res$prob_shift[1:G] == 0))
  # with prob_shift=0, prob_conserved=1: BFDR for conservation is 0 → all flagged
  expect_equal(sum(res$flag_cons), G)
  expect_equal(sum(res$flag_shift), 0L)
})

# ---- Known-answer: phase difference = P/2 -----------------------------------

test_that("phase_infer: all phase diffs = P/2 → all genes phase-shifted", {
  G <- 20; K <- 200; P <- 24
  phi1 <- make_fixed_phi(G, K, 0)
  phi2 <- make_fixed_phi(G, K, 12)   # P/2 = 12h apart
  # phase_diff = phi1 - phi2 = -12; wrapped to [-12, 12] → -12; |−12| >= 2 always
  res <- run_phase_infer(phi1, phi2, shift = 2, bfdr_alpha = 0.25)

  expect_true(all(res$prob_shift[1:G] == 1))
  expect_equal(sum(res$flag_shift), G)
  expect_equal(sum(res$flag_cons), 0L)
})

# ---- Known-answer: no maintained genes --------------------------------------

test_that("phase_infer: no maintained genes → no flags, no error", {
  G <- 10; K <- 100; P <- 24
  phi1 <- make_fixed_phi(G, K, 6)
  phi2 <- make_fixed_phi(G, K, 8)
  status <- rep("Non-rhythmic", G)
  suppressMessages(
    res <- phase_infer(phi1, phi2, status, P = P, shift = 2, compute_hdi = FALSE)
  )
  expect_equal(sum(res$flag_shift),        0L)
  expect_equal(sum(res$flag_cons),         0L)
  expect_equal(sum(res$flag_undetermined), 0L)
})

# ---- Return structure --------------------------------------------------------

test_that("phase_infer returns required fields", {
  G <- 5; K <- 100; P <- 24
  phi1 <- make_fixed_phi(G, K, 3)
  phi2 <- make_fixed_phi(G, K, 5)
  res  <- run_phase_infer(phi1, phi2)
  expected <- c("flag_shift", "flag_cons", "flag_undetermined",
                "prob_shift", "prob_conserved",
                "peak1", "peak2",
                "BFDR_shift", "BFDR_cons")
  expect_true(all(expected %in% names(res)))
})

# ---- Shift threshold sensitivity --------------------------------------------

test_that("phase_infer: shift=0 → all genes phase-shifted (|diff|>=0 always)", {
  G <- 10; K <- 100; P <- 24
  phi <- make_fixed_phi(G, K, 6)  # diff = 0, but |0| >= 0 is TRUE
  res <- run_phase_infer(phi, phi, shift = 0, bfdr_alpha = 0.25)
  expect_true(all(res$prob_shift[1:G] == 1))
})

test_that("phase_infer: larger shift yields fewer phase-shifted genes", {
  G <- 30; K <- 200; P <- 24
  set.seed(42)
  # Phase diffs uniformly distributed in [-12, 12]
  phi1 <- matrix(runif(G * K, 0, 24), G, K)
  phi2 <- matrix(runif(G * K, 0, 24), G, K)
  res2 <- run_phase_infer(phi1, phi2, shift = 2,  bfdr_alpha = 0.20)
  res6 <- run_phase_infer(phi1, phi2, shift = 6,  bfdr_alpha = 0.20)
  # stricter threshold → fewer (or equal) genes called shifted
  expect_lte(sum(res6$flag_shift), sum(res2$flag_shift))
})

# ---- prob_shift is in [0, 1] ------------------------------------------------

test_that("phase_infer: prob_shift values are in [0, 1]", {
  G <- 15; K <- 150; P <- 24
  set.seed(7)
  phi1 <- matrix(runif(G * K, 0, 24), G, K)
  phi2 <- matrix(runif(G * K, 0, 24), G, K)
  res  <- run_phase_infer(phi1, phi2, shift = 2)
  maintained <- !is.na(res$prob_shift)
  expect_true(all(res$prob_shift[maintained] >= 0))
  expect_true(all(res$prob_shift[maintained] <= 1))
})

# ---- prob_shift + prob_conserved = 1 ----------------------------------------

test_that("phase_infer: prob_shift + prob_conserved = 1 for maintained genes", {
  G <- 10; K <- 100; P <- 24
  set.seed(3)
  phi1 <- matrix(runif(G * K, 0, 24), G, K)
  phi2 <- matrix(runif(G * K, 0, 24), G, K)
  res  <- run_phase_infer(phi1, phi2, shift = 3)
  maintained <- !is.na(res$prob_shift)
  expect_equal(res$prob_shift[maintained] + res$prob_conserved[maintained],
               rep(1, sum(maintained)), tolerance = 1e-10)
})
