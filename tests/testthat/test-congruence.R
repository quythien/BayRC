# tests/testthat/test-congruence.R
# Tests for congruence() (R/internal.R) and compute_iteration_jaccard() (R/permutation_sim.R)

# Helper: build a minimal rho list from a matrix
make_rho_list <- function(mat) {
  list(rho = mat)
}

# Helper: build a binary rho matrix (genes x iterations) from per-gene prob
make_rho_matrix <- function(probs, n_iter = 100, seed = 1) {
  set.seed(seed)
  G <- length(probs)
  mat <- matrix(rbinom(G * n_iter, 1, rep(probs, n_iter)), nrow = G)
  mat
}

# ---- congruence() -----------------------------------------------------------

test_that("congruence returns required list elements", {
  mat1 <- make_rho_list(make_rho_matrix(rep(0.8, 10)))
  mat2 <- make_rho_list(make_rho_matrix(rep(0.8, 10), seed = 2))
  res  <- congruence(mat1, mat2)
  expect_true(all(c("congruence_index", "gain_index", "loss_index",
                    "gain_loss_ratio", "n_genes") %in% names(res)))
})

test_that("congruence is 1 when both conditions are identical (p=1)", {
  probs <- rep(1.0, 20)
  mat   <- make_rho_list(matrix(1L, nrow = 20, ncol = 50))
  res   <- congruence(mat, mat)
  expect_equal(res$congruence_index, 1.0, tolerance = 1e-10)
  expect_equal(res$gain_index,  0.0, tolerance = 1e-10)
  expect_equal(res$loss_index,  0.0, tolerance = 1e-10)
})

test_that("congruence is 0 when conditions are perfectly complementary", {
  # A all rhythmic, B all arrhythmic
  rho_A <- make_rho_list(matrix(1L, nrow = 10, ncol = 50))
  rho_B <- make_rho_list(matrix(0L, nrow = 10, ncol = 50))
  res   <- congruence(rho_A, rho_B)
  expect_equal(res$congruence_index, 0.0, tolerance = 1e-10)
  expect_equal(res$loss_index, 1.0, tolerance = 1e-10)
  expect_equal(res$gain_index, 0.0, tolerance = 1e-10)
})

test_that("congruence handles all-zero rho (no rhythmic genes)", {
  rho_Z <- make_rho_list(matrix(0L, nrow = 10, ncol = 50))
  res   <- congruence(rho_Z, rho_Z)
  expect_equal(res$congruence_index, 0.0)
  expect_true(is.na(res$gain_loss_ratio))
})

test_that("congruence indices sum to 1", {
  set.seed(42)
  probs1 <- runif(30, 0.2, 0.8)
  probs2 <- runif(30, 0.2, 0.8)
  mat1   <- make_rho_list(make_rho_matrix(probs1, n_iter = 200))
  mat2   <- make_rho_list(make_rho_matrix(probs2, n_iter = 200, seed = 99))
  res    <- congruence(mat1, mat2)
  # congruence + gain + loss should sum to 1 (they partition the expected union)
  total  <- res$congruence_index + res$gain_index + res$loss_index
  expect_equal(total, 1.0, tolerance = 1e-10)
})

test_that("congruence gain_index > 0 when B is more rhythmic than A", {
  rho_A <- make_rho_list(matrix(0L, nrow = 20, ncol = 100))
  rho_A$rho[1:5, ] <- 1L   # only 5 genes rhythmic in A
  rho_B <- make_rho_list(matrix(0L, nrow = 20, ncol = 100))
  rho_B$rho[1:15, ] <- 1L  # 15 genes rhythmic in B
  res <- congruence(rho_A, rho_B)
  expect_gt(res$gain_index, res$loss_index)
})

# ---- compute_iteration_jaccard() --------------------------------------------

test_that("compute_iteration_jaccard returns correct structure", {
  rho_A <- make_rho_matrix(rep(0.7, 20), n_iter = 50)
  rho_B <- make_rho_matrix(rep(0.7, 20), n_iter = 50, seed = 99)
  res   <- compute_iteration_jaccard(rho_A, rho_B)
  expect_true(all(c("jaccard_obs", "jaccard_exp", "jaccard_adj",
                    "mean_jaccard_obs", "mean_jaccard_adj") %in% names(res)))
})

test_that("compute_iteration_jaccard obs Jaccard in [0,1]", {
  rho_A <- make_rho_matrix(rep(0.5, 15), n_iter = 80)
  rho_B <- make_rho_matrix(rep(0.5, 15), n_iter = 80, seed = 7)
  res   <- compute_iteration_jaccard(rho_A, rho_B)
  valid <- res$jaccard_obs[!is.na(res$jaccard_obs)]
  expect_true(all(valid >= 0 & valid <= 1))
})

test_that("compute_iteration_jaccard mean_jaccard_obs == 1 for identical all-1 matrices", {
  mat <- matrix(1L, nrow = 10, ncol = 50)
  res <- compute_iteration_jaccard(mat, mat)
  expect_equal(res$mean_jaccard_obs, 1.0, tolerance = 1e-10)
})

test_that("compute_iteration_jaccard mean_jaccard_obs == 0 for complementary matrices", {
  mat_A <- matrix(1L, nrow = 10, ncol = 50)
  mat_B <- matrix(0L, nrow = 10, ncol = 50)
  res   <- compute_iteration_jaccard(mat_A, mat_B)
  expect_equal(res$mean_jaccard_obs, 0.0, tolerance = 1e-10)
})

# ---- compute_gain_loss_iterations() -----------------------------------------

test_that("compute_gain_loss_iterations returns correct structure", {
  rho_A <- make_rho_matrix(rep(0.6, 15), n_iter = 50)
  rho_B <- make_rho_matrix(rep(0.6, 15), n_iter = 50, seed = 5)
  res   <- compute_gain_loss_iterations(rho_A, rho_B)
  expect_true(all(c("gain_index", "loss_index", "gain_loss_ratio",
                    "expected_conserved", "expected_gained", "expected_lost") %in% names(res)))
})

test_that("compute_gain_loss_iterations indices are non-negative", {
  rho_A <- make_rho_matrix(rep(0.4, 20), n_iter = 60)
  rho_B <- make_rho_matrix(rep(0.6, 20), n_iter = 60, seed = 11)
  res   <- compute_gain_loss_iterations(rho_A, rho_B)
  expect_gte(res$gain_index, 0)
  expect_gte(res$loss_index, 0)
})

test_that("compute_gain_loss_iterations gain_index dominates for gain scenario", {
  mat_A <- matrix(0L, nrow = 30, ncol = 100)
  mat_B <- matrix(1L, nrow = 30, ncol = 100)
  res   <- compute_gain_loss_iterations(mat_A, mat_B)
  expect_gt(res$gain_index, res$loss_index)
})

test_that("compute_gain_loss_iterations loss_index dominates for loss scenario", {
  mat_A <- matrix(1L, nrow = 30, ncol = 100)
  mat_B <- matrix(0L, nrow = 30, ncol = 100)
  res   <- compute_gain_loss_iterations(mat_A, mat_B)
  expect_gt(res$loss_index, res$gain_index)
})
