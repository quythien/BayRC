# tests/testthat/test-global.R
# Tests for p_conservation_global() in R/global.R

# Helper: build a minimal rho list
make_rho_list <- function(mat) list(rho = mat)

# Helper: build a binary rho matrix
make_rho_matrix <- function(probs, n_iter = 100, seed = 1) {
  set.seed(seed)
  G <- length(probs)
  matrix(rbinom(G * n_iter, 1, rep(probs, n_iter)), nrow = G)
}

# ---- p_conservation_global --------------------------------------------------

test_that("p_conservation_global returns correct structure", {
  set.seed(1)
  G       <- 20
  n_iter  <- 50
  rho_A   <- make_rho_list(make_rho_matrix(rep(0.7, G), n_iter = n_iter))
  rho_B   <- make_rho_list(make_rho_matrix(rep(0.7, G), n_iter = n_iter, seed = 2))

  obs     <- congruence(rho_A, rho_B)
  perm    <- perm_conservation_global(rho_A, rho_B, B = 20,
                                      parallel = FALSE, delta = 3, units = "hours")
  result  <- p_conservation_global(obs, perm)

  expect_true(all(c("congruence", "gain_index", "loss_index", "ratio") %in% names(result)))
})

test_that("p_conservation_global congruence p-value is in (0, 1]", {
  set.seed(2)
  G      <- 20
  rho_A  <- make_rho_list(make_rho_matrix(rep(0.6, G), n_iter = 50))
  rho_B  <- make_rho_list(make_rho_matrix(rep(0.6, G), n_iter = 50, seed = 3))

  obs    <- congruence(rho_A, rho_B)
  perm   <- perm_conservation_global(rho_A, rho_B, B = 20,
                                     parallel = FALSE, delta = 3, units = "hours")
  result <- p_conservation_global(obs, perm)
  expect_gte(result$congruence, 0)
  expect_lte(result$congruence, 1)
})

test_that("p_conservation_global p-values are NA when obs is NA", {
  obs_na <- list(
    congruence_index = NA,
    gain_loss_ratio  = NA,
    gain_index       = 0,
    loss_index       = 0
  )
  perm <- matrix(runif(20), ncol = 2,
                 dimnames = list(NULL, c("congruence_index", "gain_loss_ratio")))
  result <- p_conservation_global(obs_na, perm)
  expect_true(is.na(result$congruence))
})

test_that("p_conservation_global ratio p-values satisfy right+left >= two_sided logic", {
  set.seed(5)
  G      <- 30
  rho_A  <- make_rho_list(make_rho_matrix(rep(0.5, G), n_iter = 80))
  rho_B  <- make_rho_list(make_rho_matrix(rep(0.8, G), n_iter = 80, seed = 6))

  obs    <- congruence(rho_A, rho_B)
  perm   <- perm_conservation_global(rho_A, rho_B, B = 30,
                                     parallel = FALSE, delta = 3, units = "hours")
  result <- p_conservation_global(obs, perm)

  if (!is.na(result$ratio$PValue_TwoSided)) {
    # two-sided p = 2*min(right, left), capped at 1
    min_one_sided <- min(result$ratio$PValue_Right, result$ratio$PValue_Left)
    expect_lte(result$ratio$PValue_TwoSided, min(1, 2 * min_one_sided + 1e-10))
  }
})
