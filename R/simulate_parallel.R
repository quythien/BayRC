################################################################################
# PARALLEL BOOTSTRAP VS DELTA METHOD COMPARISON (LINUX VERSION)
################################################################################


################################################################################
# HELPER FUNCTIONS
################################################################################

#' Generate synthetic posterior rhythmicity probability vectors
#'
#' @description
#' Simulates posterior rhythmicity probabilities (p_A, p_B) for a given
#' number of genes under a specified conservation level and skewness of the
#' marginal Beta distribution.  Used for simulation studies comparing
#' variance estimation methods.
#'
#' @param n Integer; number of genes.
#' @param conservation Numeric in \[0, 1\]; correlation between p_A and
#'   p_B; higher values simulate more conserved rhythmicity (default 0.5).
#' @param skewness Character; shape of the p_A Beta distribution:
#'   \code{"symmetric"} (Beta(2,2)), \code{"moderate_rhythmic"}
#'   (Beta(3,2)), \code{"high_rhythmic"} (Beta(5,3)), or
#'   \code{"low_rhythmic"} (Beta(2,5)).
#'
#' @return A data.frame with columns \code{p_A} and \code{p_B} (length n).
#'
#' @export
generate_gene_data_full_spectrum <- function(n, conservation = 0.5, skewness = "symmetric") {
  if (skewness == "symmetric") {
    p_A <- rbeta(n, 2, 2)
  } else if (skewness == "moderate_rhythmic") {
    p_A <- rbeta(n, 3, 2)
  } else if (skewness == "high_rhythmic") {
    p_A <- rbeta(n, 5, 3)
  } else if (skewness == "low_rhythmic") {
    p_A <- rbeta(n, 2, 5)
  } else {
    stop("Unknown skewness option")
  }
  
  p_A <- p_A * 0.9 + 0.05
  
  if (conservation > 0) {
    noise_sd <- 0.15
    lower <- 0.1; upper <- 0.9
    p_B <- pmin(pmax(p_A + (1 - conservation) * rnorm(n, 0, noise_sd), lower), upper)
  } else {
    p_B <- sample(p_A) 
  }
  list(p_A = p_A, p_B = p_B)
}

#' Compute probabilistic congruence index from probability vectors
#'
#' @description
#' Computes the fuzzy Jaccard congruence index
#' sum(p_A * p_B) / sum(p_A + p_B - p_A * p_B) from two numeric vectors
#' of posterior rhythmicity probabilities.
#'
#' @param p_A,p_B Numeric vectors of posterior rhythmicity probabilities
#'   (length G, values in \[0, 1\]).
#'
#' @return Numeric scalar; congruence index in \[0, 1\].
#'
#' @export
calculate_congruence <- function(p_A, p_B) {
  sum(p_A * p_B) / sum(p_A + p_B - p_A * p_B)
}

#' Permutation test for congruence index significance
#'
#' @description
#' Computes the observed congruence index between p_A and p_B, then
#' builds a null distribution by shuffling p_B labels B times and
#' returns a one-sided p-value.
#'
#' @param p_A,p_B Numeric vectors; posterior rhythmicity probabilities.
#' @param B Integer; number of permutations (default 1000).
#'
#' @return A list with \code{observed}, \code{p_value}, and \code{null_dist}.
#'
#' @export
permutation_test <- function(p_A, p_B, B = 1000) {
  n <- length(p_A)
  obs <- calculate_congruence(p_A, p_B)
  nulls <- replicate(B, calculate_congruence(p_A, p_B[sample(n)]))
  (sum(nulls >= obs) + 1) / (B + 1)
}

################################################################################
# BOOTSTRAP VARIANCE ESTIMATION
################################################################################

#' Bootstrap estimate of permutation variance for the congruence index
#'
#' @description
#' Estimates the variance of the permutation null distribution of the
#' congruence index by bootstrapping, for use in analytical approximation
#' testing methods.
#'
#' @param p_A,p_B Numeric vectors; posterior rhythmicity probabilities.
#' @param n_boot Integer; bootstrap replicates (default 500).
#'
#' @return A list with \code{variance} (bootstrap estimate of the
#'   permutation variance) and \code{se} (standard error).
#'
#' @export
bootstrap_permutation_variance <- function(p_A, p_B, n_boot = 500) {
  n <- length(p_A)
  
  sum_pA <- sum(p_A)
  sum_pB <- sum(p_B)
  mean_pB <- mean(p_B)
  
  mu_AB_null <- sum_pA * mean_pB
  mu_Union_null <- sum_pA + sum_pB - mu_AB_null
  expected_ratio <- mu_AB_null / mu_Union_null
  
  boot_ratios <- replicate(n_boot, {
    perm_idx <- sample(n)
    AB_sum <- sum(p_A * p_B[perm_idx])
    Union_sum <- sum(p_A + p_B[perm_idx] - p_A * p_B[perm_idx])
    AB_sum / Union_sum
  })
  
  var_ratio <- var(boot_ratios)
  
  list(
    expected_ratio = expected_ratio,
    var_ratio = var_ratio,
    sd_ratio = sqrt(var_ratio),
    boot_ratios = boot_ratios
  )
}

################################################################################
# DELTA METHOD VARIANCE
################################################################################

#' Delta-method estimate of permutation variance for the congruence index
#'
#' @description
#' Uses the delta method to analytically approximate the variance of the
#' congruence index under gene-label permutation of p_B.
#'
#' @param p_A,p_B Numeric vectors; posterior rhythmicity probabilities.
#'
#' @return A list with \code{variance} (delta-method estimate) and
#'   \code{gradient} (partial derivatives for diagnostics).
#'
#' @export
delta_variance <- function(p_A, p_B) {
  n <- length(p_A)
  AB_vec <- p_A * p_B
  Union_vec <- p_A + p_B - AB_vec
  
  mu_AB <- mean(AB_vec)
  mu_Union <- mean(Union_vec)
  
  grad <- c(1/mu_Union, -mu_AB/(mu_Union^2))
  Sigma <- cov(cbind(AB_vec, Union_vec))
  
  var_T <- as.numeric((1/n) * t(grad) %*% Sigma %*% grad)
  
  list(expected_ratio = mu_AB/mu_Union,
       var_ratio = var_T,
       sd_ratio = sqrt(var_T))
}

################################################################################
# GAUSSIAN TESTS
################################################################################

#' Gaussian significance test using bootstrap variance estimate
#'
#' @description
#' Tests whether the congruence index differs significantly from the
#' null expectation using a Gaussian approximation with bootstrap-estimated
#' variance.
#'
#' @param p_A,p_B Numeric vectors; posterior rhythmicity probabilities.
#' @param n_boot Integer; bootstrap replicates for variance estimation.
#'
#' @return A list with \code{z_score}, \code{p_value}, and
#'   \code{variance_est}.
#'
#' @export
gaussian_test_bootstrap <- function(p_A, p_B, n_boot = 500) {
  obs <- calculate_congruence(p_A, p_B)
  params <- bootstrap_permutation_variance(p_A, p_B, n_boot)
  
  if (params$var_ratio <= 0 || !is.finite(params$var_ratio)) return(NA)
  
  z <- (obs - params$expected_ratio) / params$sd_ratio
  pnorm(z, lower.tail = FALSE)
}

#' Gaussian significance test using delta-method variance estimate
#'
#' @description
#' Tests congruence index significance using a Gaussian approximation with
#' delta-method variance, avoiding the computational cost of bootstrapping.
#'
#' @param p_A,p_B Numeric vectors; posterior rhythmicity probabilities.
#'
#' @return A list with \code{z_score}, \code{p_value}, and
#'   \code{variance_est}.
#'
#' @export
gaussian_test_delta <- function(p_A, p_B) {
  obs <- calculate_congruence(p_A, p_B)
  params <- delta_variance(p_A, p_B)
  
  if (params$var_ratio <= 0 || !is.finite(params$var_ratio)) return(NA)
  
  z <- (obs - params$expected_ratio) / params$sd_ratio
  pnorm(z, lower.tail = FALSE)
}

################################################################################
# EDGEWORTH CORRECTION
################################################################################

edgeworth_correction_fixed <- function(z, skew, kurt, n) {
  phi <- dnorm(z); Phi <- pnorm(z)
  term1 <- (skew/(6*sqrt(n)))*(z^2-1)*phi
  term2 <- (kurt/(24*n))*(z^3-3*z)*phi
  term3 <- (skew^2/(72*n))*(z^5-10*z^3+15*z)*phi
  Phi_corrected <- Phi - term1 - term2 - term3
  pmax(0.0001, pmin(0.9999, Phi_corrected))
}

#' Edgeworth expansion significance test using bootstrap variance
#'
#' @description
#' Applies an Edgeworth series correction to the Gaussian approximation to
#' better handle skewness in the permutation null distribution of the
#' congruence index, using bootstrap-estimated variance and skewness.
#'
#' @param p_A,p_B Numeric vectors; posterior rhythmicity probabilities.
#' @param n_boot Integer; bootstrap replicates.
#'
#' @return A list with \code{z_score}, \code{p_value_gaussian},
#'   \code{p_value_edgeworth}, and \code{variance_est}.
#'
#' @export
edgeworth_test_bootstrap <- function(p_A, p_B, n_boot = 500) {
  n <- length(p_A)
  obs <- calculate_congruence(p_A, p_B)
  params <- bootstrap_permutation_variance(p_A, p_B, n_boot)
  
  if (params$var_ratio <= 0 || !is.finite(params$var_ratio)) return(NA)
  
  z <- (obs - params$expected_ratio) / params$sd_ratio
  
  if (n <= 30 && length(params$boot_ratios) >= 100) {
    ratio_std <- (params$boot_ratios - mean(params$boot_ratios)) / sd(params$boot_ratios)
    skew <- mean(ratio_std^3)
    kurt <- mean(ratio_std^4) - 3
    
    if (abs(skew) > 0.3) {
      return(1 - edgeworth_correction_fixed(z, skew, kurt, n))
    }
  }
  
  pnorm(z, lower.tail = FALSE)
}

#' Edgeworth expansion significance test using delta-method variance
#'
#' @description
#' Applies the Edgeworth correction using delta-method variance and skewness
#' estimates, avoiding the need for bootstrap replicates.
#'
#' @param p_A,p_B Numeric vectors; posterior rhythmicity probabilities.
#'
#' @return A list with \code{z_score}, \code{p_value_gaussian},
#'   \code{p_value_edgeworth}, and \code{variance_est}.
#'
#' @export
edgeworth_test_delta <- function(p_A, p_B) {
  n <- length(p_A)
  obs <- calculate_congruence(p_A, p_B)
  params <- delta_variance(p_A, p_B)
  
  if (params$var_ratio <= 0 || !is.finite(params$var_ratio)) return(NA)
  
  z <- (obs - params$expected_ratio) / params$sd_ratio
  
  ratios <- p_A * p_B / (p_A + p_B - p_A * p_B)
  ratio_std <- (ratios - mean(ratios)) / sd(ratios)
  skew <- mean(ratio_std^3)
  kurt <- mean(ratio_std^4) - 3
  
  if (n <= 30 && abs(skew) > 0.3) {
    return(1 - edgeworth_correction_fixed(z, skew, kurt, n))
  }
  
  pnorm(z, lower.tail = FALSE)
}

################################################################################
# VARIANCE COMPARISON FUNCTION
################################################################################

#' Compare all variance estimation methods for the congruence index
#'
#' @description
#' Applies all four variance estimation approaches (bootstrap Gaussian,
#' delta Gaussian, bootstrap Edgeworth, delta Edgeworth) and returns a
#' data.frame comparing p-values, z-scores, and variance estimates.
#'
#' @param p_A,p_B Numeric vectors; posterior rhythmicity probabilities.
#' @param n_boot Integer; bootstrap replicates (default 500).
#'
#' @return Data.frame with one row per method and columns for method name,
#'   z-score, p-value, and variance.
#'
#' @export
compare_variances <- function(p_A, p_B, n_boot = 500) {
  boot_params <- bootstrap_permutation_variance(p_A, p_B, n_boot)
  delta_params <- delta_variance(p_A, p_B)
  
  list(
    var_boot = boot_params$var_ratio,
    var_delta = delta_params$var_ratio,
    ratio = boot_params$var_ratio / delta_params$var_ratio,
    abs_diff = abs(boot_params$var_ratio - delta_params$var_ratio),
    rel_diff = abs(boot_params$var_ratio - delta_params$var_ratio) / delta_params$var_ratio
  )
}

################################################################################
# SINGLE SIMULATION RUN (for parallelization)
################################################################################

#' Run a single simulation replicate for method comparison
#'
#' @description
#' Generates synthetic p_A/p_B data via \code{generate_gene_data_full_spectrum}
#' and runs the permutation test and all four variance-based tests,
#' returning a data.frame of results for one simulation replicate.  Designed
#' for parallelisation over many replicates.
#'
#' @param sim_id Integer; replicate identifier (used as random seed offset).
#' @param n Integer; number of genes.
#' @param cons Numeric; conservation level (passed to
#'   \code{generate_gene_data_full_spectrum}).
#' @param skew Character; skewness label.
#' @param B_perm Integer; permutations for the permutation test.
#' @param B_boot Integer; bootstrap replicates for variance estimation.
#'
#' @return Data.frame with columns for sim_id, method, z_score, p_value,
#'   observed_congruence, and variance.
#'
#' @export
run_single_simulation <- function(sim_id, n, cons, skew, B_perm, B_boot) {
  dat <- generate_gene_data_full_spectrum(n, cons, skew)
  
  # Time bootstrap methods
  t1 <- system.time({
    p_gauss_boot <- gaussian_test_bootstrap(dat$p_A, dat$p_B, n_boot = B_boot)
  })
  
  t2 <- system.time({
    p_edge_boot <- edgeworth_test_bootstrap(dat$p_A, dat$p_B, n_boot = B_boot)
  })
  
  # Time delta methods
  t3 <- system.time({
    p_gauss_delta <- gaussian_test_delta(dat$p_A, dat$p_B)
  })
  
  t4 <- system.time({
    p_edge_delta <- edgeworth_test_delta(dat$p_A, dat$p_B)
  })
  
  # Permutation test (gold standard)
  p_perm <- permutation_test(dat$p_A, dat$p_B, B = B_perm)
  
  # Variance comparison
  var_comp <- compare_variances(dat$p_A, dat$p_B, n_boot = B_boot)
  
  list(
    p_perm = p_perm,
    p_gauss_boot = p_gauss_boot,
    p_edge_boot = p_edge_boot,
    p_gauss_delta = p_gauss_delta,
    p_edge_delta = p_edge_delta,
    time_gauss_boot = t1[3],
    time_edge_boot = t2[3],
    time_gauss_delta = t3[3],
    time_edge_delta = t4[3],
    var_comp = var_comp
  )
}
