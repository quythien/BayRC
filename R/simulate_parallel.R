################################################################################
# PARALLEL BOOTSTRAP VS DELTA METHOD COMPARISON (LINUX VERSION)
################################################################################


################################################################################
# HELPER FUNCTIONS
################################################################################

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

calculate_congruence <- function(p_A, p_B) {
  sum(p_A * p_B) / sum(p_A + p_B - p_A * p_B)
}

permutation_test <- function(p_A, p_B, B = 1000) {
  n <- length(p_A)
  obs <- calculate_congruence(p_A, p_B)
  nulls <- replicate(B, calculate_congruence(p_A, p_B[sample(n)]))
  (sum(nulls >= obs) + 1) / (B + 1)
}

################################################################################
# BOOTSTRAP VARIANCE ESTIMATION
################################################################################

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

gaussian_test_bootstrap <- function(p_A, p_B, n_boot = 500) {
  obs <- calculate_congruence(p_A, p_B)
  params <- bootstrap_permutation_variance(p_A, p_B, n_boot)
  
  if (params$var_ratio <= 0 || !is.finite(params$var_ratio)) return(NA)
  
  z <- (obs - params$expected_ratio) / params$sd_ratio
  pnorm(z, lower.tail = FALSE)
}

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
