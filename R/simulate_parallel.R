################################################################################
# PARALLEL BOOTSTRAP VS DELTA METHOD COMPARISON (LINUX VERSION)
################################################################################

library(reshape2)
library(moments)
library(ggplot2)
library(parallel)

set.seed(123)

# Detect number of cores
n_cores <- 40 # Leave one core free
cat("Detected", detectCores(), "cores. Using", n_cores, "cores for parallel processing.\n\n")

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

################################################################################
# SIMULATION PARAMETERS
################################################################################

pathway_sizes <- c(5, 10, 15, 30, 50, 100)
conservation_levels <- c(0, 0.3, 0.5, 0.8, 1)
skewness_levels <- c("symmetric", "moderate_rhythmic", "high_rhythmic", "low_rhythmic")
skewness_labels <- c("Symmetric Beta(2,2)", 
                     "Moderate Rhythmic Beta(3,2)", 
                     "High Rhythmic Beta(5,3)", 
                     "Low Rhythmic Beta(2,5)")

n_sim <- 1000
B_perm <- 1000
B_boot <- 500

################################################################################
# RUN PARALLEL SIMULATIONS
################################################################################

cat("=== PARALLEL BOOTSTRAP VS DELTA METHOD COMPARISON ===\n")
cat("Simulations per condition:", n_sim, "\n")
cat("Bootstrap samples:", B_boot, "\n")
cat("Permutation test samples:", B_perm, "\n")
cat("Using", n_cores, "parallel cores\n\n")

results_list <- list()
variance_comparison_list <- list()
timing_list <- list()

idx <- 1
start_time <- Sys.time()

for (skew in skewness_levels) {
  for (n in pathway_sizes) {
    for (cons in conservation_levels) {
      cat(sprintf("Running: %-20s n=%3d cons=%.1f ... ", skew, n, cons))
      flush.console()
      
      condition_start <- Sys.time()
      
      # PARALLEL EXECUTION
      sim_results <- mclapply(1:n_sim, function(i) {
        run_single_simulation(i, n, cons, skew, B_perm, B_boot)
      }, mc.cores = n_cores, mc.set.seed = TRUE)
      
      condition_time <- difftime(Sys.time(), condition_start, units = "secs")
      cat(sprintf("Done (%.1fs)\n", condition_time))
      
      # Extract results
      p_perm <- sapply(sim_results, `[[`, "p_perm")
      p_gauss_boot <- sapply(sim_results, `[[`, "p_gauss_boot")
      p_edge_boot <- sapply(sim_results, `[[`, "p_edge_boot")
      p_gauss_delta <- sapply(sim_results, `[[`, "p_gauss_delta")
      p_edge_delta <- sapply(sim_results, `[[`, "p_edge_delta")
      
      # Timing
      timing_list[[idx]] <- data.frame(
        skew = skew, n = n, cons = cons,
        time_gauss_boot = mean(sapply(sim_results, `[[`, "time_gauss_boot")),
        time_edge_boot = mean(sapply(sim_results, `[[`, "time_edge_boot")),
        time_gauss_delta = mean(sapply(sim_results, `[[`, "time_gauss_delta")),
        time_edge_delta = mean(sapply(sim_results, `[[`, "time_edge_delta"))
      )
      
      # Variance comparison
      var_comps <- lapply(sim_results, `[[`, "var_comp")
      variance_comparison_list[[idx]] <- data.frame(
        skew = skew, n = n, cons = cons,
        var_boot = sapply(var_comps, `[[`, "var_boot"),
        var_delta = sapply(var_comps, `[[`, "var_delta"),
        ratio = sapply(var_comps, `[[`, "ratio"),
        abs_diff = sapply(var_comps, `[[`, "abs_diff"),
        rel_diff = sapply(var_comps, `[[`, "rel_diff")
      )
      
      results_list[[idx]] <- list(
        skew = skew, n = n, cons = cons,
        p_perm = p_perm,
        p_gauss_boot = p_gauss_boot,
        p_edge_boot = p_edge_boot,
        p_gauss_delta = p_gauss_delta,
        p_edge_delta = p_edge_delta
      )
      idx <- idx + 1
    }
  }
}

elapsed <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("\n=== SIMULATION COMPLETE ===\n"))
cat(sprintf("Total time: %.1f minutes\n", elapsed))
cat(sprintf("Speedup estimate: %.1f× faster than sequential\n\n", n_cores * 0.8))

################################################################################
# ANALYSIS
################################################################################

timing_df <- do.call(rbind, timing_list)
variance_df <- do.call(rbind, variance_comparison_list)

summary_df <- data.frame()

for (res in results_list) {
  valid <- !is.na(res$p_gauss_boot) & !is.na(res$p_edge_boot) & 
    !is.na(res$p_gauss_delta) & !is.na(res$p_edge_delta)
  
  p_perm <- res$p_perm[valid]
  p_gauss_boot <- res$p_gauss_boot[valid]
  p_edge_boot <- res$p_edge_boot[valid]
  p_gauss_delta <- res$p_gauss_delta[valid]
  p_edge_delta <- res$p_edge_delta[valid]
  
  if (length(p_perm) == 0) next
  
  if (res$cons == 0) {
    type1_perm <- mean(p_perm < 0.05)
    type1_gauss_boot <- mean(p_gauss_boot < 0.05)
    type1_edge_boot <- mean(p_edge_boot < 0.05)
    type1_gauss_delta <- mean(p_gauss_delta < 0.05)
    type1_edge_delta <- mean(p_edge_delta < 0.05)
    power_perm <- power_gauss_boot <- power_edge_boot <- NA
    power_gauss_delta <- power_edge_delta <- NA
  } else {
    type1_perm <- type1_gauss_boot <- type1_edge_boot <- NA
    type1_gauss_delta <- type1_edge_delta <- NA
    power_perm <- mean(p_perm < 0.05)
    power_gauss_boot <- mean(p_gauss_boot < 0.05)
    power_edge_boot <- mean(p_edge_boot < 0.05)
    power_gauss_delta <- mean(p_gauss_delta < 0.05)
    power_edge_delta <- mean(p_edge_delta < 0.05)
  }
  
  mae_gauss_boot <- mean(abs(p_perm - p_gauss_boot))
  mae_edge_boot <- mean(abs(p_perm - p_edge_boot))
  mae_gauss_delta <- mean(abs(p_perm - p_gauss_delta))
  mae_edge_delta <- mean(abs(p_perm - p_edge_delta))
  
  summary_df <- rbind(summary_df, data.frame(
    skew = as.character(res$skew),
    n = as.numeric(res$n),
    cons = as.numeric(res$cons),
    type1_perm = type1_perm,
    type1_gauss_boot = type1_gauss_boot,
    type1_edge_boot = type1_edge_boot,
    type1_gauss_delta = type1_gauss_delta,
    type1_edge_delta = type1_edge_delta,
    power_perm = power_perm,
    power_gauss_boot = power_gauss_boot,
    power_edge_boot = power_edge_boot,
    power_gauss_delta = power_gauss_delta,
    power_edge_delta = power_edge_delta,
    mae_gauss_boot = mae_gauss_boot,
    mae_edge_boot = mae_edge_boot,
    mae_gauss_delta = mae_gauss_delta,
    mae_edge_delta = mae_edge_delta,
    stringsAsFactors = FALSE
  ))
}

summary_df <- merge(summary_df, timing_df, by = c("skew", "n", "cons"))

################################################################################
# COMPREHENSIVE SUMMARY
################################################################################

cat("\n========================================\n")
cat("===  TYPE I ERROR RATE COMPARISON  ===\n")
cat("========================================\n")

type1_df <- summary_df[summary_df$cons == 0, ]

cat("\nMean Type I Error (Target = 0.05):\n")
cat(sprintf("  Permutation:        %.4f\n", mean(type1_df$type1_perm, na.rm = TRUE)))
cat(sprintf("  Gauss (Bootstrap):  %.4f  [Deviation: %.4f]\n", 
            mean(type1_df$type1_gauss_boot, na.rm = TRUE),
            mean(abs(type1_df$type1_gauss_boot - 0.05), na.rm = TRUE)))
cat(sprintf("  Edge (Bootstrap):   %.4f  [Deviation: %.4f]\n", 
            mean(type1_df$type1_edge_boot, na.rm = TRUE),
            mean(abs(type1_df$type1_edge_boot - 0.05), na.rm = TRUE)))
cat(sprintf("  Gauss (Delta):      %.4f  [Deviation: %.4f]\n", 
            mean(type1_df$type1_gauss_delta, na.rm = TRUE),
            mean(abs(type1_df$type1_gauss_delta - 0.05), na.rm = TRUE)))
cat(sprintf("  Edge (Delta):       %.4f  [Deviation: %.4f]\n", 
            mean(type1_df$type1_edge_delta, na.rm = TRUE),
            mean(abs(type1_df$type1_edge_delta - 0.05), na.rm = TRUE)))

cat("\n========================================\n")
cat("===  COMPUTATIONAL EFFICIENCY  ===\n")
cat("========================================\n")

cat("\nAverage Computation Time (seconds per test):\n")
cat(sprintf("  Gauss (Bootstrap):  %.4f\n", mean(timing_df$time_gauss_boot)))
cat(sprintf("  Edge (Bootstrap):   %.4f\n", mean(timing_df$time_edge_boot)))
cat(sprintf("  Gauss (Delta):      %.4f\n", mean(timing_df$time_gauss_delta)))
cat(sprintf("  Edge (Delta):       %.4f\n", mean(timing_df$time_edge_delta)))

cat("\nSpeedup Factor (Bootstrap / Delta):\n")
cat(sprintf("  Gaussian:  %.1f×\n", mean(timing_df$time_gauss_boot) / mean(timing_df$time_gauss_delta)))
cat(sprintf("  Edgeworth: %.1f×\n", mean(timing_df$time_edge_boot) / mean(timing_df$time_edge_delta)))

cat("\n========================================\n")
cat("===  VARIANCE ESTIMATE COMPARISON  ===\n")
cat("========================================\n")

cat("\nVariance Ratio (Bootstrap / Delta):\n")
cat(sprintf("  Mean:   %.4f\n", mean(variance_df$ratio, na.rm = TRUE)))
cat(sprintf("  Median: %.4f\n", median(variance_df$ratio, na.rm = TRUE)))
cat(sprintf("  SD:     %.4f\n", sd(variance_df$ratio, na.rm = TRUE)))

cat("\nRelative Difference (|Boot - Delta| / Delta):\n")
cat(sprintf("  Mean:   %.4f (%.1f%%)\n", 
            mean(variance_df$rel_diff, na.rm = TRUE),
            mean(variance_df$rel_diff, na.rm = TRUE) * 100))
cat(sprintf("  Median: %.4f (%.1f%%)\n", 
            median(variance_df$rel_diff, na.rm = TRUE),
            median(variance_df$rel_diff, na.rm = TRUE) * 100))

cat("\n========================================\n")
cat("===  RECOMMENDATION  ===\n")
cat("========================================\n")

type1_dev_boot <- mean(abs(type1_df$type1_edge_boot - 0.05), na.rm = TRUE)
type1_dev_delta <- mean(abs(type1_df$type1_edge_delta - 0.05), na.rm = TRUE)
speedup <- mean(timing_df$time_edge_boot) / mean(timing_df$time_edge_delta)

cat("\nBased on simulation results:\n\n")
if (type1_dev_delta <= 0.01 && type1_dev_boot <= 0.01) {
  cat("✓ BOTH methods provide adequate Type I error control\n")
  cat(sprintf("✓ Delta method is %.0f× FASTER\n", speedup))
  cat("\n>>> RECOMMENDATION: Use DELTA METHOD <<<\n")
  cat("    (Faster, equally accurate)\n")
} else if (type1_dev_boot < type1_dev_delta * 0.8) {
  cat("✓ Bootstrap provides significantly better Type I error control\n")
  cat(sprintf("  (Bootstrap deviation: %.4f vs Delta: %.4f)\n", type1_dev_boot, type1_dev_delta))
  cat("\n>>> RECOMMENDATION: Use BOOTSTRAP METHOD <<<\n")
  cat("    (Better statistical properties justify computational cost)\n")
} else {
  cat("✓ Methods have similar Type I error control\n")
  cat(sprintf("✓ Delta method is %.0f× FASTER\n", speedup))
  cat("\n>>> RECOMMENDATION: Use DELTA METHOD for most applications <<<\n")
  cat("    Consider bootstrap only for small pathways (n < 15) with high skewness\n")
}

################################################################################
# SAVE RESULTS
################################################################################

save(summary_df, timing_df, variance_df, 
     file = "parallel_bootstrap_vs_delta_results.RData")

cat("\n========================================\n")
cat("Results saved to: parallel_bootstrap_vs_delta_results.RData\n")
cat("Load later with: load('parallel_bootstrap_vs_delta_results.RData')\n")
cat("========================================\n")