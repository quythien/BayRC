####################################### INTRA-COHORT ANALYSIS ###################
#── Clear environment ───────────────────────────────────────────────────────────
rm(list=ls())

#── Paths ───────────────────────────────────────────────────────────────────────
current_gtex <- "/home/qtp1/Projects/Collaborative"
current_wd   <- "/home/qtp1/Projects/Circadian"
current_aging <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging"

#── Load data ───────────────────────────────────────────────────────────────────
BA11 <- readRDS(file.path(current_aging, "data/BA11_data.rds"))
BA47 <- readRDS(file.path(current_aging, "data/BA47_data.rds"))

#── Load packages ───────────────────────────────────────────────────────────────
library(parallel)
library(edgeR)
library(Rcpp)
library(dplyr)

#── Source dependencies ─────────────────────────────────────────────────────────
WD <- "Kyle/Circadian-analysis-main/R/v1"
setwd(file.path(current_wd, WD))
scripts <- list.files("R", full.names=TRUE)
sapply(scripts, source)

#══════════════════════════════════════════════════════════════════════════════
# CORE FUNCTIONS
#══════════════════════════════════════════════════════════════════════════════

#── Run MCMC and get posterior probabilities ────────────────────────────────────
run_MCMC_get_posteriors <- function(data, tod, thin = 1, n.iter = 2500, 
                                    n.burn = 500, P = 24, seed = 1) {
  
  # Store original gene names
  original_genes <- rownames(data)
  
  if (is.null(original_genes)) {
    stop("Input data MUST have rownames (gene identifiers)!")
  }
  
  if (any(duplicated(original_genes))) {
    warning("Duplicate gene names detected! Making unique...")
    original_genes <- make.unique(original_genes)
    rownames(data) <- original_genes
  }
  
  # Prepare data list
  dat_list <- list(data = as.data.frame(data),
                   time = tod,
                   gname = original_genes)
  
  # Initialize MCMC
  a.init <- CBt_init_single(Data.list = dat_list, P = P,
                            FitCosinor = TRUE,
                            mu_M = 0, sigma_M = 10,
                            mu_A = 1, sigma_A = 10,
                            seed = seed)
  
  # Run MCMC
  CB.res <- CB_MCMC_single_rj_slice(Data.list = dat_list,
                                    Init.value = a.init,
                                    P = P, iteration = n.iter,
                                    thin = thin, n.burn = n.burn,
                                    seed = seed,
                                    p_rhythmic = rep(0.2, nrow(data)),
                                    rj.p.stay = 0.5,
                                    A_prior = "trunc_Normal_OLS_condi",
                                    mu_A = 1, sigma_A = 10^2, A.min = 0,
                                    A_wb_beta2 = 2,
                                    A_gm_shape = 1.99, A_gm_rate = 0.5,
                                    rj.phi = TRUE, rj.A = TRUE,
                                    mu_M = 0, sigma_M = 10^2,
                                    sigma_prior_v = 2, sigma_prior_s = 0)
  
  # Ensure gene order is preserved
  if (is.null(rownames(CB.res$rho))) {
    rownames(CB.res$rho) <- original_genes
  }
  CB.res$rho <- CB.res$rho[original_genes, , drop = FALSE]
  
  # Calculate posteriors
  posteriors <- rowMeans(CB.res$rho)
  names(posteriors) <- original_genes
  
  # Verification
  if (!identical(names(posteriors), original_genes)) {
    stop("CRITICAL ERROR: Gene order was not preserved by MCMC!")
  }
  
  return(posteriors)
}

#── Null-adjusted Spearman correlation ──────────────────────────────────────────
spearman <- function(p_A, p_B, n_perm = 1000) {
  
  # Verify gene names match
  if (!identical(names(p_A), names(p_B))) {
    stop("Gene names must match between p_A and p_B!")
  }
  
  # Observed correlation
  obs <- cor(p_A, p_B, method = "spearman")
  
  # Permutation null distribution
  null_dist <- sapply(1:n_perm, function(i) {
    cor(p_A, sample(p_B), method = "spearman")
  })
  
  # Calculate null statistics
  null_mean <- mean(null_dist, na.rm = TRUE)
  null_sd <- sd(null_dist, na.rm = TRUE)
  
  # Adjusted correlation
  if (abs(1 - null_mean) > 0.001) {
    adjusted <- (obs - null_mean) / (1 - null_mean)
    null_adjusted <- (null_dist - null_mean) / (1 - null_mean)
  } else {
    adjusted <- obs
    null_adjusted <- null_dist
  }
  
  # Empirical p-value (one-sided)
  p_value <- (1 + sum(null_adjusted >= adjusted, na.rm = TRUE)) / (1 + n_perm)
  
  return(list(
    observed = obs,
    null_mean = null_mean,
    null_sd = null_sd,
    adjusted = adjusted,
    p_value = p_value,
    z_score = (obs - null_mean) / null_sd
  ))
}

#── Random split function ───────────────────────────────────────────────────────
random_split_dataset <- function(expr_data, tod_vector, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n_samples <- ncol(expr_data)
  if (length(tod_vector) != n_samples) {
    stop("Expression data and TOD vector must have same length!")
  }
  
  if (is.null(rownames(expr_data))) {
    stop("Expression data must have rownames (gene names)!")
  }
  
  # Random split
  shuffled_indices <- sample(1:n_samples)
  split_point <- floor(n_samples / 2)
  indices_1 <- shuffled_indices[1:split_point]
  indices_2 <- shuffled_indices[(split_point + 1):n_samples]
  
  subset_1 <- list(
    expr = expr_data[, indices_1, drop = FALSE],
    tod = tod_vector[indices_1]
  )
  
  subset_2 <- list(
    expr = expr_data[, indices_2, drop = FALSE],
    tod = tod_vector[indices_2]
  )
  
  return(list(subset1 = subset_1, subset2 = subset_2))
}

#── Prepare brain data ──────────────────────────────────────────────────────────
prepare_brain_data <- function(brain_data, region_name) {
  expr_ids <- sub(paste0(".*", region_name, "-([0-9]+)\\.CEL\\.gz"), "\\1", 
                  colnames(brain_data$expr))
  expr_ids <- as.integer(expr_ids)
  names(expr_ids) <- colnames(brain_data$expr)
  
  tod_map <- setNames(brain_data$pheno$TOD, brain_data$pheno$ID)
  
  list(
    all_expr = brain_data$expr,
    all_tod = tod_map[as.character(expr_ids)]
  )
}

#══════════════════════════════════════════════════════════════════════════════
# MAIN ANALYSIS FUNCTION
#══════════════════════════════════════════════════════════════════════════════

run_single_repetition <- function(rep_num, ba11_data, ba47_data, 
                                  n.iter = 2000, n.burn = 500, thin = 2, 
                                  n_perm = 1000, n_cores_per_rep = 8) {  # ← ADDED PARAMETER
  
  cat("\n", rep("=", 70), "\n", sep="")
  cat(sprintf("REPETITION %d (using %d cores)\n", rep_num, n_cores_per_rep))  # ← MODIFIED
  cat(rep("=", 70), "\n\n", sep="")
  
  # Create random splits
  set.seed(rep_num)
  split_test1 <- random_split_dataset(ba11_data$all_expr, ba11_data$all_tod, seed = rep_num)
  split_test2 <- random_split_dataset(ba47_data$all_expr, ba47_data$all_tod, seed = rep_num)
  
  # Define 8 MCMC jobs
  mcmc_jobs <- list(
    # BA11: subset1, seed=1
    list(dataset = "BA11", expr = split_test1$subset1$expr, tod = split_test1$subset1$tod, 
         n.iter = n.iter, n.burn = n.burn, thin = thin, seed = 1),
    # BA11: subset1, seed=10001
    list(dataset = "BA11", expr = split_test1$subset1$expr, tod = split_test1$subset1$tod, 
         n.iter = n.iter, n.burn = n.burn, thin = thin, seed = 10001),
    # BA11: subset2, seed=1
    list(dataset = "BA11", expr = split_test1$subset2$expr, tod = split_test1$subset2$tod, 
         n.iter = n.iter, n.burn = n.burn, thin = thin, seed = 1),
    # BA11: subset2, seed=10001
    list(dataset = "BA11", expr = split_test1$subset2$expr, tod = split_test1$subset2$tod, 
         n.iter = n.iter, n.burn = n.burn, thin = thin, seed = 10001),
    # BA47: subset1, seed=1
    list(dataset = "BA47", expr = split_test2$subset1$expr, tod = split_test2$subset1$tod, 
         n.iter = n.iter, n.burn = n.burn, thin = thin, seed = 1),
    # BA47: subset1, seed=10001
    list(dataset = "BA47", expr = split_test2$subset1$expr, tod = split_test2$subset1$tod, 
         n.iter = n.iter, n.burn = n.burn, thin = thin, seed = 10001),
    # BA47: subset2, seed=1
    list(dataset = "BA47", expr = split_test2$subset2$expr, tod = split_test2$subset2$tod, 
         n.iter = n.iter, n.burn = n.burn, thin = thin, seed = 1),
    # BA47: subset2, seed=10001
    list(dataset = "BA47", expr = split_test2$subset2$expr, tod = split_test2$subset2$tod, 
         n.iter = n.iter, n.burn = n.burn, thin = thin, seed = 10001)
  )
  
  # Run all 8 jobs in parallel
  cat(sprintf("Running 8 MCMC jobs in parallel using %d cores...\n", n_cores_per_rep))  # ← MODIFIED
  start_time <- Sys.time()
  
  results <- mclapply(seq_along(mcmc_jobs), function(idx) {
    job <- mcmc_jobs[[idx]]
    tryCatch({
      cat(sprintf("[Job %d] Starting: %s subset, seed=%d\n", 
                  idx, job$dataset, job$seed))
      
      result <- run_MCMC_get_posteriors(
        data = job$expr,
        tod = job$tod,
        n.iter = job$n.iter,
        n.burn = job$n.burn, 
        thin = job$thin, 
        seed = job$seed
      )
      
      cat(sprintf("[Job %d] ✓ Completed: %d genes\n", idx, length(result)))
      return(result)
      
    }, error = function(e) {
      cat(sprintf("[Job %d] ✗ FAILED with error:\n", idx))
      print(e)
      cat(sprintf("  Dataset: %s, Seed: %d, N samples: %d\n",
                  job$dataset, job$seed, ncol(job$expr)))
      return(structure(paste("Error:", e$message), class = "try-error"))
    })
  }, mc.cores = n_cores_per_rep)
  
  # Check for errors after mclapply
  cat("\n=== Checking for failed MCMC jobs ===\n")
  failed_jobs <- sapply(results, inherits, "try-error")
  if (any(failed_jobs)) {
    cat(sprintf("WARNING: %d jobs failed!\n", sum(failed_jobs)))
    cat(sprintf("Failed job indices: %s\n", paste(which(failed_jobs), collapse=", ")))
    for (i in which(failed_jobs)) {
      cat(sprintf("\nJob %d error:\n", i))
      print(results[[i]])
    }
    stop("Cannot continue with failed MCMC jobs!")
  }
  cat("✓ All MCMC jobs completed successfully\n")
  
  
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("✓ Completed in %.2f minutes\n\n", elapsed))
  
  # Calculate concordances with null adjustment
  cat("Calculating concordances...\n")
  
  # BA11 - MCMC reliability (same data, different seeds)
  ba11_mcmc_s1 <- spearman(results[[1]], results[[2]], n_perm = n_perm)
  ba11_mcmc_s2 <- spearman(results[[3]], results[[4]], n_perm = n_perm)
  
  # BA11 - Biological reproducibility (different data)
  ba11_bio_1 <- spearman(results[[1]], results[[3]], n_perm = n_perm)
  ba11_bio_2 <- spearman(results[[2]], results[[4]], n_perm = n_perm)
  ba11_bio_3 <- spearman(results[[1]], results[[4]], n_perm = n_perm)
  ba11_bio_4 <- spearman(results[[2]], results[[3]], n_perm = n_perm)
  
  # BA47 - MCMC reliability
  ba47_mcmc_s1 <- spearman(results[[5]], results[[6]], n_perm = n_perm)
  ba47_mcmc_s2 <- spearman(results[[7]], results[[8]], n_perm = n_perm)
  
  # BA47 - Biological reproducibility
  ba47_bio_1 <- spearman(results[[5]], results[[7]], n_perm = n_perm)
  ba47_bio_2 <- spearman(results[[6]], results[[8]], n_perm = n_perm)
  ba47_bio_3 <- spearman(results[[5]], results[[8]], n_perm = n_perm)
  ba47_bio_4 <- spearman(results[[6]], results[[7]], n_perm = n_perm)
  
  # Print summary
  cat("\n--- Repetition Summary ---\n")
  cat(sprintf("BA11 MCMC Reliability:         %.4f\n", 
              mean(c(ba11_mcmc_s1$adjusted, ba11_mcmc_s2$adjusted))))
  cat(sprintf("BA11 Biological Reproducibility: %.4f\n", 
              mean(c(ba11_bio_1$adjusted, ba11_bio_2$adjusted, 
                     ba11_bio_3$adjusted, ba11_bio_4$adjusted))))
  cat(sprintf("BA47 MCMC Reliability:         %.4f\n", 
              mean(c(ba47_mcmc_s1$adjusted, ba47_mcmc_s2$adjusted))))
  cat(sprintf("BA47 Biological Reproducibility: %.4f\n", 
              mean(c(ba47_bio_1$adjusted, ba47_bio_2$adjusted, 
                     ba47_bio_3$adjusted, ba47_bio_4$adjusted))))
  
  # Return results
  list(
    repetition = rep_num,
    ba11 = list(
      mcmc_reliability = list(subset1 = ba11_mcmc_s1, subset2 = ba11_mcmc_s2),
      biological = list(seed1 = ba11_bio_1, seed10001 = ba11_bio_2,
                        cross1 = ba11_bio_3, cross2 = ba11_bio_4)
    ),
    ba47 = list(
      mcmc_reliability = list(subset1 = ba47_mcmc_s1, subset2 = ba47_mcmc_s2),
      biological = list(seed1 = ba47_bio_1, seed10001 = ba47_bio_2,
                        cross1 = ba47_bio_3, cross2 = ba47_bio_4)
    )
  )
}

#══════════════════════════════════════════════════════════════════════════════
# RUN ANALYSIS
#══════════════════════════════════════════════════════════════════════════════

# Prepare data
BA11_data <- prepare_brain_data(BA11, "BA11")
BA47_data <- prepare_brain_data(BA47, "BA47")

# Configuration
N_REPS <- 10
N_CORES_TOTAL <- 40  # Total cores available
N_REPS_PARALLEL <- 10  # Number of repetitions to run in parallel
N_CORES_PER_REP <- N_CORES_TOTAL / N_REPS_PARALLEL  
# Run repetitions in parallel
cat("\n", rep("=", 70), "\n", sep="")
cat(sprintf("STARTING %d REPETITIONS\n", N_REPS))
cat(sprintf("Configuration: %d repetitions in parallel, %d cores each\n", 
            N_REPS_PARALLEL, N_CORES_PER_REP))
cat(sprintf("Total cores used: %d\n", N_CORES_TOTAL))
cat(rep("=", 70), "\n", sep="")

overall_start <- Sys.time()

all_repetitions <- mclapply(1:N_REPS, function(i) {
  tryCatch({
    run_single_repetition(
      rep_num = i,
      ba11_data = BA11_data,
      ba47_data = BA47_data,
      n.iter = 1500,  
      n.burn = 500,
      thin = 2,
      n_perm = 1000,
      n_cores_per_rep = N_CORES_PER_REP
    )
  }, error = function(e) {
    cat("\n", rep("!", 70), "\n", sep="")
    cat(sprintf("✗✗✗ REPETITION %d FAILED ✗✗✗\n", i))
    cat(rep("!", 70), "\n", sep="")
    cat("Error message:\n")
    print(e)
    cat("\n")
    return(structure(paste("Repetition", i, "failed:", e$message), 
                     class = "try-error"))
  })
}, mc.cores = N_REPS_PARALLEL)


overall_elapsed <- difftime(Sys.time(), overall_start, units = "mins")

cat("\n", rep("=", 70), "\n", sep="")
cat(sprintf("✓ ALL REPETITIONS COMPLETED in %.2f minutes\n", overall_elapsed)) # 1522.02 minutes 
cat(sprintf("Average time per repetition: %.2f minutes\n", overall_elapsed / N_REPS)) # Average time per repetition: 152.20 minutes
cat(rep("=", 70), "\n\n", sep="")

#══════════════════════════════════════════════════════════════════════════════
# AGGREGATE RESULTS
#══════════════════════════════════════════════════════════════════════════════

# Check for failures
failed_reps <- sapply(all_repetitions, function(x) {
  inherits(x, "try-error") || is.null(x)
})

n_successful <- sum(!failed_reps)
n_failed <- sum(failed_reps)
all_repetitions_new <- all_repetitions[!failed_reps]


cat("\n", rep("=", 70), "\n", sep="")
cat("AGGREGATED RESULTS ACROSS 10 REPETITIONS\n")
cat(rep("=", 70), "\n\n", sep="")

# Extract metrics
extract_metric <- function(all_reps, dataset, test_type, metric) {
  sapply(all_reps, function(rep) {
    if (test_type == "mcmc") {
      mean(c(rep[[dataset]]$mcmc_reliability$subset1[[metric]],
             rep[[dataset]]$mcmc_reliability$subset2[[metric]]))
    } else {
      mean(c(rep[[dataset]]$biological$seed1[[metric]],
             rep[[dataset]]$biological$seed10001[[metric]],
             rep[[dataset]]$biological$cross1[[metric]],
             rep[[dataset]]$biological$cross2[[metric]]))
    }
  })
}

# BA11 metrics
ba11_mcmc_adj <- extract_metric(all_repetitions_new, "ba11", "mcmc", "adjusted")
ba11_bio_adj <- extract_metric(all_repetitions_new, "ba11", "bio", "adjusted")

# BA47 metrics
ba47_mcmc_adj <- extract_metric(all_repetitions_new, "ba47", "mcmc", "adjusted")
ba47_bio_adj <- extract_metric(all_repetitions_new, "ba47", "bio", "adjusted")

# Print results
cat("--- BA11 ---\n")
cat(sprintf("MCMC Reliability:         %.4f ± %.4f\n",
            mean(ba11_mcmc_adj), sd(ba11_mcmc_adj)))
cat(sprintf("Biological Reproducibility: %.4f ± %.4f [%.4f, %.4f]\n\n",
            mean(ba11_bio_adj), sd(ba11_bio_adj),
            quantile(ba11_bio_adj, 0.025), quantile(ba11_bio_adj, 0.975)))

cat("--- BA47 ---\n")
cat(sprintf("MCMC Reliability:         %.4f ± %.4f\n",
            mean(ba47_mcmc_adj), sd(ba47_mcmc_adj)))
cat(sprintf("Biological Reproducibility: %.4f ± %.4f [%.4f, %.4f]\n\n",
            mean(ba47_bio_adj), sd(ba47_bio_adj),
            quantile(ba47_bio_adj, 0.025), quantile(ba47_bio_adj, 0.975)))

# Create summary table
final_summary <- data.frame(
  Dataset = rep(c("BA11", "BA47"), each = 2),
  Test_Type = rep(c("MCMC_Reliability", "Biological_Reproducibility"), 2),
  Mean = c(mean(ba11_mcmc_adj), mean(ba11_bio_adj),
           mean(ba47_mcmc_adj), mean(ba47_bio_adj)),
  SD = c(sd(ba11_mcmc_adj), sd(ba11_bio_adj),
         sd(ba47_mcmc_adj), sd(ba47_bio_adj)),
  CI_Lower = c(quantile(ba11_mcmc_adj, 0.025), quantile(ba11_bio_adj, 0.025),
               quantile(ba47_mcmc_adj, 0.025), quantile(ba47_bio_adj, 0.025)),
  CI_Upper = c(quantile(ba11_mcmc_adj, 0.975), quantile(ba11_bio_adj, 0.975),
               quantile(ba47_mcmc_adj, 0.975), quantile(ba47_bio_adj, 0.975))
)

print(final_summary, row.names = FALSE)

# Save results
output_dir <- file.path(current_aging, "results", "intra_cohort_10reps")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(all_repetitions, file.path(output_dir, "all_repetitions.rds"))
write.csv(final_summary, file.path(output_dir, "summary.csv"), row.names = FALSE)

cat("\n✓ Results saved to:", output_dir, "\n")

#══════════════════════════════════════════════════════════════════════════════
# INTERPRETATION GUIDE
#══════════════════════════════════════════════════════════════════════════════

cat("\n", rep("=", 70), "\n", sep="")
cat("INTERPRETATION\n")
cat(rep("=", 70), "\n", sep="")
cat("MCMC Reliability:\n")
cat("  >0.7  = Good convergence\n")
cat("  <0.5  = Need more iterations\n\n")
cat("Biological Reproducibility:\n")
cat("  >0.5  = Strong, reproducible rhythms\n")
cat("  0.3-0.5 = Moderate rhythms\n")
cat("  <0.3  = Weak rhythms or high variability\n")
cat(rep("=", 70), "\n\n", sep="")