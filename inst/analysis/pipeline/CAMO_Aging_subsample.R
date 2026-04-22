
####################################### only intra ###################
#── Clear environment ───────────────────────────────────────────────────────────────
rm(list=ls())

#── Paths ─────────────────────────────────────────────────────────────────────────
current_gtex <- "/home/qtp1/Projects/Collaborative"
current_wd   <- "/home/qtp1/Projects/Circadian"
current_aging <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging"


current_gtex <- "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative"
current_wd   <- "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian"
current_aging <- "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Collaborative/Paper/Congruence/PNAS_aging"



#── Load BA11 and BA47 data ──────────────────────────────────────────────────────
BA11 <- readRDS(file.path(current_aging, "data/BA11_data.rds"))
BA47 <- readRDS(file.path(current_aging, "data/BA47_data.rds"))
#local 
BA11<- readRDS("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/Paper/Congruence/PNAS_aging/data/BA11_data.rds")
BA47<- readRDS("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/Paper/Congruence/PNAS_aging/data/BA47_data.rds")
COMBINED <- readRDS("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/Paper/Congruence/PNAS_aging/data/combined_data.rds")

#--Process combined data ---
COMBINED <- readRDS(file.path(current_aging, "data/combined_data.rds"))

prepare_combined_data <- function(combined_data) {
  cat("Debugging combined data preparation...\n")
  
  expr_sample_names <- colnames(combined_data$expr)
  
  if(!"sample_name" %in% colnames(combined_data$pheno)) {
    stop("sample_name column not found in phenotype data")
  }
  
  cat("Aligning phenotype data to expression matrix order...\n")
  pheno_order <- match(expr_sample_names, combined_data$pheno$sample_name)
  
  if(any(is.na(pheno_order))) {
    missing_samples <- sum(is.na(pheno_order))
    cat("WARNING:", missing_samples, "expression samples not found in phenotype data\n")
    valid_expr_samples <- !is.na(pheno_order)
    expr_sample_names <- expr_sample_names[valid_expr_samples]
    combined_data$expr <- combined_data$expr[, valid_expr_samples]
    pheno_order <- pheno_order[valid_expr_samples]
  }
  
  pheno_data <- combined_data$pheno[pheno_order, ]
  
  # Handle TOD columns
  if("TOD.x" %in% colnames(pheno_data)) {
    pheno_data$tod <- pheno_data$TOD.x
  } else if("TOD.y" %in% colnames(pheno_data)) {
    pheno_data$tod <- pheno_data$TOD.y  
  } else if("TOD" %in% colnames(pheno_data)) {
    pheno_data$tod <- pheno_data$TOD
  } else {
    stop("Cannot find any TOD information in phenotype data")
  }
  
  # Handle age group
  if("AgeGroup" %in% colnames(pheno_data)) {
    pheno_data$age_group_final <- pheno_data$AgeGroup
  } else if("age_group" %in% colnames(pheno_data)) {
    pheno_data$age_group_final <- pheno_data$age_group
  } else {
    stop("Cannot find age group information in phenotype data")
  }
  
  complete_samples <- !is.na(pheno_data$age_group_final) & 
    !is.na(pheno_data$tod) &
    pheno_data$age_group_final %in% c("younger", "older")
  
  pheno_clean_final <- pheno_data[complete_samples, ]
  expr_clean <- combined_data$expr[, complete_samples]
  
  younger_samples <- pheno_clean_final$age_group_final == "younger"
  older_samples <- pheno_clean_final$age_group_final == "older"
  
  tod_all <- pheno_clean_final$tod
  tod_younger <- pheno_clean_final$tod[younger_samples]
  tod_older <- pheno_clean_final$tod[older_samples]
  
  result <- list(
    all_expr = expr_clean,
    younger_expr = expr_clean[, younger_samples],
    older_expr = expr_clean[, older_samples],
    all_tod = tod_all,
    tod_younger = tod_younger,
    tod_older = tod_older,
    sample_info = pheno_clean_final
  )
  
  return(result)
}

COMBINED_data <- prepare_combined_data(COMBINED)

#── Packages ─────────────────────────────────────────────────────────────────────
library(parallel); library(edgeR);   library(Rcpp)
require(DESeq2)
require(limma)
require('minpack.lm')
require(doParallel)
require(biomaRt)
require(devtools)
require(AWFisher)
require(ggplot2)
require(Rcpp)
require(dplyr)
require(pROC)
require(edgeR)
require(tidyr)

source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/pathwaySelect.R"))
source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/KEGG_module.R"))
source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/KEGG_module_topology_plot.R"))

source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/multi_ACS_ADS_pathway.R"))
source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/multi_ACS_ADS_global.R"))
source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/internal.R"))
sourceCpp(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/ACS.cpp"))

load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/kegg.pathway.list_hsa.RData"))
load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/kegg.pathway.list_cel_GeneNames.RData"))
load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/hw_orth.RData"))
load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/human.pathway.list.RData"))
load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/go.pathway.list_hsa.RData"))

source("/home/qtp1/Projects/Pipeline/one_cosinor_OLS_new.R")

#── Source R scripts ─────────────────────────────────────────────────────────────
WD <- "Kyle/Circadian-analysis-main/R/v1"
setwd(file.path(current_wd, WD))
scripts <- list.files("R", full.names=TRUE)
sapply(scripts, source)

#── MCMC runner function ─────────────────────────────────────────────────────────
run_MCMC_get_posteriors <- function(data, tod, n.iter=2500, n.burn=500, P=24, seed=1) {
  dat_list <- list(data = as.data.frame(data),
                   time = tod,
                   gname = rownames(data))
  
  a.init <- CBt_init_single(Data.list=dat_list, P=P,
                            FitCosinor=TRUE,
                            mu_M=0, sigma_M=10,
                            mu_A=1, sigma_A=10,
                            seed=seed)
  
  CB.res <- CB_MCMC_single_rj_slice(Data.list=dat_list,
                                    Init.value=a.init,
                                    P=P, iteration=n.iter,
                                    thin=1, n.burn=n.burn,
                                    seed=seed,
                                    p_rhythmic=rep(0.2,nrow(data)),
                                    rj.p.stay=0.5,
                                    A_prior="trunc_Normal_OLS_condi",
                                    mu_A=1, sigma_A=10^2, A.min=0,
                                    A_wb_beta2=2,
                                    A_gm_shape=1.99, A_gm_rate=0.5,
                                    rj.phi=TRUE, rj.A=TRUE,
                                    mu_M=0, sigma_M=10^2,
                                    sigma_prior_v=2, sigma_prior_s=0)
  # 1. If rho lost rownames, recover them
  if (is.null(rownames(CB.res$rho))) {
    rownames(CB.res$rho) <- attr(CB.res$rho, "dimnames")[[1]]
  }
  
  # 2. Reorder rho matrix to match input gene order
  CB.res$rho <- CB.res$rho[rownames(data), , drop = FALSE]
  posteriors <- rowMeans(CB.res$rho)
  names(posteriors) <- rownames(data)
  
  return(posteriors)
}

#── Adjusted Spearman correlation function ───────────────────────────────────
spearman <- function(p_A, p_B, n_perm = 1000) {
  obs <- cor(p_A, p_B, method = "spearman")
  
  # Permutation
  null_dist <- sapply(1:n_perm, function(i) {
    cor(p_A, sample(p_B), method = "spearman")
  })
  
  # Calculate null statistics FIRST
  null_mean <- mean(null_dist, na.rm = TRUE)
  null_sd <- sd(null_dist, na.rm = TRUE)
  
  # Then calculate adjusted values
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
    method = "Spearman",
    formula = "cor(rank(pA), rank(pB))",
    description = "Rank correlation (monotonic)",
    observed = obs,
    null_mean = null_mean,
    null_sd = null_sd,
    adjusted = adjusted,
    p_value = p_value,
    z_score = (obs - null_mean) / null_sd,
    null_dist = null_dist
  ))
}



#── Data preparation function ────────────────────────────────────────────────────
prepare_brain_data <- function(brain_data, region_name) {
  expr_ids <- sub(paste0(".*", region_name, "-([0-9]+)\\.CEL\\.gz"), "\\1", 
                  colnames(brain_data$expr))
  expr_ids <- as.integer(expr_ids)
  names(expr_ids) <- colnames(brain_data$expr)
  
  younger_ids <- brain_data$pheno$ID[brain_data$pheno$AgeGroup == "younger"]
  older_ids <- brain_data$pheno$ID[brain_data$pheno$AgeGroup == "older"]
  
  younger_expr <- brain_data$expr[, expr_ids %in% younger_ids]
  older_expr <- brain_data$expr[, expr_ids %in% older_ids]
  
  tod_map <- setNames(brain_data$pheno$TOD, brain_data$pheno$ID)
  
  ids_younger <- sub(paste0(".*", region_name, "-([0-9]+)\\.CEL\\.gz"), "\\1", 
                     colnames(younger_expr))
  tod_younger <- tod_map[ids_younger]
  
  ids_older <- sub(paste0(".*", region_name, "-([0-9]+)\\.CEL\\.gz"), "\\1", 
                   colnames(older_expr))
  tod_older <- tod_map[ids_older]
  
  list(
    all_expr = brain_data$expr,
    younger_expr = younger_expr,
    older_expr = older_expr,
    all_tod = tod_map[as.character(expr_ids)],
    tod_younger = tod_younger,
    tod_older = tod_older
  )
}

BA11_data <- prepare_brain_data(BA11, "BA11")
BA47_data <- prepare_brain_data(BA47, "BA47")

#── Random subsetting function ───────────────────────────────────────────────────
random_split_dataset <- function(expr_data, tod_vector, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n_samples <- ncol(expr_data)
  if (length(tod_vector) != n_samples) {
    stop("Expression data and TOD vector must have same length!")
  }
  
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

#══════════════════════════════════════════════════════════════════════════════
# INTRA-COHORT HETEROGENEITY ANALYSIS (PARALLELIZED REPETITIONS)
#══════════════════════════════════════════════════════════════════════════════

run_intra_cohort_analysis <- function(dataset, dataset_name, 
                                      n_repetitions = 100,
                                      base_seed = 1000,
                                      save_dir = NULL,
                                      n.iter = 2500,
                                      n.burn = 500,
                                      n_cores = 4,
                                      n_perm = 1000) {  # <--- ADD permutation parameter
  
  if (is.null(save_dir)) {
    save_dir <- file.path(current_aging, "results", "intra_cohort", dataset_name)
  }
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("\n===============================================\n")
  cat("INTRA-COHORT ANALYSIS:", dataset_name, "\n")
  cat("Repetitions:", n_repetitions, "\n")
  cat("Permutations per repetition:", n_perm, "\n")
  cat("Cores:", n_cores, "\n")
  cat("===============================================\n\n")
  
  # Define function for single repetition
  # Define function for single repetition
  process_single_rep <- function(i) {
    rep_seed <- base_seed + i
    
    # Randomly split the dataset
    split_data <- random_split_dataset(dataset$expr, dataset$tod, seed = rep_seed)
    
    tryCatch({
      #----------------------------------------------------------
      # Run MCMC on subset 1
      #----------------------------------------------------------
      p_subset1 <- run_MCMC_get_posteriors(
        data = split_data$subset1$expr,
        tod  = split_data$subset1$tod,
        n.iter = n.iter,
        n.burn = n.burn,
        seed = rep_seed
      )
      
      # Immediately clean memory
      gc(verbose = FALSE)
      
      #----------------------------------------------------------
      # Run MCMC on subset 2
      #----------------------------------------------------------
      p_subset2 <- run_MCMC_get_posteriors(
        data = split_data$subset2$expr,
        tod  = split_data$subset2$tod,
        n.iter = n.iter,
        n.burn = n.burn,
        seed = rep_seed + 10000
      )
      
      #----------------------------------------------------------
      # Compute concordance and free heavy objects
      #----------------------------------------------------------
      concordance_result <- spearman(p_subset1, p_subset2, n_perm = n_perm)
      rm(p_subset1, p_subset2, split_data)
      gc(verbose = FALSE)
      
      #----------------------------------------------------------
      # Save lightweight result to disk (optional but recommended)
      #----------------------------------------------------------
      out <- list(
        repetition = i,
        seed = rep_seed,
        observed = concordance_result$observed,
        adjusted = concordance_result$adjusted,
        p_value = concordance_result$p_value,
        null_mean = concordance_result$null_mean,
        null_sd = concordance_result$null_sd,
        z_score = concordance_result$z_score,
        status = "success"
      )
      
      saveRDS(out, file.path(save_dir, sprintf("rep_%03d.RDS", i)))  # optional
      
      return(out)
      
    }, error = function(e) {
      rm(list = ls()); gc()
      return(list(
        repetition = i,
        seed = rep_seed,
        observed = NA,
        adjusted = NA,
        p_value = NA,
        null_mean = NA,
        null_sd = NA,
        z_score = NA,
        status = "error",
        error_message = conditionMessage(e)
      ))
    })
  }
  
  
  # PARALLELIZE the repetitions
  cat("Starting parallel processing...\n")
  start_time <- Sys.time()
  
  concordance_results <- mclapply(1:n_repetitions, 
                                  process_single_rep, 
                                  mc.cores = n_cores)
  
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  
  cat(sprintf("\n✓ Parallel processing complete! Time: %.2f minutes\n", elapsed))
  
  # Convert to data frame
  results_df <- do.call(rbind, lapply(concordance_results, function(x) {
    data.frame(
      repetition = x$repetition,
      seed = x$seed,
      observed = x$observed,
      adjusted = x$adjusted,
      p_value = x$p_value,
      null_mean = x$null_mean,
      null_sd = x$null_sd,
      z_score = x$z_score,
      status = x$status,
      stringsAsFactors = FALSE
    )
  }))
  
  # Calculate summary
  successful <- results_df %>% filter(status == "success")
  
  summary_stats <- list(
    dataset = dataset_name,
    type = "intra_cohort",
    n_repetitions = nrow(successful),
    # Observed correlations
    observed_mean = mean(successful$observed, na.rm = TRUE),
    observed_sd = sd(successful$observed, na.rm = TRUE),
    observed_median = median(successful$observed, na.rm = TRUE),
    # Adjusted correlations
    adjusted_mean = mean(successful$adjusted, na.rm = TRUE),
    adjusted_sd = sd(successful$adjusted, na.rm = TRUE),
    adjusted_median = median(successful$adjusted, na.rm = TRUE),
    adjusted_min = min(successful$adjusted, na.rm = TRUE),
    adjusted_max = max(successful$adjusted, na.rm = TRUE),
    adjusted_ci_lower = quantile(successful$adjusted, 0.025, na.rm = TRUE),
    adjusted_ci_upper = quantile(successful$adjusted, 0.975, na.rm = TRUE),
    # P-values
    mean_p_value = mean(successful$p_value, na.rm = TRUE),
    # Null distribution stats
    mean_null_mean = mean(successful$null_mean, na.rm = TRUE),
    mean_null_sd = mean(successful$null_sd, na.rm = TRUE)
  )
  
  cat("\n--- INTRA-COHORT SUMMARY ---\n")
  cat("Observed Correlation:\n")
  cat(sprintf("  Mean: %.4f ± %.4f\n", summary_stats$observed_mean, summary_stats$observed_sd))
  cat("\nAdjusted Correlation:\n")
  cat(sprintf("  Mean: %.4f ± %.4f\n", summary_stats$adjusted_mean, summary_stats$adjusted_sd))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n", summary_stats$adjusted_ci_lower, summary_stats$adjusted_ci_upper))
  cat(sprintf("\nMean p-value: %.4f\n", summary_stats$mean_p_value))
  cat(sprintf("Time per repetition: %.2f minutes\n\n", as.numeric(elapsed) / n_repetitions))
  
  # Save results
  saveRDS(list(results = results_df, summary = summary_stats),
          file.path(save_dir, "analysis.RDS"))
  write.csv(results_df, file.path(save_dir, "results.csv"), row.names = FALSE)
  
  return(list(results = results_df, summary = summary_stats))
}

#══════════════════════════════════════════════════════════════════════════════
# DEFINE DATASETS
#══════════════════════════════════════════════════════════════════════════════

all_datasets <- list(
  BA11_all = list(
    expr = BA11_data$all_expr,
    tod = BA11_data$all_tod
  )
  ,BA47_all = list(
    expr = BA47_data$all_expr,
    tod = BA47_data$all_tod
  )
  , COMBINED_younger = list(
    expr = COMBINED_data$younger_expr,
    tod = COMBINED_data$tod_younger
  ),
  COMBINED_older = list(
    expr = COMBINED_data$older_expr,
    tod = COMBINED_data$tod_older
  )
)

#══════════════════════════════════════════════════════════════════════════════
# RUN INTRA-COHORT ANALYSIS (PARALLELIZED REPETITIONS)
#══════════════════════════════════════════════════════════════════════════════
# Set number of cores for repetitions
n_cores <- 10  

intra_results <- list()

# Process each dataset sequentially, but parallelize repetitions within each
for (dataset_name in names(all_datasets)) {
  cat("\n")
  cat("=" %>% rep(78) %>% paste0(collapse=""), "\n")
  cat("PROCESSING:", dataset_name, "\n")
  cat("=" %>% rep(78) %>% paste0(collapse=""), "\n")
  
  intra_results[[dataset_name]] <- run_intra_cohort_analysis(
    dataset = all_datasets[[dataset_name]],
    dataset_name = dataset_name,
    n_repetitions = 10,  
    base_seed = 1,
    n.iter = 2500,
    n.burn = 500,
    n_cores = n_cores  
  )
}

#══════════════════════════════════════════════════════════════════════════════
# CREATE SUMMARY
#══════════════════════════════════════════════════════════════════════════════

# Create summary table
intra_summary <- do.call(rbind, lapply(names(intra_results), function(name) {
  s <- intra_results[[name]]$summary
  data.frame(
    dataset = name,
    n_reps = s$n_repetitions,
    observed_mean = s$observed_mean,
    observed_sd = s$observed_sd,
    adjusted_mean = s$adjusted_mean,
    adjusted_sd = s$adjusted_sd,
    ci_lower = s$adjusted_ci_lower,
    ci_upper = s$adjusted_ci_upper,
    mean_p_value = s$mean_p_value,
    stringsAsFactors = FALSE
  )
}))

print(intra_summary, row.names = FALSE)

# Save summary
output_dir <- file.path(current_aging, "results", "intra_cohort_analysis")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(intra_summary,
          file.path(output_dir, "intra_cohort_summary.csv"),
          row.names = FALSE)

