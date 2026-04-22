################################################################################
# CORE FUNCTIONS FOR PROBABILISTIC concordance ANALYSIS
################################################################################
# Helper function: Compute probabilistic concordance index
# compute_probabilistic_concordance <- function(p_A, p_B) {
#   intersection_probs <- p_A * p_B
#   expected_intersection <- sum(intersection_probs)
#   expected_union <- sum(p_A) + sum(p_B) - expected_intersection
#   
#   if (expected_union == 0) return(0)
#   return(expected_intersection / expected_union)
# }
# Helper function: Compute normalized probabilistic concordance index (C_prod_star)
# compute_probabilistic_concordance <- function(p_A, p_B) {
#   # Inner helper: raw fuzzy-Jaccard using product form
#   C_prod <- function(x, y) {
#     num <- sum(x * y)
#     den <- sum(x) + sum(y) - num
#     if (den == 0) return(0)
#     num / den
#   }
#   
#   # Compute pairwise concordances
#   cAB <- C_prod(p_A, p_B)
#   cAA <- C_prod(p_A, p_A)
#   cBB <- C_prod(p_B, p_B)
#   
#   # Normalize so identical = 1
#   if (cAA == 0 || cBB == 0) {
#     return(NA_real_)
#   }
#   
#   normalized_concordance <- cAB / sqrt(cAA * cBB)
#   return(normalized_concordance)
# }
# 
# compute_probabilistic_concordance <- function(p_A, p_B) {
#   valid <- complete.cases(p_A, p_B)
#   if (sum(valid) < 3) return(NA_real_)
#   
#   # Spearman correlation as scalar
#   rho <- suppressWarnings(cor(p_A[valid], p_B[valid], method = "spearman"))
#   return(as.numeric(rho))   # <— critical fix
# }
# 

# NEW: Compute iteration-level Jaccard index
compute_iteration_jaccard <- function(rho_A, rho_B) {
  # rho_A, rho_B: binary matrices (genes × iterations)
  # Returns: list with jaccard values and confusion matrices
  
  n_iter <- ncol(rho_A)
  G <- nrow(rho_A)
  
  jaccard_obs <- numeric(n_iter)
  jaccard_exp <- numeric(n_iter)
  jaccard_adj <- numeric(n_iter)
  
  for (i in 1:n_iter) {
    A <- rho_A[, i]
    B <- rho_B[, i]
    
    # Confusion matrix elements
    a <- sum(A == 1 & B == 1)  # both rhythmic
    b <- sum(A == 1 & B == 0)  # loss (A only)
    c <- sum(A == 0 & B == 1)  # gain (B only)
    
    # Observed Jaccard
    union_size <- a + b + c
    if (union_size == 0) {
      jaccard_obs[i] <- NA_real_
      jaccard_exp[i] <- NA_real_
      jaccard_adj[i] <- NA_real_
    } else {
      jaccard_obs[i] <- a / union_size
      
      # Expected Jaccard under independence
      n_A <- sum(A)
      n_B <- sum(B)
      exp_intersection <- (n_A * n_B) / G
      exp_union <- n_A + n_B - exp_intersection
      
      if (exp_union > 0) {
        jaccard_exp[i] <- exp_intersection / exp_union
      } else {
        jaccard_exp[i] <- 0
      }
      
      # Null-adjusted Jaccard
      if (abs(jaccard_exp[i]) < 1) {
        jaccard_adj[i] <- (jaccard_obs[i] - jaccard_exp[i]) / (1 - jaccard_exp[i])
      } else {
        jaccard_adj[i] <- NA_real_
      }
    }
  }
  
  return(list(
    jaccard_obs = jaccard_obs,
    jaccard_exp = jaccard_exp,
    jaccard_adj = jaccard_adj,
    mean_jaccard_obs = mean(jaccard_obs, na.rm = TRUE),
    mean_jaccard_adj = mean(jaccard_adj, na.rm = TRUE)
  ))
}

# Helper function: Compute gain and loss indices
# compute_gain_loss <- function(p_A, p_B) {
#   intersection_probs <- p_A * p_B
#   expected_intersection <- sum(intersection_probs)
#   expected_union <- sum(p_A) + sum(p_B) - expected_intersection
#   
#   if (expected_union == 0) {
#     return(list(
#       gain_index = 0,
#       loss_index = 0,
#       gain_loss_ratio = NA_real_
#     ))
#   }
#   
#   expected_loss <- sum(p_A) - expected_intersection
#   expected_gain <- sum(p_B) - expected_intersection
#   
#   union_inv <- 1 / expected_union
#   loss_index <- expected_loss * union_inv
#   gain_index <- expected_gain * union_inv
#   
#   gain_loss_ratio <- if (loss_index == 0) {
#     if (gain_index == 0) {
#       NA_real_
#     } else {
#       Inf
#     }
#   } else {
#     gain_index / loss_index
#   }
#   
#   return(list(
#     gain_index = gain_index,
#     loss_index = loss_index,
#     gain_loss_ratio = gain_loss_ratio
#   ))
# }

# Compute gain/loss across iterations
compute_gain_loss_iterations <- function(rho_A, rho_B) {
  # rho_A, rho_B: binary matrices (genes × iterations)
  
  n_iter <- ncol(rho_A)
  
  exp_conserved <- numeric(n_iter)
  exp_gained <- numeric(n_iter)
  exp_lost <- numeric(n_iter)
  
  for (i in 1:n_iter) {
    A <- rho_A[, i]
    B <- rho_B[, i]
    
    exp_conserved[i] <- sum(A == 1 & B == 1)
    exp_lost[i] <- sum(A == 1 & B == 0)
    exp_gained[i] <- sum(A == 0 & B == 1)
  }
  
  # Average across iterations
  mean_conserved <- mean(exp_conserved)
  mean_gained <- mean(exp_gained)
  mean_lost <- mean(exp_lost)
  
  # Calculate indices (matching your original logic)
  total_union <- mean_conserved + mean_gained + mean_lost
  
  if (total_union == 0) {
    return(list(
      expected_conserved = mean_conserved,
      expected_gained = mean_gained,
      expected_lost = mean_lost,
      gain_index = 0,
      loss_index = 0,
      gain_loss_ratio = NA_real_
    ))
  }
  
  gain_index <- mean_gained / total_union
  loss_index <- mean_lost / total_union
  
  gain_loss_ratio <- if (loss_index == 0) {
    if (gain_index == 0) NA_real_ else Inf
  } else {
    gain_index / loss_index
  }
  
  return(list(
    expected_conserved = mean_conserved,
    expected_gained = mean_gained,
    expected_lost = mean_lost,
    gain_index = gain_index,
    loss_index = loss_index,
    gain_loss_ratio = gain_loss_ratio
  ))
}


################################################################################
# MAIN ANALYSIS FUNCTION FOR SINGLE PATHWAY
################################################################################
analyze_pathway_concordance <- function(rho_A, rho_B, idx = NULL, 
                                        n_perm = 1000, n_boot = 500,
                                        alpha = 0.05, use_cpp = FALSE, is_global = FALSE,
                                        compute_pvalue = TRUE,
                                        compute_ci = TRUE) {
  
  # Subset if pathway indices provided
  if (!is.null(idx)) {
    rho_A_subset <- rho_A[idx, , drop = FALSE]
    rho_B_subset <- rho_B[idx, , drop = FALSE]
  } else {
    rho_A_subset <- rho_A
    rho_B_subset <- rho_B
  }
  
  n_genes <- nrow(rho_A_subset)
  n_iter <- ncol(rho_A_subset)
  
  # 1. Calculate observed concordance (iteration-level Jaccard)
  if (use_cpp) {
    obs_result <- compute_iteration_jaccard_cpp(rho_A_subset, rho_B_subset)
  } else {
    obs_result <- compute_iteration_jaccard(rho_A_subset, rho_B_subset)
  }
  obs_concordance <- obs_result$mean_jaccard_obs
  adjusted_concordance <- obs_result$mean_jaccard_adj
  
  # 2. Calculate gain/loss indices
  if (use_cpp) {
    gain_loss <- compute_gain_loss_iterations_cpp(rho_A_subset, rho_B_subset)
  } else {
    gain_loss <- compute_gain_loss_iterations(rho_A_subset, rho_B_subset)
  }
  
  # 3. PERMUTATION TEST (optional)
  if (compute_pvalue) {
    if (use_cpp) {
      null_adj_scores <- permutation_jaccard_cpp(rho_A_subset, rho_B_subset, 
                                                 n_perm = n_perm, seed = 12345)
      
      # Gain/loss permutation (only for global)
      if (is_global) {
        gainloss_p_value <- gainloss_pvalue_iterations_cpp(rho_A_subset, rho_B_subset,
                                                           n_perm = n_perm, seed = 12345)
      } else {
        gainloss_p_value <- NA_real_
      }
      
    } else {
      null_adj_scores <- numeric(n_perm)
      gain_loss_perm <- if (is_global) numeric(n_perm) else rep(NA_real_, n_perm)
      
      for (b in 1:n_perm) {
        perm_idx <- sample(n_genes)
        rho_B_perm <- rho_B_subset[perm_idx, , drop = FALSE]
        
        perm_result <- compute_iteration_jaccard(rho_A_subset, rho_B_perm)
        null_adj_scores[b] <- perm_result$mean_jaccard_adj
        
        if (is_global) {
          gl_perm <- compute_gain_loss_iterations(rho_A_subset, rho_B_perm)
          gain_loss_perm[b] <- gl_perm$gain_loss_ratio
        }
      }
      
      # Calculate gain/loss p-value (only for R version)
      if (is_global) {
        obs_ratio <- gain_loss$gain_loss_ratio
        obs_dev <- abs(log(obs_ratio))
        perm_dev <- abs(log(gain_loss_perm))
        gainloss_p_value <- (1 + sum(perm_dev >= obs_dev, na.rm = TRUE)) / (1 + n_perm)
      } else {
        gainloss_p_value <- NA_real_
      }
    }
    
    # Calculate null statistics
    null_mean <- mean(null_adj_scores, na.rm = TRUE)
    null_sd <- sd(null_adj_scores, na.rm = TRUE)
    
    # Calculate empirical p-value
    if (!is.na(adjusted_concordance) && !all(is.na(null_adj_scores))) {
      p_value <- (1 + sum(null_adj_scores >= adjusted_concordance, na.rm = TRUE)) / (1 + n_perm)
    } else {
      p_value <- NA_real_
    }
    
    # Z-score
    z_score <- if (null_sd > 0) {
      (adjusted_concordance - null_mean) / null_sd
    } else {
      NA_real_
    }
    
  } else {
    # Skip permutation test
    null_adj_scores <- numeric(0)
    p_value <- NA_real_
    gainloss_p_value <- NA_real_
    z_score <- NA_real_
    null_mean <- NA_real_
    null_sd <- NA_real_
  }
  
  # 4. BOOTSTRAP (optional)
  if (compute_ci) {
    if (use_cpp) {
      boot_result <- bootstrap_jaccard_cpp(rho_A_subset, rho_B_subset, 
                                           n_boot = n_boot, seed = 12345)
      boot_raw <- boot_result$boot_obs
      boot_adj <- boot_result$boot_adj
      
    } else {
      boot_raw <- numeric(n_boot)
      boot_adj <- numeric(n_boot)
      
      for (b in 1:n_boot) {
        boot_idx <- sample(n_genes, n_genes, replace = TRUE)
        rho_A_boot <- rho_A_subset[boot_idx, , drop = FALSE]
        rho_B_boot <- rho_B_subset[boot_idx, , drop = FALSE]
        
        boot_result <- compute_iteration_jaccard(rho_A_boot, rho_B_boot)
        boot_raw[b] <- boot_result$mean_jaccard_obs
        boot_adj[b] <- boot_result$mean_jaccard_adj
      }
    }
    
    # Remove NAs
    boot_raw <- boot_raw[!is.na(boot_raw)]
    boot_adj <- boot_adj[!is.na(boot_adj)]
    
    # Confidence intervals
    ci_raw <- quantile(boot_raw, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
    ci_adj <- quantile(boot_adj, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
    
  } else {
    # Skip bootstrap
    boot_raw <- numeric(0)
    boot_adj <- numeric(0)
    ci_raw <- c(NA_real_, NA_real_)
    ci_adj <- c(NA_real_, NA_real_)
  }
  
  return(list(
    # Observed values
    concordance = obs_concordance,
    adjusted_concordance = adjusted_concordance,
    expected_conserved = gain_loss$expected_conserved,
    expected_gained = gain_loss$expected_gained,
    expected_lost = gain_loss$expected_lost,
    gain_index = gain_loss$gain_index,
    loss_index = gain_loss$loss_index,
    gain_loss_ratio = gain_loss$gain_loss_ratio,
    gain_loss_p_value = gainloss_p_value,
    
    # Statistical inference
    p_value = p_value,
    z_score = z_score,
    null_mean = null_mean,
    null_sd = null_sd,
    
    # Confidence intervals
    ci_raw = ci_raw,
    ci_adj = ci_adj,
    
    # Additional info
    n_genes = n_genes,
    n_iter = n_iter,
    n_perm = if (compute_pvalue) n_perm else 0,
    n_boot = length(boot_raw),
    
    # Full distributions (optional, for diagnostics)
    null_distribution = null_adj_scores,
    boot_distribution_raw = boot_raw,
    boot_distribution_adj = boot_adj
  ))
}

################################################################################
# WRAPPER FOR MULTIPLE PATHWAYS
################################################################################

bootstrap_conservation_pathway_within <- function(dat1, dat2, pathway.list,
                                                  n_perm = 1000, n_boot = 500,
                                                  alpha = 0.05, use_cpp = FALSE,
                                                  compute_pvalue = TRUE,
                                                  compute_ci = TRUE) {
  
  # Extract FULL binary matrices (not averaged!)
  rho_A <- dat1$rho  # genes × iterations
  rho_B <- dat2$rho  # genes × iterations
  genes <- attr(dat1$rho, "symbols")
  
  # Determine if global or pathway-specific
  is_global <- is.character(pathway.list) && length(pathway.list) == 1 && 
    tolower(pathway.list) == "global"
  
  if (is_global) {
    K <- 1
    pathway_names <- "Global"
  } else {
    K <- length(pathway.list)
    pathway_names <- names(pathway.list)
  }
  
  # Initialize storage
  results <- list(
    conservation_scores = rep(NA_real_, K),
    adjusted_scores = rep(NA_real_, K),
    ci_lower_raw = rep(NA_real_, K),
    ci_upper_raw = rep(NA_real_, K),
    ci_lower_adj = rep(NA_real_, K),
    ci_upper_adj = rep(NA_real_, K),
    p_values = rep(NA_real_, K),
    z_scores = rep(NA_real_, K),
    null_means = rep(NA_real_, K),
    null_sds = rep(NA_real_, K),
    gain_indices = rep(NA_real_, K),
    loss_indices = rep(NA_real_, K),
    gain_loss_ratios = rep(NA_real_, K),
    gain_loss_pvalues = rep(NA_real_, K),
    pathway_sizes = rep(NA_real_, K)
  )
  
  # Store full distributions if needed
  boot_samples_matrix <- matrix(NA, nrow = K, ncol = if (compute_ci) n_boot else 1)
  null_samples_matrix <- matrix(NA, nrow = K, ncol = if (compute_pvalue) n_perm else 1)
  
  # Analyze each pathway
  for (k in 1:K) {
    cat(sprintf("Analyzing pathway %d/%d: %s\n", k, K, pathway_names[k]))
    
    # Get pathway indices
    if (is_global) {
      idx <- seq_along(genes)
    } else {
      idx <- which(genes %in% pathway.list[[k]])
    }
    
    if (length(idx) == 0) {
      warning(sprintf("Pathway %s has no genes in dataset", pathway_names[k]))
      next
    }
    
    results$pathway_sizes[k] <- length(idx)
    
    # Run analysis
    pathway_result <- analyze_pathway_concordance(
      rho_A, rho_B, idx = idx,  
      n_perm = n_perm, n_boot = n_boot,
      alpha = alpha, use_cpp = use_cpp, is_global = is_global,
      compute_pvalue = compute_pvalue,
      compute_ci = compute_ci
    )
    
    # Store results
    results$conservation_scores[k] <- pathway_result$concordance
    results$adjusted_scores[k] <- pathway_result$adjusted_concordance
    results$ci_lower_raw[k] <- pathway_result$ci_raw[1]
    results$ci_upper_raw[k] <- pathway_result$ci_raw[2]
    results$ci_lower_adj[k] <- pathway_result$ci_adj[1]
    results$ci_upper_adj[k] <- pathway_result$ci_adj[2]
    results$p_values[k] <- pathway_result$p_value
    results$z_scores[k] <- pathway_result$z_score
    results$null_means[k] <- pathway_result$null_mean
    results$null_sds[k] <- pathway_result$null_sd
    results$gain_indices[k] <- pathway_result$gain_index
    results$loss_indices[k] <- pathway_result$loss_index
    results$gain_loss_ratios[k] <- pathway_result$gain_loss_ratio
    results$gain_loss_pvalues[k] <- pathway_result$gain_loss_p_value
    
    # Store distributions
    if (compute_ci && length(pathway_result$boot_distribution_raw) > 0) {
      boot_samples_matrix[k, 1:length(pathway_result$boot_distribution_raw)] <- 
        pathway_result$boot_distribution_raw
    }
    if (compute_pvalue && length(pathway_result$null_distribution) > 0) {
      null_samples_matrix[k, 1:length(pathway_result$null_distribution)] <- 
        pathway_result$null_distribution
    }
  }
  
  # Convert to named matrices
  create_named_matrix <- function(data, rows, col_name) {
    matrix(data, ncol = 1, dimnames = list(rows, col_name))
  }
  
  return(list(
    conservation_scores = create_named_matrix(results$conservation_scores, pathway_names, "concordance"),
    adjusted_scores = create_named_matrix(results$adjusted_scores, pathway_names, "adjusted_concordance"),
    ci_lower_raw = create_named_matrix(results$ci_lower_raw, pathway_names, "ci_lower_raw"),
    ci_upper_raw = create_named_matrix(results$ci_upper_raw, pathway_names, "ci_upper_raw"),
    ci_lower_adj = create_named_matrix(results$ci_lower_adj, pathway_names, "ci_lower_adj"),
    ci_upper_adj = create_named_matrix(results$ci_upper_adj, pathway_names, "ci_upper_adj"),
    p_values = create_named_matrix(results$p_values, pathway_names, "p_value"),
    z_scores = create_named_matrix(results$z_scores, pathway_names, "z_score"),
    null_means = create_named_matrix(results$null_means, pathway_names, "null_mean"),
    null_sds = create_named_matrix(results$null_sds, pathway_names, "null_sd"),
    gain_indices = create_named_matrix(results$gain_indices, pathway_names, "gain_index"),
    loss_indices = create_named_matrix(results$loss_indices, pathway_names, "loss_index"),
    gain_loss_ratios = create_named_matrix(results$gain_loss_ratios, pathway_names, "gain_loss_ratio"),
    gain_loss_pvalues = create_named_matrix(results$gain_loss_pvalues,
                                            pathway_names,
                                            "gain_loss_p_value"),
    pathway_sizes = create_named_matrix(results$pathway_sizes, pathway_names, "n_genes"),
    boot_samples = boot_samples_matrix,
    null_samples = null_samples_matrix
  ))
}





################################################################################
# MULTI-DATASET COMPARISON
################################################################################

multi_conservation <- function(mcmc.merge.list, dataset.names,
                               select.pathway.list,
                               n_perm = 1000, n_boot = 500,
                               alpha = 0.05, min.p = 5,
                               qvalue_threshold = NULL,
                               output.dir = "Conservation_Pathway_Bootstrap",
                               use_cpp = FALSE,
                               compute_pvalue = TRUE,
                               compute_ci = TRUE) {
  
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required")
  }
  
  is_global <- is.character(select.pathway.list) && 
    length(select.pathway.list) == 1 && 
    tolower(select.pathway.list) == "global"
  
  names(mcmc.merge.list) <- dataset.names
  M <- length(mcmc.merge.list)
  P <- choose(M, 2)
  data_genes <- attr(mcmc.merge.list[[1]]$rho, "symbols")
  
  # Filter pathways by size
  if (!is_global) {
    pathway_valid <- sapply(select.pathway.list, function(p) {
      length(intersect(data_genes, p)) >= min.p
    })
    select.pathway.list <- select.pathway.list[pathway_valid]
    K <- length(select.pathway.list)
    
    if (K == 0) {
      stop("No pathways remain after filtering by min.p size.")
    }
  } else {
    K <- 1
  }
  
  pathway_names <- if (is_global) "Global" else names(select.pathway.list)
  combinations <- combn(dataset.names, m = 2)
  pair_names <- apply(combinations, 2, function(x) paste(x, collapse = "&"))
  
  # Initialize result matrices
  matrix_names <- c("Conservation", "Adjusted", "CI_lower_raw", "CI_upper_raw",
                    "CI_lower_adj", "CI_upper_adj", "PValue", "ZScore",
                    "Null_mean", "Null_sd", "Gain_index", "Loss_index", 
                    "Gain_loss_ratio", "Gain_loss_p_value", "Pathway_size")
  
  results <- lapply(setNames(matrix_names, matrix_names), function(x) {
    matrix(NA, nrow = K, ncol = P, dimnames = list(pathway_names, pair_names))
  })
  
  # Process each pair
  for (p_idx in 1:P) {
    i_name <- combinations[1, p_idx]
    j_name <- combinations[2, p_idx]
    
    cat("\n=== Processing pair", p_idx, "of", P, ":", i_name, "vs", j_name, "===\n")
    
    dat1 <- mcmc.merge.list[[i_name]]
    dat2 <- mcmc.merge.list[[j_name]]
    
    # Run analysis
    pair_results <- bootstrap_conservation_pathway_within(
      dat1, dat2, select.pathway.list,
      n_perm = n_perm, n_boot = n_boot,
      alpha = alpha, use_cpp = use_cpp,
      compute_pvalue = compute_pvalue,
      compute_ci = compute_ci
    )
    
    # Store results
    results$Conservation[, p_idx] <- pair_results$conservation_scores
    results$Adjusted[, p_idx] <- pair_results$adjusted_scores
    results$CI_lower_raw[, p_idx] <- pair_results$ci_lower_raw
    results$CI_upper_raw[, p_idx] <- pair_results$ci_upper_raw
    results$CI_lower_adj[, p_idx] <- pair_results$ci_lower_adj
    results$CI_upper_adj[, p_idx] <- pair_results$ci_upper_adj
    results$PValue[, p_idx] <- pair_results$p_values
    results$ZScore[, p_idx] <- pair_results$z_scores
    results$Null_mean[, p_idx] <- pair_results$null_means
    results$Null_sd[, p_idx] <- pair_results$null_sds
    results$Gain_index[, p_idx] <- pair_results$gain_indices
    results$Loss_index[, p_idx] <- pair_results$loss_indices
    results$Gain_loss_ratio[, p_idx] <- pair_results$gain_loss_ratios
    results$Gain_loss_p_value[, p_idx] <- pair_results$gain_loss_pvalues
    results$Pathway_size[, p_idx] <- pair_results$pathway_sizes
  }
  
  # FDR correction
  if (!is_global && compute_pvalue) {
    all_pvalues <- as.vector(results$PValue)
    valid_mask <- !is.na(all_pvalues)
    
    qvalues <- matrix(NA, nrow = K, ncol = P, dimnames = list(pathway_names, pair_names))
    
    if (sum(valid_mask) > 0) {
      qvalues[valid_mask] <- p.adjust(all_pvalues[valid_mask], method = "fdr")
    }
    
    results$QValue <- qvalues
  }
  
  # Create output dataframe
  final_df_list <- list(Pathway = pathway_names)
  
  for (p_idx in 1:P) {
    pair_name <- gsub("&", "_vs_", pair_names[p_idx])
    
    final_df_list[[paste0(pair_name, "_Concordance")]] <- results$Conservation[, p_idx]
    final_df_list[[paste0(pair_name, "_AdjustedConcordance")]] <- results$Adjusted[, p_idx]
    
    if (compute_ci) {
      final_df_list[[paste0(pair_name, "_CI_Lower_Raw")]] <- results$CI_lower_raw[, p_idx]
      final_df_list[[paste0(pair_name, "_CI_Upper_Raw")]] <- results$CI_upper_raw[, p_idx]
      final_df_list[[paste0(pair_name, "_CI_Lower_Adj")]] <- results$CI_lower_adj[, p_idx]
      final_df_list[[paste0(pair_name, "_CI_Upper_Adj")]] <- results$CI_upper_adj[, p_idx]
    }
    
    if (compute_pvalue) {
      final_df_list[[paste0(pair_name, "_PValue")]] <- results$PValue[, p_idx]
      if (!is_global) {
        final_df_list[[paste0(pair_name, "_QValue")]] <- results$QValue[, p_idx]
      }
      final_df_list[[paste0(pair_name, "_ZScore")]] <- results$ZScore[, p_idx]
    }
    
    final_df_list[[paste0(pair_name, "_GainIndex")]] <- results$Gain_index[, p_idx]
    final_df_list[[paste0(pair_name, "_LossIndex")]] <- results$Loss_index[, p_idx]
    final_df_list[[paste0(pair_name, "_GainLossRatio")]] <- results$Gain_loss_ratio[, p_idx]
    
    if (compute_pvalue && is_global) {
      final_df_list[[paste0(pair_name, "_GainLossPValue")]] <- results$Gain_loss_p_value[, p_idx]
    }
    
    final_df_list[[paste0(pair_name, "_PathwaySize")]] <- results$Pathway_size[, p_idx]
  }
  
  final_df <- as.data.frame(final_df_list)
  
  # Save to Excel
  if (!dir.exists(output.dir)) {
    dir.create(output.dir, recursive = TRUE)
  }
  
  wb <- openxlsx::createWorkbook()
  
  # Add results sheet
  openxlsx::addWorksheet(wb, sheetName = "Results")
  openxlsx::writeData(wb, sheet = "Results", x = final_df)
  
  header_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#E6E6FA")
  openxlsx::addStyle(wb, sheet = "Results", style = header_style, 
                     rows = 1, cols = 1:ncol(final_df), gridExpand = TRUE)
  openxlsx::setColWidths(wb, sheet = "Results", cols = 1:ncol(final_df), widths = "auto")
  
  # Add documentation sheet
  doc_df <- data.frame(
    Column = c("Pathway", "Concordance", "AdjustedConcordance", 
               "CI_Lower_Raw", "CI_Upper_Raw", "CI_Lower_Adj", "CI_Upper_Adj",
               "PValue", "QValue", "ZScore", "GainIndex", "LossIndex", 
               "GainLossRatio", "GainLossPValue", "PathwaySize"),
    Description = c(
      "Pathway or gene set name",
      "Raw concordance index: intersection/union of posterior probabilities",
      "Adjusted concordance: (observed - expected)/(1 - expected)",
      "Lower bound of 95% bootstrap CI for raw concordance",
      "Upper bound of 95% bootstrap CI for raw concordance",
      "Lower bound of 95% bootstrap CI for adjusted concordance",
      "Upper bound of 95% bootstrap CI for adjusted concordance",
      "Empirical one-sided p-value from permutation test",
      "FDR-corrected q-value (Benjamini-Hochberg)",
      "Standardized effect size: (observed - expected)/SD(null)",
      "Probabilistic gain index: genes rhythmic in B but not A",
      "Probabilistic loss index: genes rhythmic in A but not B",
      "Ratio of gain to loss indices",
      "P-value for gain/loss ratio deviation from 1 (global only)",
      "Number of genes in pathway"
    ),
    stringsAsFactors = FALSE
  )
  
  openxlsx::addWorksheet(wb, sheetName = "Column_Definitions")
  openxlsx::writeData(wb, sheet = "Column_Definitions", x = doc_df)
  openxlsx::addStyle(wb, sheet = "Column_Definitions", style = header_style,
                     rows = 1, cols = 1:ncol(doc_df), gridExpand = TRUE)
  openxlsx::setColWidths(wb, sheet = "Column_Definitions", 
                         cols = 1:ncol(doc_df), widths = "auto")
  
  excel_filename <- file.path(output.dir, "Conservation_Results.xlsx")
  openxlsx::saveWorkbook(wb, excel_filename, overwrite = TRUE)
  
  cat("\nResults saved to:", excel_filename, "\n")
  
  # Filter results if threshold provided
  if (!is.null(qvalue_threshold) && !is_global && compute_pvalue) {
    qvalue_cols <- grep("_QValue$", colnames(final_df))
    if (length(qvalue_cols) > 0) {
      significant_rows_mask <- apply(final_df[, qvalue_cols, drop = FALSE], 1, 
                                     function(row) any(row < qvalue_threshold, na.rm = TRUE))
      filtered_df <- final_df[significant_rows_mask, ]
      cat("\nReturning", nrow(filtered_df), "of", nrow(final_df),
          "pathways with at least one QValue <", qvalue_threshold, "\n")
      return(filtered_df)
    }
  }
  
  return(final_df)
}