
#' Permutation null distribution for global congruence index
#'
#' @title Global permutation test for circadian congruence
#'
#' @description
#' Generates a permutation null distribution for the congruence index and
#' related metrics (gain, loss, gain/loss ratio) by randomly shuffling the
#' row order of \code{dat2$rho} while keeping \code{dat1$rho} fixed.
#' Supports parallelisation via \code{parallel::mclapply} (Unix only) and
#' multi-round execution to reduce peak memory use.
#'
#' @param dat1 Named list; MCMC output for condition 1 with a \code{rho}
#'   matrix (G x K binary).
#' @param dat2 Named list; MCMC output for condition 2 with a \code{rho}
#'   matrix.
#' @param delta Numeric; passed to \code{congruence} (default 3; currently
#'   unused in the probabilistic score).
#' @param units Character; \code{"hours"} or \code{"radians"} (default
#'   \code{"hours"}).
#' @param B Integer; total number of permutations (default 1000).
#' @param ncores Integer or \code{NULL}; number of cores for parallel runs;
#'   \code{NULL} auto-detects (up to 4).
#' @param parallel Logical or \code{"auto"}; whether to parallelise
#'   (default \code{"auto"}: enabled on Unix when B >= 100).
#' @param rounds Integer; split B permutations into this many rounds to
#'   reduce memory use (default 1).
#' @param save_intermediate Logical; save each round to \code{intermediate_dir}
#'   (default \code{FALSE}).
#' @param intermediate_dir Character; directory for intermediate saves.
#'
#' @return A B x 4 numeric matrix with columns \code{congruence_index},
#'   \code{gain_index}, \code{loss_index}, \code{gain_loss_ratio}.
#'
#' @export
# Global permutation function with rounds support
perm_conservation_global <- function(dat1, dat2, delta = 3, units = "hours", B = 1000,
                                     ncores = NULL, parallel = "auto",
                                     rounds = 1, save_intermediate = FALSE,
                                     intermediate_dir = NULL) {
  
  # Auto-detect best approach - convert to logical
  if (parallel == "auto") {
    parallel <- .Platform$OS.type != "windows" && B >= 100
  }
  
  # Ensure parallel is logical
  parallel <- as.logical(parallel)
  
  # Set default ncores if using parallel
  if (is.null(ncores)) {
    ncores <- min(parallel::detectCores(logical = FALSE), 4)
  }
  
  # Calculate permutations per round
  B_per_round <- ceiling(B / rounds)
  total_B <- B_per_round * rounds
  
  if (rounds > 1) {
    cat("Splitting", B, "permutations into", rounds, "rounds of", B_per_round, "each\n")
    cat("Total permutations:", total_B, "\n")
  }
  
  # Setup intermediate saving if requested
  if (save_intermediate && !is.null(intermediate_dir)) {
    if (!dir.exists(intermediate_dir)) {
      dir.create(intermediate_dir, recursive = TRUE)
    }
  }
  
  G <- nrow(dat1$rho)
  # Updated metrics (removed Cp and Dp, added gain_loss_ratio)
  metrics <- c("congruence_index", "gain_index", "loss_index", "gain_loss_ratio")
  
  # Create output array for ALL metrics
  out <- array(NA, dim = c(total_B, length(metrics)), 
               dimnames = list(NULL, metrics))
  
  # Process each round
  for (round_num in 1:rounds) {
    if (rounds > 1) {
      cat("\n--- Round", round_num, "of", rounds, "---\n")
    }
    
    # Calculate start and end indices for this round
    start_idx <- (round_num - 1) * B_per_round + 1
    end_idx <- min(round_num * B_per_round, total_B)
    current_B <- end_idx - start_idx + 1
    
    if (rounds > 1) {
      cat("Processing permutations", start_idx, "to", end_idx, "(", current_B, "permutations)\n")
    }
    
    # Run permutations for this round
    round_results <- run_single_round_global(
      dat1, dat2, G, metrics, delta, units, current_B, parallel, ncores
    )
    
    # Store results
    out[start_idx:end_idx, ] <- round_results
    
    # Save intermediate results if requested
    if (save_intermediate && !is.null(intermediate_dir)) {
      round_file <- file.path(intermediate_dir, paste0("round_", round_num, "_results.rds"))
      saveRDS(round_results, round_file)
      if (rounds > 1) {
        cat("Intermediate results saved to:", round_file, "\n")
      }
    }
    
    # Memory cleanup
    gc()
    
    if (rounds > 1) {
      cat("Round", round_num, "completed\n")
    }
  }
  
  # If we generated more permutations than requested, trim to exact B
  if (total_B > B) {
    cat("Trimming results to exactly", B, "permutations\n")
    out <- out[1:B, , drop = FALSE]
  }
  
  return(out)
}

#' Run one round of global permutations
#'
#' @description
#' Executes \code{B_round} permutations of the congruence index and
#' related metrics for a single round.  Called by
#' \code{perm_conservation_global}.
#'
#' @param dat1,dat2 MCMC output lists with \code{rho} matrices.
#' @param G Integer; number of genes.
#' @param metrics Character vector; metric names to compute.
#' @param delta,units Passed to \code{congruence}.
#' @param B_round Integer; permutations in this round.
#' @param parallel Logical; whether to use \code{mclapply}.
#' @param ncores Integer; number of cores.
#'
#' @return B_round x length(metrics) numeric matrix.
#'
#' @keywords internal
# Helper function to run a single round of global permutations
run_single_round_global <- function(dat1, dat2, G, metrics, delta, units, B_round,
                                    parallel, ncores) {
  
  # Create output array for this round
  round_out <- array(NA, dim = c(B_round, length(metrics)), 
                     dimnames = list(NULL, metrics))
  
  if (parallel && .Platform$OS.type != "windows") {
    
    # PARALLEL VERSION (Mac/Linux)
    library(parallel)
    
    perm_results <- mclapply(1:B_round, function(b) {
      perm_indices <- sample(1:G, G, replace = FALSE)
      
      # Only rho needed for congruence calculation
      dat1_perm <- list(
        rho = dat1$rho
      )
      dat2_perm <- list(
        rho = dat2$rho[perm_indices, , drop = FALSE]
      )
      
      perm_result <- congruence(dat1_perm, dat2_perm, delta = delta, units = units)
      
      result_vector <- numeric(length(metrics))
      names(result_vector) <- metrics
      for(metric in metrics) {
        result_vector[metric] <- perm_result[[metric]]
      }
      return(result_vector)
    }, mc.cores = ncores)
    
    # Combine parallel results
    for (b in 1:B_round) {
      round_out[b, ] <- perm_results[[b]]
    }
    
  } else {
    
    # SEQUENTIAL VERSION
    cat("  Using sequential processing for this round\n")
    
    # Setup progress tracking for this round
    progress_interval <- max(1, floor(B_round/10))
    
    for (b in 1:B_round) {
      # Show progress at regular intervals
      if (b %% progress_interval == 0 || b == B_round) {
        percentage <- round(100 * b / B_round, 1)
        cat("    Round progress: Permutation", b, "of", B_round, "(", percentage, "%)\n")
        flush.console()
      }
      
      perm_indices <- sample(1:G, G, replace = FALSE)
      
      # Only rho needed for congruence calculation
      dat1_perm <- list(
        rho = dat1$rho
      )
      dat2_perm <- list(
        rho = dat2$rho[perm_indices, , drop = FALSE]
      )
      
      perm_result <- congruence(dat1_perm, dat2_perm, delta = delta, units = units)
      
      for(metric in metrics) {
        round_out[b, metric] <- perm_result[[metric]]
      }
    }
  }
  
  return(round_out)
}

#' Compute permutation p-values for global congruence scores
#'
#' @description
#' Derives one-sided (right-tailed for congruence; right/left/two-sided for
#' the gain/loss ratio) permutation p-values by comparing observed scores
#' against the null distribution from \code{perm_conservation_global}.
#'
#' @param observed_scores Named list; observed values of \code{congruence_index},
#'   \code{gain_index}, \code{loss_index}, and \code{gain_loss_ratio} from
#'   \code{congruence}.
#' @param perm_results B x 4 matrix from \code{perm_conservation_global}.
#'
#' @return A list with elements \code{congruence} (right-tailed p-value),
#'   \code{gain_index}, \code{loss_index} (observed values), and
#'   \code{ratio} (list of right/left/two-sided p-values for the ratio).
#'
#' @export
# Global p-value calculation
# For congruence_index: right-sided test only
# For gain_loss_ratio: right-sided, left-sided, and two-sided tests
p_conservation_global <- function(observed_scores, perm_results) {
  
  # Initialize p-value lists
  # For congruence_index: only right-sided
  p_congruence <- NA
  
  # For gain_loss_ratio: three types of tests
  p_ratio <- list(
    PValue_Right = NA,
    PValue_Left = NA,
    PValue_TwoSided = NA
  )
  
  # Congruence index - right-sided test (higher is better)
  if(!is.na(observed_scores$congruence_index)) {
    null_scores <- perm_results[, "congruence_index"]
    null_scores <- null_scores[!is.na(null_scores)]
    
    if(length(null_scores) > 0) {
      p_congruence <- (sum(null_scores >= observed_scores$congruence_index) + 1) / 
        (length(null_scores) + 1)
    }
  }
  
  # Gain/loss ratio - three types of tests
  if(!is.na(observed_scores$gain_loss_ratio)) {
    null_scores <- perm_results[, "gain_loss_ratio"]
    null_scores <- null_scores[!is.na(null_scores)]
    
    if(length(null_scores) > 0) {
      obs_val <- observed_scores$gain_loss_ratio
      
      # Right-sided: tests if observed ratio is higher than expected
      p_right <- (sum(null_scores >= obs_val) + 1) / (length(null_scores) + 1)
      
      # Left-sided: tests if observed ratio is lower than expected
      p_left <- (sum(null_scores <= obs_val) + 1) / (length(null_scores) + 1)
      
      # Two-sided: tests if observed ratio differs from expected (either direction)
      p_two <- 2 * min(p_right, p_left)
      p_two <- min(p_two, 1)  # Cap at 1
      
      p_ratio$PValue_Right <- p_right
      p_ratio$PValue_Left <- p_left
      p_ratio$PValue_TwoSided <- p_two
    }
  }
  
  return(list(
    congruence = p_congruence,
    gain_index = observed_scores$gain_index,
    loss_index = observed_scores$loss_index,
    ratio = p_ratio
  ))
}

#' Pairwise global conservation analysis across multiple conditions
#'
#' @title Multi-condition global circadian conservation analysis
#'
#' @description
#' Runs all pairwise comparisons among a list of MCMC outputs: for each
#' pair computes the observed congruence, gain, and loss indices; performs
#' permutation tests via \code{perm_conservation_global}; applies BH
#' multiple-testing correction; and writes an Excel summary with one sheet
#' per metric.
#'
#' @param mcmc.merge.list Named list of MCMC output lists, one per condition.
#' @param dataset.names Character vector; labels for each condition.
#' @param delta Numeric; phase-shift threshold (default 3).
#' @param units Character; \code{"hours"} (default).
#' @param B Integer; permutations per pair (default 1000).
#' @param ncores Integer or \code{NULL}; parallel cores.
#' @param parallel Logical or \code{"auto"}; parallelisation strategy.
#' @param rounds Integer; permutation rounds per pair.
#' @param save_intermediate Logical; save intermediate permutation results.
#' @param intermediate_dir Character; directory for intermediate saves.
#' @param output.dir Character; directory for the Excel output file
#'   (default \code{"Conservation_Global"}).
#'
#' @return Named list of square matrices (one per metric): congruence scores
#'   and p-values, gain and loss indices, ratio scores and p-values.
#'   Also writes an Excel file to \code{output.dir}.
#'
#' @export
# Main global conservation function with rounds support and Excel output
multi_conservation_global <- function(mcmc.merge.list, dataset.names,
                                      delta = 3, units = "hours", B = 1000,
                                      ncores = NULL, parallel = "auto",
                                      rounds = 1, save_intermediate = FALSE,
                                      intermediate_dir = "intermediate_results",
                                      output.dir = "Conservation_Global") {
  
  # Check if required package is installed
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required but not installed. Please install it with: install.packages('openxlsx')")
  }
  
  # Auto-detect best approach - convert to logical
  if (parallel == "auto") {
    parallel <- .Platform$OS.type != "windows" && B >= 100
  }
  
  # Ensure parallel is logical
  parallel <- as.logical(parallel)
  
  # Set default ncores if using parallel
  if (is.null(ncores)) {
    ncores <- min(parallel::detectCores(logical = FALSE), 4)
  }
  
  # Inform user about computational strategy
  if (parallel && .Platform$OS.type != "windows") {
    cat("Using parallel processing with", ncores, "cores for permutations\n")
  } else if (.Platform$OS.type == "windows") {
    cat("Using sequential processing (Windows detected)\n")
  } else {
    cat("Using sequential processing (parallel=FALSE or B<100)\n")
  }
  
  if (rounds > 1) {
    cat("Using", rounds, "rounds for permutations\n")
    if (save_intermediate) {
      cat("Intermediate results will be saved to:", intermediate_dir, "\n")
    }
  }
  
  names(mcmc.merge.list) <- dataset.names
  M <- length(mcmc.merge.list)
  
  # Initialize results matrices
  # Congruence index
  Congruence_scores <- matrix(NA, M, M)
  Congruence_pvalues <- matrix(NA, M, M)
  Congruence_qvalues <- matrix(NA, M, M)
  
  # Gain/loss components and ratio
  Gain_index <- matrix(NA, M, M)
  Loss_index <- matrix(NA, M, M)
  Ratio_scores <- matrix(NA, M, M)
  Ratio_pvalues_right <- matrix(NA, M, M)
  Ratio_qvalues_right <- matrix(NA, M, M)
  Ratio_pvalues_left <- matrix(NA, M, M)
  Ratio_qvalues_left <- matrix(NA, M, M)
  Ratio_pvalues_two <- matrix(NA, M, M)
  Ratio_qvalues_two <- matrix(NA, M, M)
  
  # Set row and column names
  rownames(Congruence_scores) <- colnames(Congruence_scores) <- dataset.names
  rownames(Congruence_pvalues) <- colnames(Congruence_pvalues) <- dataset.names
  rownames(Congruence_qvalues) <- colnames(Congruence_qvalues) <- dataset.names
  rownames(Gain_index) <- colnames(Gain_index) <- dataset.names
  rownames(Loss_index) <- colnames(Loss_index) <- dataset.names
  rownames(Ratio_scores) <- colnames(Ratio_scores) <- dataset.names
  rownames(Ratio_pvalues_right) <- colnames(Ratio_pvalues_right) <- dataset.names
  rownames(Ratio_qvalues_right) <- colnames(Ratio_qvalues_right) <- dataset.names
  rownames(Ratio_pvalues_left) <- colnames(Ratio_pvalues_left) <- dataset.names
  rownames(Ratio_qvalues_left) <- colnames(Ratio_qvalues_left) <- dataset.names
  rownames(Ratio_pvalues_two) <- colnames(Ratio_pvalues_two) <- dataset.names
  rownames(Ratio_qvalues_two) <- colnames(Ratio_qvalues_two) <- dataset.names
  
  # Set diagonals (comparing dataset to itself)
  diag(Congruence_scores) <- 1
  diag(Congruence_pvalues) <- 0
  diag(Congruence_qvalues) <- 0
  diag(Gain_index) <- 0
  diag(Loss_index) <- 0
  diag(Ratio_scores) <- NA  # 0/0 is undefined
  diag(Ratio_pvalues_right) <- 0
  diag(Ratio_qvalues_right) <- 0
  diag(Ratio_pvalues_left) <- 0
  diag(Ratio_qvalues_left) <- 0
  diag(Ratio_pvalues_two) <- 0
  diag(Ratio_qvalues_two) <- 0
  
  total_pairs <- choose(M, 2)
  current_pair <- 1
  
  for (i in 1:(M-1)) {
    for (j in (i+1):M) {
      cat("\n=== Processing pair", current_pair, "of", total_pairs, ":", 
          dataset.names[i], "vs", dataset.names[j], "===\n")
      
      dat1 <- mcmc.merge.list[[i]]
      dat2 <- mcmc.merge.list[[j]]
      
      # Calculate observed conservation scores
      conservation_result <- congruence(dat1, dat2, delta = delta, units = units)
      
      # Calculate permutation results with rounds
      cat("Running global permutations (B =", B, ")...\n")
      
      # Create pair-specific intermediate directory if saving intermediate results
      pair_intermediate_dir <- NULL
      if (save_intermediate) {
        pair_intermediate_dir <- file.path(intermediate_dir, paste0("pair_", current_pair, "_", 
                                                                    dataset.names[i], "_vs_", dataset.names[j]))
      }
      
      perm_results <- perm_conservation_global(
        dat1, dat2, delta = delta, units = units, B = B,
        ncores = ncores, parallel = parallel,
        rounds = rounds, save_intermediate = save_intermediate,
        intermediate_dir = pair_intermediate_dir
      )
      
      # Calculate p-values
      p_values <- p_conservation_global(conservation_result, perm_results)
      
      # Store congruence results (symmetric matrix)
      Congruence_scores[i, j] <- Congruence_scores[j, i] <- conservation_result$congruence_index
      Congruence_pvalues[i, j] <- Congruence_pvalues[j, i] <- p_values$congruence
      
      # Store gain/loss ratio results (symmetric matrix)
      Gain_index[i, j] <- Gain_index[j, i] <- p_values$gain_index
      Loss_index[i, j] <- Loss_index[j, i] <- p_values$loss_index
      Ratio_scores[i, j] <- Ratio_scores[j, i] <- conservation_result$gain_loss_ratio
      Ratio_pvalues_right[i, j] <- Ratio_pvalues_right[j, i] <- p_values$ratio$PValue_Right
      Ratio_pvalues_left[i, j] <- Ratio_pvalues_left[j, i] <- p_values$ratio$PValue_Left
      Ratio_pvalues_two[i, j] <- Ratio_pvalues_two[j, i] <- p_values$ratio$PValue_TwoSided
      
      current_pair <- current_pair + 1
    }
  }
  
  # Apply FDR correction AFTER all p-values are calculated
  cat("\nApplying FDR correction across all", total_pairs, "tests...\n")
  
  # Congruence index: single FDR correction
  p_flat <- Congruence_pvalues[upper.tri(Congruence_pvalues)]
  q_flat <- p.adjust(p_flat, method = "fdr")
  Congruence_qvalues[upper.tri(Congruence_qvalues)] <- q_flat
  Congruence_qvalues[lower.tri(Congruence_qvalues)] <- t(Congruence_qvalues)[lower.tri(Congruence_qvalues)]
  
  # Gain/loss ratio: separate FDR correction for each test type
  p_flat_right <- Ratio_pvalues_right[upper.tri(Ratio_pvalues_right)]
  q_flat_right <- p.adjust(p_flat_right, method = "fdr")
  Ratio_qvalues_right[upper.tri(Ratio_qvalues_right)] <- q_flat_right
  Ratio_qvalues_right[lower.tri(Ratio_qvalues_right)] <- t(Ratio_qvalues_right)[lower.tri(Ratio_qvalues_right)]
  
  p_flat_left <- Ratio_pvalues_left[upper.tri(Ratio_pvalues_left)]
  q_flat_left <- p.adjust(p_flat_left, method = "fdr")
  Ratio_qvalues_left[upper.tri(Ratio_qvalues_left)] <- q_flat_left
  Ratio_qvalues_left[lower.tri(Ratio_qvalues_left)] <- t(Ratio_qvalues_left)[lower.tri(Ratio_qvalues_left)]
  
  p_flat_two <- Ratio_pvalues_two[upper.tri(Ratio_pvalues_two)]
  q_flat_two <- p.adjust(p_flat_two, method = "fdr")
  Ratio_qvalues_two[upper.tri(Ratio_qvalues_two)] <- q_flat_two
  Ratio_qvalues_two[lower.tri(Ratio_qvalues_two)] <- t(Ratio_qvalues_two)[lower.tri(Ratio_qvalues_two)]
  
  # Create output directory
  if (!file.exists(output.dir)) dir.create(output.dir)
  
  # Create Excel workbook with 2 sheets
  cat("Creating Excel workbook with 2 sheets...\n")
  
  wb <- openxlsx::createWorkbook()
  
  # ========== SHEET 1: CONGRUENCE INDEX ==========
  openxlsx::addWorksheet(wb, sheetName = "congruence_index")
  
  result_df_cong <- data.frame(
    Dataset = dataset.names,
    stringsAsFactors = FALSE
  )
  
  for (j in 1:M) {
    result_df_cong[[paste0(dataset.names[j], "_Score")]] <- Congruence_scores[, j]
    result_df_cong[[paste0(dataset.names[j], "_PValue")]] <- Congruence_pvalues[, j]
    result_df_cong[[paste0(dataset.names[j], "_PValue_Adj")]] <- Congruence_qvalues[, j]
  }
  
  openxlsx::writeData(wb, sheet = "congruence_index", x = result_df_cong, startRow = 1, startCol = 1)
  
  # Format congruence sheet
  header_style <- openxlsx::createStyle(textDecoration = "bold", 
                                        fgFill = "#E6E6FA",
                                        border = "TopBottomLeftRight")
  openxlsx::addStyle(wb, sheet = "congruence_index", style = header_style, 
                     rows = 1, cols = 1:ncol(result_df_cong), gridExpand = TRUE)
  
  pvalue_cols <- grep("PValue", names(result_df_cong))
  if (length(pvalue_cols) > 0) {
    pvalue_style <- openxlsx::createStyle(numFmt = "0.0000")
    openxlsx::addStyle(wb, sheet = "congruence_index", style = pvalue_style, 
                       rows = 2:(nrow(result_df_cong) + 1), cols = pvalue_cols, gridExpand = TRUE)
  }
  
  openxlsx::setColWidths(wb, sheet = "congruence_index", cols = 1:ncol(result_df_cong), widths = "auto")
  
  # ========== SHEET 2: GAIN/LOSS RATIO ==========
  openxlsx::addWorksheet(wb, sheetName = "gain_loss_ratio")
  
  result_df_ratio <- data.frame(
    Dataset = dataset.names,
    stringsAsFactors = FALSE
  )
  
  for (j in 1:M) {
    result_df_ratio[[paste0(dataset.names[j], "_Gain_Index")]] <- Gain_index[, j]
    result_df_ratio[[paste0(dataset.names[j], "_Loss_Index")]] <- Loss_index[, j]
    result_df_ratio[[paste0(dataset.names[j], "_Gain_Loss_Ratio")]] <- Ratio_scores[, j]
    result_df_ratio[[paste0(dataset.names[j], "_PValue_Right")]] <- Ratio_pvalues_right[, j]
    result_df_ratio[[paste0(dataset.names[j], "_QValue_Right")]] <- Ratio_qvalues_right[, j]
    result_df_ratio[[paste0(dataset.names[j], "_PValue_Left")]] <- Ratio_pvalues_left[, j]
    result_df_ratio[[paste0(dataset.names[j], "_QValue_Left")]] <- Ratio_qvalues_left[, j]
    result_df_ratio[[paste0(dataset.names[j], "_PValue_TwoSided")]] <- Ratio_pvalues_two[, j]
    result_df_ratio[[paste0(dataset.names[j], "_QValue_TwoSided")]] <- Ratio_qvalues_two[, j]
  }
  
  openxlsx::writeData(wb, sheet = "gain_loss_ratio", x = result_df_ratio, startRow = 1, startCol = 1)
  
  # Format ratio sheet
  openxlsx::addStyle(wb, sheet = "gain_loss_ratio", style = header_style, 
                     rows = 1, cols = 1:ncol(result_df_ratio), gridExpand = TRUE)
  
  pvalue_cols_ratio <- grep("PValue|QValue", names(result_df_ratio))
  if (length(pvalue_cols_ratio) > 0) {
    pvalue_style <- openxlsx::createStyle(numFmt = "0.0000")
    openxlsx::addStyle(wb, sheet = "gain_loss_ratio", style = pvalue_style, 
                       rows = 2:(nrow(result_df_ratio) + 1), cols = pvalue_cols_ratio, gridExpand = TRUE)
  }
  
  openxlsx::setColWidths(wb, sheet = "gain_loss_ratio", cols = 1:ncol(result_df_ratio), widths = "auto")
  
  # Save the workbook
  excel_filename <- paste0(output.dir, "/Conservation_Results_Global_", M, "_datasets.xlsx")
  openxlsx::saveWorkbook(wb, file = excel_filename, overwrite = TRUE)
  
  # Also save CSV files for compatibility
  write.csv(Congruence_scores, 
            file = paste0(output.dir, "/Conservation_global_congruence_index_", M, ".csv"))
  write.csv(Congruence_pvalues, 
            file = paste0(output.dir, "/PValues_global_congruence_index_", M, ".csv"))
  write.csv(Congruence_qvalues, 
            file = paste0(output.dir, "/PValues_adj_global_congruence_index_", M, ".csv"))
  
  write.csv(Ratio_scores, 
            file = paste0(output.dir, "/Conservation_global_gain_loss_ratio_", M, ".csv"))
  write.csv(Ratio_pvalues_right, 
            file = paste0(output.dir, "/PValues_global_gain_loss_ratio_right_", M, ".csv"))
  write.csv(Ratio_qvalues_right, 
            file = paste0(output.dir, "/PValues_adj_global_gain_loss_ratio_right_", M, ".csv"))
  
  cat("\nResults saved to:", excel_filename, "\n")
  cat("Excel file contains 2 sheets:\n")
  cat("  1. congruence_index - conservation scores\n")
  cat("  2. gain_loss_ratio - gain/loss indices and ratio with multiple test types\n")
  cat("CSV files also created for compatibility.\n")
  
  return(list(
    Congruence_scores = Congruence_scores,
    Congruence_pvalues = Congruence_pvalues,
    Congruence_qvalues = Congruence_qvalues,
    Gain_index = Gain_index,
    Loss_index = Loss_index,
    Ratio_scores = Ratio_scores,
    Ratio_pvalues_right = Ratio_pvalues_right,
    Ratio_qvalues_right = Ratio_qvalues_right,
    Ratio_pvalues_left = Ratio_pvalues_left,
    Ratio_qvalues_left = Ratio_qvalues_left,
    Ratio_pvalues_two = Ratio_pvalues_two,
    Ratio_qvalues_two = Ratio_qvalues_two,
    delta = delta,
    units = units,
    B = B,
    rounds = rounds,
    ncores_used = ncores,
    parallel_used = parallel,
    output_file = excel_filename
  ))
}
#' Compute scalar Jaccard concordance summaries from posterior rho matrices
#'
#' @description
#' Wraps \code{compute_iteration_jaccard} to return mean observed, adjusted,
#' and null Jaccard values as a compact list for pairwise condition summaries.
#'
#' @param rho_A G x K binary matrix; posterior rho samples for condition A.
#' @param rho_B G x K binary matrix; posterior rho samples for condition B.
#'
#' @return List with elements \code{jaccard_obs}, \code{jaccard_adj}, and
#'   \code{jaccard_null} (all numeric scalars).
#'
#' @export
# Minimal concordance summary from posterior ρ matrices
# Wraps compute_iteration_jaccard to return scalar summaries per condition pair.
compute_concordance_minimal <- function(rho_A, rho_B) {
  res <- compute_iteration_jaccard(rho_A, rho_B)
  list(
    jaccard_obs  = mean(res$jaccard_obs,  na.rm = TRUE),
    jaccard_adj  = mean(res$jaccard_adj,  na.rm = TRUE),
    jaccard_null = mean(res$jaccard_exp,  na.rm = TRUE)
  )
}
