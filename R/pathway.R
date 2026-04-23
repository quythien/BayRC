
#' Compute observed congruence scores for each pathway
#'
#' @title Pathway-level circadian conservation scores
#'
#' @description
#' For each pathway in \code{select.pathway.list}, restricts the MCMC
#' posterior rho matrices to genes in that pathway and computes the
#' probabilistic congruence, gain, and loss indices using
#' \code{congruence}.
#'
#' @param dat1 Named list; MCMC output for condition 1 with a \code{rho}
#'   matrix and \code{attr(rho, "symbols")} set.
#' @param dat2 Named list; MCMC output for condition 2.
#' @param select.pathway.list Named list of character vectors; pathway
#'   gene sets (e.g. from MSigDB or KEGG).
#' @param delta Numeric; phase-shift threshold passed to \code{congruence}
#'   (default 3).
#' @param units Character; \code{"hours"} (default).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{scores}{K x 4 matrix of observed pathway congruence,
#'       gain, loss, and ratio scores.}
#'     \item{pathway_sizes}{Named integer vector; number of data genes
#'       in each pathway.}
#'     \item{expected_unions}{Named numeric vector; expected union size for
#'       each pathway.}
#'   }
#'
#' @export
# Observed pathway conservation scores with pathway metadata
conservation_pathway_circadian <- function(dat1, dat2, select.pathway.list,
                                           delta = 3, units = "hours") {
  
  pathway_names <- names(select.pathway.list)
  data_genes <- attr(dat1$rho, "symbols")
  K <- length(select.pathway.list)
  
  # Metrics to calculate
  metrics <- c("congruence_index", "gain_index", "loss_index", "gain_loss_ratio")
  # Initialize output for metrics
  conservation_scores <- matrix(NA, nrow = K, ncol = length(metrics))
  rownames(conservation_scores) <- pathway_names
  colnames(conservation_scores) <- metrics
  
  # Initialize pathway metadata
  pathway_sizes <- numeric(K)
  expected_unions <- numeric(K)
  names(pathway_sizes) <- pathway_names
  names(expected_unions) <- pathway_names
  
  for (k in 1:K) {
    pathway_genes <- select.pathway.list[[k]]
    overlap_genes <- intersect(data_genes, pathway_genes)
    
    if (length(overlap_genes) == 0) {
      conservation_scores[k, ] <- NA
      pathway_sizes[k] <- 0
      expected_unions[k] <- NA
      next
    }
    
    # Store pathway size
    pathway_sizes[k] <- length(overlap_genes)
    
    # Get indices for pathway genes
    gene_indices <- match(overlap_genes, data_genes)
    
    # Create matrices for this pathway (only rho needed)
    matrix1 <- list(
      rho = dat1$rho[gene_indices, , drop = FALSE]
    )
    matrix2 <- list(
      rho = dat2$rho[gene_indices, , drop = FALSE]
    )
    
    # Calculate conservation score
    cons_result <- congruence(matrix1, matrix2, delta = delta, units = units)
    
    # Store all metrics
    conservation_scores[k, "congruence_index"] <- cons_result$congruence_index
    conservation_scores[k, "gain_index"] <- cons_result$gain_index
    conservation_scores[k, "loss_index"] <- cons_result$loss_index
    conservation_scores[k, "gain_loss_ratio"] <- cons_result$gain_loss_ratio
    
    # Calculate and store union information
    # Recreate the union calculation from congruence function
    p_A <- rowMeans(matrix1$rho)
    p_B <- rowMeans(matrix2$rho)
    intersection_probs <- p_A * p_B
    expected_intersection <- sum(intersection_probs)
    expected_union <- sum(p_A) + sum(p_B) - expected_intersection
    
    expected_unions[k] <- expected_union
  }
  
  return(list(
    scores = conservation_scores,
    pathway_sizes = pathway_sizes,
    expected_unions = expected_unions
  ))
}

#' Generate permutation null distribution for pathway conservation scores
#'
#' @description
#' Runs gene-label permutations within each pathway to build a null
#' distribution for the congruence, gain, and loss indices.  Calls
#' \code{run_single_round} internally and supports parallel execution and
#' multi-round splits to manage memory.
#'
#' @param dat1,dat2 Named MCMC output lists with \code{rho} matrices and
#'   gene symbols set via \code{match_symbols}.
#' @param select.pathway.list Named list of character vectors; pathway
#'   gene sets.
#' @param delta Numeric; phase threshold (default 3).
#' @param units Character; \code{"hours"} (default).
#' @param B Integer; number of permutations (default 1000).
#' @param ncores Integer or \code{NULL}; parallel cores.
#' @param parallel Logical or \code{"auto"}; parallelisation strategy.
#' @param rounds Integer; rounds to split B permutations (default 1).
#' @param save_intermediate Logical; save each round's output to disk.
#' @param intermediate_dir Character; directory for intermediate saves.
#'
#' @return A named list of K arrays (one per pathway), each B x 4
#'   containing permuted congruence, gain, loss, and ratio scores.
#'
#' @export
# Permutation for pathway analysis with improved progress tracking and rounds support
perm_pathway_circadian <- function(dat1, dat2, select.pathway.list,
                                   delta = 3, units = "hours", B = 1000,
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
  total_B <- B_per_round * rounds  # Might be slightly larger than B
  
  cat("Splitting", B, "permutations into", rounds, "rounds of", B_per_round, "each\n")
  cat("Total permutations:", total_B, "\n")
  
  # Setup intermediate saving if requested
  if (save_intermediate && !is.null(intermediate_dir)) {
    if (!dir.exists(intermediate_dir)) {
      dir.create(intermediate_dir, recursive = TRUE)
    }
  }
  
  pathway_names <- names(select.pathway.list)
  data_genes <- attr(dat1$rho, "symbols")
  G <- nrow(dat1$rho)
  K <- length(select.pathway.list)
  
  # Get pathway gene indices for BOTH datasets (same genes, different order for B)
  pathway_indices_A <- vector("list", K)
  pathway_indices_B <- vector("list", K) 
  pathway_sizes <- numeric(K)
  
  for (k in 1:K) {
    pathway_genes <- select.pathway.list[[k]]
    overlap_genes <- intersect(data_genes, pathway_genes)
    gene_indices <- match(overlap_genes, data_genes)
    
    pathway_indices_A[[k]] <- gene_indices  # Fixed order for A
    pathway_indices_B[[k]] <- gene_indices  # Will be permuted for B
    pathway_sizes[k] <- length(overlap_genes)
  }
  
  # Metrics to store in permutation results
  metrics <- c("congruence_index", "gain_index", "loss_index", "gain_loss_ratio")
  
  # Create output array for ALL metrics (using total_B for consistent dimensions)
  out <- array(NA, dim = c(total_B, K, length(metrics)), 
               dimnames = list(NULL, pathway_names, metrics))
  
  # Process each round
  for (round_num in 1:rounds) {
    cat("\n--- Round", round_num, "of", rounds, "---\n")
    
    # Calculate start and end indices for this round
    start_idx <- (round_num - 1) * B_per_round + 1
    end_idx <- min(round_num * B_per_round, total_B)
    current_B <- end_idx - start_idx + 1
    
    cat("Processing permutations", start_idx, "to", end_idx, "(", current_B, "permutations)\n")
    
    # Run permutations for this round
    round_results <- run_single_round(
      dat1, dat2, pathway_indices_A, pathway_indices_B, pathway_sizes,
      K, metrics, delta, units, current_B, parallel, ncores
    )
    
    # Store results
    out[start_idx:end_idx, , ] <- round_results
    
    # Save intermediate results if requested
    if (save_intermediate && !is.null(intermediate_dir)) {
      round_file <- file.path(intermediate_dir, paste0("round_", round_num, "_results.rds"))
      saveRDS(round_results, round_file)
      cat("Intermediate results saved to:", round_file, "\n")
    }
    
    # Memory cleanup
    gc()
    
    cat("Round", round_num, "completed\n")
  }
  
  # If we generated more permutations than requested, trim to exact B
  if (total_B > B) {
    cat("Trimming results to exactly", B, "permutations\n")
    out <- out[1:B, , , drop = FALSE]
  }
  
  return(out)
}

#' Run one round of pathway-level permutations
#'
#' @description
#' Executes \code{B_round} gene-label permutations for all pathways in
#' one round, returning a B_round x K x 4 array of permuted scores.
#' Called by \code{perm_pathway_circadian}.
#'
#' @param dat1,dat2 MCMC output lists with \code{rho} matrices.
#' @param pathway_indices_A,pathway_indices_B Lists of integer vectors;
#'   row indices into the rho matrices for each pathway's genes.
#' @param pathway_sizes Integer vector; number of genes per pathway.
#' @param K Integer; number of pathways.
#' @param metrics Character vector; metric names.
#' @param delta,units Passed to \code{congruence}.
#' @param B_round Integer; permutations in this round.
#' @param parallel,ncores Parallelisation controls.
#'
#' @return B_round x K x length(metrics) numeric array.
#'
#' @keywords internal
# Helper function to run a single round of permutations
run_single_round <- function(dat1, dat2, pathway_indices_A, pathway_indices_B,
                             pathway_sizes, K, metrics, delta, units, B_round,
                             parallel, ncores) {
  
  # Create output array for this round
  round_out <- array(NA, dim = c(B_round, K, length(metrics)), 
                     dimnames = list(NULL, NULL, metrics))
  
  if (parallel && .Platform$OS.type != "windows") {
    
    # PARALLEL VERSION (Mac/Linux)
    library(parallel)
    
    perm_results <- mclapply(1:B_round, function(b) {
      
      # Initialize results for this permutation
      perm_b_results <- array(NA, dim = c(1, K, length(metrics)),
                              dimnames = list(NULL, NULL, metrics))
      
      for (k in 1:K) {
        if (pathway_sizes[k] > 0) {
          # FIXED indices for dataset A
          fixed_indices_A <- pathway_indices_A[[k]]
          
          # PERMUTED indices for dataset B
          permuted_indices_B <- sample(pathway_indices_B[[k]], pathway_sizes[k], replace = FALSE)
          
          # Only rho needed for congruence calculation
          matrix1 <- list(
            rho = dat1$rho[fixed_indices_A, , drop = FALSE]
          )
          matrix2 <- list(
            rho = dat2$rho[permuted_indices_B, , drop = FALSE]
          )
          
          # Calculate all metrics (including gain_loss_ratio from congruence)
          cons_result <- congruence(matrix1, matrix2, delta = delta, units = units)
          
          # Store all metrics
          perm_b_results[1, k, "congruence_index"] <- cons_result$congruence_index
          perm_b_results[1, k, "gain_index"] <- cons_result$gain_index
          perm_b_results[1, k, "loss_index"] <- cons_result$loss_index
          perm_b_results[1, k, "gain_loss_ratio"] <- cons_result$gain_loss_ratio
        }
      }
      
      return(perm_b_results)
      
    }, mc.cores = ncores)
    
    # Combine parallel results
    for (b in 1:B_round) {
      round_out[b, , ] <- perm_results[[b]][1, , ]
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
      
      for (k in 1:K) {
        if (pathway_sizes[k] > 0) {
          # FIXED indices for dataset A
          fixed_indices_A <- pathway_indices_A[[k]]
          
          # PERMUTED indices for dataset B
          permuted_indices_B <- sample(pathway_indices_B[[k]], pathway_sizes[k], replace = FALSE)
          
          # Only rho needed for congruence calculation
          matrix1 <- list(
            rho = dat1$rho[fixed_indices_A, , drop = FALSE]
          )
          matrix2 <- list(
            rho = dat2$rho[permuted_indices_B, , drop = FALSE]
          )
          
          # Calculate all metrics (including gain_loss_ratio from congruence)
          cons_result <- congruence(matrix1, matrix2, delta = delta, units = units)
          
          # Store all metrics
          round_out[b, k, "congruence_index"] <- cons_result$congruence_index
          round_out[b, k, "gain_index"] <- cons_result$gain_index
          round_out[b, k, "loss_index"] <- cons_result$loss_index
          round_out[b, k, "gain_loss_ratio"] <- cons_result$gain_loss_ratio
        }
      }
    }
  }
  
  return(round_out)
}

#' Compute permutation p-values for pathway conservation scores
#'
#' @description
#' Derives right-tailed (congruence) and three-sided (gain/loss ratio)
#' permutation p-values for each pathway by comparing observed scores from
#' \code{conservation_pathway_circadian} against the null distribution from
#' \code{perm_pathway_circadian}.
#'
#' @param observed_scores K x 4 matrix; observed pathway scores.
#' @param perm_results Named list of K arrays (B x 4); permutation null
#'   distributions.
#'
#' @return A data.frame with one row per pathway and columns for observed
#'   scores, p-values, and q-values (BH-adjusted).
#'
#' @export
# P-value calculation for pathway analysis
# For congruence_index: right-sided test only
# For gain_loss_ratio:


p_conservation_pathway <- function(observed_scores, perm_results) {
  K <- nrow(observed_scores)
  pathway_names <- rownames(observed_scores)
  
  # Initialize p-value matrices
  # For congruence_index: only right-sided
  p_congruence <- matrix(NA, nrow = K, ncol = 1)
  rownames(p_congruence) <- pathway_names
  colnames(p_congruence) <- "PValue"
  
  # For gain_loss_ratio: three types of tests
  p_ratio <- matrix(NA, nrow = K, ncol = 3)
  rownames(p_ratio) <- pathway_names
  colnames(p_ratio) <- c("PValue_Right", "PValue_Left", "PValue_TwoSided")
  
  
  for(k in 1:K) {
    # Congruence index - right-sided test (higher is better)
    if(!is.na(observed_scores[k, "congruence_index"])) {
      null_scores <- perm_results[, k, "congruence_index"]
      null_scores <- null_scores[!is.na(null_scores)]
      
      if(length(null_scores) > 0) {
        p_congruence[k, "PValue"] <- (sum(null_scores >= observed_scores[k, "congruence_index"]) + 1) / 
          (length(null_scores) + 1)
      }
    }
    
    # Gain/loss ratio - three types of tests 
    if(!is.na(observed_scores[k, "gain_loss_ratio"])) {
      null_scores <- perm_results[, k, "gain_loss_ratio"]
      null_scores <- null_scores[!is.na(null_scores)]
      
      if(length(null_scores) > 0) {
        obs_val <- observed_scores[k, "gain_loss_ratio"]
        
        # Right-sided: tests if GLR > 1 (more gain than loss)
        p_right <- (sum(null_scores >= obs_val) + 1) / (length(null_scores) + 1)
    
        # Left-sided: tests if GLR < 1 (more loss than gain)
        p_left <- (sum(null_scores <= obs_val) + 1) / (length(null_scores) + 1)
        
        
        # Two-sided: tests if GLR ≠ 1 (asymmetric in either direction)
        p_two <- 2 * min(p_right, p_left)
        p_two <- min(p_two, 1)
        
        p_ratio[k, "PValue_Right"] <- p_right
        p_ratio[k, "PValue_Left"] <- p_left
        p_ratio[k, "PValue_TwoSided"] <- p_two
        
      }
    }
    
  }
  return(list(
    congruence = p_congruence,
    ratio = p_ratio
  ))
}

# Main pathway function with rounds support
#' Pairwise pathway conservation analysis across multiple conditions
#'
#' @title Multi-condition pathway circadian conservation analysis
#'
#' @description
#' Runs all pairwise pathway conservation comparisons among a list of MCMC
#' outputs: computes observed pathway congruence scores, generates
#' permutation null distributions, applies BH correction, and writes
#' results to an Excel file.
#'
#' @param mcmc.merge.list Named list of MCMC output lists, one per condition.
#' @param dataset.names Character vector; condition labels.
#' @param select.pathway.list Named list of character vectors; pathway gene
#'   sets.
#' @param delta Numeric; phase threshold (default 3).
#' @param units Character; \code{"hours"} (default).
#' @param B Integer; permutations per pair (default 1000).
#' @param ncores Integer or \code{NULL}; parallel cores.
#' @param parallel Logical or \code{"auto"}.
#' @param rounds Integer; rounds per pair.
#' @param save_intermediate Logical; save intermediate results.
#' @param intermediate_dir Character; directory for intermediate saves.
#' @param output.dir Character; directory for Excel output
#'   (default \code{"Conservation_Pathway"}).
#'
#' @return Named list of pairwise pathway result data.frames.  Also writes
#'   an Excel file to \code{output.dir}.
#'
#' @export
multi_conservation_pathway <- function(mcmc.merge.list, dataset.names,
                                       select.pathway.list,
                                       delta = 3, units = "hours", B = 1000,
                                       ncores = NULL, parallel = "auto",
                                       rounds = 1, save_intermediate = FALSE,
                                       intermediate_dir = "intermediate_results",
                                       output.dir = "Conservation_Pathway") {
  
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
  P <- choose(M, 2)
  K <- length(select.pathway.list)
  
  # Storage for results
  # Congruence index
  Congruence_scores <- matrix(NA, nrow = K, ncol = P)
  Congruence_pvalues <- matrix(NA, nrow = K, ncol = P)
  Congruence_qvalues <- matrix(NA, nrow = K, ncol = P)
  
  # Gain/loss ratio and components
  Gain_index <- matrix(NA, nrow = K, ncol = P)
  Loss_index <- matrix(NA, nrow = K, ncol = P)
  Ratio_scores <- matrix(NA, nrow = K, ncol = P)
  Ratio_pvalues_right <- matrix(NA, nrow = K, ncol = P)
  Ratio_qvalues_right <- matrix(NA, nrow = K, ncol = P)
  Ratio_pvalues_left <- matrix(NA, nrow = K, ncol = P)
  Ratio_qvalues_left <- matrix(NA, nrow = K, ncol = P)
  Ratio_pvalues_two <- matrix(NA, nrow = K, ncol = P)
  Ratio_qvalues_two <- matrix(NA, nrow = K, ncol = P)
  
  # Pathway metadata (same for all comparisons)
  Pathway_sizes <- matrix(NA, nrow = K, ncol = P)
  Expected_unions <- matrix(NA, nrow = K, ncol = P)
  
  # Set row names
  pathway_names <- names(select.pathway.list)
  rownames(Congruence_scores) <- pathway_names
  rownames(Congruence_pvalues) <- pathway_names
  rownames(Congruence_qvalues) <- pathway_names
  rownames(Gain_index) <- pathway_names
  rownames(Loss_index) <- pathway_names
  rownames(Ratio_scores) <- pathway_names
  rownames(Ratio_pvalues_right) <- pathway_names
  rownames(Ratio_qvalues_right) <- pathway_names
  rownames(Ratio_pvalues_left) <- pathway_names
  rownames(Ratio_qvalues_left) <- pathway_names
  rownames(Ratio_pvalues_two) <- pathway_names
  rownames(Ratio_qvalues_two) <- pathway_names
  rownames(Pathway_sizes) <- pathway_names
  rownames(Expected_unions) <- pathway_names
  
  # Generate column names
  combinations <- combn(dataset.names, m = 2)
  pair_names <- apply(combinations, 2, FUN = function(x) paste(x, collapse = "&"))
  
  colnames(Congruence_scores) <- pair_names
  colnames(Congruence_pvalues) <- pair_names
  colnames(Congruence_qvalues) <- pair_names
  colnames(Gain_index) <- pair_names
  colnames(Loss_index) <- pair_names
  colnames(Ratio_scores) <- pair_names
  colnames(Ratio_pvalues_right) <- pair_names
  colnames(Ratio_qvalues_right) <- pair_names
  colnames(Ratio_pvalues_left) <- pair_names
  colnames(Ratio_qvalues_left) <- pair_names
  colnames(Ratio_pvalues_two) <- pair_names
  colnames(Ratio_qvalues_two) <- pair_names
  colnames(Pathway_sizes) <- pair_names
  colnames(Expected_unions) <- pair_names
  
  # Process each pair
  pair_idx <- 1
  total_pairs <- P
  
  for (i in 1:(M-1)) {
    for (j in (i+1):M) {
      cat("\n=== Processing pair", pair_idx, "of", total_pairs, ":", 
          dataset.names[i], "vs", dataset.names[j], "===\n")
      
      dat1 <- mcmc.merge.list[[i]]
      dat2 <- mcmc.merge.list[[j]]
      
      # Calculate observed conservation scores with metadata
      obs_results <- conservation_pathway_circadian(dat1, dat2, select.pathway.list,
                                                    delta = delta, units = units)
      
      conservation_scores <- obs_results$scores
      
      # Store metadata
      Pathway_sizes[, pair_idx] <- obs_results$pathway_sizes
      Expected_unions[, pair_idx] <- obs_results$expected_unions
      
      # Calculate permutation results with rounds
      cat("Running pathway permutations (B =", B, ")...\n")
      
      # Create pair-specific intermediate directory if saving intermediate results
      pair_intermediate_dir <- NULL
      if (save_intermediate) {
        pair_intermediate_dir <- file.path(intermediate_dir, paste0("pair_", pair_idx, "_", 
                                                                    gsub("&", "_vs_", pair_names[pair_idx])))
      }
      
      perm_results <- perm_pathway_circadian(
        dat1, dat2, select.pathway.list,
        delta = delta, units = units, B = B,
        ncores = ncores, parallel = parallel,
        rounds = rounds, save_intermediate = save_intermediate,
        intermediate_dir = pair_intermediate_dir
      )
      
      # Calculate p-values
      p_values_all <- p_conservation_pathway(conservation_scores, perm_results)
      
      # Store congruence results
      Congruence_scores[, pair_idx] <- conservation_scores[, "congruence_index"]
      Congruence_pvalues[, pair_idx] <- p_values_all$congruence[, "PValue"]
      
      # Store gain/loss ratio results
      Gain_index[, pair_idx] <- conservation_scores[, "gain_index"]
      Loss_index[, pair_idx] <- conservation_scores[, "loss_index"]
      Ratio_scores[, pair_idx] <- conservation_scores[, "gain_loss_ratio"]
      Ratio_pvalues_right[, pair_idx] <- p_values_all$ratio[, "PValue_Right"]
      Ratio_pvalues_left[, pair_idx] <- p_values_all$ratio[, "PValue_Left"]
      Ratio_pvalues_two[, pair_idx] <- p_values_all$ratio[, "PValue_TwoSided"]
      
      pair_idx <- pair_idx + 1
    }
  }
  
  # FDR CORRECTION
  cat("\nApplying FDR correction...\n")
  
  # Congruence index: single FDR correction
  p_flat <- as.vector(Congruence_pvalues)
  q_flat <- p.adjust(p_flat, method = "fdr")
  Congruence_qvalues <- matrix(q_flat, nrow = K, ncol = P)
  rownames(Congruence_qvalues) <- pathway_names
  colnames(Congruence_qvalues) <- pair_names
  
  # Gain/loss ratio: separate FDR correction for each test type
  p_flat_right <- as.vector(Ratio_pvalues_right)
  q_flat_right <- p.adjust(p_flat_right, method = "fdr")
  Ratio_qvalues_right <- matrix(q_flat_right, nrow = K, ncol = P)
  rownames(Ratio_qvalues_right) <- pathway_names
  colnames(Ratio_qvalues_right) <- pair_names
  
  p_flat_left <- as.vector(Ratio_pvalues_left)
  q_flat_left <- p.adjust(p_flat_left, method = "fdr")
  Ratio_qvalues_left <- matrix(q_flat_left, nrow = K, ncol = P)
  rownames(Ratio_qvalues_left) <- pathway_names
  colnames(Ratio_qvalues_left) <- pair_names
  
  p_flat_two <- as.vector(Ratio_pvalues_two)
  q_flat_two <- p.adjust(p_flat_two, method = "fdr")
  Ratio_qvalues_two <- matrix(q_flat_two, nrow = K, ncol = P)
  rownames(Ratio_qvalues_two) <- pathway_names
  colnames(Ratio_qvalues_two) <- pair_names
  
  # Create output directory
  if (!file.exists(output.dir)) dir.create(output.dir)
  
  # Create Excel workbook with 2 sheets
  cat("Creating Excel workbook with 2 sheets...\n")
  
  wb <- openxlsx::createWorkbook()
  
  # ========== SHEET 1: CONGRUENCE INDEX ==========
  openxlsx::addWorksheet(wb, sheetName = "congruence_index")
  
  result_df_cong <- data.frame(
    Pathway = pathway_names,
    stringsAsFactors = FALSE
  )
  
  for (pair_idx in 1:P) {
    pair_name <- pair_names[pair_idx]
    clean_pair_name <- gsub("&", "_vs_", pair_name)
    
    result_df_cong[[paste0(clean_pair_name, "_Pathway_Size")]] <- Pathway_sizes[, pair_idx]
    result_df_cong[[paste0(clean_pair_name, "_Expected_Union")]] <- Expected_unions[, pair_idx]
    result_df_cong[[paste0(clean_pair_name, "_Score")]] <- Congruence_scores[, pair_idx]
    result_df_cong[[paste0(clean_pair_name, "_PValue")]] <- Congruence_pvalues[, pair_idx]
    result_df_cong[[paste0(clean_pair_name, "_PValue_Adj")]] <- Congruence_qvalues[, pair_idx]
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
    Pathway = pathway_names,
    stringsAsFactors = FALSE
  )
  
  for (pair_idx in 1:P) {
    pair_name <- pair_names[pair_idx]
    clean_pair_name <- gsub("&", "_vs_", pair_name)
    
    result_df_ratio[[paste0(clean_pair_name, "_Pathway_Size")]] <- Pathway_sizes[, pair_idx]
    result_df_ratio[[paste0(clean_pair_name, "_Gain_Index")]] <- Gain_index[, pair_idx]
    result_df_ratio[[paste0(clean_pair_name, "_Loss_Index")]] <- Loss_index[, pair_idx]
    result_df_ratio[[paste0(clean_pair_name, "_Gain_Loss_Ratio")]] <- Ratio_scores[, pair_idx]
    result_df_ratio[[paste0(clean_pair_name, "_PValue_Right")]] <- Ratio_pvalues_right[, pair_idx]
    result_df_ratio[[paste0(clean_pair_name, "_QValue_Right")]] <- Ratio_qvalues_right[, pair_idx]
    result_df_ratio[[paste0(clean_pair_name, "_PValue_Left")]] <- Ratio_pvalues_left[, pair_idx]
    result_df_ratio[[paste0(clean_pair_name, "_QValue_Left")]] <- Ratio_qvalues_left[, pair_idx]
    result_df_ratio[[paste0(clean_pair_name, "_PValue_TwoSided")]] <- Ratio_pvalues_two[, pair_idx]
    result_df_ratio[[paste0(clean_pair_name, "_QValue_TwoSided")]] <- Ratio_qvalues_two[, pair_idx]
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
  excel_filename <- paste0(output.dir, "/Conservation_Results_", M, "_datasets.xlsx")
  openxlsx::saveWorkbook(wb, file = excel_filename, overwrite = TRUE)
  
  cat("\nResults saved to:", excel_filename, "\n")
  cat("Excel file contains 2 sheets:\n")
  cat("  1. congruence_index - conservation scores with pathway metadata\n")
  cat("  2. gain_loss_ratio - gain/loss indices and ratio with multiple test types\n")
  
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
    Pathway_sizes = Pathway_sizes,
    Expected_unions = Expected_unions,
    delta = delta,
    units = units,
    B = B,
    rounds = rounds,
    ncores_used = ncores,
    parallel_used = parallel,
    output_file = excel_filename
  ))
}
