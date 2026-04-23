#' Select circadian-enriched pathways using FGSEA
#'
#' @title FGSEA-based pathway selection for circadian transitions
#'
#' @description
#' Applies fast gene-set enrichment analysis (FGSEA) to a posterior
#' gene-ranking metric derived from two BayRC MCMC outputs.  Genes are
#' ranked by their log-odds-ratio of gain vs loss, or by conservation,
#' gain, loss, or union probability; pathways are tested for enrichment
#' toward either end of this ranking.  Pathway-level effect sizes and
#' probabilistic gain/loss/conservation indices are added to the output.
#'
#' @param mcmc.merge.list Named list of exactly 2 MCMC output lists; the
#'   first is condition A (reference) and the second is condition B.
#'   Each must have \code{attr(rho, "symbols")} set.
#' @param pathway.list Named list of character vectors; pathway gene sets.
#' @param dataset.names Character vector of length 2 or \code{NULL};
#'   labels for conditions A and B (default auto-detected from list names).
#' @param ranking.method Character; gene-ranking metric: one of
#'   \code{"log_odds_ratio"} (default), \code{"union"}, \code{"conserved"},
#'   \code{"gain"}, or \code{"loss"}.
#' @param score_type Character; fgsea score type: \code{"std"} (two-sided),
#'   \code{"pos"} (gain only), or \code{"neg"} (loss only, inverts ranking).
#' @param pathwaysize.lower.cut Integer; minimum pathway size (default 10).
#' @param pathwaysize.upper.cut Integer; maximum pathway size (default 200).
#' @param qvalue.cut Numeric; adjusted p-value significance cutoff
#'   (default 0.05).
#' @param epsilon Numeric; small constant added to probabilities to avoid
#'   log(0) (default 0.001).
#' @param n_top_genes Integer; number of top genes to report per category
#'   per pathway (default 5).
#' @param nperm Integer; fgsea permutations (default 10000).
#' @param nproc Integer; fgsea parallel processes (default 1).
#' @param seed Integer; random seed (default 12345).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{results}{Data.frame of fgsea results (one row per pathway)
#'       augmented with \code{Gain_Index}, \code{Loss_Index},
#'       \code{Conserved_Index}, \code{Gain_Loss_Ratio_Arithmetic},
#'       expected counts, and top gene columns.}
#'     \item{gene_rankings}{Data.frame of per-gene ranking metrics and
#'       posterior probabilities.}
#'     \item{parameters}{List of analysis settings.}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' res <- pathSelect(list(human = human_res, mouse = mouse_res),
#'                   pathway.list = msig_pathways,
#'                   ranking.method = "log_odds_ratio")
#' head(res$results)
#' }
pathSelect <- function(mcmc.merge.list,
                       pathway.list,
                       dataset.names = NULL,
                       ranking.method = c("log_odds_ratio", "union", "conserved", "gain", "loss"),
                       score_type = c("std", "pos", "neg"),
                       pathwaysize.lower.cut = 10,
                       pathwaysize.upper.cut = 200,
                       qvalue.cut = 0.05,
                       epsilon = 0.001,
                       n_top_genes = 5,
                       nperm = 10000,
                       nproc = 1,
                       seed = 12345) {
  
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required. Install with: BiocManager::install('fgsea')")
  }
  
  score_type <- match.arg(score_type)
  ranking.method <- match.arg(ranking.method)
  set.seed(seed)
  
  # Validate input
  if (length(mcmc.merge.list) != 2) {
    stop("mcmc.merge.list must contain exactly 2 datasets")
  }
  
  if (is.null(dataset.names)) {
    dataset.names <- names(mcmc.merge.list)
    if (is.null(dataset.names)) {
      dataset.names <- c("Condition_A", "Condition_B")
    }
  }
  
  # Extract data
  dat1 <- mcmc.merge.list[[1]]
  dat2 <- mcmc.merge.list[[2]]
  
  genes1 <- attr(dat1$rho, "symbols")
  genes2 <- attr(dat2$rho, "symbols")
  
  if (!identical(genes1, genes2)) {
    stop("Gene sets must be identical across datasets")
  }
  
  gene_names <- genes1
  n_genes <- length(gene_names)
  
  # Filter pathways by size
  pathway_sizes <- sapply(pathway.list, length)
  size_filter <- pathway_sizes >= pathwaysize.lower.cut & 
    pathway_sizes <= pathwaysize.upper.cut
  filtered_pathways <- pathway.list[size_filter]
  
  if (length(filtered_pathways) == 0) {
    stop("No pathways pass size filters")
  }
  
  # Determine if we need to invert ranking for loss-only tests
  # ONLY invert for log_odds_ratio, NOT for probability-based rankings!
  invert_ranking <- (score_type == "neg") && (ranking.method == "log_odds_ratio")
  
  cat("\n=== FGSEA PATHWAY ANALYSIS ===\n")
  cat("Comparing:", dataset.names[1], "vs", dataset.names[2], "\n")
  cat("Ranking method:", ranking.method, "\n")
  cat("Test type:",
      switch(score_type,
             "std" = "Two-sided (gain or loss)",
             "pos" = "One-sided (gain only)",
             "neg" = "One-sided (loss only)"), "\n")
  if (invert_ranking) {
    cat("Note: Rankings inverted for loss-only testing\n")
  }
  cat("Pathways:", length(filtered_pathways), "of", length(pathway.list), "\n")
  cat("Genes:", n_genes, "\n")
  cat("Permutations:", nperm, "\n\n")
  
  # ===== CALCULATE POSTERIOR PROBABILITIES =====
  cat("Calculating posterior probabilities...\n")
  
  # Calculate point estimates
  P_A <- rowMeans(dat1$rho)
  P_B <- rowMeans(dat2$rho)
  names(P_A) <- gene_names
  names(P_B) <- gene_names
  
  # Calculate transition probabilities (continuous, no thresholding)
  P_gain <- P_B * (1 - P_A)
  P_loss <- P_A * (1 - P_B)
  P_cons <- P_A * P_B
  P_union <- 1 - (1 - P_A) * (1 - P_B)
  
  # Summary statistics
  cat("Probabilistic rhythmic gene summary:\n")
  cat("  Expected rhythmic in", dataset.names[1], ":", round(sum(P_A), 1), "\n")
  cat("  Expected rhythmic in", dataset.names[2], ":", round(sum(P_B), 1), "\n")
  cat("  Expected gain:", round(sum(P_gain), 1), "\n")
  cat("  Expected loss:", round(sum(P_loss), 1), "\n")
  cat("  Expected conserved:", round(sum(P_cons), 1), "\n")
  cat("  Expected union:", round(sum(P_union), 1), "\n\n")
  
  # ===== HELPER FUNCTIONS =====
  
  # Helper: Calculate ranking metric with optional inversion
  calculate_ranking_metric <- function(P_A, P_B, method, epsilon = 0.001, invert = FALSE) {
    if (method == "log_odds_ratio") {
      loss_signal <- P_A * (1 - P_B) + epsilon
      gain_signal <- P_B * (1 - P_A) + epsilon
      metric <- log(gain_signal) - log(loss_signal)
      
      # If testing for loss only, invert the ranking
      if (invert) {
        metric <- -metric
      }
    } else if (method == "conserved") {
      # Conservation probability
      metric <- P_A * P_B + epsilon
    } else if (method == "union") {
      # Fisher-like: -2 * log(P(non-rhythmic in both))
      metric <- -2 * log((1 - P_A) * (1 - P_B) + epsilon)
    } else if (method == "gain") {
      # Gain probability
      metric <- P_B * (1 - P_A) + epsilon
    } else if (method == "loss") {
      # Loss probability
      metric <- P_A * (1 - P_B) + epsilon
    }
    return(metric)
  }
  
  # Helper: Calculate pathway effect size (GEOMETRIC MEAN)
  calculate_pathway_effect_size <- function(P_A, P_B, pathway_genes, method, epsilon = 0.001) {
    P_A_pathway <- P_A[pathway_genes]
    P_B_pathway <- P_B[pathway_genes]
    
    if (method == "log_odds_ratio") {
      loss_signal <- P_A_pathway * (1 - P_B_pathway) + epsilon
      gain_signal <- P_B_pathway * (1 - P_A_pathway) + epsilon
      log_OR <- log(gain_signal) - log(loss_signal)
      effect_size <- exp(mean(log_OR, na.rm = TRUE))
      effect_name <- "Geometric_Mean_OR"
    } else if (method == "conserved") {
      # Geometric mean of conservation probabilities
      conserved_probs <- P_A_pathway * P_B_pathway
      effect_size <- exp(mean(log(conserved_probs + epsilon), na.rm = TRUE))
      effect_name <- "Geometric_Mean_Conserved"
    } else if (method == "union") {
      # Geometric mean of union probabilities
      union_probs <- 1 - (1 - P_A_pathway) * (1 - P_B_pathway)
      effect_size <- exp(mean(log(union_probs + epsilon), na.rm = TRUE))
      effect_name <- "Geometric_Mean_Union"
    } else if (method == "gain") {
      # Geometric mean of gain probabilities
      gain_probs <- P_B_pathway * (1 - P_A_pathway)
      effect_size <- exp(mean(log(gain_probs + epsilon), na.rm = TRUE))
      effect_name <- "Geometric_Mean_Gain"
    } else if (method == "loss") {
      # Geometric mean of loss probabilities
      loss_probs <- P_A_pathway * (1 - P_B_pathway)
      effect_size <- exp(mean(log(loss_probs + epsilon), na.rm = TRUE))
      effect_name <- "Geometric_Mean_Loss"
    }
    
    return(list(value = effect_size, name = effect_name))
  }
  
  # Helper: Calculate gain/loss indices using CONTINUOUS probabilities
  calculate_gain_loss_indices <- function(pathway_genes, P_A, P_B, P_gain, P_loss, P_cons) {
    # Get indices for genes in this pathway
    pathway_idx <- match(pathway_genes, gene_names)
    pathway_idx <- pathway_idx[!is.na(pathway_idx)]
    
    if (length(pathway_idx) == 0) {
      return(list(
        gain_index = 0,
        loss_index = 0,
        congruence_index = 0,
        gain_loss_ratio_arithmetic = NA_real_,
        expected_n_gain = 0,
        expected_n_loss = 0,
        expected_n_conserved = 0
      ))
    }
    
    # Get probabilities for pathway genes
    P_A_pathway <- P_A[pathway_idx]
    P_B_pathway <- P_B[pathway_idx]
    P_gain_pathway <- P_gain[pathway_idx]
    P_loss_pathway <- P_loss[pathway_idx]
    P_cons_pathway <- P_cons[pathway_idx]
    
    # Calculate expected counts (continuous)
    expected_gain <- sum(P_gain_pathway)
    expected_loss <- sum(P_loss_pathway)
    expected_cons <- sum(P_cons_pathway)
    
    # Union for normalization
    expected_union <- sum(P_A_pathway) + sum(P_B_pathway) - expected_cons
    
    if (expected_union == 0) {
      return(list(
        gain_index = 0,
        loss_index = 0,
        congruence_index = 0,
        gain_loss_ratio_arithmetic = NA_real_,
        expected_n_gain = 0,
        expected_n_loss = 0,
        expected_n_conserved = 0
      ))
    }
    
    # Indices (normalized by union)
    gain_index <- expected_gain / expected_union
    loss_index <- expected_loss / expected_union
    congruence_index <- expected_cons / expected_union
    
    # Arithmetic ratio
    gain_loss_ratio_arithmetic <- if (expected_loss == 0) {
      if (expected_gain == 0) {
        NA_real_
      } else {
        1e10  # Gain only, no loss
      }
    } else {
      expected_gain / expected_loss
    }
    
    return(list(
      gain_index = gain_index,
      loss_index = loss_index,
      congruence_index = congruence_index,
      gain_loss_ratio_arithmetic = gain_loss_ratio_arithmetic,
      expected_n_gain = expected_gain,
      expected_n_loss = expected_loss,
      expected_n_conserved = expected_cons
    ))
  }
  
  # Helper: Get top genes by ranking transition probabilities (no thresholding)
  get_top_genes_by_category <- function(pathway_genes, P_A, P_B, 
                                        P_gain, P_loss, P_cons, n_top = 5) {
    # Get indices for genes in this pathway
    pathway_idx <- match(pathway_genes, gene_names)
    pathway_idx <- pathway_idx[!is.na(pathway_idx)]
    
    if (length(pathway_idx) == 0) {
      return(list(
        top_gain_genes = character(0),
        top_loss_genes = character(0),
        top_conserved_genes = character(0)
      ))
    }
    
    genes_in_pathway <- gene_names[pathway_idx]
    P_gain_pathway <- P_gain[pathway_idx]
    P_loss_pathway <- P_loss[pathway_idx]
    P_cons_pathway <- P_cons[pathway_idx]
    
    # Rank by gain probability
    gain_order <- order(P_gain_pathway, decreasing = TRUE)
    top_gain <- genes_in_pathway[gain_order][1:min(n_top, length(genes_in_pathway))]
    
    # Rank by loss probability
    loss_order <- order(P_loss_pathway, decreasing = TRUE)
    top_loss <- genes_in_pathway[loss_order][1:min(n_top, length(genes_in_pathway))]
    
    # Rank by conservation probability
    cons_order <- order(P_cons_pathway, decreasing = TRUE)
    top_cons <- genes_in_pathway[cons_order][1:min(n_top, length(genes_in_pathway))]
    
    return(list(
      top_gain_genes = top_gain,
      top_loss_genes = top_loss,
      top_conserved_genes = top_cons
    ))
  }
  
  # ===== MAIN ANALYSIS =====
  
  # Calculate ranking metric with optional inversion
  ranking_metric <- calculate_ranking_metric(P_A, P_B, ranking.method, epsilon, invert_ranking)
  names(ranking_metric) <- gene_names
  
  # Add reproducible noise to break ties
  set.seed(seed)
  noise_sd <- sd(ranking_metric, na.rm = TRUE) * 1e-6
  noise <- rnorm(n_genes, mean = 0, sd = noise_sd)
  ranking_metric <- ranking_metric + noise
  
  cat("Ranking statistics:\n")
  cat("  Range: [", round(min(ranking_metric), 4), ",", round(max(ranking_metric), 4), "]\n")
  if (ranking.method %in% c("union")) {
    cat("  Interpretation: Fisher-like statistic (-2*log(prob))\n")
  }
  cat("\n")
  
  # Run fgsea
  cat("Running fgsea...\n")
  if (ranking.method %in% c("union", "conserved", "gain", "loss")) {
    fgsea_scoreType <- "pos"
  } else if (ranking.method == "log_odds_ratio") {
    fgsea_scoreType <- switch(score_type,
                              "std" = "std",
                              "pos" = "pos",
                              "neg" = "pos") # 'neg' uses inverted ranking with 'pos' test
  }
  
  fgsea_results <- fgsea::fgsea(
    pathways = filtered_pathways,
    stats = ranking_metric,
    minSize = pathwaysize.lower.cut,
    maxSize = pathwaysize.upper.cut,
    nPermSimple = nperm,
    scoreType = fgsea_scoreType,
    nproc = nproc
  )
  
  # Convert to data.frame and remove unwanted columns
  cat("Adding custom metrics...\n")
  
  fgsea_results <- as.data.frame(fgsea_results)
  fgsea_results$nMoreExtreme <- NULL
  fgsea_results$leadingEdge <- NULL
  
  n_pathways <- nrow(fgsea_results)
  
  # Initialize columns for custom metrics
  pathway_effects <- numeric(n_pathways)
  gain_indices <- numeric(n_pathways)
  loss_indices <- numeric(n_pathways)
  conserved_indices <- numeric(n_pathways)
  gain_loss_ratios_arithmetic <- numeric(n_pathways)
  expected_n_gain <- numeric(n_pathways)
  expected_n_loss <- numeric(n_pathways)
  expected_n_conserved <- numeric(n_pathways)
  
  top_gain_genes_list <- vector("list", n_pathways)
  top_loss_genes_list <- vector("list", n_pathways)
  top_conserved_genes_list <- vector("list", n_pathways)
  
  for (k in 1:n_pathways) {
    pathway_name <- fgsea_results$pathway[k]
    pathway_genes <- intersect(filtered_pathways[[pathway_name]], gene_names)
    
    # Calculate pathway effect (GEOMETRIC MEAN)
    pathway_effect <- calculate_pathway_effect_size(P_A, P_B, pathway_genes, ranking.method, epsilon)
    pathway_effects[k] <- pathway_effect$value
    
    if (k == 1) effect_name <- pathway_effect$name
    
    # Calculate gain/loss indices using CONTINUOUS probabilities
    gl_indices <- calculate_gain_loss_indices(
      pathway_genes, P_A, P_B, P_gain, P_loss, P_cons
    )
    
    gain_indices[k] <- gl_indices$gain_index
    loss_indices[k] <- gl_indices$loss_index
    conserved_indices[k] <- gl_indices$congruence_index
    gain_loss_ratios_arithmetic[k] <- gl_indices$gain_loss_ratio_arithmetic
    expected_n_gain[k] <- gl_indices$expected_n_gain
    expected_n_loss[k] <- gl_indices$expected_n_loss
    expected_n_conserved[k] <- gl_indices$expected_n_conserved
    
    # Get top genes by category (ranked by transition probabilities)
    top_genes <- get_top_genes_by_category(
      pathway_genes, P_A, P_B, P_gain, P_loss, P_cons, n_top_genes
    )
    
    top_gain_genes_list[[k]] <- top_genes$top_gain_genes
    top_loss_genes_list[[k]] <- top_genes$top_loss_genes
    top_conserved_genes_list[[k]] <- top_genes$top_conserved_genes
  }
  
  # Add columns to fgsea results
  fgsea_results[[effect_name]] <- pathway_effects
  fgsea_results$Gain_Index <- gain_indices
  fgsea_results$Loss_Index <- loss_indices
  fgsea_results$Conserved_Index <- conserved_indices
  fgsea_results$Gain_Loss_Ratio_Arithmetic <- gain_loss_ratios_arithmetic
  fgsea_results$Expected_N_Gain <- expected_n_gain
  fgsea_results$Expected_N_Loss <- expected_n_loss
  fgsea_results$Expected_N_Conserved <- expected_n_conserved
  fgsea_results$Top_Gain_Genes <- sapply(top_gain_genes_list, function(x) paste(x, collapse = ", "))
  fgsea_results$Top_Loss_Genes <- sapply(top_loss_genes_list, function(x) paste(x, collapse = ", "))
  fgsea_results$Top_Conserved_Genes <- sapply(top_conserved_genes_list, function(x) paste(x, collapse = ", "))
  
  # Add method-specific columns
  if (ranking.method == "log_odds_ratio") {
    fgsea_results$Odds_Ratio_Geometric <- pathway_effects
    fgsea_results$Log_OR_Geometric <- log(pathway_effects)
  } else if (ranking.method == "conserved") {
    fgsea_results$Log_Conserved_Geometric <- log(pathway_effects)
  } else if (ranking.method == "union") {
    fgsea_results$Log_Union_Geometric <- log(pathway_effects)
  } else if (ranking.method == "gain") {
    fgsea_results$Log_Gain_Geometric <- log(pathway_effects)
  } else if (ranking.method == "loss") {
    fgsea_results$Log_Loss_Geometric <- log(pathway_effects)
  }
  
  # Add significance flag
  fgsea_results$Significant <- fgsea_results$padj < qvalue.cut
  
  # Sort by p-value
  fgsea_results <- fgsea_results[order(fgsea_results$pval), ]
  rownames(fgsea_results) <- NULL
  
  # Summary
  cat("\n=== RESULTS SUMMARY ===\n")
  cat("Significant pathways (Q <", qvalue.cut, "):", sum(fgsea_results$Significant, na.rm = TRUE), "\n")
  
  if (sum(fgsea_results$Significant, na.rm = TRUE) > 0) {
    cat("\nTop significant pathways:\n")
    sig_results <- head(fgsea_results[fgsea_results$Significant, ], 5)
    for (i in 1:nrow(sig_results)) {
      cat(sprintf("  %s: ES=%.3f, NES=%.3f, Q=%.4f\n",
                  sig_results$pathway[i],
                  sig_results$ES[i],
                  sig_results$NES[i],
                  sig_results$padj[i]))
      cat(sprintf("           %s=%.4f\n",
                  effect_name,
                  sig_results[[effect_name]][i]))
      cat(sprintf("           Expected: Gain=%.1f (%.3f), Loss=%.1f (%.3f), Conserved=%.1f (%.3f)\n",
                  sig_results$Expected_N_Gain[i],
                  sig_results$Gain_Index[i],
                  sig_results$Expected_N_Loss[i],
                  sig_results$Loss_Index[i],
                  sig_results$Expected_N_Conserved[i],
                  sig_results$Conserved_Index[i]))
    }
  }
  
  # Gene rankings dataframe with probabilistic metrics
  gene_rankings <- data.frame(
    gene = gene_names,
    ranking_metric = ranking_metric,
    P_A = P_A,
    P_B = P_B,
    P_gain = P_gain,
    P_loss = P_loss,
    P_conserved = P_cons,
    P_union = P_union,
    stringsAsFactors = FALSE
  )
  
  # Sort by ranking metric
  gene_rankings <- gene_rankings[order(-gene_rankings$ranking_metric), ]
  
  if (ranking.method == "log_odds_ratio") {
    # Store the ACTUAL log OR 
    actual_log_OR <- calculate_ranking_metric(
      gene_rankings$P_A, 
      gene_rankings$P_B, 
      "log_odds_ratio", 
      epsilon, 
      invert = FALSE
    )
    gene_rankings$log_OR <- actual_log_OR
    gene_rankings$OR <- exp(actual_log_OR)
  }
  
  return(list(
    results = fgsea_results,
    gene_rankings = gene_rankings,
    parameters = list(
      ranking_method = ranking.method,
      dataset_A = dataset.names[1],
      dataset_B = dataset.names[2],
      n_genes = n_genes,
      n_pathways_tested = nrow(fgsea_results),
      npermutations = nperm,
      qvalue_cutoff = qvalue.cut,
      score_type = score_type,
      ranking_inverted = invert_ranking,
      expected_rhythmic_A = sum(P_A),
      expected_rhythmic_B = sum(P_B),
      expected_gain = sum(P_gain),
      expected_loss = sum(P_loss),
      expected_conserved = sum(P_cons),
      expected_union = sum(P_union)
    )
  ))
}