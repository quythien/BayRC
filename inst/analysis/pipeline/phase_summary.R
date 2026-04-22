compare_phase_shift_thresholds <- function(
    phi_matrix1,
    phi_matrix2,
    gain_loss_status,
    shift_thresholds = c(4, 3 , 2, 1 , 0.5),
    bfdr_grid = seq(0.25, 0.05, by = -0.05),
    P = 24
) {
  require(circular)
  
  # Maintained genes only
  maintained_idx <- which(gain_loss_status == "Maintained")
  n_maintained <- length(maintained_idx)
  
  if (n_maintained == 0) {
    cat("No maintained genes.\n")
    return(NULL)
  }
  
  # Compute phdiff matrix once
  phase_diff_matrix <- ((phi_matrix1 - phi_matrix2 + P/2) %% P) - P/2
  
  # Initialize results
  results <- data.frame()
  
  for (bfdr_alpha in bfdr_grid) {
    
    for (shift_val in shift_thresholds) {
      
      prob_shift <- prob_conserved <- numeric(n_maintained)
      
      for (k in seq_along(maintained_idx)) {
        i <- maintained_idx[k]
        vec <- phase_diff_matrix[i, ]
        prob_shift[k] <- mean(abs(vec) >= shift_val)
        prob_conserved[k] <- 1 - prob_shift[k]
      }
      
      # ---------------------------
      # BFDR — Shifted
      # ---------------------------
      ord_shift <- order(-prob_shift)
      BFDR_shift_sorted <- cumsum(1 - prob_shift[ord_shift]) / seq_along(prob_shift)
      idx_shift <- ord_shift[BFDR_shift_sorted <= bfdr_alpha]
      
      # ---------------------------
      # BFDR — Conserved
      # ---------------------------
      ord_cons <- order(-prob_conserved)
      BFDR_cons_sorted <- cumsum(1 - prob_conserved[ord_cons]) / seq_along(prob_conserved)
      idx_cons <- ord_cons[BFDR_cons_sorted <= bfdr_alpha]
      
      # ---------------------------
      # Exclusive Categorization
      # ---------------------------
      overlap <- intersect(idx_shift, idx_cons)
      
      # Remove overlap → Undetermined
      if (length(overlap) > 0) {
        idx_shift <- setdiff(idx_shift, overlap)
        idx_cons <- setdiff(idx_cons, overlap)
      }
      
      idx_union <- union(idx_shift, idx_cons)
      idx_undetermined <- setdiff(seq_len(n_maintained), idx_union)
      
      # Counts
      N_Shifted <- length(idx_shift)
      N_Conserved <- length(idx_cons)
      N_Undetermined <- length(idx_undetermined)
      
      # Append to results
      results <- rbind(
        results,
        data.frame(
          BFDR_alpha = bfdr_alpha,
          Shift_Threshold = shift_val,
          N_Maintained = n_maintained,
          N_Shifted = N_Shifted,
          N_Conserved = N_Conserved,
          N_Undetermined = N_Undetermined,
          stringsAsFactors = FALSE
        )
      )
    }
  }
  
  return(results)
}

human_LUNG <- list(
  rho = mcmc_data_human$LUN,
  phi = mcmc_phi_human$LUN
)

baboon_LUNG <- list(
  rho = mcmc_data_baboon$LUN,
  phi = mcmc_phi_baboon$LUN
)


pA <- rowMeans(human_LUNG$rho)
pB <- rowMeans(baboon_LUNG$rho)

# Step 1: classify gain/loss first
trans_lun <- transition_classify(pA, pB, bfdr_alpha = 0.2)
circadian_genes <- c(
  "BHLHE40","BHLHE41","BMAL1","ARNTL","BTRC","CLOCK","CREB1",
  "CRY1","CRY2","CSNK1D","CSNK1E","CUL1","DBP","FBXL3","FBXW11",
  "NFIL3","NPAS2","NR1D1","NR1D2","PER1","PER2","PER3",
  "RORA","RORB","RORC"
)
gain_genes  <- names(trans_lun$gain_genes)
loss_genes  <- names(trans_lun$loss_genes)
cons_genes  <- names(trans_lun$cons_genes)
circ_gain  <- intersect(circadian_genes, gain_genes)
circ_loss  <- intersect(circadian_genes, loss_genes)
circ_cons  <- intersect(circadian_genes, cons_genes)

list(
  circadian_gain = circ_gain,
  circadian_loss = circ_loss,
  circadian_conserved = circ_cons
)

phase_lung <- compare_phase_shift_thresholds(
  phi_matrix1 = human_LUNG$phi,
  phi_matrix2 = baboon_LUNG$phi,
  gain_loss_status = trans_lun$gain_loss_status
)



######
compare_phase_shift_thresholds_circadian <- function(
    phi_matrix1,
    phi_matrix2,
    gain_loss_status,
    circadian_genes,
    gene_names,
    shift_thresholds = c(4, 3, 2, 1, 0.5),
    bfdr_grid = seq(0.25, 0.05, by = -0.05),
    P = 24
) {
  
  # Maintained genes
  maintained_idx <- which(gain_loss_status == "Maintained")
  maintained_genes <- gene_names[maintained_idx]
  
  # Restrict to circadian genes
  circ_idx <- maintained_idx[maintained_genes %in% circadian_genes]
  n_circ <- length(circ_idx)
  
  if (n_circ == 0) {
    cat("No maintained circadian genes.\n")
    return(NULL)
  }
  
  # Phase difference matrix
  phase_diff_matrix <- ((phi_matrix1 - phi_matrix2 + P/2) %% P) - P/2
  
  # Store rows safely
  res_list <- list()
  rr <- 1
  
  for (bfdr_alpha in bfdr_grid) {
    for (shift_val in shift_thresholds) {
      
      prob_cons <- prob_ahead <- prob_behind <- numeric(n_circ)
      
      for (k in seq_along(circ_idx)) {
        i <- circ_idx[k]
        vec <- phase_diff_matrix[i, ]
        prob_ahead[k]  <- mean(vec >=  shift_val)
        prob_behind[k] <- mean(vec <= -shift_val)
        prob_cons[k]   <- 1 - (prob_ahead[k] + prob_behind[k])
      }
      
      bfdr_apply <- function(p) {
        ord <- order(-p)
        bfdr <- cumsum(1 - p[ord]) / seq_along(p)
        ord[bfdr <= bfdr_alpha]
      }
      
      idx_ahead  <- bfdr_apply(prob_ahead)
      idx_behind <- bfdr_apply(prob_behind)
      idx_cons   <- bfdr_apply(prob_cons)
      
      # Enforce exclusivity
      overlap <- Reduce(intersect, list(idx_ahead, idx_behind, idx_cons))
      if (length(overlap) > 0) {
        idx_ahead  <- setdiff(idx_ahead, overlap)
        idx_behind <- setdiff(idx_behind, overlap)
        idx_cons   <- setdiff(idx_cons, overlap)
      }
      
      idx_union <- union(union(idx_ahead, idx_behind), idx_cons)
      idx_undet <- setdiff(seq_len(n_circ), idx_union)
      
      res_list[[rr]] <- data.frame(
        BFDR_alpha = bfdr_alpha,
        Shift_Threshold = shift_val,
        N_Circadian = n_circ,
        N_Conserved = length(idx_cons),
        N_Shifted_Ahead = length(idx_ahead),
        N_Shifted_Behind = length(idx_behind),
        N_Undetermined = length(idx_undet)
      )
      rr <- rr + 1
    }
  }
  
  do.call(rbind, res_list)
}

options(width = 200)

phase_lung_circ <- compare_phase_shift_thresholds_circadian(
  phi_matrix1 = human_LUNG$phi,
  phi_matrix2 = baboon_LUNG$phi,
  gain_loss_status = trans_lun$gain_loss_status,
  circadian_genes = circadian_genes,
  gene_names = rownames(human_LUNG$phi)
)

#############



baboon_HIP <- list(
  rho = mcmc_data_baboon$HIP,
  phi = mcmc_phi_baboon$HIP
)

baboon_SCN <- list(
  rho = mcmc_data_baboon$SCN,
  phi = mcmc_phi_baboon$SCN
)


pA <- rowMeans(baboon_SCN$rho)
pB <- rowMeans(baboon_HIP$rho)

# Step 1: classify gain/loss first
trans_scn <- transition_classify(pA, pB, bfdr_alpha = 0.25)
circadian_genes <- c(
  "BHLHE40","BHLHE41","BMAL1","ARNTL","BTRC","CLOCK","CREB1",
  "CRY1","CRY2","CSNK1D","CSNK1E","CUL1","DBP","FBXL3","FBXW11",
  "NFIL3","NPAS2","NR1D1","NR1D2","PER1","PER2","PER3",
  "RORA","RORB","RORC"
)
gain_genes  <- names(trans_scn$gain_genes)
loss_genes  <- names(trans_scn$loss_genes)
cons_genes  <- names(trans_scn$cons_genes)
circ_gain  <- intersect(circadian_genes, gain_genes)
circ_loss  <- intersect(circadian_genes, loss_genes)
circ_cons  <- intersect(circadian_genes, cons_genes)

list(
  circadian_gain = circ_gain,
  circadian_loss = circ_loss,
  circadian_conserved = circ_cons
)






# Step 2: run the exclusive threshold comparison
phase_baboon <- compare_phase_shift_thresholds(
  phi_matrix1 = baboon_SCN$phi,
  phi_matrix2 = baboon_HIP$phi,
  gain_loss_status = trans_scn$gain_loss_status
)




# View table
phase_lung
phase_baboon

#### Pairwise: 

compare_all_cross_species_pairs <- function(
    mcmc_data_human,
    mcmc_phi_human,
    mcmc_data_baboon,
    mcmc_phi_baboon,
    bfdr_alpha = 0.25,
    shift_threshold = 2,  # hours
    bfdr_phase = 0.25
) {
  
  human_tissues <- names(mcmc_data_human)
  baboon_tissues <- names(mcmc_data_baboon)
  
  # Find matching tissue names
  common_tissues <- intersect(human_tissues, baboon_tissues)
  
  results_summary <- list()
  
  for (tissue in common_tissues) {
    
    cat("\n=== Analyzing:", tissue, "===\n")
    
    # Extract data
    human_data <- list(
      rho = mcmc_data_human[[tissue]],
      phi = mcmc_phi_human[[tissue]]
    )
    
    baboon_data <- list(
      rho = mcmc_data_baboon[[tissue]],
      phi = mcmc_phi_baboon[[tissue]]
    )
    
    # Classify gain/loss
    pA <- rowMeans(human_data$rho)
    pB <- rowMeans(baboon_data$rho)
    
    trans <- transition_classify(pA, pB, bfdr_alpha = bfdr_alpha)
    
    # Compare phase shifts
    phase_results <- compare_phase_shift_thresholds(
      phi_matrix1 = human_data$phi,
      phi_matrix2 = baboon_data$phi,
      gain_loss_status = trans$gain_loss_status,
      shift_thresholds = shift_threshold,
      bfdr_grid = bfdr_phase
    )
    
    # Store results
    if (!is.null(phase_results)) {
      results_summary[[tissue]] <- list(
        transition = trans,
        phase = phase_results,
        n_maintained = sum(trans$gain_loss_status == "Maintained"),
        n_gained = sum(trans$gain_loss_status == "Gain"),
        n_lost = sum(trans$gain_loss_status == "Loss")
      )
    }
  }
  
  return(results_summary)
}

# Run comparison
all_pairs <- compare_all_cross_species_pairs(
  mcmc_data_human,
  mcmc_phi_human,
  mcmc_data_baboon,
  mcmc_phi_baboon
)

# Create summary table
summary_df <- do.call(rbind, lapply(names(all_pairs), function(tissue) {
  res <- all_pairs[[tissue]]$phase
  data.frame(
    Tissue = tissue,
    N_Maintained = res$N_Maintained,
    N_Conserved = res$N_Conserved,
    N_Shifted = res$N_Shifted,
    N_Undetermined = res$N_Undetermined,
    Pct_Conserved = round(100 * res$N_Conserved / res$N_Maintained, 1),
    Pct_Shifted = round(100 * res$N_Shifted / res$N_Maintained, 1)
  )
}))

# Sort by different criteria
cat("\n=== Tissues with HIGHEST Phase Conservation ===\n")
print(summary_df[order(-summary_df$N_Conserved), ])

cat("\n=== Tissues with HIGHEST Phase Shift ===\n")
print(summary_df[order(-summary_df$N_Shifted), ])

cat("\n=== Tissues with HIGHEST Percentage Phase Conservation ===\n")
print(summary_df[order(-summary_df$Pct_Conserved), ])