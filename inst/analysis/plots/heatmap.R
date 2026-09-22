# plot_heatmap(): ComplexHeatmap of rhythmicity and phase status for a pathway
# in two conditions, sourced by the exploratory scripts.
library(ComplexHeatmap)
library(circlize)

plot_heatmap <- function(data1, data2,
                                    pathway_genes,
                                    pathway_name,
                                    phase_results,
                                    transition_results = NULL,
                                    group_names = c("Group1", "Group2"),
                                    save_path = NULL,
                                    n_bins = 24,
                                    versions = c("full", "rhythmic_only", "both")) {

  versions <- match.arg(versions)
  
  if(!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Install ComplexHeatmap: BiocManager::install('ComplexHeatmap')")
  }
  if(!requireNamespace("circular", quietly = TRUE)) {
    stop("Install circular: install.packages('circular')")
  }
  
  # ==========================================================================
  # 1. DATA PREPARATION
  # ==========================================================================
  
  all_genes <- rownames(data1$rho)
  overlap_genes <- intersect(all_genes, pathway_genes)
  n_genes <- length(overlap_genes)
  
  cat("\nPathway:", pathway_name, "\n")
  cat("Genes:", n_genes, "\n")
  
  rho_1 <- rowMeans(data1$rho[overlap_genes, ], na.rm = TRUE)
  rho_2 <- rowMeans(data2$rho[overlap_genes, ], na.rm = TRUE)
  
  # rhythmicity status comes from transition_classify(), which is required
  if (!is.null(transition_results)) {
    all_gene_names <- names(transition_results$gain_loss_status)
    match_idx_trans <- match(overlap_genes, all_gene_names)
    
    status_vec <- transition_results$gain_loss_status[match_idx_trans]
    status_vec[is.na(status_vec)] <- "Non-rhythmic"
    
    # NA marks a gene rhythmic in neither condition
    concordance <- rep(NA_character_, n_genes)
    concordance[status_vec == "Maintained"] <- "Conserved"
    concordance[status_vec == "Gain"] <- paste0("Gain in ", group_names[2])  # rhythmic in Group2 only
    concordance[status_vec == "Loss"] <- paste0("Loss in ", group_names[2])  # rhythmic in Group1 only
    
  } else {
    stop("transition_results is required. Please provide output from transition_classify()")
  }
  
  # phase status from the phase_infer() flags; NA when neither flag is set
  phase_gene_names <- names(phase_results$peak1)
  match_idx_phase <- match(overlap_genes, phase_gene_names)
  
  phase_status <- rep(NA_character_, n_genes)
  
  for(i in 1:n_genes) {
    idx <- match_idx_phase[i]
    
    if(!is.na(idx)) {
      if(phase_results$flag_shift[idx]) {
        phase_status[i] <- "Shifted"
      } else if(phase_results$flag_cons[idx]) {
        phase_status[i] <- "Conserved"
      }
    }
  }

  # ==========================================================================
  # DETERMINE VERSIONS TO GENERATE
  # ==========================================================================

  versions_to_run <- if (versions == "both") c("full", "rhythmic_only") else versions

  # the full gene set, restored for each version
  orig_overlap_genes <- overlap_genes
  orig_n_genes <- n_genes
  orig_rho_1 <- rho_1
  orig_rho_2 <- rho_2
  orig_concordance <- concordance
  orig_phase_status <- phase_status
  orig_match_idx_phase <- match_idx_phase

  for (current_version in versions_to_run) {

    if (current_version == "rhythmic_only") {
      # genes rhythmic in at least one condition
      rhythmic_idx <- which(!is.na(orig_concordance))
      if (length(rhythmic_idx) == 0) {
        cat("  No rhythmic genes found, skipping rhythmic_only version\n")
        next
      }
      overlap_genes <- orig_overlap_genes[rhythmic_idx]
      n_genes <- length(overlap_genes)
      rho_1 <- orig_rho_1[rhythmic_idx]
      rho_2 <- orig_rho_2[rhythmic_idx]
      concordance <- orig_concordance[rhythmic_idx]
      phase_status <- orig_phase_status[rhythmic_idx]
      match_idx_phase <- orig_match_idx_phase[rhythmic_idx]
      cat("  [rhythmic_only] Filtered to", n_genes, "rhythmic genes\n")
    } else {
      overlap_genes <- orig_overlap_genes
      n_genes <- orig_n_genes
      rho_1 <- orig_rho_1
      rho_2 <- orig_rho_2
      concordance <- orig_concordance
      phase_status <- orig_phase_status
      match_idx_phase <- orig_match_idx_phase
      cat("  [full] Using all", n_genes, "genes\n")
    }

  # ==========================================================================
  # 2. GENE ORDERING
  # ==========================================================================

  deltaPhi_for_priority <- phase_results$deltaPhi.Est[match_idx_phase]

  # Group2 peak on the -6 to 18 display range, for sorting within a class
  peak2_for_sort <- phase_results$peak2[match_idx_phase]
  peak2_for_sort[peak2_for_sort > 18] <- peak2_for_sort[peak2_for_sort > 18] - 24
  peak2_for_sort[is.na(peak2_for_sort)] <- 99  # NA last

  priority <- rep(0, n_genes)
  max_rho <- pmax(rho_1, rho_2)
  # deltaPhi is Group1 - Group2, so a negative value means Group2 peaks later
  priority[!is.na(phase_status) & phase_status == "Shifted" & !is.na(deltaPhi_for_priority) & deltaPhi_for_priority < 0] <- 5.5   # shifted later
  priority[!is.na(phase_status) & phase_status == "Shifted" & (is.na(deltaPhi_for_priority) | deltaPhi_for_priority >= 0)] <- 5.0  # shifted earlier
  priority[!is.na(phase_status) & phase_status == "Conserved"] <- 4
  priority[!is.na(concordance) & concordance == "Conserved" & is.na(phase_status)] <- 3  # maintained, phase undetermined
  priority[!is.na(concordance) & grepl("^Gain", concordance)] <- 2.5
  priority[!is.na(concordance) & grepl("^Loss", concordance)] <- 2.0
  priority[is.na(concordance)] <- 1

  # priority descending, then Group2 peak ascending
  order_idx <- order(priority, peak2_for_sort, decreasing = c(TRUE, FALSE))
  genes_ord <- overlap_genes[order_idx]
  rho_1_ord <- rho_1[order_idx]
  rho_2_ord <- rho_2[order_idx]
  phase_status_ord <- phase_status[order_idx]
  concordance_ord  <- concordance[order_idx]
  
  loss_label <- paste0("Loss in ", group_names[2])
  gain_label <- paste0("Gain in ", group_names[2])

  # every class appears in the legend, Gain before Loss
  concordance_ord <- factor(
    concordance_ord,
    levels = c("Conserved", gain_label, loss_label)
  )
  
  phi_1 <- data1$phi[genes_ord, , drop = FALSE]
  phi_2 <- data2$phi[genes_ord, , drop = FALSE]
  if(colnames(phi_1)[1] == "phi.store") phi_1 <- phi_1[, -1]
  if(colnames(phi_2)[1] == "phi.store") phi_2 <- phi_2[, -1]
  
  # where deltaPhi.Est is NA it is taken from the peaks, as Group1 - Group2
  # wrapped to [-12, 12]; anything still missing is drawn as 0
  match_idx_delta <- match(genes_ord, phase_gene_names)
  deltaPhi_ord <- phase_results$deltaPhi.Est[match_idx_delta]

  peak1_ord <- phase_results$peak1[match_idx_delta]
  peak2_ord <- phase_results$peak2[match_idx_delta]

  for(i in 1:n_genes) {
    if(is.na(deltaPhi_ord[i]) && !is.na(peak1_ord[i]) && !is.na(peak2_ord[i])) {
      delta <- peak1_ord[i] - peak2_ord[i]

      if(delta > 12) {
        delta <- delta - 24
      } else if(delta < -12) {
        delta <- delta + 24
      }

      deltaPhi_ord[i] <- delta
    }
  }
  
  deltaPhi_ord[is.na(deltaPhi_ord)] <- 0
  
  # ==========================================================================
  # 3. LEFT ANNOTATIONS (with complete legends)
  # ==========================================================================
  
  # both phase levels appear in the legend
  phase_status_factor <- factor(
    phase_status_ord,
    levels = c("Shifted", "Conserved")
  )

  loss_label <- paste0("Loss in ", group_names[2])
  gain_label <- paste0("Gain in ", group_names[2])

  conc_colors <- c("#FFA500", "#4169E1", "#9370DB")
  names(conc_colors) <- c("Conserved", loss_label, gain_label)
  
  phase_colors <- c("Shifted" = "#E63946", "Conserved" = "#06A77D")
  
  left_ha <- rowAnnotation(
    `Rhythmicity Status`  = concordance_ord,
    `Phase Status` = phase_status_factor,
    col = list(`Rhythmicity Status` = conc_colors,
               `Phase Status` = phase_colors),
    
    annotation_name_side = "bottom",
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
    na_col = "white"
  )
  
  # ==========================================================================
  # 4. MAIN HEATMAP (ρ values)
  # ==========================================================================
  
  heatmap_mat <- cbind(rho_1_ord, rho_2_ord)
  colnames(heatmap_mat) <- group_names
  rownames(heatmap_mat) <- genes_ord
  
  ht_main <- Heatmap(
    heatmap_mat,
    name = "Rhythmicity_Prob",
    col = colorRamp2(c(0, 0.5, 0.7, 0.9, 1), c("#fff5f0", "#fee0d2", "#fcae91", "#fb6a4a", "#ef3b2c")),
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    
    column_names_side = "bottom",
    column_names_centered = TRUE,
    column_names_gp = gpar(fontsize = 9, fontface = "bold"),
    column_title_side = "top",
    column_title = NULL,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    
    left_annotation = left_ha,
    width = unit(3, "cm"),
    border = TRUE,
    
    top_annotation = HeatmapAnnotation(
      space = anno_empty(border = FALSE, height = unit(5, "mm")),
      show_annotation_name = FALSE
    ),
    
    heatmap_legend_param = list(
      title = expression(Pr(rho == 1)),
      title_gp = gpar(fontsize = 10, fontface = "bold")
    )
  )
  
  # ==========================================================================
  # 5. PHASE HISTOGRAM MATRICES (with context-aware intensity)
  # ==========================================================================
  
  # posterior phase histograms on -6 to 18, each row scaled to a maximum of 1
  # where the gene is rhythmic in this condition and 0.2 where it is not
  make_histogram_matrix <- function(phi_mat, from = -6, to = 18, n_bins, 
                                    is_group1, concordance_status, group_num) {
    breaks <- seq(from, to, length.out = n_bins + 1)
    hist_mat <- matrix(0, nrow = nrow(phi_mat), ncol = n_bins)
    
    for(i in 1:nrow(phi_mat)) {
      x_vals <- as.numeric(phi_mat[i, ])
      x_vals <- ifelse(x_vals > 18, x_vals - 24, x_vals)
      hist_counts <- hist(x_vals, breaks = breaks, plot = FALSE)$counts
      
      if(max(hist_counts) > 0) {
        conc_status <- concordance_status[i]
        
        if(is.na(conc_status)) {
          intensity <- 0.2
        } else if(conc_status == "Conserved") {
          intensity <- 1.0
        } else if(group_num == 1) {
          if(grepl("Loss", conc_status)) {
            intensity <- 1.0
          } else if(grepl("Gain", conc_status)) {
            intensity <- 0.2
          } else {
            intensity <- 0.2
          }
        } else {
          if(grepl("Gain", conc_status)) {
            intensity <- 1.0
          } else if(grepl("Loss", conc_status)) {
            intensity <- 0.2
          } else {
            intensity <- 0.2
          }
        }
        
        hist_mat[i, ] <- (hist_counts / max(hist_counts)) * intensity
      }
    }
    
    # only the bins nearest -6, 0, 6, 12 and 18 are labelled
    bin_centers <- round(breaks[-1] - diff(breaks)/2, 1)
    col_names <- rep("", n_bins)

    target_labels <- c(-6, 0, 6, 12, 18)
    for(label in target_labels) {
      closest_idx <- which.min(abs(bin_centers - label))
      if (length(closest_idx) > 0) {
        col_names[closest_idx] <- as.character(label)
      }
    }
    
    colnames(hist_mat) <- col_names
    rownames(hist_mat) <- rownames(phi_mat)
    return(hist_mat)
  }
  
  hist_mat_1 <- make_histogram_matrix(phi_1, n_bins = n_bins, 
                                      concordance_status = concordance_ord, 
                                      group_num = 1)
  hist_mat_2 <- make_histogram_matrix(phi_2, n_bins = n_bins, 
                                      concordance_status = concordance_ord, 
                                      group_num = 2)
  
  ht_phase1 <- Heatmap(
    hist_mat_1,
    name = paste0("Phase_", group_names[1]),
    col = colorRamp2(c(0, 0.5, 1), c("white", "#6baed6", "#08519c")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 8),
    column_names_centered = TRUE,
    column_title_side = "bottom",
    column_title = NULL,
    show_heatmap_legend = FALSE,
    width = unit(4, "cm"),
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5)
  )
  
  ht_phase2 <- Heatmap(
    hist_mat_2,
    name = paste0("Phase_", group_names[2]),
    col = colorRamp2(c(0, 0.5, 1), c("white", "#fc9272", "#a50f15")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 8),
    column_names_centered = TRUE,
    column_title_side = "bottom",
    column_title = NULL,
    show_heatmap_legend = FALSE,
    width = unit(4, "cm"),
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5)
  )
  
  # ==========================================================================
  # 6. DELTA PHI BARPLOT (symmetric around 0, calculated for ALL genes)
  # ==========================================================================
  
  # grey for genes without a phase class; shifted later red, earlier blue
  delta_colors <- rep("#E0E0E0", n_genes)
  delta_colors[!is.na(phase_status_ord) & phase_status_ord == "Shifted" & !is.na(deltaPhi_ord) & deltaPhi_ord < 0] <- "#E63946"
  delta_colors[!is.na(phase_status_ord) & phase_status_ord == "Shifted" & (is.na(deltaPhi_ord) | deltaPhi_ord >= 0)] <- "#4361EE"
  delta_colors[!is.na(phase_status_ord) & phase_status_ord == "Conserved"] <- "#06A77D"

  # bars show Group2 - Group1, positive when Group2 peaks later
  deltaPhi_display <- -deltaPhi_ord

  # symmetric axis, at least +/- 12 h
  max_abs_delta <- max(abs(deltaPhi_display), na.rm = TRUE)
  axis_limit <- max(12, ceiling(max_abs_delta))

  delta_peak_legend <- Legend(
    title = "Delta Peak",
    labels = c(paste0(group_names[2], " later"), paste0(group_names[2], " earlier"), "Conserved", "Non-classified"),
    legend_gp = gpar(fill = c("#E63946", "#4361EE", "#06A77D", "#E0E0E0")),
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )

  right_ha <- rowAnnotation(
    `Delta Peak (hours)` = anno_barplot(
      deltaPhi_display,
      baseline = 0,
      bar_width = 1,
      gp = gpar(fill = delta_colors, col = NA),
      axis_param = list(
        at = c(-axis_limit, -axis_limit/2, 0, axis_limit/2, axis_limit),
        labels = c(paste0("-", axis_limit), paste0("-", axis_limit/2), "0",
                   as.character(axis_limit/2), as.character(axis_limit)),
        side = "bottom",
        gp = gpar(fontsize = 8)
      ),
      ylim = c(-axis_limit, axis_limit),
      width = unit(3, "cm")
    ),
    annotation_name_side = "bottom",
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
  )
  
  # ==========================================================================
  # 7. GENE NAMES
  # ==========================================================================
  
  gene_ha <- rowAnnotation(
    Genes = anno_text(
      genes_ord,
      gp = gpar(fontsize = 7),
      just = "left",
      width = max_text_width(genes_ord, gp = gpar(fontsize = 7)) + unit(2, "mm")
    ),
    show_annotation_name = FALSE
  )
  
  # ==========================================================================
  # 8. COMBINE + DRAW
  # ==========================================================================
  
  ht_list <- ht_main + ht_phase1 + ht_phase2 + right_ha + gene_ha
  
  if(!is.null(save_path)) {
    if(!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
    version_suffix <- if (current_version == "rhythmic_only") "_rhythmic_only" else ""
    filename <- paste0(save_path, "/", gsub("[^A-Za-z0-9]", "_", pathway_name), "_integrated", version_suffix, ".pdf")
    fig_height <- 4 + (n_genes * 0.15)
    fig_height <- max(6, min(fig_height, 20))
    cat("Saving:", filename, "\n")
    pdf(filename, width = 10, height = fig_height)

    draw(ht_list,
         heatmap_legend_side = "left",
         annotation_legend_side = "left",
         annotation_legend_list = list(delta_peak_legend),
         merge_legend = TRUE,
         padding = unit(c(2, 2, 2, 2), "mm"))

    title_text <- pathway_name

    grid.text(
      title_text,
      x = unit(0.5, "npc"),
      y = unit(1, "npc") - unit(3, "mm"),
      gp = gpar(fontsize = 14, fontface = "bold")
    )
    
    # Add peak time titles
    decorate_heatmap_body(paste0("Phase_", group_names[1]), {
      grid.text(
        paste0(group_names[1], " Peak Time"),
        x = unit(0.5, "npc"),
        y = unit(0, "npc") - unit(9, "mm"),
        gp = gpar(fontsize = 9, fontface = "bold")
      )
    })
    
    decorate_heatmap_body(paste0("Phase_", group_names[2]), {
      grid.text(
        paste0(group_names[2], " Peak Time"),
        x = unit(0.5, "npc"),
        y = unit(0, "npc") - unit(9, "mm"),
        gp = gpar(fontsize = 9, fontface = "bold")
      )
    })
    
    dev.off()
    cat("✓ Saved\n")
  } else {
    draw(ht_list,
         heatmap_legend_side = "left",
         annotation_legend_side = "left",
         annotation_legend_list = list(delta_peak_legend),
         merge_legend = TRUE,
         padding = unit(c(2, 2, 2, 2), "mm"))

    grid.text(
      pathway_name,
      x = unit(0.5, "npc"),
      y = unit(1, "npc") - unit(3, "mm"),
      gp = gpar(fontsize = 14, fontface = "bold")
    )
    
    # Add peak time titles
    decorate_heatmap_body(paste0("Phase_", group_names[1]), {
      grid.text(
        paste0(group_names[1], " Peak Time"),
        x = unit(0.5, "npc"),
        y = unit(0, "npc") - unit(7, "mm"),
        gp = gpar(fontsize = 9, fontface = "bold")
      )
    })
    
    decorate_heatmap_body(paste0("Phase_", group_names[2]), {
      grid.text(
        paste0(group_names[2], " Peak Time"),
        x = unit(0.5, "npc"),
        y = unit(0, "npc") - unit(7, "mm"),
        gp = gpar(fontsize = 9, fontface = "bold")
      )
    })
  }
  
  n_shifted <- sum(!is.na(phase_status_ord) & phase_status_ord == "Shifted")
  n_conserved <- sum(!is.na(phase_status_ord) & phase_status_ord == "Conserved")
  n_gain <- sum(!is.na(concordance_ord) & grepl("^Gain", concordance_ord))
  n_loss <- sum(!is.na(concordance_ord) & grepl("^Loss", concordance_ord))

  cat("\n=== SUMMARY [", toupper(current_version), "] ===\n")
  cat("Total genes:    ", n_genes, "\n")
  cat("Phase Shifted:  ", n_shifted, "\n")
  cat("Phase Conserved:", n_conserved, "\n")
  cat("Gain in", group_names[2], ":", n_gain, "\n")
  cat("Loss in", group_names[2], ":", n_loss, "\n\n")

  }  # versions

  return(invisible(list(genes = genes_ord, deltaPhi = deltaPhi_ord)))
}

