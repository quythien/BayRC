# BayRC Pathway Heatmap Visualization
# Requires: ComplexHeatmap (Bioconductor), circular

#' Integrated pathway heatmap with rhythmicity and phase information
#'
#' @title Pathway-level expression heatmap
#'
#' @description
#' Produces a \code{ComplexHeatmap} visualisation for a given pathway,
#' showing posterior rhythmicity probabilities and phase information for
#' both conditions side by side.  Genes can be coloured by gain/loss/
#' maintained status from \code{transition_classify} output.  The plot
#' can show all genes, only rhythmic genes, or both (two panels).
#'
#' @param data1 Named list; MCMC output for condition 1 with \code{rho}
#'   and \code{phi} matrices.
#' @param data2 Named list; MCMC output for condition 2.
#' @param pathway_genes Character vector; gene names in the pathway.
#' @param pathway_name Character; pathway label used in the plot title.
#' @param phase_results Output from \code{phase_infer} or a similar list
#'   with per-gene phase metrics.
#' @param transition_results Output from \code{transition_classify} or
#'   \code{NULL} (default \code{NULL}); used to colour genes by
#'   gain/loss/maintained status.
#' @param group_names Character vector of length 2; condition labels
#'   (default \code{c("Group1", "Group2")}).
#' @param save_path Character or \code{NULL}; file path for saving the
#'   plot PNG; if \code{NULL} the plot is drawn but not saved.
#' @param n_bins Integer; number of bins for the phase colour wheel
#'   (default 24).
#' @param col_main Colour mapping for the posterior-probability panel, as
#'   returned by \code{circlize::colorRamp2}; defaults to the package's own
#'   five-stop white-to-red ramp over \[0, 1\].
#' @param col_phase1 Colour mapping for the first phase panel; defaults to a
#'   three-stop white-to-blue ramp over \[0, 1\].
#' @param col_phase2 Colour mapping for the second phase panel; defaults to a
#'   three-stop white-to-red ramp over \[0, 1\].
#' @param legend_names Character vector of length 2; the condition names used
#'   in the legend text, which defaults to \code{group_names}. One legend
#'   shared by two panels needs wording that fits both, while the column and
#'   peak-time annotations keep the names in \code{group_names}.
#' @param extra_legends List of \code{ComplexHeatmap::Legend} objects packed
#'   into the shared legend file alongside the heatmap's own, so a figure whose
#'   panels are not all heatmaps still has one legend region.
#' @param legend_max_width Numeric centimetres or \code{NULL}; caps the width of
#'   the shared legend so it wraps onto further rows. Use it when the strip sits
#'   under one panel rather than the whole row, where scaling to fit would
#'   shrink the type instead.
#' @param canvas_width Numeric inches, \code{NULL} or \code{"fit"}; the page
#'   the heatmap is drawn on. The default narrows the page for a pathway with
#'   few genes, so the fixed-width blocks fill it rather than sitting in
#'   margin, and \code{"fit"} sizes the page to the drawn content.
#' @param data3 Named list or \code{NULL}; a second condition to compare
#'   against \code{data2}. When given, the posterior block gains a column, the
#'   peak-time blocks gain one for it and a second offset block is drawn.
#' @param phase_results3 Output from \code{phase_infer} for \code{data3}
#'   against \code{data2}; required when \code{data3} is given.
#' @param transition_results3 Output from \code{transition_classify} for
#'   \code{data3} against \code{data2}.
#' @param col_phase3 Colour mapping for the third phase panel.
#' @param block_width Numeric centimetres; width of each peak-time block.
#' @param delta_width Numeric centimetres; width of each offset block.
#' @param title_size Numeric; type size of the pathway title.
#' @param row_order Character vector or \code{NULL}; the genes to draw and the
#'   order to draw them in, replacing the ranking the function would otherwise
#'   apply. Genes absent from the data are dropped.
#' @param font_scale Numeric multiplier on the type drawn inside the heatmap:
#'   gene names, tick labels and the title. A figure that places the heatmap
#'   beside a scatter scales the two panels differently, so raise this until the
#'   rendered type matches its neighbour.
#' @param legend_path Character or \code{NULL}; when given, the legends are
#'   also packed horizontally and written on their own to
#'   \code{<legend_path>.pdf}, for a figure whose panels share one legend.
#' @param show_title Logical; draw \code{pathway_name} above the heatmap
#'   (default \code{TRUE}). Set \code{FALSE} for a panel whose figure carries
#'   one title over several panels.
#' @param show_legend Logical; draw the heatmap and annotation legends
#'   (default \code{TRUE}). Set \code{FALSE} for a panel that shares the
#'   legend of another panel in the same figure.
#' @param legend_side Character; which side the legends sit on, passed to
#'   \code{ComplexHeatmap::draw} (default \code{"left"}).
#' @param versions Character; which version to produce: \code{"full"},
#'   \code{"rhythmic_only"}, or \code{"both"} (default \code{"full"}).
#'
#' @return Called for side effects; invisibly returns the heatmap object.
#'
#' @export
plot_heatmap <- function(data1, data2, data3 = NULL,
                          pathway_genes,
                          pathway_name,
                          phase_results,
                          phase_results3 = NULL,
                          transition_results = NULL,
                          transition_results3 = NULL,
                          group_names = c("Group1", "Group2"),
                          save_path = NULL,
                          n_bins = 24,
                          col_main = circlize::colorRamp2(c(0, 0.5, 0.7, 0.9, 1),
                                                c("#fff5f0", "#fee0d2", "#fcae91", "#fb6a4a", "#ef3b2c")),
                          col_phase1 = circlize::colorRamp2(c(0, 0.5, 1), c("white", "#6baed6", "#08519c")),
                          col_phase2 = circlize::colorRamp2(c(0, 0.5, 1), c("white", "#fc9272", "#a50f15")),
                          col_phase3 = circlize::colorRamp2(c(0, 0.5, 1), c("white", "#a1d99b", "#006d2c")),
                          legend_names = group_names,
                          legend_path = NULL,
                          extra_legends = list(),
                          legend_max_width = NULL,
                          canvas_width = NULL,
                          title_size = 22,
                          block_width = 4,
                          delta_width = 5,
                          row_order = NULL,
                          font_scale = 1,
                          show_title = TRUE,
                          show_legend = TRUE,
                          legend_side = "left",
                          versions = c("full", "rhythmic_only", "both")) {

  versions <- match.arg(versions)
  fs <- function(size) size * font_scale
  # the region names sit under blocks of a fixed width and the title sits a
  # fixed distance above the body, so both run into their neighbours if they
  # take the full scaling
  fs_block <- function(size) size * min(font_scale, 1.1)
  
  if(!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Install ComplexHeatmap: BiocManager::install('ComplexHeatmap')")
  }
  if(!requireNamespace("circular", quietly = TRUE)) {
    stop("Install circular: install.packages('circular')")
  }
  # requireNamespace() above only confirms the packages are installed; the
  # heatmap calls below (rowAnnotation, Heatmap, colorRamp2, etc.) use bare
  # unqualified names, so these must also be attached to the search path.
  suppressPackageStartupMessages(library(ComplexHeatmap))
  if (requireNamespace("circlize", quietly = TRUE)) {
    suppressPackageStartupMessages(library(circlize))
  }
  
  # ==========================================================================
  # 1. DATA PREPARATION
  # ==========================================================================
  
  all_genes <- rownames(data1$rho)
  if (!is.null(data3)) all_genes <- intersect(all_genes, rownames(data3$rho))
  overlap_genes <- if (is.null(row_order))
    intersect(all_genes, pathway_genes) else row_order[row_order %in% all_genes]
  n_genes <- length(overlap_genes)
  
  cat("\nPathway:", pathway_name, "\n")
  cat("Genes:", n_genes, "\n")
  
  rho_1 <- rowMeans(data1$rho[overlap_genes, ], na.rm = TRUE)
  rho_2 <- rowMeans(data2$rho[overlap_genes, ], na.rm = TRUE)
  rho_3 <- if (is.null(data3)) NULL
           else rowMeans(data3$rho[overlap_genes, ], na.rm = TRUE)
  
  # Use transition_classify results if provided
  if (!is.null(transition_results)) {
    # Get gain/loss status from transition_classify
    all_gene_names <- names(transition_results$gain_loss_status)
    match_idx_trans <- match(overlap_genes, all_gene_names)
    
    status_vec <- transition_results$gain_loss_status[match_idx_trans]
    status_vec[is.na(status_vec)] <- "Non-rhythmic"
    
    # Map to concordance categories
    concordance <- rep(NA_character_, n_genes)
    concordance[status_vec == "Maintained"] <- "Conserved"
    # gain and loss are the comparator's, against the reference in data1
    concordance[status_vec == "Gain"] <-
      if (is.null(data3)) paste0("Gain in ", legend_names[2]) else "Gain"
    concordance[status_vec == "Loss"] <-
      if (is.null(data3)) paste0("Loss in ", legend_names[2]) else "Loss"
    
  } else {
    stop("transition_results is required. Please provide output from transition_classify()")
  }
  
  # Get phase results - directly use flag_shift and flag_cons
  phase_gene_names <- names(phase_results$peak1)
  match_idx_phase <- match(overlap_genes, phase_gene_names)
  
  phase_status <- rep(NA_character_, n_genes)
  
  for(i in 1:n_genes) {
    idx <- match_idx_phase[i]
    
    if(!is.na(idx)) {
      # Directly use the flags from phase_infer
      if(phase_results$flag_shift[idx]) {
        phase_status[i] <- "Shifted"
      } else if(phase_results$flag_cons[idx]) {
        phase_status[i] <- "Conserved"
      }
      # If neither flag is TRUE, phase_status remains NA (undetermined)
    }
  }

  # ==========================================================================
  # DETERMINE VERSIONS TO GENERATE
  # ==========================================================================

  versions_to_run <- if (versions == "both") c("full", "rhythmic_only") else versions

  # Store original data for looping
  orig_overlap_genes <- overlap_genes
  orig_n_genes <- n_genes
  orig_rho_1 <- rho_1
  orig_rho_2 <- rho_2
  orig_rho_3 <- rho_3
  orig_concordance <- concordance
  orig_phase_status <- phase_status
  orig_match_idx_phase <- match_idx_phase

  for (current_version in versions_to_run) {

    # Filter genes based on version
    if (current_version == "rhythmic_only") {
      # Keep only genes that are rhythmic in at least one condition (concordance is not NA)
      rhythmic_idx <- which(!is.na(orig_concordance))
      if (length(rhythmic_idx) == 0) {
        cat("  No rhythmic genes found, skipping rhythmic_only version\n")
        next
      }
      overlap_genes <- orig_overlap_genes[rhythmic_idx]
      n_genes <- length(overlap_genes)
      rho_1 <- orig_rho_1[rhythmic_idx]
      rho_2 <- orig_rho_2[rhythmic_idx]
      rho_3 <- if (is.null(orig_rho_3)) NULL else orig_rho_3[rhythmic_idx]
      concordance <- orig_concordance[rhythmic_idx]
      phase_status <- orig_phase_status[rhythmic_idx]
      match_idx_phase <- orig_match_idx_phase[rhythmic_idx]
      cat("  [rhythmic_only] Filtered to", n_genes, "rhythmic genes\n")
    } else {
      # Use all genes
      overlap_genes <- orig_overlap_genes
      n_genes <- orig_n_genes
      rho_1 <- orig_rho_1
      rho_2 <- orig_rho_2
      rho_3 <- orig_rho_3
      concordance <- orig_concordance
      phase_status <- orig_phase_status
      match_idx_phase <- orig_match_idx_phase
      cat("  [full] Using all", n_genes, "genes\n")
    }

  # ==========================================================================
  # 2. GENE ORDERING
  # ==========================================================================

  # Get deltaPhi for shifted later/earlier classification
  deltaPhi_for_priority <- phase_results$deltaPhi.Est[match_idx_phase]

  # Get peak2 (Group2 peak time) for sorting within categories
  peak2_for_sort <- phase_results$peak2[match_idx_phase]
  # Convert to -6 to 18 display range for proper visual sorting
  # Peaks > 18 should wrap to negative (e.g., 20 -> -4)
  peak2_for_sort[peak2_for_sort > 18] <- peak2_for_sort[peak2_for_sort > 18] - 24
  peak2_for_sort[is.na(peak2_for_sort)] <- 99  # Put NA at end

  priority <- rep(0, n_genes)
  max_rho <- pmax(rho_1, rho_2)
  # Split "Shifted" into later and earlier
  # NOTE: deltaPhi from phase_infer is (Group1 - Group2), so:
  #   - Negative deltaPhi = peak moved LATER in Group2 = "Shifted later"
  #   - Positive deltaPhi = peak moved EARLIER in Group2 = "Shifted earlier"
  priority[!is.na(phase_status) & phase_status == "Shifted" & !is.na(deltaPhi_for_priority) & deltaPhi_for_priority < 0] <- 5.5   # Shifted later first (negative deltaPhi)
  priority[!is.na(phase_status) & phase_status == "Shifted" & (is.na(deltaPhi_for_priority) | deltaPhi_for_priority >= 0)] <- 5.0  # Shifted earlier second (positive deltaPhi)
  priority[!is.na(phase_status) & phase_status == "Conserved"] <- 4
  priority[!is.na(concordance) & concordance == "Conserved" & is.na(phase_status)] <- 3  # Rhythmically conserved, phase undetermined
  priority[!is.na(concordance) & grepl("^Gain", concordance)] <- 2.5  # Gain grouped first
  priority[!is.na(concordance) & grepl("^Loss", concordance)] <- 2.0  # Loss grouped second
  priority[is.na(concordance)] <- 1

  # Sort by priority (descending), then by peak2 (ascending) to show time trend
  # a caller that supplied row_order has already fixed the order
  order_idx <- if (is.null(row_order))
    order(priority, peak2_for_sort, decreasing = c(TRUE, FALSE)) else
    seq_along(overlap_genes)
  genes_ord <- overlap_genes[order_idx]
  rho_1_ord <- rho_1[order_idx]
  rho_2_ord <- rho_2[order_idx]
  rho_3_ord <- if (is.null(rho_3)) NULL else rho_3[order_idx]
  phase_status_ord <- phase_status[order_idx]
  concordance_ord  <- concordance[order_idx]
  
  loss_label <- if (is.null(data3)) paste0("Loss in ", legend_names[2]) else "Loss"
  gain_label <- if (is.null(data3)) paste0("Gain in ", legend_names[2]) else "Gain"

  # Force all categories to appear in legend (Gain before Loss)
  concordance_ord <- factor(
    concordance_ord,
    levels = c("Conserved", gain_label, loss_label)
  )
  
  phi_1 <- data1$phi[genes_ord, , drop = FALSE]
  phi_2 <- data2$phi[genes_ord, , drop = FALSE]
  if(colnames(phi_1)[1] == "phi.store") phi_1 <- phi_1[, -1]
  if(colnames(phi_2)[1] == "phi.store") phi_2 <- phi_2[, -1]
  
  # Get deltaPhi from phase_results where available
  match_idx_delta <- match(genes_ord, phase_gene_names)
  deltaPhi_ord <- phase_results$deltaPhi.Est[match_idx_delta]
  
  # For genes where deltaPhi is NA, calculate it from peak1 and peak2
  peak1_ord <- phase_results$peak1[match_idx_delta]
  peak2_ord <- phase_results$peak2[match_idx_delta]
  
  for(i in 1:n_genes) {
    if(is.na(deltaPhi_ord[i]) && !is.na(peak1_ord[i]) && !is.na(peak2_ord[i])) {
      # Calculate circular difference: Group1 - Group2 (to match phase_infer convention)
      delta <- peak1_ord[i] - peak2_ord[i]

      # Wrap to [-12, 12] range for circular difference
      if(delta > 12) {
        delta <- delta - 24
      } else if(delta < -12) {
        delta <- delta + 24
      }

      deltaPhi_ord[i] <- delta
    }
  }
  
  # Set any remaining NAs to 0
  deltaPhi_ord[is.na(deltaPhi_ord)] <- 0
  
  # ==========================================================================
  # 3. LEFT ANNOTATIONS (with complete legends)
  # ==========================================================================
  
  # Define all possible levels for phase status (ensures both appear in legend)
  phase_status_factor <- factor(
    phase_status_ord,
    levels = c("Shifted", "Conserved")
  )
  
  # Define colors for all concordance states (dynamically named)
  loss_label <- if (is.null(data3)) paste0("Loss in ", legend_names[2]) else "Loss"
  gain_label <- if (is.null(data3)) paste0("Gain in ", legend_names[2]) else "Gain"

  conc_colors <- c("#FFA500", "#4169E1", "#9370DB")
  names(conc_colors) <- c("Conserved", loss_label, gain_label)
  
  phase_colors <- c("Shifted" = "#E63946", "Conserved" = "#06A77D")
  
  phase_status_factor3 <- if (is.null(phase_results3)) NULL else {
    i3 <- match(genes_ord, names(phase_results3$peak1))
    v <- rep(NA_character_, n_genes)
    v[phase_results3$flag_cons[i3]  %in% TRUE] <- "Conserved"
    v[phase_results3$flag_shift[i3] %in% TRUE] <- "Shifted"
    factor(v, levels = c("Shifted", "Conserved"))
  }
  conc_ord3 <- if (is.null(transition_results3)) NULL else {
    j3 <- match(genes_ord, names(transition_results3$gain_loss_status))
    v <- transition_results3$gain_loss_status[j3]
    out <- rep(NA_character_, n_genes)
    out[v == "Maintained"] <- "Conserved"
    out[v == "Gain"] <- gain_label
    out[v == "Loss"] <- loss_label
    factor(out, levels = c("Conserved", gain_label, loss_label))
  }
  # the strip names carry the comparator, so they are spliced in as a named list
  left_ha <- if (!is.null(phase_status_factor3)) {
    strip_names <- c(paste0("Rhythm: ", legend_names[2]),
                     paste0("Phase: ", legend_names[2]),
                     paste0("Rhythm: ", legend_names[3]),
                     paste0("Phase: ", legend_names[3]))
    strips <- setNames(list(
      concordance_ord,
      phase_status_factor,
      if (is.null(conc_ord3)) concordance_ord else conc_ord3,
      phase_status_factor3), strip_names)
    strip_cols <- setNames(list(conc_colors, phase_colors, conc_colors, phase_colors),
                           strip_names)
    strip_legends <- setNames(list(
      list(title = "Rhythmicity status",
           title_gp = gpar(fontsize = 18, fontface = "bold"),
           labels_gp = gpar(fontsize = 16)),
      list(title = "Phase status",
           title_gp = gpar(fontsize = 18, fontface = "bold"),
           labels_gp = gpar(fontsize = 16))), strip_names[1:2])
    # both comparators share one key, so only the first pair of strips shows it
    do.call(rowAnnotation, c(strips, list(
      col = strip_cols,
      annotation_legend_param = strip_legends,
      show_legend = c(TRUE, TRUE, FALSE, FALSE),
      # the two pairs repeat one colour key, so each strip names its comparison
      show_annotation_name = TRUE,
      annotation_name_side = "bottom",
      annotation_name_rot = 90,
      annotation_name_gp = gpar(fontsize = fs(14), fontface = "bold"),
      # a rotated name is as wide as its type is tall, so the strips carry it
      simple_anno_size = unit(5, "mm"),
      gap = unit(3.5, "mm"),
      na_col = "white")))
  } else rowAnnotation(
    `Rhythmicity Status`  = concordance_ord,
    `Phase Status` = phase_status_factor,
    col = list(`Rhythmicity Status` = conc_colors,
               `Phase Status` = phase_colors),
    annotation_legend_param = list(
      `Rhythmicity Status` = list(title_gp = gpar(fontsize = 18, fontface = "bold"),
                                  labels_gp = gpar(fontsize = 16)),
      `Phase Status` = list(title_gp = gpar(fontsize = 18, fontface = "bold"),
                            labels_gp = gpar(fontsize = 16))),
    
    show_annotation_name = FALSE,
    simple_anno_size = unit(3, "mm"),
    na_col = "white"
  )
  
  # ==========================================================================
  # 4. MAIN HEATMAP (rho values)
  # ==========================================================================
  
  heatmap_mat <- if (is.null(rho_3_ord)) cbind(rho_1_ord, rho_2_ord)
                 else cbind(rho_1_ord, rho_2_ord, rho_3_ord)
  colnames(heatmap_mat) <- group_names[seq_len(ncol(heatmap_mat))]
  rownames(heatmap_mat) <- genes_ord
  
  ht_main <- Heatmap(
    heatmap_mat,
    name = "Rhythmicity_Prob",
    col = col_main,
    
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = !is.null(data3),
    column_names_rot = 90,
    
    column_names_side = "bottom",
    column_names_centered = TRUE,
    column_names_gp = gpar(fontsize = fs(15), fontface = "bold"),
    column_title_side = "top",
    column_title = NULL,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    
    left_annotation = left_ha,
    width = unit(1.5 * ncol(heatmap_mat), "cm"),
    border = TRUE,
    
    top_annotation = HeatmapAnnotation(
      space = anno_empty(border = FALSE, height = unit(5, "mm")),
      show_annotation_name = FALSE
    ),
    
    heatmap_legend_param = list(
      title = expression(Pr(rho == 1)),
      title_gp = gpar(fontsize = 18, fontface = "bold"),
      labels_gp = gpar(fontsize = 16)
    )
  )
  
  # ==========================================================================
  # 5. PHASE HISTOGRAM MATRICES (with context-aware intensity)
  # ==========================================================================
  
  # every peak-time block spans this window, and its axis is drawn against it
  phase_from <- -6
  phase_to   <- 18

  # Create histogram matrix for visualization with condition-specific weighting
  make_histogram_matrix <- function(phi_mat, from = phase_from, to = phase_to, n_bins,
                                    is_group1, concordance_status, group_num) {
    breaks <- seq(from, to, length.out = n_bins + 1)
    hist_mat <- matrix(0, nrow = nrow(phi_mat), ncol = n_bins)
    
    for(i in 1:nrow(phi_mat)) {
      x_vals <- as.numeric(phi_mat[i, ])
      x_vals <- ifelse(x_vals > 18, x_vals - 24, x_vals)
      hist_counts <- hist(x_vals, breaks = breaks, plot = FALSE)$counts
      
      if(max(hist_counts) > 0) {
        conc_status <- concordance_status[i]
        
        # Determine intensity based on rhythmicity in this condition
        if(is.na(conc_status)) {
          # Non-rhythmic in both: low intensity
          intensity <- 0.2
        } else if(conc_status == "Conserved") {
          # Conserved: full intensity in both conditions
          intensity <- 1.0
        } else if(group_num == 1) {
          # Group 1 histogram
          if(grepl("Loss", conc_status)) {
            # Loss in Group2 means rhythmic in Group1: full intensity
            intensity <- 1.0
          } else if(grepl("Gain", conc_status)) {
            # Gain in Group2 means non-rhythmic in Group1: low intensity
            intensity <- 0.2
          } else {
            intensity <- 0.2
          }
        } else {
          # Group 2 histogram
          if(grepl("Gain", conc_status)) {
            # Gain in Group2 means rhythmic in Group2: full intensity
            intensity <- 1.0
          } else if(grepl("Loss", conc_status)) {
            # Loss in Group2 means non-rhythmic in Group2: low intensity
            intensity <- 0.2
          } else {
            intensity <- 0.2
          }
        }
        
        hist_mat[i, ] <- (hist_counts / max(hist_counts)) * intensity
      }
    }
    
    rownames(hist_mat) <- rownames(phi_mat)
    return(hist_mat)
  }
  
  hist_mat_1 <- make_histogram_matrix(phi_1, n_bins = n_bins, 
                                      concordance_status = concordance_ord, 
                                      group_num = 1)
  hist_mat_2 <- make_histogram_matrix(phi_2, n_bins = n_bins, 
                                      concordance_status = concordance_ord, 
                                      group_num = 2)
  if (!is.null(data3)) {
    phi_3 <- data3$phi[genes_ord, , drop = FALSE]
    if (colnames(phi_3)[1] == "phi.store") phi_3 <- phi_3[, -1]
    hist_mat_3 <- make_histogram_matrix(phi_3, n_bins = n_bins,
                                        concordance_status = if (is.null(conc_ord3))
                                          concordance_ord else conc_ord3,
                                        group_num = 2)
  }
  
  ht_phase1 <- Heatmap(
    hist_mat_1,
    name = paste0("Phase_", group_names[1]),
    col = col_phase1,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    # the peak-hour axis is drawn under the body once the layout is fixed
    show_column_names = FALSE,
    column_title = NULL,
    show_heatmap_legend = FALSE,
    width = unit(block_width, "cm"),
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5)
  )
  
  ht_phase2 <- Heatmap(
    hist_mat_2,
    name = paste0("Phase_", group_names[2]),
    col = col_phase2,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    # the peak-hour axis is drawn under the body once the layout is fixed
    show_column_names = FALSE,
    column_title = NULL,
    show_heatmap_legend = FALSE,
    width = unit(block_width, "cm"),
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5)
  )
  
  ht_phase3 <- if (is.null(data3)) NULL else Heatmap(
    hist_mat_3,
    name = paste0("Phase_", group_names[3]),
    col = col_phase3,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    # the peak-hour axis is drawn under the body once the layout is fixed
    show_column_names = FALSE,
    column_title = NULL,
    show_heatmap_legend = FALSE,
    width = unit(block_width, "cm"),
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5)
  )

  # ==========================================================================
  # 6. DELTA PHI BARPLOT (symmetric around 0, calculated for ALL genes)
  # ==========================================================================
  
  # Set gray for genes without phase status, colors for those with phase classification
  delta_colors <- rep("#E0E0E0", n_genes)  # Medium gray for non-classified
  # Split shifted into later (red) and earlier (blue) based on deltaPhi
  # NOTE: deltaPhi from phase_infer is (Group1 - Group2)
  # We DISPLAY as (Group2 - Group1), so:
  #   - Negative internal deltaPhi -> Positive display -> Group2 peaks LATER = "Shifted later" (red)
  #   - Positive internal deltaPhi -> Negative display -> Group2 peaks EARLIER = "Shifted earlier" (blue)
  delta_colors[!is.na(phase_status_ord) & phase_status_ord == "Shifted" & !is.na(deltaPhi_ord) & deltaPhi_ord < 0] <- "#E63946"   # Shifted later (red)
  delta_colors[!is.na(phase_status_ord) & phase_status_ord == "Shifted" & (is.na(deltaPhi_ord) | deltaPhi_ord >= 0)] <- "#4361EE"  # Shifted earlier (blue)
  delta_colors[!is.na(phase_status_ord) & phase_status_ord == "Conserved"] <- "#06A77D"

  # Flip deltaPhi for display: show as (Group2 - Group1)
  # Positive = Group2 is later, Negative = Group2 is earlier
  deltaPhi_display <- -deltaPhi_ord

  # Determine symmetric axis limits
  max_abs_delta <- max(abs(deltaPhi_display), na.rm = TRUE)
  axis_limit <- max(12, ceiling(max_abs_delta))

  if (!is.null(phase_results3)) {
    m3 <- match(genes_ord, names(phase_results3$peak1))
    delta3_ord <- phase_results3$deltaPhi.Est[m3]

    # the offset comes from the peaks themselves where the comparison left none,
    # as it does for the first comparator
    peak1_3 <- phase_results3$peak1[m3]
    peak2_3 <- phase_results3$peak2[m3]
    for (i in 1:n_genes) {
      if (is.na(delta3_ord[i]) && !is.na(peak1_3[i]) && !is.na(peak2_3[i])) {
        delta <- peak1_3[i] - peak2_3[i]
        if (delta > 12) {
          delta <- delta - 24
        } else if (delta < -12) {
          delta <- delta + 24
        }
        delta3_ord[i] <- delta
      }
    }
    delta3_ord[is.na(delta3_ord)] <- 0
    status3 <- rep(NA_character_, n_genes)
    status3[phase_results3$flag_cons[m3]  %in% TRUE] <- "Conserved"
    status3[phase_results3$flag_shift[m3] %in% TRUE] <- "Shifted"
    delta3_colors <- rep("#E0E0E0", n_genes)
    delta3_colors[!is.na(status3) & status3 == "Shifted" &
                  !is.na(delta3_ord) & delta3_ord < 0]  <- "#E63946"
    delta3_colors[!is.na(status3) & status3 == "Shifted" &
                  (is.na(delta3_ord) | delta3_ord >= 0)] <- "#4361EE"
    delta3_colors[!is.na(status3) & status3 == "Conserved"] <- "#06A77D"
    delta3_display <- -delta3_ord
    axis_limit <- max(axis_limit, ceiling(max(abs(delta3_display), na.rm = TRUE)))
  }

  # one key serves both offset blocks, so with two comparators it names the reference
  delta_labels <- if (is.null(data3))
    c(paste0(legend_names[2], " later"), paste0(legend_names[2], " earlier"))
  else c(paste0("Later than ", legend_names[1]), paste0("Earlier than ", legend_names[1]))
  delta_labels <- c(delta_labels, "Within the window", "Not classified")

  # Create legend for Delta Peak colors
  delta_peak_legend <- Legend(
    title = "Delta Peak",
    labels = delta_labels,
    legend_gp = gpar(fill = c("#E63946", "#4361EE", "#06A77D", "#E0E0E0")),
    title_gp = gpar(fontsize = 18, fontface = "bold"),
    labels_gp = gpar(fontsize = 16)
  )

  right_ha <- rowAnnotation(
    `Delta Peak (hours)` = anno_barplot(
      deltaPhi_display,
      baseline = 0,
      bar_width = 1,
      gp = gpar(fill = delta_colors, col = NA),
      # the offset axis is drawn under the body along with the peak-hour axes
      axis = FALSE,
      ylim = c(-axis_limit, axis_limit),
      width = unit(delta_width, "cm")
    ),
    show_annotation_name = FALSE
  )
  
  right_ha3 <- if (is.null(phase_results3)) NULL else rowAnnotation(
    `Delta Peak 2 (hours)` = anno_barplot(
      delta3_display,
      baseline = 0,
      bar_width = 1,
      gp = gpar(fill = delta3_colors, col = NA),
      axis = FALSE,
      ylim = c(-axis_limit, axis_limit),
      width = unit(delta_width, "cm")
    ),
    show_annotation_name = FALSE
  )

  # ==========================================================================
  # 7. GENE NAMES
  # ==========================================================================
  
  gene_ha <- rowAnnotation(
    Genes = anno_text(
      genes_ord,
      gp = gpar(fontsize = fs(15)),
      just = "left",
      width = max_text_width(genes_ord, gp = gpar(fontsize = fs(15))) + unit(2, "mm")
    ),
    show_annotation_name = FALSE
  )
  
  # ==========================================================================
  # 8. COMBINE + DRAW
  # ==========================================================================
  
  ht_list <- ht_main + ht_phase1 + ht_phase2
  if (!is.null(ht_phase3)) ht_list <- ht_list + ht_phase3
  ht_list <- ht_list + right_ha
  if (!is.null(right_ha3)) ht_list <- ht_list + right_ha3
  ht_list <- ht_list + gene_ha

  # ==========================================================================
  # 9. TYPE AROUND THE BODY: axes, block names, title and the room they take
  # ==========================================================================

  pt_mm <- 25.4 / 72
  axis_fs <- fs(15)
  name_gp <- gpar(fontsize = fs_block(15), fontface = "bold")
  # the axis numerals hang just under the body and the block names hang a line
  # under them, so a two-line name grows downward and never meets the numerals
  axis_y <- unit(0, "npc") - unit(1.5, "mm")
  name_y <- axis_y - unit(axis_fs * 1.2, "bigpts") - unit(2.5, "mm")

  offset_names <- if (is.null(data3)) "Delta peak (h)" else
    sprintf("%s - %s (h)", legend_names[2], legend_names[1])
  if (!is.null(phase_results3))
    offset_names <- c(offset_names, sprintf("%s - %s (h)", legend_names[3], legend_names[1]))

  pdf(NULL)
  mm_wide <- function(x, gp) convertWidth(max_text_width(x, gp = gp), "mm", valueOnly = TRUE)
  # a block name wider than its block spills into the gaps either side of it
  spill <- max(mm_wide(group_names, name_gp) - block_width * 10,
               mm_wide(offset_names, name_gp) - delta_width * 10)
  # a three-condition panel names its columns under the body in rotated type,
  # and that room already holds the axes and block names
  reserved <- if (is.null(data3)) 0 else
    mm_wide(group_names, gpar(fontsize = fs(15), fontface = "bold"))
  dev.off()
  block_gap <- max(4, spill + 3)
  ht_gaps <- if (is.null(data3)) unit(c(1.2, block_gap, block_gap, 1.2), "mm") else
    unit(c(4, block_gap, block_gap, block_gap, if (!is.null(right_ha3)) block_gap, 4), "mm")

  n_name_lines <- max(lengths(strsplit(c(group_names, offset_names), "\n", fixed = TRUE)))
  hang_mm <- 7 + (axis_fs + n_name_lines * fs_block(15)) * 1.2 * pt_mm
  below_mm <- max(4, hang_mm - reserved)
  # the title is centred 1.6 lines below the page edge and the body begins 2.4
  # lines down, so the title keeps its headroom when the page is scaled
  title_mm <- if (show_title) fs_block(title_size) * 2.4 * pt_mm else 2

  # every peak-time block spans phase_from to phase_to, so its labels sit at
  # fixed fractions of the block; each end label is justified inward and stays
  # inside the block at any scale
  draw_axis <- function(ticks, labels, at) {
    grid.segments(x0 = unit(ticks, "npc"), x1 = unit(ticks, "npc"),
                  y0 = unit(0, "npc"), y1 = unit(0, "npc") - unit(1, "mm"))
    grid.text(labels, x = unit(at, "npc"), y = axis_y, hjust = at, vjust = 1,
              gp = gpar(fontsize = axis_fs))
  }
  phase_at  <- function(h) (h - phase_from) / (phase_to - phase_from)
  offset_at <- function(h) (h + axis_limit) / (2 * axis_limit)

  decorate_all <- function() {
    if (show_title) {
      upViewport(0)
      grid.text(pathway_name, x = unit(0.5, "npc"),
                y = unit(1, "npc") - unit(fs_block(title_size) * 1.6, "bigpts"),
                gp = gpar(fontsize = fs_block(title_size), fontface = "bold"))
    }
    for (i in seq_len(if (is.null(data3)) 2 else 3))
      decorate_heatmap_body(paste0("Phase_", group_names[i]), {
        draw_axis(phase_at(seq(phase_from, phase_to, by = 6)),
                  c(phase_from, 6, phase_to), phase_at(c(phase_from, 6, phase_to)))
        grid.text(group_names[i], x = unit(0.5, "npc"), y = name_y, vjust = 1,
                  gp = name_gp)
      })
    for (j in seq_along(offset_names))
      decorate_annotation(c("Delta Peak (hours)", "Delta Peak 2 (hours)")[j], {
        draw_axis(offset_at(seq(-axis_limit, axis_limit, by = 6)),
                  c(-axis_limit, 0, axis_limit), offset_at(c(-axis_limit, 0, axis_limit)))
        grid.text(offset_names[j], x = unit(0.5, "npc"), y = name_y, vjust = 1,
                  gp = name_gp)
      })
  }

  # A figure whose panels share one legend draws them with show_legend = FALSE
  # and places this file beneath the pair. The strip is placed at its own width
  # rather than a panel's, so its type is set smaller than the in-panel text.
  if (!is.null(legend_path) && current_version == versions_to_run[1]) {
    base_legends <- list(
      Legend(title = "Rhythmicity Status", labels = names(conc_colors),
             legend_gp = gpar(fill = conc_colors),
             title_gp = gpar(fontsize = 10, fontface = "bold"),
             labels_gp = gpar(fontsize = 8)),
      Legend(title = "Phase Status", labels = names(phase_colors),
             legend_gp = gpar(fill = phase_colors),
             title_gp = gpar(fontsize = 10, fontface = "bold"),
             labels_gp = gpar(fontsize = 8)),
      Legend(title = expression(Pr(rho == 1)), col_fun = col_main, at = c(0, 0.5, 1),
             title_gp = gpar(fontsize = 10, fontface = "bold"),
             labels_gp = gpar(fontsize = 8), direction = "horizontal"),
      Legend(title = "Delta Peak",
             labels = delta_labels,
             legend_gp = gpar(fill = c("#E63946", "#4361EE", "#06A77D", "#E0E0E0")),
             title_gp = gpar(fontsize = 10, fontface = "bold"),
             labels_gp = gpar(fontsize = 8)))
    pack_args <- list(direction = "horizontal", gap = unit(6, "mm"))
    # a strip placed under one panel rather than the whole row has to wrap, or
    # it is scaled down to fit and its type shrinks with it
    if (!is.null(legend_max_width))
      pack_args$max_width <- unit(legend_max_width, "cm")
    shared <- do.call(packLegend, c(base_legends, extra_legends, pack_args))
    pdf(NULL)
    lw <- convertWidth(grobWidth(shared@grob), "in", valueOnly = TRUE)
    lh <- convertHeight(grobHeight(shared@grob), "in", valueOnly = TRUE)
    dev.off()
    pdf(paste0(legend_path, ".pdf"), width = lw + 0.2, height = lh + 0.2)
    grid.newpage()
    draw(shared)
    dev.off()
    cat("Saving:", paste0(legend_path, ".pdf"), "\n")
  }

  draw_list <- function() draw(
         ht_list,
         heatmap_legend_side = legend_side,
         annotation_legend_side = legend_side,
         annotation_legend_list = list(delta_peak_legend),
         show_heatmap_legend = show_legend,
         show_annotation_legend = show_legend,
         merge_legend = TRUE,
         ht_gap = ht_gaps,
         # the type under the body hangs into the bottom margin and the title
         # sits in the top one
         padding = unit(c(below_mm, 2, title_mm, 2), "mm"))

  if(!is.null(save_path)) {
    if(!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
    # Add version suffix to filename
    version_suffix <- if (current_version == "rhythmic_only") "_rhythmic_only" else ""
    filename <- paste0(save_path, "/", gsub("[^A-Za-z0-9]", "_", pathway_name), "_integrated", version_suffix, ".pdf")
    fig_height <- 4 + (n_genes * 0.24)
    fig_height <- max(6, min(fig_height, 26))
    # legends laid out in a row below the heatmap need their own band
    if (show_legend && legend_side == "bottom") fig_height <- fig_height + 1.2
    # the margins that hold the title and the type under the body are added to
    # the page, so the body keeps its height
    fig_height <- fig_height + (below_mm + title_mm) / 25.4
    cat("Saving:", filename, "\n")
    # the blocks are a fixed 16 cm, so a short heatmap on a 10 in canvas is
    # mostly margin and renders small beside a taller panel. Narrowing the
    # canvas for a few-gene pathway lets the blocks fill it and brings the
    # panel's aspect closer to square.
    fig_width <- if (identical(canvas_width, "fit")) {
      # every component has an absolute width, so a draw on a null device
      # measures the page that holds them with no margin to spare
      pdf(NULL, width = 30, height = fig_height)
      drawn <- draw_list()
      w <- convertWidth(drawn@ht_list_param$width, "in", valueOnly = TRUE)
      dev.off()
      w + 0.1
    } else if (is.null(canvas_width))
      max(7.5, min(10, 4.5 + fig_height * 0.45)) else canvas_width
    pdf(filename, width = fig_width, height = fig_height)
    draw_list()
    decorate_all()
    dev.off()
    cat("Saved\n")
  } else {
    draw_list()
    decorate_all()
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

  }  # End of for loop over versions

  return(invisible(list(genes = genes_ord, deltaPhi = deltaPhi_ord)))
}
