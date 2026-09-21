# Baboon putamen versus visual cortex: rhythmic
# transitions, phase inference, pathway enrichment and the Figure 3C, 4B and 5B
# panels. Putamen is condition A here and in Baboon_PUT_SUN.R.

library(BayRC)
library(dplyr)

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
source(file.path(analysis.dir, "plots", "theme_bayrc.R"))
source(file.path(analysis.dir, "plots", "peak_concordance.R"))
source(file.path(analysis.dir, "plots", "palette_concordance.R"))
source(file.path(analysis.dir, "plots", "shared_legend.R"))
source(file.path(analysis.dir, "pipeline", "run_record.R"))
source(file.path(analysis.dir, "pipeline", "plot_cache.R"))

# --replot redraws every figure from the tables a full run left behind
replot <- "--replot" %in% commandArgs(trailingOnly = TRUE)

# frozen analysis parameters
bfdr_alpha     <- 0.25
shift          <- 2
stage1_q       <- 0.20
stage2_q       <- 0.05
nperm          <- 10000
min_measured   <- 15
panel_pathways <- "KEGG Parkinson disease"
# Figure 4A and 4B sit side by side, so both panels fix the same colour and size
# scales; Baboon_PUT_SUN.R repeats these two values.
q_limits       <- c(0.05, 0.0005)
size_limits    <- c(0, 120)

fig.dir <- file.path(BAYRC_FIGURE_DIR, "baboon_PUT_VIC")
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)

# --replot reads the plot cache instead of the draws; a full run writes it
cache.file <- file.path(fig.dir, "plot_data.rds")
if (replot) {
  cache <- read_plot_cache(cache.file)
  put <- cache$panelA; vic <- cache$panelB
  measured <- cache$measured; pA <- cache$pA; pB <- cache$pB
  trans <- cache$trans; status <- trans$gain_loss_status
  phase <- cache$phase; maintained <- cache$maintained
  phase_class <- cache$phase_class
} else {
  load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
  load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))

  # direction is PUT to VIC, so a loss is a gene that stops being rhythmic in VIC
  put <- list(rho = mcmc_data_baboon$PUT, phi = mcmc_phi_baboon$PUT)
  vic <- list(rho = mcmc_data_baboon$VIC, phi = mcmc_phi_baboon$VIC)
  measured <- rownames(put$rho)

  pA <- rowMeans(put$rho)
  pB <- rowMeans(vic$rho)
  trans <- transition_classify(pA, pB, bfdr_alpha = bfdr_alpha)
  status <- trans$gain_loss_status

  phase <- phase_infer(phi_matrix1 = put$phi, phi_matrix2 = vic$phi,
                       gain_loss_status = status, bfdr_alpha = bfdr_alpha,
                       shift = shift, P = 24, compute_hdi = TRUE)

  maintained <- names(status)[status == "Maintained"]
  phase_class <- rep("Undetermined", length(status))
  names(phase_class) <- names(status)
  phase_class[phase$flag_cons]  <- "Phase-conserved"
  phase_class[phase$flag_shift] <- "Phase-shifted"
}

# signed difference, VIC minus PUT, on the circle
delta <- ((phase$peak2 - phase$peak1 + 12) %% 24) - 12

clock_genes <- c("BHLHE40", "BHLHE41", "BMAL1", "BTRC", "CLOCK", "CREB1",
                 "CRY1", "CRY2", "CSNK1D", "CSNK1E", "CUL1", "DBP", "FBXL3",
                 "FBXW11", "NFIL3", "NPAS2", "NR1D1", "NR1D2", "PER1", "PER2",
                 "PER3", "RORA", "RORB", "RORC")

# Figure 3C, condition A on x as in the other two panels
p <- peak_concordance_plot(
  peak_x = phase$peak1[maintained], peak_y = phase$peak2[maintained],
  phase_class = phase_class, label_genes = clock_genes,
  title = "Putamen versus visual cortex",
  xlab = "Peak Hour - Putamen (ZT)",
  ylab = "Peak Hour - Visual cortex (ZT)", window = shift)
bayrc_save(p + theme(legend.position = "none"), file.path(fig.dir, "Baboon_PUT_VIC_Peak_Concordance"),
           width = 5, height = 4.4)

# gene sets cut to the measured genes before the enrichment sees them
kegg <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))
kegg <- lapply(kegg, function(g) replace(g, g == "ARNTL", "BMAL1"))
kegg <- lapply(kegg, intersect, measured)
kegg <- kegg[lengths(kegg) >= min_measured]

select_pathways <- function(plist, method)
  pathSelect(mcmc.merge.list = list(PUT = put, VIC = vic),
             pathway.list = plist, dataset.names = c("PUT", "VIC"),
             ranking.method = method, score_type = "pos", qvalue.cut = stage2_q,
             pathwaysize.lower.cut = min_measured,
             pathwaysize.upper.cut = length(measured),
             nperm = nperm, nproc = 1)$results

# stage 1 keeps the rhythmically active pathways and stage 2 tests the
# transitions within those; --replot reads both from the tables a full run wrote
if (replot) {
  source_run <- require_run_record(fig.dir, c("stage1_union.csv", "stage2_significant.csv",
                                              "plot_data.rds"))
  union_res <- read.csv(file.path(fig.dir, "stage1_union.csv"))
  sig <- read.csv(file.path(fig.dir, "stage2_significant.csv"))
  active <- union_res$pathway[union_res$q < stage1_q]
  stage2 <- sig
  stage2_full <- NULL
} else {
  source_run <- NULL
  union_res <- select_pathways(kegg, "union")
  union_res$q <- p.adjust(union_res$pval, "BH")
  active <- union_res$pathway[union_res$q < stage1_q]

  # pathSelect names its effect column after the ranking method, so the three
  # runs are cut to the shared columns before they are stacked
  stage2_cols <- c("pathway", "size", "pval", "Expected_N_Gain",
                   "Expected_N_Loss", "Expected_N_Conserved")
  stage2_full <- lapply(c("gain", "loss", "conserved"), function(m) {
    if (!length(active)) return(NULL)
    r <- select_pathways(kegg[active], m)
    r$q <- p.adjust(r$pval, "BH")
    r$direction <- m
    r
  })
  names(stage2_full) <- c("gain", "loss", "conserved")
  stage2 <- do.call(rbind, lapply(stage2_full, function(r)
    if (is.null(r)) NULL else r[, c(stage2_cols, "q", "direction")]))
  if (is.null(stage2))
    stage2 <- data.frame(pathway = character(), size = integer(),
                         pval = numeric(), Expected_N_Gain = numeric(),
                         Expected_N_Loss = numeric(),
                         Expected_N_Conserved = numeric(), q = numeric(),
                         direction = character())
  sig <- stage2[stage2$q < stage2_q, ]
}

# Figure 4: a dot per pathway and transition that clears the stage-2 cut, sized
# by the expected gene count for that transition and shaded by -log10(q)
wrap_label <- function(x, width = 26)
  vapply(x, function(s) paste(strwrap(s, width), collapse = "\n"), character(1),
         USE.NAMES = FALSE)

if (nrow(sig)) {
# a pathway that cleared stage 1 but has no enriched transition would draw an
# empty row, so the panel keeps only the pathways with a dot
plot4 <- stage2[stage2$pathway %in% sig$pathway, ]
plot4$n_expected <- with(plot4,
  ifelse(direction == "gain", Expected_N_Gain,
  ifelse(direction == "loss", Expected_N_Loss, Expected_N_Conserved)))
plot4$direction <- factor(plot4$direction, levels = c("gain", "loss", "conserved"))
plot4$label <- wrap_label(sub(" - multiple diseases", "",
                              sub("^KEGG ", "", plot4$pathway)))

# best q across the transitions puts the strongest pathway at the top
best_q <- tapply(plot4$q, plot4$label, min)
plot4$label <- factor(plot4$label, levels = names(sort(best_q, decreasing = TRUE)))

q_breaks <- c(0.05, 0.01, 0.001)
fig4 <- ggplot(plot4[plot4$q < stage2_q, ],
               aes(x = direction, y = label, size = n_expected,
                   colour = -log10(q))) +
  geom_point() +
  scale_colour_gradientn(colours = enrichment_colors, name = "q",
                         limits = -log10(q_limits), oob = scales::squish,
                         breaks = -log10(q_breaks),
                         labels = format(q_breaks, drop0trailing = TRUE)) +
  scale_size_continuous(name = "expected gene count", range = c(2.5, 9),
                        limits = size_limits) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  labs(title = "Putamen versus visual cortex",
       x = "Transition", y = NULL) +
  theme_bayrc(base_size = 12) +
  theme(axis.text.y = element_text(size = 10))
# Figure 4A and 4B share one legend, written here and placed beneath the pair
bayrc_save(fig4 + theme(legend.position = "none"),
           file.path(fig.dir, "PUT_VIC_transition_enrichment"),
           width = 5.4, height = 3.9)
save_plot_legend(fig4, file.path(fig.dir, "transition_enrichment_legend"))
}

selected <- unique(sig$pathway)

# pathway concordance metrics behind the enrichment table, and the tables a
# later --replot reads back
if (!replot) {
  metrics <- if (!length(selected)) NULL else multi_conservation(
    mcmc.merge.list = list(PUT = put, VIC = vic),
    dataset.names = c("PUT", "VIC"), select.pathway.list = kegg[selected],
    n_perm = 1000, n_boot = 1000,
    output.dir = file.path(fig.dir, "multiconservation"), use_cpp = TRUE)

  write.csv(union_res[order(union_res$pval), c("pathway", "size", "pval", "q")],
            file.path(fig.dir, "stage1_union.csv"), row.names = FALSE)
  write.csv(sig[order(sig$direction, sig$pval),
                c("pathway", "direction", "size", "pval", "q",
                  "Expected_N_Gain", "Expected_N_Loss", "Expected_N_Conserved")],
            file.path(fig.dir, "stage2_significant.csv"), row.names = FALSE)
  if (!is.null(metrics))
    write.csv(metrics, file.path(fig.dir, "pathway_metrics.csv"), row.names = FALSE)
  write_plot_cache(cache.file, dataA = put, dataB = vic,
                   panel_genes = unlist(kegg[panel_pathways], use.names = FALSE),
                   measured = measured, pA = pA, pB = pB, trans = trans,
                   phase = phase, maintained = maintained,
                   phase_class = phase_class)
}

# Figure 5B, and the legend both Figure 5 panels share
for (pw in panel_pathways) {
  if (!pw %in% names(kegg)) stop("pathway not in the gene set list: ", pw)
  plot_heatmap(data1 = put, data2 = vic, pathway_genes = kegg[[pw]],
               pathway_name = pw, phase_results = phase,
               transition_results = trans,
               group_names = c("Putamen", "Visual cortex"),
               legend_names = c("Putamen", "the compared region"),
               versions = "both", save_path = fig.dir,
               show_title = FALSE, show_legend = FALSE,
               legend_path = file.path(fig.dir, "parkinson_heatmap_legend"),
               # the strip sits under panel A alone, so it wraps to that width
               legend_max_width = 11)
}

write_run_record(file.path(fig.dir, "run_record.txt"), "applications/Baboon_PUT_VIC.R",
                 list(bfdr_alpha = bfdr_alpha, shift = shift,
                      stage1_q = stage1_q, stage2_q = stage2_q, nperm = nperm,
                      min_measured = min_measured,
                      pathway_list = "kegg_pathway_list_hsa.rds",
                      panels = panel_pathways),
                 repo = analysis.dir,
                 replot_of = if (is.null(source_run)) NULL else source_run$run_at)

shifted <- names(phase$flag_shift)[phase$flag_shift]
conserved <- names(phase$flag_cons)[phase$flag_cons]
within <- circ_diff(phase$peak1[maintained], phase$peak2[maintained]) <= shift

cat("\n== PUT vs VIC, the numbers the paper quotes ==\n")
cat("rhythmic in PUT:", sum(pA >= trans$tau_rhythmic_A),
    " rhythmic in VIC:", sum(pB >= trans$tau_rhythmic_B), "\n")
cat("gain:", sum(status == "Gain"), " loss:", sum(status == "Loss"),
    " maintained:", length(maintained), "\n")
cat(sprintf("within +/-%g h: %d of %d (%.1f%%)\n", shift, sum(within),
            length(maintained), 100 * mean(within)))
cat("phase-shifted:", length(shifted), " phase-conserved:", length(conserved),
    " undetermined:", length(maintained) - length(shifted) - length(conserved), "\n")
# the offset over the whole maintained set, not only the genes past the window
cat(sprintf("offset over maintained: mean %+.2f h, median %+.2f, SD %.2f, %.1f%% positive\n",
            mean(delta[maintained]), median(delta[maintained]),
            sd(delta[maintained]), 100 * mean(delta[maintained] > 0)))
if (length(shifted))
  cat(sprintf("  within the shifted class alone: mean %+.2f h, SD %.2f\n",
              mean(delta[shifted]), sd(delta[shifted])))
cat("pathways tested:", length(kegg), " stage 1 active:", length(active),
    " stage 2 significant:", length(selected), "\n")
if (nrow(sig))
  cat(sprintf("Figure 4 dot range: q %.4f to %.4f, expected genes %.1f to %.1f\n",
              min(sig$q), max(sig$q), min(plot4$n_expected[plot4$q < stage2_q]),
              max(plot4$n_expected[plot4$q < stage2_q])))
cat("\n")
print(union_res[union_res$q < stage1_q, c("pathway", "size", "pval", "q")],
      row.names = FALSE)
print(sig[order(sig$direction, sig$pval), c("pathway", "direction", "pval", "q")],
      row.names = FALSE)
# the gain and loss NES say how the active pathways sit against background on
# the two transitions that carry no enrichment
for (m in if (is.null(stage2_full)) character() else c("gain", "loss")) {
  r <- stage2_full[[m]]
  if (is.null(r)) next
  cat("\n", m, "enrichment across the stage-1 active set\n")
  print(r[order(r$NES), c("pathway", "NES", "pval", "q")], row.names = FALSE)
}
for (pw in panel_pathways) {
  g <- kegg[[pw]]
  s <- status[g]
  cat("\n", pw, ": measured", length(g), "| gain", sum(s == "Gain"),
      "| loss", sum(s == "Loss"), "| maintained", sum(s == "Maintained"), "\n")
  msh <- intersect(g, shifted); mco <- intersect(g, conserved)
  cat("  phase-shifted", length(msh), "| phase-conserved", length(mco), "\n")
  if (length(msh)) print(round(sort(delta[msh]), 2))
  if (length(mco)) print(round(delta[mco], 2))
}
cat("\nfigures:", fig.dir, "\n")
