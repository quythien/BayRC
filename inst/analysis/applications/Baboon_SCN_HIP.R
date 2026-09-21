# Baboon suprachiasmatic nucleus versus hippocampus: rhythmic transitions,
# phase inference, pathway enrichment and the Figure 3A panel.

library(BayRC)
library(dplyr)

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
source(file.path(analysis.dir, "plots", "theme_bayrc.R"))
source(file.path(analysis.dir, "plots", "peak_concordance.R"))
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
panel_pathways <- "KEGG Circadian rhythm"

fig.dir <- file.path(BAYRC_FIGURE_DIR, "baboon_SCN_HIP")
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)

# --replot reads the plot cache instead of the draws; a full run writes it
cache.file <- file.path(fig.dir, "plot_data.rds")
if (replot) {
  cache <- read_plot_cache(cache.file)
  scn <- cache$panelA; hip <- cache$panelB
  measured <- cache$measured; pA <- cache$pA; pB <- cache$pB
  trans <- cache$trans; status <- trans$gain_loss_status
  phase <- cache$phase; maintained <- cache$maintained
  phase_class <- cache$phase_class
} else {
  load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
  load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))

  # direction is SCN to HIP
  scn <- list(rho = mcmc_data_baboon$SCN, phi = mcmc_phi_baboon$SCN)
  hip <- list(rho = mcmc_data_baboon$HIP, phi = mcmc_phi_baboon$HIP)
  measured <- rownames(scn$rho)

  pA <- rowMeans(scn$rho)
  pB <- rowMeans(hip$rho)
  trans <- transition_classify(pA, pB, bfdr_alpha = bfdr_alpha)
  status <- trans$gain_loss_status

  phase <- phase_infer(phi_matrix1 = scn$phi, phi_matrix2 = hip$phi,
                       gain_loss_status = status, bfdr_alpha = bfdr_alpha,
                       shift = shift, P = 24, compute_hdi = TRUE)

  maintained <- names(status)[status == "Maintained"]
  phase_class <- rep("Undetermined", length(status))
  names(phase_class) <- names(status)
  phase_class[phase$flag_cons]  <- "Phase-conserved"
  phase_class[phase$flag_shift] <- "Phase-shifted"
}
delta <- ((phase$peak2 - phase$peak1 + 12) %% 24) - 12

clock_genes <- c("BHLHE40", "BHLHE41", "BMAL1", "BTRC", "CLOCK", "CREB1",
                 "CRY1", "CRY2", "CSNK1D", "CSNK1E", "CUL1", "DBP", "FBXL3",
                 "FBXW11", "NFIL3", "NPAS2", "NR1D1", "NR1D2", "PER1", "PER2",
                 "PER3", "RORA", "RORB", "RORC")

# Figure 3A
p <- peak_concordance_plot(
  peak_x = phase$peak1[maintained], peak_y = phase$peak2[maintained],
  phase_class = phase_class, label_genes = clock_genes,
  title = "SCN versus hippocampus",
  xlab = "Peak Hour - Suprachiasmatic nucleus (ZT)",
  ylab = "Peak Hour - Hippocampus (ZT)", window = shift)
bayrc_save(p + theme(legend.position = "none"), file.path(fig.dir, "Baboon_SCN_HIP_Peak_Concordance"),
           width = 5, height = 4.4)
# the legend the three Figure 3 panels share
save_plot_legend(p, file.path(fig.dir, "phase_class_legend"))

# gene sets cut to the measured genes before the enrichment sees them
kegg <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))
kegg <- lapply(kegg, function(g) replace(g, g == "ARNTL", "BMAL1"))
kegg <- lapply(kegg, intersect, measured)
kegg <- kegg[lengths(kegg) >= min_measured]

select_pathways <- function(plist, method)
  pathSelect(mcmc.merge.list = list(SCN = scn, HIP = hip), pathway.list = plist,
             dataset.names = c("SCN", "HIP"), ranking.method = method,
             score_type = "pos", qvalue.cut = stage2_q,
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
} else {
  source_run <- NULL
  union_res <- select_pathways(kegg, "union")
  union_res$q <- p.adjust(union_res$pval, "BH")
  active <- union_res$pathway[union_res$q < stage1_q]

  # pathSelect names its effect column after the ranking method, so the three
  # runs are cut to the shared columns before they are stacked
  stage2_cols <- c("pathway", "size", "pval", "Expected_N_Gain",
                   "Expected_N_Loss", "Expected_N_Conserved")
  stage2 <- do.call(rbind, lapply(c("gain", "loss", "conserved"), function(m) {
    if (!length(active)) return(NULL)
    r <- select_pathways(kegg[active], m)[, stage2_cols]
    r$q <- p.adjust(r$pval, "BH")
    r$direction <- m
    r
  }))
  if (is.null(stage2))
    stage2 <- data.frame(pathway = character(), size = integer(),
                         pval = numeric(), Expected_N_Gain = numeric(),
                         Expected_N_Loss = numeric(),
                         Expected_N_Conserved = numeric(), q = numeric(),
                         direction = character())
  sig <- stage2[stage2$q < stage2_q, ]
}

# the tables a later --replot reads back
if (!replot) {
  write.csv(union_res[order(union_res$pval), c("pathway", "size", "pval", "q")],
            file.path(fig.dir, "stage1_union.csv"), row.names = FALSE)
  write.csv(sig[order(sig$direction, sig$pval),
                c("pathway", "direction", "size", "pval", "q")],
            file.path(fig.dir, "stage2_significant.csv"), row.names = FALSE)
  write_plot_cache(cache.file, dataA = scn, dataB = hip,
                   panel_genes = unlist(kegg[panel_pathways], use.names = FALSE),
                   measured = measured, pA = pA, pB = pB, trans = trans,
                   phase = phase, maintained = maintained,
                   phase_class = phase_class)
}

for (pw in panel_pathways) {
  if (!pw %in% names(kegg)) stop("pathway not in the gene set list: ", pw)
  plot_heatmap(data1 = scn, data2 = hip, pathway_genes = kegg[[pw]],
               pathway_name = pw, phase_results = phase,
               transition_results = trans,
               group_names = c("Baboon SCN", "Baboon HIP"),
               versions = "both", save_path = fig.dir)
}

write_run_record(file.path(fig.dir, "run_record.txt"), "applications/Baboon_SCN_HIP.R",
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

cat("\n== SCN vs HIP, the numbers the paper quotes ==\n")
cat("rhythmic in SCN:", sum(pA >= trans$tau_rhythmic_A),
    " rhythmic in HIP:", sum(pB >= trans$tau_rhythmic_B), "\n")
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
if (length(conserved))
  cat(sprintf("mean |delta| among phase-conserved: %.2f h\n",
              mean(abs(delta[conserved]))))
cat("pathways tested:", length(kegg), " stage 1 active:", length(active),
    " stage 2 significant:", length(unique(sig$pathway)), "\n\n")
print(union_res[union_res$q < stage1_q, c("pathway", "size", "pval", "q")],
      row.names = FALSE)
if (nrow(sig))
  print(sig[order(sig$direction, sig$pval), c("pathway", "direction", "pval", "q")],
        row.names = FALSE)
cat("\nfigures:", fig.dir, "\n")

# Figure 3D: how the four clock genes rhythmic in SCN hold their SCN phase
# across the rest of the atlas. Drawn from the table scn_clock_reference.R
# writes, so this script does not repeat that scan.
clock_ref <- file.path(BAYRC_OUTPUT_DIR, "scn_clock_reference.csv")
if (!file.exists(clock_ref)) {
  message("no scn_clock_reference.csv under ", BAYRC_OUTPUT_DIR,
          "; skipping the clock reference panel")
} else {
  ref <- read.csv(clock_ref, stringsAsFactors = FALSE)
  # the canonical clock, in loop order. NPAS2 and RORC are absent from the
  # atlas. Genes that are not rhythmic in SCN draw an empty column, which is
  # itself the reading.
  clock_order <- c("CLOCK", "BMAL1", "PER1", "PER2", "CRY1", "CRY2",
                   "NR1D1", "NR1D2", "DBP")
  ref <- ref[ref$gene %in% clock_order, ]
  ref$dphi[ref$status != "Maintained"] <- NA

  # tissues ordered by how much of the clock they hold with SCN
  held <- tapply(ref$status == "Maintained", ref$tissue, sum)
  ref$tissue <- factor(ref$tissue, levels = names(sort(held)))
  ref$gene <- factor(ref$gene, levels = clock_order)

  # the differences run -6.4 to +3.9 h, so capping at 6 spends the whole ramp
  # on the range the data occupy and keeps one hour visibly distinct
  lim <- 6
  fig3d <- ggplot(ref, aes(x = gene, y = tissue, fill = dphi)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_gradientn(colours = bayrc_div(256), limits = c(-lim, lim),
                         oob = scales::squish,
                         na.value = "#bdbab3", name = "Δφ from SCN (h)",
                         breaks = seq(-6, 6, by = 2)) +
    labs(title = "Clock genes against SCN", x = NULL, y = NULL) +
    theme_bayrc(base_size = 10) +
    theme(axis.text.x = element_text(face = "italic"),
          panel.grid = element_blank())
  bayrc_save(fig3d, file.path(fig.dir, "SCN_clock_reference"),
             width = 9, height = 4.2)

  cat("\nclock genes rhythmic in SCN, against", nlevels(ref$tissue), "tissues\n")
  print(with(ref[ref$status == "Maintained", ],
             data.frame(gene = levels(gene),
                        maintained = as.vector(table(gene)),
                        phase_conserved = as.vector(tapply(conserved, gene, sum)),
                        median_dphi = round(as.vector(tapply(dphi, gene, median)), 2))),
        row.names = FALSE)
}
