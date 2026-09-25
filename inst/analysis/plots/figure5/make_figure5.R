# Figure 5: the circadian pathway in lung across three species.
#
# Human is the reference group; baboon and mouse are comparator groups. Both
# comparisons use the frozen parameters of the case study the paper reports.

library(BayRC)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Paths come from config.R; override any of them with the matching env var.
this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
# the shared look and the scatter live with the case-study scripts
source(file.path(analysis.dir, "plots", "theme_bayrc.R"))
source(file.path(analysis.dir, "plots", "peak_concordance.R"))

# the panels go beside the other figure panels
here <- Sys.getenv("FIG5_DIR", unset = file.path(BAYRC_FIGURE_DIR, "figure5"))
dir.create(here, recursive = TRUE, showWarnings = FALSE)
res  <- BAYRC_RESULT_DIR
# the mouse summaries sit beside the human and baboon ones, under summary/hm
if (!file.exists(file.path(res, "summary", "hm", "mcmc_rho_BF3.RData")))
  stop("no mouse summary under ", file.path(res, "summary", "hm"),
       "; run mouse/build_mouse_summaries.R first")

bfdr_alpha <- 0.25
shift      <- 2
tissue     <- "LUN"

grab <- function(dir, sp) {
  r <- new.env(); load(file.path(res, "summary", dir, "mcmc_rho_BF3.RData"), envir = r)
  p <- new.env(); load(file.path(res, "summary", dir, "phi", "mcmc_phi_BF3.RData"), envir = p)
  list(rho = get(paste0("mcmc_data_", sp), envir = r)[[tissue]],
       phi = get(paste0("mcmc_phi_", sp),  envir = p)[[tissue]])
}
bab <- grab("hb", "baboon"); hum_b <- grab("hb", "human")
mou <- grab("hm", "mouse");  hum_m <- grab("hm", "human")

# both comparisons are restricted to the genes measured in both
universe <- intersect(rownames(bab$rho), rownames(mou$rho))
cat("genes measured in both comparisons:", length(universe), "\n")
sub <- function(x) lapply(x, function(m) m[universe, , drop = FALSE])
bab <- sub(bab); mou <- sub(mou); hum <- sub(hum_b); hum_m <- sub(hum_m)

kegg <- readRDS(system.file("extdata", "kegg_pathway_list_hsa.rds", package = "BayRC"))
kegg <- lapply(kegg, function(g) replace(g, g == "ARNTL", "BMAL1"))
genes <- intersect(kegg[["KEGG Circadian rhythm"]], universe)
cat("pathway genes drawn:", length(genes), "\n")

# the reference goes first, so gain and loss are the comparator's
pair <- function(ref, cmp) {
  tr <- transition_classify(rowMeans(ref$rho), rowMeans(cmp$rho),
                            bfdr_alpha = bfdr_alpha)
  ph <- phase_infer(phi_matrix1 = ref$phi, phi_matrix2 = cmp$phi,
                    gain_loss_status = tr$gain_loss_status,
                    bfdr_alpha = bfdr_alpha, shift = shift, P = 24,
                    compute_hdi = TRUE)
  list(trans = tr, phase = ph)
}
hb <- pair(hum, bab)
hm <- pair(hum_m, mou)

# rows run by which species carry the rhythm, all three first and none last
called <- function(x) {
  p <- rowMeans(x$rho)
  bfdr_from_posterior(p, alpha = bfdr_alpha)$rhythmic_genes[match(genes, names(p))]
}
r1 <- called(hum); r2 <- called(bab); r3 <- called(mou)
block <- ifelse( r1 &  r2 &  r3, 1,
         ifelse( r1 &  r2 & !r3, 2,
         ifelse( r1 & !r2 &  r3, 3,
         ifelse(!r1 &  r2 &  r3, 4,
         ifelse( r1 & !r2 & !r3, 5,
         ifelse(!r1 &  r2 & !r3, 6,
         ifelse(!r1 & !r2 &  r3, 7, 8)))))))
tie <- rowMeans(hum$rho[genes, , drop = FALSE])
row_order <- genes[order(block, -tie)]
cat("genes per block:", paste(table(factor(block, levels = 1:8)), collapse = " "), "\n")

# Panels A and B: the reference on x in both, so the two read the same way
clock <- c("BMAL1", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2",
           "NR1D1", "NR1D2", "RORA", "DBP", "NFIL3")
scatter <- function(res, cmp_name, file, show_key) {
  st <- res$trans$gain_loss_status
  keep <- names(st)[st == "Maintained"]
  cls <- rep("Undetermined", length(st)); names(cls) <- names(st)
  lev <- c("Phase-conserved", "Phase-shifted", "Undetermined")
  cls[res$phase$flag_cons]  <- "Phase-conserved"
  cls[res$phase$flag_shift] <- "Phase-shifted"
  cls <- factor(cls, levels = lev)
  p <- peak_concordance_plot(
    peak_x = res$phase$peak1[keep], peak_y = res$phase$peak2[keep],
    phase_class = cls, label_genes = clock,
    title = sprintf("Human versus %s lung", tolower(cmp_name)),
    xlab = "Peak Hour - Human lung (ZT)",
    ylab = sprintf("Peak Hour - %s lung (ZT)", cmp_name),
    window = shift)
  # the pair shares one key, carried by the first panel
  p <- p + if (show_key) {
    theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
          legend.title = element_text(size = 16), legend.text = element_text(size = 15))
  } else {
    theme(legend.position = "none")
  }
  # canvas sized so type prints at 5 pt on a text-width page
  cairo_pdf(file.path(here, file), width = 5.1, height = 4.8)
  print(p)
  dev.off()
  cat(sprintf("  %s: %d maintained; clock genes present: %s\n", cmp_name,
              length(keep), paste(intersect(clock, keep), collapse = " ")))
}
scatter(hb, "Baboon", "Fig5A_human_baboon.pdf", TRUE)
scatter(hm, "Mouse",  "Fig5B_human_mouse.pdf", FALSE)

plot_heatmap(
  data1 = hum, data2 = bab, data3 = mou,
  pathway_genes = genes, pathway_name = "KEGG Circadian rhythm",
  phase_results = hb$phase, phase_results3 = hm$phase,
  transition_results = hb$trans, transition_results3 = hm$trans,
  group_names = c("Human lung\n(Reference)", "Baboon lung", "Mouse lung"),
  legend_names = c("Human", "Baboon", "Mouse"),
  row_order = row_order, canvas_width = "fit",
  # blocks as wide as their titles, key beneath; the panel stays near 1.2:1
  block_width = 3.2, delta_width = 5.6, font_scale = 1.15,
  show_legend = TRUE, legend_side = "bottom",
  versions = "full", save_path = file.path(here, "lung3"))
# assemble_figure5.py reads the three panels under one naming scheme
file.copy(file.path(here, "lung3", "KEGG_Circadian_rhythm_integrated.pdf"),
          file.path(here, "Fig5C_heatmap.pdf"), overwrite = TRUE)
cat("done\n")
