## Redraw the Figure 2 panels from the stored concordance matrices.
##
## The matrices are written by plots/heatmap_baboon.R and
## plots/heatmap_circadian_pairs.R, so this needs neither the MCMC output nor
## the Rcpp concordance routine and is the quick path when only the drawing
## changes.
##
## Usage:
##   Rscript replot_figure2.R [outdir]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
here <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
analysis.dir <- here
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
source(file.path(analysis.dir, "plots", "palette_concordance.R"))
source(file.path(analysis.dir, "plots", "shared_legend.R"))
source(file.path(analysis.dir, "pipeline", "run_record.R"))

suppressPackageStartupMessages({library(pheatmap); library(BayRC)})

args   <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else
  file.path(BAYRC_FIGURE_DIR, "figure2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# the concordance matrices are written under the analysis output root
fig_dir <- BAYRC_OUTPUT_DIR

panels <- list(
  Fig2A_genomewide = list(
    csv   = file.path(fig_dir, "all_plots", "Baboon_Concordance_Matrix.csv"),
    title = "Genome-wide rhythmicity"),
  Fig2B_circadian = list(
    csv   = file.path(fig_dir, "heatmap_circadian_pairs", "within_baboon",
                      "pairwise_concordance_baboon_circadian_Matrix.csv"),
    title = "Circadian pathway")
)

off_max <- setNames(numeric(length(panels)), names(panels))

for (stem in names(panels)) {
  p <- panels[[stem]]
  if (!file.exists(p$csv)) {
    stop("missing concordance matrix: ", p$csv)
  }
  m <- as.matrix(read.csv(p$csv, row.names = 1, check.names = FALSE))
  diag(m) <- 1
  d <- as.dist(1 - m)

  # the assembler gives each panel half the figure width, so a wide canvas is
  # scaled down twice over and the tissue codes stop being readable in print
  cairo_pdf(file.path(outdir, paste0(stem, ".pdf")), width = 7, height = 7)
  pheatmap(m,
           cluster_rows  = hclust(d, method = "ward.D2"),
           cluster_cols  = hclust(d, method = "ward.D2"),
           color         = concordance_colors,
           breaks        = concordance_breaks,
           border_color  = NA,
           main          = p$title,
           fontsize      = 14,
           legend        = FALSE)
  dev.off()

  off <- m[row(m) != col(m)]
  off_max[stem] <- max(off)
  cat(sprintf("%s  %d tissues  off-diagonal range %.3f to %.3f, %d pairs above %.2f\n",
              stem, nrow(m), min(off), max(off), sum(off > concordance_max),
              concordance_max))
}

# Both panels are drawn on this one scale, so the colour bar is written once
# and the assembler places it under the pair.
concordance_fun <- circlize::colorRamp2(
  seq(0, concordance_max, length.out = length(concordance_colors)),
  concordance_colors)

# the circadian panel runs past the cap, so its top label declares that the
# darkest cells are clamped rather than reading as exactly 0.5
bar_labels <- format(concordance_legend, digits = 2)
if (any(off_max > concordance_max))
  bar_labels[length(bar_labels)] <- sprintf("≥ %.1f", concordance_max)

concordance_bar <- ComplexHeatmap::Legend(
  col_fun = concordance_fun, title = "Adjusted c-score",
  at = concordance_legend, labels = bar_labels,
  direction = "horizontal", legend_width = unit(7, "cm"),
  title_position = "lefttop", title_gp = gpar(fontsize = 10, fontface = "bold"),
  labels_gp = gpar(fontsize = 9))
# the top label sits past the end of the bar, so the file needs room for it
save_legend_grob(concordance_bar@grob, file.path(outdir, "Fig2_concordance_legend"),
                 pad = 0.35)

write_run_record(file.path(outdir, "run_record.txt"), "plots/replot_figure2.R",
                 c(list(colour_cap = concordance_max,
                        clustering = "ward.D2 on 1 - concordance"),
                   setNames(lapply(names(panels), function(n)
                     sprintf("%s (written %s)", panels[[n]]$csv,
                             format(file.info(panels[[n]]$csv)$mtime,
                                    "%Y-%m-%d %H:%M"))), names(panels))),
                 repo = analysis.dir)
