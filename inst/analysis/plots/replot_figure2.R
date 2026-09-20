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
source(file.path(here, "palette_concordance.R"))

suppressPackageStartupMessages(library(pheatmap))

args   <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else
  file.path(Sys.getenv("BAYRC_FIGURE_DIR", unset = "."), "fig2_replot")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

fig_dir <- Sys.getenv("BAYRC_FIGURE_DIR", unset = ".")

panels <- list(
  Fig2A_genomewide = list(
    csv   = file.path(fig_dir, "all_plots", "Baboon_Concordance_Matrix.csv"),
    title = "Baboon Genome-wide Rhythmicity Concordance Heatmap"),
  Fig2B_circadian = list(
    csv   = file.path(fig_dir, "heatmap_circadian_pairs", "within_baboon",
                      "pairwise_concordance_baboon_circadian_Matrix.csv"),
    title = "Baboon Circadian Pathway Concordance Heatmap")
)

for (stem in names(panels)) {
  p <- panels[[stem]]
  if (!file.exists(p$csv)) {
    stop("missing concordance matrix: ", p$csv)
  }
  m <- as.matrix(read.csv(p$csv, row.names = 1, check.names = FALSE))
  diag(m) <- 1
  d <- as.dist(1 - m)

  cairo_pdf(file.path(outdir, paste0(stem, ".pdf")), width = 9, height = 8)
  pheatmap(m,
           cluster_rows  = hclust(d, method = "ward.D2"),
           cluster_cols  = hclust(d, method = "ward.D2"),
           color         = concordance_colors,
           breaks        = concordance_breaks,
           border_color  = NA,
           main          = p$title,
           legend_breaks = concordance_legend,
           legend_labels = format(concordance_legend, digits = 2))
  dev.off()

  off <- m[row(m) != col(m)]
  cat(sprintf("%s  off-diagonal range %.3f to %.3f\n",
              stem, min(off), max(off)))
}
