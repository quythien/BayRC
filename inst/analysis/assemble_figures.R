################################################################################
# Collect the paper's figure panels and assemble the numbered figures.
#
# The analysis scripts each drop their panels somewhere under BAYRC_OUTPUT_DIR,
# under names that describe the comparison rather than the figure. Until now
# the step from those panels to Figure_2 ... Figure_6 was done by hand, which
# is why no script in the repository produces a file called Figure_*.pdf. This
# script does that step: it copies every panel to paper/subfigures/ under a
# stable name and then merges the ones that are simple multi-panel stacks.
#
# Usage:
#   Rscript assemble_figures.R [paper.dir]
#
# paper.dir defaults to <BAYRC_OUTPUT_DIR>/../paper.
#
# Figures 3 and 4 place two panels side by side. R has no dependency-free way
# to do that with finished PDFs, so the script merges them into a two-page PDF
# and says so; the side-by-side placement is the one remaining manual step.
################################################################################

# Paths come from config.R; override any of them with the matching env var.
this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(if (is.na(this.file)) getwd() else dirname(normalizePath(this.file)),
                 "config.R"))

args <- commandArgs(trailingOnly = TRUE)
paper.dir <- if (length(args) >= 1) args[1] else
  file.path(dirname(BAYRC_OUTPUT_DIR), "paper")

sub.dir <- file.path(paper.dir, "subfigures")
fig.dir <- file.path(paper.dir, "figures")
dir.create(sub.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)

fig.out <- file.path(BAYRC_OUTPUT_DIR, "figure")

# Where each panel is written, and what it is called once collected ------------
# The left-hand name is the file the generating script writes; the right-hand
# name is what this script calls it in paper/subfigures/.

panels <- list(
  # Figure 2: genome-wide and KEGG-circadian concordance heatmaps
  list(from = file.path(BAYRC_FIGURE_DIR,
                        "Baboon_Concordance_Heatmap_0.25_ward.D2.pdf"),
       to = "F2A_baboon_genomewide_concordance.pdf",
       script = "plots/heatmap_baboon.R"),
  list(from = file.path(fig.out, "heatmap_circadian_pairs", "within_baboon",
                        "Within_Baboon_Circadian_Heatmap_ward.D2.pdf"),
       to = "F2B_baboon_circadian_concordance.pdf",
       script = "plots/heatmap_circadian_pairs.R within_baboon"),

  # Figure 3: within-species phase concordance scatters
  list(from = file.path(fig.out, "baboon_brain",
                        "Baboon_SCN_HIP_Peak_Concordance_0.25_2h_new.pdf"),
       to = "F3A_baboon_SCN_HIP_phase_concordance.pdf",
       script = "Baboon_SCN_HIP.R"),
  list(from = file.path(fig.out, "baboon_brain",
                        "Baboon_SUN_PUT_Peak_Concordance_0.25_2h_new.pdf"),
       to = "F3B_baboon_SUN_PUT_phase_concordance.pdf",
       script = "Baboon_SUN_PUT.R"),

  # Figure 4: SUN-PUT enrichment dotplots
  list(from = file.path(fig.out, "baboon_brain",
                        "SUN_PUT_shifted_GO_BP_dotplot.pdf"),
       to = "F4A_SUN_PUT_shifted_GO_BP_dotplot.pdf",
       script = "plot_enrich_SUN_PUT.R"),
  list(from = file.path(fig.out, "baboon_brain",
                        "SUN_PUT_shifted_KEGG_dotplot.pdf"),
       to = "F4B_SUN_PUT_shifted_KEGG_dotplot.pdf",
       script = "plot_enrich_SUN_PUT.R"),

  # Figure 5: SUN-PUT heatmaps
  list(from = file.path(BAYRC_OUTPUT_DIR, "heatmap_baboon_SUN_PUT",
                        "heatmap_SUN_PUT.pdf"),
       to = "F5_SUN_PUT_heatmap.pdf",
       script = "Baboon_SUN_PUT.R (heatmap section)"),

  # Figure 6: cross-species lung
  list(from = file.path(fig.out, "baboon_human_lung",
                        "Baboon_Human_LUN_Peak_Concordance_0.25_2h.pdf"),
       to = "F6_baboon_human_LUN_phase_concordance.pdf",
       script = "Baboon_Human_LUN.R")
)

# Which collected panels make up each numbered figure -------------------------

figures <- list(
  Figure_2 = c("F2A_baboon_genomewide_concordance.pdf",
               "F2B_baboon_circadian_concordance.pdf"),
  Figure_3 = c("F3A_baboon_SCN_HIP_phase_concordance.pdf",
               "F3B_baboon_SUN_PUT_phase_concordance.pdf"),
  Figure_4 = c("F4A_SUN_PUT_shifted_GO_BP_dotplot.pdf",
               "F4B_SUN_PUT_shifted_KEGG_dotplot.pdf"),
  Figure_5 = "F5_SUN_PUT_heatmap.pdf",
  Figure_6 = "F6_baboon_human_LUN_phase_concordance.pdf"
)

# Collect ---------------------------------------------------------------------

missing <- character(0)
for (p in panels) {
  dest <- file.path(sub.dir, p$to)
  if (file.exists(p$from)) {
    file.copy(p$from, dest, overwrite = TRUE)
    cat("collected", p$to, "\n")
  } else {
    missing <- c(missing, sprintf("%s  (run %s)", p$to, p$script))
  }
}

# Merge -----------------------------------------------------------------------
# qpdf and pdftools are the usual ways to concatenate PDFs from R; pdfunite
# from poppler does the same job from the shell. Whichever is present is used.

merge_pdf <- function(inputs, output) {
  if (requireNamespace("qpdf", quietly = TRUE)) {
    qpdf::pdf_combine(inputs, output); return(TRUE)
  }
  if (requireNamespace("pdftools", quietly = TRUE)) {
    pdftools::pdf_combine(inputs, output); return(TRUE)
  }
  if (nzchar(Sys.which("pdfunite"))) {
    system2("pdfunite", c(inputs, output)); return(file.exists(output))
  }
  FALSE
}

no.merger <- character(0)
for (nm in names(figures)) {
  want <- file.path(sub.dir, figures[[nm]])
  have <- want[file.exists(want)]
  if (!length(have)) next
  out <- file.path(fig.dir, paste0(nm, ".pdf"))
  if (length(have) == 1L) {
    file.copy(have, out, overwrite = TRUE)
    cat("wrote", basename(out), "\n")
  } else if (merge_pdf(have, out)) {
    cat("wrote", basename(out), "as a", length(have),
        "page PDF; panels still need side-by-side placement\n")
  } else {
    no.merger <- c(no.merger, nm)
  }
  if (length(have) < length(want))
    cat("  ", nm, "is incomplete:", length(have), "of", length(want), "panels\n")
}

cat("\nFigure 1 is a hand-drawn flowchart and has no generating script.\n")

if (length(missing)) {
  cat("\nPanels not found:\n")
  cat(paste0("  ", missing, collapse = "\n"), "\n")
}
if (length(no.merger)) {
  cat("\nNo PDF merger available (install qpdf or pdftools, or put pdfunite",
      "on the PATH). Multi-panel figures left unassembled:",
      paste(no.merger, collapse = ", "), "\n")
}
cat("\npanels:", sub.dir, "\nfigures:", fig.dir, "\n")
