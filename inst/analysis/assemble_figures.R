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
dry.run <- "--dry-run" %in% args
args <- setdiff(args, "--dry-run")
paper.dir <- if (length(args) >= 1) args[1] else
  file.path(dirname(BAYRC_OUTPUT_DIR), "paper")

sub.dir <- file.path(paper.dir, "subfigures")
fig.dir <- file.path(paper.dir, "figures")
if (!dry.run) {
  dir.create(sub.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
}

fig.out <- file.path(BAYRC_OUTPUT_DIR, "figure")

# Where each panel is written, and what it is called once collected ------------
# The left-hand name is the file the generating script writes; the right-hand
# name is what this script calls it in paper/subfigures/.

# Panels are resolved by directory and pattern rather than by an exact file
# name, so a panel still resolves when the generating script changes the cap in
# a file name, or when a KEGG release spells a pathway differently.

find_panel <- function(dir, pattern, recursive = FALSE) {
  if (!dir.exists(dir)) return(NA_character_)
  hits <- list.files(dir, pattern = pattern, full.names = TRUE,
                     recursive = recursive)
  if (!length(hits)) return(NA_character_)
  ## a stale panel from an earlier run can match too, so take the newest
  hits <- hits[order(file.info(hits)$mtime, decreasing = TRUE)]
  if (length(hits) > 1) {
    cat("  ", length(hits), "files match", pattern, "under", dir, "\n")
    cat("     taking the newest:", basename(hits[1]), "\n")
    cat("     also matched:", paste(basename(hits[-1]), collapse = ", "), "\n")
  }
  hits[1]
}

brain.dir <- file.path(fig.out, "baboon_brain")
lung.dir  <- file.path(fig.out, "baboon_human_lung")
sunput.heatmaps <- file.path(BAYRC_OUTPUT_DIR, "heatmap_baboon_SUN_PUT")
circadian.dir <- file.path(BAYRC_OUTPUT_DIR, "heatmap_circadian_pairs",
                           "within_baboon")

panels <- list(
  # Figure 2: genome-wide and KEGG-circadian concordance heatmaps
  list(dir = BAYRC_FIGURE_DIR,
       pattern = "^Baboon_Concordance_Heatmap_0[.]5_ward[.]D2[.]pdf$",
       to = "F2A_baboon_genomewide_concordance.pdf",
       script = "plots/heatmap_baboon.R"),
  list(dir = circadian.dir,
       pattern = "_Heatmap_ward[.]D2[.]pdf$",
       to = "F2B_baboon_circadian_concordance.pdf",
       script = "plots/heatmap_circadian_pairs.R within_baboon"),

  # Figure 3: within-species phase concordance scatters
  list(dir = brain.dir,
       pattern = "^Baboon_SCN_HIP_Peak_Concordance_.*_2h_new[.]pdf$",
       to = "F3A_baboon_SCN_HIP_phase_concordance.pdf",
       script = "Baboon_SCN_HIP.R"),
  list(dir = brain.dir,
       pattern = "^Baboon_SUN_PUT_Peak_Concordance_.*_2h_new[.]pdf$",
       to = "F3B_baboon_SUN_PUT_phase_concordance.pdf",
       script = "Baboon_SUN_PUT.R"),

  # Figure 4: SUN-PUT enrichment dotplots
  list(dir = brain.dir,
       pattern = "^SUN_PUT_shifted_GO_BP_dotplot[.]pdf$",
       to = "F4A_SUN_PUT_shifted_GO_BP_dotplot.pdf",
       script = "plot_enrich_SUN_PUT.R"),
  list(dir = brain.dir,
       pattern = "^SUN_PUT_shifted_KEGG_dotplot[.]pdf$",
       to = "F4B_SUN_PUT_shifted_KEGG_dotplot.pdf",
       script = "plot_enrich_SUN_PUT.R"),

  # Figure 5: SUN-PUT pathway heatmaps. plot_heatmap() writes one file per
  # pathway into its own sub-directory, and also a rhythmic-only version;
  # the patterns take the full version of the two pathways the figure shows.
  list(dir = sunput.heatmaps, recursive = TRUE,
       pattern = "Parkinson.*_integrated[.]pdf$",
       to = "F5A_SUN_PUT_KEGG_Parkinson_heatmap.pdf",
       script = "Baboon_SUN_PUT.R (heatmap section)"),
  list(dir = sunput.heatmaps, recursive = TRUE,
       pattern = "Oxidative_phosphorylation_integrated[.]pdf$",
       to = "F5B_SUN_PUT_KEGG_OxPhos_heatmap.pdf",
       script = "Baboon_SUN_PUT.R (heatmap section)"),

  # Figure 6: cross-species lung
  list(dir = lung.dir,
       pattern = "^Baboon_Human_LUN_Peak_Concordance_.*_2h[.]pdf$",
       to = "F6_baboon_human_LUN_phase_concordance.pdf",
       script = "Baboon_Human_LUN.R")
)

for (i in seq_along(panels))
  panels[[i]]$from <- find_panel(panels[[i]]$dir, panels[[i]]$pattern,
                                 isTRUE(panels[[i]]$recursive))

# Which collected panels make up each numbered figure -------------------------

figures <- list(
  Figure_2 = c("F2A_baboon_genomewide_concordance.pdf",
               "F2B_baboon_circadian_concordance.pdf"),
  Figure_3 = c("F3A_baboon_SCN_HIP_phase_concordance.pdf",
               "F3B_baboon_SUN_PUT_phase_concordance.pdf"),
  Figure_4 = c("F4A_SUN_PUT_shifted_GO_BP_dotplot.pdf",
               "F4B_SUN_PUT_shifted_KEGG_dotplot.pdf"),
  Figure_5 = c("F5A_SUN_PUT_KEGG_Parkinson_heatmap.pdf",
               "F5B_SUN_PUT_KEGG_OxPhos_heatmap.pdf"),
  Figure_6 = "F6_baboon_human_LUN_phase_concordance.pdf"
)

# Collect ---------------------------------------------------------------------

if (dry.run) {
  cat("resolved panels\n\n")
  for (nm in names(figures)) {
    cat(nm, "\n")
    for (want in figures[[nm]]) {
      p <- Filter(function(x) x$to == want, panels)[[1]]
      cat("  ", want, "\n      ",
          if (is.na(p$from)) paste0("not found: ", p$pattern, " under ", p$dir,
                                    "  (run ", p$script, ")") else p$from,
          "\n")
    }
  }
  cat("\nFigure 1 is a hand-drawn flowchart and has no generating script.\n")
  cat("\nnothing written; drop --dry-run to collect and merge\n")
  quit(save = "no")
}

missing <- character(0)
for (p in panels) {
  dest <- file.path(sub.dir, p$to)
  if (!is.na(p$from) && file.exists(p$from)) {
    file.copy(p$from, dest, overwrite = TRUE)
    cat("collected", p$to, "\n")
  } else {
    missing <- c(missing, sprintf("%s  (no %s under %s; run %s)",
                                  p$to, p$pattern, p$dir, p$script))
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
