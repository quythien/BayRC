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
# this script only arranges files that other scripts produced
bayrc.needs.summary <- FALSE
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

fig2.dir   <- file.path(BAYRC_FIGURE_DIR, "figure2")
scnhip.dir <- file.path(BAYRC_FIGURE_DIR, "baboon_SCN_HIP")
sunput.dir <- file.path(BAYRC_FIGURE_DIR, "baboon_SUN_PUT")
lung.dir   <- file.path(BAYRC_FIGURE_DIR, "baboon_human_LUN")

panels <- list(
  # Figure 2: genome-wide and KEGG-circadian concordance heatmaps
  list(dir = fig2.dir, pattern = "^Fig2A_genomewide[.]pdf$",
       to = "F2A_baboon_genomewide_concordance.pdf",
       script = "plots/heatmap_baboon.R then plots/replot_figure2.R"),
  list(dir = fig2.dir, pattern = "^Fig2B_circadian[.]pdf$",
       to = "F2B_baboon_circadian_concordance.pdf",
       script = "plots/heatmap_circadian_pairs.R within_baboon then plots/replot_figure2.R"),

  # Figure 3: within-species phase concordance scatters
  list(dir = scnhip.dir, pattern = "^Baboon_SCN_HIP_Peak_Concordance[.]pdf$",
       to = "F3A_baboon_SCN_HIP_phase_concordance.pdf",
       script = "applications/Baboon_SCN_HIP.R"),
  list(dir = sunput.dir, pattern = "^Baboon_SUN_PUT_Peak_Concordance[.]pdf$",
       to = "F3B_baboon_SUN_PUT_phase_concordance.pdf",
       script = "applications/Baboon_SUN_PUT.R"),

  # Figure 4: SUN-PUT pathway transition enrichment, a single panel
  list(dir = sunput.dir, pattern = "^SUN_PUT_transition_enrichment[.]pdf$",
       to = "F4_SUN_PUT_transition_enrichment.pdf",
       script = "applications/Baboon_SUN_PUT.R"),

  # Figure 5: SUN-PUT pathway heatmaps, the rhythmic-only version
  list(dir = sunput.dir,
       pattern = "^KEGG_Parkinson_disease_integrated_rhythmic_only[.]pdf$",
       to = "F5A_SUN_PUT_KEGG_Parkinson_heatmap.pdf",
       script = "applications/Baboon_SUN_PUT.R"),
  list(dir = sunput.dir,
       pattern = "^KEGG_Oxidative_phosphorylation_integrated_rhythmic_only[.]pdf$",
       to = "F5B_SUN_PUT_KEGG_OxPhos_heatmap.pdf",
       script = "applications/Baboon_SUN_PUT.R"),

  # Figure 6: cross-species lung
  list(dir = lung.dir, pattern = "^Baboon_Human_LUN_Peak_Concordance[.]pdf$",
       to = "F6A_baboon_human_LUN_phase_concordance.pdf",
       script = "applications/Baboon_Human_LUN.R"),
  list(dir = lung.dir,
       pattern = "^KEGG_Circadian_rhythm_integrated_rhythmic_only[.]pdf$",
       to = "F6B_baboon_human_LUN_circadian_heatmap.pdf",
       script = "applications/Baboon_Human_LUN.R")
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
  Figure_4 = "F4_SUN_PUT_transition_enrichment.pdf",
  Figure_5 = c("F5A_SUN_PUT_KEGG_Parkinson_heatmap.pdf",
               "F5B_SUN_PUT_KEGG_OxPhos_heatmap.pdf"),
  Figure_6 = c("F6A_baboon_human_LUN_phase_concordance.pdf",
               "F6B_baboon_human_LUN_circadian_heatmap.pdf")
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
# Panels are placed side by side and lettered, which pdflatex does through
# graphicx. Each panel is scaled to the same width, so the page is as tall as
# the tallest scaled panel.

page_size <- function(pdf) {
  info <- system2("pdfinfo", shQuote(pdf), stdout = TRUE, stderr = FALSE)
  d <- sub(".*: *", "", grep("^Page size", info, value = TRUE)[1])
  as.numeric(strsplit(sub(" pts.*", "", d), " x ")[[1]])
}

side_by_side <- function(inputs, output, labels = LETTERS[seq_along(inputs)],
                         panel.width = 324, gutter = 9, margin = 9,
                         label.space = 22) {
  if (!nzchar(Sys.which("pdflatex")) || !nzchar(Sys.which("pdfinfo")))
    return(FALSE)
  size <- lapply(inputs, page_size)
  scaled.h <- vapply(size, function(d) panel.width * d[2] / d[1], numeric(1))
  paper.w <- length(inputs) * panel.width +
             (length(inputs) - 1) * gutter + 2 * margin
  paper.h <- max(scaled.h) + label.space + 2 * margin

  panel <- function(i) sprintf(
    "\\begin{minipage}[t]{%.1fbp}\\raggedright\\textbf{\\sffamily\\large %s}\\\\[2bp]\n\\includegraphics[width=%.1fbp]{%s}\\end{minipage}",
    panel.width, labels[i], panel.width, inputs[i])

  tex <- c("\\documentclass[11pt]{article}",
    sprintf("\\usepackage[paperwidth=%.1fbp,paperheight=%.1fbp,margin=%.1fbp]{geometry}",
            paper.w, paper.h, margin),
    "\\usepackage{graphicx}", "\\pagestyle{empty}",
    "\\setlength{\\parindent}{0pt}",
    "\\begin{document}\\noindent",
    paste(vapply(seq_along(inputs), panel, character(1)),
          collapse = sprintf("\\hspace{%.1fbp}\n", gutter)),
    "\\end{document}")

  work <- file.path(tempdir(), "assemble")
  dir.create(work, showWarnings = FALSE)
  writeLines(tex, file.path(work, "fig.tex"))
  system2("pdflatex", c("-interaction=batchmode", "-halt-on-error",
                        "-output-directory", shQuote(work),
                        shQuote(file.path(work, "fig.tex"))),
          stdout = FALSE, stderr = FALSE)
  built <- file.path(work, "fig.pdf")
  if (!file.exists(built)) return(FALSE)
  file.copy(built, output, overwrite = TRUE)
}

unassembled <- character(0)
for (nm in names(figures)) {
  want <- file.path(sub.dir, figures[[nm]])
  have <- want[file.exists(want)]
  if (!length(have)) next
  out <- file.path(fig.dir, paste0(nm, ".pdf"))
  if (length(have) == 1L) {
    file.copy(have, out, overwrite = TRUE)
    cat("wrote", basename(out), "\n")
  } else if (isTRUE(side_by_side(have, out))) {
    cat("wrote", basename(out), "from", length(have), "panels side by side\n")
  } else {
    unassembled <- c(unassembled, nm)
  }
  if (length(have) < length(want))
    cat("  ", nm, "is incomplete:", length(have), "of", length(want), "panels\n")
}

cat("\nFigure 1 is a hand-drawn flowchart and has no generating script.\n")

if (length(missing)) {
  cat("\nPanels not found:\n")
  cat(paste0("  ", missing, collapse = "\n"), "\n")
}
if (length(unassembled)) {
  cat("\nSide-by-side assembly needs pdflatex and pdfinfo on the PATH.",
      "Left unassembled:", paste(unassembled, collapse = ", "), "\n")
}
cat("\npanels:", sub.dir, "\nfigures:", fig.dir, "\n")
