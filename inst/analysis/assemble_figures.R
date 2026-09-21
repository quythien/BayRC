################################################################################
# Collect the paper's figure panels and assemble the numbered figures.
#
# The analysis scripts each drop their panels somewhere under BAYRC_OUTPUT_DIR,
# under names that describe the comparison rather than the figure. This script
# takes the step from those panels to Figure_2 ... Figure_6: it copies every
# panel to paper/subfigures/ under a stable name, then letters and merges them
# into the numbered figure.
#
# Usage:
#   Rscript assemble_figures.R [paper.dir]
#
# paper.dir defaults to <BAYRC_OUTPUT_DIR>/../paper.
#
# Panels run along a row, or down the page for the figures named in
# stacked.figures. Both layouts go through pdflatex and graphicx, so the whole
# assembly is reproducible from the panels the analysis scripts write.
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
putsun.dir <- file.path(BAYRC_FIGURE_DIR, "baboon_PUT_SUN")
putvic.dir <- file.path(BAYRC_FIGURE_DIR, "baboon_PUT_VIC")
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
  list(dir = putsun.dir, pattern = "^Baboon_PUT_SUN_Peak_Concordance[.]pdf$",
       to = "F3B_baboon_PUT_SUN_phase_concordance.pdf",
       script = "applications/Baboon_PUT_SUN.R"),
  list(dir = putvic.dir, pattern = "^Baboon_PUT_VIC_Peak_Concordance[.]pdf$",
       to = "F3C_baboon_PUT_VIC_phase_concordance.pdf",
       script = "applications/Baboon_PUT_VIC.R"),

  # Figure 4: pathway transition enrichment in the two putamen circuits
  list(dir = putsun.dir, pattern = "^PUT_SUN_transition_enrichment[.]pdf$",
       to = "F4A_PUT_SUN_transition_enrichment.pdf",
       script = "applications/Baboon_PUT_SUN.R"),
  list(dir = putvic.dir, pattern = "^PUT_VIC_transition_enrichment[.]pdf$",
       to = "F4B_PUT_VIC_transition_enrichment.pdf",
       script = "applications/Baboon_PUT_VIC.R"),
  list(dir = putvic.dir, pattern = "^transition_enrichment_legend[.]pdf$",
       to = "F4L_transition_enrichment_legend.pdf",
       script = "applications/Baboon_PUT_VIC.R"),

  # Figure 5: KEGG Parkinson disease in both circuits, the rhythmic-only version
  list(dir = putsun.dir,
       pattern = "^KEGG_Parkinson_disease_integrated_rhythmic_only[.]pdf$",
       to = "F5A_PUT_SUN_KEGG_Parkinson_heatmap.pdf",
       script = "applications/Baboon_PUT_SUN.R"),
  list(dir = putvic.dir,
       pattern = "^KEGG_Parkinson_disease_integrated_rhythmic_only[.]pdf$",
       to = "F5B_PUT_VIC_KEGG_Parkinson_heatmap.pdf",
       script = "applications/Baboon_PUT_VIC.R"),
  list(dir = putvic.dir, pattern = "^parkinson_heatmap_legend[.]pdf$",
       to = "F5L_parkinson_heatmap_legend.pdf",
       script = "applications/Baboon_PUT_VIC.R"),

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
               "F3B_baboon_PUT_SUN_phase_concordance.pdf",
               "F3C_baboon_PUT_VIC_phase_concordance.pdf"),
  Figure_4 = c("F4A_PUT_SUN_transition_enrichment.pdf",
               "F4B_PUT_VIC_transition_enrichment.pdf"),
  Figure_5 = c("F5A_PUT_SUN_KEGG_Parkinson_heatmap.pdf",
               "F5B_PUT_VIC_KEGG_Parkinson_heatmap.pdf"),
  Figure_5_row = c("F5A_PUT_SUN_KEGG_Parkinson_heatmap.pdf",
                   "F5B_PUT_VIC_KEGG_Parkinson_heatmap.pdf"),
  Figure_6 = c("F6A_baboon_human_LUN_phase_concordance.pdf",
               "F6B_baboon_human_LUN_circadian_heatmap.pdf")
)

# Figure 5's two heatmaps are each as wide as the text block, so they stack.
# Figure_5_row holds the same two panels along a row for comparison; every
# other multi-panel figure runs its panels along a row.
stacked.figures <- "Figure_5"

# A figure listed here is drawn with its panels carrying no legend of their own
# and this one placed under the row. The panels' scales have to match for that
# to be right, which for Figure 4 is what q_limits and size_limits fix.
shared.legends <- list(Figure_4 = "F4L_transition_enrichment_legend.pdf",
                       Figure_5 = "F5L_parkinson_heatmap_legend.pdf",
                       Figure_5_row = "F5L_parkinson_heatmap_legend.pdf")

# Collect ---------------------------------------------------------------------

if (dry.run) {
  cat("resolved panels\n\n")
  for (nm in names(figures)) {
    cat(nm, "\n")
    for (want in c(figures[[nm]], shared.legends[[nm]])) {
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
# Panels are lettered and placed by pdflatex through graphicx. Every panel is
# scaled to one width, so a row is as tall as its tallest panel and a stack is
# as tall as its panels together.

page_size <- function(pdf) {
  info <- system2("pdfinfo", shQuote(pdf), stdout = TRUE, stderr = FALSE)
  d <- sub(".*: *", "", grep("^Page size", info, value = TRUE)[1])
  as.numeric(strsplit(sub(" pts.*", "", d), " x ")[[1]])
}

side_by_side <- function(inputs, output, labels = LETTERS[seq_along(inputs)],
                         panel.width = 324, gutter = 9, margin = 9,
                         label.space = 22, below = NA_character_) {
  if (!nzchar(Sys.which("pdflatex")) || !nzchar(Sys.which("pdfinfo")))
    return(FALSE)
  size <- lapply(inputs, page_size)
  scaled.h <- vapply(size, function(d) panel.width * d[2] / d[1], numeric(1))
  paper.w <- length(inputs) * panel.width +
             (length(inputs) - 1) * gutter + 2 * margin
  paper.h <- max(scaled.h) + label.space + 2 * margin

  # a legend the panels share sits centred under the row, at its own width
  strip <- ""
  if (!is.na(below)) {
    d <- page_size(below)
    strip.w <- min(d[1], paper.w - 2 * margin)
    paper.h <- paper.h + strip.w * d[2] / d[1] + gutter
    strip <- sprintf("\n\n\\vspace{%.1fbp}\\centerline{\\includegraphics[width=%.1fbp]{%s}}",
                     gutter, strip.w, below)
  }

  panel <- function(i) sprintf(
    "\\begin{minipage}[t]{%.1fbp}\\raggedright\\textbf{\\sffamily\\large %s}\\\\[2bp]\n\\includegraphics[width=%.1fbp]{%s}\\end{minipage}",
    panel.width, labels[i], panel.width, inputs[i])

  tex <- c("\\documentclass[11pt]{article}",
    sprintf("\\usepackage[paperwidth=%.1fbp,paperheight=%.1fbp,margin=%.1fbp]{geometry}",
            paper.w, paper.h, margin),
    "\\usepackage{graphicx}", "\\pagestyle{empty}",
    "\\setlength{\\parindent}{0pt}",
    "\\begin{document}\\noindent",
    paste0(paste(vapply(seq_along(inputs), panel, character(1)),
                 collapse = sprintf("\\hspace{%.1fbp}\n", gutter)), strip),
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

stacked <- function(inputs, output, labels = LETTERS[seq_along(inputs)],
                    panel.width = 468, gap = 14, margin = 9,
                    label.space = 22, below = NA_character_) {
  if (!nzchar(Sys.which("pdflatex")) || !nzchar(Sys.which("pdfinfo")))
    return(FALSE)
  size <- lapply(inputs, page_size)
  scaled.h <- vapply(size, function(d) panel.width * d[2] / d[1], numeric(1))
  paper.w <- panel.width + 2 * margin
  paper.h <- sum(scaled.h) + length(inputs) * label.space +
             (length(inputs) - 1) * gap + 2 * margin

  # a legend the panels share sits centred under the stack, at its own width
  strip <- ""
  if (!is.na(below)) {
    d <- page_size(below)
    strip.w <- min(d[1], panel.width)
    paper.h <- paper.h + strip.w * d[2] / d[1] + gap
    strip <- sprintf("\\\\[%.1fbp]\n\\centerline{\\includegraphics[width=%.1fbp]{%s}}",
                     gap, strip.w, below)
  }

  panel <- function(i) sprintf(
    "\\textbf{\\sffamily\\large %s}\\\\[2bp]\n\\includegraphics[width=%.1fbp]{%s}",
    labels[i], panel.width, inputs[i])

  tex <- c("\\documentclass[11pt]{article}",
    sprintf("\\usepackage[paperwidth=%.1fbp,paperheight=%.1fbp,margin=%.1fbp]{geometry}",
            paper.w, paper.h, margin),
    "\\usepackage{graphicx}", "\\pagestyle{empty}",
    "\\setlength{\\parindent}{0pt}",
    "\\begin{document}\\noindent",
    paste0(paste(vapply(seq_along(inputs), panel, character(1)),
                 collapse = sprintf("\\\\[%.1fbp]\n", gap)), strip),
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
  down <- nm %in% stacked.figures
  legend <- file.path(sub.dir, shared.legends[[nm]])
  legend <- if (length(legend) && file.exists(legend)) legend else NA_character_
  merge_panels <- if (down) function(...) stacked(..., below = legend) else
    function(...) side_by_side(..., below = legend)
  if (length(have) == 1L && is.na(legend)) {
    file.copy(have, out, overwrite = TRUE)
    cat("wrote", basename(out), "\n")
  } else if (isTRUE(merge_panels(have, out))) {
    cat("wrote", basename(out), "from", length(have),
        if (down) "panels stacked" else "panels side by side",
        if (is.na(legend)) "\n" else "under a shared legend\n")
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
