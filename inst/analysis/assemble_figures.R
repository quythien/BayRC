################################################################################
# Collect the paper's figure panels and assemble the numbered figures.
#
# The analysis scripts each drop their panels somewhere under BAYRC_OUTPUT_DIR,
# under names that describe the comparison rather than the figure. This script
# takes the step from those panels to Figure_2 and Figure_4: it copies every
# panel to paper/subfigures/ under a stable name, then letters and merges them
# into the numbered figure. Figures 2 and 6 are laid out by their own scripts,
# named in external.figures, and only have their panels collected here.
#
# Usage:
#   Rscript assemble_figures.R [paper.dir]
#
# paper.dir defaults to <BAYRC_OUTPUT_DIR>/../paper.
#
# Panels run along a row, or down the page for the figures named in
# stacked.figures. Both layouts are built with pdflatex and graphicx.
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

# Panels are matched by directory and pattern rather than exact file name.

find_panel <- function(dir, pattern, recursive = FALSE) {
  if (!dir.exists(dir)) return(NA_character_)
  hits <- list.files(dir, pattern = pattern, full.names = TRUE,
                     recursive = recursive)
  if (!length(hits)) return(NA_character_)
  ## with several matches, take the newest
  hits <- hits[order(file.info(hits)$mtime, decreasing = TRUE)]
  if (length(hits) > 1) {
    cat("  ", length(hits), "files match", pattern, "under", dir, "\n")
    cat("     taking the newest:", basename(hits[1]), "\n")
    cat("     also matched:", paste(basename(hits[-1]), collapse = ", "), "\n")
  }
  hits[1]
}

fig6.dir   <- file.path(BAYRC_FIGURE_DIR, "figure6")
scnhip.dir <- file.path(BAYRC_FIGURE_DIR, "baboon_SCN_HIP")
putsun.dir <- file.path(BAYRC_FIGURE_DIR, "baboon_PUT_SUN")
putvic.dir <- file.path(BAYRC_FIGURE_DIR, "baboon_PUT_VIC")
lung.dir   <- file.path(BAYRC_FIGURE_DIR, "baboon_human_LUN")

panels <- list(
  # Figure 6: genome-wide and KEGG-circadian concordance heatmaps
  list(dir = fig6.dir, pattern = "^Fig6A_genomewide[.]pdf$",
       to = "F6A_baboon_genomewide_concordance.pdf",
       script = "plots/heatmap_baboon.R then plots/replot_figure6.R"),
  list(dir = fig6.dir, pattern = "^Fig6B_circadian[.]pdf$",
       to = "F6B_baboon_circadian_concordance.pdf",
       script = "plots/heatmap_circadian_pairs.R within_baboon then plots/replot_figure6.R"),
  list(dir = fig6.dir, pattern = "^Fig6C_circadian_membership[.]pdf$",
       to = "F6C_circadian_membership.pdf",
       script = "plots/figure6_panelC.R"),
  list(dir = fig6.dir, pattern = "^Fig6_concordance_legend[.]pdf$",
       to = "F6L_concordance_legend.pdf",
       script = "plots/replot_figure6.R"),

  # Figure 2: within-species phase concordance scatters
  list(dir = scnhip.dir, pattern = "^Baboon_SCN_HIP_Peak_Concordance[.]pdf$",
       to = "F2E_baboon_SCN_HIP_phase_concordance.pdf",
       script = "applications/Baboon_SCN_HIP.R"),
  list(dir = scnhip.dir, pattern = "^phase_class_legend[.]pdf$",
       to = "F2L_phase_class_legend.pdf",
       script = "applications/Baboon_SCN_HIP.R"),
  list(dir = putsun.dir, pattern = "^Baboon_PUT_SUN_Peak_Concordance[.]pdf$",
       to = "F2C_baboon_PUT_SUN_phase_concordance.pdf",
       script = "applications/Baboon_PUT_SUN.R"),
  list(dir = putvic.dir, pattern = "^Baboon_PUT_VIC_Peak_Concordance[.]pdf$",
       to = "F2D_baboon_PUT_VIC_phase_concordance.pdf",
       script = "applications/Baboon_PUT_VIC.R"),
  list(dir = scnhip.dir, pattern = "^SCN_clock_reference[.]pdf$",
       to = "SCN_clock_reference.pdf",
       script = "applications/Baboon_SCN_HIP.R"),

  # Figure 3 is one plot covering both circuits, written straight to the paper
  # figures by plots/pathway_transition_panels.R, so it has no panels here.

  # Figure 4: KEGG Parkinson disease in both circuits, the rhythmic-only version
  list(dir = putsun.dir,
       pattern = "^KEGG_Parkinson_disease_integrated_rhythmic_only[.]pdf$",
       to = "F4A_PUT_SUN_KEGG_Parkinson_heatmap.pdf",
       script = "applications/Baboon_PUT_SUN.R"),
  list(dir = putvic.dir,
       pattern = "^KEGG_Parkinson_disease_integrated_rhythmic_only[.]pdf$",
       to = "F4B_PUT_VIC_KEGG_Parkinson_heatmap.pdf",
       script = "applications/Baboon_PUT_VIC.R"),
  list(dir = putvic.dir, pattern = "^parkinson_heatmap_legend[.]pdf$",
       to = "F4L_parkinson_heatmap_legend.pdf",
       script = "applications/Baboon_PUT_VIC.R"),

  # Figure 5: cross-species lung, restricted to the common three-species set.
  list(dir = file.path(BAYRC_FIGURE_DIR, "figure5"),
       pattern = "^Fig5A_human_baboon[.]pdf$",
       to = "F5A_human_baboon_phase_concordance.pdf",
       script = "plots/figure5/make_figure5.R"),
  list(dir = file.path(BAYRC_FIGURE_DIR, "figure5"),
       pattern = "^Fig5B_human_mouse[.]pdf$",
       to = "F5B_human_mouse_phase_concordance.pdf",
       script = "plots/figure5/make_figure5.R"),
  list(dir = file.path(BAYRC_FIGURE_DIR, "figure5"),
       pattern = "^Fig5C_heatmap[.]pdf$",
       to = "F5C_lung_circadian_heatmap.pdf",
       script = "plots/figure5/make_figure5.R")
)

for (i in seq_along(panels))
  panels[[i]]$from <- find_panel(panels[[i]]$dir, panels[[i]]$pattern,
                                 isTRUE(panels[[i]]$recursive))

# Which collected panels make up each numbered figure -------------------------

figures <- list(
  Figure_6 = c("F6A_baboon_genomewide_concordance.pdf",
               "F6B_baboon_circadian_concordance.pdf",
               "F6C_circadian_membership.pdf"),
  Figure_2 = c("F2C_baboon_PUT_SUN_phase_concordance.pdf",
               "F2D_baboon_PUT_VIC_phase_concordance.pdf",
               "F2E_baboon_SCN_HIP_phase_concordance.pdf"),
  Figure_4 = c("F4A_PUT_SUN_KEGG_Parkinson_heatmap.pdf",
               "F4B_PUT_VIC_KEGG_Parkinson_heatmap.pdf"),
  Figure_4_row = c("F4A_PUT_SUN_KEGG_Parkinson_heatmap.pdf",
                   "F4B_PUT_VIC_KEGG_Parkinson_heatmap.pdf"),
  Figure_5 = c("F5A_human_baboon_phase_concordance.pdf",
               "F5B_human_mouse_phase_concordance.pdf",
               "F5C_lung_circadian_heatmap.pdf")
)

# Figures laid out by their own script; this one only collects their panels.
external.figures <- c(Figure_6 = "plots/figure6/assemble_figure6_nature.py",
                      Figure_5 = "plots/figure5/assemble_figure5.py")

# Figure 4's heatmaps are each as wide as the text block, so they stack;
# Figure_4_row lays the same two panels along a row.
stacked.figures <- "Figure_4"

# A figure listed here wraps its panels into rows of this many.
figure.columns <- list(Figure_6 = 2, Figure_2 = 3)
# a legend that belongs to one panel rather than the row sits under that panel
legend.under <- list(Figure_5 = 2L, Figure_4_row = 1L)
# a legend that belongs to the panels in one row sits directly beneath that row
# the colour bar belongs to the two heatmaps, not to the membership panel below
legend.after.row <- list(Figure_6 = 1L, Figure_2 = 1L)

# caption beside each panel letter, under a title the panels share
panel.captions <- list(
  Figure_4     = c("Putamen versus substantia nigra",
                   "Putamen versus visual cortex"),
  Figure_4_row = c("Putamen versus substantia nigra",
                   "Putamen versus visual cortex"))

figure.titles <- list(Figure_4 = "KEGG Parkinson disease",
                      Figure_4_row = "KEGG Parkinson disease",
                      Figure_5 = "Cross-species lung")

shared.legends <- list(Figure_6 = "F6L_concordance_legend.pdf",
                       Figure_2 = "F2L_phase_class_legend.pdf",
                       Figure_4 = "F4L_parkinson_heatmap_legend.pdf",
                       Figure_4_row = "F4L_parkinson_heatmap_legend.pdf")

# Collect ---------------------------------------------------------------------

if (dry.run) {
  cat("resolved panels\n\n")
  for (nm in names(figures)) {
    cat(nm, if (nm %in% names(external.figures))
      paste0("  (panels only; laid out by ", external.figures[[nm]], ")"), "\n")
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

# Point size that prints at label.pt once the page is scaled to the text width
label_size <- function(paper.w, label.pt = 8.3, textwidth = 488.5)
  label.pt * paper.w / textwidth

side_by_side <- function(inputs, output, labels = LETTERS[seq_along(inputs)],
                         panel.width = 324, gutter = 9, margin = 9,
                         label.space = 28, below = NA_character_,
                         below.under = NA_integer_,
                         below.after = NA_integer_,
                         title = NA_character_, ncol = NULL) {
  if (!nzchar(Sys.which("pdflatex")) || !nzchar(Sys.which("pdfinfo")))
    return(FALSE)
  if (is.null(ncol)) ncol <- length(inputs)
  size <- lapply(inputs, page_size)
  # panels fill each row in turn; a row with fewer panels spreads them over the
  # same content width, and a row is as tall as its tallest panel
  row.of <- ceiling(seq_along(inputs) / ncol)
  content.w <- min(ncol, length(inputs)) * panel.width +
               (min(ncol, length(inputs)) - 1) * gutter
  paper.w <- content.w + 2 * margin
  width.of <- vapply(row.of, function(r) {
    k <- sum(row.of == r); (content.w - (k - 1) * gutter) / k }, numeric(1))
  scaled.h <- vapply(seq_along(inputs), function(i)
    width.of[i] * size[[i]][2] / size[[i]][1], numeric(1))
  paper.h <- sum(tapply(scaled.h, row.of, max)) +
             length(unique(row.of)) * label.space +
             (length(unique(row.of)) - 1) * gutter + 2 * margin
  head <- ""
  if (!is.na(title)) {
    paper.h <- paper.h + 24
    head <- sprintf(paste0("\\centerline{\\sffamily\\bfseries",
                             sprintf("\\fontsize{%.1f}{%.1f}\\selectfont",
                                     label_size(paper.w, 10.5),
                                     label_size(paper.w, 10.5) * 1.2),
                             " %s}\\vspace{8bp}\n\n\\noindent "), title)
  }

  # a legend the panels share sits centred under the row at its own width, or
  # under one named panel when below.under gives that panel's index
  strip <- ""
  if (!is.na(below)) {
    d <- page_size(below)
    if (length(below.under) == 1 && !is.na(below.under)) {
      strip.w <- width.of[below.under]
      lead <- (below.under - 1) * (width.of[below.under] + gutter)
      strip <- sprintf("\n\n\\vspace{%.1fbp}\\noindent\\hspace*{%.1fbp}\\includegraphics[width=%.1fbp]{%s}",
                       gutter, lead, strip.w, below)
    } else {
      strip.w <- min(d[1], paper.w - 2 * margin)
      strip <- sprintf("\n\n\\vspace{%.1fbp}\\centerline{\\includegraphics[width=%.1fbp]{%s}}",
                       gutter, strip.w, below)
    }
    paper.h <- paper.h + strip.w * d[2] / d[1] + gutter
  }

  # A legend after a given row sits inline in a content-width minipage, since
  # \centerline would end the paragraph and drop the later rows.
  mid.strip <- character(0)
  if (length(below.after) == 1 && !is.na(below.after) && !is.na(below)) {
    d <- page_size(below)
    strip.w <- min(d[1], content.w)
    mid.strip <- sprintf(
      "\\begin{minipage}{%.1fbp}\\centering\\includegraphics[width=%.1fbp]{%s}\\end{minipage}",
      content.w, strip.w, below)
    strip <- ""
    # the strip is a row of its own, so it takes one more separator
    paper.h <- paper.h + gutter
  } else {
    below.after <- length(unique(row.of))
  }

  lab.pt <- label_size(paper.w)
  panel <- function(i) sprintf(
    paste0("\\begin{minipage}[t]{%.1fbp}\\raggedright\\textbf{\\sffamily",
           sprintf("\\fontsize{%.1f}{%.1f}\\selectfont", lab.pt, lab.pt * 1.2),
           " %s}\\\\[2bp]\n\\includegraphics[width=%.1fbp]{%s}\\end{minipage}"),
    width.of[i], labels[i], width.of[i], inputs[i])

  tex <- c("\\documentclass[11pt]{article}",
    sprintf("\\usepackage[paperwidth=%.1fbp,paperheight=%.1fbp,margin=%.1fbp]{geometry}",
            paper.w, paper.h, margin),
    "\\usepackage{graphicx}", "\\pagestyle{empty}",
    "\\setlength{\\parindent}{0pt}",
    "\\begin{document}\\noindent",
    paste0(head, paste(append(
             vapply(split(seq_along(inputs), row.of), function(ix)
               paste(vapply(ix, panel, character(1)),
                     collapse = sprintf("\\hspace{%.1fbp}\n", gutter)),
               character(1)),
             # a legend belonging to the panels above it goes between the rows
             # rather than at the foot of the figure
             values = mid.strip, after = below.after),
             collapse = sprintf("\\\\[%.1fbp]\n", gutter)), strip),
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
                    label.space = 28, below = NA_character_,
                    title = NA_character_) {
  if (!nzchar(Sys.which("pdflatex")) || !nzchar(Sys.which("pdfinfo")))
    return(FALSE)
  size <- lapply(inputs, page_size)
  scaled.h <- vapply(size, function(d) panel.width * d[2] / d[1], numeric(1))
  paper.w <- panel.width + 2 * margin
  paper.h <- sum(scaled.h) + length(inputs) * label.space +
             (length(inputs) - 1) * gap + 2 * margin
  head <- ""
  if (!is.na(title)) {
    paper.h <- paper.h + 24
    head <- sprintf(paste0("\\centerline{\\sffamily\\bfseries",
                             sprintf("\\fontsize{%.1f}{%.1f}\\selectfont",
                                     label_size(paper.w, 10.5),
                                     label_size(paper.w, 10.5) * 1.2),
                             " %s}\\vspace{8bp}\n\n\\noindent "), title)
  }

  # a legend the panels share sits centred under the stack, at its own width
  strip <- ""
  if (!is.na(below)) {
    d <- page_size(below)
    strip.w <- min(d[1], panel.width)
    paper.h <- paper.h + strip.w * d[2] / d[1] + gap
    strip <- sprintf("\\\\[%.1fbp]\n\\centerline{\\includegraphics[width=%.1fbp]{%s}}",
                     gap, strip.w, below)
  }

  lab.pt <- label_size(paper.w)
  panel <- function(i) sprintf(
    paste0("\\textbf{\\sffamily",
           sprintf("\\fontsize{%.1f}{%.1f}\\selectfont", lab.pt, lab.pt * 1.2),
           " %s}\\\\[2bp]\n\\includegraphics[width=%.1fbp]{%s}"),
    labels[i], panel.width, inputs[i])

  tex <- c("\\documentclass[11pt]{article}",
    sprintf("\\usepackage[paperwidth=%.1fbp,paperheight=%.1fbp,margin=%.1fbp]{geometry}",
            paper.w, paper.h, margin),
    "\\usepackage{graphicx}", "\\pagestyle{empty}",
    "\\setlength{\\parindent}{0pt}",
    "\\begin{document}\\noindent",
    paste0(head, paste(vapply(seq_along(inputs), panel, character(1)),
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
  if (nm %in% names(external.figures)) {
    cat(nm, "is laid out by", external.figures[[nm]], "\n")
    next
  }
  want <- file.path(sub.dir, figures[[nm]])
  have <- want[file.exists(want)]
  if (!length(have)) next
  if (length(have) != length(want)) {
    unassembled <- c(unassembled, nm)
    next
  }
  out <- file.path(fig.dir, paste0(nm, ".pdf"))
  down <- nm %in% stacked.figures
  legend <- file.path(sub.dir, shared.legends[[nm]])
  legend <- if (length(legend) && file.exists(legend)) legend else NA_character_
  title <- if (is.null(figure.titles[[nm]])) NA_character_ else figure.titles[[nm]]
  caps <- panel.captions[[nm]]
  # the letter keeps the left margin and takes no width, so the caption beside
  # it centres over the panel rather than starting where the letter ends
  lab <- if (is.null(caps)) LETTERS[seq_along(have)] else
    mapply(function(l, c)
             sprintf("\\makebox[0pt][l]{%s}\\makebox[\\linewidth]{%s}", l, c),
           LETTERS[seq_along(have)], caps[seq_along(have)], USE.NAMES = FALSE)
  if (nm == "Figure_2") lab <- LETTERS[3:5]
  merge_panels <- if (down)
    function(...) stacked(..., labels = lab, below = legend, title = title) else
    function(...) side_by_side(..., labels = lab, below = legend, title = title,
                              ncol = figure.columns[[nm]],
                              below.under = legend.under[[nm]],
                              below.after = legend.after.row[[nm]])
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
