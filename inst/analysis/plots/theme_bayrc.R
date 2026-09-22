## Shared look for the paper figures. Source this at the top of every panel
## script so Figures 2-6 agree on palette, type and spacing.

library(ggplot2)
library(grDevices)

## Nimbus Sans is the Helvetica clone and is the only installed family with a
## real bold face; plain "sans" falls back to DejaVu, which is wider.
bayrc_family <- local({
  have <- tryCatch(any(grepl("Nimbus Sans", system("fc-list :family", intern = TRUE))),
                   error = function(e) FALSE)
  if (have) "Nimbus Sans" else "sans"
})

## the plain pdf() and postscript() devices need the family registered; Nimbus
## Sans takes Helvetica's metrics
if (bayrc_family == "Nimbus Sans" &&
    !"Nimbus Sans" %in% names(grDevices::pdfFonts())) {
  grDevices::pdfFonts("Nimbus Sans" = grDevices::Type1Font(
    "Nimbus Sans", grDevices::pdfFonts()$Helvetica$metrics))
  grDevices::postscriptFonts("Nimbus Sans" = grDevices::Type1Font(
    "Nimbus Sans", grDevices::postscriptFonts()$Helvetica$metrics))
}

bayrc_ink   <- "#0b0b0b"   # all text, axis lines and ticks
bayrc_ink2  <- "#52514e"   # mid grey, for marks that sit behind text
bayrc_ink3  <- "#7a7974"   # reference lines
bayrc_grid  <- "#e4e3de"   # gridlines

## two-series categorical, then the ordered set used for nested levels
bayrc_pair   <- c(all = "#9dbfe0", highlight = "#eaa37c")
bayrc_levels <- c("#5b7fb3", "#6fb3a0", "#dcb662", "#cf95b4")

## sequential ramp for concordance-style heatmaps on [0, 1]
bayrc_seq <- function(n = 256)
  colorRampPalette(c("#f7f9fc", "#c9dcef", "#8fb4d9", "#5b7fb3", "#33507a"))(n)

## diverging ramp for signed quantities such as a phase shift
bayrc_div <- function(n = 256)
  colorRampPalette(c("#5b7fb3", "#a8c3de", "#f4f2ee", "#eaba86", "#d1854f"))(n)

base_size <- 12

theme_bayrc <- function(base_size = 12) {
  theme_minimal(base_size = base_size, base_family = bayrc_family) +
    theme(
      plot.background   = element_rect(fill = NA, colour = NA),
      panel.background  = element_rect(fill = NA, colour = NA),
      panel.grid        = element_blank(),
      axis.line         = element_line(colour = bayrc_ink, linewidth = 0.5),
      axis.ticks        = element_line(colour = bayrc_ink, linewidth = 0.5),
      axis.ticks.length = unit(3.5, "pt"),
      axis.text         = element_text(colour = bayrc_ink, size = base_size * 0.85),
      axis.title        = element_text(colour = bayrc_ink, size = base_size),
      plot.title        = element_text(colour = bayrc_ink, size = base_size * 1.05,
                                       face = "bold", hjust = 0.5),
      plot.subtitle     = element_text(colour = bayrc_ink, size = base_size * 0.85,
                                       hjust = 0.5),
      plot.tag          = element_text(colour = bayrc_ink, size = base_size * 1.3,
                                       face = "bold"),
      legend.background = element_rect(fill = NA, colour = NA),
      legend.key        = element_rect(fill = NA, colour = NA),
      legend.title      = element_text(colour = bayrc_ink, size = base_size * 0.85),
      legend.text       = element_text(colour = bayrc_ink, size = base_size * 0.8),
      strip.text        = element_text(colour = bayrc_ink, size = base_size,
                                       face = "bold")
    )
}

## pheatmap and ComplexHeatmap do not read the ggplot theme; pass these instead
bayrc_heat_args <- function(diverging = FALSE, n = 256) {
  list(color = if (diverging) bayrc_div(n) else bayrc_seq(n),
       border_color = NA,
       fontfamily = bayrc_family,
       fontsize = base_size * 0.85)
}

## one writer for every panel, so sizes and embedding stay consistent
bayrc_save <- function(plot, file, width = 7.1, height = 5.2, dpi = 600) {
  cairo_pdf(paste0(file, ".pdf"), width = width, height = height,
            family = bayrc_family, bg = "transparent")
  print(plot); dev.off()
  ggsave(paste0(file, ".png"), plot, width = width, height = height,
         dpi = dpi, bg = "transparent")
  invisible(file)
}
