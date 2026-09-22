# Writes a plot's legend on its own, for figures whose panels share one legend.

library(ggplot2)
library(grid)
library(gtable)

# Writes a legend grob to its own PDF at the grob's size, measured on a null device.
save_legend_grob <- function(grob, path, pad = 0.15, right = 0) {
  pdf(NULL)
  w <- convertWidth(grobWidth(grob), "in", valueOnly = TRUE)
  h <- convertHeight(grobHeight(grob), "in", valueOnly = TRUE)
  dev.off()
  # cairo, as for the panels
  cairo_pdf(paste0(path, ".pdf"), width = max(w, 1) + pad + right,
            height = max(h, 0.4) + pad)
  grid.newpage()
  # left-anchored, so an end label overhanging the bar falls into `right`
  pushViewport(viewport(x = unit(pad / 2, "in"), width = unit(max(w, 1), "in"),
                        just = "left"))
  grid.draw(grob)
  popViewport()
  dev.off()
  cat("Saving:", paste0(path, ".pdf"), "\n")
}

# Writes the legend of a ggplot, laid out horizontally on the given side.
save_plot_legend <- function(plot, path, side = "bottom") {
  pdf(NULL)
  g <- ggplotGrob(plot + theme(legend.position = side,
                               legend.box = "horizontal",
                               legend.direction = "horizontal"))
  dev.off()
  box <- g$grobs[[which(g$layout$name == paste0("guide-box-", side))]]
  # the outermost break label overhangs the colour bar, so the box gets a margin
  save_legend_grob(gtable::gtable_add_padding(box, unit(c(1, 4, 1, 4), "mm")), path)
}
