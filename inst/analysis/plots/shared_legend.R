# Writes a plot's legend on its own, for figures whose panels share one legend.

library(ggplot2)
library(grid)
library(gtable)

# Writes a finished legend grob to its own PDF at the size it asks for. The
# measuring needs a device, and this one keeps it off the default that would
# otherwise be opened in the working directory.
save_legend_grob <- function(grob, path, pad = 0.15) {
  pdf(NULL)
  w <- convertWidth(grobWidth(grob), "in", valueOnly = TRUE)
  h <- convertHeight(grobHeight(grob), "in", valueOnly = TRUE)
  dev.off()
  # cairo carries the glyphs the labels use, and matches the device the panels
  # themselves are drawn on
  cairo_pdf(paste0(path, ".pdf"), width = max(w, 1) + pad,
            height = max(h, 0.4) + pad)
  grid.newpage()
  grid.draw(grob)
  dev.off()
  cat("Saving:", paste0(path, ".pdf"), "\n")
}

# The panels of such a figure are drawn with no legend of their own, so the
# legend has to come from a copy of one panel laid out with its guides at the
# bottom. Only the box on the side the guides were sent to holds them.
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
