# Writes a plot's legend on its own, for figures whose panels share one legend.

library(ggplot2)
library(grid)
library(gtable)

# The panels of such a figure are drawn with no legend of their own, so the
# legend has to come from a copy of one panel laid out with its guides at the
# bottom. Only the box on the side the guides were sent to holds them.
save_plot_legend <- function(plot, path, side = "bottom") {
  g <- ggplotGrob(plot + theme(legend.position = side,
                               legend.box = "horizontal",
                               legend.direction = "horizontal"))
  box <- g$grobs[[which(g$layout$name == paste0("guide-box-", side))]]
  # the outermost break label overhangs the colour bar, so the box gets a margin
  box <- gtable::gtable_add_padding(box, unit(c(1, 4, 1, 4), "mm"))
  # measuring the box needs a device open, and the sizes then set the real one
  pdf(NULL)
  w <- convertWidth(sum(box$widths), "in", valueOnly = TRUE)
  h <- convertHeight(sum(box$heights), "in", valueOnly = TRUE)
  dev.off()
  pdf(paste0(path, ".pdf"), width = max(w, 1), height = max(h, 0.4))
  grid.newpage()
  grid.draw(box)
  dev.off()
  cat("Saving:", paste0(path, ".pdf"), "\n")
}
