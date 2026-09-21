# Writes a plot's legend on its own, for figures whose panels share one legend.

library(ggplot2)
library(grid)

# The panels of such a figure are drawn with no legend of their own, so the
# legend has to come from a copy of one panel laid out with its guides at the
# bottom. Only the box on the side the guides were sent to holds them.
save_plot_legend <- function(plot, path, side = "bottom") {
  g <- ggplotGrob(plot + theme(legend.position = side,
                               legend.box = "horizontal",
                               legend.direction = "horizontal"))
  box <- g$grobs[[which(g$layout$name == paste0("guide-box-", side))]]
  w <- convertWidth(sum(box$widths), "in", valueOnly = TRUE)
  h <- convertHeight(sum(box$heights), "in", valueOnly = TRUE)
  pdf(paste0(path, ".pdf"), width = max(w, 1), height = max(h, 0.4))
  grid.newpage()
  grid.draw(box)
  dev.off()
  cat("Saving:", paste0(path, ".pdf"), "\n")
}
