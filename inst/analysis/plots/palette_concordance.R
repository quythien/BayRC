## Palette and scale for the concordance heatmaps (Figure 6).
##
## A GnBu ramp, with the scale capped at 0.5 so the off-diagonal structure
## stays visible.

concordance_colors <- colorRampPalette(c(
  "#DFF2F4", "#C6E9DC", "#B3E2C6", "#82CEC1",
  "#47AACC", "#3394C2", "#2282B9", "#1372B1"
))(200)

concordance_max    <- 0.5
concordance_breaks <- seq(0, concordance_max, length.out = 201)
concordance_legend <- seq(0, concordance_max, length.out = 6)

## Ramp for the pathway transition dotplot (Figure 3), light at the stage-2 q
## cut and dark at the strongest enrichment; lightness falls from L* 83 to 41.
enrichment_colors <- colorRampPalette(c(
  "#A8D8DC", "#7FC4A8", "#4FA89C", "#C97A4E", "#A03E5C"
))(256)
