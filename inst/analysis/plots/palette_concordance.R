## Palette and scale for the concordance heatmaps (Figure 2).
##
## The ramp is GnBu and the scale is capped at 0.5 so the off-diagonal
## structure stays visible; the stops below were read off the published
## figure so a rebuild reproduces it exactly.

concordance_colors <- colorRampPalette(c(
  "#DFF2F4", "#C6E9DC", "#B3E2C6", "#82CEC1",
  "#47AACC", "#3394C2", "#2282B9", "#1372B1"
))(200)

concordance_max    <- 0.5
concordance_breaks <- seq(0, concordance_max, length.out = 201)
concordance_legend <- seq(0, concordance_max, length.out = 6)

## Ramp for the pathway transition dotplot (Figure 4). The light end is the
## stage-2 q cut and the dark end the strongest enrichment; the first three
## stops share the concordance family, and lightness falls from L* 83 to 41 so
## the teal-to-warm hue change still reads as ordered.
enrichment_colors <- colorRampPalette(c(
  "#A8D8DC", "#7FC4A8", "#4FA89C", "#C97A4E", "#A03E5C"
))(256)
