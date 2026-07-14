# Builds two pathway heatmaps from the bundled quickstart-scale OMF
# (omental fat) vs THR (thyroid) posterior (2,500 iterations, from
# inst/analysis/quickstart_baboon_OMF_THR.R). Assumes mcmc_OMF, mcmc_THR,
# kegg, trans, and phase are already built in the environment (see that
# script for the exact steps).
#
# Two pathways chosen from the loss-direction pathSelect() results, both
# real KEGG hits (padj < 0.20 among 220 testable pathways, nperm = 500).
# Gene-level gain/loss/maintained counts below are at BFDR alpha = 0.20
# (pathSelect's own significance test does not depend on alpha; only
# these per-gene counts do):
#  - KEGG Long-term depression: strongest statistics (padj = 0.0065),
#    4 gain, 4 loss, 4 maintained, 4 non-rhythmic of 16 matched genes.
#  - KEGG GnRH signaling pathway: padj = 0.0258, second-strongest hit,
#    3 gain, 2 loss, 6 maintained, 11 non-rhythmic of 22 matched genes.

suppressMessages(library(ggplot2))

pathways <- list(
  "KEGG Long-term depression" = "man/figures/pathway_heatmap_demo_ltd.png",
  "KEGG GnRH signaling pathway" = "man/figures/pathway_heatmap_demo_gnrh.png"
)

for (pathway_name in names(pathways)) {
  png(pathways[[pathway_name]], width = 2000, height = 2400, res = 200, type = "cairo")
  plot_heatmap(
    data1 = mcmc_OMF, data2 = mcmc_THR,
    pathway_genes = kegg[[pathway_name]],
    pathway_name  = pathway_name,
    phase_results = phase, transition_results = trans,
    group_names = c("OMF", "THR")
  )
  dev.off()
  cat("saved", pathways[[pathway_name]], "\n")
}
