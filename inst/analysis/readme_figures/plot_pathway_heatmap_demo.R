# Builds two pathway heatmaps from the bundled quickstart-scale OMF (omental
# fat) vs THR (thyroid) posterior (2,500 iterations). Assumes mcmc_OMF,
# mcmc_THR, kegg, trans and phase are in the session, built as in
# exploratory/quickstart_baboon_OMF_THR.R.
#
# Two loss-direction pathSelect() hits (padj < 0.20 among 298 testable
# pathways, nperm = 500); gene counts are at BFDR alpha = 0.20:
#  - KEGG Long-term depression: Q_loss = 0.0282,
#    4 gain, 2 maintained, 10 non-rhythmic of 16 matched genes.
#  - KEGG GnRH signaling pathway: Q_loss = 0.1620,
#    4 gain, 5 maintained, 13 non-rhythmic of 22 matched genes.

suppressMessages(library(ggplot2))

pathways <- list(
  "KEGG Long-term depression" = "man/figures/pathway_heatmap_demo_ltd.png",
  "KEGG GnRH signaling pathway" = "man/figures/pathway_heatmap_demo_gnrh.png"
)

for (pathway_name in names(pathways)) {
  # the blocks and the left legend are fixed widths, so the canvas holds both
  png(pathways[[pathway_name]], width = 2700, height = 2400, res = 200, type = "cairo")
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
