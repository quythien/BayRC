# Builds two pathway heatmaps from the real manuscript-scale OMF (omental
# fat) vs THR (thyroid) posterior (2,001 iterations, mcmc_rho_BF3.RData /
# mcmc_phi_BF3.RData; not bundled with the package). Load those files and
# build mcmc_OMF, mcmc_THR via match_symbols() first (see
# omf_thr_full_pipeline.R for the exact steps), or load a cached copy with
# those objects, kegg, trans, and phase already built.
#
# Two pathways chosen from the real Stage 2 loss/conserved-direction
# pathSelect() results, both real KEGG hits (padj < 0.20 among 298
# testable pathways, nperm = 500). Gene-level gain/loss/maintained counts
# below are at BFDR alpha = 0.20 (pathSelect's own significance test does
# not depend on alpha; only these per-gene counts do):
#  - KEGG Long-term depression: strongest statistics (padj = 0.0018) and
#    a balanced gene-level mix (4 gain, 5 loss, 6 maintained, 4
#    non-rhythmic of 19 matched genes).
#  - KEGG Circadian rhythm: the weaker of the two statistically
#    (padj = 0.12, still under the 0.20 cutoff) but the clearest possible
#    thematic fit for a circadian-genomics package (2 gain, 2 loss, 9
#    maintained, 10 non-rhythmic of 23 matched genes).

suppressMessages(library(ggplot2))

pathways <- list(
  "KEGG Long-term depression" = "man/figures/pathway_heatmap_demo_ltd.png",
  "KEGG Circadian rhythm" = "man/figures/pathway_heatmap_demo_circadian.png"
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
