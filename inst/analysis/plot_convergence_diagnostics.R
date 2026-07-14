# Builds man/figures/convergence_trace_demo.png from the real ILE posterior
# and mcmc_diagnostics() output. Run quickstart_baboon_ILE_KIC.R first (or
# load a cached mcmc_ILE), then compute diag_ILE <- mcmc_diagnostics(mcmc_ILE)
# before sourcing this script.
#
# Picks one gene with high phi ESS (fast-mixing, tight trace) and one with
# low phi ESS (slow-mixing, sticky trace) to show what good vs. poor mixing
# actually looks like, rather than showing only a flattering example.

suppressMessages(library(ggplot2))

ess_phi <- diag_ILE$ess_phi
gene_hi <- names(sort(ess_phi, decreasing = TRUE))[1]
gene_lo <- names(sort(ess_phi))[1]

make_trace <- function(gene, label) {
  i <- match(gene, rownames(mcmc_ILE$phi))
  df <- data.frame(iter = seq_len(ncol(mcmc_ILE$phi)), phi = mcmc_ILE$phi[i, ],
                    rho = mcmc_ILE$rho[i, ])
  ggplot(df, aes(x = iter, y = phi, color = factor(rho))) +
    geom_line(aes(group = 1), color = "gray70", linewidth = 0.2) +
    geom_point(size = 0.4, alpha = 0.6) +
    scale_color_manual(values = c("0" = "gray70", "1" = "#4C79A6"), guide = "none") +
    scale_y_continuous(limits = c(0, 24), breaks = seq(0, 24, 6)) +
    labs(title = sprintf("%s (ESS phi = %.0f)", gene, ess_phi[gene]),
         subtitle = label, x = "retained MCMC iteration", y = "phi (h)") +
    theme_minimal(base_size = 10)
}

p1 <- make_trace(gene_hi, "fastest-mixing gene in ILE")
p2 <- make_trace(gene_lo, "slowest-mixing gene in ILE")

png("man/figures/convergence_trace_demo.png", width = 1800, height = 900, res = 200, type = "cairo")
gridExtra::grid.arrange(p1, p2, ncol = 2,
  top = grid::textGrob("Phase trace plots: fastest vs. slowest-mixing gene (baboon ILE)",
                        gp = grid::gpar(fontsize = 12, fontface = "bold")))
dev.off()
cat("saved man/figures/convergence_trace_demo.png\n")
cat("high ESS gene:", gene_hi, "=", ess_phi[gene_hi], "\n")
cat("low ESS gene:", gene_lo, "=", ess_phi[gene_lo], "\n")
