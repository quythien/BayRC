# Builds man/figures/circular_hdi_demo.png from the real manuscript-scale
# OMF posterior (mcmc_rho_BF3.RData / mcmc_phi_BF3.RData, 2,001 iterations;
# not bundled with the package). Load those two files and build mcmc_OMF via
# match_symbols() first (see omf_thr_full_pipeline.R for the exact steps),
# or load a cached copy with mcmc_OMF already built.

suppressMessages(library(ggplot2))

genes <- c("BMAL1", "NR1D1", "DBP", "FAM76B")

bf_OMF <- summarize_bay(mcmc_OMF$rho, BF = 3, p_rhythmic = 0.2)

panels <- list()
for (g in genes) {
  i <- match(g, mcmc_OMF$gname)
  rho <- mcmc_OMF$rho[i, ]
  phi_rhy <- mcmc_OMF$phi[i, rho == 1]
  p_rhy <- mean(rho)
  bf <- bf_OMF[g, "BayesF"]
  hdi <- circular_HDI(phi_rhy, credMass = 0.95, P = 24)
  med <- circular_median(phi_rhy, P = 24)
  panels[[g]] <- list(phi = phi_rhy, hdi = hdi, med = med, p_rhy = p_rhy, bf = bf)
}

make_panel <- function(g) {
  p <- panels[[g]]
  df <- data.frame(phi = p$phi)

  breaks <- seq(0, 24, by = 2)
  counts <- table(cut(df$phi, breaks = breaks, right = FALSE, include.lowest = TRUE))
  max_h  <- max(counts)
  ring_r <- max_h * 1.25

  if (p$hdi$upper >= p$hdi$lower) {
    arc_df <- data.frame(x = seq(p$hdi$lower, p$hdi$upper, length.out = 100), seg = 1)
  } else {
    # Wrap-around HDI: two disjoint segments (lower->24 and 0->upper) need
    # distinct `seg` groups, or geom_path connects them straight across the
    # circle instead of leaving the non-HDI majority of the circle empty.
    arc_df <- rbind(data.frame(x = seq(p$hdi$lower, 24, length.out = 60), seg = 1),
                     data.frame(x = seq(0, p$hdi$upper, length.out = 60), seg = 2))
  }
  arc_df$y <- ring_r

  ggplot(df, aes(x = phi)) +
    geom_histogram(breaks = breaks, fill = "#4C79A6", color = "white",
                    boundary = 0, closed = "left") +
    geom_path(data = arc_df, aes(x = x, y = y, group = seg), color = "#C0392B", linewidth = 3,
              lineend = "round") +
    geom_point(data = data.frame(x = p$med, y = ring_r), aes(x = x, y = y),
               color = "black", size = 2, shape = 18) +
    scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 6),
                        labels = c("ZT0", "ZT6", "ZT12", "ZT18", "ZT24")) +
    scale_y_continuous(limits = c(0, ring_r * 1.15)) +
    coord_polar(start = 0) +
    labs(title = sprintf("%s  (P[rhythmic] = %.2f, BF %s)", g, p$p_rhy,
                          if (p$bf >= 1e6) "> 1e6 (unbounded)" else sprintf("= %.0f", p$bf)),
         subtitle = sprintf("median = %.1fh, 95%% circular HDI = [%.1f, %.1f]h  (red arc; black = median)",
                             p$med, p$hdi$lower, p$hdi$upper),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
          plot.subtitle = element_text(size = 8, color = "#C0392B", hjust = 0.5))
}

plots <- lapply(genes, make_panel)

png("man/figures/circular_hdi_demo.png", width = 1600, height = 1550, res = 200, type = "cairo")
gridExtra::grid.arrange(grobs = plots, ncol = 2,
  top = grid::textGrob("Posterior phase distribution and 95% circular HDI (baboon OMF, manuscript-scale posterior, 4 genes)",
                        gp = grid::gpar(fontsize = 12, fontface = "bold")))
dev.off()
cat("saved man/figures/circular_hdi_demo.png\n")
for (g in genes) cat(g, ": HDI [", round(panels[[g]]$hdi$lower, 2), ",", round(panels[[g]]$hdi$upper, 2),
                     "]  median=", round(panels[[g]]$med, 2), "  P(rhy)=", round(panels[[g]]$p_rhy, 2), "\n")
