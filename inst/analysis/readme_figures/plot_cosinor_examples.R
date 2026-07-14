# Builds man/figures/cosinor_examples.png: raw expression scatter + OLS
# cosinor fit (BayRC's one_cosinor_OLS(), the classical frequentist
# baseline) for one example gene per transition/phase category, baboon
# OMF vs THR, on the bundled quickstart-scale posterior (2,500
# iterations, from inst/analysis/quickstart_baboon_OMF_THR.R). Assumes
# mcmc_OMF, mcmc_THR, trans, phase are already built in the environment
# (see quickstart_baboon_OMF_THR.R for the exact steps).
#
# trans$gain_loss_status and phase$flag_* are already named by gene
# symbol (inherited from rownames(mcmc_OMF$rho)/rownames(mcmc_OMF$phi),
# set by match_symbols()); do not try to re-name them from mcmc_OMF$gname,
# which does not exist on a match_symbols()-processed object and is NULL.

suppressMessages(library(ggplot2))

baboon <- readRDS(system.file("extdata", "baboon_OMF_THR_GSE98965.rds", package = "BayRC"))
zt <- baboon$zt

# A handful of gene symbols in this dataset are Excel date-autocorrupted
# (e.g. "10-Sep" for SEPT10, "2-Mar" for MARCH2); skip those for the
# illustrative plot so panel titles are readable gene symbols.
is_clean_symbol <- function(x) !grepl("^[0-9]+-[A-Za-z]{3}$", x)

status <- trans$gain_loss_status

# Rank candidates within a category by effect size (not just the first
# alphabetical match) so the example is visually representative of the
# category: strongest gain/loss posterior support, tightest phase
# agreement for "conserved", largest phase difference for "shifted".
pick_gene <- function(candidate_names, rank_by = NULL, decreasing = TRUE) {
  candidate_names <- candidate_names[!is.na(candidate_names)]
  candidate_names <- candidate_names[is_clean_symbol(candidate_names)]
  candidate_names <- intersect(candidate_names, rownames(baboon$expr_OMF))
  if (!is.null(rank_by)) {
    ord <- order(rank_by[candidate_names], decreasing = decreasing, na.last = NA)
    candidate_names <- candidate_names[ord]
  }
  candidate_names[1]
}

# CRY2 (a core clock gene) and ARNTL (BMAL1) are used directly for Gain
# and Phase-shifted, since both land in real, well-supported categories
# here (CRY2: p_gain = 0.87; ARNTL: phase-shifted, deltaPhi = -3.1h) and
# are far more recognizable than an arbitrary top-ranked gene.
gain_gene  <- "CRY2"
loss_gene  <- pick_gene(names(status)[status == "Loss"], rank_by = trans$p_loss)
cons_gene  <- pick_gene(names(phase$flag_cons)[phase$flag_cons %in% TRUE],
                        rank_by = abs(phase$deltaPhi.Est), decreasing = FALSE)
shift_gene <- "ARNTL"
undet_gene <- pick_gene(names(phase$flag_undetermined)[phase$flag_undetermined %in% TRUE],
                        rank_by = abs(phase$prob_shift - 0.5), decreasing = FALSE)

genes <- c(Gain = gain_gene, Loss = loss_gene, "Phase-conserved" = cons_gene,
           "Phase-shifted" = shift_gene, Undetermined = undet_gene)
cat("Selected genes:\n"); print(genes)

fit_panel <- function(gene, category, tissue, expr_mat) {
  y <- log2(expr_mat[gene, ] + 1)
  fit <- BayRC:::one_cosinor_OLS(tod = zt, y = y, period = 24)
  time_seq <- seq(0, 24, length.out = 200)
  omega <- 2 * pi / 24
  fitted_y <- fit$M + fit$A * cos(omega * (time_seq - fit$phase))

  df_pts <- data.frame(zt = zt, y = y)
  df_fit <- data.frame(zt = time_seq, y = fitted_y)

  ggplot() +
    geom_point(data = df_pts, aes(x = zt, y = y), color = "#377EB8", size = 2, alpha = 0.7) +
    geom_line(data = df_fit, aes(x = zt, y = y), color = "#E41A1C", linewidth = 1) +
    scale_x_continuous(breaks = seq(0, 24, 6)) +
    labs(title = sprintf("%s: %s (%s)", category, gene, tissue),
         subtitle = sprintf("peak = %.1fh, R² = %.2f, p = %.2e", fit$peak %% 24, fit$R2, fit$pvalue),
         x = "ZT (h)", y = "log2(expr + 1)") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 8, color = "#E41A1C"))
}

plots <- list()
for (cat in names(genes)) {
  g <- genes[[cat]]
  plots[[paste0(cat, "_OMF")]] <- fit_panel(g, cat, "OMF", baboon$expr_OMF)
  plots[[paste0(cat, "_THR")]] <- fit_panel(g, cat, "THR", baboon$expr_THR)
}

png("man/figures/cosinor_examples.png", width = 1500, height = 2200, res = 200, type = "cairo")
gridExtra::grid.arrange(grobs = plots, ncol = 2,
  top = grid::textGrob("Cosinor fit examples by transition/phase category (baboon OMF vs THR)",
                        gp = grid::gpar(fontsize = 12, fontface = "bold")))
dev.off()
cat("saved man/figures/cosinor_examples.png\n")
