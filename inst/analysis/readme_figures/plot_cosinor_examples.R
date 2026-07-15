# Builds man/figures/cosinor_examples.png: raw expression scatter + OLS
# cosinor fit (BayRC's one_cosinor_OLS()/Cosinor_fit(), the classical
# frequentist baseline) for one example gene per transition/phase
# category (Gain, Loss, Phase-conserved, Phase-shifted), baboon OMF vs
# THR, on the bundled quickstart-scale posterior (2,500 iterations, from
# inst/analysis/quickstart_baboon_OMF_THR.R). Assumes mcmc_OMF, mcmc_THR,
# bf_OMF, bf_THR, trans, phase are already built in the environment (see
# quickstart_baboon_OMF_THR.R for the exact steps).
#
# trans$gain_loss_status and phase$flag_* are already named by gene
# symbol (inherited from rownames(mcmc_OMF$rho)/rownames(mcmc_OMF$phi),
# set by match_symbols()); do not try to re-name them from mcmc_OMF$gname,
# which does not exist on a match_symbols()-processed object and is NULL.
#
# p/q come from Cosinor_fit() run genome-wide on each tissue (BH-adjusted
# across all 5,066 genes), not from a single-gene fit, so q is a real,
# not fabricated, multiple-testing-adjusted value. BF comes from the
# RJMCMC posterior (bf_OMF/bf_THR), a separate estimate from a different
# model; the classical fit is for visualization only.

suppressMessages(library(ggplot2))

baboon <- readRDS(system.file("extdata", "baboon_OMF_THR_GSE98965.rds", package = "BayRC"))
zt <- baboon$zt

is_clean_symbol <- function(x) !grepl("^[0-9]+-[A-Za-z]{3}$", x)

status <- trans$gain_loss_status
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
# here and are far more recognizable than an arbitrary top-ranked gene.
gain_gene  <- "CRY2"
loss_gene  <- pick_gene(names(status)[status == "Loss"], rank_by = trans$p_loss)
cons_gene  <- pick_gene(names(phase$flag_cons)[phase$flag_cons %in% TRUE],
                        rank_by = abs(phase$deltaPhi.Est), decreasing = FALSE)
shift_gene <- "ARNTL"

genes <- c(Gain = gain_gene, Loss = loss_gene, "Phase-conserved" = cons_gene,
           "Phase-shifted" = shift_gene)
cat("Selected genes:\n"); print(genes)

# Genome-wide classical cosinor fit, both tissues, for real BH q-values.
x_OMF <- list(data = as.data.frame(log2(baboon$expr_OMF + 1)), time = zt, gname = baboon$gene_symbol)
x_THR <- list(data = as.data.frame(log2(baboon$expr_THR + 1)), time = zt, gname = baboon$gene_symbol)
rhythm_OMF <- Cosinor_fit(x_OMF)$rhythm
rhythm_THR <- Cosinor_fit(x_THR)$rhythm
rownames(rhythm_OMF) <- rhythm_OMF$gname
rownames(rhythm_THR) <- rhythm_THR$gname

fmt_bf <- function(bf) {
  if (bf >= 1e6) "> 1e6" else if (bf < 1) sprintf("= %.3f", bf) else sprintf("= %.0f", bf)
}

time_seq <- seq(0, 24, length.out = 200)
omega <- 2 * pi / 24

pts_list <- list(); fit_list <- list(); stat_list <- list(); peak_list <- list()
for (cat in names(genes)) {
  g <- genes[[cat]]
  for (tissue in c("OMF", "THR")) {
    expr_mat <- if (tissue == "OMF") baboon$expr_OMF else baboon$expr_THR
    rhythm   <- if (tissue == "OMF") rhythm_OMF else rhythm_THR
    bf_tab   <- if (tissue == "OMF") bf_OMF else bf_THR

    y <- log2(expr_mat[g, ] + 1)
    r <- rhythm[g, ]
    fitted_y <- r$M + r$A * cos(omega * (time_seq - r$phase))
    bf <- bf_tab[g, "BayesF"]
    peak <- r$peak %% 24

    row_label <- sprintf("%s: %s", cat, g)
    pts_list[[paste(cat, tissue)]] <- data.frame(row_label = row_label, tissue = tissue, zt = zt, y = y)
    fit_list[[paste(cat, tissue)]] <- data.frame(row_label = row_label, tissue = tissue, zt = time_seq, y = fitted_y)
    peak_list[[paste(cat, tissue)]] <- data.frame(row_label = row_label, tissue = tissue, peak = peak)
    stat_list[[paste(cat, tissue)]] <- data.frame(
      row_label = row_label, tissue = tissue,
      label = sprintf("peak = %.1fh\np = %.2e, q = %.2e\nBF %s",
                       peak, r$pvalue, r$qvalue, fmt_bf(bf))
    )
  }
}

row_order <- sprintf("%s: %s", names(genes), genes)
pts_df  <- do.call(rbind, pts_list)
fit_df  <- do.call(rbind, fit_list)
stat_df <- do.call(rbind, stat_list)
peak_df <- do.call(rbind, peak_list)
pts_df$row_label  <- factor(pts_df$row_label,  levels = row_order)
fit_df$row_label  <- factor(fit_df$row_label,  levels = row_order)
stat_df$row_label <- factor(stat_df$row_label, levels = row_order)
peak_df$row_label <- factor(peak_df$row_label, levels = row_order)
pts_df$tissue  <- factor(pts_df$tissue,  levels = c("OMF", "THR"))
fit_df$tissue  <- factor(fit_df$tissue,  levels = c("OMF", "THR"))
stat_df$tissue <- factor(stat_df$tissue, levels = c("OMF", "THR"))
peak_df$tissue <- factor(peak_df$tissue, levels = c("OMF", "THR"))

p <- ggplot() +
  geom_point(data = pts_df, aes(x = zt, y = y), color = "#377EB8", size = 2, alpha = 0.7) +
  geom_line(data = fit_df, aes(x = zt, y = y), color = "#E41A1C", linewidth = 1) +
  geom_vline(data = peak_df, aes(xintercept = peak),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_text(data = stat_df, aes(x = 12, y = Inf, label = label),
            vjust = 1.15, size = 2.8, color = "black", lineheight = 0.95) +
  scale_x_continuous(breaks = seq(0, 24, 6)) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.42))) +
  facet_grid(rows = vars(row_label), cols = vars(tissue), scales = "free_y", switch = "y") +
  labs(title = "Cosinor fit examples by transition/phase category, baboon OMF vs THR",
       x = "ZT (h)", y = "log2(expr + 1)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        strip.text = element_text(face = "bold", size = 10),
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        strip.placement = "outside",
        strip.background.x = element_rect(fill = "grey92", color = NA),
        strip.background.y = element_blank(),
        panel.spacing.y = unit(1.1, "lines"),
        panel.spacing.x = unit(1.3, "lines"),
        panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.4),
        panel.grid.minor = element_blank())

png("man/figures/cosinor_examples.png", width = 1750, height = 2200, res = 200, type = "cairo")
print(p)
dev.off()
cat("saved man/figures/cosinor_examples.png\n")
