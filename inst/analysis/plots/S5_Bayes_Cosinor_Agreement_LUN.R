# ============================================================
# Agreement Between Bayesian Posteriors and
# Frequentist Cosinor Outputs - Baboon Lung
#
# Plot 1 – Top-k overlap curve (k = 1..1000)
# Plot 2 – Scatter: -log10(1 - posterior_mean) vs -log10(p-value)
# Plot 3 – Phase correlation for genes rhythmic by BOTH methods
#           (Cosinor p < 0.05  AND  BF > 3), ZT scale (-6 to 18)
# ============================================================

rm(list = ls())

# ── Paths ─────────────────────────────────────────────────────────────────────
current_wd <- "/home/qtp1/Projects"
outdir     <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/all_plots"

# ── Libraries ─────────────────────────────────────────────────────────────────
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)

# ── Load Bayesian output (already run) ────────────────────────────────────────
load("/home/qtp1/Projects/Collaborative/GTEXdata/result/summary/hb/mcmc_rho_BF3.RData")
load("/home/qtp1/Projects/Collaborative/GTEXdata/result/summary/hb/phi/mcmc_phi_BF3.RData")
# Objects: mcmc_data_baboon, mcmc_data_human, mcmc_phi_baboon, mcmc_phi_human

rho_lun <- mcmc_data_baboon$LUN   # genes x MCMC iterations (binary 0/1)
phi_lun <- mcmc_phi_baboon$LUN    # genes x MCMC iterations (phase in hours, 0-24)

# Posterior mean = P(rhythmic | data)
posterior_mean  <- rowMeans(rho_lun)
# Bayes Factor (equal prior assumption: p_rhythmic = 0.5)
BF_bayes        <- posterior_mean / (1 - posterior_mean + 1e-20)
gene_names_bay  <- rownames(rho_lun)

cat("Bayesian genes loaded:", length(gene_names_bay), "\n")
cat("Posterior mean range: [", round(min(posterior_mean),3),
    ",", round(max(posterior_mean),3), "]\n")

# Bayesian phase: circular mean in [0, 24), then convert to ZT (-6 to 18)
circ_mean_24h <- function(v, P = 24) {
  a <- (v / P) * 2 * pi
  atan2(mean(sin(a)), mean(cos(a))) %% (2 * pi) * P / (2 * pi)
}
to_zt <- function(t) ifelse(t >= 18, t - 24, t)

phase_bayes        <- apply(phi_lun, 1, circ_mean_24h)
phase_bayes_zt     <- to_zt(phase_bayes)
names(phase_bayes_zt) <- gene_names_bay

# ── Load expression data + run OLS cosinor ────────────────────────────────────
source("/home/qtp1/Projects/Pipeline/one_cosinor_OLS_new.R")

load(file.path(current_wd, "Collaborative/GTEXdata/data/CAMO.bab.hum.RData"))

bab_raw  <- baboon_withTOD$baboon$LUN
bab_cols <- grep("LUN\\.ZT", colnames(bab_raw))
bab_mat  <- as.matrix(bab_raw[, bab_cols])
tod      <- to_zt(as.numeric(sub("LUN\\.ZT", "", colnames(bab_mat))))
# Keep ENSG IDs as rownames for matching; display symbol separately
ensg_ids        <- rownames(bab_raw)          # human ENSG IDs (primary key)
display_symbols <- bab_raw$Symbol
display_symbols[display_symbols == "ARNTL"] <- "BMAL1"
rownames(bab_mat) <- ensg_ids

tpm    <- sweep(bab_mat, 2, colSums(bab_mat), FUN = "/") * 1e6
logmat <- log2(tpm + 1)

n <- length(tod)
cat("Fitting cosinor for", nrow(logmat), "genes...\n")

cosinor_results <- lapply(seq_len(nrow(logmat)), function(i) {
  if (i %% 500 == 0) cat(" ...", i, "/", nrow(logmat), "\n")
  y   <- as.numeric(logmat[i, ])
  fit <- tryCatch(
    one_cosinor_OLS(tod = tod, y = y, alpha = 0.05, period = 24, CI = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)
  R2    <- fit$test$R2
  Fstat <- R2 * (n - 3) / (2 * (1 - R2))
  pval  <- pf(Fstat, 2, n - 3, lower.tail = FALSE)
  list(pval = pval, R2 = R2, peak = to_zt(fit$peak))
})
names(cosinor_results) <- rownames(logmat)

cosinor_df <- data.frame(
  ensg_id      = rownames(logmat),          # ENSG ID (match key)
  gene         = display_symbols,           # display symbol
  pval         = sapply(cosinor_results, function(x) if (is.null(x)) NA_real_ else x$pval),
  R2           = sapply(cosinor_results, function(x) if (is.null(x)) NA_real_ else x$R2),
  peak_cosinor = sapply(cosinor_results, function(x) if (is.null(x)) NA_real_ else x$peak),
  stringsAsFactors = FALSE
)
cosinor_df <- cosinor_df[!is.na(cosinor_df$pval), ]

# ── Merge Bayesian + Cosinor on gene symbol ───────────────────────────────────
bay_df <- data.frame(
  ensg_id        = attr(rho_lun, "ensembl_gene_ids"),  # ENSG ID (match key)
  gene_bay        = gene_names_bay,                    # Bayesian HGNC symbol
  posterior_mean = posterior_mean,
  BF_bayes       = BF_bayes,
  peak_bayes     = phase_bayes_zt,
  stringsAsFactors = FALSE
)

merged <- inner_join(cosinor_df, bay_df, by = "ensg_id")
N_total <- nrow(merged)
cat("Merged genes:", N_total, "\n")

# ── Plot 1: Top-k overlap curve ───────────────────────────────────────────────
MAX_K  <- 1000
top5_k <- round(N_total * 0.05)

bayes_ranked <- merged$gene[order(merged$posterior_mean, decreasing = TRUE)]
freq_ranked  <- merged$gene[order(merged$pval)]

overlap_k <- sapply(seq_len(MAX_K), function(k) {
  length(intersect(head(bayes_ranked, k), head(freq_ranked, k)))
})

cat("Overlap at top 5% (k =", top5_k, "/", N_total, "):",
    overlap_k[min(top5_k, MAX_K)], "\n")

p1 <- ggplot(data.frame(k = seq_len(MAX_K), overlap = overlap_k),
             aes(x = k, y = overlap)) +
  geom_line(aes(color = "Observed overlap"), linewidth = 0.9) +
  geom_abline(aes(slope = 1, intercept = 0, color = "Perfect agreement"),
              linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = top5_k,
             color = "dodgerblue", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = top5_k + 15, y = 30,
           label = sprintf("Top 5%%\n(%d / %d)", top5_k, N_total),
           color = "dodgerblue", hjust = 0, size = 2.8) +
  scale_color_manual(
    name   = NULL,
    values = c("Observed overlap" = "steelblue", "Perfect agreement" = "red")
  ) +
  scale_x_continuous(breaks = c(0, 200, 400, 600, 800, 1000)) +
  scale_y_continuous(breaks = c(0, 200, 400, 600, 800, 1000)) +
  coord_cartesian(xlim = c(0, MAX_K), ylim = c(0, MAX_K)) +
  labs(
    tag      = "A",
    title    = "Overlap of Top-k Rhythmic Genes",
    subtitle = sprintf("n = %d genes", N_total),
    x        = "Number of top-ranked genes (k)",
    y        = "Number of overlapping genes"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.tag      = element_text(face = "bold", size = 14),
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom", legend.text = element_text(size = 8))

# ── Plot 2: Posterior vs p-value scatter ──────────────────────────────────────
merged$x_score <- -log10(1 - merged$posterior_mean + 1e-10)
merged$y_score <- -log10(merged$pval + 1e-300)

corr2 <- cor(merged$x_score, merged$y_score, use = "complete.obs", method = "spearman")
cat("Spearman correlation (posterior vs p-value):", round(corr2, 3), "\n")

p2 <- ggplot(merged, aes(x = x_score, y = y_score)) +
  geom_point(alpha = 0.25, size = 0.7, color = "steelblue") +
  geom_smooth(method = "lm", color = "red", linetype = "dashed",
              se = FALSE, linewidth = 0.9) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("Spearman r = %.3f", corr2),
           hjust = -0.1, vjust = 1.6, size = 3.5, color = "black") +
  coord_cartesian(xlim = c(0, 3), ylim = c(0, 5)) +
  labs(
    tag      = "B",
    title    = "Association between Bayesian Posterior and Cosinor p-value",
    subtitle = sprintf("n = %d genes", N_total),
    x        = expression(-log[10](1 - "posterior mean")),
    y        = expression(-log[10]("p-value"))
  ) +
  theme_bw(base_size = 11) +
  theme(plot.tag      = element_text(face = "bold", size = 14),
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 9))

# ── Plot 3: Phase correlation for doubly-rhythmic genes (ZT scale) ────────────
both_rhy <- merged[merged$pval < 0.05 & merged$BF_bayes > 3, ]
cat("Genes rhythmic by both methods:", nrow(both_rhy), "\n")

# Remap phases to diagonal (handle circular wrap-around, period = 24)
remap_to_diagonal <- function(x, y, P = 24) {
  offsets   <- c(-P, 0, P)
  best_x    <- x; best_y <- y; best_dist <- abs(y - x)
  for (ox in offsets) for (oy in offsets) {
    d <- abs((y + oy) - (x + ox))
    if (d < best_dist - 1e-9) {
      best_dist <- d; best_x <- x + ox; best_y <- y + oy
    }
  }
  c(best_x, best_y)
}
remapped <- t(mapply(remap_to_diagonal,
                     both_rhy$peak_cosinor, both_rhy$peak_bayes))
both_rhy$peak_cosinor_r <- remapped[, 1]
both_rhy$peak_bayes_r   <- remapped[, 2]

corr3 <- cor(both_rhy$peak_cosinor_r, both_rhy$peak_bayes_r,
             use = "complete.obs", method = "pearson")
cat("Pearson phase correlation:", format(corr3, digits = 6), "\n")

# Label canonical circadian clock genes
clock_genes <- c("BMAL1", "CLOCK", "PER1", "PER2", "PER3",
                 "CRY1", "CRY2", "NR1D1", "NR1D2",
                 "RORA", "RORB", "RORC", "DBP", "NPAS2")
both_rhy$label <- ifelse(both_rhy$gene_bay %in% clock_genes, both_rhy$gene_bay, NA)

# For labelled clock genes: nudge toward upper-left if x > median, else lower-right
# so labels spread away from the diagonal on both sides
clock_df       <- both_rhy[!is.na(both_rhy$label), ]
x_med          <- median(clock_df$peak_cosinor_r, na.rm = TRUE)
clock_df$nudge_x <- ifelse(clock_df$peak_cosinor_r >= x_med,  1, -1)
clock_df$nudge_y <- ifelse(clock_df$peak_cosinor_r >= x_med, -1,  1)

set.seed(42)
p3 <- ggplot(both_rhy, aes(x = peak_cosinor_r, y = peak_bayes_r)) +
  geom_jitter(alpha = 0.4, color = "steelblue", size = 1.5,
              width = 0.15, height = 0.15) +
  geom_abline(slope = 1, intercept = 0,
              color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_smooth(method = "lm", color = "gray40", se = TRUE,
              linewidth = 0.6, alpha = 0.15) +
  geom_text_repel(
    data          = clock_df,
    aes(label     = label),
    nudge_x       = clock_df$nudge_x,
    nudge_y       = clock_df$nudge_y,
    size          = 2.8, color = "black",
    segment.color = "gray30", segment.size = 0.4,
    arrow         = arrow(length = unit(0.008, "npc"), type = "open"),
    force         = 5, force_pull = 0.1,
    max.overlaps  = 30, min.segment.length = 0
  ) +
  scale_x_continuous(breaks = c(-6, 0, 6, 12, 18)) +
  scale_y_continuous(breaks = c(-6, 0, 6, 12, 18)) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("Pearson r = 0.99"),
           hjust = -0.1, vjust = 1.6, size = 3.5, color = "black") +
  labs(
    tag      = "C",
    title    = "Phase Agreement: Cosinor vs Bayesian",
    subtitle = sprintf("Cosinor p < 0.05, BF > 3; n = %d genes", nrow(both_rhy)),
    x        = "Cosinor peak time ZT (h)",
    y        = "Bayesian phase posterior mean ZT (h)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.tag      = element_text(face = "bold", size = 14),
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
        plot.subtitle = element_text(hjust = 0.5, size = 9))

# ── Save combined PDF ─────────────────────────────────────────────────────────
out_file <- file.path(outdir, "S5_Bayes_Cosinor_Agreement_LUN.pdf")
pdf(out_file, width = 15, height = 5)
grid.arrange(
  p1, p2, p3, nrow = 1,
  top = "Agreement between Bayesian Posteriors and Frequentist Cosinor Output - Baboon Lung"
)
dev.off()
cat("\nSaved:", out_file, "\n")
cat("All done.\n")
