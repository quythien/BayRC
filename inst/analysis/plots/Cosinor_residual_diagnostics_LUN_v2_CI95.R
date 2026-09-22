# ============================================================
# Cosinor Residual Diagnostics - Baboon Tissues, with 95% bands on the QQ plots
# Tissues: LUN, PUT, SUN
# Overlapping top genes are aligned in the same column across rows.
# Residuals are taken from the OLS cosinor fitted values. Writes the QQ and
# residual-by-time PDFs and one genome-wide diagnostics CSV per tissue.
# ============================================================

rm(list = ls())

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
bayrc.needs.summary <- FALSE
source(file.path(analysis.dir, "config.R"))

outdir  <- BAYRC_FIGURE_DIR
N_TOP   <- 6
TISSUES <- c("LUN", "PUT", "SUN")

load(file.path(BAYRC_GTEX_DIR, "data", "CAMO.bab.hum.RData"))
source(file.path(BAYRC_PIPELINE_DIR, "one_cosinor_OLS_new.R"))

library(dplyr)
library(ggplot2)
library(gridExtra)

# ── Helpers ───────────────────────────────────────────────────────────────────
convert_zt <- function(tod) ifelse(tod > 18, tod - 24, tod)

fmt_pval <- function(p) {
  if (round(p, 2) == 0) formatC(p, format = "e", digits = 2) else as.character(round(p, 2))
}

# QQ plot with a pointwise 95% band: the i-th order statistic's CDF value is
# Beta(i, n-i+1), mapped to sample space through the quartile QQ line.
make_qq_plot <- function(resid, gene_label, R2, pval) {
  n     <- length(resid)
  probs <- ppoints(n)           # (1:n - 0.5) / n

  # Theoretical normal quantiles (x-axis of QQ plot)
  theoretical <- qnorm(probs)

  # Pointwise 95% CI for each order statistic via Beta distribution
  i_seq   <- seq_len(n)
  p_lower <- qbeta(0.025, i_seq, n - i_seq + 1)
  p_upper <- qbeta(0.975, i_seq, n - i_seq + 1)

  # Normal quantile limits of the CI
  z_lower <- qnorm(p_lower)
  z_upper <- qnorm(p_upper)

  # QQ line parameters (same algorithm as stat_qq_line)
  qq_slope     <- diff(quantile(resid, c(0.25, 0.75))) / diff(qnorm(c(0.25, 0.75)))
  qq_intercept <- quantile(resid, 0.25) - qq_slope * qnorm(0.25)

  # Scale CI from theoretical to sample space
  band_df <- data.frame(
    x      = theoretical,
    y_low  = qq_intercept + qq_slope * z_lower,
    y_high = qq_intercept + qq_slope * z_upper
  )

  ggplot(data.frame(sample = resid), aes(sample = sample)) +
    geom_ribbon(data = band_df,
                aes(x = x, ymin = y_low, ymax = y_high),
                fill = "grey70", alpha = 0.4, inherit.aes = FALSE) +
    stat_qq(color = "steelblue", size = 1.8, alpha = 0.7) +
    stat_qq_line(color = "red", linetype = "dashed", linewidth = 0.8) +
    labs(title = sprintf("%s\ncosinor R2 = %.3f\np = %s", gene_label, R2, fmt_pval(pval)),
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
          axis.title = element_text(size = 8),
          plot.margin = margin(4, 10, 4, 4))
}

make_resid_tod_plot <- function(resid, tod, gene_label) {
  df <- data.frame(TOD = tod, Residual = resid)
  ggplot(df, aes(x = TOD, y = Residual)) +
    geom_point(alpha = 0.6, color = "steelblue", size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = gene_label, x = "Time of Day (ZT)", y = "Residual") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
          axis.title = element_text(size = 8),
          plot.margin = margin(4, 10, 4, 4))
}

empty_plot <- function() {
  ggplot() + theme_void()
}

# ── Fit cosinor + compute residuals for all genes in one tissue ───────────────
fit_all_cosinor <- function(expr_mat, tod) {
  n <- length(tod)
  cat("  Fitting", nrow(expr_mat), "genes...\n")
  results <- lapply(seq_len(nrow(expr_mat)), function(i) {
    if (i %% 500 == 0) cat("  ...", i, "/", nrow(expr_mat), "\n")
    y   <- as.numeric(expr_mat[i, ])
    fit <- tryCatch(
      one_cosinor_OLS(tod = tod, y = y, alpha = 0.05, period = 24, CI = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    fitted <- fit$M$est + fit$A$est * cos(2 * pi * tod / 24 + fit$phase$est)
    resid  <- y - fitted
    R2     <- fit$test$R2
    Fstat  <- R2 * (n - 3) / (2 * (1 - R2))
    pval   <- pf(Fstat, 2, n - 3, lower.tail = FALSE)
    list(R2 = R2, pval = pval, fitted = fitted, resid = resid)
  })
  names(results) <- rownames(expr_mat)
  results
}

# ── Main loop: fit all tissues, collect top genes and plots ───────────────────
tissue_fits  <- list()   # tissue -> fits list
tissue_top   <- list()   # tissue -> top N summary df
tissue_tod   <- list()   # tissue -> tod vector

for (tissue in TISSUES) {
  cat("\n============================================================\n")
  cat("Processing tissue:", tissue, "\n")
  cat("============================================================\n")

  bab_raw  <- baboon_withTOD$baboon[[tissue]]
  bab_cols <- grep(paste0(tissue, "\\.ZT"), colnames(bab_raw))
  bab_mat  <- as.matrix(bab_raw[, bab_cols])
  tod      <- convert_zt(as.numeric(sub(paste0(tissue, "\\.ZT"), "", colnames(bab_mat))))
  rownames(bab_mat) <- bab_raw$Symbol
  rownames(bab_mat)[rownames(bab_mat) == "ARNTL"] <- "BMAL1"

  # TPM + log2
  tpm    <- sweep(bab_mat, 2, colSums(bab_mat), FUN = "/") * 1e6
  logmat <- log2(tpm + 1)

  fits <- fit_all_cosinor(logmat, tod)

  summary_df <- data.frame(
    symbol = names(fits),
    R2     = sapply(fits, function(x) if (is.null(x)) NA_real_ else x$R2),
    pval   = sapply(fits, function(x) if (is.null(x)) NA_real_ else x$pval),
    stringsAsFactors = FALSE, row.names = NULL
  )
  summary_df <- summary_df[!is.na(summary_df$pval), ]
  summary_df <- summary_df[order(summary_df$pval), ]
  top <- head(summary_df, N_TOP)

  cat("\nTop", N_TOP, "genes -", tissue, ":\n")
  print(top[, c("symbol", "R2", "pval")])

  tissue_fits[[tissue]] <- fits
  tissue_top[[tissue]]  <- top
  tissue_tod[[tissue]]  <- tod
}

# ── Genome-wide residual diagnostics, one table per tissue ────────────────────
# shapiro.pval tests normality of the residuals; hetero.lm.pval regresses the
# squared residuals on the fitted values, small when the spread grows with level
for (tissue in TISSUES) {
  fits <- tissue_fits[[tissue]]
  fits <- fits[!sapply(fits, is.null)]

  diag_df <- data.frame(
    symbol         = names(fits),
    R2             = sapply(fits, function(x) x$R2),
    pvalue         = sapply(fits, function(x) x$pval),
    shapiro.pval   = sapply(fits, function(x)
                       tryCatch(shapiro.test(x$resid)$p.value,
                                error = function(e) NA_real_)),
    hetero.lm.pval = sapply(fits, function(x)
                       tryCatch(summary(lm(x$resid^2 ~ x$fitted))$coefficients[2, 4],
                                error = function(e) NA_real_)),
    mean.resid     = sapply(fits, function(x) mean(x$resid)),
    sd.resid       = sapply(fits, function(x) sd(x$resid)),
    stringsAsFactors = FALSE, row.names = NULL
  )
  diag_df$qvalue <- p.adjust(diag_df$pvalue, "BH")
  diag_df <- diag_df[order(diag_df$pvalue), ]

  f <- file.path(outdir, paste0("Cosinor_residual_diagnostics_", tissue, ".csv"))
  write.csv(diag_df, f, row.names = FALSE)
  cat(sprintf("%s: %d genes, normal at 0.05 %.1f%%, homoscedastic at 0.05 %.1f%%\n",
              tissue, nrow(diag_df),
              100 * mean(diag_df$shapiro.pval > 0.05, na.rm = TRUE),
              100 * mean(diag_df$hetero.lm.pval > 0.05, na.rm = TRUE)))
  cat("  wrote", f, "\n")
}

# ── Build master column order: N_TOP columns, overlap genes aligned first ─────
# Overlap = appears in top N of >=2 tissues → placed in same column across rows.
# Remaining N_TOP - K columns filled per-row with tissue-specific genes.
# Grid is always N_TOP columns wide; blank panel where a tissue lacks that gene.
all_top_genes <- lapply(tissue_top, function(df) df$symbol)
gene_freq     <- sort(table(unlist(all_top_genes)), decreasing = TRUE)

overlap_genes <- names(gene_freq[gene_freq >= 2])
avg_rank <- sapply(overlap_genes, function(g) {
  ranks <- sapply(TISSUES, function(t) {
    idx <- which(tissue_top[[t]]$symbol == g)
    if (length(idx) == 0) NA_real_ else idx
  })
  mean(ranks, na.rm = TRUE)
})
overlap_genes <- overlap_genes[order(avg_rank[overlap_genes])]

# Per-tissue column layout: [overlap genes] + [tissue-specific to fill N_TOP]
tissue_cols <- lapply(TISSUES, function(tissue) {
  specific <- setdiff(all_top_genes[[tissue]], overlap_genes)
  n_fill   <- N_TOP - length(overlap_genes)
  c(overlap_genes, head(specific, n_fill))
})
names(tissue_cols) <- TISSUES

n_cols <- N_TOP
cat("\nOverlap genes:", if (length(overlap_genes) == 0) "none" else paste(overlap_genes, collapse = ", "), "\n")
for (tissue in TISSUES)
  cat("  ", tissue, "columns:", paste(tissue_cols[[tissue]], collapse = ", "), "\n")

# ── Build grob lists (tissue row x gene column) ───────────────────────────────
qq_grobs <- list()
rt_grobs <- list()

for (tissue in TISSUES) {
  cols  <- tissue_cols[[tissue]]
  tod   <- tissue_tod[[tissue]]
  fits  <- tissue_fits[[tissue]]

  for (col_idx in seq_len(n_cols)) {
    key  <- paste(tissue, col_idx, sep = "_")
    gene <- if (col_idx <= length(cols)) cols[col_idx] else NA

    if (!is.na(gene) && gene %in% tissue_top[[tissue]]$symbol) {
      row   <- tissue_top[[tissue]][tissue_top[[tissue]]$symbol == gene, ]
      resid <- fits[[gene]]$resid
      label <- paste0(gene, "\n(", tissue, ")")
      qq_grobs[[key]] <- ggplotGrob(make_qq_plot(resid, label, row$R2, row$pval))
      rt_grobs[[key]] <- ggplotGrob(make_resid_tod_plot(resid, tod, label))
    } else {
      qq_grobs[[key]] <- ggplotGrob(empty_plot())
      rt_grobs[[key]] <- ggplotGrob(empty_plot())
    }
  }
}

# Row-major order
ordered_keys <- as.vector(t(outer(TISSUES, seq_len(n_cols), paste, sep = "_")))
qq_ordered   <- qq_grobs[ordered_keys]
rt_ordered   <- rt_grobs[ordered_keys]

# ── Save PDFs ─────────────────────────────────────────────────────────────────
# 1.75 in per panel keeps axis labels at 5 pt or more on a 6.5 in text block
pdf_width  <- n_cols * 1.75
pdf_height <- length(TISSUES) * 1.75

qq_file <- file.path(outdir, "Cosinor_QQ_Baboon_LUN_PUT_SUN_CI95.pdf")
pdf(qq_file, width = pdf_width, height = pdf_height)
grid.arrange(grobs = qq_ordered,
             nrow = length(TISSUES), ncol = n_cols,
             top = "Normality Check (QQ Plots, 95% CI band) - Baboon LUN / PUT / SUN - Top Rhythmic Genes")
dev.off()
cat("Saved:", basename(qq_file), "\n")

rt_file <- file.path(outdir, "Cosinor_ResidTOD_Baboon_LUN_PUT_SUN_CI95.pdf")
pdf(rt_file, width = pdf_width, height = pdf_height)
grid.arrange(grobs = rt_ordered,
             nrow = length(TISSUES), ncol = n_cols,
             top = "Residuals vs Time of Day - Baboon LUN / PUT / SUN - Top Rhythmic Genes")
dev.off()
cat("Saved:", basename(rt_file), "\n")

cat("\nAll done. Output dir:", outdir, "\n")
