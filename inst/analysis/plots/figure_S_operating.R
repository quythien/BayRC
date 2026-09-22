## Supplementary figure: false discovery rate, power and type I error of BFDR
## across signal strength and sample size.
##
## Reads the swept calibration runs that pipeline/bfdr_calibration.R writes as
## bfdr_calibration_A<A>_s1_n<n>_seed<s>.rds, scores each against its known
## truth, and draws three panels:
##   A  realised against nominal FDR, one line per A/sigma, faceted by n
##   B  power against A/sigma at each nominal level, one line per n
##   C  type I error against A/sigma, one line per n
## Realised FDR pools false calls over calls across replicates; power and type I
## error are replicate means. Also writes the pooled table beside the figure.

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
bayrc.needs.summary <- FALSE
source(file.path(analysis.dir, "config.R"))
source(file.path(analysis.dir, "plots", "theme_bayrc.R"))
suppressPackageStartupMessages({ library(BayRC); library(ggplot2) })

cal.dir <- file.path(BAYRC_OUTPUT_DIR, "calibration")
fs <- list.files(cal.dir, "^bfdr_calibration_A[0-9.]+_s1_n[0-9]+_seed[0-9]+[.]rds$",
                 full.names = TRUE)
if (!length(fs)) stop("no swept calibration runs under ", cal.dir)
alphas <- c(0.05, 0.10, 0.25)

key <- regmatches(basename(fs), regexec("A([0-9.]+)_s1_n([0-9]+)_seed([0-9]+)", basename(fs)))
runs <- data.frame(file = fs,
                   A = as.numeric(sapply(key, `[`, 2)),
                   n = as.integer(sapply(key, `[`, 3)))

rows <- list()
for (i in seq_len(nrow(runs))) {
  x <- readRDS(runs$file[i]); tr <- x$truth
  for (a in alphas) {
    call <- bfdr_from_posterior(x$post, alpha = a)$rhythmic_genes
    rows[[length(rows) + 1]] <- data.frame(
      A = runs$A[i], n = runs$n[i], alpha = a,
      calls = sum(call), false = sum(call & tr == 0),
      power = mean(call[tr == 1]), type1 = mean(call[tr == 0]))
  }
}
d <- do.call(rbind, rows)
pooled <- do.call(rbind, lapply(split(d, list(d$A, d$n, d$alpha), drop = TRUE),
  function(x) data.frame(A = x$A[1], n = x$n[1], alpha = x$alpha[1],
                         replicates = nrow(x),
                         fdr = sum(x$false) / max(1, sum(x$calls)),
                         power = mean(x$power), type1 = mean(x$type1))))
pooled <- pooled[order(pooled$n, pooled$A, pooled$alpha), ]
write.csv(pooled, file.path(cal.dir, "bfdr_operating_characteristics.csv"), row.names = FALSE)

pooled$snr  <- factor(sprintf("A/σ = %g", pooled$A), levels = sprintf("A/σ = %g", sort(unique(pooled$A))))
pooled$nlab <- factor(sprintf("n = %d", pooled$n), levels = sprintf("n = %d", sort(unique(pooled$n))))
pooled$alab <- factor(sprintf("α = %g", pooled$alpha), levels = sprintf("α = %g", alphas))

## panel letters sit at the top left of each panel, as in Figures 2-6; signal
## strength is ordered so it takes the sequential ramp, sample size the levels
letter <- theme(plot.title = element_text(hjust = 0, face = "bold"),
                plot.title.position = "plot")
snr_cols <- setNames(bayrc_seq(length(levels(pooled$snr)) + 2)[-(1:2)], levels(pooled$snr))
n_cols   <- setNames(bayrc_levels[seq_along(levels(pooled$nlab))], levels(pooled$nlab))

pA <- ggplot(pooled, aes(alpha, fdr, colour = snr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = bayrc_ink) +
  geom_line() + geom_point(size = 1.6) +
  facet_wrap(~ nlab, nrow = 1) +
  scale_x_continuous(breaks = alphas) +
  scale_colour_manual(values = snr_cols) +
  labs(x = "Nominal FDR", y = "Realised FDR", colour = NULL, title = "A") +
  theme_bayrc() + letter

pB <- ggplot(pooled, aes(A, power, colour = nlab)) +
  geom_line() + geom_point(size = 1.6) +
  facet_wrap(~ alab, nrow = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = n_cols) +
  labs(x = "A/σ", y = "Power", colour = NULL, title = "B") +
  theme_bayrc() + letter

pC <- ggplot(pooled, aes(A, type1, colour = nlab)) +
  geom_hline(aes(yintercept = alpha), linetype = "dashed", colour = bayrc_ink) +
  geom_line() + geom_point(size = 1.6) +
  facet_wrap(~ alab, nrow = 1) +
  scale_colour_manual(values = n_cols) +
  labs(x = "A/σ", y = "Type I error", colour = NULL, title = "C") +
  theme_bayrc() + letter

out <- file.path(BAYRC_FIGURE_DIR, "Figure_S_operating.pdf")
cairo_pdf(out, width = 9.5, height = 9, family = bayrc_family)
gridExtra::grid.arrange(pA, pB, pC, ncol = 1)
invisible(dev.off())
cat("wrote", out, "from", nrow(runs), "runs\n")
print(format(pooled[, c("n", "A", "alpha", "fdr", "power", "type1")], digits = 3), row.names = FALSE)
