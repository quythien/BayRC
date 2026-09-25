## Supplementary figure: ranking accuracy, power, type I error and false
## discovery rate of BayRC across signal strength and sample size.
##
## Reads the swept calibration runs that pipeline/bfdr_calibration.R writes as
## bfdr_calibration_A<A>_s1_n<n>_seed<s>.rds, regenerates each simulated data set
## from its seed to fit an ordinary least-squares cosinor alongside, and draws
## five panels:
##   A  AUC of the BayRC posterior and the cosinor F-test p-value, by A/sigma and n
##   B  BFDR power against A/sigma at each nominal level, one line per n
##   C  BFDR type I error against A/sigma, one line per n
##   D  realised against nominal FDR for BFDR, faceted by n, one line per A/sigma
##   E  the same for the cosinor F-test with Benjamini-Hochberg
## Realised FDR is the replicate mean of false calls over max(calls, 1), as FDR
## is defined; AUC, power and type I error are replicate means too. Every point
## carries a 95% interval over the available replicates, the mean plus and minus
## 1.96 standard errors. Also writes the pooled table beside the figure.

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
# Saved-summary mode redraws the figure without fitting or rescoring simulations.
# Usage: Rscript figure_S_operating.R --summaries <CSV directory> <figure directory>
args <- commandArgs(trailingOnly = TRUE)
summary_mode <- length(args) > 0 && args[1] == "--summaries"
if (summary_mode) {
  stopifnot(length(args) == 3)
  cal.dir <- normalizePath(args[2])
  BAYRC_FIGURE_DIR <- normalizePath(args[3])
} else {
  bayrc.needs.summary <- FALSE
  source(file.path(analysis.dir, "config.R"))
  cal.dir <- file.path(BAYRC_OUTPUT_DIR, "calibration")
}
source(file.path(analysis.dir, "plots", "theme_bayrc.R"))
suppressPackageStartupMessages(library(ggplot2))
if (!summary_mode) suppressPackageStartupMessages(library(BayRC))

if (!summary_mode) {
fs <- list.files(cal.dir, "^bfdr_calibration_A[0-9.]+_s1_n[0-9]+_seed[0-9]+[.]rds$",
                 full.names = TRUE)
if (!length(fs)) stop("no swept calibration runs under ", cal.dir)
alphas <- c(0.05, 0.10, 0.15, 0.20, 0.25)
P <- 24; omega <- 2 * pi / P

key <- regmatches(basename(fs), regexec("A([0-9.]+)_s1_n([0-9]+)_seed([0-9]+)", basename(fs)))
runs <- data.frame(file = fs,
                   A    = as.numeric(sapply(key, `[`, 2)),
                   n    = as.integer(sapply(key, `[`, 3)),
                   seed = as.integer(sapply(key, `[`, 4)))

## the swept generator of bfdr_calibration.R, in the same draw order
simulate <- function(G, n, A, seed) {
  tod <- seq(0, P, length.out = n + 1)[-(n + 1)]
  set.seed(1000 + seed)
  truth <- rbinom(G, 1, 0.2)
  runif(G, 0.3, 1.2)
  phase <- runif(G, 0, P)
  rnorm(G, 5, 1)
  Amp <- ifelse(truth == 1, A, 0)
  Y <- 5 + Amp * cos(omega * (matrix(tod, G, n, byrow = TRUE) - phase)) +
       matrix(rnorm(G * n, 0, 1), G, n)
  list(Y = Y, tod = tod, truth = truth)
}

## cosinor F-test p-value for every row of Y
cosinor_p <- function(Y, tod) {
  X    <- cbind(1, cos(omega * tod), sin(omega * tod))
  hat  <- X %*% solve(crossprod(X), t(X))
  rss1 <- rowSums((Y - Y %*% t(hat))^2)
  rss0 <- rowSums((Y - rowMeans(Y))^2)
  N    <- length(tod)
  pf(((rss0 - rss1) / 2) / (rss1 / (N - 3)), 2, N - 3, lower.tail = FALSE)
}

## area under the ROC curve from ranks, ties shared
auc <- function(score, truth) {
  r <- rank(score); n1 <- sum(truth == 1); n0 <- sum(truth == 0)
  (sum(r[truth == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

rows <- list(); aucs <- list()
for (i in seq_len(nrow(runs))) {
  x <- readRDS(runs$file[i]); tr <- x$truth
  s <- simulate(length(tr), runs$n[i], runs$A[i], runs$seed[i])
  if (!identical(as.integer(s$truth), as.integer(tr)))
    stop(basename(runs$file[i]), ": regenerated truth does not match the stored truth")
  aucs[[i]] <- data.frame(A = runs$A[i], n = runs$n[i],
                          BayRC = auc(x$post, tr),
                          Cosinor = auc(-cosinor_p(s$Y, s$tod), tr))
  q <- p.adjust(cosinor_p(s$Y, s$tod), "BH")
  for (a in alphas) {
    call <- bfdr_from_posterior(x$post, alpha = a)$rhythmic_genes
    bh   <- q <= a
    rows[[length(rows) + 1]] <- data.frame(
      A = runs$A[i], n = runs$n[i], alpha = a,
      calls = sum(call), false = sum(call & tr == 0),
      power = mean(call[tr == 1]), type1 = mean(call[tr == 0]),
      bh_calls = sum(bh), bh_false = sum(bh & tr == 0),
      bh_power = mean(bh[tr == 1]))
  }
}
d <- do.call(rbind, rows)
se <- function(v) if (length(v) > 1) sd(v) / sqrt(length(v)) else NA_real_
pooled <- do.call(rbind, lapply(split(d, list(d$A, d$n, d$alpha), drop = TRUE),
  function(x) data.frame(A = x$A[1], n = x$n[1], alpha = x$alpha[1],
                         replicates = nrow(x),
                         fdr = mean(x$false / pmax(x$calls, 1)),
                         fdr_se = se(x$false / pmax(x$calls, 1)),
                         power = mean(x$power), power_se = se(x$power),
                         type1 = mean(x$type1), type1_se = se(x$type1),
                         bh_fdr = mean(x$bh_false / pmax(x$bh_calls, 1)),
                         bh_fdr_se = se(x$bh_false / pmax(x$bh_calls, 1)),
                         bh_power = mean(x$bh_power))))
pooled <- pooled[order(pooled$n, pooled$A, pooled$alpha), ]
write.csv(pooled, file.path(cal.dir, "bfdr_operating_characteristics.csv"), row.names = FALSE)

a <- do.call(rbind, aucs)
auc_pooled <- aggregate(cbind(BayRC, Cosinor) ~ A + n, data = a, FUN = mean)
auc_se <- aggregate(cbind(BayRC, Cosinor) ~ A + n, data = a, FUN = se)
names(auc_se)[3:4] <- c("BayRC_se", "Cosinor_se")
auc_pooled <- merge(auc_pooled, auc_se, by = c("A", "n"))
auc_pooled <- auc_pooled[order(auc_pooled$n, auc_pooled$A), ]
write.csv(auc_pooled, file.path(cal.dir, "bayrc_cosinor_auc.csv"), row.names = FALSE)
 } else {
  pooled <- read.csv(file.path(cal.dir, "bfdr_operating_characteristics.csv"))
  auc_pooled <- read.csv(file.path(cal.dir, "bayrc_cosinor_auc.csv"))
  alphas <- sort(unique(pooled$alpha))
}
auc_long <- rbind(data.frame(auc_pooled[, c("A", "n")], method = "BayRC",
                             auc = auc_pooled$BayRC, se = auc_pooled$BayRC_se),
                  data.frame(auc_pooled[, c("A", "n")], method = "Cosinor",
                             auc = auc_pooled$Cosinor, se = auc_pooled$Cosinor_se))

pooled$snr  <- factor(sprintf("A/σ = %g", pooled$A), levels = sprintf("A/σ = %g", sort(unique(pooled$A))))
pooled$nlab <- factor(sprintf("n = %d", pooled$n), levels = sprintf("n = %d", sort(unique(pooled$n))))
pooled$alab <- factor(sprintf("α = %g", pooled$alpha), levels = sprintf("α = %g", alphas))
auc_long$nlab <- factor(sprintf("n = %d", auc_long$n), levels = levels(pooled$nlab))

## panel letters at the top left; signal strength on the sequential ramp,
## sample size on the level colours
letter <- theme(plot.title = element_text(hjust = 0, face = "bold"),
                plot.title.position = "plot")
snr_cols <- setNames(bayrc_seq(length(levels(pooled$snr)) + 2)[-(1:2)], levels(pooled$snr))
n_cols   <- setNames(bayrc_levels[seq_along(levels(pooled$nlab))], levels(pooled$nlab))
method_cols <- c(BayRC = bayrc_levels[1], Cosinor = bayrc_ink3)

## the two curves differ by at most 0.007, so they are offset along A/sigma to
## keep both visible; cosinor is drawn dashed with open points
nudge <- position_dodge(width = 0.13)
pA <- ggplot(auc_long, aes(A, auc, colour = method, linetype = method, shape = method)) +
  geom_errorbar(aes(ymin = auc - 1.96 * se, ymax = auc + 1.96 * se), width = 0.08,
                linetype = "solid", position = nudge) +
  geom_line(position = nudge) + geom_point(size = 1.8, stroke = 0.7, position = nudge) +
  facet_wrap(~ nlab, nrow = 1) +
  scale_y_continuous(limits = c(0.5, 1)) +
  scale_colour_manual(values = method_cols) +
  scale_linetype_manual(values = c(BayRC = "solid", Cosinor = "22")) +
  scale_shape_manual(values = c(BayRC = 16, Cosinor = 1)) +
  labs(x = "A/σ", y = "AUC", colour = NULL, linetype = NULL, shape = NULL, title = "A") +
  theme_bayrc() + letter

pB <- ggplot(pooled, aes(A, power, colour = nlab)) +
  geom_errorbar(aes(ymin = power - 1.96 * power_se, ymax = power + 1.96 * power_se), width = 0.08) +
  geom_line() + geom_point(size = 1.6) +
  facet_wrap(~ alab, nrow = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = n_cols) +
  labs(x = "A/σ", y = "Power", colour = NULL, title = "B") +
  theme_bayrc() + letter + theme(strip.text.x = element_text(size = 9))

pC <- ggplot(pooled, aes(A, type1, colour = nlab)) +
  geom_hline(aes(yintercept = alpha), linetype = "dashed", colour = bayrc_ink) +
  geom_errorbar(aes(ymin = type1 - 1.96 * type1_se, ymax = type1 + 1.96 * type1_se), width = 0.08) +
  geom_line() + geom_point(size = 1.6) +
  facet_wrap(~ alab, nrow = 1) +
  scale_colour_manual(values = n_cols) +
  labs(x = "A/σ", y = "Type I error", colour = NULL, title = "C") +
  theme_bayrc() + letter + theme(strip.text.x = element_text(size = 9))

## the amplitude curves run close together at the larger sample sizes, so they
## are offset along the nominal level to keep them apart
nudge_d <- position_dodge(width = 0.014)
fdr_panel_data <- droplevels(pooled[pooled$A != 0.5, ])
pD <- ggplot(fdr_panel_data, aes(alpha, fdr, colour = snr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = bayrc_ink) +
  geom_errorbar(aes(ymin = fdr - 1.96 * fdr_se, ymax = fdr + 1.96 * fdr_se), width = 0.006,
                position = nudge_d) +
  geom_line(position = nudge_d) + geom_point(size = 1.6, position = nudge_d) +
  facet_wrap(~ nlab, nrow = 1) +
  # every alpha is plotted; labelling every other one keeps the facets legible
  scale_x_continuous(breaks = alphas[c(TRUE, FALSE)]) +
  scale_colour_manual(values = snr_cols) +
  labs(x = "Nominal FDR", y = "Realised FDR", colour = NULL, title = "D") +
  theme_bayrc() + letter +
  # the outer tick labels of neighbouring facets would otherwise touch
  theme(panel.spacing.x = grid::unit(5, "mm"))

pE <- ggplot(fdr_panel_data, aes(alpha, bh_fdr, colour = snr)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = bayrc_ink) +
  geom_errorbar(aes(ymin = bh_fdr - 1.96 * bh_fdr_se, ymax = bh_fdr + 1.96 * bh_fdr_se),
                width = 0.006, position = nudge_d) +
  geom_line(position = nudge_d) + geom_point(size = 1.6, position = nudge_d) +
  facet_wrap(~ nlab, nrow = 1) +
  scale_x_continuous(breaks = alphas[c(TRUE, FALSE)]) +
  scale_colour_manual(values = snr_cols) +
  labs(x = "Nominal FDR", y = "Realised FDR", colour = NULL, title = "E") +
  theme_bayrc() + letter +
  theme(panel.spacing.x = grid::unit(5, "mm"))

## the three legends are different widths, so patchwork stacks the panels and
## lines their plotting areas up down the page
stacked <- patchwork::wrap_plots(pA, pB, pC, pD, pE, ncol = 1)

out <- file.path(BAYRC_FIGURE_DIR, "Figure_S_operating.pdf")
if (Sys.info()[["sysname"]] == "Darwin") {
  quartz(type = "pdf", file = out, width = 9.5, height = 15, family = bayrc_family)
} else {
  cairo_pdf(out, width = 9.5, height = 15, family = bayrc_family)
}
print(stacked)
invisible(dev.off())
cat("wrote", out, "with replicate counts", paste(sort(unique(pooled$replicates)), collapse = ", "), "per setting\n")
print(format(auc_pooled, digits = 3), row.names = FALSE)
print(format(pooled[, c("n", "A", "alpha", "fdr", "power", "type1")], digits = 3), row.names = FALSE)
