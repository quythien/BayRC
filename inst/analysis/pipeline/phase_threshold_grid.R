## Sweep shift and bfdr_alpha to see which drives the undetermined class.
##
## Usage: Rscript phase_threshold_grid.R [tissueA] [tissueB] [out.csv]

suppressPackageStartupMessages(library(BayRC))

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

args <- commandArgs(trailingOnly = TRUE)
tA   <- if (length(args) >= 1) args[1] else "SUN"
tB   <- if (length(args) >= 2) args[2] else "PUT"
out  <- if (length(args) >= 3) args[3] else
  file.path(BAYRC_OUTPUT_DIR, sprintf("phase_threshold_grid_%s_%s.csv", tA, tB))

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))
rho <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
phi <- get(grep("^mcmc_phi",  ls(), value = TRUE)[1])

shifts <- c(1, 1.5, 2, 2.5, 3, 4)
alphas <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)

res <- list()
for (a in alphas) {
  tr <- transition_classify(rowMeans(rho[[tA]]), rowMeans(rho[[tB]]),
                            bfdr_alpha = a)
  n_maint <- sum(tr$gain_loss_status == "Maintained")
  for (s in shifts) {
    ph <- phase_infer(phi_matrix1 = phi[[tA]], phi_matrix2 = phi[[tB]],
                      gain_loss_status = tr$gain_loss_status,
                      bfdr_alpha = a, shift = s, P = 24, compute_hdi = FALSE)
    res[[length(res) + 1]] <- data.frame(
      alpha = a, shift = s, maintained = n_maint,
      shifted   = sum(ph$flag_shift, na.rm = TRUE),
      conserved = sum(ph$flag_cons,  na.rm = TRUE),
      undetermined = n_maint - sum(ph$flag_shift, na.rm = TRUE) -
                     sum(ph$flag_cons, na.rm = TRUE))
  }
}

d <- do.call(rbind, res)
d$pct_undet <- round(100 * d$undetermined / d$maintained, 1)
write.csv(d, out, row.names = FALSE)

cat(sprintf("\n%s vs %s\n\n", tA, tB))
cat("percent undetermined, rows = BFDR alpha, cols = shift (h)\n\n")
print(round(xtabs(pct_undet ~ alpha + shift, data = d), 1))
cat("\nmaintained genes by alpha:\n")
print(unique(d[, c("alpha", "maintained")]), row.names = FALSE)
cat("\nwrote", out, "\n")
