## Signed circular phase difference SUN - PUT over all maintained genes, and
## over the shifted subset, to separate the offset from the shift threshold.
##
## Usage: Rscript sunput_global_offset.R [shift]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
suppressPackageStartupMessages(library(BayRC))

args  <- commandArgs(trailingOnly = TRUE)
shift <- if (length(args) >= 1) as.numeric(args[1]) else 2

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))
rho <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
phi <- get(grep("^mcmc_phi",  ls(), value = TRUE)[1])

tr <- transition_classify(rowMeans(rho[["PUT"]]), rowMeans(rho[["SUN"]]),
                          bfdr_alpha = 0.25)
ph <- phase_infer(phi_matrix1 = phi[["PUT"]], phi_matrix2 = phi[["SUN"]],
                  gain_loss_status = tr$gain_loss_status,
                  bfdr_alpha = 0.25, shift = shift, P = 24, compute_hdi = FALSE)

d <- (ph$peak2 - ph$peak1) %% 24
d <- ifelse(d > 12, d - 24, d)
maint <- tr$gain_loss_status == "Maintained"

report <- function(x, label) {
  x <- x[is.finite(x)]
  cat(sprintf("%-22s n = %3d  mean %+5.2f  median %+5.2f  sd %4.2f  %%pos %5.1f\n",
              label, length(x), mean(x), median(x), sd(x), 100 * mean(x > 0)))
}
report(d[maint], "all maintained")
report(d[maint & ph$flag_shift %in% TRUE], "shifted only")
report(d[maint & ph$flag_cons  %in% TRUE], "phase-conserved only")
report(abs(d[maint]), "all maintained |diff|")

cat("\nfraction of maintained genes with |diff| above the threshold:",
    sprintf("%.3f\n", mean(abs(d[maint]) > shift, na.rm = TRUE)))
cat("one-sample test that the maintained offset is zero:\n")
print(t.test(d[maint]))
