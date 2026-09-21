## Core clock genes in every tissue compared against the suprachiasmatic
## nucleus: transition status and phase difference relative to the master clock.

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
suppressPackageStartupMessages(library(BayRC))

bfdr_alpha <- 0.25; shift <- 2
clock <- c("CLOCK", "BMAL1", "ARNTL", "NPAS2", "PER1", "PER2", "PER3",
           "CRY1", "CRY2", "NR1D1", "NR1D2", "RORC", "DBP", "TEF", "HLF")

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))
rho <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
phi <- get(grep("^mcmc_phi",  ls(), value = TRUE)[1])

genes <- rownames(rho[["SCN"]])
clock <- clock[clock %in% genes]
cat("clock genes present:", paste(clock, collapse = " "), "\n")
cat("posterior rhythmicity in SCN:\n")
print(round(rowMeans(rho[["SCN"]])[clock], 2))

others <- setdiff(names(rho), "SCN")
res <- do.call(rbind, lapply(others, function(t) {
  tr <- transition_classify(rowMeans(rho[["SCN"]]), rowMeans(rho[[t]]),
                            bfdr_alpha = bfdr_alpha)
  ph <- phase_infer(phi_matrix1 = phi[["SCN"]], phi_matrix2 = phi[[t]],
                    gain_loss_status = tr$gain_loss_status,
                    bfdr_alpha = bfdr_alpha, shift = shift, P = 24,
                    compute_hdi = FALSE)
  d <- (ph$peak2 - ph$peak1) %% 24; d <- ifelse(d > 12, d - 24, d)
  i <- match(clock, genes)
  data.frame(tissue = t, gene = clock,
             status = tr$gain_loss_status[i],
             shifted = ph$flag_shift[i], conserved = ph$flag_cons[i],
             dphi = round(d[i], 2))
}))
write.csv(res, file.path(BAYRC_OUTPUT_DIR, "scn_clock_reference.csv"),
          row.names = FALSE)

cat("\nmaintained against SCN, out of", length(others), "tissues:\n")
m <- tapply(res$status == "Maintained", res$gene, sum)
print(sort(m, decreasing = TRUE))
cat("\nof those maintained, how many are phase-conserved within 2 h:\n")
c2 <- tapply(res$conserved %in% TRUE, res$gene, sum)
print(c2[names(sort(m, decreasing = TRUE))])
cat("\nmedian phase difference from SCN among maintained:\n")
print(round(tapply(res$dphi[res$status == "Maintained"],
                   res$gene[res$status == "Maintained"], median, na.rm = TRUE), 2))
