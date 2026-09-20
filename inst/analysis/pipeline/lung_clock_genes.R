## Per-gene detail for the core clock genes in the baboon-human lung pair:
## where each one sits relative to the maintained threshold, and how the phase
## posterior classifies it at two BFDR levels.
##
## Sign convention: delta = human peak - baboon peak, on the circle, positive
## when human peaks later. Reported as the posterior median.
##
## Usage: Rscript lung_clock_genes.R [alpha1] [alpha2] [shift]

suppressPackageStartupMessages(library(BayRC))

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

args <- commandArgs(trailingOnly = TRUE)
a1    <- if (length(args) >= 1) as.numeric(args[1]) else 0.25
a2    <- if (length(args) >= 2) as.numeric(args[2]) else 0.30
shift <- if (length(args) >= 3) as.numeric(args[3]) else 2

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))

rhoA <- mcmc_data_baboon$LUN; rhoB <- mcmc_data_human$LUN
phiA <- mcmc_phi_baboon$LUN;  phiB <- mcmc_phi_human$LUN
pA <- rowMeans(rhoA); pB <- rowMeans(rhoB)

## the panel the paper's clock claim is drawn from
clock <- c("BMAL1", "ARNTL", "CLOCK", "NPAS2", "PER1", "PER2", "PER3",
           "CRY1", "CRY2", "NR1D1", "NR1D2", "RORA", "RORB", "RORC",
           "DBP", "TEF", "HLF", "NFIL3", "BHLHE40", "BHLHE41", "CIART",
           "CSNK1D", "CSNK1E", "FBXL3")
clock <- intersect(clock, names(pA))

## signed phase difference, human minus baboon
dphi <- ((phiB - phiA + 12) %% 24) - 12
med  <- apply(dphi[clock, , drop = FALSE], 1, median)
absd <- abs(med)

## lowest bfdr_alpha at which each gene is called maintained
alphas <- seq(0.01, 0.95, by = 0.01)
first <- setNames(rep(NA_real_, length(clock)), clock)
for (a in alphas) {
  st <- capture.output(tr <- transition_classify(pA, pB, bfdr_alpha = a))
  m <- names(tr$gain_loss_status)[tr$gain_loss_status == "Maintained"]
  hit <- clock[clock %in% m & is.na(first[clock])]
  if (length(hit)) first[hit] <- a
  if (!anyNA(first)) break
}

at_alpha <- function(a) {
  st <- capture.output(tr <- transition_classify(pA, pB, bfdr_alpha = a))
  ph <- phase_infer(phi_matrix1 = phiA, phi_matrix2 = phiB,
                    gain_loss_status = tr$gain_loss_status,
                    bfdr_alpha = a, shift = shift, P = 24, compute_hdi = FALSE)
  cls <- rep("not maintained", length(clock)); names(cls) <- clock
  maint <- clock[tr$gain_loss_status[clock] == "Maintained"]
  cls[maint] <- "undetermined"
  cls[clock[ph$flag_cons[clock]]]  <- "phase-conserved"
  cls[clock[ph$flag_shift[clock]]] <- "phase-shifted"
  data.frame(maintained = clock %in% maint,
             p_cons = round(ph$prob_conserved[clock], 3),
             p_shift = round(ph$prob_shift[clock], 3),
             class = cls, row.names = clock)
}

r1 <- at_alpha(a1); r2 <- at_alpha(a2)

out <- data.frame(
  gene = clock,
  pA_baboon = round(pA[clock], 3),
  pB_human  = round(pB[clock], 3),
  delta_h   = round(med, 2),
  abs_delta = round(absd, 2),
  within_2h = absd <= shift,
  first_alpha_maintained = first[clock],
  maint_a1 = r1$maintained, p_cons_a1 = r1$p_cons,
  p_shift_a1 = r1$p_shift, class_a1 = r1$class,
  maint_a2 = r2$maintained, p_cons_a2 = r2$p_cons,
  p_shift_a2 = r2$p_shift, class_a2 = r2$class,
  row.names = NULL)
out <- out[order(out$first_alpha_maintained, out$abs_delta), ]

out.dir <- file.path(BAYRC_OUTPUT_DIR, "param_search")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(out.dir, sprintf("lung_clock_genes_a%g_a%g.csv", a1, a2)),
          row.names = FALSE)

cat("\nsign convention: delta = human peak - baboon peak (positive = human later)\n")
cat("alpha columns: a1 =", a1, " a2 =", a2, " | phase window =", shift, "h\n\n")
print(out, row.names = FALSE)

seven <- c("CRY1", "CRY2", "PER2", "DBP", "NR1D1", "NR1D2", "BMAL1")
s <- out[out$gene %in% seven, ]
cat("\nthe seven genes the paper names:\n")
cat("  within", shift, "h window:", sum(s$within_2h), "of", nrow(s), "\n")
cat("  maintained at alpha", a1, ":", sum(s$maint_a1), "\n")
cat("  maintained at alpha", a2, ":", sum(s$maint_a2), "\n")
