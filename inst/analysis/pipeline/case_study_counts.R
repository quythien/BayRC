## Recompute every gene-level count the Applications section reports, for the
## three case studies, at one choice of bfdr_alpha and phase window.
##
## Usage: Rscript case_study_counts.R [bfdr_alpha] [shift]

suppressPackageStartupMessages(library(BayRC))

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

args  <- commandArgs(trailingOnly = TRUE)
alpha <- if (length(args) >= 1) as.numeric(args[1]) else 0.25
shift <- if (length(args) >= 2) as.numeric(args[2]) else 2

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))

kegg354 <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))
kegg.env <- new.env()
load(file.path(BAYRC_PATHWAY_DIR, "kegg.pathway.list_hsa.RData"), envir = kegg.env)
kegg229 <- kegg.env$kegg.pathway.list_hsa

circ_diff <- function(a, b) { d <- abs(a - b); ifelse(d > 12, 24 - d, d) }

glr <- function(pA, pB) {
  n_union <- sum(pA + pB - pA * pB)
  gain <- sum((1 - pA) * pB) / n_union
  loss <- sum(pA * (1 - pB)) / n_union
  c(gain = gain, loss = loss, cons = sum(pA * pB) / n_union, glr = gain / loss)
}

report <- function(label, A, B, phiA, phiB) {
  cat("\n\n########", label, " alpha =", alpha, " shift =", shift, "h ########\n")
  pA <- rowMeans(A); pB <- rowMeans(B)
  cat("rhythmic in A (tau_A):", sum(pA >= bfdr_from_posterior(pA, alpha)$threshold),
      sprintf("(%.1f%%)", 100 * mean(pA >= bfdr_from_posterior(pA, alpha)$threshold)), "\n")
  cat("rhythmic in B (tau_B):", sum(pB >= bfdr_from_posterior(pB, alpha)$threshold),
      sprintf("(%.1f%%)", 100 * mean(pB >= bfdr_from_posterior(pB, alpha)$threshold)), "\n")
  tr <- transition_classify(pA, pB, bfdr_alpha = alpha)
  st <- tr$gain_loss_status
  g <- glr(pA, pB)
  cat(sprintf("expected proportions  gain %.3f  loss %.3f  cons %.3f  GLR %.3f\n",
              g["gain"], g["loss"], g["cons"], g["glr"]))
  ph <- phase_infer(phi_matrix1 = phiA, phi_matrix2 = phiB,
                    gain_loss_status = st, bfdr_alpha = alpha,
                    shift = shift, P = 24, compute_hdi = FALSE)
  maint <- names(st)[st == "Maintained"]
  d <- circ_diff(ph$peak1[maint], ph$peak2[maint])
  cat(sprintf("maintained %d | within +/-2h %d (%.1f%%) | within +/-1h %d (%.1f%%)\n",
              length(maint), sum(d <= 2), 100 * mean(d <= 2),
              sum(d <= 1), 100 * mean(d <= 1)))
  sh <- names(ph$flag_shift)[ph$flag_shift]
  co <- names(ph$flag_cons)[ph$flag_cons]
  cat(sprintf("phase-shifted %d (%.1f%% of maintained) | phase-conserved %d (%.1f%%) | undetermined %d\n",
              length(sh), 100 * length(sh) / length(maint),
              length(co), 100 * length(co) / length(maint),
              length(maint) - length(sh) - length(co)))
  if (length(co)) cat(sprintf("mean |delta| among phase-conserved %.2f h\n",
                              mean(circ_diff(ph$peak1[co], ph$peak2[co]))))
  invisible(list(status = st, phase = ph, pA = pA, pB = pB))
}

r1 <- report("SCN vs HIP (baboon)", mcmc_data_baboon$SCN, mcmc_data_baboon$HIP,
             mcmc_phi_baboon$SCN, mcmc_phi_baboon$HIP)
r2 <- report("PUT vs SUN (baboon)", mcmc_data_baboon$PUT, mcmc_data_baboon$SUN,
             mcmc_phi_baboon$PUT, mcmc_phi_baboon$SUN)
r3 <- report("baboon LUN vs human LUN", mcmc_data_baboon$LUN, mcmc_data_human$LUN,
             mcmc_phi_baboon$LUN, mcmc_phi_human$LUN)

## Parkinson pathway breakdown for SUN-PUT, under both KEGG releases
for (tag in c("kegg229", "kegg354")) {
  L <- get(tag)
  nm <- grep("Parkinson", names(L), value = TRUE)
  if (!length(nm)) next
  gset <- intersect(L[[nm]], names(r2$status))
  st <- r2$status[gset]
  sh <- names(r2$phase$flag_shift)[r2$phase$flag_shift]
  co <- names(r2$phase$flag_cons)[r2$phase$flag_cons]
  dphi <- ((r2$phase$peak2[gset] - r2$phase$peak1[gset] + 12) %% 24) - 12
  msh <- intersect(gset, sh); mco <- intersect(gset, co)
  cat("\n--", tag, "|", nm, "--\n")
  cat(sprintf("expressed %d | maintained %d | loss %d | gain %d\n",
              length(gset), sum(st == "Maintained"), sum(st == "Loss"),
              sum(st == "Gain")))
  cat(sprintf("phase-shifted %d | phase-conserved %d | undetermined %d\n",
              length(msh), length(mco),
              sum(st == "Maintained") - length(msh) - length(mco)))
  if (length(msh)) {
    cat(sprintf("shifted mean delta (SUN - PUT) %+.2f h, SD %.2f, range %.1f to %.1f, same sign %.0f%%\n",
                mean(dphi[msh]), sd(dphi[msh]), min(dphi[msh]), max(dphi[msh]),
                100 * max(mean(dphi[msh] > 0), mean(dphi[msh] < 0))))
    print(round(sort(dphi[msh]), 2))
  }
  if (length(mco)) { cat("phase-conserved:\n"); print(round(dphi[mco], 2)) }
}

## core clock genes in cross-species lung
clock <- c("CRY1", "CRY2", "PER1", "PER2", "PER3", "DBP", "NR1D1", "NR1D2",
           "ARNTL", "BMAL1", "CLOCK", "NPAS2", "RORA", "BHLHE40", "BHLHE41",
           "NFIL3", "TEF", "HLF", "CIART")
maint3 <- names(r3$status)[r3$status == "Maintained"]
cc <- intersect(clock, maint3)
cat("\n-- cross-species lung, conserved clock genes --\n")
if (length(cc)) {
  d <- circ_diff(r3$phase$peak1[cc], r3$phase$peak2[cc])
  print(data.frame(gene = cc, abs_delta = round(d, 2),
                   within2h = d <= 2), row.names = FALSE)
}
