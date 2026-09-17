################################################################################
# Signed phase offset for the Figure 3 panels, both samplers.
#
# The published 3B caption reads a "systematic +3 h offset visible as
# displacement above the diagonal". This recomputes that offset so it can be
# checked against the fixed sampler.
#
# Sign convention, taken from the case-study scripts rather than assumed:
#
#   3A  Baboon_SCN_HIP.R  phi_matrix1 = HIP, phi_matrix2 = SCN
#                         plot x = SCN, y = Hippocampus
#   3B  Baboon_SUN_PUT.R  phi_matrix1 = PUT, phi_matrix2 = SUN
#                         plot x = Substantia nigra, y = Putamen
#
# phase_infer() builds its difference as phi_matrix1 - phi_matrix2 wrapped to
# [-12, 12), so in both panels matrix1 is the y axis and
#
#   deltaPhi > 0  ==  the gene sits ABOVE the diagonal.
#
# deltaPhi.Est is only filled in when compute_hdi = TRUE, which is why it
# reads NA at the default; it is requested explicitly here.
#
# Usage:
#   Rscript phase_offset_check.R
################################################################################

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(if (is.na(this.file)) dirname(getwd()) else
                   dirname(dirname(normalizePath(this.file))), "config.R"))

suppressPackageStartupMessages(library(BayRC))

BFDR_ALPHA <- 0.25
SHIFT      <- 2
P          <- 24

# the scripts' own helper, so the peak columns are on the same scale as the plot
to_zt <- function(t_cos) ifelse(t_cos >= 18, t_cos - 24, t_cos)

arms <- c(
  baseline = Sys.getenv("BAYRC_BASELINE_SUMMARY",
    unset = file.path(BAYRC_GTEX_DIR, "result", "summary", "hb")),
  fixed    = Sys.getenv("BAYRC_FIXED_SUMMARY",
    unset = file.path(BAYRC_GTEX_DIR, "result_fixed", "summary", "hb")))

# panel, y-axis tissue (matrix1), x-axis tissue (matrix2)
panels <- list(
  list(panel = "3A", y = "SCN_HIP_y_HIP", t1 = "HIP", t2 = "SCN"),
  list(panel = "3B", y = "SUN_PUT_y_PUT", t1 = "PUT", t2 = "SUN"))

one_arm <- function(dir) {
  e <- new.env()
  load(file.path(dir, "mcmc_rho_BF3.RData"), e)
  load(file.path(dir, "phi", "mcmc_phi_BF3.RData"), e)
  out <- lapply(panels, function(p) {
    rho1 <- e$mcmc_data_baboon[[p$t1]]; phi1 <- e$mcmc_phi_baboon[[p$t1]]
    rho2 <- e$mcmc_data_baboon[[p$t2]]; phi2 <- e$mcmc_phi_baboon[[p$t2]]

    trans <- transition_classify(rowMeans(rho1), rowMeans(rho2),
                                 bfdr_alpha = BFDR_ALPHA)
    inner <- phase_infer(phi_matrix1 = phi1, phi_matrix2 = phi2,
                         gain_loss_status = trans$gain_loss_status,
                         bfdr_alpha = BFDR_ALPHA, shift = SHIFT, P = P,
                         compute_hdi = TRUE)

    keep <- names(trans$gain_loss_status)[trans$gain_loss_status == "Maintained"]
    # The published "within +/-2 h" count comes from the difference of the two
    # peak estimates, not from deltaPhi.Est: for 3B the peak difference gives
    # 116 of 553 (21.0%), which is the number in the caption, while
    # deltaPhi.Est gives 118 (21.3%). They agree exactly for 3A.
    dp <- ((inner$peak1[keep] - inner$peak2[keep] + P/2) %% P) - P/2
    d  <- inner$deltaPhi.Est[keep]
    within <- abs(dp) < SHIFT

    data.frame(
      panel = p$panel, pair = paste0(p$t2, " vs ", p$t1),
      Rc = length(keep),
      within_2h = sum(within, na.rm = TRUE),
      pct_within = round(100 * mean(within, na.rm = TRUE), 1),
      median_delta = round(stats::median(d, na.rm = TRUE), 2),
      q25 = round(stats::quantile(d, 0.25, na.rm = TRUE), 2),
      q75 = round(stats::quantile(d, 0.75, na.rm = TRUE), 2),
      IQR = round(stats::IQR(d, na.rm = TRUE), 2),
      median_peakdiff = round(stats::median(dp, na.rm = TRUE), 2),
      # share of maintained genes sitting above the plotted diagonal
      frac_above = round(mean(dp > 0, na.rm = TRUE), 3),
      n_na = sum(is.na(dp)),
      stringsAsFactors = FALSE)
  })
  rm(e); invisible(gc())
  do.call(rbind, out)
}

res <- do.call(rbind, lapply(names(arms), function(a) {
  cat("--", a, "\n")
  r <- one_arm(arms[[a]])
  r$arm <- a
  r
}))
res <- res[, c("arm", "panel", "pair", "Rc", "within_2h", "pct_within",
               "median_delta", "q25", "q75", "IQR", "median_peakdiff",
               "frac_above", "n_na")]

options(width = 200)
cat("\nsigned phase offset over the maintained set\n")
cat("delta = phi(y axis) - phi(x axis), wrapped to [-12, 12); positive is above the diagonal\n")
cat("3A y = HIP, x = SCN     3B y = PUT, x = SUN\n\n")
print(res[order(res$panel, res$arm), ], row.names = FALSE)

out.file <- file.path(BAYRC_OUTPUT_DIR, "phase_offset_check.csv")
dir.create(dirname(out.file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, out.file, row.names = FALSE)
cat("\nwrote", out.file, "\n")
