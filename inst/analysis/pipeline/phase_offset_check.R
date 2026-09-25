################################################################################
# Signed phase offset for the Figure 2 panels, for the previous run and the
# current run. The 3B caption reports an offset of about 3 h above the diagonal.
#
# Sign convention, as in the case-study scripts:
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
# deltaPhi.Est is filled in only when compute_hdi = TRUE, so it is requested here.
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

runs <- c(
  previous = Sys.getenv("BAYRC_BASELINE_SUMMARY",
    unset = file.path(BAYRC_GTEX_DIR, "result", "summary", "hb")),
  current  = Sys.getenv("BAYRC_FIXED_SUMMARY",
    unset = file.path(BAYRC_GTEX_DIR, "result_fixed", "summary", "hb")))

# panel, y-axis tissue (matrix1), x-axis tissue (matrix2)
panels <- list(
  list(panel = "3A", y = "SCN_HIP_y_HIP", t1 = "HIP", t2 = "SCN"),
  list(panel = "3B", y = "SUN_PUT_y_PUT", t1 = "PUT", t2 = "SUN"))

one_run <- function(dir) {
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
    # the caption's "within +/-2 h" count uses the peak difference, not deltaPhi.Est
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

res <- do.call(rbind, lapply(names(runs), function(a) {
  cat("--", a, "\n")
  r <- one_run(runs[[a]])
  r$run <- a
  r
}))
res <- res[, c("run", "panel", "pair", "Rc", "within_2h", "pct_within",
               "median_delta", "q25", "q75", "IQR", "median_peakdiff",
               "frac_above", "n_na")]

options(width = 200)
cat("\nsigned phase offset over the maintained set\n")
cat("delta = phi(y axis) - phi(x axis), wrapped to [-12, 12); positive is above the diagonal\n")
cat("3A y = HIP, x = SCN     3B y = PUT, x = SUN\n\n")
print(res[order(res$panel, res$run), ], row.names = FALSE)

out.file <- file.path(BAYRC_OUTPUT_DIR, "phase_offset_check.csv")
dir.create(dirname(out.file), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, out.file, row.names = FALSE)
cat("\nwrote", out.file, "\n")
