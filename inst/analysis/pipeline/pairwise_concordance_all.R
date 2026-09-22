################################################################################
# All pairwise genome-wide concordance, for one set of tissues.
#
# The pair loop runs in parallel; the permutation and bootstrap loops inside
# multi_conservation() stay serial.
#
# Usage:
#   Rscript pairwise_concordance_all.R <mode> [cores]
#
#   within_baboon   choose(26, 2) = 325 pairs
#   within_human    choose(26, 2) = 325 pairs
#   cross_matched   26 pairs, baboon T against human T, as the paper reports;
#                   its stored table covers 23 of the 26, leaving out MUA, SCN
#                   and VIC
#   cross_species   26 x 26 = 676 pairs, every baboon-human combination, as a
#                   reference for the matched diagonal
#
# Writes <mode>_pairwise_concordance.csv and .rds to BAYRC_OUTPUT_DIR.
################################################################################

# Paths come from inst/analysis/config.R; override with the matching env var.
this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(if (is.na(this.file)) dirname(getwd()) else
                   dirname(dirname(normalizePath(this.file))), "config.R"))

suppressPackageStartupMessages({
  library(BayRC)
  library(parallel)
})

args  <- commandArgs(trailingOnly = TRUE)
mode  <- if (length(args) >= 1) args[1] else "within_baboon"
cores <- if (length(args) >= 2) as.integer(args[2]) else
  max(1L, parallel::detectCores() - 4L)
# "infer" adds the permutation p-value and bootstrap CI; Figure 2 uses neither.
infer <- length(args) >= 3 && args[3] == "infer"

N_PERM <- 1000
N_BOOT <- 1000

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))

# multi_conservation() reads the gene names off attr(rho, "symbols")
as_dataset <- function(rho, phi) {
  if (is.null(attr(rho, "symbols"))) attr(rho, "symbols") <- rownames(rho)
  list(rho = rho, phi = phi)
}

baboon <- setNames(lapply(names(mcmc_data_baboon), function(t)
  as_dataset(mcmc_data_baboon[[t]], mcmc_phi_baboon[[t]])), names(mcmc_data_baboon))
human <- setNames(lapply(names(mcmc_data_human), function(t)
  as_dataset(mcmc_data_human[[t]], mcmc_phi_human[[t]])), names(mcmc_data_human))
rm(mcmc_data_baboon, mcmc_data_human, mcmc_phi_baboon, mcmc_phi_human)
invisible(gc())

pairs <- switch(mode,
  within_baboon = {
    cb <- combn(names(baboon), 2, simplify = FALSE)
    lapply(cb, function(p) list(a = p[1], b = p[2],
                                da = baboon[[p[1]]], db = baboon[[p[2]]],
                                la = paste0("bab_", p[1]), lb = paste0("bab_", p[2])))
  },
  within_human = {
    cb <- combn(names(human), 2, simplify = FALSE)
    lapply(cb, function(p) list(a = p[1], b = p[2],
                                da = human[[p[1]]], db = human[[p[2]]],
                                la = paste0("hum_", p[1]), lb = paste0("hum_", p[2])))
  },
  cross_matched = {
    tis <- intersect(names(baboon), names(human))
    lapply(tis, function(t)
      list(a = t, b = t, da = baboon[[t]], db = human[[t]],
           la = paste0("bab_", t), lb = paste0("hum_", t)))
  },
  cross_species = {
    grid <- expand.grid(b = names(baboon), h = names(human),
                        stringsAsFactors = FALSE)
    lapply(seq_len(nrow(grid)), function(i) {
      b <- grid$b[i]; h <- grid$h[i]
      list(a = b, b = h, da = baboon[[b]], db = human[[h]],
           la = paste0("bab_", b), lb = paste0("hum_", h))
    })
  },
  stop("mode must be within_baboon, within_human, cross_matched or cross_species"))

# each pair holds its own two chains; drop the ~4 GB species lists before forking
if (mode == "within_baboon") rm(human)
if (mode == "within_human")  rm(baboon)
invisible(gc())

cat(mode, ":", length(pairs), "pairs on", cores, "cores",
    if (infer) "with permutation and bootstrap inference" else
               "adjusted concordance only", "\n")

blank_row <- function(p, err = NA_character_)
  data.frame(tissue1 = p$a, tissue2 = p$b, raw = NA_real_, adjusted = NA_real_,
             gain = NA_real_, loss = NA_real_, gain_loss_ratio = NA_real_,
             ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_,
             z_score = NA_real_, null_mean = NA_real_, null_sd = NA_real_,
             error = err, stringsAsFactors = FALSE)

one_pair <- function(p) {
  res <- try(multi_conservation(
    mcmc.merge.list     = list(p$da, p$db),
    dataset.names       = c(p$la, p$lb),
    select.pathway.list = "global",
    n_perm = N_PERM, n_boot = N_BOOT,
    output.dir = tempdir(), save_output = FALSE,
    use_cpp = TRUE,
    compute_pvalue = infer, compute_ci = infer), silent = TRUE)
  if (inherits(res, "try-error")) return(blank_row(p, as.character(res)))

  # a data.frame without inference, a list with it; columns are
  # "<a>_vs_<b>_<quantity>", so match on the suffix
  r <- if (is.data.frame(res)) res else res[[1]]
  if (!is.data.frame(r)) return(blank_row(p, "unexpected result shape"))
  grab <- function(suffix) {
    j <- grep(paste0("_", suffix, "$"), names(r))
    if (!length(j)) NA_real_ else as.numeric(r[1, j[1]])
  }
  out <- blank_row(p)
  out$raw             <- grab("Concordance")
  out$adjusted        <- grab("AdjustedConcordance")
  out$gain            <- grab("GainIndex")
  out$loss            <- grab("LossIndex")
  out$gain_loss_ratio <- grab("GainLossRatio")
  out$ci_lower        <- grab("CI_lower_adj")
  out$ci_upper        <- grab("CI_upper_adj")
  out$p_value         <- grab("PValue")
  out$z_score         <- grab("ZScore")
  out$null_mean       <- grab("Null_mean")
  out$null_sd         <- grab("Null_sd")
  out
}

t0 <- Sys.time()
out <- mclapply(pairs, one_pair, mc.cores = cores, mc.preschedule = FALSE)
cat("elapsed", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    "min\n")

# mclapply returns the error object when a child dies; keep one-row data.frames
bad <- !vapply(out, function(x) is.data.frame(x) && nrow(x) == 1L, logical(1))
if (any(bad)) {
  cat("WARNING:", sum(bad), "pairs failed in the worker\n")
  out <- out[!bad]
}
tab <- do.call(rbind, out)
tab <- tab[order(-tab$adjusted), ]

dir.create(BAYRC_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(BAYRC_OUTPUT_DIR, paste0(mode, "_pairwise_concordance"))
utils::write.csv(tab, paste0(stem, ".csv"), row.names = FALSE)
saveRDS(tab, paste0(stem, ".rds"))

cat("\n", mode, ": ", sum(!is.na(tab$adjusted)), " of ", length(pairs),
    " pairs succeeded\n", sep = "")
cat("adjusted concordance  median ",
    round(stats::median(tab$adjusted, na.rm = TRUE), 4),
    "  range ", paste(round(range(tab$adjusted, na.rm = TRUE), 4),
                      collapse = " to "), "\n", sep = "")
print(utils::head(tab[, c("tissue1", "tissue2", "adjusted", "p_value")], 10),
      row.names = FALSE)
cat("\nwrote", paste0(stem, ".csv"), "\n")
