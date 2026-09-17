################################################################################
# All pairwise genome-wide concordance, for one set of tissues.
#
# Replaces the serial double loop in pairwise_Con.R and the one in
# plots/heatmap_baboon.R. Each pair is independent, so the OUTER pair loop is
# what gets parallelised; the permutation and bootstrap loops inside
# multi_conservation() stay serial, because nesting parallelism there would
# oversubscribe the node and does not help -- there are far more pairs than
# cores.
#
# Usage:
#   Rscript pairwise_concordance_all.R <mode> [cores]
#
#   within_baboon   choose(26, 2) = 325 pairs
#   within_human    choose(26, 2) = 325 pairs
#   cross_matched   26 pairs, baboon T against human T -- this is what the
#                   paper reports; its stored table covers 23 of the 26,
#                   leaving out MUA, SCN and VIC
#   cross_species   26 x 26 = 676 pairs, every baboon-human combination. Not
#                   what the paper reports; useful as a null for how special
#                   the matched diagonal is.
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
# "infer" adds the permutation p-value and bootstrap CI. Figure 2 does not use
# them -- plots/heatmap_baboon.R calls multi_conservation() with
# compute_pvalue = FALSE and compute_ci = FALSE and plots the adjusted
# concordance alone -- and they cost about three orders of magnitude more time
# per pair, so they are off unless asked for.
infer <- length(args) >= 3 && args[3] == "infer"

N_PERM <- 1000
N_BOOT <- 1000

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))

# multi_conservation() reads the gene names off attr(rho, "symbols"); the
# summary builder sets it, but a chain read straight from an RDS may not have.
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

# Each pair already holds references to the two chains it needs, so the
# species-level lists can go. They are ~4 GB each, and every forked worker
# that touches them makes the kernel copy those pages -- with 18 workers that
# was enough to fill the node.
if (mode == "within_baboon") rm(human)
if (mode == "within_human")  rm(baboon)
invisible(gc())

cat(mode, ":", length(pairs), "pairs on", cores, "cores",
    if (infer) "with permutation and bootstrap inference" else
               "adjusted concordance only", "\n")

one_pair <- function(p) {
  res <- try(multi_conservation(
    mcmc.merge.list     = list(p$da, p$db),
    dataset.names       = c(p$la, p$lb),
    select.pathway.list = "global",
    n_perm = N_PERM, n_boot = N_BOOT,
    output.dir = tempdir(), save_output = FALSE,
    use_cpp = TRUE,
    compute_pvalue = infer, compute_ci = infer), silent = TRUE)
  if (inherits(res, "try-error"))
    return(data.frame(tissue1 = p$a, tissue2 = p$b, raw = NA_real_,
                      adjusted = NA_real_, ci_lower = NA_real_,
                      ci_upper = NA_real_, p_value = NA_real_,
                      z_score = NA_real_, null_mean = NA_real_,
                      null_sd = NA_real_, error = as.character(res),
                      stringsAsFactors = FALSE))
  r <- res[[1]]
  data.frame(tissue1 = p$a, tissue2 = p$b,
             raw = r[1, "Conservation"], adjusted = r[1, "Adjusted"],
             ci_lower = r[1, "CI_lower_adj"], ci_upper = r[1, "CI_upper_adj"],
             p_value = r[1, "PValue"], z_score = r[1, "ZScore"],
             null_mean = r[1, "Null_mean"], null_sd = r[1, "Null_sd"],
             error = NA_character_, stringsAsFactors = FALSE)
}

t0 <- Sys.time()
out <- mclapply(pairs, one_pair, mc.cores = cores, mc.preschedule = FALSE)
cat("elapsed", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
    "min\n")

# mclapply returns the condition object rather than throwing when a child dies,
# so anything that is not a one-row data.frame is a failed pair.
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
