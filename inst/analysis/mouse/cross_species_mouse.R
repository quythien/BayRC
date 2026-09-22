## Human against mouse across the six tissues the GSE54651 atlas shares with
## GTEx, on the frozen parameters the rest of the paper uses.
##
## Baboon against human is computed on the same rows for reference, so the two
## cross-species comparisons can be read side by side and the section can be
## put on the tissue that carries it best rather than on lung by default.

args   <- commandArgs(trailingOnly = TRUE)
hmdir  <- args[1]
hbdir  <- args[2]
outdir <- args[3]

suppressPackageStartupMessages(library(BayRC))

bfdr_alpha <- 0.25
shift      <- 2

## circular difference in hours, taking the shorter way round, as in
## plots/peak_concordance.R
circ_diff <- function(a, b, P = 24) {
  d <- abs(a - b)
  ifelse(d > P / 2, P - d, d)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

load(file.path(hmdir, "mcmc_rho_BF3.RData"))
load(file.path(hmdir, "phi", "mcmc_phi_BF3.RData"))

tissues <- names(mcmc_data_mouse)
genes   <- rownames(mcmc_data_mouse[[1]])

clock_genes <- c("BMAL1", "CLOCK", "PER1", "PER2", "CRY1", "CRY2",
                 "NR1D1", "NR1D2", "RORA", "DBP", "NFIL3")

## condition A is human throughout, so both comparisons are human-referenced
one_pair <- function(tis, comparison, A, B, nameA, nameB) {
  pA <- rowMeans(A$rho)
  pB <- rowMeans(B$rho)
  trans <- transition_classify(pA, pB, bfdr_alpha = bfdr_alpha)
  status <- trans$gain_loss_status

  phase <- phase_infer(phi_matrix1 = A$phi, phi_matrix2 = B$phi,
                       gain_loss_status = status, bfdr_alpha = bfdr_alpha,
                       shift = shift, P = 24, compute_hdi = TRUE)

  maintained <- names(status)[status == "Maintained"]
  within <- if (length(maintained))
    circ_diff(phase$peak1[maintained], phase$peak2[maintained]) <= shift else
    logical(0)

  conc <- multi_conservation(
    mcmc.merge.list     = list(A, B),
    dataset.names       = c(nameA, nameB),
    select.pathway.list = "global",
    output.dir = tempdir(), save_output = FALSE, use_cpp = TRUE,
    compute_pvalue = FALSE, compute_ci = FALSE)
  r <- if (is.data.frame(conc)) conc else conc[[1]]
  grab <- function(suffix) {
    j <- grep(paste0("_", suffix, "$"), names(r))
    if (!length(j)) NA_real_ else as.numeric(r[1, j[1]])
  }

  data.frame(
    comparison        = comparison,
    tissue            = tis,
    genes             = length(pA),
    mean_post_A       = round(mean(pA), 4),
    mean_post_B       = round(mean(pB), 4),
    rhythmic_A        = sum(pA >= trans$tau_rhythmic_A),
    rhythmic_B        = sum(pB >= trans$tau_rhythmic_B),
    gain              = trans$n_gain,
    loss              = trans$n_loss,
    maintained        = length(maintained),
    within_window     = sum(within),
    within_pct        = if (length(maintained))
                          round(100 * mean(within), 1) else NA_real_,
    phase_conserved   = sum(phase$flag_cons),
    phase_shifted     = sum(phase$flag_shift),
    concordance_raw   = round(grab("Concordance"), 4),
    concordance_adj   = round(grab("AdjustedConcordance"), 4),
    stringsAsFactors  = FALSE)
}

cat("condition A is human, bfdr_alpha", bfdr_alpha, "shift", shift, "h\n\n")

hm <- do.call(rbind, lapply(tissues, function(tis)
  one_pair(tis, "human_vs_mouse",
           list(rho = mcmc_data_human[[tis]], phi = mcmc_phi_human[[tis]]),
           list(rho = mcmc_data_mouse[[tis]], phi = mcmc_phi_mouse[[tis]]),
           paste0("hum_", tis), paste0("mou_", tis))))

## Baboon on the same rows, read from the human/baboon summary and cut to the
## three-way gene set.
e <- new.env()
load(file.path(hbdir, "mcmc_rho_BF3.RData"), envir = e)
load(file.path(hbdir, "phi", "mcmc_phi_BF3.RData"), envir = e)

cut_to <- function(m) {
  idx <- match(genes, rownames(m))
  out <- m[idx, , drop = FALSE]
  attr(out, "symbols") <- genes
  out
}

hb <- do.call(rbind, lapply(tissues, function(tis)
  one_pair(tis, "human_vs_baboon",
           list(rho = mcmc_data_human[[tis]], phi = mcmc_phi_human[[tis]]),
           list(rho = cut_to(e$mcmc_data_baboon[[tis]]),
                phi = cut_to(e$mcmc_phi_baboon[[tis]])),
           paste0("hum_", tis), paste0("bab_", tis))))
rm(e); invisible(gc())

tab <- rbind(hm, hb)
utils::write.csv(tab, file.path(outdir, "cross_species_six_tissues.csv"),
                 row.names = FALSE)

cat("\n=== human against mouse, ranked by adjusted concordance ===\n")
print(hm[order(-hm$concordance_adj), ], row.names = FALSE)
cat("\n=== human against baboon, same rows, for reference ===\n")
print(hb[order(-hb$concordance_adj), ], row.names = FALSE)

## Posterior rhythmicity of the core clock. The mouse matrix is keyed on human
## identifiers, so the rows carry human symbols.
clock <- do.call(rbind, lapply(tissues, function(tis) {
  g <- intersect(clock_genes, rownames(mcmc_data_mouse[[tis]]))
  data.frame(tissue = tis, gene = g,
             human = round(rowMeans(mcmc_data_human[[tis]])[g], 3),
             mouse = round(rowMeans(mcmc_data_mouse[[tis]])[g], 3),
             stringsAsFactors = FALSE)
}))
utils::write.csv(clock, file.path(outdir, "cross_species_clock.csv"),
                 row.names = FALSE)

cat("\n=== core clock posterior rhythmicity, mouse ===\n")
wide <- reshape(clock[, c("tissue", "gene", "mouse")], idvar = "gene",
                timevar = "tissue", direction = "wide")
names(wide) <- sub("^mouse\\.", "", names(wide))
print(wide, row.names = FALSE)

cat("\n=== core clock posterior rhythmicity, human ===\n")
wideh <- reshape(clock[, c("tissue", "gene", "human")], idvar = "gene",
                 timevar = "tissue", direction = "wide")
names(wideh) <- sub("^human\\.", "", names(wideh))
print(wideh, row.names = FALSE)

record <- c(
  paste("script       ", "src/cross_species_mouse.R"),
  paste("mode         ", "full run"),
  paste("BayRC version", as.character(utils::packageVersion("BayRC"))),
  paste("summaries    ", hmdir),
  paste("baboon source", hbdir),
  paste("rho written  ", format(file.info(file.path(hmdir,
        "mcmc_rho_BF3.RData"))$mtime, "%Y-%m-%d %H:%M")),
  paste("run at       ", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "parameters",
  sprintf("  %-14s %s", "bfdr_alpha", bfdr_alpha),
  sprintf("  %-14s %s", "shift", shift),
  sprintf("  %-14s %s", "P", 24),
  sprintf("  %-14s %s", "genes", length(genes)),
  sprintf("  %-14s %s", "direction", "condition A human in both comparisons"),
  sprintf("  %-14s %s", "tissues", paste(tissues, collapse = ", ")))
writeLines(record, file.path(outdir, "run_record.txt"))

cat("\nwrote", file.path(outdir, "cross_species_six_tissues.csv"), "\n")
