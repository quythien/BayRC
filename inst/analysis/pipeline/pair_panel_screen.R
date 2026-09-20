## Screen every baboon tissue pair for a pathway that would make a readable
## Figure 5 panel: enough maintained genes, split across both phase-shifted and
## phase-conserved. The panel metric is min(n_shift, n_cons) within a pathway.
##
## No enrichment is run here; this is the cheap pass over all pairs.
##
## Usage: Rscript pair_panel_screen.R [--time-one] [ncores]

suppressPackageStartupMessages({library(BayRC); library(parallel)})

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

BFDR_ALPHA <- 0.25
SHIFT <- 2
MIN_MEASURED <- 15

args <- commandArgs(trailingOnly = TRUE)
time.one <- "--time-one" %in% args
ncores <- suppressWarnings(as.integer(setdiff(args, "--time-one")[1]))
if (is.na(ncores)) ncores <- max(1, detectCores() - 4)

out.dir <- file.path(BAYRC_OUTPUT_DIR, "pair_panel_screen")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))
tissues <- names(mcmc_data_baboon)
measured <- rownames(mcmc_data_baboon[[1]])

kegg <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))
kegg <- lapply(kegg, function(g) replace(g, g == "ARNTL", "BMAL1"))
kegg <- lapply(kegg, intersect, measured)
kegg <- kegg[lengths(kegg) >= MIN_MEASURED]

pairs <- combn(tissues, 2, simplify = FALSE)
if (time.one) pairs <- pairs[1]

screen_pair <- function(pr) {
  tA <- pr[1]; tB <- pr[2]
  pA <- rowMeans(mcmc_data_baboon[[tA]]); pB <- rowMeans(mcmc_data_baboon[[tB]])
  invisible(capture.output(
    tr <- transition_classify(pA, pB, bfdr_alpha = BFDR_ALPHA)))
  st <- tr$gain_loss_status
  invisible(capture.output(
    ph <- phase_infer(phi_matrix1 = mcmc_phi_baboon[[tA]],
                      phi_matrix2 = mcmc_phi_baboon[[tB]],
                      gain_loss_status = st, bfdr_alpha = BFDR_ALPHA,
                      shift = SHIFT, P = 24, compute_hdi = FALSE)))
  shifted <- names(ph$flag_shift)[ph$flag_shift]
  conserv <- names(ph$flag_cons)[ph$flag_cons]
  do.call(rbind, lapply(names(kegg), function(nm) {
    g <- kegg[[nm]]; s <- st[g]
    data.frame(pair = paste(tA, tB, sep = "-"), pathway = nm,
               measured = length(g),
               n_maint = sum(s == "Maintained"), n_loss = sum(s == "Loss"),
               n_gain = sum(s == "Gain"),
               n_shift = sum(g %in% shifted), n_cons = sum(g %in% conserv),
               stringsAsFactors = FALSE)
  }))
}

t0 <- Sys.time()
res <- if (time.one) lapply(pairs, screen_pair) else
  mclapply(pairs, screen_pair, mc.cores = ncores, mc.preschedule = FALSE)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

if (time.one) {
  cat(sprintf("\none pair took %.1f s\n", elapsed))
  cat(sprintf("325 pairs serial: %.1f min; on %d cores: %.1f min\n",
              325 * elapsed / 60, ncores, 325 * elapsed / 60 / ncores))
  quit(save = "no")
}

bad <- !vapply(res, is.data.frame, logical(1))
if (any(bad)) cat("failed pairs:", sum(bad), "\n")
tab <- do.call(rbind, res[!bad])
tab$panel <- pmin(tab$n_shift, tab$n_cons)
saveRDS(tab, file.path(out.dir, "pair_pathway_counts.rds"))

best <- do.call(rbind, lapply(split(tab, tab$pair), function(d) {
  a <- d[which.max(d$panel), ]; b <- d[which.max(d$n_maint), ]
  data.frame(pair = a$pair,
             panel_pathway = a$pathway, panel = a$panel,
             panel_maint = a$n_maint, panel_shift = a$n_shift,
             panel_cons = a$n_cons, panel_loss = a$n_loss,
             densest_pathway = b$pathway, densest_maint = b$n_maint,
             densest_shift = b$n_shift, densest_cons = b$n_cons)
}))
best <- best[order(-best$panel, -best$panel_maint), ]
write.csv(best, file.path(out.dir, "pair_panel_best.csv"), row.names = FALSE)

cat(sprintf("\n%d pairs screened in %.1f min on %d cores\n",
            nrow(best), elapsed / 60, ncores))
cat("\n=== top 20 by min(shift, conserved) within a pathway ===\n")
print(head(best, 20), row.names = FALSE)
