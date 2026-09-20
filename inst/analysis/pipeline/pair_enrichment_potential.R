## Rank all tissue pairs on the two things that predict a usable pathway
## result: enough conserved genes to fill a panel, and enough contrast in the
## Stage 1 ranking for a pathway to separate from background.

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
rho <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
p   <- lapply(rho, rowMeans)

rank.csv <- file.path(BAYRC_OUTPUT_DIR, "pair_maintained_rank.csv")
mr <- read.csv(rank.csv)

sat <- mapply(function(a, b) {
  S <- -2 * log(pmax((1 - p[[a]]) * (1 - p[[b]]), 1e-300))
  100 * mean(S > 10)
}, mr$A, mr$B)

mr$sat <- round(sat, 1)
## pairs that enriched in the scan sat between roughly 2 and 12 percent
mr$band <- mr$sat >= 2 & mr$sat <= 12
out <- mr[mr$band & mr$maint >= 400, ]
out <- out[order(-out$maint), ]
write.csv(out, file.path(BAYRC_OUTPUT_DIR, "pair_enrichment_potential.csv"),
          row.names = FALSE)

cat("pairs in the contrast band with at least 400 conserved genes:",
    nrow(out), "of", nrow(mr), "\n\n")
print(head(out[, c("A", "B", "gain", "loss", "maint", "sat")], 30),
      row.names = FALSE)
