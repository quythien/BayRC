## Rank every tissue pair by the genes that would appear in a pathway heatmap:
## maintained, and maintained plus the gain and loss classes.
##
## Usage: Rscript pair_maintained_rank.R [bfdr_alpha]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
suppressPackageStartupMessages(library(BayRC))

args  <- commandArgs(trailingOnly = TRUE)
alpha <- if (length(args) >= 1) as.numeric(args[1]) else 0.25

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
rho <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
tis <- names(rho)
p   <- lapply(rho, rowMeans)
cat("tissues:", length(tis), "\n")

cmb <- combn(tis, 2)
res <- do.call(rbind, apply(cmb, 2, function(z) {
  tr <- transition_classify(p[[z[1]]], p[[z[2]]], bfdr_alpha = alpha)
  s  <- tr$gain_loss_status
  data.frame(A = z[1], B = z[2],
             gain = sum(s == "Gain"), loss = sum(s == "Loss"),
             maint = sum(s == "Maintained"),
             panel = sum(s %in% c("Gain", "Loss", "Maintained")))
}))
res <- res[order(-res$maint), ]
write.csv(res, file.path(BAYRC_OUTPUT_DIR, "pair_maintained_rank.csv"),
          row.names = FALSE)

cat("\ntop 20 pairs by maintained genes:\n")
print(head(res, 20), row.names = FALSE)
cat("\nwhere the SUN pairs sit:\n")
sr <- res[res$A == "SUN" | res$B == "SUN", ]
sr$rank <- match(rownames(sr), rownames(res))
print(head(sr[, c("A","B","gain","loss","maint","panel","rank")], 8), row.names = FALSE)
