## Compare thresholded rhythmic counts against the threshold-free posterior
## expectation for each tissue, to separate a real signal difference from the
## effect of the BFDR cut.
##
## Usage: Rscript sunput_rhythmic_budget.R [tissueA] [tissueB]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
suppressPackageStartupMessages(library(BayRC))

args <- commandArgs(trailingOnly = TRUE)
tA <- if (length(args) >= 1) args[1] else "PUT"
tB <- if (length(args) >= 2) args[2] else "SUN"

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
rho <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
p <- list(rowMeans(rho[[tA]]), rowMeans(rho[[tB]]))
names(p) <- c(tA, tB)

for (a in c(0.10, 0.25)) {
  tr <- transition_classify(p[[1]], p[[2]], bfdr_alpha = a)
  cat(sprintf("\nbfdr_alpha %.2f: gain %d  loss %d  maintained %d\n", a,
              sum(tr$gain_loss_status == "Gain"),
              sum(tr$gain_loss_status == "Loss"),
              sum(tr$gain_loss_status == "Maintained")))
}

cat("\n", sprintf("%-10s %10s %10s %8s\n", "", "expected", "at p>0.5", "ratio"))
for (t in names(p)) {
  cat(sprintf("%-10s %10.0f %10d %8.2f\n", t, sum(p[[t]]), sum(p[[t]] > 0.5),
              sum(p[[t]] > 0.5) / sum(p[[t]])))
}
cat(sprintf("\n%-10s %10.2f %10.2f\n", "B/A ratio",
            sum(p[[2]]) / sum(p[[1]]),
            sum(p[[2]] > 0.5) / sum(p[[1]] > 0.5)))

cat("\nposterior probability distribution, genes per band:\n")
br <- c(0, .1, .3, .5, .7, .9, 1)
tab <- sapply(p, function(x) table(cut(x, br, include.lowest = TRUE)))
print(tab)
cat("\nmean posterior probability:\n"); print(round(sapply(p, mean), 3))
