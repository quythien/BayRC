# Posterior phase summaries for every measured member of the KEGG circadian
# pathway drawn in Figure 2C: the circular mean peak, the 95% circular HDI and
# the resultant length of each tissue-gene cell called rhythmic. The table is
# written beside the Figure 2 panels as pathway_phase_summary.csv, which
# explore_phase_groups.R reads to draw panel D.
suppressPackageStartupMessages(library(BayRC))

# Paths come from config.R; override any of them with the matching env var.
this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

fig2.dir <- file.path(BAYRC_FIGURE_DIR, "figure2")
membership <- file.path(fig2.dir, "Fig2C_circadian_membership.csv")
if (!file.exists(membership))
  stop("no ", membership, "; run plots/figure2_panelC.R first")
cells <- read.csv(membership)

message("Loading phase posterior")
env <- new.env()
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"), envir = env)
phi <- env$mcmc_phi_baboon

cells$peak <- cells$lower <- cells$upper <- cells$resultant <- NA_real_
for (i in seq_len(nrow(cells))) {
  if (!cells$called[i]) next
  x <- phi[[cells$tissue[i]]][cells$gene[i], ]
  x <- x[is.finite(x)] %% 24
  z <- mean(exp(2i * pi * x / 24))
  ci <- circular_HDI(x, credMass = .95, P = 24)
  cells$peak[i] <- (Arg(z) %% (2 * pi)) * 24 / (2 * pi)
  cells$lower[i] <- ci$lower
  cells$upper[i] <- ci$upper
  cells$resultant[i] <- Mod(z)
}

out <- file.path(fig2.dir, "pathway_phase_summary.csv")
write.csv(cells, out, row.names = FALSE)
message("Exported ", nrow(cells), " cells, ", sum(cells$called), " rhythmic, to ", out)
