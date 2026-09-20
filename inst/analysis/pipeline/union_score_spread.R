## Spread of the Stage 1 union score for a set of tissue pairs. The union test
## ranks genes by evidence of rhythmicity in either condition, so a pair where
## nearly every gene is rhythmic has little to rank on.
##
## Usage: Rscript union_score_spread.R A1 B1 A2 B2 ...

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
rho <- get(grep("^mcmc_data", ls(), value = TRUE)[1])

args  <- commandArgs(trailingOnly = TRUE)
pairs <- split(args, ceiling(seq_along(args) / 2))

cat(sprintf("%-10s %8s %8s %8s %8s %8s %8s\n", "pair", "medS", "sdS",
            "IQR", "%sat", "rhyA", "rhyB"))
for (p in pairs) {
  pA <- rowMeans(rho[[p[1]]]); pB <- rowMeans(rho[[p[2]]])
  S  <- -2 * log(pmax((1 - pA) * (1 - pB), 1e-300))
  cat(sprintf("%-10s %8.2f %8.2f %8.2f %7.1f%% %8.0f %8.0f\n",
              paste(p, collapse = "-"), median(S), sd(S),
              IQR(S), 100 * mean(S > 10), sum(pA), sum(pB)))
}
cat("\nmedS/sdS: centre and spread of the Stage 1 ranking statistic\n")
cat("%sat: genes whose union score exceeds 10, i.e. near-certain rhythmicity\n")
