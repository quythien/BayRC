## Agreement between the Bayesian Bayes factor and the cosinor F-test in baboon
## lung, across raw and multiplicity-adjusted cosinor thresholds.

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
rho <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
post <- rowMeans(rho[["LUN"]])
BF   <- post / (1 - post + 1e-20)

source(file.path(BAYRC_PIPELINE_DIR, "one_cosinor_OLS_new.R"))
dat <- readRDS(file.path(BAYRC_DATA_DIR, "GTEXdata", "Baboon_processed.rds"))
mat <- dat$LUN$exprs; tod <- dat$LUN$tod
logmat <- log2(mat + 1)
logmat <- logmat[rownames(logmat) %in% names(BF), , drop = FALSE]

pv <- vapply(seq_len(nrow(logmat)), function(i) {
  y <- as.numeric(logmat[i, ])
  f <- try(one_cosinor_OLS(tod = tod, y = y, alpha = 0.05, period = 24,
                           CI = FALSE), silent = TRUE)
  if (inherits(f, "try-error")) return(NA_real_)
  n <- length(y); R2 <- f$R2
  pf(R2 / (1 - R2) * (n - 3) / 2, 2, n - 3, lower.tail = FALSE)
}, numeric(1))
names(pv) <- rownames(logmat)

g   <- intersect(names(BF), names(pv)[!is.na(pv)])
bay <- BF[g] > 3
p   <- pv[g]
q   <- p.adjust(p, "BH")
bon <- p.adjust(p, "bonferroni")

rules <- list("p < 0.05" = p < 0.05, "p < 0.01" = p < 0.01,
              "BH q < 0.05" = q < 0.05, "BH q < 0.01" = q < 0.01,
              "Bonferroni < 0.05" = bon < 0.05)

cat(sprintf("genes compared: %d   Bayesian BF > 3: %d\n\n", length(g), sum(bay)))
cat(sprintf("%-18s %7s %7s %9s %9s %8s\n", "cosinor rule", "cosinor",
            "both", "bayes only", "cos only", "neither"))
for (nm in names(rules)) {
  cos <- rules[[nm]]
  cat(sprintf("%-18s %7d %7d %9d %9d %8d\n", nm, sum(cos), sum(bay & cos),
              sum(bay & !cos), sum(!bay & cos), sum(!bay & !cos)))
}
