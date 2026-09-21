## Agreement between the Bayesian Bayes factor and the cosinor F-test in baboon
## lung, across raw and multiplicity-adjusted cosinor thresholds. Follows the
## same fit as plots/S5_Bayes_Cosinor_Agreement_LUN.R.

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

BAYRC_PIPELINE_DIR <- Sys.getenv("BAYRC_PIPELINE_DIR",
                        unset = file.path(dirname(BAYRC_DATA_DIR), "Pipeline"))
source(file.path(BAYRC_PIPELINE_DIR, "one_cosinor_OLS_new.R"))

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
rho  <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
post <- rowMeans(rho[["LUN"]])
BF   <- post / (1 - post + 1e-20)

to_zt <- function(t_cos) ifelse(t_cos >= 18, t_cos - 24, t_cos)

load(file.path(BAYRC_GTEX_DIR, "data", "CAMO.bab.hum.RData"))
bab_raw  <- baboon_withTOD$baboon$LUN
bab_cols <- grep("LUN\\.ZT", colnames(bab_raw))
bab_mat  <- as.matrix(bab_raw[, bab_cols])
tod      <- to_zt(as.numeric(sub("LUN\\.ZT", "", colnames(bab_mat))))
symbols  <- bab_raw$Symbol
symbols[symbols == "ARNTL"] <- "BMAL1"

tpm    <- sweep(bab_mat, 2, colSums(bab_mat), FUN = "/") * 1e6
logmat <- log2(tpm + 1)
n      <- length(tod)

pv <- vapply(seq_len(nrow(logmat)), function(i) {
  fit <- tryCatch(one_cosinor_OLS(tod = tod, y = as.numeric(logmat[i, ]),
                                  alpha = 0.05, period = 24, CI = FALSE),
                  error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  R2 <- fit$test$R2
  pf(R2 * (n - 3) / (2 * (1 - R2)), 2, n - 3, lower.tail = FALSE)
}, numeric(1))
names(pv) <- symbols

pv  <- pv[!is.na(pv) & !duplicated(names(pv))]
g   <- intersect(names(BF), names(pv))
bay <- BF[g] > 3
p   <- pv[g]

rules <- list("p < 0.05"          = p < 0.05,
              "p < 0.01"          = p < 0.01,
              "BH q < 0.05"       = p.adjust(p, "BH") < 0.05,
              "BH q < 0.01"       = p.adjust(p, "BH") < 0.01,
              "Bonferroni < 0.05" = p.adjust(p, "bonferroni") < 0.05)

cat(sprintf("\ngenes compared: %d   Bayesian BF > 3: %d\n\n", length(g), sum(bay)))
cat(sprintf("%-18s %8s %7s %11s %10s %8s\n", "cosinor rule", "cosinor",
            "both", "bayes only", "cos only", "neither"))
for (nm in names(rules)) {
  cos <- rules[[nm]]
  cat(sprintf("%-18s %8d %7d %11d %10d %8d\n", nm, sum(cos), sum(bay & cos),
              sum(bay & !cos), sum(!bay & cos), sum(!bay & !cos)))
}
