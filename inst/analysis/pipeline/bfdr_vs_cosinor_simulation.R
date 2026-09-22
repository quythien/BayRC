## BFDR against cosinor regression on the calibration data sets.
##
## Regenerates each data set that bfdr_calibration.R simulated, with the same
## generator and seed, and fits an ordinary least-squares cosinor with an F-test
## per gene. The Bayesian posteriors are the ones bfdr_calibration.R stored, so
## no sampler is re-run. Both rules are then scored against the known truth:
##
##   BFDR at alpha          realised false discovery rate and power
##   cosinor + BH at alpha  the same two, at the same nominal level
##   cosinor p < alpha      type I error on the non-rhythmic genes, and power
##
## Realised FDR is pooled over replicates as total false calls over total calls;
## power and type I error are means over replicates. Run after
## bfdr_calibration.R has written its ten seeds.

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
bayrc.needs.summary <- FALSE
source(file.path(analysis.dir, "config.R"))
suppressPackageStartupMessages(library(BayRC))

cal.dir <- file.path(BAYRC_OUTPUT_DIR, "calibration")
seeds   <- as.integer(sub(".*seed([0-9]+)[.]rds$", "\\1",
                      list.files(cal.dir, "^bfdr_calibration_seed[0-9]+[.]rds$")))
alphas  <- c(0.05, 0.10, 0.20, 0.25)
P <- 24; omega <- 2 * pi / P; G <- 2000L

## the calibration run's own sampling times
tod <- tryCatch({
  e <- new.env()
  load(file.path(BAYRC_GTEX_DIR, "data", "CAMO.bab.hum.RData"), envir = e)
  m <- e$baboon_withTOD$baboon$LUN
  cc <- grep("LUN[.]ZT", colnames(m))
  as.numeric(sub("LUN[.]ZT", "", colnames(m)[cc]))
}, error = function(e) seq(0, P, length.out = 13)[-13])
N <- length(tod)

## one design matrix serves every gene
X    <- cbind(1, cos(omega * tod), sin(omega * tod))
hat  <- X %*% solve(crossprod(X), t(X))
f_p  <- function(Y) {
  rss1 <- rowSums((Y - Y %*% t(hat))^2)
  rss0 <- rowSums((Y - rowMeans(Y))^2)
  Fst  <- ((rss0 - rss1) / 2) / (rss1 / (N - 3))
  pf(Fst, 2, N - 3, lower.tail = FALSE)
}

score <- function(call, truth) c(calls = sum(call), false = sum(call & truth == 0),
                                 power = mean(call[truth == 1]))

rows <- list()
for (s in seeds) {
  cal <- readRDS(file.path(cal.dir, sprintf("bfdr_calibration_seed%d.rds", s)))

  ## identical generator, identical draw order
  set.seed(1000 + s)
  truth <- rbinom(G, 1, 0.2)
  Amp   <- ifelse(truth == 1, runif(G, 0.3, 1.2), 0)
  phase <- runif(G, 0, P)
  M     <- rnorm(G, 5, 1)
  Y <- M + Amp * cos(omega * (matrix(tod, G, N, byrow = TRUE) - phase)) +
       matrix(rnorm(G * N, 0, 0.5), G, N)
  if (!identical(as.integer(truth), as.integer(cal$truth)))
    stop("seed ", s, ": regenerated truth does not match the stored truth")

  p <- f_p(Y)
  for (a in alphas) {
    bcall <- bfdr_from_posterior(cal$post, alpha = a)$rhythmic_genes
    b  <- score(bcall, truth)
    bh <- score(p.adjust(p, "BH") <= a, truth)
    un <- score(p <= a, truth)
    rows[[length(rows) + 1]] <- data.frame(
      seed = s, alpha = a,
      bfdr_calls = b[["calls"]], bfdr_false = b[["false"]], bfdr_power = b[["power"]],
      bfdr_type1 = mean(bcall[truth == 0]),
      bh_calls   = bh[["calls"]], bh_false  = bh[["false"]], bh_power  = bh[["power"]],
      raw_power  = un[["power"]], raw_type1 = mean((p <= a)[truth == 0]))
  }
}
d <- do.call(rbind, rows)

pooled <- do.call(rbind, lapply(split(d, d$alpha), function(x) data.frame(
  alpha        = x$alpha[1],
  bfdr_fdr     = sum(x$bfdr_false) / sum(x$bfdr_calls),
  bfdr_power   = mean(x$bfdr_power),
  bfdr_type1   = mean(x$bfdr_type1),
  bh_fdr       = sum(x$bh_false) / max(1, sum(x$bh_calls)),
  bh_power     = mean(x$bh_power),
  raw_type1    = mean(x$raw_type1),
  raw_power    = mean(x$raw_power))))

cat(sprintf("%d replicates of %d genes at n = %d; regenerated truth matched in every one\n\n",
            length(seeds), G, N))
print(format(pooled, digits = 3), row.names = FALSE)
write.csv(pooled, file.path(cal.dir, "bfdr_vs_cosinor_summary.csv"), row.names = FALSE)
