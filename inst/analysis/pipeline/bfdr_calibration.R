## Realised versus nominal BFDR for the current sampler, from simulated data
## with known rhythmic status.
##
## The design matches the baboon atlas: 12 samples at the observed times, and
## the sampler is called with the priors and bounds CAMO_h_b.R uses, so the
## posteriors carry the same granularity as the real run (2000 kept draws).
##
## Usage: Rscript bfdr_calibration.R <seed> [G] [outfile]

suppressPackageStartupMessages(library(BayRC))

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

args <- commandArgs(trailingOnly = TRUE)
seed <- if (length(args) >= 1) as.integer(args[1]) else 1L
G    <- if (length(args) >= 2) as.integer(args[2]) else 2000L
out.dir <- file.path(BAYRC_OUTPUT_DIR, "calibration")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
outfile <- if (length(args) >= 3) args[3] else
  file.path(out.dir, sprintf("bfdr_calibration_seed%d.rds", seed))

P <- 24; omega <- 2 * pi / P

## observed sampling times, falling back to an even grid
tod <- tryCatch({
  e <- new.env()
  load(file.path(BAYRC_GTEX_DIR, "data", "CAMO.bab.hum.RData"), envir = e)
  m <- e$baboon_withTOD$baboon$LUN
  cc <- grep("LUN[.]ZT", colnames(m))
  as.numeric(sub("LUN[.]ZT", "", colnames(m)[cc]))
}, error = function(e) seq(0, P, length.out = 13)[-13])
N <- length(tod)

## truth at the rate the sampler is told to expect
set.seed(1000 + seed)
truth <- rbinom(G, 1, 0.2)
Amp   <- ifelse(truth == 1, runif(G, 0.3, 1.2), 0)
phase <- runif(G, 0, P)
M     <- rnorm(G, 5, 1)
sigma <- 0.5
Y <- M + Amp * cos(omega * (matrix(tod, G, N, byrow = TRUE) - phase)) +
     matrix(rnorm(G * N, 0, sigma), G, N)
rownames(Y) <- paste0("G", seq_len(G))

dat <- list(data = as.data.frame(Y), time = tod, gname = rownames(Y))
init <- CBt_init_single(Data.list = dat, P = P, FitCosinor = TRUE,
                        mu_M = 0, sigma_M = 10, mu_A = 1, sigma_A = 10,
                        seed = seed)

t0 <- Sys.time()
res <- CB_MCMC_single_rj_slice(
  Data.list = dat, Init.value = init, P = P,
  iteration = 2500, thin = 1, n.burn = 500, seed = seed,
  p_rhythmic = rep(0.2, G), rj.p.stay = 0.5,
  A_prior = "trunc_Normal_OLS_condi",
  mu_A = 1, sigma_A = 10^2, A.min = 0,
  rj.phi = TRUE, rj.A = TRUE,
  mu_M = 0, sigma_M = 10^2,
  sigma_prior_v = 2, sigma_prior_s = 0,
  diagnostics = FALSE)
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

post <- rowMeans(res$rho)
nominal <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40)
rows <- lapply(nominal, function(a) {
  tau <- bfdr_from_posterior(post, alpha = a)$threshold
  called <- post >= tau
  TP <- sum(called & truth == 1); FP <- sum(called & truth == 0)
  data.frame(seed = seed, G = G, N = N, nominal = a,
             n_called = sum(called), TP = TP, FP = FP,
             realised = if (TP + FP > 0) FP / (TP + FP) else NA_real_,
             power = TP / sum(truth), threshold = tau)
})
d <- do.call(rbind, rows)

saveRDS(list(grid = d, post = post, truth = truth, amp = Amp,
             kept = ncol(res$rho), minutes = mins,
             med_halfrange = median(apply(Y, 1, function(x) (max(x) - min(x)) / 2))),
        outfile)

cat(sprintf("seed %d | G %d | N %d | kept %d | true rhythmic %d | %.1f min\n",
            seed, G, N, ncol(res$rho), sum(truth), mins))
print(d, row.names = FALSE)
