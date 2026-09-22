## Realised versus nominal BFDR for the current sampler, from simulated data
## with known rhythmic status.
##
## The design matches the baboon atlas: 12 samples at the observed times, and
## the sampler is called with the priors and bounds CAMO_h_b.R uses, so the
## posteriors carry the same granularity as the real run (2000 kept draws).
##
## Usage: Rscript bfdr_calibration.R <seed> [G] [outfile] [N] [A] [sigma]
##
## With no optional arguments this is the atlas design: the observed lung times,
## amplitudes drawn on [0.3, 1.2] and residual SD 0.5. N replaces the times with N
## samples evenly spaced over one period; A gives every rhythmic gene that one
## amplitude, and sigma sets the residual SD, so A / sigma can be swept.

suppressPackageStartupMessages(library(BayRC))

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

args <- commandArgs(trailingOnly = TRUE)

## --summary pools the replicates already on disk instead of running another
if ("--summary" %in% args) {
  out.dir <- file.path(BAYRC_OUTPUT_DIR, "calibration")
  f <- list.files(out.dir, "^bfdr_calibration_seed.*[.]rds$", full.names = TRUE)
  if (!length(f)) stop("no replicates under ", out.dir)
  reps <- lapply(f, readRDS)
  g <- do.call(rbind, lapply(reps, `[[`, "grid"))
  ## pooled over replicates: total false calls over total calls
  pooled <- do.call(rbind, lapply(split(g, g$nominal), function(d)
    data.frame(nominal = d$nominal[1],
               replicates = nrow(d),
               mean_called = round(mean(d$n_called), 1),
               realised = round(sum(d$FP) / sum(d$TP + d$FP), 4),
               se = round(sd(d$realised) / sqrt(nrow(d)), 4),
               min = round(min(d$realised), 3),
               max = round(max(d$realised), 3),
               ratio = round((sum(d$FP) / sum(d$TP + d$FP)) / d$nominal[1], 2),
               power = round(mean(d$power), 3))))
  cat(sprintf("\n%d replicates, G = %d, N = %d, %d kept draws\n",
              length(reps), reps[[1]]$grid$G[1], reps[[1]]$grid$N[1],
              reps[[1]]$kept))
  cat(sprintf("median simulated (max-min)/2 = %.2f\n",
              median(sapply(reps, `[[`, "med_halfrange"))))
  print(pooled, row.names = FALSE)
  write.csv(pooled, file.path(out.dir, "bfdr_calibration_summary.csv"),
            row.names = FALSE)
  quit(save = "no")
}

seed <- if (length(args) >= 1) as.integer(args[1]) else 1L
G    <- if (length(args) >= 2) as.integer(args[2]) else 2000L
N_arg <- if (length(args) >= 4 && args[4] != "-") as.integer(args[4]) else NA_integer_
A_arg <- if (length(args) >= 5 && args[5] != "-") as.numeric(args[5]) else NA_real_
S_arg <- if (length(args) >= 6 && args[6] != "-") as.numeric(args[6]) else NA_real_
out.dir <- file.path(BAYRC_OUTPUT_DIR, "calibration")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
outfile <- if (length(args) >= 3 && nzchar(args[3]) && args[3] != "-") args[3] else
  file.path(out.dir, if (is.na(N_arg) && is.na(A_arg)) sprintf("bfdr_calibration_seed%d.rds", seed)
                     else sprintf("bfdr_calibration_A%g_s%g_n%d_seed%d.rds",
                                  A_arg, if (is.na(S_arg)) 0.5 else S_arg,
                                  if (is.na(N_arg)) 12L else N_arg, seed))

P <- 24; omega <- 2 * pi / P

## observed sampling times, falling back to an even grid
tod <- tryCatch({
  e <- new.env()
  load(file.path(BAYRC_GTEX_DIR, "data", "CAMO.bab.hum.RData"), envir = e)
  m <- e$baboon_withTOD$baboon$LUN
  cc <- grep("LUN[.]ZT", colnames(m))
  as.numeric(sub("LUN[.]ZT", "", colnames(m)[cc]))
}, error = function(e) seq(0, P, length.out = 13)[-13])
if (!is.na(N_arg)) tod <- seq(0, P, length.out = N_arg + 1)[-(N_arg + 1)]
N <- length(tod)

## truth at the rate the sampler is told to expect
set.seed(1000 + seed)
truth <- rbinom(G, 1, 0.2)
Amp   <- ifelse(truth == 1, runif(G, 0.3, 1.2), 0)
if (!is.na(A_arg)) Amp <- ifelse(truth == 1, A_arg, 0)
phase <- runif(G, 0, P)
M     <- rnorm(G, 5, 1)
if (!is.na(A_arg)) M <- rep(5, G)   # the swept design fixes the mesor
sigma <- if (is.na(S_arg)) 0.5 else S_arg
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
