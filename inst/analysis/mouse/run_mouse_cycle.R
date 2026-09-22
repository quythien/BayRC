## One mouse tissue of the GSE54651 atlas, restricted to a single circadian
## cycle. The eight samples are four phases (CT 22, 4, 10, 16) measured twice,
## so cycle 1 is CT 22 to 40 and cycle 2 is CT 46 to 64. Each arm therefore fits
## four observations, which leaves one residual degree of freedom; the arms are
## run to compare the two cycles against each other and against the pooled fit.

args    <- commandArgs(trailingOnly = TRUE)
tissue  <- args[1]
cycle   <- as.integer(args[2])
outroot <- args[3]
libp    <- args[4]
seed    <- if (length(args) >= 5) as.integer(args[5]) else 1L
n.iter  <- if (length(args) >= 6) as.integer(args[6]) else 2500L
n.burn  <- if (length(args) >= 7) as.integer(args[7]) else 500L

.libPaths(c(libp, .libPaths()))
suppressPackageStartupMessages(library(BayRC))

# CAMO.mouse.hum.RData sits under BAYRC_GTEX_DIR/data, beside CAMO.bab.hum.RData
datafile <- file.path(Sys.getenv("BAYRC_GTEX_DIR"), "data", "CAMO.mouse.hum.RData")
if (!nzchar(Sys.getenv("BAYRC_GTEX_DIR")) || !file.exists(datafile))
  stop("set BAYRC_GTEX_DIR to the directory holding data/CAMO.mouse.hum.RData")
load(datafile)

Y   <- as.matrix(mice$count_clean[[tissue]])
tod <- as.numeric(mice$tod[[tissue]])
stopifnot(ncol(Y) == length(tod))

keep <- if (cycle == 1L) tod < 46 else tod >= 46
Y   <- Y[, keep, drop = FALSE]
tod <- tod[keep]
stopifnot(ncol(Y) == 4L, length(unique(tod %% 24)) == 4L)

cat("job: mouse", tissue, "cycle", cycle, "on", Sys.info()[["nodename"]], "\n")
cat("genes:", nrow(Y), " samples:", ncol(Y),
    " tod:", paste(tod, collapse = ", "), "\n")

dat.input <- list(data = as.data.frame(Y), time = tod, gname = rownames(Y))
a.init <- CBt_init_single(Data.list = dat.input, P = 24, FitCosinor = TRUE,
                          mu_M = 0, sigma_M = 10, mu_A = 1, sigma_A = 10,
                          seed = seed)
t0 <- Sys.time()
CB.res <- CB_MCMC_single_rj_slice(
  Data.list = dat.input, Init.value = a.init, P = 24,
  iteration = n.iter, thin = 1, n.burn = n.burn, seed = seed,
  diagnostics = FALSE, p_rhythmic = rep(0.2, nrow(Y)), rj.p.stay = 0.5,
  A_prior = "trunc_Normal_OLS_condi", mu_A = 1, sigma_A = 10^2, A.min = 0,
  rj.phi = TRUE, rj.A = TRUE, mu_M = 0, sigma_M = 10^2,
  sigma_prior_v = 2, sigma_prior_s = 0)
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

if (n.iter - n.burn == 2000L) stopifnot(ncol(CB.res$rho) == 2001L)
stopifnot(sum(is.na(CB.res$rho)) == 0L)

save_dir <- file.path(outroot, "mice_cycle", paste0("cycle", cycle), tissue)
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
save_file <- file.path(save_dir, paste0("mice_", tissue, "_cycle", cycle, "_bay_", seed, ".RDS"))
saveRDS(CB.res, file = save_file)

writeLines(c(
  paste("script        src/run_mouse_cycle.R"),
  paste("tissue       ", tissue),
  paste("cycle        ", cycle, "of 2"),
  paste("input        ", datafile, "written",
        format(file.info(datafile)$mtime, "%Y-%m-%d")),
  paste("samples      ", ncol(Y), "at CT", paste(tod, collapse = " ")),
  paste("genes        ", nrow(Y)),
  paste("iteration    ", n.iter, " n.burn", n.burn, " thin 1  seed", seed),
  paste("A_prior       trunc_Normal_OLS_condi   p_rhythmic 0.2"),
  paste("run at       ", format(Sys.time(), "%Y-%m-%d %H:%M"))),
  file.path(save_dir, "run_record.txt"))

cat(sprintf("DONE mouse %s cycle %d G=%d N=%d ncol=%d NA=%d mins=%.1f meanPost=%.4f\n",
    tissue, cycle, nrow(Y), ncol(Y), ncol(CB.res$rho),
    sum(is.na(CB.res$rho)), mins, mean(rowMeans(CB.res$rho))))
