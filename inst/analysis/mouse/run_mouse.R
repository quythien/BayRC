## Fixed-arm BayRC MCMC for one mouse tissue of the GSE54651 atlas.
##
## The mouse companion to run_fixed.R. Same sampler, same settings, same
## thin = 1 so the stored chain has 2001 columns, and output goes under
## result_fixed/ beside the human and baboon runs.

args    <- commandArgs(trailingOnly = TRUE)
tissue  <- args[1]
outroot <- args[2]
libp    <- args[3]
seed    <- if (length(args) >= 4) as.integer(args[4]) else 1L
n.iter  <- if (length(args) >= 5) as.integer(args[5]) else 2500L
n.burn  <- if (length(args) >= 6) as.integer(args[6]) else 500L

.libPaths(c(libp, .libPaths()))
suppressPackageStartupMessages(library(BayRC))

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
bayrc.needs.summary <- FALSE
source(file.path(analysis.dir, "config.R"))

# CAMO.mouse.hum.RData sits under BAYRC_GTEX_DIR/data, beside CAMO.bab.hum.RData
datafile <- file.path(BAYRC_GTEX_DIR, "data", "CAMO.mouse.hum.RData")
if (!file.exists(datafile))
  stop("set BAYRC_GTEX_DIR to the directory holding data/CAMO.mouse.hum.RData")
load(datafile)

## mice$count_clean is log2 CPM keyed on human Ensembl identifiers, so it is
## used as-is, the same way run_fixed.R uses gtex$CPM.large.clean.
Y   <- as.matrix(mice$count_clean[[tissue]])
tod <- mice$tod[[tissue]]
stopifnot(min(Y) < 0)
stopifnot(ncol(Y) == length(tod))

cat("job: mouse", tissue, "on", Sys.info()[["nodename"]], "\n")
cat("input:", datafile, "written",
    format(file.info(datafile)$mtime, "%Y-%m-%d"), "\n")
cat("genes:", nrow(Y), " samples:", ncol(Y),
    " tod:", paste(tod, collapse = ", "), "\n")

dat.input <- list(data = as.data.frame(Y), time = tod, gname = rownames(Y))

a.init <- CBt_init_single(Data.list = dat.input, P = 24, FitCosinor = TRUE,
                          mu_M = 0, sigma_M = 10, mu_A = 1, sigma_A = 10,
                          seed = seed)

t0 <- Sys.time()
CB.res <- CB_MCMC_single_rj_slice(
  Data.list  = dat.input,
  Init.value = a.init,
  P          = 24,
  iteration  = n.iter,
  thin       = 1,
  n.burn     = n.burn,
  seed       = seed,
  diagnostics = FALSE,
  p_rhythmic = rep(0.2, nrow(Y)),
  rj.p.stay  = 0.5,
  A_prior    = "trunc_Normal_OLS_condi",
  mu_A = 1, sigma_A = 10^2, A.min = 0,
  rj.phi = TRUE, rj.A = TRUE,
  mu_M = 0, sigma_M = 10^2,
  sigma_prior_v = 2, sigma_prior_s = 0
)
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

if (n.iter - n.burn == 2000L) stopifnot(ncol(CB.res$rho) == 2001L)
stopifnot(sum(is.na(CB.res$rho)) == 0L)

save_dir <- file.path(outroot, "mice", tissue)
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
save_file <- file.path(save_dir, paste0("mice_", tissue, "_bay_", seed, ".RDS"))
saveRDS(CB.res, file = save_file)
if (!file.exists(save_file)) stop("Error: File not saved correctly")

## The record beside the output names the script, the input with its date and
## the settings, so every number taken from this chain can be traced back.
record <- c(
  paste("script       ", "src/run_mouse.R"),
  paste("mode         ", "full run"),
  paste("BayRC version", as.character(utils::packageVersion("BayRC"))),
  paste("species      ", "mouse"),
  paste("tissue       ", tissue),
  paste("input        ", datafile),
  paste("input written", format(file.info(datafile)$mtime, "%Y-%m-%d %H:%M")),
  paste("output       ", save_file),
  paste("run at       ", format(Sys.time(), "%Y-%m-%d %H:%M")),
  paste("node         ", Sys.info()[["nodename"]]),
  paste("minutes      ", sprintf("%.1f", mins)),
  "parameters",
  sprintf("  %-14s %s", "genes", nrow(Y)),
  sprintf("  %-14s %s", "samples", ncol(Y)),
  sprintf("  %-14s %s", "tod", paste(tod, collapse = ", ")),
  sprintf("  %-14s %s", "transform", "none, mice$count_clean is log2 CPM"),
  sprintf("  %-14s %s", "iteration", n.iter),
  sprintf("  %-14s %s", "n.burn", n.burn),
  sprintf("  %-14s %s", "thin", 1),
  sprintf("  %-14s %s", "seed", seed),
  sprintf("  %-14s %s", "P", 24),
  sprintf("  %-14s %s", "p_rhythmic", 0.2),
  sprintf("  %-14s %s", "rj.p.stay", 0.5),
  sprintf("  %-14s %s", "A_prior", "trunc_Normal_OLS_condi"),
  sprintf("  %-14s %s", "A.max", "per gene, (max-min)/2"),
  sprintf("  %-14s %s", "mu_A", 1),
  sprintf("  %-14s %s", "sigma_A", "10^2"),
  sprintf("  %-14s %s", "mu_M", 0),
  sprintf("  %-14s %s", "sigma_M", "10^2"),
  sprintf("  %-14s %s", "sigma_prior_v", 2),
  sprintf("  %-14s %s", "sigma_prior_s", 0))
writeLines(record, file.path(save_dir, "run_record.txt"))

cat(sprintf("DONE mouse %s G=%d N=%d ncol=%d NA=%d mins=%.1f meanPost=%.4f\n",
            tissue, nrow(Y), ncol(Y), ncol(CB.res$rho),
            sum(is.na(CB.res$rho)), mins, mean(rowMeans(CB.res$rho))))
