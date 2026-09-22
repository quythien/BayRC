## BayRC MCMC for one human or baboon tissue of the matched atlas.
##
## One tissue per process, so a driver can fan the 52 runs out across cores.
## thin = 1, so a 2,500 iteration run with a 500 iteration burn-in stores
## 2,001 posterior samples per gene. Output goes to
## <outroot>/<species>/<tissue>, which the summary step reads.
##
##   Rscript run_fixed.R <species> <tissue> <outroot> <libpath> [seed] [iter] [burn]
##
## mouse/run_mouse.R is the mouse counterpart and takes the same settings.

args    <- commandArgs(trailingOnly = TRUE)
species <- args[1]                       # "human" | "baboon"
tissue  <- args[2]
outroot <- args[3]
libp    <- args[4]
seed    <- if (length(args) >= 5) as.integer(args[5]) else 1L
n.iter  <- if (length(args) >= 6) as.integer(args[6]) else 2500L
n.burn  <- if (length(args) >= 7) as.integer(args[7]) else 500L

if (length(args) < 4)
  stop("usage: run_fixed.R <species> <tissue> <outroot> <libpath> [seed] [iter] [burn]")
if (!species %in% c("human", "baboon"))
  stop("species must be 'human' or 'baboon', not ", species)

.libPaths(c(libp, .libPaths()))
suppressPackageStartupMessages(library(BayRC))

# CAMO.bab.hum.RData sits under BAYRC_GTEX_DIR/data
datafile <- file.path(Sys.getenv("BAYRC_GTEX_DIR"), "data", "CAMO.bab.hum.RData")
if (!nzchar(Sys.getenv("BAYRC_GTEX_DIR")) || !file.exists(datafile))
  stop("set BAYRC_GTEX_DIR to the directory holding data/CAMO.bab.hum.RData")
load(datafile)

if (species == "human") {
  ## gtex$CPM.large.clean is already on the log scale, so it is used as-is.
  Y   <- as.matrix(gtex$CPM.large.clean[[tissue]])
  tod <- gtex$tod[[tissue]]
  stopifnot(min(Y) < 0)
  prefix <- "gtex"
} else {
  ## the first two baboon columns carry identifiers rather than samples
  raw <- baboon_withTOD$baboon[[tissue]][, -(1:2)]
  Y   <- log2(as.matrix(raw) + 1)
  tod <- baboon_withTOD$tod[[tissue]]
  stopifnot(min(Y) >= 0)
  prefix <- "baboon"
}
stopifnot(ncol(Y) == length(tod))

cat("job:", species, tissue, "on", Sys.info()[["nodename"]], "\n")
cat("input:", datafile, "written",
    format(file.info(datafile)$mtime, "%Y-%m-%d"), "\n")
cat("genes:", nrow(Y), " samples:", ncol(Y), "\n")

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

save_dir <- file.path(outroot, species, tissue)
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
save_file <- file.path(save_dir, paste0(prefix, "_", tissue, "_bay_", seed, ".RDS"))
saveRDS(CB.res, file = save_file)
if (!file.exists(save_file)) stop("chain was not written to ", save_file)

cat(sprintf("DONE %s %s G=%d N=%d ncol=%d NA=%d mins=%.1f meanPost=%.4f\n",
            species, tissue, nrow(Y), ncol(Y), ncol(CB.res$rho),
            sum(is.na(CB.res$rho)), mins, mean(rowMeans(CB.res$rho))))
