## Extract the SUN and PUT rho/phi slices into a small cache the parameter
## search can reload without reading the full genome-wide summaries.
##
## Usage: Rscript sunput_cache.R [tissueA] [tissueB]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

args <- commandArgs(trailingOnly = TRUE)
tA <- if (length(args) >= 1) args[1] else "PUT"
tB <- if (length(args) >= 2) args[2] else "SUN"

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))

cache <- list(
  A = list(rho = mcmc_data_baboon[[tA]], phi = mcmc_phi_baboon[[tA]]),
  B = list(rho = mcmc_data_baboon[[tB]], phi = mcmc_phi_baboon[[tB]]),
  names = c(tA, tB),
  summary.dir = BAYRC_SUMMARY_DIR)

out.dir <- file.path(BAYRC_OUTPUT_DIR, "param_search")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out.dir, sprintf("cache_%s_%s.rds", tA, tB))
saveRDS(cache, out)
cat("wrote", out, "\n")
cat("genes:", nrow(cache$A$rho), " iterations:", ncol(cache$A$rho), "\n")
