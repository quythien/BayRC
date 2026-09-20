# inst/analysis/config.R
# Central path configuration for all BayRC analysis scripts.
#
# Usage: source("config.R") at the top of each analysis script,
# then use the variables below instead of hardcoded paths.
#
# Override any path by setting the corresponding environment variable
# before launching R, e.g.:
#   export BAYRC_DATA_DIR=/path/to/your/data
#   Rscript inst/analysis/Baboon_SUN_PUT.R

BAYRC_DATA_DIR    <- Sys.getenv("BAYRC_DATA_DIR",
                       unset = "/home/qtp1/Projects/Collaborative")
BAYRC_WD_DIR      <- Sys.getenv("BAYRC_WD_DIR",
                       unset = "/home/qtp1/Projects/Circadian")
BAYRC_AGING_DIR   <- Sys.getenv("BAYRC_AGING_DIR",
                       unset = file.path(BAYRC_DATA_DIR,
                                         "Paper/Congruence/PNAS_aging"))
BAYRC_PACKAGE_DIR <- Sys.getenv("BAYRC_PACKAGE_DIR",
                       unset = file.path(BAYRC_WD_DIR,
                                         "Kyle/Circadian-analysis-main/R/v1/BayRC"))
BAYRC_PATHWAY_DIR <- Sys.getenv("BAYRC_PATHWAY_DIR",
                       unset = file.path(BAYRC_WD_DIR,
                                         "Kyle/Circadian-analysis-main/R/pathway_data"))
# The helper sources the analysis scripts pull in live beside the package, in
# R/v1/R/Thien, not inside it. Several scripts used to point at
# BAYRC_PACKAGE_DIR/Thien, which has never existed.
BAYRC_THIEN_DIR   <- Sys.getenv("BAYRC_THIEN_DIR",
                       unset = file.path(dirname(BAYRC_PACKAGE_DIR),
                                         "R", "Thien"))
BAYRC_GTEX_DIR    <- Sys.getenv("BAYRC_GTEX_DIR",
                       unset = file.path(BAYRC_DATA_DIR, "GTEXdata"))

# Per-tissue MCMC output, and the rho/phi summaries built from it by
# pipeline/summarize_rho_phi.R. Point BAYRC_RESULT_DIR at result_fixed/ to run
# the downstream analysis on the corrected sampler instead of the 2025 run.
BAYRC_RESULT_DIR  <- Sys.getenv("BAYRC_RESULT_DIR",
                       unset = file.path(BAYRC_GTEX_DIR, "result"))
BAYRC_SUMMARY_DIR <- Sys.getenv("BAYRC_SUMMARY_DIR",
                       unset = file.path(BAYRC_RESULT_DIR, "summary", "hb"))

# Where figures, tables and intermediate analysis output are written. The
# analysis scripts used to reassign output.dir several times per file; they now
# take it from here and build sub-directories underneath.
BAYRC_OUTPUT_DIR  <- Sys.getenv("BAYRC_OUTPUT_DIR",
                       unset = file.path(BAYRC_AGING_DIR, "results", "baboon",
                                         "output_final"))
BAYRC_FIGURE_DIR  <- Sys.getenv("BAYRC_FIGURE_DIR",
                       unset = file.path(BAYRC_AGING_DIR, "all_plots"))

# Validate that critical directories exist and warn if not
.check_dir <- function(path, name) {
  if (!dir.exists(path))
    warning("config.R: ", name, " does not exist: ", path,
            "\n  Set env var ", name, " to override.")
}
.check_dir(BAYRC_DATA_DIR,    "BAYRC_DATA_DIR")
.check_dir(BAYRC_WD_DIR,      "BAYRC_WD_DIR")
.check_dir(BAYRC_PATHWAY_DIR, "BAYRC_PATHWAY_DIR")
.check_dir(BAYRC_RESULT_DIR,  "BAYRC_RESULT_DIR")
rm(.check_dir)

# Output directories are created rather than warned about, since a fresh
# checkout will not have them yet.
for (.d in c(BAYRC_OUTPUT_DIR, BAYRC_FIGURE_DIR))
  dir.create(.d, recursive = TRUE, showWarnings = FALSE)
rm(.d)

# The summary directory decides which MCMC run every downstream number comes
# from, and nothing further down the pipeline names it. Report it and stop if
# the rho summary is not there, so a run against the wrong tree cannot pass
# unnoticed.
# A script that only arranges existing files sets bayrc.needs.summary <- FALSE
# before sourcing this, since it reports no numbers of its own.
.rho <- file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData")
.needs <- !exists("bayrc.needs.summary") || isTRUE(bayrc.needs.summary)
if (.needs && !file.exists(.rho))
  stop("config.R: no mcmc_rho_BF3.RData under BAYRC_SUMMARY_DIR: ",
       BAYRC_SUMMARY_DIR,
       "\n  Set BAYRC_RESULT_DIR (or BAYRC_SUMMARY_DIR) to the run you mean.")

message("BayRC config loaded. Override paths via environment variables (see config.R).")
message("  summaries: ", BAYRC_SUMMARY_DIR)
if (file.exists(.rho))
  message("  written:   ", format(file.info(.rho)$mtime, "%Y-%m-%d %H:%M"))
message("  output:    ", BAYRC_OUTPUT_DIR)
rm(.rho, .needs)
