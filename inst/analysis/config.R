# inst/analysis/config.R
# Central path configuration for all BayRC analysis scripts.
#
# Usage: source("config.R") at the top of each analysis script,
# then use the variables below instead of hardcoded paths.
#
# Override any path by setting the corresponding environment variable
# before launching R, e.g.:
#   export BAYRC_DATA_DIR=/path/to/your/data
#   Rscript inst/analysis/applications/Baboon_PUT_SUN.R

# Every default is derived from where this file sits: inst/analysis inside the
# package checkout, which in turn sits at Kyle/Circadian-analysis-main/R/v1/BayRC
# under the Circadian project directory, with the collaborative data beside that
# directory. A checkout placed elsewhere sets the variables instead.
.config.file <- local({
  for (f in rev(sys.frames()))
    if (is.character(f$ofile)) return(normalizePath(f$ofile))
  file.path(getwd(), "config.R")
})
.up <- function(path, n) {
  for (i in seq_len(n)) path <- dirname(path)
  path
}
BAYRC_PACKAGE_DIR <- Sys.getenv("BAYRC_PACKAGE_DIR",
                       unset = .up(.config.file, 3))
BAYRC_WD_DIR      <- Sys.getenv("BAYRC_WD_DIR",
                       unset = .up(BAYRC_PACKAGE_DIR, 5))
BAYRC_DATA_DIR    <- Sys.getenv("BAYRC_DATA_DIR",
                       unset = file.path(dirname(BAYRC_WD_DIR), "Collaborative"))
BAYRC_AGING_DIR   <- Sys.getenv("BAYRC_AGING_DIR",
                       unset = file.path(BAYRC_DATA_DIR,
                                         "Paper/Congruence/PNAS_aging"))
# The pathway lists are read from R/pathway_data beside the checkout where that
# directory exists, and otherwise from the copies the package ships in extdata.
BAYRC_PATHWAY_DIR <- Sys.getenv("BAYRC_PATHWAY_DIR", unset = "")
if (!nzchar(BAYRC_PATHWAY_DIR)) {
  .pathway <- c(file.path(BAYRC_WD_DIR, "Kyle/Circadian-analysis-main/R/pathway_data"),
                file.path(BAYRC_PACKAGE_DIR, "inst", "extdata"),
                system.file("extdata", package = "BayRC"))
  BAYRC_PATHWAY_DIR <- .pathway[c(which(dir.exists(.pathway)), 1)[1]]
}
# The helper sources the analysis scripts pull in sit beside the package in some
# checkouts and inside it in others, so both are tried before either is assumed.
BAYRC_THIEN_DIR   <- Sys.getenv("BAYRC_THIEN_DIR", unset = "")
if (!nzchar(BAYRC_THIEN_DIR)) {
  .thien <- c(file.path(dirname(BAYRC_PACKAGE_DIR), "R", "Thien"),
              file.path(BAYRC_PACKAGE_DIR, "Thien"))
  BAYRC_THIEN_DIR <- if (any(dir.exists(.thien)))
    .thien[which(dir.exists(.thien))[1]] else .thien[1]
}
BAYRC_GTEX_DIR    <- Sys.getenv("BAYRC_GTEX_DIR",
                       unset = file.path(BAYRC_DATA_DIR, "GTEXdata"))
# one_cosinor_OLS_new.R lives beside the collaborative data, not in the package
BAYRC_PIPELINE_DIR <- Sys.getenv("BAYRC_PIPELINE_DIR",
                       unset = file.path(dirname(BAYRC_DATA_DIR), "Pipeline"))

# Per-tissue MCMC output, and the rho/phi summaries built from it by
# pipeline/summarize_rho_phi.R. The default is the 2025 run; point
# BAYRC_RESULT_DIR at result_fixed/ for the current run.
BAYRC_RESULT_DIR  <- Sys.getenv("BAYRC_RESULT_DIR",
                       unset = file.path(BAYRC_GTEX_DIR, "result"))
BAYRC_SUMMARY_DIR <- Sys.getenv("BAYRC_SUMMARY_DIR",
                       unset = file.path(BAYRC_RESULT_DIR, "summary", "hb"))

# Where figures, tables and intermediate output are written; scripts build
# sub-directories underneath.
BAYRC_OUTPUT_DIR  <- Sys.getenv("BAYRC_OUTPUT_DIR",
                       unset = file.path(BAYRC_AGING_DIR, "results", "baboon",
                                         "output_final"))
BAYRC_FIGURE_DIR  <- Sys.getenv("BAYRC_FIGURE_DIR",
                       unset = file.path(BAYRC_AGING_DIR, "all_plots"))

# A file outside the repository is reached through the variable that locates
# it, so a missing one names what to set.
bayrc_file <- function(dir, ...) {
  path <- file.path(dir, ...)
  if (!file.exists(path))
    stop("config.R: ", path, " does not exist.\n  Set env var ",
         deparse(substitute(dir)), " to the directory holding ",
         paste(..., collapse = "/"), ".", call. = FALSE)
  path
}

# Warn about missing input directories
.check_dir <- function(path, name) {
  if (!dir.exists(path))
    warning("config.R: ", name, " does not exist: ", path,
            "\n  Set env var ", name, " to override.")
}
.check_dir(BAYRC_DATA_DIR,    "BAYRC_DATA_DIR")
.check_dir(BAYRC_WD_DIR,      "BAYRC_WD_DIR")
.check_dir(BAYRC_PATHWAY_DIR, "BAYRC_PATHWAY_DIR")
.check_dir(BAYRC_RESULT_DIR,  "BAYRC_RESULT_DIR")
rm(.check_dir, .config.file, .up)
suppressWarnings(rm(.pathway, .thien))

# Output directories are created rather than warned about, since a fresh
# checkout will not have them yet.
for (.d in c(BAYRC_OUTPUT_DIR, BAYRC_FIGURE_DIR))
  dir.create(.d, recursive = TRUE, showWarnings = FALSE)
rm(.d)

# Report the summary directory every downstream number comes from, and stop if
# its rho summary is missing. A script that only arranges existing files sets
# bayrc.needs.summary <- FALSE before sourcing this.
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
