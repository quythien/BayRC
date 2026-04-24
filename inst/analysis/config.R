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
BAYRC_THIEN_DIR   <- Sys.getenv("BAYRC_THIEN_DIR",
                       unset = file.path(BAYRC_PACKAGE_DIR, "Thien"))
BAYRC_GTEX_DIR    <- Sys.getenv("BAYRC_GTEX_DIR",
                       unset = file.path(BAYRC_DATA_DIR, "GTEXdata"))

# Validate that critical directories exist and warn if not
.check_dir <- function(path, name) {
  if (!dir.exists(path))
    warning("config.R: ", name, " does not exist: ", path,
            "\n  Set env var ", name, " to override.")
}
.check_dir(BAYRC_DATA_DIR,    "BAYRC_DATA_DIR")
.check_dir(BAYRC_WD_DIR,      "BAYRC_WD_DIR")
.check_dir(BAYRC_PATHWAY_DIR, "BAYRC_PATHWAY_DIR")
rm(.check_dir)

message("BayRC config loaded. Override paths via environment variables (see config.R).")
