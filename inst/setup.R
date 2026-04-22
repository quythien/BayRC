# BayRC Setup Script
# Run this once to install all required and suggested packages before using BayRC.

# ── CRAN packages ─────────────────────────────────────────────────────────────
cran_pkgs <- c(
  # Required (Imports)
  "Rcpp",
  "circular",
  "ggplot2",
  "dplyr",
  "openxlsx",
  # Suggested — optional but needed for full functionality
  "parallel",
  "tibble",
  "purrr",
  "RColorBrewer",
  "edgeR",
  "knitr",
  "rmarkdown",
  "testthat"
)

install_if_missing <- function(pkgs) {
  to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(to_install) > 0) {
    message("Installing: ", paste(to_install, collapse = ", "))
    install.packages(to_install)
  } else {
    message("All CRAN packages already installed.")
  }
}

install_if_missing(cran_pkgs)

# ── Bioconductor packages ─────────────────────────────────────────────────────
bioc_pkgs <- c(
  "fgsea",       # pathway enrichment (used by pathSelect)
  "DESeq2",      # optional: differential expression pre-processing
  "biomaRt",     # optional: gene annotation and ortholog mapping
  "KEGGREST"     # optional: KEGG pathway retrieval
)

bioc_missing <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(bioc_missing) > 0) {
  message("Installing Bioconductor packages: ", paste(bioc_missing, collapse = ", "))
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(bioc_missing)
} else {
  message("All Bioconductor packages already installed.")
}

# ── Development tools (for loading BayRC from source) ─────────────────────────
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

# ── Load BayRC ────────────────────────────────────────────────────────────────
# Option 1: installed package
# library(BayRC)

# Option 2: load from source during development
BAYRC_ROOT <- "/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/BayRC"
devtools::load_all(BAYRC_ROOT)
message("BayRC loaded successfully.")
