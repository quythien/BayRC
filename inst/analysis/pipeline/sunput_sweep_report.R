## Print the pathways selected by named cells of the sunput_enrich_sweep grid.
##
## Usage: Rscript sunput_sweep_report.R ["list|size|stage1|qcut" ...]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

hits <- readRDS(file.path(BAYRC_OUTPUT_DIR, "param_search",
                          "enrich_sweep_hits.rds"))

cells <- commandArgs(trailingOnly = TRUE)
if (!length(cells)) cells <- c("kegg354|raw5-500|q<0.05|0.2",
                               "kegg354|raw5-500|p<0.05|0.2",
                               "kegg229|raw5-500|p<0.05|0.2")

for (k in cells) {
  if (is.null(hits[[k]])) {
    cat("\nno such cell:", k, "\n")
    next
  }
  cat("\n=====", k, "=====\n")
  for (rk in c("gain", "loss", "conserved")) {
    p <- hits[[k]][[rk]]
    cat(sprintf("-- %s (%d)\n", rk, length(p)))
    if (length(p)) cat(paste0("   ", p, collapse = "\n"), "\n")
  }
}
