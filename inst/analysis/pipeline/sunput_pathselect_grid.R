## Run the union, gain, loss and conservation enrichment once per KEGG release
## with the size filter opened up, so the stage-1 and stage-2 cut-offs can be
## swept afterwards from the stored tables.
##
## Usage: Rscript sunput_pathselect_grid.R [nperm]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
source(file.path(BAYRC_THIEN_DIR, "pathwaySelect.R"))

suppressPackageStartupMessages({library(dplyr); library(fgsea)})

args <- commandArgs(trailingOnly = TRUE)
nperm <- if (length(args) >= 1) as.integer(args[1]) else 10000

out.dir <- file.path(BAYRC_OUTPUT_DIR, "param_search")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
cache <- readRDS(file.path(out.dir, "cache_PUT_SUN.rds"))
datA <- cache$A; datB <- cache$B

kegg.env <- new.env()
load(file.path(BAYRC_PATHWAY_DIR, "kegg.pathway.list_hsa.RData"), envir = kegg.env)
lists <- list(kegg229 = kegg.env$kegg.pathway.list_hsa,
              kegg354 = readRDS(file.path(BAYRC_PATHWAY_DIR,
                                          "kegg_pathway_list_hsa.rds")))

res <- list()
for (ln in names(lists)) {
  for (rk in c("union", "gain", "loss", "conserved")) {
    cat("\n####", ln, rk, "####\n")
    r <- pathSelect(mcmc.merge.list = list(PUT = datA, SUN = datB),
                    pathway.list = lists[[ln]],
                    dataset.names = c("PUT", "SUN"),
                    ranking.method = rk, score_type = "pos",
                    qvalue.cut = 0.20,
                    pathwaysize.lower.cut = 1,
                    pathwaysize.upper.cut = 5000,
                    nperm = nperm, nproc = 1)
    tab <- r$results
    tab$list <- ln; tab$ranking <- rk
    tab$raw_size <- sapply(lists[[ln]][tab$pathway], length)
    tab$measured <- sapply(lists[[ln]][tab$pathway],
                           function(p) sum(p %in% rownames(datA$rho)))
    res[[paste(ln, rk, sep = "_")]] <- tab
  }
}

saveRDS(res, file.path(out.dir, sprintf("pathselect_raw_nperm%d.rds", nperm)))
all <- do.call(rbind, lapply(res, function(x)
  x[, c("list", "ranking", "pathway", "size", "raw_size", "measured",
        "pval", "padj", "NES")]))
write.csv(all, file.path(out.dir, sprintf("pathselect_all_nperm%d.csv", nperm)),
          row.names = FALSE)
cat("\nwrote", out.dir, "\n")
