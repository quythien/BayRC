## Over-representation of the phase-shifted gene set across the grid cells the
## phase sweep stored, with and without restricting the background to the genes
## the atlas measures.
##
## Usage: Rscript sunput_shifted_ora_grid.R

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(dplyr)
})

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

out.dir <- file.path(BAYRC_OUTPUT_DIR, "param_search")
cells <- readRDS(file.path(out.dir, "phase_grid_genes_PUT_SUN.rds"))
cache <- readRDS(file.path(out.dir, "cache_PUT_SUN.rds"))
measured <- rownames(cache$A$rho)

to_entrez <- function(g) {
  m <- suppressWarnings(suppressMessages(
    bitr(g, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)))
  unique(m$ENTREZID)
}
bg_all <- to_entrez(measured)

key <- "Oxidative phosphorylation|Parkinson|Alzheimer|Huntington|neurodegener|Proteasome|respirat|mitochond"

rows <- list(); tops <- list()
for (cn in names(cells)) {
  g <- to_entrez(cells[[cn]]$gene)
  if (!length(g)) next
  for (bg in c("annotated", "measured")) {
    universe <- if (bg == "measured") bg_all else NULL
    for (db in c("GO_BP", "KEGG")) {
      e <- try({
        if (db == "GO_BP")
          enrichGO(g, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                   pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                   universe = universe, readable = TRUE)
        else
          enrichKEGG(g, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1,
                     universe = universe)
      }, silent = TRUE)
      if (inherits(e, "try-error") || is.null(e)) next
      df <- as.data.frame(e)
      if (!nrow(df)) next
      rows[[length(rows) + 1]] <- data.frame(
        cell = cn, background = bg, db = db, genes = length(g),
        n_q05 = sum(df$p.adjust < 0.05), n_q20 = sum(df$p.adjust < 0.20),
        min_q = signif(min(df$p.adjust), 3),
        key_q05 = sum(df$p.adjust < 0.05 & grepl(key, df$Description,
                                                 ignore.case = TRUE)))
      tops[[paste(cn, bg, db)]] <- head(
        df[order(df$p.adjust), c("Description", "GeneRatio", "pvalue", "p.adjust")], 12)
    }
  }
}

d <- do.call(rbind, rows)
write.csv(d, file.path(out.dir, "shifted_ora_grid.csv"), row.names = FALSE)
saveRDS(tops, file.path(out.dir, "shifted_ora_tops.rds"))
print(d, row.names = FALSE)
for (nm in names(tops)) {
  cat("\n====", nm, "====\n"); print(tops[[nm]], row.names = FALSE)
}
