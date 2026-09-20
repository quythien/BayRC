## Single-stage enrichment for one tissue pair: test every pathway for gain,
## loss and conservation directly, with no Stage 1 screen.
##
## Usage: Rscript pair_single_stage_enrichment.R tissueA tissueB [nperm]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
suppressPackageStartupMessages({library(BayRC); library(dplyr); library(fgsea)})

args  <- commandArgs(trailingOnly = TRUE)
tA <- args[1]; tB <- args[2]
nperm <- if (length(args) >= 3) as.integer(args[3]) else 10000
min_measured <- 15

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
rho <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
genes <- rownames(rho[[tA]])

datA <- list(rho = rho[[tA]]); datB <- list(rho = rho[[tB]])
attr(datA$rho, "symbols") <- genes; attr(datB$rho, "symbols") <- genes

kegg <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))
kegg <- lapply(kegg, function(p) intersect(p, genes))
kegg <- kegg[sapply(kegg, length) >= min_measured]
cat("pathways tested:", length(kegg), "\n")

out <- list()
for (rk in c("gain", "loss", "conserved")) {
  r <- pathSelect(mcmc.merge.list = setNames(list(datA, datB), c(tA, tB)),
                  pathway.list = kegg, dataset.names = c(tA, tB),
                  ranking.method = rk, score_type = "pos", qvalue.cut = 0.20,
                  pathwaysize.lower.cut = min_measured,
                  pathwaysize.upper.cut = length(genes),
                  nperm = nperm, nproc = 1)
  t <- r$results
  t$q <- p.adjust(t$pval, "BH")
  t$direction <- rk
  out[[rk]] <- t[, c("pathway", "size", "direction", "pval", "q")]
}
d <- do.call(rbind, out)
d <- d[order(d$q), ]
write.csv(d, file.path(BAYRC_OUTPUT_DIR,
          sprintf("single_stage_%s_%s.csv", tA, tB)), row.names = FALSE)

for (cut in c(0.05, 0.20)) {
  s <- d[d$q < cut, ]
  cat(sprintf("\n=== single stage, q < %.2f: %d rows over %d pathways ===\n",
              cut, nrow(s), length(unique(s$pathway))))
  if (nrow(s)) print(head(s, 20), row.names = FALSE)
}
