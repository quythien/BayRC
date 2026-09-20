## Full two-stage enrichment for one tissue pair, keeping every row of both
## stages rather than only the significant ones.
##
## Usage: Rscript pair_two_stage_full.R tissueA tissueB [nperm]

suppressPackageStartupMessages({library(BayRC); library(dplyr)})

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

STAGE1_Q <- 0.20; STAGE2_Q <- 0.05; MIN_MEASURED <- 15
args <- commandArgs(trailingOnly = TRUE)
tA <- args[1]; tB <- args[2]
nperm <- if (length(args) >= 3) as.integer(args[3]) else 10000

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
rho <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
measured <- rownames(rho[[tA]])
datA <- list(rho = rho[[tA]]); datB <- list(rho = rho[[tB]])
attr(datA$rho, "symbols") <- measured; attr(datB$rho, "symbols") <- measured

kegg <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))
kegg <- lapply(kegg, intersect, measured)
kegg <- kegg[sapply(kegg, length) >= MIN_MEASURED]

run <- function(sets, rk) {
  r <- pathSelect(mcmc.merge.list = setNames(list(datA, datB), c(tA, tB)),
                  pathway.list = sets, dataset.names = c(tA, tB),
                  ranking.method = rk, score_type = "pos", qvalue.cut = 0.20,
                  pathwaysize.lower.cut = MIN_MEASURED,
                  pathwaysize.upper.cut = length(measured),
                  nperm = nperm, nproc = 1)
  t <- r$results; t$q <- p.adjust(t$pval, "BH"); t$direction <- rk
  t[, c("pathway", "size", "direction", "NES", "pval", "q")]
}

s1 <- run(kegg, "union")
s1 <- s1[order(s1$q), ]
active <- s1$pathway[s1$q < STAGE1_Q]

s2 <- do.call(rbind, lapply(c("gain", "loss", "conserved"),
                            function(rk) run(kegg[active], rk)))
s2 <- s2[order(s2$q), ]

out <- file.path(BAYRC_OUTPUT_DIR, sprintf("two_stage_%s_%s", tA, tB))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write.csv(s1, file.path(out, "stage1_all.csv"), row.names = FALSE)
write.csv(s2, file.path(out, "stage2_all.csv"), row.names = FALSE)

cat(sprintf("\n===== %s vs %s =====\n", tA, tB))
cat(sprintf("\nSTAGE 1: %d tested, %d active at q < %.2f\n",
            nrow(s1), length(active), STAGE1_Q))
print(head(s1[, c("pathway", "size", "NES", "pval", "q")], 15), row.names = FALSE)
cat(sprintf("\nSTAGE 2: all %d rows over the %d active pathways\n",
            nrow(s2), length(active)))
print(s2, row.names = FALSE)
cat(sprintf("\nsignificant at q < %.2f: %d rows\n", STAGE2_Q, sum(s2$q < STAGE2_Q)))
cat("wrote", out, "\n")
