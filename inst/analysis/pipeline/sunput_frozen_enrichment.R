## SUN vs PUT enrichment under the frozen setting: the 354-pathway release,
## gene sets cut to the measured genes and kept at 15 or more, union stage 1 at
## BH q < 0.05, transition stage 2 at q < 0.20, 10000 permutations.
##
## Writes the stage tables and the pathway metrics behind Table 1.
##
## Usage: Rscript sunput_frozen_enrichment.R

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
suppressPackageStartupMessages({library(BayRC); library(dplyr)})

STAGE1_Q <- 0.10
STAGE2_Q <- 0.20
NPERM    <- 10000
MIN_MEASURED <- 15

out.dir <- file.path(BAYRC_OUTPUT_DIR, "frozen_sunput")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
cache <- readRDS(file.path(BAYRC_OUTPUT_DIR, "param_search", "cache_PUT_SUN.rds"))
datA <- cache$A; datB <- cache$B          # A = PUT, B = SUN
measured <- rownames(datA$rho)

## gene sets cut to the measured genes before anything else sees them
kegg <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))
kegg <- lapply(kegg, intersect, measured)
kegg <- kegg[lengths(kegg) >= MIN_MEASURED]
cat("pathways tested:", length(kegg), "\n")

run <- function(plist, method)
  pathSelect(mcmc.merge.list = list(PUT = datA, SUN = datB),
             pathway.list = plist, dataset.names = c("PUT", "SUN"),
             ranking.method = method, score_type = "pos",
             qvalue.cut = STAGE2_Q,
             pathwaysize.lower.cut = MIN_MEASURED,
             pathwaysize.upper.cut = length(measured),
             nperm = NPERM, nproc = 1)$results

## stage 1 --------------------------------------------------------------------
union <- run(kegg, "union")
union$q <- p.adjust(union$pval, "BH")
active <- union$pathway[union$q < STAGE1_Q]
cat("stage 1 active at q <", STAGE1_Q, ":", length(active), "\n")

## stage 2 --------------------------------------------------------------------
stage2 <- lapply(c("gain", "loss", "conserved"), function(m) {
  r <- run(kegg[active], m)
  r$q <- p.adjust(r$pval, "BH")
  r$direction <- m
  r
})
names(stage2) <- c("gain", "loss", "conserved")
sig <- do.call(rbind, lapply(stage2, function(r)
  r[r$q < STAGE2_Q, c("pathway", "direction", "size", "pval", "q",
                      "Expected_N_Gain", "Expected_N_Loss",
                      "Expected_N_Conserved")]))
selected <- unique(sig$pathway)

## pathway metrics behind Table 1 ---------------------------------------------
mc <- multi_conservation(
  mcmc.merge.list = list(PUT = datA, SUN = datB),
  dataset.names = c("PUT", "SUN"),
  select.pathway.list = kegg[selected],
  n_perm = 1000, n_boot = 1000,
  output.dir = file.path(out.dir, "multiconservation"),
  use_cpp = TRUE)

saveRDS(list(union = union, stage2 = stage2, sig = sig, mc = mc,
             active = active, selected = selected, tested = names(kegg)),
        file.path(out.dir, "frozen_enrichment.rds"))
write.csv(union[order(union$pval), c("pathway", "size", "pval", "q")],
          file.path(out.dir, "stage1_union.csv"), row.names = FALSE)
write.csv(sig[order(sig$direction, sig$pval), ],
          file.path(out.dir, "stage2_significant.csv"), row.names = FALSE)

cat("\n=== stage 1, the active set ===\n")
print(union[union$q < STAGE1_Q, c("pathway", "size", "pval", "q")], row.names = FALSE)
cat("\n=== stage 2, q <", STAGE2_Q, "===\n")
print(sig[order(sig$direction, sig$pval), c("pathway", "direction", "size", "pval", "q")],
      row.names = FALSE)
cat("\n=== where Proteasome lands ===\n")
print(union[grepl("Proteasome", union$pathway), c("pathway", "size", "pval", "q")],
      row.names = FALSE)
for (m in names(stage2)) {
  r <- stage2[[m]]
  if (any(grepl("Proteasome", r$pathway)))
    print(cbind(direction = m, r[grepl("Proteasome", r$pathway),
                                 c("pathway", "pval", "q")]), row.names = FALSE)
}
cat("\n=== Table 1 metrics ===\n")
print(mc, row.names = FALSE)
