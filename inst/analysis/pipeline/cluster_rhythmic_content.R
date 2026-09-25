## What separates the two tissue clusters in the Figure 6B concordance heatmap.
##
## The concordance score compares which genes are called rhythmic in a pair, so
## a cluster can be tight either because its tissues agree on timing or because
## they simply carry more rhythmic content. This reports both, for the split the
## panel's own dendrogram gives.
##
## Usage: Rscript cluster_rhythmic_content.R

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
suppressPackageStartupMessages(library(BayRC))

bfdr_alpha <- 0.25
pathway    <- "KEGG Circadian rhythm"

conc <- file.path(BAYRC_OUTPUT_DIR, "heatmap_circadian_pairs", "within_baboon",
                  "pairwise_concordance_baboon_circadian_Matrix.csv")
m <- as.matrix(read.csv(conc, row.names = 1, check.names = FALSE))
diag(m) <- 1
groups <- cutree(hclust(as.dist(1 - m), method = "ward.D2"), k = 2)
within <- vapply(sort(unique(groups)), function(g) {
  s <- m[groups == g, groups == g, drop = FALSE]
  mean(s[row(s) != col(s)])
}, numeric(1))
tight <- names(groups)[groups == which.max(within)]

cat("tight cluster, mean within-cluster c-score", sprintf("%.3f\n", max(within)))
cat("  ", paste(sort(tight), collapse = " "), "\n")
cat("remaining, mean within-group c-score", sprintf("%.3f\n", min(within)))
cat("  ", paste(sort(setdiff(names(groups), tight)), collapse = " "), "\n")

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
rho <- get(grep("^mcmc_data", ls(), value = TRUE)[1])

kegg <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))
kegg <- lapply(kegg, function(g) replace(g, g == "ARNTL", "BMAL1"))
path_genes <- kegg[[pathway]]

tissues <- intersect(names(groups), names(rho))
res <- do.call(rbind, lapply(tissues, function(t) {
  p   <- rowMeans(rho[[t]])
  tau <- bfdr_from_posterior(p, alpha = bfdr_alpha)$threshold
  pg  <- intersect(path_genes, names(p))
  data.frame(tissue = t,
             group = if (t %in% tight) "cluster" else "rest",
             pathway_rhythmic = sum(p[pg] >= tau),
             pathway_measured = length(pg),
             genome_percent = round(100 * mean(p >= tau), 1))
}))
res$pathway_percent <- round(100 * res$pathway_rhythmic / res$pathway_measured, 1)
res <- res[order(res$group, -res$pathway_percent), ]

print(res, row.names = FALSE)
write.csv(res, file.path(BAYRC_OUTPUT_DIR, "cluster_rhythmic_content.csv"),
          row.names = FALSE)

cat("\ngroup means\n")
for (v in c("pathway_percent", "genome_percent"))
  cat(sprintf("  %-16s cluster %5.1f%%   rest %5.1f%%\n", v,
              mean(res[[v]][res$group == "cluster"]),
              mean(res[[v]][res$group == "rest"])))

# whether the concordance a tissue shows is explained by how much it carries
mean_c <- vapply(tissues, function(t) mean(m[t, setdiff(tissues, t)]), numeric(1))
cat(sprintf("\ncorrelation of mean c-score with pathway rhythmic percent: %.3f\n",
            cor(mean_c[res$tissue], res$pathway_percent)))
cat(sprintf("correlation of mean c-score with genome-wide percent:      %.3f\n",
            cor(mean_c[res$tissue], res$genome_percent)))
cat("\nwrote", file.path(BAYRC_OUTPUT_DIR, "cluster_rhythmic_content.csv"), "\n")
