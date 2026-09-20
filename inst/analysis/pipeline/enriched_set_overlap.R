## Gene overlap among an arbitrary set of pathways, on measured genes only.
##
## Usage: Rscript enriched_set_overlap.R "Pathway 1" "Pathway 2" ...

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
rho   <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
genes <- rownames(rho[[1]])
kegg  <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))

hits <- commandArgs(trailingOnly = TRUE)
sets <- lapply(kegg[hits], function(p) intersect(p, genes))
names(sets) <- substr(sub("^KEGG ", "", sub(" - multiple diseases", "", hits)), 1, 22)

n <- length(sets)
jac <- matrix(0, n, n, dimnames = list(names(sets), names(sets)))
for (i in seq_len(n)) for (j in seq_len(n))
  jac[i, j] <- length(intersect(sets[[i]], sets[[j]])) /
               length(union(sets[[i]], sets[[j]]))

cat("measured genes per pathway:\n"); print(sapply(sets, length))
cat("\nJaccard:\n"); print(round(jac, 2))
cat("\nsum of set sizes:", sum(sapply(sets, length)),
    "  distinct genes:", length(unique(unlist(sets))), "\n")
cat("redundancy:", round(1 - length(unique(unlist(sets))) /
    sum(sapply(sets, length)), 2), "\n")
