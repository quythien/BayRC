## Gene overlap among the enriched SUN-PUT pathways, to show how much of the
## signal is shared rather than independent.

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
rho   <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
genes <- rownames(rho[["SUN"]])
kegg  <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))

hits <- c("KEGG Oxidative phosphorylation", "KEGG Thyroid hormone synthesis",
          "KEGG Retrograde endocannabinoid signaling", "KEGG Parkinson disease",
          "KEGG Prion disease", "KEGG Amyotrophic lateral sclerosis",
          "KEGG Pathways of neurodegeneration - multiple diseases")
sets <- lapply(kegg[hits], function(p) intersect(p, genes))
names(sets) <- sub("^KEGG ", "", sub(" - multiple diseases", "", hits))

n <- length(sets)
shared <- jac <- matrix(0, n, n, dimnames = list(names(sets), names(sets)))
for (i in seq_len(n)) for (j in seq_len(n)) {
  shared[i, j] <- length(intersect(sets[[i]], sets[[j]]))
  jac[i, j]    <- shared[i, j] / length(union(sets[[i]], sets[[j]]))
}

cat("measured genes per pathway:\n"); print(sapply(sets, length))
cat("\nshared measured genes:\n"); print(shared)
cat("\nJaccard:\n"); print(round(jac, 2))

ox <- sets[["Oxidative phosphorylation"]]
cat("\nfraction of each pathway that is also in Oxidative phosphorylation:\n")
print(round(sapply(sets, function(s) length(intersect(s, ox)) / length(s)), 2))
cat("\ngenes in the union of all seven:", length(unique(unlist(sets))), "\n")
cat("sum of individual sizes:", sum(sapply(sets, length)), "\n")
