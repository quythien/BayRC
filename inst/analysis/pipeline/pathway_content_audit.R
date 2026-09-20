## How much of each tested pathway is electron transport chain or proteasome
## content, and where the neurodegeneration maps sit by size.

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
rho   <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
genes <- rownames(rho[[1]])
kegg  <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))

meas <- sapply(kegg, function(p) sum(genes %in% p))
k    <- kegg[meas >= 15]
m    <- meas[meas >= 15]

ox <- intersect(kegg[["KEGG Oxidative phosphorylation"]], genes)
pr <- intersect(kegg[["KEGG Proteasome"]], genes)
nox <- sapply(k, function(p) length(intersect(p, ox)))
npr <- sapply(k, function(p) length(intersect(p, pr)))

pat <- "Parkinson|Alzheimer|Huntington|Prion|lateral sclerosis|neurodegener"
neuro <- grep(pat, names(k), ignore.case = TRUE)

cat("pathways tested:", length(k), "   median measured size:", median(m), "\n\n")
cat("the neurodegeneration maps:\n")
d <- data.frame(pathway = substr(sub("^KEGG ", "", names(k)[neuro]), 1, 36),
                measured = m[neuro], size_rank = rank(-m)[neuro],
                ETC = nox[neuro], proteasome = npr[neuro])
d$pct_ETC_prot <- round(100 * (d$ETC + d$proteasome) / d$measured)
print(d[order(d$size_rank), ], row.names = FALSE)

cat("\npathways carrying >= 20 of the", length(ox), "measured ETC genes:\n")
print(sort(nox[nox >= 20], decreasing = TRUE))
cat("\nof those, how many are neurodegeneration maps:",
    sum(grepl(pat, names(nox)[nox >= 20], ignore.case = TRUE)), "\n")
