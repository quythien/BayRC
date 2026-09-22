## Genome-wide concordance and its permutation p-value for the three baboon
## pairs the case studies quote: SCN against hippocampus, and putamen against
## substantia nigra and against visual cortex.
##
## 500 permutations of gene labels, so the smallest p-value the test can return
## is 1/501 = 0.002. Writes one row per pair to genome_concordance_pvalues.csv
## under BAYRC_OUTPUT_DIR.

suppressPackageStartupMessages(library(BayRC))
this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
source(file.path(analysis.dir, "config.R"))

hb <- new.env()
load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"), envir = hb)

pairs <- list(c("SCN", "HIP"), c("PUT", "SUN"), c("PUT", "VIC"))
n_perm <- 500

rows <- lapply(pairs, function(p) {
  a <- hb$mcmc_data_baboon[[p[1]]]; b <- hb$mcmc_data_baboon[[p[2]]]
  g <- intersect(rownames(a), rownames(b))
  ## subsetting drops the symbols attribute the concordance code reads genes from
  as_mcmc <- function(m) { m <- m[g, ]; attr(m, "symbols") <- g; list(rho = m) }
  out <- file.path(BAYRC_OUTPUT_DIR, "genome_concordance", paste(p, collapse = "_"))
  res <- multi_conservation(
    mcmc.merge.list = list(as_mcmc(a), as_mcmc(b)),
    dataset.names   = paste0("baboon_", p),
    select.pathway.list = "global",
    n_perm = n_perm, n_boot = 100,
    output.dir = out, use_cpp = TRUE,
    compute_pvalue = TRUE, compute_ci = TRUE)
  r <- as.data.frame(res)
  pick <- function(pat) r[[grep(pat, names(r), value = TRUE)[1]]][1]
  data.frame(pair = paste(p, collapse = "-"), genes = length(g),
             adjusted_concordance = pick("AdjustedConcordance$"),
             p_value = pick("PValue$"), permutations = n_perm)
})

tab <- do.call(rbind, rows)
print(tab, row.names = FALSE)
write.csv(tab, file.path(BAYRC_OUTPUT_DIR, "genome_concordance_pvalues.csv"),
          row.names = FALSE)
