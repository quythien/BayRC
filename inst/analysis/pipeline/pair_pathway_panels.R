## Draw integrated pathway heatmaps for an arbitrary tissue pair and pathway
## list, into an exploration directory separate from the case studies.
##
## Usage: Rscript pair_pathway_panels.R tissueA tissueB "Pathway 1" "Pathway 2" ...

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
suppressPackageStartupMessages(library(BayRC))

args <- commandArgs(trailingOnly = TRUE)
tA <- args[1]; tB <- args[2]; pws <- args[-(1:2)]
bfdr_alpha <- 0.25; shift <- 2

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))
rho  <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
phi  <- get(grep("^mcmc_phi",  ls(), value = TRUE)[1])
kegg <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))

datA <- list(rho = rho[[tA]], phi = phi[[tA]])
datB <- list(rho = rho[[tB]], phi = phi[[tB]])
attr(datA$rho, "symbols") <- rownames(rho[[tA]])
attr(datB$rho, "symbols") <- rownames(rho[[tB]])

trans <- transition_classify(rowMeans(rho[[tA]]), rowMeans(rho[[tB]]),
                             bfdr_alpha = bfdr_alpha)
phase <- phase_infer(phi_matrix1 = phi[[tA]], phi_matrix2 = phi[[tB]],
                     gain_loss_status = trans$gain_loss_status,
                     bfdr_alpha = bfdr_alpha, shift = shift, P = 24,
                     compute_hdi = FALSE)

fig.dir <- file.path(BAYRC_FIGURE_DIR, "explore", paste(tA, tB, sep = "_"))
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)

for (pw in pws) {
  if (!pw %in% names(kegg)) { cat("not in list:", pw, "\n"); next }
  plot_heatmap(data1 = datA, data2 = datB, pathway_genes = kegg[[pw]],
               pathway_name = pw, phase_results = phase,
               transition_results = trans,
               group_names = c(paste("Baboon", tA), paste("Baboon", tB)),
               versions = "both", save_path = fig.dir)
}
cat("\nwrote", fig.dir, "\n")
