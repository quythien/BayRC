## Count the genes each candidate Figure 4 pathway puts into the integrated
## heatmap, by transition class and by phase call, at a range of bfdr_alpha.
##
## Usage: Rscript sunput_panel_density.R [shift]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
suppressPackageStartupMessages(library(BayRC))

args  <- commandArgs(trailingOnly = TRUE)
shift <- if (length(args) >= 1) as.numeric(args[1]) else 2

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))
rho <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
phi <- get(grep("^mcmc_phi",  ls(), value = TRUE)[1])

kegg <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))
cand <- c("KEGG Oxidative phosphorylation", "KEGG Parkinson disease",
          "KEGG Prion disease", "KEGG Thyroid hormone synthesis",
          "KEGG Pathways of neurodegeneration - multiple diseases",
          "KEGG Protein processing in endoplasmic reticulum",
          "KEGG Proteasome")

genes <- rownames(rho[["SUN"]])
pA <- rowMeans(rho[["PUT"]]); pB <- rowMeans(rho[["SUN"]])

rows <- list()
for (a in c(0.20, 0.25, 0.30)) {
  tr <- transition_classify(pA, pB, bfdr_alpha = a)
  ph <- phase_infer(phi_matrix1 = phi[["PUT"]], phi_matrix2 = phi[["SUN"]],
                    gain_loss_status = tr$gain_loss_status,
                    bfdr_alpha = a, shift = shift, P = 24, compute_hdi = FALSE)
  for (p in cand) {
    i <- which(genes %in% kegg[[p]])
    st <- tr$gain_loss_status[i]
    rows[[length(rows) + 1]] <- data.frame(
      pathway = substr(sub("^KEGG ", "", p), 1, 34), alpha = a,
      measured = length(i),
      gain = sum(st == "Gain"), loss = sum(st == "Loss"),
      maint = sum(st == "Maintained"),
      panel = sum(st %in% c("Gain", "Loss", "Maintained")),
      shifted = sum(ph$flag_shift[i], na.rm = TRUE),
      cons    = sum(ph$flag_cons[i],  na.rm = TRUE))
  }
}
d <- do.call(rbind, rows)
write.csv(d, file.path(BAYRC_OUTPUT_DIR, "param_search",
                       sprintf("panel_density_shift%g.csv", shift)),
          row.names = FALSE)
for (a in unique(d$alpha)) {
  cat(sprintf("\n=== bfdr_alpha %.2f, shift %g h ===\n", a, shift))
  print(d[d$alpha == a, -2], row.names = FALSE)
}
