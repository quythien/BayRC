# Baboon suprachiasmatic nucleus versus hippocampus: rhythmic transitions,
# phase inference, pathway enrichment and the Figure 3A panel.

library(BayRC)
library(dplyr)

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
source(file.path(analysis.dir, "plots", "theme_bayrc.R"))
source(file.path(analysis.dir, "plots", "peak_concordance.R"))
source(file.path(analysis.dir, "pipeline", "run_record.R"))

# frozen analysis parameters
bfdr_alpha     <- 0.25
shift          <- 2
stage1_q       <- 0.10
stage2_q       <- 0.20
nperm          <- 10000
min_measured   <- 15
panel_pathways <- "KEGG Circadian rhythm"

fig.dir <- file.path(BAYRC_FIGURE_DIR, "baboon_SCN_HIP")
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))

# direction is SCN to HIP
scn <- list(rho = mcmc_data_baboon$SCN, phi = mcmc_phi_baboon$SCN)
hip <- list(rho = mcmc_data_baboon$HIP, phi = mcmc_phi_baboon$HIP)
measured <- rownames(scn$rho)

pA <- rowMeans(scn$rho)
pB <- rowMeans(hip$rho)
trans <- transition_classify(pA, pB, bfdr_alpha = bfdr_alpha)
status <- trans$gain_loss_status

phase <- phase_infer(phi_matrix1 = scn$phi, phi_matrix2 = hip$phi,
                     gain_loss_status = status, bfdr_alpha = bfdr_alpha,
                     shift = shift, P = 24, compute_hdi = TRUE)

maintained <- names(status)[status == "Maintained"]
phase_class <- rep("Undetermined", length(status))
names(phase_class) <- names(status)
phase_class[phase$flag_cons]  <- "Phase-conserved"
phase_class[phase$flag_shift] <- "Phase-shifted"
delta <- ((phase$peak2 - phase$peak1 + 12) %% 24) - 12

clock_genes <- c("BHLHE40", "BHLHE41", "BMAL1", "BTRC", "CLOCK", "CREB1",
                 "CRY1", "CRY2", "CSNK1D", "CSNK1E", "CUL1", "DBP", "FBXL3",
                 "FBXW11", "NFIL3", "NPAS2", "NR1D1", "NR1D2", "PER1", "PER2",
                 "PER3", "RORA", "RORB", "RORC")

# Figure 3A
p <- peak_concordance_plot(
  peak_x = phase$peak1[maintained], peak_y = phase$peak2[maintained],
  phase_class = phase_class, label_genes = clock_genes,
  title = "Circadian Peak Concordance: Baboon SCN versus Hippocampus",
  xlab = "Peak Hour - Suprachiasmatic nucleus (ZT)",
  ylab = "Peak Hour - Hippocampus (ZT)", window = shift)
bayrc_save(p, file.path(fig.dir, "Baboon_SCN_HIP_Peak_Concordance"),
           width = 9, height = 8)

# gene sets cut to the measured genes before the enrichment sees them
kegg <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))
kegg <- lapply(kegg, function(g) replace(g, g == "ARNTL", "BMAL1"))
kegg <- lapply(kegg, intersect, measured)
kegg <- kegg[lengths(kegg) >= min_measured]

select_pathways <- function(plist, method)
  pathSelect(mcmc.merge.list = list(SCN = scn, HIP = hip), pathway.list = plist,
             dataset.names = c("SCN", "HIP"), ranking.method = method,
             score_type = "pos", qvalue.cut = stage2_q,
             pathwaysize.lower.cut = min_measured,
             pathwaysize.upper.cut = length(measured),
             nperm = nperm, nproc = 1)$results

union_res <- select_pathways(kegg, "union")
union_res$q <- p.adjust(union_res$pval, "BH")
active <- union_res$pathway[union_res$q < stage1_q]

# pathSelect names its effect column after the ranking method, so the three
# runs are cut to the shared columns before they are stacked
stage2_cols <- c("pathway", "size", "pval", "Expected_N_Gain",
                 "Expected_N_Loss", "Expected_N_Conserved")
stage2 <- do.call(rbind, lapply(c("gain", "loss", "conserved"), function(m) {
  if (!length(active)) return(NULL)
  r <- select_pathways(kegg[active], m)[, stage2_cols]
  r$q <- p.adjust(r$pval, "BH")
  r$direction <- m
  r
}))
if (is.null(stage2))
  stage2 <- data.frame(pathway = character(), size = integer(),
                       pval = numeric(), Expected_N_Gain = numeric(),
                       Expected_N_Loss = numeric(),
                       Expected_N_Conserved = numeric(), q = numeric(),
                       direction = character())
sig <- stage2[stage2$q < stage2_q, ]

write.csv(union_res[order(union_res$pval), c("pathway", "size", "pval", "q")],
          file.path(fig.dir, "stage1_union.csv"), row.names = FALSE)
write.csv(sig[order(sig$direction, sig$pval),
              c("pathway", "direction", "size", "pval", "q")],
          file.path(fig.dir, "stage2_significant.csv"), row.names = FALSE)

for (pw in panel_pathways) {
  if (!pw %in% names(kegg)) stop("pathway not in the gene set list: ", pw)
  plot_heatmap(data1 = scn, data2 = hip, pathway_genes = kegg[[pw]],
               pathway_name = pw, phase_results = phase,
               transition_results = trans,
               group_names = c("Baboon SCN", "Baboon HIP"),
               versions = "both", save_path = fig.dir)
}

write_run_record(file.path(fig.dir, "run_record.txt"), "applications/Baboon_SCN_HIP.R",
                 list(bfdr_alpha = bfdr_alpha, shift = shift,
                      stage1_q = stage1_q, stage2_q = stage2_q, nperm = nperm,
                      min_measured = min_measured,
                      pathway_list = "kegg_pathway_list_hsa.rds",
                      panels = panel_pathways),
                 repo = analysis.dir)

shifted <- names(phase$flag_shift)[phase$flag_shift]
conserved <- names(phase$flag_cons)[phase$flag_cons]
within <- circ_diff(phase$peak1[maintained], phase$peak2[maintained]) <= shift

cat("\n== SCN vs HIP, the numbers the paper quotes ==\n")
cat("rhythmic in SCN:", sum(pA >= trans$tau_rhythmic_A),
    " rhythmic in HIP:", sum(pB >= trans$tau_rhythmic_B), "\n")
cat("gain:", sum(status == "Gain"), " loss:", sum(status == "Loss"),
    " maintained:", length(maintained), "\n")
cat(sprintf("within +/-%g h: %d of %d (%.1f%%)\n", shift, sum(within),
            length(maintained), 100 * mean(within)))
cat("phase-shifted:", length(shifted), " phase-conserved:", length(conserved),
    " undetermined:", length(maintained) - length(shifted) - length(conserved), "\n")
if (length(conserved))
  cat(sprintf("mean |delta| among phase-conserved: %.2f h\n",
              mean(abs(delta[conserved]))))
cat("pathways tested:", length(kegg), " stage 1 active:", length(active),
    " stage 2 significant:", length(unique(sig$pathway)), "\n\n")
print(union_res[union_res$q < stage1_q, c("pathway", "size", "pval", "q")],
      row.names = FALSE)
if (nrow(sig))
  print(sig[order(sig$direction, sig$pval), c("pathway", "direction", "pval", "q")],
        row.names = FALSE)
cat("\nfigures:", fig.dir, "\n")
