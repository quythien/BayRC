## Sweep bfdr_alpha and the phase window for SUN vs PUT and record how the
## mitochondrial and proteasome gene sets that carry the pathway heatmaps
## behave in each cell.
##
## One alpha per invocation keeps the sweep parallel; the cells are collected
## afterwards from the per-alpha files.
##
## Usage: Rscript sunput_phase_grid.R [alpha] [cache.rds]

suppressPackageStartupMessages(library(BayRC))

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

out.dir <- file.path(BAYRC_OUTPUT_DIR, "param_search")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
alpha.arg <- if (length(args) >= 1) as.numeric(args[1]) else NA
cache.file <- if (length(args) >= 2) args[2] else
  file.path(out.dir, "cache_PUT_SUN.rds")
cache <- readRDS(cache.file)
datA <- cache$A; datB <- cache$B
nmA <- cache$names[1]; nmB <- cache$names[2]

## marker sets, taken from both KEGG releases so the grid is not tied to one
kegg.env <- new.env()
load(file.path(BAYRC_PATHWAY_DIR, "kegg.pathway.list_hsa.RData"), envir = kegg.env)
kegg229 <- kegg.env$kegg.pathway.list_hsa
kegg354 <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))
genes <- rownames(datA$rho)

pick <- function(lst, nm) if (nm %in% names(lst)) intersect(lst[[nm]], genes) else character(0)
marker <- list(
  OxPhos    = union(pick(kegg229, "KEGG Oxidative phosphorylation"),
                    pick(kegg354, "KEGG Oxidative phosphorylation")),
  Proteasome = union(pick(kegg229, "KEGG Proteasome"),
                     pick(kegg354, "KEGG Proteasome")),
  Parkinson = union(pick(kegg229, "KEGG Parkinson's disease"),
                    pick(kegg354, "KEGG Parkinson disease")))

pA <- rowMeans(datA$rho); pB <- rowMeans(datB$rho)

## signed median phase difference, B minus A, on the circle
dphi <- ((datB$phi - datA$phi + 12) %% 24) - 12
med  <- apply(dphi, 1, median)
rm(dphi)

shifts <- c(1, 1.5, 2, 2.5, 3, 4)
alphas <- c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40)
if (!is.na(alpha.arg)) alphas <- alpha.arg
tag <- if (is.na(alpha.arg)) "" else sprintf("_a%g", alpha.arg)

rows <- list(); gene.tab <- list()
for (a in alphas) {
  tr <- transition_classify(pA, pB, bfdr_alpha = a)
  st <- tr$gain_loss_status
  n_maint <- sum(st == "Maintained"); n_gain <- sum(st == "Gain")
  n_loss <- sum(st == "Loss")
  for (s in shifts) {
    ph <- phase_infer(phi_matrix1 = datA$phi, phi_matrix2 = datB$phi,
                      gain_loss_status = st, bfdr_alpha = a, shift = s,
                      P = 24, compute_hdi = FALSE)
    shifted <- names(ph$flag_shift)[ph$flag_shift]
    cons    <- names(ph$flag_cons)[ph$flag_cons]
    row <- data.frame(alpha = a, shift = s, gain = n_gain, loss = n_loss,
                      maintained = n_maint, shifted = length(shifted),
                      conserved = length(cons),
                      undetermined = n_maint - length(shifted) - length(cons))
    for (mk in names(marker)) {
      g <- marker[[mk]]
      m_maint <- intersect(g, names(st)[st == "Maintained"])
      m_sh <- intersect(g, shifted); m_co <- intersect(g, cons)
      row[[paste0(mk, "_maint")]] <- length(m_maint)
      row[[paste0(mk, "_shift")]] <- length(m_sh)
      row[[paste0(mk, "_cons")]]  <- length(m_co)
      row[[paste0(mk, "_mean_d")]] <- if (length(m_sh)) round(mean(med[m_sh]), 2) else NA
      row[[paste0(mk, "_same_sign")]] <- if (length(m_sh))
        round(max(mean(med[m_sh] > 0), mean(med[m_sh] < 0)), 3) else NA
    }
    rows[[length(rows) + 1]] <- row
    if (s %in% c(1.5, 2, 2.5))
      gene.tab[[sprintf("a%g_s%g", a, s)]] <- data.frame(
        gene = shifted, delta = round(med[shifted], 2))
  }
}

d <- do.call(rbind, rows)
d$pct_undet <- round(100 * d$undetermined / d$maintained, 1)
write.csv(d, file.path(out.dir, sprintf("phase_grid_%s_%s%s.csv", nmA, nmB, tag)),
          row.names = FALSE)
saveRDS(gene.tab, file.path(out.dir,
        sprintf("phase_grid_genes_%s_%s%s.rds", nmA, nmB, tag)))

cat("\n", nmA, "vs", nmB, "\n\n")
print(d[, c("alpha", "shift", "maintained", "shifted", "conserved",
            "undetermined", "pct_undet")], row.names = FALSE)
cat("\nmarker sets\n")
print(d[, c("alpha", "shift", "OxPhos_maint", "OxPhos_shift", "OxPhos_cons",
            "OxPhos_mean_d", "OxPhos_same_sign", "Proteasome_shift",
            "Parkinson_maint", "Parkinson_shift", "Parkinson_cons",
            "Parkinson_mean_d", "Parkinson_same_sign")], row.names = FALSE)
