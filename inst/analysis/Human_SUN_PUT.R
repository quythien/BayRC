# Human Substantia Nigra vs Putamen


rm(list=ls())

current_gtex <- "/home/qtp1/Projects/Collaborative"
current_wd <- "/home/qtp1/Projects/Circadian"
current_aging <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging"

outdir = "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/all_plots"

load(file = "/home/qtp1/Projects/Collaborative/GTEXdata/result/summary/hb/mcmc_rho_BF3.RData")
load(file.path(current_gtex, "GTEXdata/result/summary/hb/phi/mcmc_phi_BF3.RData"))

# Objects from these:
# "mcmc_data_baboon" "mcmc_data_human"  "mcmc_phi_baboon"  "mcmc_phi_human"

library(KEGGREST)
library(parallel)

options(width = 10000)

#── Packages ─────────────────────────────────────────────────────────────────────
library(parallel); library(edgeR);   library(Rcpp)
require(DESeq2)
require(limma)
require('minpack.lm')
require(doParallel)
require(biomaRt)
require(devtools)
require(AWFisher)
require(ggplot2)
require(Rcpp)
require(dplyr)
require(pROC)
require(edgeR)

thien_dir <- file.path(current_wd, "Kyle/Circadian-analysis-main/R/v1/BayRC/Thien")
if (file.exists(file.path(thien_dir, "pathwaySelect.R"))) {
  source(file.path(thien_dir, "pathwaySelect.R"))
}
if (file.exists(file.path(thien_dir, "multi_pathway.R"))) {
  source(file.path(thien_dir, "multi_pathway.R"))
}
if (file.exists(file.path(thien_dir, "multi_global.R"))) {
  source(file.path(thien_dir, "multi_global.R"))
}
if (file.exists(file.path(thien_dir, "Permutation_Sim.R"))) {
  source(file.path(thien_dir, "Permutation_Sim.R"))
}
if (file.exists(file.path(thien_dir, "internal.R"))) {
  source(file.path(thien_dir, "internal.R"))
}
if (file.exists(file.path(thien_dir, "congruence.cpp"))) {
  Rcpp::sourceCpp(file.path(thien_dir, "congruence.cpp"))
}
if (file.exists(file.path(thien_dir, "permutation_functions.cpp"))) {
  Rcpp::sourceCpp(file.path(thien_dir, "permutation_functions.cpp"))
}
sourceCpp(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Wei/ACS.cpp"))

load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/hw_orth.RData"))
load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/human.pathway.list.RData"))
load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/go.pathway.list_hsa.RData"))

kegg.pathway.list_hsa <- readRDS("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/pathway_data/kegg_pathway_list_hsa.rds")

#── Source R scripts ─────────────────────────────────────────────────────────────
WD <- "Kyle/Circadian-analysis-main/R/v1"
setwd(file.path(current_wd, WD))
scripts <- list.files("BayRC", pattern="\\.R$", full.names=TRUE)
sapply(scripts, source)

# Reset wd
setwd(current_aging)

#---------------------------------------------------------------------------------
to_zt <- function(t_cos) ifelse(t_cos >= 18, t_cos - 24, t_cos)

#---------------------------------------------------------------------------------
# BRAIN OBJECTS (Human Substantia Nigra vs Putamen)
#---------------------------------------------------------------------------------

human_PUT <- list(
  rho = mcmc_data_human$PUT,
  phi = mcmc_phi_human$PUT
)

human_SUN <- list(
  rho = mcmc_data_human$SUN,
  phi = mcmc_phi_human$SUN
)

#---------------------------------------------------------------------------------
# Global concordance score
output.dir <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/results/human/output_final"
if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
}
cat("Directory ready:", output.dir, "\n")

# Analyze ALL genes together
t_start <- Sys.time()

results_global1 <- multi_conservation(
  mcmc.merge.list = list(human_PUT, human_SUN),
  dataset.names = c("human_PUT", "human_SUN"),
  select.pathway.list = "global",
  n_perm = 500,
  n_boot = 100,
  output.dir = file.path(output.dir, "PUT_SUN_concordance"),
  use_cpp = TRUE,
  compute_pvalue = TRUE,
  compute_ci = TRUE
)

t_end <- Sys.time()
elapsed_time <- difftime(t_end, t_start, units = "mins")
print(elapsed_time)

#---------------------------------------------------------------------------------
# Outer (rhythmic) and inner (phase) analysis
#---------------------------------------------------------------------------------

pA    <- rowMeans(human_PUT$rho)
pB    <- rowMeans(human_SUN$rho)

# BFDR: Rhythmic biomarkers
rhy   <- detect_rhy(human_PUT, human_SUN, 0.25)
rhy_summary <- data.frame(
  Category = c("Rhythmic in PUT",
               "Rhythmic in SUN"),
  Count = c(rhy$n_rhythmic_A,
            rhy$n_rhythmic_B)
)
rhy_summary

# BFDR:
trans_outer <- transition_classify(pA, pB, bfdr_alpha = 0.25)

phase_inner <- phase_infer(
  phi_matrix1 = human_PUT$phi,
  phi_matrix2 = human_SUN$phi,
  gain_loss_status = trans_outer$gain_loss_status,
  bfdr_alpha = 0.25,
  shift = 2,
  P = 24, compute_hdi = TRUE
)

# === PHASE INFERENCE SUMMARY ===
cat("\n=== PHASE INFERENCE SUMMARY ===\n")
cat("Maintained genes:", length(which(trans_outer$gain_loss_status == "Maintained")), "\n")
cat("BFDR α = 0.25  | shift threshold = 2 h\n")
cat("Significant phase-shifted genes:", sum(phase_inner$flag_shift, na.rm = TRUE), "\n")
cat("Significant phase-conserved genes:", sum(phase_inner$flag_cons, na.rm = TRUE), "\n")

#---------------------------------------------------------------------------------
# CHECK PEAK TIMES FOR PARKINSON DISEASE GENES
#---------------------------------------------------------------------------------
cat("\n=== PARKINSON DISEASE PATHWAY ANALYSIS ===\n")

pd_genes <- kegg.pathway.list_hsa[["KEGG Parkinson disease"]]
cat("Total Parkinson pathway genes:", length(pd_genes), "\n")

# Get peak times for Parkinson genes
pd_peak1 <- phase_inner$peak1[pd_genes]
pd_peak2 <- phase_inner$peak2[pd_genes]

# Create dataframe
pd_df <- data.frame(
  Gene = pd_genes,
  PUT = pd_peak1,
  SUN = pd_peak2
)

# Complete cases
pd_complete <- pd_df[complete.cases(pd_df), ]
cat("Genes with complete data:", nrow(pd_complete), "\n")

# Calculate difference
pd_complete$Diff <- pd_complete$SUN - pd_complete$PUT

cat("\n=== SUMMARY STATISTICS (SUN - PUT) ===\n")
print(summary(pd_complete$Diff))

cat("\nMean PUT peak:", mean(pd_complete$PUT), "\n")
cat("Mean SUN peak:", mean(pd_complete$SUN), "\n")
cat("Mean difference (SUN - PUT):", mean(pd_complete$Diff), "\n")

# Count genes where SUN is earlier vs later
n_sun_earlier <- sum(pd_complete$Diff < 0, na.rm = TRUE)
n_sun_later <- sum(pd_complete$Diff > 0, na.rm = TRUE)

cat("\nSUN earlier than PUT:", n_sun_earlier, "genes\n")
cat("SUN later than PUT:", n_sun_later, "genes\n")

# Show first 20 genes
cat("\n=== FIRST 20 PARKINSON GENES ===\n")
print(head(pd_complete, 20))

# Check NDUFB5 specifically
cat("\n=== NDUFB5 SPECIFICALLY ===\n")
if ("NDUFB5" %in% pd_complete$Gene) {
  ndufb5_data <- pd_complete[pd_complete$Gene == "NDUFB5", ]
  print(ndufb5_data)
  if (ndufb5_data$SUN > ndufb5_data$PUT) {
    cat("=> SUN peaks", ndufb5_data$SUN - ndufb5_data$PUT, "hours LATER than PUT\n")
  } else {
    cat("=> SUN peaks", ndufb5_data$PUT - ndufb5_data$SUN, "hours EARLIER than PUT\n")
  }
} else {
  cat("NDUFB5 not found in complete data\n")
}

#---------------------------------------------------------------------------------
# CHECK ALL RHYTHMIC GENES
#---------------------------------------------------------------------------------
cat("\n=== ALL RHYTHMIC GENES ANALYSIS ===\n")

all_genes <- names(phase_inner$peak1)
all_df <- data.frame(
  Gene = all_genes,
  PUT = phase_inner$peak1,
  SUN = phase_inner$peak2
)

all_complete <- all_df[complete.cases(all_df), ]
cat("Total genes with complete data:", nrow(all_complete), "\n")

all_complete$Diff <- all_complete$SUN - all_complete$PUT

cat("\n=== SUMMARY STATISTICS - ALL GENES (SUN - PUT) ===\n")
print(summary(all_complete$Diff))

cat("\nMean PUT peak (all genes):", mean(all_complete$PUT), "\n")
cat("Mean SUN peak (all genes):", mean(all_complete$SUN), "\n")
cat("Mean difference (SUN - PUT):", mean(all_complete$Diff), "\n")

n_sun_earlier_all <- sum(all_complete$Diff < 0, na.rm = TRUE)
n_sun_later_all <- sum(all_complete$Diff > 0, na.rm = TRUE)

cat("\nSUN earlier than PUT:", n_sun_earlier_all, "genes\n")
cat("SUN later than PUT:", n_sun_later_all, "genes\n")
