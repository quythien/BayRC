#── Paths ─────────────────────────────────────────────────────────────────────────
current_gtex <- "/home/qtp1/Projects/Collaborative"
current_wd   <- "/home/qtp1/Projects/Circadian"
current_aging <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging"
output = "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/output"
base_result_dir <- file.path(current_aging, "results", "brain_regions")
options(width = Sys.getenv("COLUMNS"))
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

source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/pathwaySelect.R"))
source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/KEGG_module.R"))
source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/KEGG_module_topology_plot.R"))

source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/multi_ACS_ADS_pathway.R"))
source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/multi_ACS_ADS_global.R"))
source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/internal.R"))
sourceCpp(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/ACS.cpp"))

load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/kegg.pathway.list_hsa.RData"))
load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/kegg.pathway.list_cel_GeneNames.RData"))
load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/hw_orth.RData"))
load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/human.pathway.list.RData"))
load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/go.pathway.list_hsa.RData"))

#── Benchmarking with cosinor results ────────────────────────────────────────────────────────────────
source("/home/qtp1/Projects/Pipeline/one_cosinor_OLS_new.R")

#── Source R scripts ─────────────────────────────────────────────────────────────
WD <- "Kyle/Circadian-analysis-main/R/v1"
setwd(file.path(current_wd, WD))
scripts <- list.files("R", full.names=TRUE)
sapply(scripts, source)



#── Real Data ─────────────────────────────────────────────────────────────
source("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/src/Thien/Permutation_Sim.R")
source("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/src/Thien/Permutation_Sim.R")
mcmc_age = readRDS("/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/data/mcmc_young_old.rds")

output.dir <- "/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/src/Thien"

if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
}

result_bootstrap <- multi_conservation_pathway_bootstrap(
  mcmc.merge.list = mcmc_age,
  dataset.names = c("younger", "older"),
  select.pathway.list = kegg.pathway.list_hsa,
  delta = 4,
  units = "hours",
  n_boot = 500,  # Add this parameter (300-500 recommended)
  output.dir = file.path(output.dir, "bootstrap_results")  # Separate folder
)
require(openxlsx)

# File paths
file1 <- "/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/src/Thien/bootstrap_results/Conservation_Bootstrap_younger_vs_older.xlsx"
file2 <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/results/brain_regions/Conservation_Pathway_Test_1e6/Conservation_Results_2_datasets.xlsx"
require(openxlsx)

# File paths
file1 <- "/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/src/Thien/bootstrap_results/Conservation_Bootstrap_younger_vs_older.xlsx"
file2 <- "/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/src/Thien/permutation/Conservation_Results_Permutation_100000_datasets.xlsx"

# Read data
bootstrap_cong <- read.xlsx(file1, sheet = "congruence_index")
perm_cong <- read.xlsx(file2, sheet = "congruence_index")

# Match pathways
common_pathways <- intersect(bootstrap_cong$Pathway, perm_cong$Pathway)

cat("Total pathways in Bootstrap:", nrow(bootstrap_cong), "\n")
cat("Total pathways in Permutation:", nrow(perm_cong), "\n")
cat("Common pathways:", length(common_pathways), "\n\n")

# Extract p-values for common pathways
bootstrap_cong <- bootstrap_cong[bootstrap_cong$Pathway %in% common_pathways, ]
perm_cong <- perm_cong[perm_cong$Pathway %in% common_pathways, ]

# Sort by pathway name for matching
bootstrap_cong <- bootstrap_cong[order(bootstrap_cong$Pathway), ]
perm_cong <- perm_cong[order(perm_cong$Pathway), ]

################################################################################
# EXTRACT P-VALUES
################################################################################

bootstrap_p <- bootstrap_cong$younger_vs_older_PValue
perm_p <- perm_cong$younger_vs_older_PValue

# Remove NAs
valid_idx <- !is.na(bootstrap_p) & !is.na(perm_p)
bootstrap_p <- bootstrap_p[valid_idx]
perm_p <- perm_p[valid_idx]
pathways_valid <- bootstrap_cong$Pathway[valid_idx]

cat("Valid pathways for comparison:", length(pathways_valid), "\n\n")

################################################################################
# CORRELATIONS
################################################################################

cor_spearman <- cor(bootstrap_p, perm_p, method = "spearman")
cor_pearson <- cor(bootstrap_p, perm_p, method = "pearson")

cat("=== CORRELATIONS ===\n")
cat("Spearman correlation:", round(cor_spearman, 4), "\n")
cat("Pearson correlation:", round(cor_pearson, 4), "\n\n")

################################################################################
# AGREEMENT ON SIGNIFICANCE
################################################################################

bootstrap_sig <- bootstrap_p < 0.05
perm_sig <- perm_p < 0.05
agreement <- sum(bootstrap_sig == perm_sig) / length(bootstrap_sig)

cat("Agreement on significance (p < 0.05):", round(agreement * 100, 1), "%\n")
cat("Both significant:", sum(bootstrap_sig & perm_sig), "\n")
cat("Bootstrap only:", sum(bootstrap_sig & !perm_sig), "\n")
cat("Permutation only:", sum(!bootstrap_sig & perm_sig), "\n")
cat("Both non-significant:", sum(!bootstrap_sig & !perm_sig), "\n\n")

# Discrepancies
discrepant <- bootstrap_sig != perm_sig
if (sum(discrepant) > 0) {
  cat("Pathways with discrepant calls (", sum(discrepant), "):\n\n")
  discrepant_df <- data.frame(
    Pathway = pathways_valid[discrepant],
    Pathway_Size = bootstrap_cong$younger_vs_older_Pathway_Size[valid_idx][discrepant],
    Bootstrap_P = round(bootstrap_p[discrepant], 6),
    Perm_P = round(perm_p[discrepant], 6),
    Bootstrap_Sig = bootstrap_sig[discrepant],
    Perm_Sig = perm_sig[discrepant]
  )
  print(discrepant_df, row.names = FALSE)
}

################################################################################
# ERROR METRICS
################################################################################

mean_abs_diff <- mean(abs(bootstrap_p - perm_p))
rmse <- sqrt(mean((bootstrap_p - perm_p)^2))

cat("\n=== ERROR METRICS ===\n")
cat("Mean absolute difference:", round(mean_abs_diff, 6), "\n")
cat("Root mean square error:", round(rmse, 6), "\n\n")

################################################################################
# VISUALIZATION
################################################################################

pdf("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/src/Thien/bootstrap_results/bootstrap_vs_permutation_comparison.pdf", 
    width = 12, height = 10)

par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# 1. Scatter plot - all p-values
plot(perm_p, bootstrap_p, 
     xlab = "Permutation P-value", 
     ylab = "Bootstrap P-value",
     main = paste0("All P-values (n=", length(bootstrap_p), ")\nSpearman r = ", round(cor_spearman, 3)),
     pch = 16, col = rgb(0, 0, 0, 0.3), cex = 0.8)
abline(0, 1, col = "red", lwd = 2)
abline(h = 0.05, lty = 2, col = "blue")
abline(v = 0.05, lty = 2, col = "blue")
legend("topleft", legend = c("Perfect agreement", "α = 0.05"), 
       col = c("red", "blue"), lty = c(1, 2), lwd = c(2, 1))

# 2. Zoomed scatter
small_p_idx <- perm_p < 0.2 & bootstrap_p < 0.2
if (sum(small_p_idx) > 0) {
  plot(perm_p[small_p_idx], bootstrap_p[small_p_idx], 
       xlab = "Permutation P-value", 
       ylab = "Bootstrap P-value",
       main = paste0("Zoomed: P-values < 0.2 (n=", sum(small_p_idx), ")"),
       pch = 16, col = rgb(0, 0, 1, 0.5), cex = 1)
  abline(0, 1, col = "red", lwd = 2)
  abline(h = 0.05, lty = 2, col = "blue")
  abline(v = 0.05, lty = 2, col = "blue")
}

# 3. Bland-Altman plot
mean_p <- (bootstrap_p + perm_p) / 2
diff_p <- bootstrap_p - perm_p
plot(mean_p, diff_p,
     xlab = "Mean of two methods",
     ylab = "Bootstrap - Permutation",
     main = "Bland-Altman Plot",
     pch = 16, col = rgb(0, 0, 0, 0.3))
abline(h = 0, col = "red", lwd = 2)
abline(h = mean(diff_p) + 1.96*sd(diff_p), lty = 2, col = "blue")
abline(h = mean(diff_p) - 1.96*sd(diff_p), lty = 2, col = "blue")
legend("topright", legend = c("Mean bias", "95% limits"), 
       col = c("red", "blue"), lty = c(1, 2))

# 4. Histogram of differences
hist(diff_p, breaks = 30,
     xlab = "Bootstrap P-value - Permutation P-value",
     main = "Distribution of Differences",
     col = "lightblue")
abline(v = 0, col = "red", lwd = 2)
abline(v = mean(diff_p), col = "darkred", lwd = 2, lty = 2)

dev.off()

cat("✓ Plot saved: bootstrap_vs_permutation_comparison.pdf\n")
