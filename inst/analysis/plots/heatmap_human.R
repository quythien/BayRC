# Heatmap for concordance of Human tissues (pairwise)
library(pheatmap)
current_gtex <- "/home/qtp1/Projects/Collaborative"
current_wd <- "/home/qtp1/Projects/Circadian"
current_aging <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging"
output.dir <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/results/baboon/output_final"

load(file = "/home/qtp1/Projects/Collaborative/GTEXdata/result/summary/hb/mcmc_rho_BF3.RData")
load(file.path(current_gtex, "GTEXdata/result/summary/hb/phi/mcmc_phi_BF3.RData"))

# Thien functions
thien_dir <- file.path(current_wd, "Kyle/Circadian-analysis-main/R/v1/BayRC/Thien")
source(file.path(thien_dir, "Permutation_Sim.R"))
Rcpp::sourceCpp(file.path(thien_dir, "congruence.cpp"))

###############################################################################
# 1. Build human MCMC objects
###############################################################################

human_tissues <- names(mcmc_data_human)

human_list <- lapply(human_tissues, function(tis) {
  list(
    rho = mcmc_data_human[[tis]],
    phi = mcmc_phi_human[[tis]]
  )
})
names(human_list) <- human_tissues

###############################################################################
# 2. Pairwise combinations
###############################################################################
pairwise <- combn(human_tissues, 2, simplify = FALSE)

###############################################################################
# 3. Run multi_conservation for each pair
###############################################################################

pair_results <- list()

for (pair in pairwise) {
  t1 <- pair[1]
  t2 <- pair[2]
  cat("\nRunning:", t1, "vs", t2, "\n")

  res <- multi_conservation(
    mcmc.merge.list = list(human_list[[t1]], human_list[[t2]]),
    dataset.names   = c(t1, t2),
    select.pathway.list = "global",
    n_perm = 1000,
    n_boot = 1000,
    output.dir = file.path(output.dir, "concordance_pairwise_human"),
    use_cpp = TRUE,
    compute_pvalue = FALSE,
    compute_ci = FALSE
  )

  pair_name <- paste0(t1, "_vs_", t2)
  pair_results[[pair_name]] <- res
}

###############################################################################
# 4. Extract Adjusted Concordance matrices for all pairs
###############################################################################

AC_list <- list()

for (nm in names(pair_results)) {
  df <- pair_results[[nm]]
  colname <- paste0(nm, "_AdjustedConcordance")
  AC_list[[nm]] <- df[[colname]]
}

saveRDS(
  pair_results,
  file = file.path(output.dir, "concordance_pairwise_human", "pairwise_concordance_human.rds")
)

###############################################################################
# 5. Build pairwise Adjusted Concordance summary heatmap
###############################################################################

ACI_mat <- matrix(NA, nrow = length(human_tissues), ncol = length(human_tissues))
rownames(ACI_mat) <- human_tissues
colnames(ACI_mat) <- human_tissues

for (name in names(AC_list)) {
  parts <- strsplit(name, "_vs_")[[1]]
  t1 <- parts[1]
  t2 <- parts[2]
  mean_ACI <- mean(AC_list[[name]], na.rm = TRUE)
  ACI_mat[t1, t2] <- mean_ACI
  ACI_mat[t2, t1] <- mean_ACI
}

diag(ACI_mat) <- 1

###############################################################################
# 6. Save PDF heatmap
###############################################################################

col_fun <- colorRampPalette(c(
  "#0011FF",
  "#3F00FF",
  "#7F00FF",
  "#FF4D4D",
  "#FF0000"
))(200)

outdir <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/all_plots"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ACI_dissim <- 1 - ACI_mat
diag(ACI_dissim) <- 0

row_dist <- as.dist(ACI_dissim)
col_dist <- as.dist(ACI_dissim)

methods <- c("ward.D2", "complete", "average")

for (method in methods) {
  pdf(file.path(outdir, paste0("Human_Concordance_Heatmap_Dissim_", method, ".pdf")),
      width = 9, height = 8)

  pheatmap(
    ACI_mat,
    cluster_rows = hclust(row_dist, method = method),
    cluster_cols = hclust(col_dist, method = method),
    color = col_fun,
    main = paste0("Human Rhythmicity Concordance Heatmap"),
    fontsize = 10,
    border_color = NA,
    legend = TRUE,
    legend_breaks = c(0, 0.25, 0.5, 0.75, 1),
    legend_labels = c("0", "0.25", "0.50", "0.75", "1.00")
  )

  dev.off()
  cat("Saved dissimilarity-based:", method, "\n")
}

write.csv(
  ACI_mat,
  file = file.path(outdir, "Human_Concordance_Matrix.csv"),
  row.names = TRUE
)
