# Standalone: circadian KEGG concordance matrices for within-baboon, within-human, and cross-species pairs
# Usage:
#   Rscript heatmap_circadian_pairs.R within_baboon
#   Rscript heatmap_circadian_pairs.R within_baboon_with_scn
#   Rscript heatmap_circadian_pairs.R within_human
#   Rscript heatmap_circadian_pairs.R cross_species

library(pheatmap)

if (!exists("mode", inherits = FALSE)) {
  mode <- "all"
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) mode <- args[1]
}

current_gtex <- "/home/qtp1/Projects/Collaborative"
current_wd <- "/home/qtp1/Projects/Circadian"
output.dir <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/results/baboon/output_final"

load(file = "/home/qtp1/Projects/Collaborative/GTEXdata/result/summary/hb/mcmc_rho_BF3.RData")
load(file.path(current_gtex, "GTEXdata/result/summary/hb/phi/mcmc_phi_BF3.RData"))

# Thien functions
thien_dir <- file.path(current_wd, "Kyle/Circadian-analysis-main/R/v1/BayRC/Thien")
source(file.path(thien_dir, "Permutation_Sim.R"))
Rcpp::sourceCpp(file.path(thien_dir, "congruence.cpp"))

# Pathway list
kegg.pathway.list_hsa <- readRDS(
  "/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/pathway_data/kegg_pathway_list_hsa.rds"
)

pathway_name <- "KEGG Circadian rhythm"
if (!pathway_name %in% names(kegg.pathway.list_hsa)) {
  stop("KEGG Circadian rhythm not found in kegg.pathway.list_hsa")
}

# Replace ARNTL with BMAL1
circadian_genes <- kegg.pathway.list_hsa[[pathway_name]]
circadian_genes <- ifelse(circadian_genes == "ARNTL", "BMAL1", circadian_genes)
kegg.pathway.list_hsa[[pathway_name]] <- circadian_genes

# Build tissue lists (rho/phi objects)
baboon_tissues <- names(mcmc_data_baboon)
human_tissues <- names(mcmc_data_human)

make_tissue_obj <- function(rho_mat, phi_mat) {
  if (is.null(attr(rho_mat, "symbols"))) {
    dn <- dimnames(rho_mat)
    if (!is.null(dn) && length(dn) > 0 && !is.null(dn[[1]])) {
      attr(rho_mat, "symbols") <- dn[[1]]
    }
  }
  list(rho = rho_mat, phi = phi_mat)
}

baboon_list <- lapply(baboon_tissues, function(tis) {
  make_tissue_obj(mcmc_data_baboon[[tis]], mcmc_phi_baboon[[tis]])
})
names(baboon_list) <- baboon_tissues

human_list <- lapply(human_tissues, function(tis) {
  make_tissue_obj(mcmc_data_human[[tis]], mcmc_phi_human[[tis]])
})
names(human_list) <- human_tissues

# Output base
base_path <- file.path(output.dir, "heatmap_circadian_pairs")
dir.create(base_path, recursive = TRUE, showWarnings = FALSE)

run_pairs_to_matrix <- function(tissues, data_list, label, out_prefix, subdir) {
  run_path <- file.path(base_path, subdir)
  dir.create(run_path, recursive = TRUE, showWarnings = FALSE)

  pairwise <- combn(tissues, 2, simplify = FALSE)
  pair_results <- list()

  for (pair in pairwise) {
    t1 <- pair[1]
    t2 <- pair[2]
    cat("\nRunning:", label, t1, "vs", t2, "\n")

    res <- multi_conservation(
      mcmc.merge.list = list(data_list[[t1]], data_list[[t2]]),
      dataset.names   = c(t1, t2),
      select.pathway.list = kegg.pathway.list_hsa[pathway_name],
      n_perm = 1000,
      n_boot = 1000,
      output.dir = file.path(run_path, "concordance_pairwise"),
      use_cpp = TRUE,
      compute_pvalue = FALSE,
      compute_ci = FALSE
    )

    pair_name <- paste0(t1, "_vs_", t2)
    pair_results[[pair_name]] <- res
  }

  # Build ACI matrix
  ACI_mat <- matrix(NA, nrow = length(tissues), ncol = length(tissues))
  rownames(ACI_mat) <- tissues
  colnames(ACI_mat) <- tissues

  for (name in names(pair_results)) {
    parts <- strsplit(name, "_vs_")[[1]]
    t1 <- parts[1]
    t2 <- parts[2]
    colname <- paste0(name, "_AdjustedConcordance")
    mean_ACI <- mean(pair_results[[name]][[colname]], na.rm = TRUE)
    ACI_mat[t1, t2] <- mean_ACI
    ACI_mat[t2, t1] <- mean_ACI
  }

  diag(ACI_mat) <- 1

  saveRDS(pair_results, file = file.path(run_path, paste0(out_prefix, ".rds")))
  write.csv(ACI_mat, file = file.path(run_path, paste0(out_prefix, "_Matrix.csv")), row.names = TRUE)

  # Heatmap
  ACI_dissim <- 1 - ACI_mat
  diag(ACI_dissim) <- 0
  row_dist <- as.dist(ACI_dissim)
  col_dist <- as.dist(ACI_dissim)

  col_fun <- colorRampPalette(c("#0011FF", "#3F00FF", "#7F00FF", "#FF4D4D", "#FF0000"))(200)
  breaksList <- seq(0, 1, length.out = 201)
  methods <- c("ward.D2", "complete", "average")

  for (method in methods) {
    pdf(file.path(run_path, paste0(out_prefix, "_Heatmap_", method, ".pdf")), width = 9, height = 8)
    pheatmap(
      ACI_mat,
      cluster_rows = hclust(row_dist, method = method),
      cluster_cols = hclust(col_dist, method = method),
      color = col_fun,
      breaks = breaksList,
      main = paste0(label, " Circadian Pathway Concordance"),
      fontsize = 10,
      border_color = NA,
      legend = TRUE,
      legend_breaks = c(0, 0.25, 0.5, 0.75, 1),
      legend_labels = c("0", "0.25", "0.50", "0.75", "1.00")
    )
    dev.off()
  }
}

run_cross_species <- function(baboon_tissues_x, human_tissues_x, subdir) {
  run_path <- file.path(base_path, subdir)
  dir.create(run_path, recursive = TRUE, showWarnings = FALSE)

  cross_results <- list()
  ACI_cross <- matrix(NA, nrow = length(baboon_tissues_x), ncol = length(human_tissues_x))
  rownames(ACI_cross) <- baboon_tissues_x
  colnames(ACI_cross) <- human_tissues_x

  for (bt in baboon_tissues_x) {
    for (ht in human_tissues_x) {
    cat("\nRunning: Cross", bt, "vs", ht, "\n")
    res <- multi_conservation(
      mcmc.merge.list = list(baboon_list[[bt]], human_list[[ht]]),
      dataset.names   = c(paste0("Baboon_", bt), paste0("Human_", ht)),
      select.pathway.list = kegg.pathway.list_hsa[pathway_name],
      n_perm = 1000,
      n_boot = 1000,
      output.dir = file.path(run_path, "concordance_cross"),
      use_cpp = TRUE,
      compute_pvalue = FALSE,
      compute_ci = FALSE
    )
    colname <- paste0("Baboon_", bt, "_vs_Human_", ht, "_AdjustedConcordance")
    val <- res[[colname]][1]
    ACI_cross[bt, ht] <- val
    cross_results[[paste0(bt, "_vs_", ht)]] <- res
    }
  }

  saveRDS(cross_results, file = file.path(run_path, "cross_species_circadian_results.rds"))
  write.csv(ACI_cross, file = file.path(run_path, "Cross_Circadian_Concordance_Matrix.csv"), row.names = TRUE)

  # Heatmap (full matrix)
  ACI_dissim <- 1 - ACI_cross
  # Use dist() on rectangular matrix for row/col clustering
  row_dist <- dist(ACI_dissim)
  col_dist <- dist(t(ACI_dissim))

  col_fun <- colorRampPalette(c("#0011FF", "#3F00FF", "#7F00FF", "#FF4D4D", "#FF0000"))(200)
  breaksList <- seq(0, 1, length.out = 201)
  methods <- c("ward.D2", "complete", "average")

  for (method in methods) {
    pdf(file.path(run_path, paste0("Cross_Circadian_Heatmap_", method, ".pdf")), width = 9, height = 8)
    pheatmap(
      ACI_cross,
      cluster_rows = hclust(row_dist, method = method),
      cluster_cols = hclust(col_dist, method = method),
      color = col_fun,
      breaks = breaksList,
      main = "Cross-species Circadian Pathway Concordance",
      fontsize = 10,
      border_color = NA,
      legend = TRUE,
      legend_breaks = c(0, 0.25, 0.5, 0.75, 1),
      legend_labels = c("0", "0.25", "0.50", "0.75", "1.00")
    )
    dev.off()
  }
}

if (mode %in% c("all", "within_baboon")) {
  cat("\n=== Within-baboon pairs (exclude SCN) ===\n")
  within_baboon_tissues <- setdiff(baboon_tissues, "SCN")
  run_pairs_to_matrix(
    within_baboon_tissues,
    baboon_list,
    "Baboon",
    "pairwise_concordance_baboon_circadian",
    "within_baboon"
  )
}

if (mode %in% c("within_baboon_with_scn")) {
  cat("\n=== Within-baboon pairs (include SCN) ===\n")
  within_baboon_tissues <- baboon_tissues
  run_pairs_to_matrix(
    within_baboon_tissues,
    baboon_list,
    "Baboon",
    "pairwise_concordance_baboon_circadian_with_scn",
    "within_baboon_with_scn"
  )
}

if (mode %in% c("all", "within_human")) {
  cat("\n=== Within-human pairs (exclude SCN) ===\n")
  within_human_tissues <- setdiff(human_tissues, "SCN")
  run_pairs_to_matrix(
    within_human_tissues,
    human_list,
    "Human",
    "pairwise_concordance_human_circadian",
    "within_human"
  )
}

if (mode %in% c("all", "cross_species")) {
  cat("\n=== Cross-species pairs (exclude SCN) ===\n")
  baboon_x <- setdiff(baboon_tissues, "SCN")
  human_x <- setdiff(human_tissues, "SCN")
  run_cross_species(baboon_x, human_x, "cross_species")
}

cat("\nDONE: circadian concordance outputs in ", base_path, "\n")
