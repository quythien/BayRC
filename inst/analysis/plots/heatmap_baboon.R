# Heatmap for concordance of Baboon
library(pheatmap)
current_gtex <- "/home/qtp1/Projects/Collaborative"
current_wd <- "/home/qtp1/Projects/Circadian"
current_aging <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging"



load(file = "/home/qtp1/Projects/Collaborative/GTEXdata/result/summary/hb/mcmc_rho_BF3.RData")
load(file.path(current_gtex, "GTEXdata/result/summary/hb/phi/mcmc_phi_BF3.RData"))  

###############################################################################
# 1. Build baboon MCMC objects
###############################################################################

baboon_tissues <- names(mcmc_data_baboon)

baboon_list <- lapply(baboon_tissues, function(tis) {
  list(
    rho = mcmc_data_baboon[[tis]],
    phi = mcmc_phi_baboon[[tis]]
  )
})
names(baboon_list) <- baboon_tissues

###############################################################################
# 2. Pairwise combinations
###############################################################################
pairwise <- combn(baboon_tissues, 2, simplify = FALSE)

###############################################################################
# 3. Run multi_conservation for each pair
###############################################################################

pair_results <- list()

for (pair in pairwise) {
  
  t1 <- pair[1]
  t2 <- pair[2]
  
  cat("\nRunning:", t1, "vs", t2, "\n")
  
  res <- multi_conservation(
    mcmc.merge.list = list(baboon_list[[t1]], baboon_list[[t2]]),
    dataset.names   = c(t1, t2),
    select.pathway.list = "global",
    n_perm = 1000,
    n_boot = 1000,
    output.dir = file.path(output.dir, "concordance_pairwise"),
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
  AC_list[[nm]] <- df[[colname]]   # this is a numeric vector (length 1 here)
}



saveRDS(
  pair_results,
  file = file.path(output.dir, "concordance_pairwise",  "pairwise_concordance_baboon.rds")
)


###############################################################################
# 5. Build pairwise Adjusted Concordance summary heatmap
###############################################################################

# Initialize empty matrix
ACI_mat <- matrix(NA, nrow = length(baboon_tissues), ncol = length(baboon_tissues))
rownames(ACI_mat) <- baboon_tissues
colnames(ACI_mat) <- baboon_tissues

# Fill symmetric matrix with mean ACI across pathways
for (name in names(AC_list)) {
  parts <- strsplit(name, "_vs_")[[1]]
  t1 <- parts[1]
  t2 <- parts[2]
  
  mean_ACI <- mean(AC_list[[name]], na.rm = TRUE)
  ACI_mat[t1, t2] <- mean_ACI
  ACI_mat[t2, t1] <- mean_ACI
}

diag(ACI_mat) <- 1  # perfect concordance with itself


###############################################################################
# 6. Save PDF heatmap 
###############################################################################

library(pheatmap)

# Your color palette
col_fun <- colorRampPalette(c(
  "#0011FF",   # electric blue
  "#3F00FF",   # violet
  "#7F00FF",   # deep purple
  "#FF4D4D",   # soft red
  "#FF0000"    # pure red for diagonal
))(200)

outdir <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/all_plots"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Convert concordance to dissimilarity
ACI_dissim <- 1 - ACI_mat
diag(ACI_dissim) <- 0  # distance to self = 0

# Create distance objects
row_dist <- as.dist(ACI_dissim)
col_dist <- as.dist(ACI_dissim)

# List of methods to test
methods <- c("ward.D2", "complete", "average")

# Generate a heatmap for each method using dissimilarity
for (method in methods) {
  pdf(file.path(outdir, paste0("Baboon_Concordance_Heatmap_Dissim_", method, ".pdf")),
      width = 9, height = 8)
  
  pheatmap(
    ACI_mat,  # still display original concordance values
    cluster_rows = hclust(row_dist, method = method),
    cluster_cols = hclust(col_dist, method = method),
    color = col_fun,
    main = paste0("Baboon Rhythmicity Concordance Heatmap"),
    fontsize = 10,
    border_color = NA,
    legend = TRUE,
    legend_breaks = c(0, 0.25, 0.5, 0.75, 1),
    legend_labels = c("0", "0.25", "0.50", "0.75", "1.00")
  )
  
  dev.off()
  cat("Saved dissimilarity-based:", method, "\n")
}


# 1. Save the full ACI matrix as CSV
write.csv(
  ACI_mat,
  file = file.path(outdir, "Baboon_Concordance_Matrix.csv"),
  row.names = TRUE
)


##### Modified version 


library(pheatmap)

# Your original color palette
col_fun <- colorRampPalette(c(
  "#0011FF",   # electric blue
  "#3F00FF",   # violet
  "#7F00FF",   # deep purple
  "#FF4D4D",   # soft red
  "#FF0000"    # pure red
))(200)

# 1. Determine a good upper limit for the color scale (e.g., 0.25 or the 95th percentile)
# This ensures the colors are used for the off-diagonal variation.
max_val <- 0.25
# Alternatively, use: max_val <- quantile(ACI_mat[row(ACI_mat) != col(ACI_mat)], 0.99)

# 2. Create breaks that focus on the 0 to max_val range
# We create 201 breaks for 200 colors.
breaksList <- seq(0, max_val, length.out = 201)

# Generate the heatmap
for (method in methods) {
  pdf(file.path(outdir, paste0("Baboon_Concordance_Heatmap_0.25_", method, ".pdf")),
      width = 9, height = 8)
  
  pheatmap(
    ACI_mat,
    cluster_rows = hclust(row_dist, method = method),
    cluster_cols = hclust(col_dist, method = method),
    color = col_fun,
    breaks = breaksList, 
    main = paste0("Baboon Rhythmicity Concordance"),
    fontsize = 10,
    border_color = NA,
    legend = TRUE,
    # Adjust legend labels to match the new scale
    legend_breaks = seq(0, max_val, length.out = 5),
    legend_labels = format(seq(0, max_val, length.out = 5), digits = 2)
  )
  
  dev.off()
}
