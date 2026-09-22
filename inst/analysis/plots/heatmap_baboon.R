# Genome-wide pairwise concordance across baboon tissues: writes the pair
# results, the concordance matrix and heatmaps under three linkage methods.
library(pheatmap)
library(BayRC)
# Paths come from inst/analysis/config.R; override any of them with the
# matching env var.
this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
source(file.path(analysis.dir, "plots", "theme_bayrc.R"))
source(file.path(analysis.dir, "plots", "palette_concordance.R"))

current_gtex <- BAYRC_DATA_DIR
current_wd <- BAYRC_WD_DIR
current_aging <- BAYRC_AGING_DIR
output.dir <- BAYRC_OUTPUT_DIR

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))

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
  AC_list[[nm]] <- df[[colname]]   # one value per pathway; one here
}

saveRDS(
  pair_results,
  file = file.path(output.dir, "concordance_pairwise",  "pairwise_concordance_baboon.rds")
)

###############################################################################
# 5. Build pairwise Adjusted Concordance summary heatmap
###############################################################################

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

diag(ACI_mat) <- 1

###############################################################################
# 6. Save PDF heatmap 
###############################################################################

library(pheatmap)

col_fun <- concordance_colors

outdir <- BAYRC_FIGURE_DIR
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# tissues are clustered on 1 - concordance
ACI_dissim <- 1 - ACI_mat
diag(ACI_dissim) <- 0

row_dist <- as.dist(ACI_dissim)
col_dist <- as.dist(ACI_dissim)

methods <- c("ward.D2", "complete", "average")

for (method in methods) {
  cairo_pdf(file.path(outdir, paste0("Baboon_Concordance_Heatmap_Dissim_", method, ".pdf")),
            width = 9, height = 8, family = bayrc_family)
  
  pheatmap(
    ACI_mat,  # cells show concordance, not the distance
    cluster_rows = hclust(row_dist, method = method),
    cluster_cols = hclust(col_dist, method = method),
    color = col_fun,
    main = paste0("Baboon Rhythmicity Concordance Heatmap"),
    fontsize = bayrc_heat_args()$fontsize,
    border_color = bayrc_heat_args()$border_color,
    fontfamily = bayrc_family,
    legend = TRUE,
    legend_breaks = c(0, 0.25, 0.5, 0.75, 1),
    legend_labels = c("0", "0.25", "0.50", "0.75", "1.00")
  )
  
  dev.off()
  cat("Saved dissimilarity-based:", method, "\n")
}

write.csv(
  ACI_mat,
  file = file.path(outdir, "Baboon_Concordance_Matrix.csv"),
  row.names = TRUE
)

##### The same heatmaps on the capped scale

library(pheatmap)

# ramp, cap and breaks come from palette_concordance.R
col_fun <- concordance_colors
max_val <- concordance_max
breaksList <- concordance_breaks

for (method in methods) {
  cairo_pdf(file.path(outdir, paste0("Baboon_Concordance_Heatmap_0.5_", method, ".pdf")),
            width = 9, height = 8, family = bayrc_family)
  
  pheatmap(
    ACI_mat,
    cluster_rows = hclust(row_dist, method = method),
    cluster_cols = hclust(col_dist, method = method),
    color = col_fun,
    breaks = breaksList, 
    main = paste0("Baboon Rhythmicity Concordance"),
    fontsize = bayrc_heat_args()$fontsize,
    border_color = bayrc_heat_args()$border_color,
    fontfamily = bayrc_family,
    legend = TRUE,
    legend_breaks = seq(0, max_val, length.out = 5),
    legend_labels = format(seq(0, max_val, length.out = 5), digits = 2)
  )
  
  dev.off()
}
