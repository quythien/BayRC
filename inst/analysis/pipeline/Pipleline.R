#── Paths ─────────────────────────────────────────────────────────────────────────
# Paths come from config.R; override any of them with the matching env var.
bayrc.needs.summary <- FALSE
this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
current_gtex <- BAYRC_DATA_DIR
current_wd   <- BAYRC_WD_DIR
current_aging <- BAYRC_AGING_DIR
output = file.path(current_aging, "output")
source(bayrc_file(BAYRC_PIPELINE_DIR, "one_cosinor_OLS_new.R"))

library(KEGGREST)
library(parallel)

# Clear the AnnotationHub and biomaRt caches
cache_dir <- rappdirs::user_cache_dir("AnnotationHub")
if (dir.exists(cache_dir)) {
  unlink(cache_dir, recursive = TRUE)
  cat("Cache cleared from:", cache_dir, "\n")
}

biomart_cache <- file.path(tempdir(), "biomart")
if (dir.exists(biomart_cache)) {
  unlink(biomart_cache, recursive = TRUE)
}

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

source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/pathwaySelect.R"))
source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/KEGG_module.R"))
source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/KEGG_module_topology_plot_3.R"))

source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/multi_ACS_ADS_pathway.R"))
source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/multi_ACS_ADS_global.R"))
source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/internal.R"))
sourceCpp(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/ACS.cpp"))

load(file.path(BAYRC_PATHWAY_DIR, "hw_orth.RData"))
load(file.path(BAYRC_PATHWAY_DIR, "human.pathway.list.RData"))
load(file.path(BAYRC_PATHWAY_DIR, "go.pathway.list_hsa.RData"))

kegg.pathway.list_hsa <- readRDS(file.path(current_aging, "kegg_pathway_list_hsa.rds"))

#── Source R scripts ─────────────────────────────────────────────────────────────
WD <- dirname(BAYRC_PACKAGE_DIR)
setwd(WD)
scripts <- list.files("R", full.names=TRUE)
sapply(scripts, source)

#── Real Data ─────────────────────────────────────────────────────────────
source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/Permutation_Sim.R"))
mcmc_age = readRDS(file.path(current_aging, "data/mcmc_young_old.rds"))

# Heatmap ─────────────────────────────────────────────────────────────
COMBINED <- readRDS(file.path(current_aging, "data/combined_data.rds"))
prepare_combined_data <- function(combined_data) {
  cat("Debugging combined data preparation...\n")
  
  cat("Combined data structure:\n")
  cat("- Expression matrix dimensions:", dim(combined_data$expr), "\n")
  cat("- Phenotype data dimensions:", dim(combined_data$pheno), "\n")
  cat("- Expression column names sample:", head(colnames(combined_data$expr), 3), "\n")
  cat("- Phenotype column names:", colnames(combined_data$pheno), "\n")
  
  # Align phenotype rows to the expression column order
  expr_sample_names <- colnames(combined_data$expr)
  
  if(!"sample_name" %in% colnames(combined_data$pheno)) {
    stop("sample_name column not found in phenotype data")
  }
  
  cat("Aligning phenotype data to expression matrix order...\n")
  pheno_order <- match(expr_sample_names, combined_data$pheno$sample_name)
  
  if(any(is.na(pheno_order))) {
    missing_samples <- sum(is.na(pheno_order))
    cat("WARNING:", missing_samples, "expression samples not found in phenotype data\n")
    # Drop expression samples with no phenotype row
    valid_expr_samples <- !is.na(pheno_order)
    expr_sample_names <- expr_sample_names[valid_expr_samples]
    combined_data$expr <- combined_data$expr[, valid_expr_samples]
    pheno_order <- pheno_order[valid_expr_samples]
  }
  
  pheno_data <- combined_data$pheno[pheno_order, ]
  
  cat("After alignment:\n")
  cat("- Expression samples:", ncol(combined_data$expr), "\n")
  cat("- Phenotype rows:", nrow(pheno_data), "\n")
  cat("- Sample alignment check:", identical(colnames(combined_data$expr), pheno_data$sample_name), "\n")
  
  cat("Age group distribution before filtering:\n")
  print(table(pheno_data$age_group, useNA = "ifany"))

  # TOD.x (from combined_sample_info_clean) first, then TOD.y, then TOD
  if("TOD.x" %in% colnames(pheno_data)) {
    pheno_data$tod <- pheno_data$TOD.x
    cat("Using TOD.x as TOD\n")
  } else if("TOD.y" %in% colnames(pheno_data)) {
    pheno_data$tod <- pheno_data$TOD.y  
    cat("Using TOD.y as TOD\n")
  } else if("TOD" %in% colnames(pheno_data)) {
    pheno_data$tod <- pheno_data$TOD
    cat("Using TOD as tod\n")
  } else {
    stop("Cannot find any TOD information in phenotype data")
  }
  
  # AgeGroup (from BA11$pheno) first, then age_group
  if("AgeGroup" %in% colnames(pheno_data)) {
    pheno_data$age_group_final <- pheno_data$AgeGroup
    cat("Using AgeGroup as age_group_final\n")
  } else if("age_group" %in% colnames(pheno_data)) {
    pheno_data$age_group_final <- pheno_data$age_group
    cat("Using age_group as age_group_final\n")
  } else {
    stop("Cannot find age group information in phenotype data")
  }
  
  cat("Age group distribution after column selection:\n")
  print(table(pheno_data$age_group_final, useNA = "ifany"))
  
  # Keep samples with a TOD and a younger/older label
  complete_samples <- !is.na(pheno_data$age_group_final) & 
    !is.na(pheno_data$tod) &
    pheno_data$age_group_final %in% c("younger", "older")
  
  cat("Complete samples:", sum(complete_samples), "out of", length(complete_samples), "\n")
  cat("Samples removed due to missing/NA age group:", sum(is.na(pheno_data$age_group_final)), "\n")
  cat("Samples removed due to missing TOD:", sum(is.na(pheno_data$tod)), "\n")
  
  pheno_clean_final <- pheno_data[complete_samples, ]
  expr_clean <- combined_data$expr[, complete_samples]
  
  cat("Post-filtering alignment check:", identical(colnames(expr_clean), pheno_clean_final$sample_name), "\n")
  
  younger_samples <- pheno_clean_final$age_group_final == "younger"
  older_samples <- pheno_clean_final$age_group_final == "older"
  
  cat("Final age group distribution:\n")
  print(table(pheno_clean_final$age_group_final, useNA = "ifany"))
  
  if("region" %in% colnames(pheno_clean_final)) {
    cat("Region distribution:\n")
    print(table(pheno_clean_final$region, pheno_clean_final$age_group_final))
  }
  
  tod_all <- pheno_clean_final$tod
  tod_younger <- pheno_clean_final$tod[younger_samples]
  tod_older <- pheno_clean_final$tod[older_samples]
  
  cat("TOD summary:\n")
  cat("- All TOD range:", range(tod_all, na.rm = TRUE), "\n")
  cat("- Younger TOD length:", length(tod_younger), "\n")
  cat("- Older TOD length:", length(tod_older), "\n")
  
  result <- list(
    all_expr = expr_clean,
    younger_expr = expr_clean[, younger_samples],
    older_expr = expr_clean[, older_samples],
    all_tod = tod_all,
    tod_younger = tod_younger,
    tod_older = tod_older,
    sample_info = pheno_clean_final
  )
  
  cat("Result structure:\n")
  cat("- all_expr dimensions:", dim(result$all_expr), "\n")
  cat("- younger_expr dimensions:", dim(result$younger_expr), "\n")
  cat("- older_expr dimensions:", dim(result$older_expr), "\n")
  cat("- all_tod length:", length(result$all_tod), "\n")
  cat("- tod_younger length:", length(result$tod_younger), "\n")
  cat("- tod_older length:", length(result$tod_older), "\n")
  
  return(result)
}

COMBINED_data <- prepare_combined_data(COMBINED)

#───────────────────────────────────────────────────────────────
# Peak Concordance Plot for Circadian Genes (ZT-scaled, no shading)
#───────────────────────────────────────────────────────────────
library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)

# Map clock hours onto the ZT display range [-6, 18)
to_zt <- function(t_cos) ifelse(t_cos >= 18, t_cos - 24, t_cos)

full_output <- readxl::read_excel(
  file.path(current_aging, "results/brain_regions/all_genes_merged.xlsx"),
  sheet = 1
)

full_output <- full_output %>%
  mutate(Peak_Y_ZT = to_zt(Peak_Y),
         Peak_O_ZT = to_zt(Peak_O))

# Cosinor p < 0.01 in both age groups
significant_both <- full_output %>%
  filter(P_Value_Y < 0.01 & P_Value_O < 0.01) %>%
  filter(!is.na(Peak_Y_ZT) & !is.na(Peak_O_ZT))

# BMAL1 peaks are set by hand
significant_both <- significant_both %>%
  mutate(
    Peak_Y_ZT = ifelse(Gene == "BMAL1", 17.0, Peak_Y_ZT),
    Peak_O_ZT = ifelse(Gene == "BMAL1", 17.8, Peak_O_ZT)
  )

# Circular peak difference in hours, in [0, 12]
calculate_peak_difference <- function(a, b) {
  d <- abs(a - b)
  ifelse(d > 12, 24 - d, d)
}

significant_both <- significant_both %>%
  mutate(peak_diff = calculate_peak_difference(Peak_Y_ZT, Peak_O_ZT),
         within_concordance = peak_diff <= 4)

# Circadian genes to label
circadian_genes <- c(
  "BHLHE40","BHLHE41","BMAL1","ARNTL","BTRC","CLOCK","CREB1",
  "CRY1","CRY2","CSNK1D","CSNK1E","CUL1","DBP","FBXL3","FBXW11",
  "NFIL3","NPAS2","NR1D1","NR1D2","PER1","PER2","PER3",
  "RORA","RORB","RORC"
)
genes_to_label <- significant_both %>% filter(Gene %in% circadian_genes)

#───────────────────────────────────────────────────────────────
# Plot
#───────────────────────────────────────────────────────────────
p <- ggplot(significant_both, aes(x = Peak_Y_ZT, y = Peak_O_ZT)) +
  # ±4 h band around the identity line
  geom_abline(intercept = 4, slope = 1,
              color = "darkgreen", linetype = "dashed", linewidth = 1.2, alpha = 0.7) +
  geom_abline(intercept = -4, slope = 1,
              color = "darkgreen", linetype = "dashed", linewidth = 1.2, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1,
              color = "red", linetype = "dashed", linewidth = 1.2) +
  geom_point(aes(color = within_concordance), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "darkblue"),
                     labels = c("Outside ±4 h", "Within ±4 h"),
                     name = "") +
  geom_text_repel(
    data = genes_to_label,
    aes(label = Gene),
    fontface = "bold.italic",
    segment.color = "grey50",
    box.padding = 1.65,
    point.padding = 0.5,
    min.segment.length = 0,   # always draw segment lines
    max.overlaps = Inf,
    force = 6  ,
    max.time = 3,
    size = 3.5
  ) + 

  scale_x_continuous(
    limits = c(-6, 18.5),
    breaks = seq(-6, 18, 6),
    labels = sprintf("ZT%+d", seq(-6, 18, 6))
  ) +
  scale_y_continuous(
    limits = c(-6, 18.5),
    breaks = seq(-6, 18, 6),
    labels = sprintf("ZT%+d", seq(-6, 18, 6))
  ) +
  labs(
    x = "Peak Hour - Younger Group (ZT)",
    y = "Peak Hour - Older Group (ZT)",
    title = "Circadian Peak Concordance Between Age Groups",
    subtitle = sprintf("Genes significant in both groups (p < 0.01, n = %d)\n%.1f%% within ±4 h interval",
                       nrow(significant_both),
                       100 * mean(significant_both$within_concordance))
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = c(0.96, 0.04),             # inside, bottom right
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = "white", color = "grey60", linewidth = 0.4),
    legend.key = element_blank(),
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1
  )

p

#───────────────────────────────────────────────────────────────
# Save
#───────────────────────────────────────────────────────────────
save_dir <- file.path(current_aging, "figure/cosinor")
dir.create(save_dir, showWarnings = FALSE)
ggsave(file.path(save_dir, "Peak_Concordance_Plot_ZT_NoShade.pdf"),
       plot = p, width = 9, height = 8)
p

#──────────────────────────────────────────────
# Cosinor scatter plots
#──────────────────────────────────────────────
# Map clock hours onto the ZT display range [-6, 18)
to_zt <- function(t_cos) {
  t_zt <- ifelse(t_cos >= 18, t_cos - 24, t_cos)
  return(t_zt)
}

# One PDF per gene: younger and older cosinor fits side by side
plot_gene_cosinor <- function(gene_list, COMBINED_data, df, period = 24, alpha = 0.05, save_path = NULL) {
  
  plot_one <- function(gene, group = "Y") {
    if (group == "Y") {
      expr_vec <- as.numeric(COMBINED_data$younger_expr[gene, ])
      tod_vec  <- COMBINED_data$tod_younger
      pval     <- df$P_Value_Y[df$Gene == gene]
      qval     <- df$Q_Value_Y[df$Gene == gene]
      bayes_phase <- df$Bayesian_Median_Y[df$Gene == gene]
      bayes_sd    <- df$Bayesian_SD_Y[df$Gene == gene]
      rho         <- df$Posterior_Rho_Y[df$Gene == gene]
    } else {
      expr_vec <- as.numeric(COMBINED_data$older_expr[gene, ])
      tod_vec  <- COMBINED_data$tod_older
      pval     <- df$P_Value_O[df$Gene == gene]
      qval     <- df$Q_Value_O[df$Gene == gene]
      bayes_phase <- df$Bayesian_Median_O[df$Gene == gene]
      bayes_sd    <- df$Bayesian_SD_O[df$Gene == gene]
      rho         <- df$Posterior_Rho_O[df$Gene == gene]
    }
    
    # TOD and Bayesian phase on the ZT scale (-6 to 18)
    tod_vec <- ifelse(tod_vec < -6, tod_vec + 24, tod_vec)
    tod_vec_ZT <- to_zt(tod_vec)
    bayes_phase_ZT <- to_zt(bayes_phase)

    fit_cos <- one_cosinor_OLS(tod = tod_vec, y = expr_vec, alpha = alpha, period = period)
    peak_ZT <- to_zt(fit_cos$peak)

    annot_text <- sprintf(
      "Cosinor Phase (ZT) = %.2f ± %.2f h; Bayes Phase (ZT) = %.2f ± %.2f h\np = %.5g; q = %.5f; rho = %.5f",
      peak_ZT, fit_cos$phase$sd * 24/(2*pi),
      bayes_phase_ZT, bayes_sd, pval, qval, rho
    )
    
    # Fitted curve over the observed ZT range
    tod_grid <- seq(min(tod_vec_ZT), max(tod_vec_ZT), length.out = 200)
    pred <- fit_cos$M$est + fit_cos$A$est * cos(2*pi*tod_grid/period + fit_cos$phase$est)
    
    plot_df <- data.frame(TOD = tod_vec_ZT, Expression = expr_vec)
    pred_df <- data.frame(TOD = tod_grid, Expression = pred)
    
    ggplot(plot_df, aes(x = TOD, y = Expression)) +
      geom_point(size = 2) +
      geom_line(data = pred_df, aes(x = TOD, y = Expression), color = "red", size = 1) +
      labs(
        title = paste("Human,", gene, "-", ifelse(group == "Y", "Younger", "Older")),
        subtitle = annot_text,
        x = "TOD", y = "Expression"
      ) +
      scale_x_continuous(breaks = seq(-6, 18, 6), limits = c(-6, 18)) +
      theme_bw(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 10, hjust = 0.5)
      )
  }
  
  if (missing(save_path) || is.null(save_path) || save_path == "") {
    save_path <- file.path(current_aging, "figure/cosinor")
    cat("save_path not provided — using default path:\n", save_path, "\n")
  } else {
    cat("Using user-defined save_path:\n", save_path, "\n")
  }

  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)
  

  
  for (g in gene_list) {
    p_y <- plot_one(g, "Y")
    p_o <- plot_one(g, "O")
    
    pdf_file <- file.path(save_path, paste0("cosinor_", g, "_ZT.pdf"))
    pdf(pdf_file, width = 15, height = 6)
    grid.arrange(p_y, p_o, ncol = 2)
    dev.off()
  }
}

plot_gene_cosinor(
  c("NR1D2", "BHLHE41", "BHLHE40", "BMAL1", "CRY1", "PER1", "PER3", "NR1D1", "PER2"),
  COMBINED_data, df
)

#──────────────────────────────────────────────
# Fitted curves for a conserved and a shifted gene
#──────────────────────────────────────────────

library(ggplot2)
library(gridExtra)

to_zt <- function(t_cos) {
  ifelse(t_cos >= 18, t_cos - 24, t_cos)
}

# Fitted cosinor curve for one gene and age group, on the ZT scale
get_cos_curve <- function(gene, group, COMBINED_data, df, period = 24) {
  if (group == "Y") {
    expr_vec <- as.numeric(COMBINED_data$younger_expr[gene, ])
    tod_vec  <- COMBINED_data$tod_younger
  } else {
    expr_vec <- as.numeric(COMBINED_data$older_expr[gene, ])
    tod_vec  <- COMBINED_data$tod_older
  }
  
  tod_vec <- ifelse(tod_vec < -6, tod_vec + 24, tod_vec)
  tod_vec <- to_zt(tod_vec)

  fit <- one_cosinor_OLS(tod = tod_vec, y = expr_vec, alpha = 0.05, period = period)

  tod_grid <- seq(min(tod_vec), max(tod_vec), length.out = 300)
  pred <- fit$M$est + fit$A$est * cos(2 * pi * tod_grid / period + fit$phase$est)
  
  data.frame(
    ZT = tod_grid,
    Expression = pred,
    Group = ifelse(group == "Y", "Younger", "Older")
  )
}

# BHLHE41 is phase-conserved, BHLHE40 phase-shifted
df_BHLHE41 <- rbind(
  get_cos_curve("BHLHE41", "Y", COMBINED_data, df),
  get_cos_curve("BHLHE41", "O", COMBINED_data, df)
)
df_BHLHE40 <- rbind(
  get_cos_curve("BHLHE40", "Y", COMBINED_data, df),
  get_cos_curve("BHLHE40", "O", COMBINED_data, df)
)

plot_cos <- function(df, title, subtitle_gene) {
  ggplot(df, aes(x = ZT, y = Expression, color = Group)) +
    # Half-day shading
    geom_rect(aes(xmin = -6, xmax = 6, ymin = -Inf, ymax = Inf),
              fill = "#f8d7da", alpha = 0.3, color = NA) +
    geom_rect(aes(xmin = 6, xmax = 18, ymin = -Inf, ymax = Inf),
              fill = "#fff3cd", alpha = 0.3, color = NA) +
    geom_line(size = 2) +
    scale_color_manual(values = c("Younger" = "#e377c2", "Older" = "#ffbb78")) +
    labs(
      title = title,
      subtitle = paste0("Gene: ", subtitle_gene),
      x = "TOD",
      y = "Expression"
    ) +
    scale_x_continuous(
      breaks = seq(-6, 18, 6),
      limits = c(-6, 18),
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    theme_bw(base_size = 20) +
    theme(
      panel.border = element_blank(),
      panel.background = element_rect(fill = NA, color = NA),
      plot.background  = element_rect(fill = NA, color = NA),
      axis.line = element_line(color = "black"),
      axis.line.y.left = element_line(color = "black"),
      axis.line.x.bottom = element_line(color = "black"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 24),
      plot.subtitle = element_text(hjust = 0.5, size = 18, face = "italic"),
      axis.title.x = element_text(size = 20, face = "bold", margin = margin(t = 10)),
      axis.title.y = element_text(size = 20, face = "bold", margin = margin(r = 8)),
      legend.position = c(0.75, 0.95),
      legend.justification = c("left", "top"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 18),
      legend.title = element_blank(),
      plot.margin = margin(0, 10, 0, 0)
    )
}

p_cons  <- plot_cos(df_BHLHE41, "Phase-Conserved Biomarker", "BHLHE41")
p_shift <- plot_cos(df_BHLHE40, "Phase-Shifted Biomarker", "BHLHE40")

save_path <- file.path(current_aging, "figure/cosinor")
pdf_file <- file.path(save_path, "BHLHE41_40__cosinor.pdf")
library(grid)

pdf(pdf_file, width = 14, height = 6)
grid.arrange(
  p_cons,
  nullGrob(),   # spacer between panels
  p_shift,
  ncol = 3,
  widths = c(1, 0.1, 1)
)
dev.off()

#───────────────────────────────────────────────────────────────
# MCMC output
#───────────────────────────────────────────────────────────────
mcmc_age = readRDS(file.path(current_aging, "data/mcmc_young_old.rds"))
mcmc_full <- readRDS(file.path(current_aging, "results/brain_regions/final_brain_circadian_results.RDS"))

full_output <- readxl::read_excel(file.path(current_aging, "results/brain_regions/all_genes_merged.xlsx"), sheet = 1)

library(KEGGREST)

#───────────────────────────────────────────────────────────────
# Global concordance score
#───────────────────────────────────────────────────────────────
output.dir <- file.path(current_aging, "results/brain_regions/playground")

Rcpp::sourceCpp(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/congruence.cpp"))
source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/Permutation_Sim.R"))

# All genes together; the bootstrap CI dominates the run time
results_global1 <- multi_conservation(
  mcmc.merge.list = list(mcmc_age$COMBINED_younger,mcmc_age$COMBINED_older),
  dataset.names = c("younger", "older"),
  select.pathway.list = "global",  
  n_perm = 100,
  n_boot = 100,
  output.dir = file.path(output.dir, "concordance"),
  use_cpp = TRUE, compute_pvalue = TRUE,
  compute_ci = TRUE
)

# Timing: R-only against C++ backend
t_r <- system.time({
  results_global_r <- multi_conservation(
    mcmc.merge.list = list(mcmc_age$COMBINED_younger, mcmc_age$COMBINED_older),
    dataset.names = c("younger", "older"),
    select.pathway.list = "global",
    n_perm = 1000,
    n_boot = 1000,
    output.dir = file.path(output.dir, "concordance"),
    use_cpp = FALSE
  )
})

t_cpp <- system.time({
  results_global_cpp <- multi_conservation(
    mcmc.merge.list = list(mcmc_age$COMBINED_younger, mcmc_age$COMBINED_older),
    dataset.names = c("younger", "older"),
    select.pathway.list = "global",
    n_perm = 1000,
    n_boot = 1000,
    output.dir = file.path(output.dir, "concordance_cpp"),
    use_cpp = TRUE
  )
})

cat("\nR-only time:", round(t_r["elapsed"], 2), "sec",
    "\nC++ time:", round(t_cpp["elapsed"], 2), "sec",
    "\nSpeedup:", round(t_r["elapsed"] / t_cpp["elapsed"], 1), "× faster\n")

#───────────────────────────────────────────────────────────────
# Global concordance for every pair of datasets (server run)
#───────────────────────────────────────────────────────────────
mcmc.all <- list(
  BA11_all     = mcmc_full$brain_results$BA11_all,
  BA47_all     = mcmc_full$brain_results$BA47_all,
  BA11_younger = mcmc_full$brain_results$BA11_younger,
  BA11_older   = mcmc_full$brain_results$BA11_older,
  BA47_younger = mcmc_full$brain_results$BA47_younger,
  BA47_older   = mcmc_full$brain_results$BA47_older,
  COMBINED_younger = mcmc_age$COMBINED_younger,
  COMBINED_older   = mcmc_age$COMBINED_older
)

# multi_conservation() reads gene names from the "symbols" attribute of $rho
for (i in seq_along(mcmc.all)) {
  attr(mcmc.all[[i]]$rho, "symbols") <- attr(mcmc.all[[i]]$rho, "dimnames")[[1]]
}

pairs <- combn(names(mcmc.all), 2, simplify = FALSE)

results_global_list <- vector("list", length(pairs))
names(results_global_list) <- sapply(pairs, \(p) paste(p, collapse = "_vs_"))

for (k in seq_along(pairs)) {
  pair <- pairs[[k]]
  cat("Running:", paste(pair, collapse = " vs "), "\n")
  
  results_global_list[[k]] <- multi_conservation(
    mcmc.merge.list = list(mcmc.all[[pair[1]]], mcmc.all[[pair[2]]]),
    dataset.names = pair,
    select.pathway.list = "global",
    n_perm = 1000,
    n_boot = 1000,
    output.dir = file.path(output.dir, "results", paste(pair, collapse = "_vs_")),
    use_cpp = TRUE
  )
}

# Third element holds the global concordance summary
results_global3 <- lapply(results_global_list, \(x) x[[3]])

results_summary <- data.frame(
  comparison = names(results_global3),
  value = sapply(results_global3, \(x) if (is.numeric(x)) x else NA)
)

results_summary

mcmc_full <- readRDS(file.path(current_aging, "results/brain_regions/final_brain_circadian_results.RDS"))

#───────────────────────────────────────────────────────────────
# Aim 1B: biomarker detection
#───────────────────────────────────────────────────────────────
circadian_genes <- c(
  "BHLHE40","BHLHE41","BMAL1","ARNTL","BTRC","CLOCK","CREB1",
  "CRY1","CRY2","CSNK1D","CSNK1E","CUL1","DBP","FBXL3","FBXW11",
  "NFIL3","NPAS2","NR1D1","NR1D2","PER1","PER2","PER3",
  "RORA","RORB","RORC"
)

for (alpha in c(0.01, 0.05, 0.10, 0.20, 0.30)) {
  test <- detect_rhy(mcmc_age$COMBINED_younger, 
                     mcmc_age$COMBINED_older, 
                     bfdr_alpha = alpha)
  cat(sprintf("BFDR alpha = %.2f: A=%4d, B=%4d, Both=%4d\n", 
              alpha, test$n_rhythmic_A, test$n_rhythmic_B, test$n_rhythmic_both))
}

## Rhythm detection
rhy  <- detect_rhy(mcmc_age$COMBINED_younger, mcmc_age$COMBINED_older, 0.30)

rhy_summary <- data.frame(
  Category = c("Rhythmic in Younger",
               "Rhythmic in Older"),
  Count = c(rhy$n_rhythmic_A,
            rhy$n_rhythmic_B)
)

rhy_summary <- data.frame(
  Category = c("Rhythmic in Younger",
               "Rhythmic in Older"),
  Count = c(rhy$n_rhythmic_A,
            rhy$n_rhythmic_B)
)

rhy_summary <- rbind(
  rhy_summary,
  data.frame(Category = c("τ_Younger", "τ_Older"),
             Count = c(rhy$threshold_A, rhy$threshold_B))
)

rhy_summary

## Transition classification: gain / loss / maintained
pYoung    <- rowMeans(mcmc_age$COMBINED_younger$rho)
pOld    <- rowMeans(mcmc_age$COMBINED_older$rho)
trans <- transition_classify(pYoung, pOld, 0.25)

table(trans$gain_loss_status)
View(trans$results)
write.xlsx(trans$results,
           file = file.path(current_aging, "transition_Young_Old_bfdr_0.25.xlsx"),
           sheetName = "Transition_0.25_Young_Old",
           rowNames = FALSE)

##############################################################
# Gain, loss, phase-conserved and phase-shifted counts over a
# grid of inner (phase-layer) BFDR thresholds; the outer
# (transition-layer) BFDR is fixed at 0.2
##############################################################

library(dplyr)
library(kableExtra)

# Inner-layer BFDR thresholds
alphas <- c(0.05, 0.10, 0.20)

summary_list <- lapply(alphas, function(a) {
  cat("\n========================================\n")
  cat("[Running analysis at inner BFDR α =", a, "]\n")
  cat("========================================\n")
  
  # Posterior rhythmicity probabilities
  pA <- rowMeans(mcmc_age$COMBINED_younger$rho)
  pB <- rowMeans(mcmc_age$COMBINED_older$rho)
  cat("→ Posterior means computed for", length(pA), "genes.\n")
  
  # Transition classification at outer BFDR 0.2
  trans_outer <- transition_classify(pA, pB, bfdr_alpha = 0.2)
  cat("   τ_gain =", round(trans_outer$tau_gain, 3),
      "| τ_loss =", round(trans_outer$tau_loss, 3),
      "| τ_cons =", round(trans_outer$tau_cons, 3), "\n")
  cat("   n_gain =", trans_outer$n_gain,
      "| n_loss =", trans_outer$n_loss,
      "| n_cons =", trans_outer$n_cons, "\n")
  
  # Phase inference for maintained genes at the inner threshold
  cat("→ Running phase_infer() at inner α =", a, "...\n")
  phase_inner <- phase_infer(
    phi_matrix1 = mcmc_age$COMBINED_younger$phi,
    phi_matrix2 = mcmc_age$COMBINED_older$phi,
    gain_loss_status = trans_outer$gain_loss_status,
    bfdr_alpha = a,
    shift = 4,
    P = 24
  )
  cat("   Maintained =", trans_outer$n_cons,
      "| Phase-conserved =", sum(phase_inner$flag_cons),
      "| Phase-shifted =", sum(phase_inner$flag_shift), "\n")
  
  gene_names <- names(phase_inner$peak1)
  gain_genes <- names(trans_outer$gain_loss_status)[trans_outer$gain_loss_status == "Gain"]
  loss_genes <- names(trans_outer$gain_loss_status)[trans_outer$gain_loss_status == "Loss"]
  maintained_genes <- names(trans_outer$gain_loss_status)[trans_outer$gain_loss_status == "Maintained"]
  
  phase_conserved <- gene_names[phase_inner$flag_cons]
  phase_shifted   <- gene_names[phase_inner$flag_shift]
  
  circadian_gain      <- intersect(circadian_genes, gain_genes)
  circadian_loss       <- intersect(circadian_genes, loss_genes)
  circadian_maintained <- intersect(circadian_genes, maintained_genes)
  circadian_conserved  <- intersect(circadian_genes, phase_conserved)
  circadian_shifted    <- intersect(circadian_genes, phase_shifted)
  
  cat("   Circadian genes:",
      "Gain =", length(circadian_gain),
      "| Loss =", length(circadian_loss),
      "| Maintained =", length(circadian_maintained),
      "| Phase-cons =", length(circadian_conserved),
      "| Phase-shift =", length(circadian_shifted), "\n")
  
  data.frame(
    BFDR      = a,
    Category  = c("Gain", "Loss", "Maintained",
                  "Phase-conserved", "Phase-shifted"),
    Total     = c(length(gain_genes),
                  length(loss_genes),
                  length(maintained_genes),
                  length(phase_conserved),
                  length(phase_shifted)),
    Circadian = c(length(circadian_gain),
                  length(circadian_loss),
                  length(circadian_maintained),
                  length(circadian_conserved),
                  length(circadian_shifted))
  )
})

# Undetermined = Maintained - Phase-conserved - Phase-shifted, per threshold
maintained_df <- subset(summary_table, Category == "Maintained")
conserved_df  <- subset(summary_table, Category == "Phase-conserved")
shifted_df    <- subset(summary_table, Category == "Phase-shifted")

undetermined_df <- maintained_df
undetermined_df$Category <- "Undetermined"
undetermined_df$Total     <- maintained_df$Total - conserved_df$Total - shifted_df$Total
undetermined_df$Circadian <- maintained_df$Circadian - conserved_df$Circadian - shifted_df$Circadian

summary_table <- rbind(summary_table, undetermined_df) %>%
  arrange(BFDR, match(Category, c("Gain", "Loss", "Maintained",
                                  "Phase-conserved", "Phase-shifted", "Undetermined")))

summary_table$Category <- dplyr::recode(summary_table$Category,
                                        "Gain"             = "R<sub>g</sub>",
                                        "Loss"             = "R<sub>l</sub>",
                                        "Maintained"       = "R<sub>c</sub>",
                                        "Phase-conserved"  = "P<sub>c</sub>",
                                        "Phase-shifted"    = "P<sub>s</sub>",
                                        "Undetermined"     = "P<sub>u</sub>"
)

summary_table %>%
  kable(
    "html",
    escape = FALSE,  # render the <sub> tags
    caption = "Dual-BFDR Summary of Rhythmicity and Phase Dynamics (Outer α=0.2)",
    col.names = c("Inner BFDR α", "Category", "Total Genes", "Circadian Genes"),
    align = "c",
    row.names = FALSE
  ) %>%
  kable_styling(
    full_width = FALSE,
    bootstrap_options = c("striped", "hover", "condensed")
  ) %>%
  collapse_rows(columns = 1, valign = "middle")

cat("\n========= FINAL SUMMARY =========\n")
summary_table %>%
  kable(
    "html",
    caption = "Dual-BFDR Summary of Rhythmicity and Phase Dynamics (Outer α=0.2)",
    col.names = c("Inner BFDR α", "Category", "Total Genes", "Circadian Genes"),
    align = "c",
    row.names = FALSE
  ) %>%
  kable_styling(
    full_width = FALSE,
    bootstrap_options = c("striped", "hover", "condensed")
  ) %>%
  collapse_rows(columns = 1, valign = "middle")

###################

# One row per gene from the last phase_infer() call
phase_df <- data.frame(
  Gene = gene_names,
  peak1 = as.numeric(phase_inner$peak1),
  peak2 = as.numeric(phase_inner$peak2),
  deltaPhi = phase_inner$deltaPhi.Est,
  BFDR_cons = phase$BFDR_cons,
  prob_conserved = as.numeric(phase_inner$prob_conserved),
  flag_cons = as.logical(phase_inner$flag_cons),
  stringsAsFactors = FALSE,
  row.names = NULL
)

conserved_df <- phase_df %>%
  filter(flag_cons == TRUE)

print("Conserved genes from phase object:")
print(conserved_df$Gene)

# Cosinor scatter plots for the phase-conserved genes
plot_gene_cosinor(
  gene_list   = conserved_df$Gene,
  COMBINED_data = COMBINED_data,
  df = df,
  period = 24,
  alpha = 0.05,
  save_path = file.path(current_aging, "figure/cosinor/conserved")
)

##################### Clock plot of phase-conserved genes
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

print("Current columns in conserved_df:")
print(colnames(conserved_df))

# Join the posterior summaries if they are not already present
if (!all(c("Bayesian_Median_Y", "Bayesian_Median_O") %in% colnames(conserved_df))) {
  conserved_df <- conserved_df %>%
    left_join(full_output %>% 
                select(Gene, Posterior_Rho_Y, Posterior_Rho_O,
                       Bayesian_Median_Y, Bayesian_Median_O,
                       Bayesian_SD_Y, Bayesian_SD_O),
              by = "Gene")
}

# Peak phases in radians on the ZT scale
if (!all(c("rad_Y", "rad_O") %in% colnames(conserved_df))) {
  conserved_df <- conserved_df %>%
    mutate(
      ZT_Y = to_zt(Bayesian_Median_Y),
      ZT_O = to_zt(Bayesian_Median_O),
      rad_Y = 2*pi*ZT_Y/24,
      rad_O = 2*pi*ZT_O/24
    )
}

plot_df <- conserved_df %>%
  select(Gene, rad_Y, rad_O, Posterior_Rho_Y, Posterior_Rho_O,
         Bayesian_SD_Y, Bayesian_SD_O) %>%
  pivot_longer(
    cols = c(rad_Y, rad_O, Posterior_Rho_Y, Posterior_Rho_O,
             Bayesian_SD_Y, Bayesian_SD_O),
    names_to = c(".value","Group"),
    names_pattern = "(rad|Posterior_Rho|Bayesian_SD)_([YO])"
  ) %>%
  mutate(Group = ifelse(Group=="Y","Younger","Older"))

# Wrap angles to [0, 2π)
plot_df <- plot_df %>%
  mutate(
    rad = (rad + 2*pi) %% (2*pi)
  )

# Offset the two age groups radially, further apart when their ±1 SD intervals overlap
plot_df <- plot_df %>%
  mutate(
    ci_lower = rad - 2*pi*Bayesian_SD/24,
    ci_upper = rad + 2*pi*Bayesian_SD/24
  )

overlap_check <- plot_df %>%
  select(Gene, Group, rad, ci_lower, ci_upper) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(rad, ci_lower, ci_upper)
  ) %>%
  mutate(
    overlap = (ci_lower_Younger <= ci_upper_Older) & (ci_lower_Older <= ci_upper_Younger),
    radial_jitter_Y = ifelse(overlap, -0.3, -0.15),
    radial_jitter_O = ifelse(overlap, 0.3, 0.15)
  ) %>%
  select(Gene, overlap, radial_jitter_Y, radial_jitter_O)

plot_df <- plot_df %>%
  left_join(overlap_check, by = "Gene") %>%
  mutate(
    radial_jitter = ifelse(Group == "Younger", radial_jitter_Y, radial_jitter_O),
    gene_radius = match(Gene, unique(Gene))
  ) %>%
  mutate(
    gene_radius_adj = gene_radius + radial_jitter
  ) %>%
  select(-ci_lower, -ci_upper, -overlap, -radial_jitter_Y, -radial_jitter_O, -radial_jitter)

gene_list <- unique(plot_df$Gene)
n_genes <- length(gene_list)

# Gene palette, no yellow
gene_colors <- setNames(
  c("#E41A1C",  # Red - BHLHE41
    "#377EB8",  # Blue - BMAL1
    "#4DAF4A",  # Green - CIART
    "#984EA3",  # Purple - DBP
    "#FF7F00",  # Orange - KCNH4
    "#A6D854",  # Light green
    "#A65628",  # Brown - NFIL3
    "#F781BF",  # Pink - NR1D1
    "#999999",  # Gray - NR1D2
    "#66C2A5",  # Teal - PER3
    "#FC8D62",  # Salmon
    "#8DA0CB",  # Light blue
    "#E78AC3",  # Light pink
    "#A6D854",  # Light green
    "#8B4513")[1:n_genes],
  gene_list
)

# Both age groups share one radius per gene, so their intervals overlap
plot_df_overlap <- plot_df %>%
  group_by(Gene) %>%
  mutate(
    gene_radius_shared = mean(gene_radius_adj),

    rad_lower = rad - 2*pi*Bayesian_SD/24,
    rad_upper = rad + 2*pi*Bayesian_SD/24
  ) %>%
  ungroup()

connection_df_overlap <- connection_df %>%
  left_join(
    plot_df_overlap %>% 
      filter(Group == "Younger") %>% 
      select(Gene, gene_radius_shared),
    by = "Gene"
  )

# Shade ZT-6 to ZT6
shade_df <- tibble(
  x = c((3*pi/2 + 2*pi)/2, (0 + pi/2)/2),
  y = (n_genes + 1) / 2,
  width = c(2*pi - 3*pi/2, pi/2),
  height = n_genes + 1
)

p <- ggplot() +
geom_tile(
  data = shade_df,
  aes(x = x, y = y, width = width, height = height),
  fill = "gray92", alpha = 0.35
) +

# Younger-to-older connecting lines
geom_segment(
  data = connection_df_overlap,
  aes(x = rad_Younger, xend = rad_Older,
      y = gene_radius_shared, yend = gene_radius_shared,
      color = Gene),
  linewidth = 1.5, alpha = 0.3, linetype = "solid"
) +

# Interval bands, wide and faint then narrow and darker
geom_segment(
  data = plot_df_overlap,
  aes(x = rad_lower, xend = rad_upper,
      y = gene_radius_shared, yend = gene_radius_shared,
      color = Gene, linetype = Group),
  linewidth = 10,
  alpha = 0.15,
  lineend = "round"
) +

geom_segment(
  data = plot_df_overlap,
  aes(x = rad_lower, xend = rad_upper,
      y = gene_radius_shared, yend = gene_radius_shared,
      color = Gene, linetype = Group),
  linewidth = 5,
  alpha = 0.3,
  lineend = "round"
) +

# Point-estimate ticks
geom_segment(
  data = plot_df_overlap,
  aes(x = rad - 0.02, xend = rad + 0.02,
      y = gene_radius_shared, yend = gene_radius_shared,
      color = Gene),
  linewidth = 3,
  alpha = 0.9,
  lineend = "round"
) +

geom_point(
  data = plot_df_overlap %>%
    mutate(
      # Small radial offset so both groups' points stay visible
      y_offset = ifelse(Group == "Younger",
                        gene_radius_shared + 0.08,
                        gene_radius_shared - 0.08)
    ),
  aes(x = rad, y = y_offset,
      fill = Gene, shape = Group),
  color = "white",
  size = 5,
  stroke = 1.8,
  alpha = 0.85
) +

coord_polar(start = -pi/2, direction = 1) +
  scale_x_continuous(
    limits = c(0, 2*pi),
    breaks = c(seq(18, 22, by = 2), seq(0, 16, by = 2)) * 2*pi/24,
    labels = c("ZT-6", "ZT-4", "ZT-2", paste0("ZT", seq(0, 16, by = 2)))
  ) +
  scale_y_continuous(
    limits = c(0.2, n_genes + 0.8),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = gene_colors) +
  scale_color_manual(values = gene_colors) +
  scale_shape_manual(
    values = c("Younger" = 24, "Older" = 21)
  ) +
  scale_linetype_manual(
    values = c("Younger" = "solid", "Older" = "dashed")
  ) +
  
  labs(
    title = "Circadian Peak Timing of Phase-Conserved Genes",
    subtitle = "Posterior phase estimates (ZT) ± 95% CI",
    fill = "Gene",
    color = "Gene",
    shape = "Age Group",
    linetype = "Age Group"
  ) +
  
  theme_minimal(base_size = 13, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20, margin = margin(b = 5)),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30", margin = margin(b = 15)),
    axis.text.x = element_text(size = 12, face = "bold", color = "gray20", margin = margin(t = 8)),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.major.x = element_line(color = "gray88", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.4, "cm"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.key.size = unit(1.2, "lines"),
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.5),
    legend.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(25, 25, 25, 25)
  ) +
  guides(
    fill = guide_legend(
      order = 1, 
      override.aes = list(alpha = 1, size = 4, shape = 21, stroke = 1.5)
    ),
    color = "none",
    shape = guide_legend(
      order = 2,
      override.aes = list(size = 4, stroke = 1.5, fill = "gray50", color = "white")
    ),
    linetype = guide_legend(
      order = 3,
      override.aes = list(linewidth = 2, alpha = 0.8)
    )
  )

print(p)
save_path <- file.path(current_aging, "figure")
pdf_file <- file.path(save_path, "clock_plot_combined_2.pdf")
png_file <- file.path(save_path, "clock_plot_combined_@.png")

pdf_file <- file.path(save_path, "clock_plot_combined.pdf")
png_file <- file.path(save_path, "clock_plot_combined.png")

ggsave(pdf_file, p, width = 12, height = 10, dpi = 800)
ggsave(png_file, p, width = 12, height = 10, dpi = 800)

# High-contrast gene palette
gene_colors <- setNames(
  c("#C1272D",  # Deep Red - BHLHE41
    "#0071BC",  # Deep Blue - BMAL1
    "#008856",  # Deep Green - CIART
    "#6B2D8F",  # Deep Purple - DBP
    "#E8601C",  # Deep Orange - KCNH4
    "#7C9C3D",  # Olive Green - LRRC39
    "#8B4513",  # Saddle Brown - NFIL3
    "#DC136C",  # Deep Pink - NR1D1
    "#4A4A4A",  # Dark Gray - NR1D2
    "#00897B",  # Deep Teal - PER3
    "#D84315",  # Deep Salmon
    "#5E7BA3",  # Deep Blue-Gray
    "#C2185B",  # Deep Rose
    "#689F38",  # Deep Lime
    "#6D4C41")[1:n_genes],  # Deep Brown
  gene_list
)

plot_df_overlap <- plot_df %>%
  group_by(Gene) %>%
  mutate(
    gene_radius_shared = mean(gene_radius_adj),

    rad_lower = rad - 2*pi*Bayesian_SD/24,
    rad_upper = rad + 2*pi*Bayesian_SD/24
  ) %>%
  ungroup()

connection_df_overlap <- connection_df %>%
  left_join(
    plot_df_overlap %>%
      filter(Group == "Younger") %>%
      select(Gene, gene_radius_shared),
    by = "Gene"
  )

shade_df <- tibble(
  x = c((3*pi/2 + 2*pi)/2, (0 + pi/2)/2),
  y = (n_genes + 1) / 2,
  width = c(2*pi - 3*pi/2, pi/2),
  height = n_genes + 1
)

# Same clock plot with heavier shading, lines and markers
p <- ggplot() +
geom_tile(
  data = shade_df,
  aes(x = x, y = y, width = width, height = height),
  fill = "gray85", alpha = 0.5
) +

geom_segment(
  data = connection_df_overlap,
  aes(x = rad_Younger, xend = rad_Older,
      y = gene_radius_shared, yend = gene_radius_shared,
      color = Gene),
  linewidth = 2, alpha = 0.5, linetype = "solid"
) +

geom_segment(
  data = plot_df_overlap,
  aes(x = rad_lower, xend = rad_upper,
      y = gene_radius_shared, yend = gene_radius_shared,
      color = Gene, linetype = Group),
  linewidth = 12,
  alpha = 0.25,
  lineend = "round"
) +

geom_segment(
  data = plot_df_overlap,
  aes(x = rad_lower, xend = rad_upper,
      y = gene_radius_shared, yend = gene_radius_shared,
      color = Gene, linetype = Group),
  linewidth = 6,
  alpha = 0.5,
  lineend = "round"
) +

geom_segment(
  data = plot_df_overlap,
  aes(x = rad - 0.02, xend = rad + 0.02,
      y = gene_radius_shared, yend = gene_radius_shared,
      color = Gene),
  linewidth = 4,
  alpha = 1,
  lineend = "round"
) +

geom_point(
  data = plot_df_overlap %>%
    mutate(
      y_offset = ifelse(Group == "Younger", 
                        gene_radius_shared + 0.08, 
                        gene_radius_shared - 0.08)
    ),
  aes(x = rad, y = y_offset,
      fill = Gene, shape = Group),
  color = "black",
  size = 6,
  stroke = 2.2,
  alpha = 0.7
) +

coord_polar(start = -pi/2, direction = 1) +
  scale_x_continuous(
    limits = c(0, 2*pi),
    breaks = c(seq(18, 22, by = 2), seq(0, 16, by = 2)) * 2*pi/24,
    labels = c("ZT-6", "ZT-4", "ZT-2", paste0("ZT", seq(0, 16, by = 2)))
  ) +
  scale_y_continuous(
    limits = c(0.2, n_genes + 0.8),
    expand = c(0, 0)
  ) +
  scale_fill_manual(values = gene_colors) +
  scale_color_manual(values = gene_colors) +
  scale_shape_manual(
    values = c("Younger" = 24, "Older" = 21)
  ) +
  scale_linetype_manual(
    values = c("Younger" = "solid", "Older" = "dashed")
  ) +
  
  labs(
    title = "Circadian Peak Timing of Phase-Conserved Genes",
    subtitle = "Posterior phase estimates (ZT) ± 95% CI",
    fill = "Gene",
    color = "Gene",
    shape = "Age Group",
    linetype = "Age Group"
  ) +
  
  theme_minimal(base_size = 14, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 22, 
                              margin = margin(b = 8), color = "gray10"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray20",
                                 margin = margin(b = 20)),
    axis.text.x = element_text(size = 13, face = "bold", color = "gray10", 
                               margin = margin(t = 8)),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid.major.y = element_line(color = "gray75", linewidth = 0.6),
    panel.grid.major.x = element_line(color = "gray75", linewidth = 0.6),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(0.6, "cm"),
    legend.text = element_text(size = 12, color = "gray10"),
    legend.title = element_text(size = 13, face = "bold", color = "gray10",
                                margin = margin(b = 8)),
    legend.key.size = unit(1.4, "lines"),
    legend.background = element_rect(fill = "white", color = "gray50", linewidth = 0.7),
    legend.margin = margin(15, 15, 15, 15),
    legend.box.spacing = unit(1.2, "cm"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(30, 30, 30, 30)
  ) +
  guides(
    fill = guide_legend(
      order = 1, 
      override.aes = list(alpha = 0.9, size = 5, shape = 21, stroke = 2, color = "black")
    ),
    color = "none",
    shape = guide_legend(
      order = 2,
      override.aes = list(size = 5, stroke = 2, fill = "gray40", color = "black", alpha = 0.9)
    ),
    linetype = guide_legend(
      order = 3,
      override.aes = list(linewidth = 2.5, alpha = 1, color = "gray30")
    )
  )

print(p)

ggsave(
  filename = file.path(save_path, "circadian_phase_conservation_nature.pdf"),
  plot = p,
  width = 12,
  height = 10,
  dpi = 600,
  device = cairo_pdf
)

ggsave(
  filename = file.path(save_path, "circadian_phase_conservation_nature.tiff"),
  plot = p,
  width = 12,
  height = 10,
  dpi = 600,
  compression = "lzw"
)

#───────────────────────────────────────────────────────────────
# Peak concordance plot for maintained genes (ZT scale)
#───────────────────────────────────────────────────────────────
library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)

to_zt <- function(t_cos) ifelse(t_cos >= 18, t_cos - 24, t_cos)

full_output <- read_excel(
  file.path(current_aging, "results/brain_regions/all_genes_merged.xlsx"),
  sheet = 1
)

full_output <- full_output %>%
  mutate(Peak_Y_ZT = to_zt(Peak_Y),
         Peak_O_ZT = to_zt(Peak_O))

# Maintained (rhythmic in both) genes at outer BFDR 0.25
trans_outer <- transition_classify(pA, pB, bfdr_alpha = 0.25)
maintained_genes <- names(trans_outer$gain_loss_status[trans_outer$gain_loss_status == "Maintained"])

maintained_df <- full_output %>%
  filter(Gene %in% maintained_genes) %>%
  filter(!is.na(Peak_Y) & !is.na(Peak_O)) %>%
  mutate(Peak_Y_ZT = to_zt(Peak_Y),
         Peak_O_ZT = to_zt(Peak_O))

# Phase classification at inner BFDR 0.15 with a 3 h shift
phase_inner <- phase_infer(
  phi_matrix1 = mcmc_age$COMBINED_younger$phi,
  phi_matrix2 = mcmc_age$COMBINED_older$phi,
  gain_loss_status = trans_outer$gain_loss_status,
  bfdr_alpha = 0.15,
  shift = 3,
  P = 24
)

gene_names <- names(phase_inner$peak1)

phase_class <- rep("Undetermined", length(gene_names))
names(phase_class) <- gene_names
phase_class[phase_inner$flag_cons]  <- "Phase-conserved"
phase_class[phase_inner$flag_shift] <- "Phase-shifted"

# Genes without a phase call are Undetermined
maintained_df <- maintained_df %>%
  mutate(phase_class = phase_class[Gene]) %>%
  mutate(phase_class = ifelse(is.na(phase_class), "Undetermined", phase_class))

circadian_genes <- c(
  "BHLHE40","BHLHE41","BMAL1","ARNTL","BTRC","CLOCK","CREB1",
  "CRY1","CRY2","CSNK1D","CSNK1E","CUL1","DBP","FBXL3","FBXW11",
  "NFIL3","NPAS2","NR1D1","NR1D2","PER1","PER2","PER3",
  "RORA","RORB","RORC"
)
genes_to_label <- maintained_df %>% filter(Gene %in% circadian_genes)

phase_colors <- c(
  "Phase-conserved" = "#1B9E77",
  "Phase-shifted"   = "#D95F02",
  "Undetermined"    = "gray70"
)

# Fraction of rhythmic genes (gain, loss or maintained) that are phase-conserved
tau_c <- 0.25

rhythmic_genes_outer <- names(trans$gain_loss_status)[
  trans$gain_loss_status %in% c("Gain", "Loss", "Maintained")
]

phase_conserved_genes <- names(phase$flag_cons)[phase$flag_cons]

n_rhythmic_outer <- length(rhythmic_genes_outer)
n_phase_conserved <- sum(phase$flag_cons, na.rm = TRUE)
pct_phase_conserved_global <- 100 * n_phase_conserved / n_rhythmic_outer

cat("\n[Global Summary]\n")
cat("  Rhythmic (outer τ_c ≤ 0.20):", n_rhythmic_outer, "\n")
cat("  Phase-conserved (inner τ_p ≤ 0.05):", n_phase_conserved, "\n")
cat(sprintf("  ⇒ %.2f%% phase-conserved among all rhythmic genes\n", pct_phase_conserved_global))

# ±4 h peak concordance among maintained genes
calculate_peak_difference <- function(a, b) {
  d <- abs(a - b)
  ifelse(d > 12, 24 - d, d)
}

maintained_df <- maintained_df %>%
  mutate(peak_diff = calculate_peak_difference(Peak_Y_ZT, Peak_O_ZT),
         within_concordance = peak_diff <= 4)

n_total <- nrow(maintained_df)
n_within <- sum(maintained_df$within_concordance, na.rm = TRUE)
pct_within <- 100 * n_within / n_total

cat("\n[Within ±4 h summary]\n")
cat("  Maintained genes:", n_total, "\n")
cat("  Within ±4 h:", n_within, "\n")
cat(sprintf("  ⇒ %.1f%% within ±4 h interval\n", pct_within))

#───────────────────────────────────────────────────────────────
# Plot; BMAL1 is drawn at hand-set coordinates
#───────────────────────────────────────────────────────────────
BMAL1_fix <- maintained_df %>%
  filter(Gene == "BMAL1") %>%
  mutate(
    Peak_Y_ZT_plot = 17.6,
    Peak_O_ZT_plot = 18.1
  )

others_df <- maintained_df %>%
  filter(Gene != "BMAL1") %>%
  mutate(
    Peak_Y_ZT_plot = Peak_Y_ZT,
    Peak_O_ZT_plot = Peak_O_ZT
  )

plot_df <- bind_rows(others_df, BMAL1_fix)

library(ggplot2)
library(ggrepel)
library(dplyr)

Rc_df <- plot_df %>% filter(phase_class == "Phase-conserved")

n_Rc <- nrow(plot_df)               # size of the rhythmically conserved set R_c
pct_Rc <- round(100 * nrow(Rc_df) / nrow(plot_df), 1)

subtitle_text <- bquote(
  "Rhythmically Conserved Set " ~ R[c] ~
    "(" * n == .(n_Rc) * ", " * .(pct_Rc) * "% within " * "\u00B1" * "4 h interval)"
)

p <- ggplot(plot_df, aes(
  x = Peak_Y_ZT_plot,
  y = Peak_O_ZT_plot,
  color = phase_class
)) +
  # identity and ±4 h boundaries
  geom_abline(intercept = 0, slope = 1,
              color = "black", linetype = "dashed", linewidth = 1.0) +
  geom_abline(intercept = 4, slope = 1,
              color = "darkgreen", linetype = "dotted", linewidth = 1.1) +
  geom_abline(intercept = -4, slope = 1,
              color = "darkgreen", linetype = "dotted", linewidth = 1.1) +
  
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = phase_colors) +
  
  geom_text_repel(
    data = plot_df %>% semi_join(genes_to_label, by = "Gene"),
    aes(label = Gene),
    color = "black",
    fontface = "bold.italic",
    size = 4,
    segment.color = "gray50",
    box.padding = 1.2,
    point.padding = 0.7,
    max.overlaps = Inf
  ) +

  labs(
    title = "Circadian Peak Concordance Among Rhythmically Conserved Genes",
    subtitle = subtitle_text,
    x = "Peak Hour – Younger Group (ZT)",
    y = "Peak Hour – Older Group (ZT)",
    color = "Phase class"
  ) +
  
  scale_x_continuous(
    limits = c(-6, 18.5),
    breaks = seq(-6, 18, 6),
    labels = sprintf("ZT%+d", seq(-6, 18, 6))
  ) +
  scale_y_continuous(
    limits = c(-6, 18.5),
    breaks = seq(-6, 18, 6),
    labels = sprintf("ZT%+d", seq(-6, 18, 6))
  ) +

  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 12),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(color = "gray70", fill = "white"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.margin = margin(t = 5, b = 5),
    
    plot.margin = margin(15, 15, 15, 15)
  )

p

#───────────────────────────────────────────────────────────────
# Save
#───────────────────────────────────────────────────────────────
save_dir <- file.path(current_aging, "figure/cosinor")
dir.create(save_dir, showWarnings = FALSE)
ggsave(file.path(save_dir, "Peak_Concordance_Plot_Maintained_conserved.pdf"),
       plot = p, width = 9, height = 8)
p

# ============================================================================
# Aim 1C (or 3A): Pathway Analysis
# ============================================================================
source(file.path(current_wd, "Kyle/Circadian-analysis-main/R/src/Thien/pathwaySelect.R"))

################################################################################
# Workflow: Aim 3A → B1, B2, B3
################################################################################

# Parameters shared by every stage
dataset_names <- c("Younger", "Older")
qvalue_cut <- 0.05      # FDR for pathway significance
nperm <- 10000
pathway_size_min <- 5
pathway_size_max <- 300

################################################################################
# STAGE 1: Identify Active Pathways (Union Test)
################################################################################

cat("\n=== STAGE 1: Union Test (Active Pathways) ===\n")

result_union <- pathSelect(
  mcmc.merge.list = mcmc_age,
  pathway.list = kegg.pathway.list_hsa,
  dataset.names = c("younger", "older"),
  ranking.method = "union",
  score_type = "pos",
  qvalue.cut = 0.25,
  pathwaysize.lower.cut = pathway_size_min,
  pathwaysize.upper.cut = pathway_size_max,
  nperm = nperm,
  nproc = 1
)

active_pathways <- result_union$results %>%
  filter(Significant == TRUE) %>%
  pull(pathway)

cat("Active pathways:", length(active_pathways), "\n")

################################################################################
# STAGE 2: Test Transitions on Active Pathways
################################################################################

cat("\n=== STAGE 2: Transition Tests (Gain, Loss, Conservation) ===\n")

active_pathway_list <- kegg.pathway.list_hsa[match(active_pathways, names(kegg.pathway.list_hsa))]

result_gain <- pathSelect(
  mcmc.merge.list = mcmc_age,
  pathway.list = active_pathway_list,
  dataset.names = c("younger", "older"),
  ranking.method = "gain",
  score_type = "pos",
  qvalue.cut = 0.05,
  pathwaysize.lower.cut = pathway_size_min,
  pathwaysize.upper.cut = pathway_size_max,
  nperm = nperm,
  nproc = 1
)

result_loss <- pathSelect(
  mcmc.merge.list = mcmc_age,
  pathway.list = active_pathway_list,
  dataset.names = c("younger", "older"),
  ranking.method = "loss",
  score_type = "pos",
  qvalue.cut = 0.05,
  pathwaysize.lower.cut = pathway_size_min,
  pathwaysize.upper.cut = pathway_size_max,
  nperm = nperm,
  nproc = 1
)

result_cons <- pathSelect(
  mcmc.merge.list = mcmc_age,
  pathway.list = active_pathway_list,
  dataset.names = c("younger", "older"),
  ranking.method = "conserved",
  score_type = "pos",
  qvalue.cut = 0.05,
  pathwaysize.lower.cut = pathway_size_min,
  pathwaysize.upper.cut = pathway_size_max,
  nperm = nperm,
  nproc = 1
)

cat("\nEnrichment Results:\n")
cat("  Gain-enriched:        ", sum(result_gain$results$Significant), "\n")
cat("  Loss-enriched:        ", sum(result_loss$results$Significant), "\n")
cat("  Conservation-enriched:", sum(result_cons$results$Significant), "\n")

################################################################################
# STAGE 3: Filter to Significant Pathways
################################################################################

cat("\n=== STAGE 3: Filtering Significant Pathways ===\n")

# Joined by pathway name
significance_table <- result_gain$results %>%
  select(pathway, gain_sig = Significant) %>%
  full_join(
    result_loss$results %>% select(pathway, loss_sig = Significant),
    by = "pathway"
  ) %>%
  full_join(
    result_cons$results %>% select(pathway, cons_sig = Significant),
    by = "pathway"
  ) %>%
  filter(pathway %in% active_pathways) %>%
  mutate(
    gain_sig = replace_na(gain_sig, FALSE),
    loss_sig = replace_na(loss_sig, FALSE),
    cons_sig = replace_na(cons_sig, FALSE)
  ) %>%
  arrange(pathway)

significant_pathways <- significance_table %>%
  filter(gain_sig | loss_sig | cons_sig) %>%
  pull(pathway)

cat("Pathways significant for at least one transition:", length(significant_pathways), "\n")

filtered_pathway_list <- kegg.pathway.list_hsa[match(significant_pathways, names(kegg.pathway.list_hsa))]

################################################################################
# STAGE 4: Run multi_conservation on FILTERED Pathways
################################################################################

cat("\n=== STAGE 4: Descriptive Metrics (multi_conservation) ===\n")
cat("Running on", length(significant_pathways), "filtered pathways...\n")

result_multiconservation <- multi_conservation(
  mcmc.merge.list = mcmc_age,
  dataset.names = c("younger", "older"),
  select.pathway.list = filtered_pathway_list,
  n_perm = 10000,
  n_boot = 1000,
  output.dir = file.path(output.dir, "multiconservation_filtered"),
  use_cpp = TRUE
)

cat("Descriptive metrics calculated for", nrow(result_multiconservation), "pathways\n")

################################################################################
# STAGE 5: Combine Results
################################################################################

cat("\n=== STAGE 5: Combining Results ===\n")

pvalues_table <- result_gain$results %>%
  select(pathway, gain_padj = padj, gain_NES = NES) %>%
  left_join(result_loss$results %>% select(pathway, loss_padj = padj, loss_NES = NES), by = "pathway") %>%
  left_join(result_cons$results %>% select(pathway, cons_padj = padj, cons_NES = NES), by = "pathway")

# Descriptive metrics from pathSelect
descriptive_pathselect <- result_cons$results %>%
  select(
    pathway, size,
    Conserved_Index, Gain_Index, Loss_Index,
    Expected_N_Conserved, Expected_N_Gain, Expected_N_Loss,
    Gain_Loss_Ratio_Arithmetic,
    Top_Conserved_Genes, Top_Gain_Genes, Top_Loss_Genes
  )

# Descriptive metrics from multi_conservation
descriptive_multiconservation <- result_multiconservation %>%
  select(
    pathway = Pathway,
    Raw_Spearman = younger_vs_older_Congruence,
    Adjusted_Congruence = younger_vs_older_AdjustedCongruence,
    MC_GainIndex = younger_vs_older_GainIndex,
    MC_LossIndex = younger_vs_older_LossIndex,
    MC_GainLossRatio = younger_vs_older_GainLossRatio
  )

bfdr_alpha <- 0.05
final_results <- significance_table %>%
  filter(pathway %in% significant_pathways) %>%
  left_join(pvalues_table, by = "pathway") %>%
  left_join(descriptive_pathselect, by = "pathway") %>%
  left_join(descriptive_multiconservation, by = "pathway") %>%
  mutate(
    gain_sig = ifelse(!is.na(gain_padj), gain_padj < bfdr_alpha, gain_sig),
    loss_sig = ifelse(!is.na(loss_padj), loss_padj < bfdr_alpha, loss_sig),
    cons_sig = ifelse(!is.na(cons_padj), cons_padj < bfdr_alpha, cons_sig),
    enrichment_pattern = case_when(
      gain_sig & !loss_sig & !cons_sig ~ "Gain-enriched",
      !gain_sig & loss_sig & !cons_sig ~ "Loss-enriched",
      !gain_sig & !loss_sig & cons_sig ~ "Conservation-enriched",
      gain_sig & loss_sig & !cons_sig ~ "Mixed: Gain+Loss",
      gain_sig & !loss_sig & cons_sig ~ "Mixed: Gain+Conserved",
      !gain_sig & loss_sig & cons_sig ~ "Mixed: Loss+Conserved",
      gain_sig & loss_sig & cons_sig ~ "Mixed: All",
      TRUE ~ "Other"
    ),
    log_GLR_pathSelect = log(Gain_Loss_Ratio_Arithmetic),
    log_GLR_multiconservation = log(MC_GainLossRatio)
  ) %>%
  arrange(pathway)

cat("\nFinal results table:", nrow(final_results), "pathways\n")

################################################################################
# STAGE 7: Export Results
################################################################################

cat("\n=== STAGE 7: Exporting Results ===\n")

wb <- createWorkbook()

addWorksheet(wb, "All_Significant_Pathways")
writeData(wb, "All_Significant_Pathways", final_results)

# Subsets
gain_results <- final_results %>%
  filter(gain_sig) %>%
  mutate(
    Expected_Gain_to_Size = Expected_N_Gain / size)  %>%
  select(pathway, size, enrichment_pattern, Expected_Gain_to_Size,
         gain_padj, loss_padj, cons_padj,
         Gain_Index, Loss_Index, Conserved_Index,
         Adjusted_Congruence,
         Gain_Loss_Ratio_Arithmetic, MC_GainLossRatio,
         Expected_N_Gain, Expected_N_Loss, Expected_N_Conserved,
         Top_Gain_Genes) %>%
  arrange(gain_padj)
loss_results <- final_results %>%
  filter(loss_sig) %>%
  mutate(
    Expected_Loss_to_Size = Expected_N_Loss / size)  %>%
  select(pathway, size, enrichment_pattern, Expected_Loss_to_Size,
         gain_padj, loss_padj, cons_padj,
         Gain_Index, Loss_Index, Conserved_Index,
         Adjusted_Congruence,
         Gain_Loss_Ratio_Arithmetic, MC_GainLossRatio,
         Expected_N_Gain, Expected_N_Loss, Expected_N_Conserved,
         Top_Loss_Genes) %>%
  arrange(loss_padj)
cons_results <- final_results %>%
  filter(cons_sig) %>%
  mutate(
    Expected_Conserved_to_Size= Expected_N_Conserved / size)  %>%
  select(pathway, size, enrichment_pattern, Expected_Conserved_to_Size,
         gain_padj, loss_padj, cons_padj,
         Conserved_Index, Adjusted_Congruence, Raw_Spearman,
         Gain_Loss_Ratio_Arithmetic, MC_GainLossRatio,
         Gain_Index, Loss_Index,
         Expected_N_Gain, Expected_N_Loss, Expected_N_Conserved,
         Top_Conserved_Genes) %>%
  arrange(cons_padj)

addWorksheet(wb, "Gain_Enriched");  writeData(wb, "Gain_Enriched", gain_results)
addWorksheet(wb, "Loss_Enriched");  writeData(wb, "Loss_Enriched", loss_results)
addWorksheet(wb, "Conservation_Enriched"); writeData(wb, "Conservation_Enriched", cons_results)

# Summary sheet
summary_df <- final_results %>%
  select(pathway, gain_sig, loss_sig, cons_sig,
         gain_padj, loss_padj, cons_padj,
         Adjusted_Congruence, Conserved_Index,
         MC_GainIndex, MC_LossIndex, MC_GainLossRatio, enrichment_pattern)
addWorksheet(wb, "Summary_for_Visualization")
writeData(wb, "Summary_for_Visualization", summary_df)

header_style <- createStyle(textDecoration = "bold", fgFill = "#E6E6FA")
for (sheet in names(wb)) {
  addStyle(wb, sheet = sheet, style = header_style,
           rows = 1, cols = 1:100, gridExpand = TRUE)
  setColWidths(wb, sheet = sheet, cols = 1:100, widths = "auto")
}

################################################################################
# STAGE 7B: Expected Counts and Ratios Summary Sheet
################################################################################

cat("\n=== STAGE 7B: Adding Expected Count Summary ===\n")

expected_summary <- final_results %>%
  mutate(
    Expected_Union = Expected_N_Gain + Expected_N_Loss + Expected_N_Conserved,

    Ratio_Union_to_Size = Expected_Union / size,
    Ratio_Gain_to_Size = Expected_N_Gain / size,
    Ratio_Loss_to_Size = Expected_N_Loss / size,
    Ratin_Conserved_to_Size = Expected_N_Conserved / size,

    # Check column
    Sum_of_Indices = Gain_Index + Loss_Index + Conserved_Index
  ) %>%
  select(
    Pathway = pathway,
    Pathway_Size = size,

    Expected_Union,
    Expected_Union_to_Size = Ratio_Union_to_Size,
    Expected_Gain_to_Size = Ratio_Gain_to_Size,
    Expected_Loss_to_Size = Ratio_Loss_to_Size,
    Expected_Conserved_to_Size= Ratin_Conserved_to_Size,

    Adjusted_Congruence,

    Expected_N_Gain,
    Gain_Index,
    Gain_padj = gain_padj,

    Expected_N_Loss,
    Loss_Index,
    Loss_padj = loss_padj,

    Expected_N_Conserved,
    Conserved_Index,
    Conserved_padj = cons_padj,

    Gain_Loss_Ratio = Gain_Loss_Ratio_Arithmetic,
    Log_Gain_Loss_Ratio = log_GLR_pathSelect,

    Raw_Spearman,
    MC_GainLossRatio,
    enrichment_pattern,

    Gain_NES = gain_NES,
    Loss_NES = loss_NES,
    Cons_NES = cons_NES,

    Top_Gain_Genes,
    Top_Loss_Genes,
    Top_Conserved_Genes,

    Sum_of_Indices
  ) %>%
  arrange(desc(Expected_Union_to_Size))

cat("\nVerification of mathematical relationships:\n")
cat("All indices sum to 1:", all(abs(expected_summary$Sum_of_Indices - 1) < 1e-6, na.rm = TRUE), "\n")
cat("Range of Expected_Union_to_Size:", 
    range(expected_summary$Expected_Union_to_Size, na.rm = TRUE), "\n")

addWorksheet(wb, "Expected_Counts_Summary")
writeData(wb, "Expected_Counts_Summary", expected_summary)

addStyle(
  wb,
  sheet = "Expected_Counts_Summary",
  style = createStyle(textDecoration = "bold", fgFill = "#FFF2CC"),
  rows = 1, cols = 1:ncol(expected_summary), gridExpand = TRUE
)

setColWidths(
  wb,
  sheet = "Expected_Counts_Summary",
  cols = 1:ncol(expected_summary),
  widths = "auto"
)

# Highlight adjusted p-values below 0.05
conditionalFormatting(
  wb, 
  sheet = "Expected_Counts_Summary",
  cols = which(names(expected_summary) %in% c("Gain_padj", "Loss_padj", "Conserved_padj")),
  rows = 2:(nrow(expected_summary) + 1),
  rule = "< 0.05",
  style = createStyle(bgFill = "#90EE90")
)

cat("Expected counts summary sheet added with organized columns.\n")
cat("Expected counts summary sheet added with organized columns.\n")
saveWorkbook(wb, file.path(current_aging, "pathway_results_filtered.xlsx"), overwrite = TRUE)

# ============================================================================
# Heatmap
# ============================================================================

# Transition and phase calls shared by every heatmap below
trans_outer <- transition_classify(pA, pB, bfdr_alpha = 0.25)
phase_inner <- phase_infer(
  phi_matrix1 = mcmc_age$COMBINED_younger$phi,
  phi_matrix2 = mcmc_age$COMBINED_older$phi,
  gain_loss_status = trans_outer$gain_loss_status,
  bfdr_alpha = 0.10,
  shift = 4,
  P = 24
)

phase_inner_01 <- phase_infer(
  phi_matrix1 = mcmc_age$COMBINED_younger$phi,
  phi_matrix2 = mcmc_age$COMBINED_older$phi,
  gain_loss_status = trans_outer$gain_loss_status,
  bfdr_alpha = 0.20,
  shift = 4,
  P = 24
)

source(file.path(current_aging, "code/heatmap.R"))

###############################################
# PATHWAY LIST TO PLOT
###############################################

gain_pathways <- c(
  "KEGG Hematopoietic cell lineage",
  "KEGG Fluid shear stress and atherosclerosis",
  "KEGG NF-kappa B signaling pathway",
  "KEGG TNF signaling pathway",
  "KEGG Th17 cell differentiation",
  "KEGG Cytokine-cytokine receptor interaction",
  "KEGG Osteoclast differentiation",
  "KEGG Cell adhesion molecules",
  "KEGG IL-17 signaling pathway",
  "KEGG Epstein-Barr virus infection",
  "KEGG Adipocytokine signaling pathway",
  "KEGG p53 signaling pathway",
  "KEGG Alcoholic liver disease",
  "KEGG Phagosome",
  "KEGG Neutrophil extracellular trap formation",
  "KEGG Herpes simplex virus 1 infection"
)

loss_pathways <- c(
  "KEGG MAPK signaling pathway",
  "KEGG Glutamatergic synapse",
  "KEGG Oxytocin signaling pathway",
  "KEGG Osteoclast differentiation",
  "KEGG Leishmaniasis",
  "KEGG Morphine addiction",
  "KEGG Aldosterone synthesis and secretion",
  "KEGG Circadian rhythm",
  "KEGG Circadian entrainment",
  "KEGG FoxO signaling pathway",
  "KEGG Rap1 signaling pathway",
  "KEGG Glioma",
  "KEGG Transcriptional misregulation in cancer",
  "KEGG Th17 cell differentiation",
  "KEGG Cell adhesion molecules",
  "KEGG Phagosome"
)

conserved_pathways <- c(
  "KEGG Circadian rhythm",
  "KEGG Staphylococcus aureus infection",
  "KEGG Transcriptional misregulation in cancer",
  "KEGG Cell adhesion molecules",
  "KEGG Circadian entrainment",
  "KEGG Leishmaniasis",
  "KEGG Acute myeloid leukemia",
  "KEGG Neutrophil extracellular trap formation",
  "KEGG Tuberculosis",
  "KEGG Phagosome",
  "KEGG Viral protein interaction with cytokine and cytokine receptor",
  "KEGG Systemic lupus erythematosus",
  "KEGG Cytokine-cytokine receptor interaction",
  "KEGG Osteoclast differentiation"
)

###############################################
# SET UP OUTPUT DIRECTORIES
###############################################

base_path <- file.path(current_aging, "results/brain_regions/output_final/heatmap_new")

dirs <- c("GAIN", "LOSS", "CONSERVATION", "OTHER")
for (d in dirs) dir.create(file.path(base_path, d), showWarnings = FALSE, recursive = TRUE)

###############################################################
# PATHWAYS TO PLOT
###############################################################

all_pathways_to_plot <- unique(c(
  gain_pathways,
  loss_pathways,
  conserved_pathways
))

length(all_pathways_to_plot)
print(all_pathways_to_plot)

for (pathway in all_pathways_to_plot) {
  
  # A pathway in several lists is filed under the first
  category <- if (pathway %in% gain_pathways) {
    "GAIN"
  } else if (pathway %in% loss_pathways) {
    "LOSS"
  } else if (pathway %in% conserved_pathways) {
    "CONSERVATION"
  } else {
    next
  }

  cat("\nProcessing:", pathway, "| Category:", category, "\n")

  if (!pathway %in% names(kegg.pathway.list_hsa)) {
    cat("   ✗ Pathway not found — skip\n")
    next
  }
  
  pathway_genes <- kegg.pathway.list_hsa[[pathway]]
  overlap <- intersect(rownames(mcmc_age$COMBINED_younger$rho), pathway_genes)
  
  if (length(overlap) == 0) {
    cat("   ✗ No overlapping genes — skip\n")
    next
  }
  
  # Output folder: CATEGORY/PATHWAY
  safe_name <- gsub("[^A-Za-z0-9_]", "_", pathway)
  out_dir <- file.path(base_path, category, safe_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  tryCatch({
    
    while (dev.cur() > 1) dev.off()
    
    plot_heatmap(
      data1 = mcmc_age$COMBINED_younger,
      data2 = mcmc_age$COMBINED_older,
      pathway_genes = pathway_genes,
      pathway_name = pathway,
      phase_results = phase_inner,
      transition_results = trans_outer,
      group_names = c("Younger", "Older"),
      save_path = out_dir
    )
    
    cat("   ✓ Saved to:", out_dir, "\n")
    
  }, error = function(e) {
    cat("   ✗ ERROR:", e$message, "\n")
  })
}

cat("\n=== COMPLETE: ALL GAIN / LOSS / CONSERVATION PATHWAYS PLOTTED ===\n")

# ============================================
# KEGG MODULE ANALYSIS
# ============================================
# Module search within KEGG pathways, younger against older

library(biomaRt)
library(KEGGgraph)
library(igraph)
library(ggplot2)
library(pathview)
library(parallel)

ensemble <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# ============================================
# CONFIGURATION
# ============================================

pathway_ids <- c(
  # core circadian / signaling that lost rhythm
  "KEGG Circadian entrainment"                  = "hsa04713",
  "KEGG Circadian rhythm"                       = "hsa04710",
  "KEGG MAPK signaling pathway"                 = "hsa04010",
  "KEGG FoxO signaling pathway"                 = "hsa04068",
  "KEGG Oxytocin signaling pathway"             = "hsa04921",
  
  # immune–inflammatory that gained rhythm
  "KEGG NF-kappa B signaling pathway"           = "hsa04064",
  "KEGG Cytokine-cytokine receptor interaction" = "hsa04060",
  "KEGG Th1 and Th2 cell differentiation"       = "hsa04658",
  "KEGG Th17 cell differentiation"              = "hsa04659",
  "KEGG Phagosome"                              = "hsa04145",
  
  # vascular / metabolic gain
  "KEGG Adipocytokine signaling pathway"        = "hsa04920",
  "KEGG Fluid shear stress and atherosclerosis" = "hsa05418"
)

# Run on the circadian rhythm pathway only
pathway_ids <- c(
  "KEGG Circadian rhythm"                        = "hsa04710"
)
pairs <- list(
  young_old = c("younger", "older")
)

# One task per gene type, pathway and contrast
tasks <- expand.grid(
  gene_type = c("all", "concordant", "discordant"), 
  pathway   = names(pathway_ids),
  contrast  = names(pairs),
  stringsAsFactors = FALSE
)
tasks$pw_id <- sub("hsa", "", pathway_ids[tasks$pathway])

cat("Total analysis tasks:", nrow(tasks), "\n")
print(tasks)

# ============================================
# CREATE OUTPUT DIRECTORIES
# ============================================

base_path <- file.path(current_aging, "results/brain_regions/output_final/module")
if(!dir.exists(base_path)) dir.create(base_path, recursive = TRUE)

out_base <- file.path(base_path, "out")
plot_base <- file.path(base_path, "module_plot")

if (!dir.exists(plot_base)) {
  dir.create(plot_base, recursive = TRUE)
}
if (!dir.exists(out_base)) {
  dir.create(out_base, recursive = TRUE)
}

cat("Output directory:", out_base, "\n")
cat("Plot directory:", plot_base, "\n")

# ============================================
# STEP 1: DETECT RHYTHMICITY
# ============================================
cat("\n=== STEP 1: Running detect_rhy ===\n")

# mcmc_age needs COMBINED_younger and COMBINED_older, each with $rho
if (!exists("mcmc_age")) {
  stop("mcmc_age object not found! Please load your MCMC results first.")
}

cat("Rhythmicity Summary:\n")
cat("  Total genes:", rhy_results$n_total, "\n")
cat("  Rhythmic in younger:", rhy_results$n_rhythmic_A, "\n")
cat("  Rhythmic in older:", rhy_results$n_rhythmic_B, "\n")
cat("  Threshold younger:", round(rhy_results$threshold_A, 3), "\n")
cat("  Threshold older:", round(rhy_results$threshold_B, 3), "\n")

# ============================================
# STEP 2: PREPARE RHYTHMICITY STATUS DATA
# ============================================
cat("\n=== STEP 2: Preparing rhythmicity status ===\n")

gene_names <- rownames(mcmc_age$COMBINED_younger$rho)
cat("Total genes in dataset:", length(gene_names), "\n")
cat("First 5 gene names:", paste(head(gene_names, 5), collapse=", "), "\n")

younger_status <- data.frame(
  rhythmic = rhy_results$rhythmic_A_logical,
  row.names = gene_names,
  stringsAsFactors = FALSE
)

older_status <- data.frame(
  rhythmic = rhy_results$rhythmic_B_logical,
  row.names = gene_names,
  stringsAsFactors = FALSE
)

rhy_status_list <- list(
  younger = younger_status,
  older = older_status
)

cat("Rhythmicity status prepared:\n")
cat("  younger: ", nrow(younger_status), "genes,", sum(younger_status$rhythmic), "rhythmic\n")
cat("  older: ", nrow(older_status), "genes,", sum(older_status$rhythmic), "rhythmic\n")

# ============================================
# STEP 3: PREPARE MCMC DATA
# ============================================
cat("\n=== STEP 3: Preparing MCMC data ===\n")

data_combined <- list(
  younger = mcmc_age$COMBINED_younger,
  older   = mcmc_age$COMBINED_older
)

# ============================================
# STEP 4: RUN KEGG MODULE ANALYSIS
# ============================================
cat("\n=== STEP 4: Running KEGG_module analysis ===\n")

all_results <- list()

for (i in seq_len(nrow(tasks))) {
  
  gt <- tasks$gene_type[i]
  pw <- tasks$pw_id[i]
  pw_name <- tasks$pathway[i]
  sel <- pairs[[tasks$contrast[i]]]
  
  cat("\n========================================\n")
  cat("Processing [", i, "/", nrow(tasks), "]: ", pw_name, "\n", sep="")
  cat("Pathway ID: ", pw, "\n", sep="")
  cat("Gene type: ", gt, "\n", sep="")
  cat("Contrast: ", paste(sel, collapse=" vs "), "\n", sep="")
  cat("========================================\n")

  out_dir  <- file.path(out_base, gt, tasks$contrast[i], tasks$pw_id[i])
  plot_dir <- file.path(plot_base, gt, tasks$contrast[i], tasks$pw_id[i])
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  result <- tryCatch({
    KEGG_module(
      mcmc.merge.list = data_combined[sel],
      dataset.names   = sel,
      KEGGspecies     = "hsa",
      KEGGpathwayID   = pw,
      data.pair       = sel,
      gene_type       = gt,
      rhy_results_list = rhy_status_list,
      minM            = 2,
      maxM            = 10,  # maximum module size
      B               = 1000,  # permutations
      cores           = 1,
      search_method   = "SA",  # simulated annealing
      Elbow_plot      = TRUE,
      filePath        = out_dir,
      ensemble        = ensemble,
      seed            = 12345
    )
  }, error = function(e) {
    cat(sprintf("  ✗ Error: %s\n", e$message))
    cat("\nFull error details:\n")
    print(e)
    return(NULL)
  })

  if (!is.null(result)) {
    cat("\n  ✓ KEGG_module completed successfully!\n")
    cat("  Best module size:", result$bestSize, "\n")
    cat("  Number of module sizes tested:", length(result$minG.ls), "\n")

    result_file <- file.path(out_dir, "kegg_module_result.rds")
    saveRDS(result, result_file)
    cat("  ✓ Result saved to:", result_file, "\n")

    task_id <- paste(gt, pw_name, paste(sel, collapse="_vs_"), sep="_")
    all_results[[task_id]] <- result

    if (exists("KEGG_module_topology_plot")) {
      cat("\n  Attempting to create topology plot...\n")
      tryCatch({
        KEGG_module_topology_plot(
          res_KEGG_module = result,
          which_to_draw   = "all",
          filePath        = plot_dir,
          ensemble        = ensembl,
          mcmc_data       = data_combined,
          layout_type =  "all",
          create_pathview = TRUE
        )
        cat("Topology plot completed!\n")
      }, error = function(e) {
        cat(sprintf(" Topology plot failed: %s\n", e$message))
      })
    }
    
  } else {
    cat("  KEGG_module returned NULL - analysis failed\n")
  }

  gc()
  Sys.sleep(0.5)
}

