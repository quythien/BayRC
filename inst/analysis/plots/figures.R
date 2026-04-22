# Load package
########################################
COMBINED <- readRDS("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/Paper/Congruence/PNAS_aging/data/combined_data.rds")

prepare_combined_data <- function(combined_data) {
  cat("Debugging combined data preparation...\n")
  
  # Check input structure
  cat("Combined data structure:\n")
  cat("- Expression matrix dimensions:", dim(combined_data$expr), "\n")
  cat("- Phenotype data dimensions:", dim(combined_data$pheno), "\n")
  cat("- Expression column names sample:", head(colnames(combined_data$expr), 3), "\n")
  cat("- Phenotype column names:", colnames(combined_data$pheno), "\n")
  
  # CRITICAL FIX: Align phenotype data to expression matrix order
  expr_sample_names <- colnames(combined_data$expr)
  
  if(!"sample_name" %in% colnames(combined_data$pheno)) {
    stop("sample_name column not found in phenotype data")
  }
  
  cat("Aligning phenotype data to expression matrix order...\n")
  # Match phenotype rows to expression column order
  pheno_order <- match(expr_sample_names, combined_data$pheno$sample_name)
  
  if(any(is.na(pheno_order))) {
    missing_samples <- sum(is.na(pheno_order))
    cat("WARNING:", missing_samples, "expression samples not found in phenotype data\n")
    # Remove missing samples from expression matrix
    valid_expr_samples <- !is.na(pheno_order)
    expr_sample_names <- expr_sample_names[valid_expr_samples]
    combined_data$expr <- combined_data$expr[, valid_expr_samples]
    pheno_order <- pheno_order[valid_expr_samples]
  }
  
  # Reorder phenotype data to match expression matrix
  pheno_data <- combined_data$pheno[pheno_order, ]
  
  cat("After alignment:\n")
  cat("- Expression samples:", ncol(combined_data$expr), "\n")
  cat("- Phenotype rows:", nrow(pheno_data), "\n")
  cat("- Sample alignment check:", identical(colnames(combined_data$expr), pheno_data$sample_name), "\n")
  
  # Handle the merge result columns (TOD.x, TOD.y) and age group NAs
  cat("Age group distribution before filtering:\n")
  print(table(pheno_data$age_group, useNA = "ifany"))
  
  # Use TOD.x (from combined_sample_info_clean) as primary, fallback to TOD.y if needed
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
  
  # Use AgeGroup (from BA11$pheno) as primary, fallback to age_group if needed
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
  
  # Remove samples with missing age group or TOD data
  # Now we can safely use positional indexing because data is aligned
  complete_samples <- !is.na(pheno_data$age_group_final) & 
    !is.na(pheno_data$tod) &
    pheno_data$age_group_final %in% c("younger", "older")
  
  cat("Complete samples:", sum(complete_samples), "out of", length(complete_samples), "\n")
  cat("Samples removed due to missing/NA age group:", sum(is.na(pheno_data$age_group_final)), "\n")
  cat("Samples removed due to missing TOD:", sum(is.na(pheno_data$tod)), "\n")
  
  # Now the indexing will work correctly because data is aligned
  pheno_clean_final <- pheno_data[complete_samples, ]
  expr_clean <- combined_data$expr[, complete_samples]
  
  # Verify alignment is maintained
  cat("Post-filtering alignment check:", identical(colnames(expr_clean), pheno_clean_final$sample_name), "\n")
  
  # Create subsets
  younger_samples <- pheno_clean_final$age_group_final == "younger"
  older_samples <- pheno_clean_final$age_group_final == "older"
  
  cat("Final age group distribution:\n")
  print(table(pheno_clean_final$age_group_final, useNA = "ifany"))
  
  if("region" %in% colnames(pheno_clean_final)) {
    cat("Region distribution:\n")
    print(table(pheno_clean_final$region, pheno_clean_final$age_group_final))
  }
  
  # Check if we have valid TOD data
  tod_all <- pheno_clean_final$tod
  tod_younger <- pheno_clean_final$tod[younger_samples]
  tod_older <- pheno_clean_final$tod[older_samples]
  
  cat("TOD summary:\n")
  cat("- All TOD range:", range(tod_all, na.rm = TRUE), "\n")
  cat("- Younger TOD length:", length(tod_younger), "\n")
  cat("- Older TOD length:", length(tod_older), "\n")
  
  # Return organized data with standardized column names
  result <- list(
    all_expr = expr_clean,
    younger_expr = expr_clean[, younger_samples],
    older_expr = expr_clean[, older_samples],
    all_tod = tod_all,
    tod_younger = tod_younger,
    tod_older = tod_older,
    sample_info = pheno_clean_final
  )
  
  # Verify result structure
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



#######################
library(readxl)
library(VennDiagram)
file_path <- "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/Paper/Congruence/PNAS_aging/results/brain_regions/all_genes_merged.xlsx"
df <- read_excel(file_path, sheet = 1)
library(readxl)
library(VennDiagram)
library(gridExtra)
library(grid)
library(readxl)
library(VennDiagram)
library(gridExtra)
library(grid)
source("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Pipeline/one_cosinor_OLS_new.R")
# Load
file_path <- "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/Paper/Congruence/PNAS_aging/results/brain_regions/all_genes_merged.xlsx"
df <- read_excel(file_path, sheet = 1)

genes <- df$Gene
thresholds <- c(0.05, 0.01, 0.001)

plots_p <- list()
plots_q <- list()

# P-values
for (t in thresholds) {
  sig_Y <- genes[df$P_Value_Y < t]
  sig_O <- genes[df$P_Value_O < t]
  
  venn <- draw.pairwise.venn(
    area1 = length(sig_Y),
    area2 = length(sig_O),
    cross.area = length(intersect(sig_Y, sig_O)),
    category = c("Young", "Old"),
    fill = c("skyblue", "orange"),
    alpha = 0.5,
    scaled = FALSE,
    ext.text = FALSE,
    cat.col = c("skyblue4", "orange4"),
    cat.cex = 1.5,
    cex = 2
  )
  
  # Convert to grob and add title above
  g <- arrangeGrob(
    grobs = venn,
    top = textGrob(paste0("p < ", t), gp = gpar(fontsize = 20, fontface = "bold"))
  )
  plots_p <- c(plots_p, list(g))
}

# Q-values
for (t in thresholds) {
  sig_Y <- genes[df$Q_Value_Y < t]
  sig_O <- genes[df$Q_Value_O < t]
  
  venn <- draw.pairwise.venn(
    area1 = length(sig_Y),
    area2 = length(sig_O),
    cross.area = length(intersect(sig_Y, sig_O)),
    category = c("Young", "Old"),
    fill = c("lightgreen", "pink"),
    alpha = 0.5,
    scaled = FALSE,
    ext.text = FALSE,
    cat.col = c("darkgreen", "red"),
    cat.cex = 1.5,
    cex = 2
  )
  
  g <- arrangeGrob(
    grobs = venn,
    top = textGrob(paste0("q < ", t), gp = gpar(fontsize = 20, fontface = "bold"))
  )
  plots_q <- c(plots_q, list(g))
}

# Arrange 2x3
pdf("venn_young_old_2x3.pdf", width = 18, height = 12)
grid.arrange(grobs = c(plots_p, plots_q), nrow = 2, ncol = 3)
dev.off()


###########################
library(ggplot2)
library(gridExtra)

plot_gene_cosinor <- function(gene_list, COMBINED_data, df, period = 24, alpha = 0.05) {
  
  plot_one <- function(gene, group = "Y") {
    # Select phenotype subset
    if (group == "Y") {
      expr_vec <- as.numeric(COMBINED_data$younger_expr[gene, ])
      tod_vec  <- COMBINED_data$tod_younger
      pval     <- df$P_Value_Y[df$Gene == gene]
      qval     <- df$Q_Value_Y[df$Gene == gene]
      bayes_phase <- df$Bayesian_Median_Y[df$Gene == gene]
      bayes_sd    <- df$Bayesian_SD_Y[df$Gene == gene]
      rho         <- df$Posterior_Rho_Y[df$Gene == gene]
      BFDR_conserve        <- df$BFDR_conserve[df$Gene == gene]
      BFDR_shift        <- df$BFDR_shift[df$Gene == gene]
    } else {
      expr_vec <- as.numeric(COMBINED_data$older_expr[gene, ])
      tod_vec  <- COMBINED_data$tod_older
      pval     <- df$P_Value_O[df$Gene == gene]
      qval     <- df$Q_Value_O[df$Gene == gene]
      bayes_phase <- df$Bayesian_Median_O[df$Gene == gene]
      bayes_sd    <- df$Bayesian_SD_O[df$Gene == gene]
      rho         <- df$Posterior_Rho_O[df$Gene == gene]
      BFDR_conserve        <- df$BFDR_conserve[df$Gene == gene]
      BFDR_shift        <- df$BFDR_shift[df$Gene == gene]
    }
    
    # Fit cosinor
    fit_cos <- one_cosinor_OLS(tod = tod_vec, y = expr_vec, alpha = alpha, period = period)
    
    # Annotation text
    annot_text <- sprintf(
      "Cosinor Phase = %.2f ± %.2f h; Bayes Phase = %.2f ± %.2f h\np = %.5g; q = %.5f; rho = %.5f",
      fit_cos$peak, fit_cos$phase$sd * 24/(2*pi),
      bayes_phase, bayes_sd, pval, qval, rho
    )

    # Predicted curve
    tod_grid <- seq(min(tod_vec), max(tod_vec), length.out = 200)
    pred <- fit_cos$M$est + fit_cos$A$est * cos(2*pi*tod_grid/period + fit_cos$phase$est)
    
    plot_df <- data.frame(TOD = tod_vec, Expression = expr_vec)
    pred_df <- data.frame(TOD = tod_grid, Expression = pred)
    
    # Make plot
    ggplot(plot_df, aes(x = TOD, y = Expression)) +
      geom_point(size = 2) +
      geom_line(data = pred_df, aes(x = TOD, y = Expression), color = "red", size = 1) +
      labs(title = paste("Human,", gene, "-", ifelse(group == "Y", "Younger", "Older")),
           subtitle = annot_text,
           x = "TOD", y = "Expression") +
      theme_bw(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(size = 10, hjust = 0.5))  # centered subtitle
  }
  
  # Build list of plots: each gene has 2 (Y, O)
  plot_list <- list()
  for (g in gene_list) {
    plot_list <- c(plot_list, list(plot_one(g, "Y"), plot_one(g, "O")))
  }
  
  # Arrange 2 per row (Young vs Old side by side)
  grid.arrange(grobs = plot_list, ncol = 2)
}

plot_gene_cosinor(c("NR1D2", "BHLHE41", "BHLHE40", "BMAL1", "CRY1", "PER1", "PER3", "NR1D1", "PER2"), COMBINED_data, df)



