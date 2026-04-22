#── Clear environment ───────────────────────────────────────────────────────────────
rm(list=ls())

#── Paths ─────────────────────────────────────────────────────────────────────────
current_gtex <- "/home/qtp1/Projects/Collaborative"
current_wd   <- "/home/qtp1/Projects/Circadian"
current_aging <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging"

#── Load BA11 and BA47 data ──────────────────────────────────────────────────────
BA11 <- readRDS(file.path(current_aging, "data/BA11_data.rds"))
BA47 <- readRDS(file.path(current_aging, "data/BA47_data.rds"))

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
# Add the cosinor method results here
source("/home/qtp1/Projects/Pipeline/one_cosinor_OLS_new.R")


#── Source R scripts ─────────────────────────────────────────────────────────────
WD <- "Kyle/Circadian-analysis-main/R/v1"
setwd(file.path(current_wd, WD))
scripts <- list.files("R", full.names=TRUE)
sapply(scripts, source)

#── MCMC runner function ─────────────────────────────────────────────────────────
run_MCMC <- function(data, tod, variable_name, save_path,
                     n.iter=2500, n.burn=500, P=24, seed=1) {
  dat_list <- list(data = as.data.frame(data),
                   time = tod,
                   gname = rownames(data))
  a.init <- CBt_init_single(Data.list=dat_list, P=P,
                            FitCosinor=TRUE,
                            mu_M=0, sigma_M=10,
                            mu_A=1, sigma_A=10,
                            seed=seed)
  CB.res <- CB_MCMC_single_rj_slice(Data.list=dat_list,
                                    Init.value=a.init,
                                    P=P, iteration=n.iter,
                                    thin=20, n.burn=n.burn,
                                    seed=seed,
                                    p_rhythmic=rep(0.2,nrow(data)),
                                    rj.p.stay=0.5,
                                    A_prior="trunc_Normal_OLS_condi",
                                    mu_A=1, sigma_A=10^2, A.min=0,
                                    A_wb_beta2=2,
                                    A_gm_shape=1.99, A_gm_rate=0.5,
                                    rj.phi=TRUE, rj.A=TRUE,
                                    mu_M=0, sigma_M=10^2,
                                    sigma_prior_v=2, sigma_prior_s=0)
  file_name <- paste0(variable_name, "_bay_", seed, ".RDS")
  file_path <- file.path(save_path, file_name)
  saveRDS(CB.res, file=file_path)
  if (!file.exists(file_path)) stop("Error: File not saved correctly")
  return(file_name)
}

#── Data preparation function ────────────────────────────────────────────────────
prepare_brain_data <- function(brain_data, region_name) {
  # Extract IDs from column names
  expr_ids <- sub(paste0(".*", region_name, "-([0-9]+)\\.CEL\\.gz"), "\\1", 
                  colnames(brain_data$expr))
  expr_ids <- as.integer(expr_ids)
  names(expr_ids) <- colnames(brain_data$expr)
  
  # Get younger and older IDs
  younger_ids <- brain_data$pheno$ID[brain_data$pheno$AgeGroup == "younger"]
  older_ids <- brain_data$pheno$ID[brain_data$pheno$AgeGroup == "older"]
  
  # Subset expression data
  younger_expr <- brain_data$expr[, expr_ids %in% younger_ids]
  older_expr <- brain_data$expr[, expr_ids %in% older_ids]
  
  # Create TOD mapping
  tod_map <- setNames(brain_data$pheno$TOD, brain_data$pheno$ID)
  
  # Get TOD for younger samples
  ids_younger <- sub(paste0(".*", region_name, "-([0-9]+)\\.CEL\\.gz"), "\\1", 
                     colnames(younger_expr))
  tod_younger <- tod_map[ids_younger]
  
  # Get TOD for older samples
  ids_older <- sub(paste0(".*", region_name, "-([0-9]+)\\.CEL\\.gz"), "\\1", 
                   colnames(older_expr))
  tod_older <- tod_map[ids_older]
  
  # Return organized data
  list(
    all_expr = brain_data$expr,
    younger_expr = younger_expr,
    older_expr = older_expr,
    all_tod = tod_map[as.character(expr_ids)],
    tod_younger = tod_younger,
    tod_older = tod_older
  )
}


#── Prepare data for both brain regions ─────────────────────────────────────────
BA11_data <- prepare_brain_data(BA11, "BA11")
BA47_data <- prepare_brain_data(BA47, "BA47")

#── Define analysis configurations ──────────────────────────────────────────────
# Create analysis groups
analysis_groups <- list(
  # Full datasets
  BA11_all = list(data = BA11_data$all_expr, tod = BA11_data$all_tod, name = "BA11_all"),
  BA47_all = list(data = BA47_data$all_expr, tod = BA47_data$all_tod, name = "BA47_all"),
  
  # Age-stratified BA11
  BA11_younger = list(data = BA11_data$younger_expr, tod = BA11_data$tod_younger, name = "BA11_younger"),
  BA11_older = list(data = BA11_data$older_expr, tod = BA11_data$tod_older, name = "BA11_older"),
  
  # Age-stratified BA47
  BA47_younger = list(data = BA47_data$younger_expr, tod = BA47_data$tod_younger, name = "BA47_younger"),
  BA47_older = list(data = BA47_data$older_expr, tod = BA47_data$tod_older, name = "BA47_older")
)

#── Set up result directories ───────────────────────────────────────────────────
base_result_dir <- file.path(current_aging, "results", "brain_regions")
num_cores <- 6

#── Processing function ──────────────────────────────────────────────────────────
process_brain_group <- function(group_info) {
  group_name <- group_info$name
  tryCatch({
    # Create directory
    save_path <- file.path(base_result_dir, group_name)
    dir.create(save_path, recursive=TRUE, showWarnings=FALSE)
    
    # Run MCMC
    result_file <- run_MCMC(data = group_info$data,
                            tod = group_info$tod,
                            variable_name = group_name,
                            save_path = save_path)
    
    cat("Completed:", group_name, "- File:", result_file, "\n")
    return(result_file)
    
  }, error=function(e){
    error_msg <- paste(Sys.time(), "- Error in", group_name, ":\n", conditionMessage(e), "\n")
    cat(error_msg, file=file.path(base_result_dir, "error_log.txt"), append=TRUE)
    cat(error_msg)
    return(e)
  })
}


#── Run MCMC analysis for all groups ────────────────────────────────────────────
cat("Starting MCMC analysis for", length(analysis_groups), "groups...\n")
results <- mclapply(analysis_groups, process_brain_group, mc.cores=num_cores)

#── Load and summarize results ──────────────────────────────────────────────────
cat("\n=== Loading and summarizing results ===\n")

# Function to load results
load_brain_results <- function() {
  result_data <- list()
  
  for (group_name in names(analysis_groups)) {
    result_path <- file.path(base_result_dir, group_name)
    rds_files <- list.files(result_path, pattern = "\\.RDS$", full.names = TRUE)
    
    if (length(rds_files) == 1) {
      result_data[[group_name]] <- readRDS(rds_files[1])
      cat("Loaded:", group_name, "\n")
    } else {
      warning(paste("Group", group_name, "has", length(rds_files), "RDS files (expected 1)"))
    }
  }
  
  return(result_data)
}

# Load results
brain_results <- load_brain_results()

# Summarize rhythmic genes
brain_summary <- lapply(brain_results, function(x) {
  summarize_bay(x$rho, BF = 1, p_rhythmic = 0.2)
})

cat("\n=== Rhythmic Gene Counts ===\n")
for (bf_threshold in c(1, 3, 5 )) {
  cat(paste("Bayes Factor >", bf_threshold, ":\n"))
  counts <- sapply(names(brain_summary), function(group) {
    sum(brain_summary[[group]]$BayesF > bf_threshold)
  })
  
  result_string <- paste0(names(counts), ": ", counts, collapse = " | ")
  cat(result_string, "\n\n")
}

# Bayes Factor > 1 :
#  BA11_all: 1388 | BA11_younger: 1370 | BA11_older: 1519 
#  BA47_all: 821  | BA47_younger: 2208 | BA47_older: 1255
# 
# Bayes Factor > 3 :
#   BA11_all: 597 | BA11_younger: 437 | BA11_older: 431 
#   BA47_all: 354 | BA47_younger: 723 | BA47_older: 383
# 
# Bayes Factor > 5 :
# BA11_all: 402 | BA11_younger: 257 | BA11_older: 264 |
# BA47_all: 247 BA47_younger: 438 | BA47_older: 235 

#── Save final results ──────────────────────────────────────────────────────────
final_results <- list(
  brain_results = brain_results,
  brain_summary = brain_summary,
  analysis_groups = analysis_groups
)

saveRDS(final_results, file.path(base_result_dir, "final_brain_circadian_results.RDS"))

#── Example of how to use the saved results ────────────────────────────────────
cat("\n=== Usage Examples ===\n")
cat("# To load results later:\n")
cat("final_results <- readRDS('path/to/final_brain_circadian_results.RDS')\n\n")

cat("# To get rhythmic genes for BA11_younger with BF > 3:\n")
cat("rhythmic_genes <- final_results$brain_summary$BA11_younger$BayesF > 3\n")
cat("gene_names <- names(rhythmic_genes)[rhythmic_genes]\n\n")

cat("# To get raw MCMC rho values for concordance analysis:\n")
cat("rho_BA11 <- final_results$brain_results$BA11_all$rho\n")
cat("rho_BA47 <- final_results$brain_results$BA47_all$rho\n")
cat("concordance_value <- concordance(rho_BA11, rho_BA47)\n\n")

cat("# To access original expression data:\n")
cat("ba11_expr <- final_results$analysis_groups$BA11_all$data\n")
cat("ba11_tod <- final_results$analysis_groups$BA11_all$tod\n")

#── Create meta-result  ────────────────────────────────────
# Simple function to get just the adjusted concordance value
get_adj_concordance <- function(data1, data2, name1, name2) {
  if (name1 %in% names(data1) && name2 %in% names(data2)) {
    result <- concordance(data1[[name1]]$rho, data2[[name2]]$rho)
    return(result$concordance_adj)
  } else {
    return(NA)
  }
}


# Define conditions
conditions <- c("BA11_all", "BA11_younger", "BA11_older", 
                "BA47_all", "BA47_younger", "BA47_older")

# Create empty matrix
n <- length(conditions)
concordance_table <- matrix(NA, nrow = n, ncol = n, 
                            dimnames = list(conditions, conditions))

# Fill diagonal with 1.0 (perfect self-concordance)
diag(concordance_table) <- 1.0

# Calculate all pairwise concordances
cat("Calculating adjusted concordances...\n")

# Manual calculation of each pair (since we know they work individually)
concordance_table["BA11_all", "BA11_younger"] <- get_adj_concordance(brain_results, brain_results, "BA11_all", "BA11_younger")
concordance_table["BA11_all", "BA11_older"] <- get_adj_concordance(brain_results, brain_results, "BA11_all", "BA11_older")
concordance_table["BA11_all", "BA47_all"] <- get_adj_concordance(brain_results, brain_results, "BA11_all", "BA47_all")
concordance_table["BA11_all", "BA47_younger"] <- get_adj_concordance(brain_results, brain_results, "BA11_all", "BA47_younger")
concordance_table["BA11_all", "BA47_older"] <- get_adj_concordance(brain_results, brain_results, "BA11_all", "BA47_older")

concordance_table["BA11_younger", "BA11_older"] <- get_adj_concordance(brain_results, brain_results, "BA11_younger", "BA11_older")
concordance_table["BA11_younger", "BA47_all"] <- get_adj_concordance(brain_results, brain_results, "BA11_younger", "BA47_all")
concordance_table["BA11_younger", "BA47_younger"] <- get_adj_concordance(brain_results, brain_results, "BA11_younger", "BA47_younger")
concordance_table["BA11_younger", "BA47_older"] <- get_adj_concordance(brain_results, brain_results, "BA11_younger", "BA47_older")

concordance_table["BA11_older", "BA47_all"] <- get_adj_concordance(brain_results, brain_results, "BA11_older", "BA47_all")
concordance_table["BA11_older", "BA47_younger"] <- get_adj_concordance(brain_results, brain_results, "BA11_older", "BA47_younger")
concordance_table["BA11_older", "BA47_older"] <- get_adj_concordance(brain_results, brain_results, "BA11_older", "BA47_older")

concordance_table["BA47_all", "BA47_younger"] <- get_adj_concordance(brain_results, brain_results, "BA47_all", "BA47_younger")
concordance_table["BA47_all", "BA47_older"] <- get_adj_concordance(brain_results, brain_results, "BA47_all", "BA47_older")

concordance_table["BA47_younger", "BA47_older"] <- get_adj_concordance(brain_results, brain_results, "BA47_younger", "BA47_older")

# Make the matrix symmetric
for(i in 1:n) {
  for(j in 1:n) {
    if(i > j && !is.na(concordance_table[j, i])) {
      concordance_table[i, j] <- concordance_table[j, i]
    }
  }
}

# Print the table
cat("\n=== ADJUSTED CONCORDANCE TABLE ===\n")
print(round(concordance_table, 4))

# Create formatted table
cat("\n=== FORMATTED TABLE ===\n")
cat("              ")
for(col in conditions) {
  cat(sprintf("%12s", col))
}
cat("\n")

for(i in 1:n) {
  cat(sprintf("%-12s  ", conditions[i]))
  for(j in 1:n) {
    if(is.na(concordance_table[i,j])) {
      cat(sprintf("%12s", "NA"))
    } else {
      cat(sprintf("%12.4f", concordance_table[i,j]))
    }
  }
  cat("\n")
}

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
upper_tri <- concordance_table[upper.tri(concordance_table)]
upper_tri <- upper_tri[!is.na(upper_tri)]

cat(sprintf("Number of pairwise comparisons: %d\n", length(upper_tri)))
cat(sprintf("Mean adjusted concordance: %.4f\n", mean(upper_tri)))
cat(sprintf("Median adjusted concordance: %.4f\n", median(upper_tri)))
cat(sprintf("Range: %.4f to %.4f\n", min(upper_tri), max(upper_tri)))
cat(sprintf("Standard deviation: %.4f\n", sd(upper_tri)))

# Find highest and lowest concordances
max_val <- max(upper_tri)
min_val <- min(upper_tri)

# Find positions of max and min values
max_indices <- which(concordance_table == max_val & upper.tri(concordance_table), arr.ind = TRUE)
min_indices <- which(concordance_table == min_val & upper.tri(concordance_table), arr.ind = TRUE)

cat(sprintf("\nHighest adjusted concordance: %.4f between %s and %s\n", 
            max_val, 
            conditions[max_indices[1,1]], 
            conditions[max_indices[1,2]]))

cat(sprintf("Lowest adjusted concordance: %.4f between %s and %s\n", 
            min_val, 
            conditions[min_indices[1,1]], 
            conditions[min_indices[1,2]]))

# Save table to CSV if desired
write.csv(concordance_table, "adjusted_concordance_table.csv", row.names = TRUE)


##################################
# Add this after loading BA11 and BA47 data
COMBINED <- readRDS(file.path(current_aging, "data/combined_data.rds"))

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


# Check the data structure
cat("Combined dataset summary:\n")
cat("Total samples:", ncol(COMBINED_data$all_expr), "\n")
cat("Younger samples:", ncol(COMBINED_data$younger_expr), "\n")
cat("Older samples:", ncol(COMBINED_data$older_expr), "\n")
cat("Sample distribution by region and age:\n")
print(table(COMBINED_data$sample_info$region, COMBINED_data$sample_info$age_group))

analysis_groups_extended <- list(
  #Combined datasets
  COMBINED_all = list(data = COMBINED_data$all_expr, tod = COMBINED_data$all_tod, name = "COMBINED_all"),
  COMBINED_younger = list(data = COMBINED_data$younger_expr, tod = COMBINED_data$tod_younger, name = "COMBINED_younger"),
  COMBINED_older = list(data = COMBINED_data$older_expr, tod = COMBINED_data$tod_older, name = "COMBINED_older")
)

cat("Starting extended MCMC analysis for", length(analysis_groups_extended), "groups...\n")
results_extended <- mclapply(analysis_groups_extended, process_brain_group, mc.cores=num_cores)


#── Load and summarize extended results ─────────────────────────────────────────
cat("\n=== Loading and summarizing extended results ===\n")

# Function to load extended results
load_extended_results <- function() {
  result_data <- list()
  
  for (group_name in names(analysis_groups_extended)) {
    result_path <- file.path(base_result_dir, group_name)
    rds_files <- list.files(result_path, pattern = "\\.RDS$", full.names = TRUE)
    
    if (length(rds_files) == 1) {
      result_data[[group_name]] <- readRDS(rds_files[1])
      cat("Loaded:", group_name, "\n")
    } else {
      warning(paste("Group", group_name, "has", length(rds_files), "RDS files (expected 1)"))
    }
  }
  
  return(result_data)
}


# Load extended results
brain_results_extended <- load_extended_results()
brain_results_extended <- c(brain_results, brain_results_extended)
# Summarize rhythmic genes for extended results
brain_summary_extended <- lapply(brain_results_extended, function(x) {
  summarize_bay(x$rho, BF = 1, p_rhythmic = 0.2)
})

cat("\n=== Extended Rhythmic Gene Counts ===\n")
for (bf_threshold in c(1, 3, 5)) {
  cat(paste("Bayes Factor >", bf_threshold, ":\n"))
  counts <- sapply(names(brain_summary_extended), function(group) {
    sum(brain_summary_extended[[group]]$BayesF > bf_threshold)
  })
  
  # Print in organized format
  original_groups <- c("BA11_all", "BA11_younger", "BA11_older", "BA47_all", "BA47_younger", "BA47_older")
  combined_groups <- c("COMBINED_all", "COMBINED_younger", "COMBINED_older")
  
  cat("Original regions:\n")
  original_string <- paste0(original_groups, ": ", counts[original_groups], collapse = " | ")
  cat(original_string, "\n")
  
  cat("Combined regions:\n")
  combined_string <- paste0(combined_groups, ": ", counts[combined_groups], collapse = " | ")
  cat(combined_string, "\n\n")
}


# Bayes Factor > 1 :
#  COMBINED_all: 1160 | COMBINED_younger: 1906 | COMBINED_older: 1710
#  BA11_all: 1388 | BA11_younger: 1370 | BA11_older: 1519 
#  BA47_all: 821  | BA47_younger: 2208 | BA47_older: 1255
# 
# Bayes Factor > 3 :
#   COMBINED_all: 510 | COMBINED_younger: 832 | COMBINED_older: 689
#   BA11_all: 597 | BA11_younger: 437 | BA11_older: 431 
#   BA47_all: 354 | BA47_younger: 723 | BA47_older: 383
# 
# Bayes Factor > 5 :
# BF > 5: COMBINED_all: 350 | COMBINED_younger: 565 | COMBINED_older: 448
# BA11_all: 402 | BA11_younger: 257 | BA11_older: 264 |
# BA47_all: 247 BA47_younger: 438 | BA47_older: 235 


# BF > 1: COMBINED_all: 1160 | COMBINED_younger: 1906 | COMBINED_older: 1710
# BF > 3: COMBINED_all: 510 | COMBINED_younger: 832 | COMBINED_older: 689
# BF > 5: COMBINED_all: 350 | COMBINED_younger: 565 | COMBINED_older: 448

#── Extended concordance analysis ───────────────────────────────────────────────
# Extended conditions including combined datasets
conditions_extended <- c("BA11_all", "BA11_younger", "BA11_older", 
                         "BA47_all", "BA47_younger", "BA47_older",
                         "COMBINED_all", "COMBINED_younger", "COMBINED_older")

# Create extended concordance matrix
n_ext <- length(conditions_extended)
concordance_table_extended <- matrix(NA, nrow = n_ext, ncol = n_ext, 
                                     dimnames = list(conditions_extended, conditions_extended))

# Fill diagonal with 1.0 (perfect self-concordance)
diag(concordance_table_extended) <- 1.0

cat("Calculating extended adjusted concordances...\n")

# Function to safely calculate concordance
safe_concordance <- function(name1, name2) {
  tryCatch({
    if (name1 %in% names(brain_results_extended) && name2 %in% names(brain_results_extended)) {
      result <- concordance(brain_results_extended[[name1]]$rho, brain_results_extended[[name2]]$rho)
      return(result$concordance_adj)
    } else {
      return(NA)
    }
  }, error = function(e) {
    cat("Error calculating concordance between", name1, "and", name2, ":", e$message, "\n")
    return(NA)
  })
}

# Calculate all pairwise concordances for extended dataset
for(i in 1:(n_ext-1)) {
  for(j in (i+1):n_ext) {
    name1 <- conditions_extended[i]
    name2 <- conditions_extended[j]
    concordance_val <- safe_concordance(name1, name2)
    concordance_table_extended[i, j] <- concordance_val
    concordance_table_extended[j, i] <- concordance_val  # Make symmetric
    
    if (!is.na(concordance_val)) {
      cat(sprintf("Concordance %s vs %s: %.4f\n", name1, name2, concordance_val))
    }
  }
}

# Print the extended concordance table
cat("\n=== EXTENDED ADJUSTED CONCORDANCE TABLE ===\n")
print(round(concordance_table_extended, 4))

# Extended summary statistics
cat("\n=== EXTENDED SUMMARY STATISTICS ===\n")
upper_tri_ext <- concordance_table_extended[upper.tri(concordance_table_extended)]
upper_tri_ext <- upper_tri_ext[!is.na(upper_tri_ext)]

cat(sprintf("Number of pairwise comparisons: %d\n", length(upper_tri_ext)))
cat(sprintf("Mean adjusted concordance: %.4f\n", mean(upper_tri_ext)))
cat(sprintf("Median adjusted concordance: %.4f\n", median(upper_tri_ext)))
cat(sprintf("Range: %.4f to %.4f\n", min(upper_tri_ext), max(upper_tri_ext)))
cat(sprintf("Standard deviation: %.4f\n", sd(upper_tri_ext)))

# Specific comparisons of interest for combined data
cat("\n=== KEY COMPARISONS WITH COMBINED DATA ===\n")
key_comparisons <- list(
  c("BA11_all", "COMBINED_all"),
  c("BA47_all", "COMBINED_all"), 
  c("BA11_younger", "COMBINED_younger"),
  c("BA47_younger", "COMBINED_younger"),
  c("BA11_older", "COMBINED_older"),
  c("BA47_older", "COMBINED_older"),
  c("COMBINED_younger", "COMBINED_older")
)

for(comparison in key_comparisons) {
  name1 <- comparison[1]
  name2 <- comparison[2]
  val <- concordance_table_extended[name1, name2]
  if(!is.na(val)) {
    cat(sprintf("%s vs %s: %.4f\n", name1, name2, val))
  }
}

#── Save extended results ───────────────────────────────────────────────────────
final_results_extended <- list(
  brain_results = brain_results_extended,
  brain_summary = brain_summary_extended,
  analysis_groups = analysis_groups_extended,
  concordance_table = concordance_table_extended,
  combined_sample_info = COMBINED_data$sample_info
)

saveRDS(final_results_extended, file.path(base_result_dir, "final_brain_circadian_results_extended.RDS"))

# Save extended concordance table
write.csv(concordance_table_extended, 
          file.path(base_result_dir, "extended_concordance_table.csv"), 
          row.names = TRUE)

#################################################################################
library(circular)

output = "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/output"
# Load the extended Bayesian results
final_results_extended <- readRDS(file.path(base_result_dir, "final_brain_circadian_results_extended.RDS"))

# Load cosinor results
com_older <- readxl::read_excel(file.path(output, "combined_cosinor_older.xlsx"), sheet = 1)
com_younger <- readxl::read_excel(file.path(output, "combined_cosinor_younger.xlsx"), sheet = 1)

# Extract Bayesian results 
bayesian_combined_younger <- final_results_extended$brain_results$COMBINED_younger
bayesian_combined_older <- final_results_extended$brain_results$COMBINED_older

# Extract summary statistics
summary_combined_younger <- final_results_extended$brain_summary$COMBINED_younger
summary_combined_older <- final_results_extended$brain_summary$COMBINED_older

rownames(bayesian_combined_younger$phi) <- rownames(summary_combined_younger)
rownames(bayesian_combined_older$phi) <- rownames(summary_combined_older)

cat("Younger phi matrix rownames (first 5):", head(rownames(bayesian_combined_younger$phi), 5), "\n")
cat("Older phi matrix rownames (first 5):", head(rownames(bayesian_combined_older$phi), 5), "\n")

add_bayesian_info <- function(cosinor_df, bayesian_phi_matrix) {
  
  # Use rownames as gene identifiers
  gene_ids <- rownames(bayesian_phi_matrix)
  
  # Match cosinor_df Gene to rownames
  matched_rows <- match(cosinor_df$Gene, gene_ids)
  
  # For each gene, compute Bayesian circular median and variance
  bayesian_median_vec <- rep(NA_real_, nrow(cosinor_df))
  bayesian_variance_vec <- rep(NA_real_, nrow(cosinor_df))
  
  for (i in seq_len(nrow(cosinor_df))) {
    idx <- matched_rows[i]
    if (!is.na(idx)) {
      posterior_samples <- bayesian_phi_matrix[idx, ]
      bayesian_median_vec[i] <- circular_median(posterior_samples)
      bayesian_variance_vec[i] <- circular_variance(posterior_samples)
    }
  }
  
  # Add new columns
  cosinor_df <- cosinor_df %>%
    mutate(
      Bayesian_Median = bayesian_median_vec,
      Bayesian_SD = sqrt(bayesian_variance_vec)
    )
  
  return(cosinor_df)
}

#── Helper function to add posterior using rownames ────────────────────────────────
add_posterior <- function(cosinor_df, summarize_df) {
  
  # Extract posterior probability using rownames
  posterior_vec <- summarize_df$RowAverage
  summarize_genes <- rownames(summarize_df)
  
  cosinor_df <- cosinor_df %>%
    mutate(
      Posterior_Rho = posterior_vec[match(Gene, summarize_genes)]
    )
  
  return(cosinor_df)
}



#── Helper function to add Rhythmicity and BayesF using rownames ───────────────────
add_summary_metrics <- function(cosinor_df, summary_result) {
  
  # Add Rhythmicity and BayesF using rownames for matching
  cosinor_df <- cosinor_df %>%
    mutate(
      Rhythmicity = summary_result$Rhythmicity[match(Gene, rownames(summary_result))],
      BayesF = summary_result$BayesF[match(Gene, rownames(summary_result))]
    )
  
  return(cosinor_df)
}

circular_median <- function(samples, P = 24) {
  samples_rad <- circular((samples / P) * (2 * pi), units = "radians")
  median_rad <- median.circular(samples_rad, type = "median")
  median_time <- (median_rad / (2 * pi)) * P
  return(median_time %% P)
}

# Your circular variance function
circular_variance <- function(samples, P = 24) {
  samples_rad <- circular((samples / P) * (2 * pi), units = "radians")
  var_circ <- 1 - rho.circular(samples_rad)  # 1 - mean resultant length
  return(var_circ)
}

#── Modified comprehensive function using rownames ─────────────────────────────────
create_comprehensive_combined <- function(cosinor_df, bayesian_result, summary_result, age_group) {
  
  cat(paste("Processing", age_group, "group...\n"))
  
  # Add original Bayesian information (circular median and variance)
  combined_df <- add_bayesian_info(cosinor_df, bayesian_result$phi)
  combined_df <- add_posterior(combined_df, summary_result)

  
  # NEW: Add Rhythmicity and BayesF from summary
  combined_df <- add_summary_metrics(combined_df, summary_result)
  
  # Add Age Group
  combined_df <- combined_df %>%
    mutate(
      Age_Group = age_group
    )
  
  return(combined_df)
}

comprehensive_younger <- create_comprehensive_combined(
  com_younger, 
  bayesian_combined_younger, 
  summary_combined_younger, 
  "Younger"
)

comprehensive_older <- create_comprehensive_combined(
  com_older, 
  bayesian_combined_older, 
  summary_combined_older, 
  "Older"
)

library(openxlsx)

# Write comprehensive younger results
write.xlsx(comprehensive_younger, 
           file = file.path(output, "comprehensive_combined_younger.xlsx"),
           sheetName = "Younger_Combined",
           rowNames = FALSE)

# Write comprehensive older results  
write.xlsx(comprehensive_older,
           file = file.path(output, "comprehensive_combined_older.xlsx"), 
           sheetName = "Older_Combined",
           rowNames = FALSE)
####################################################################

library(ggplot2)
library(dplyr)
library(viridis)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
# pathway_genes <- c("CRY1", "CRY2", "ARNTL", "PER1", "BHLHE41", "BHLHE40", "DBP",
#                    "PER3", "PER2", "NR1D1", "NR1D2", "BMAL1", "CLOCK")
# Choose pathway
pathway_name <- "KEGG Caffeine metabolism"
pathway_genes <- kegg.pathway.list_hsa[[pathway_name]] 

# Get younger and older data
younger_data <- analysis_groups_extended$COMBINED_younger
older_data <- analysis_groups_extended$COMBINED_older

# Find overlap between pathway genes and expression data
all_genes <- rownames(younger_data$data)
overlap_genes <- intersect(all_genes, pathway_genes)
cat("Found", length(overlap_genes), "pathway genes in expression data\n")

# Extract expression data for pathway genes
younger_expr <- younger_data$data[overlap_genes, , drop = FALSE]
older_expr <- older_data$data[overlap_genes, , drop = FALSE]

# Get TOD information
younger_tod <- younger_data$tod
older_tod <- older_data$tod

# Remove duplicates and use even TOD sampling
get_unique_samples <- function(tod_vector, expr_matrix) {
  unique_tod_idx <- !duplicated(tod_vector)
  return(list(
    tod = tod_vector[unique_tod_idx],
    expr = expr_matrix[, unique_tod_idx, drop = FALSE]
  ))
}

select_even_tod <- function(tod_vector, n_samples) {
  tod_order <- order(tod_vector)
  step_size <- length(tod_vector) / n_samples
  selected_positions <- round(seq(1, length(tod_vector), by = step_size))[1:n_samples]
  return(tod_order[selected_positions])
}

# Get unique samples for each age group
younger_unique <- get_unique_samples(younger_tod, younger_expr)
older_unique <- get_unique_samples(older_tod, older_expr)

cat("Unique younger samples:", length(younger_unique$tod), "\n")
cat("Unique older samples:", length(older_unique$tod), "\n")

# Determine equal sample size
n_samples <- min(31, length(younger_unique$tod), length(older_unique$tod))
cat("Using", n_samples, "samples per age group\n")

# Even sampling across TOD
younger_subset_idx <- select_even_tod(younger_unique$tod, n_samples)
older_subset_idx <- select_even_tod(older_unique$tod, n_samples)

# Extract final subset
younger_expr_final <- younger_unique$expr[, younger_subset_idx, drop = FALSE]
older_expr_final <- older_unique$expr[, older_subset_idx, drop = FALSE]
younger_tod_final <- younger_unique$tod[younger_subset_idx]
older_tod_final <- older_unique$tod[older_subset_idx]

# Convert TOD from current range to 0-24 hours
convert_tod_to_24h <- function(tod_vector) {
  # Add 6 to shift the range, then use modulo to wrap around
  converted <- (tod_vector + 6) %% 24
  return(converted)
}

# Apply conversion to both age groups
younger_tod_24h <- convert_tod_to_24h(younger_tod_final)
older_tod_24h <- convert_tod_to_24h(older_tod_final)

# Order by the converted TOD values
younger_order <- order(younger_tod_24h)
older_order <- order(older_tod_24h)

younger_expr_ordered <- younger_expr_final[, younger_order, drop = FALSE]
older_expr_ordered <- older_expr_final[, older_order, drop = FALSE]

# First create combined data WITHOUT gap for clustering
combined_expr_for_clustering <- cbind(younger_expr_ordered, older_expr_ordered)

# Do the clustering on data without NAs
hclust_result <- hclust(dist(combined_expr_for_clustering, method = "euclidean"), method = "complete")

# Create a separator column of NAs
separator_col <- matrix(NA, nrow = nrow(younger_expr_ordered), ncol = 3)
colnames(separator_col) <- paste0("separator_", 1:3)

# Now combine with separator for display
combined_expr_ordered <- cbind(younger_expr_ordered, separator_col, older_expr_ordered)

# Update TOD and age group vectors to include separator
separator_tod <- rep(NA, 3)
separator_age <- rep("Separator", 3)

combined_tod_ordered <- c(younger_tod_24h[younger_order], separator_tod, older_tod_24h[older_order])
combined_age_ordered <- c(rep("Younger", ncol(younger_expr_ordered)), 
                          separator_age,
                          rep("Older", ncol(older_expr_ordered)))

# Create annotation that treats separator samples differently
sample_annotation <- data.frame(
  TOD = ifelse(is.na(combined_tod_ordered), NA, combined_tod_ordered),
  AgeGroup = ifelse(combined_age_ordered == "Separator", NA, combined_age_ordered)
)
rownames(sample_annotation) <- colnames(combined_expr_ordered)

# Only include actual age groups in annotation colors
ann_colors <- list(
  AgeGroup = c("Younger" = "#2E86AB", "Older" = "#A23B72"),
  TOD = colorRampPalette(c("#440154", "#31688e", "#35b779", "#fde725"))(24)
)

# Save as PDF
pdf(
  file = paste0("/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/results/brain_regions/plot/", 
                pathway_name, "_heatmap_0to24.pdf"),
  width = 14,
  height = 8
)

pheatmap(combined_expr_ordered,
         annotation_col = sample_annotation,
         annotation_colors = ann_colors,
         cluster_rows = hclust_result,
         cluster_cols = FALSE,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_colnames = FALSE,
         main = paste("Expression across TOD:", pathway_name, "(0-24h Scale)"),
         fontsize = 10,
         fontsize_row = 8,
         na_col = "transparent",
         border_color = NA)

dev.off()

cat("Heatmap saved successfully!\n")



###########


