################################################################################
# ALL PAIRWISE TISSUE COMPARISONS WITHIN human
################################################################################

library(dplyr)
library(ggplot2)

# Load your data
load("/home/qtp1/Projects/Collaborative/GTEXdata/result/summary/hb/mcmc_rho_BF3.RData")

################################################################################
# SETUP OUTPUT DIRECTORY
################################################################################

output_dir <- "/home/qtp1/Projects/Collaborative/GTEXdata/result/summary/hb/human_pairwise_concordance"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\nOutput directory:", output_dir, "\n")

################################################################################
# PREPARE DATA
################################################################################

prepare_mcmc_data <- function(mcmc_data_list) {
  cat("Preparing", length(mcmc_data_list), "datasets...\n")
  
  for (i in seq_along(mcmc_data_list)) {
    dataset_name <- names(mcmc_data_list)[i]
    
    if (is.null(attr(mcmc_data_list[[i]], "symbols"))) {
      if (!is.null(attr(mcmc_data_list[[i]], "dimnames")[[1]])) {
        attr(mcmc_data_list[[i]], "symbols") <- attr(mcmc_data_list[[i]], "dimnames")[[1]]
      } else if (!is.null(rownames(mcmc_data_list[[i]]))) {
        attr(mcmc_data_list[[i]], "symbols") <- rownames(mcmc_data_list[[i]])
      }
    }
  }
  
  return(mcmc_data_list)
}

create_mcmc_list_format <- function(rho_matrix) {
  if (is.list(rho_matrix) && !is.null(rho_matrix$rho)) {
    return(rho_matrix)
  } else {
    return(list(rho = rho_matrix))
  }
}

cat("\nPreparing human data:\n")
mcmc_data_human_prepared <- prepare_mcmc_data(mcmc_data_human)

################################################################################
# GET ALL TISSUE NAMES
################################################################################

all_tissues <- names(mcmc_data_human_prepared)
n_tissues <- length(all_tissues)
n_comparisons <- choose(n_tissues, 2)

cat("\n")
cat("="  %>% rep(80) %>% paste0(collapse=""), "\n")
cat("human PAIRWISE TISSUE CONCORDANCE ANALYSIS\n")
cat("="  %>% rep(80) %>% paste0(collapse=""), "\n")
cat("\n")
cat(sprintf("Number of tissues: %d\n", n_tissues))
cat(sprintf("Number of pairwise comparisons: %d\n", n_comparisons))
cat("\n")
cat("Tissues:\n")
print(all_tissues)
cat("\n")

################################################################################
# RUN ALL PAIRWISE COMPARISONS
################################################################################

cat("\n")
cat("Starting pairwise comparisons...\n")
cat("This will take a while (estimated ", round(n_comparisons * 30 / 60), " minutes)\n\n")

# Storage
results_df <- data.frame(
  tissue1 = character(),
  tissue2 = character(),
  raw_concordance = numeric(),
  adjusted_concordance = numeric(),
  ci_lower = numeric(),
  ci_upper = numeric(),
  null_mean = numeric(),
  null_sd = numeric(),
  p_value = numeric(),
  z_score = numeric(),
  stringsAsFactors = FALSE
)

# Track progress
counter <- 0
start_time <- Sys.time()

for (i in 1:(n_tissues - 1)) {
  for (j in (i + 1):n_tissues) {
    
    counter <- counter + 1
    tissue1 <- all_tissues[i]
    tissue2 <- all_tissues[j]
    
    # Estimate time remaining
    if (counter > 1) {
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      avg_time <- elapsed / (counter - 1)
      remaining <- avg_time * (n_comparisons - counter + 1)
      eta <- format(Sys.time() + remaining, "%H:%M:%S")
      cat(sprintf("[%3d/%3d] %s vs %s (ETA: %s)... ", 
                  counter, n_comparisons, tissue1, tissue2, eta))
    } else {
      cat(sprintf("[%3d/%3d] %s vs %s... ", 
                  counter, n_comparisons, tissue1, tissue2))
    }
    
    tryCatch({
      # Prepare data
      mcmc_list <- list(
        create_mcmc_list_format(mcmc_data_human_prepared[[tissue1]]),
        create_mcmc_list_format(mcmc_data_human_prepared[[tissue2]])
      )
      
      # Run analysis
      result <- multi_conservation(
        mcmc.merge.list = mcmc_list,
        dataset.names = c(tissue1, tissue2),
        select.pathway.list = "global",
        n_perm = 1000,
        n_boot = 1000,
        output.dir = tempdir(),
        use_cpp = TRUE
      )
      
      # Extract results (assuming result is a data frame with one row)
      results_df <- rbind(results_df, data.frame(
        tissue1 = tissue1,
        tissue2 = tissue2,
        raw_concordance = result[1, 2],
        adjusted_concordance = result[1, 3],
        ci_lower = result[1, 6],
        ci_upper = result[1, 7],
        null_mean = result[1, 10],
        null_sd = result[1, 11],
        p_value = result[1, 8],
        z_score = result[1, 9],
        stringsAsFactors = FALSE
      ))
      
      cat("✓\n")
      
      # Save intermediate results every 25 comparisons
      if (counter %% 25 == 0) {
        saveRDS(results_df, file.path(output_dir, "intermediate_results.rds"))
      }
      
    }, error = function(e) {
      cat(sprintf("✗ ERROR: %s\n", e$message))
      
      results_df <- rbind(results_df, data.frame(
        tissue1 = tissue1,
        tissue2 = tissue2,
        raw_concordance = NA,
        adjusted_concordance = NA,
        ci_lower = NA,
        ci_upper = NA,
        null_mean = NA,
        null_sd = NA,
        p_value = NA,
        z_score = NA,
        stringsAsFactors = FALSE
      ))
    })
    
  }
}

total_time <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("\nTotal time: %.1f minutes\n", total_time))

################################################################################
# CLEAN AND SORT RESULTS
################################################################################

cat("\n")
cat("="  %>% rep(80) %>% paste0(collapse=""), "\n")
cat("PROCESSING RESULTS\n")
cat("="  %>% rep(80) %>% paste0(collapse=""), "\n")
cat("\n")

# Remove failed comparisons
results_clean <- results_df[complete.cases(results_df), ]

cat(sprintf("Successful comparisons: %d / %d\n", nrow(results_clean), nrow(results_df)))

if (nrow(results_clean) < nrow(results_df)) {
  failed <- results_df[!complete.cases(results_df), ]
  cat("\nFailed comparisons:\n")
  print(failed[, c("tissue1", "tissue2")])
}

# Sort by adjusted concordance
results_sorted <- results_clean %>%
  arrange(desc(adjusted_concordance))

################################################################################
# SAVE RESULTS
################################################################################

cat("\n")
cat("Saving results...\n")

# 1. Full results CSV
csv_file <- file.path(output_dir, "human_all_pairwise_results.csv")
write.csv(results_sorted, csv_file, row.names = FALSE)
cat(sprintf("  ✓ Saved: %s\n", csv_file))

# 2. RDS file for easy loading in R
rds_file <- file.path(output_dir, "human_all_pairwise_results.rds")
saveRDS(results_sorted, rds_file)
cat(sprintf("  ✓ Saved: %s\n", rds_file))

# 3. Summary statistics
summary_stats <- data.frame(
  metric = c("Raw Concordance", "", "", "",
             "Adjusted Concordance", "", "", "",
             "Null Mean", "", "", "",
             "P-value", "", "", ""),
  statistic = rep(c("Min", "Max", "Mean", "Median"), 4),
  value = c(
    min(results_clean$raw_concordance),
    max(results_clean$raw_concordance),
    mean(results_clean$raw_concordance),
    median(results_clean$raw_concordance),
    min(results_clean$adjusted_concordance),
    max(results_clean$adjusted_concordance),
    mean(results_clean$adjusted_concordance),
    median(results_clean$adjusted_concordance),
    min(results_clean$null_mean),
    max(results_clean$null_mean),
    mean(results_clean$null_mean),
    median(results_clean$null_mean),
    min(results_clean$p_value),
    max(results_clean$p_value),
    mean(results_clean$p_value),
    median(results_clean$p_value)
  )
)

summary_file <- file.path(output_dir, "summary_statistics.csv")
write.csv(summary_stats, summary_file, row.names = FALSE)
cat(sprintf("  ✓ Saved: %s\n", summary_file))

# 4. Top and bottom concordance pairs
n_top <- min(20, nrow(results_clean))

top_pairs <- results_sorted[1:n_top, c("tissue1", "tissue2", "adjusted_concordance", "p_value")]
top_file <- file.path(output_dir, "top_20_concordant_pairs.csv")
write.csv(top_pairs, top_file, row.names = FALSE)
cat(sprintf("  ✓ Saved: %s\n", top_file))

bottom_pairs <- tail(results_sorted, n_top)[, c("tissue1", "tissue2", "adjusted_concordance", "p_value")]
bottom_file <- file.path(output_dir, "bottom_20_concordant_pairs.csv")
write.csv(bottom_pairs, bottom_file, row.names = FALSE)
cat(sprintf("  ✓ Saved: %s\n", bottom_file))

################################################################################
# PRINT SUMMARY
################################################################################

cat("\n")
cat("="  %>% rep(80) %>% paste0(collapse=""), "\n")
cat("SUMMARY STATISTICS\n")
cat("="  %>% rep(80) %>% paste0(collapse=""), "\n")
cat("\n")

print(summary_stats)

cat("\n")
cat("TOP 10 MOST CONCORDANT TISSUE PAIRS:\n")
cat("-"  %>% rep(80) %>% paste0(collapse=""), "\n")
print(head(results_sorted[, c("tissue1", "tissue2", "adjusted_concordance", "p_value")], 20))

cat("\n")
cat("BOTTOM 10 LEAST CONCORDANT TISSUE PAIRS:\n")
cat("-"  %>% rep(80) %>% paste0(collapse=""), "\n")
print(tail(results_sorted[, c("tissue1", "tissue2", "adjusted_concordance", "p_value")], 10))

################################################################################
# VISUALIZATIONS
################################################################################

cat("\n")
cat("Creating visualizations...\n")

# 1. Histogram of adjusted concordance
p1 <- ggplot(results_clean, aes(x = adjusted_concordance)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = median(results_clean$adjusted_concordance),
             color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = mean(results_clean$adjusted_concordance),
             color = "blue", linetype = "dashed", size = 1) +
  labs(title = "Distribution of Adjusted Concordance (All human Tissue Pairs)",
       subtitle = sprintf("Red = Median (%.3f), Blue = Mean (%.3f)",
                          median(results_clean$adjusted_concordance),
                          mean(results_clean$adjusted_concordance)),
       x = "Adjusted Concordance Index",
       y = "Count") +
  theme_minimal()

ggsave(file.path(output_dir, "histogram_adjusted_concordance.png"), 
       p1, width = 10, height = 6, dpi = 300)
cat("  ✓ Saved histogram\n")

# 2. Raw vs Adjusted
p2 <- ggplot(results_clean, aes(x = raw_concordance, y = adjusted_concordance)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(title = "Raw vs Adjusted Concordance",
       x = "Raw Concordance",
       y = "Adjusted Concordance") +
  theme_minimal()

ggsave(file.path(output_dir, "raw_vs_adjusted.png"), 
       p2, width = 8, height = 6, dpi = 300)
cat("  ✓ Saved scatter plot\n")

# 3. Null mean distribution
p3 <- ggplot(results_clean, aes(x = null_mean)) +
  geom_histogram(bins = 30, fill = "coral", color = "black", alpha = 0.7) +
  geom_vline(xintercept = median(results_clean$null_mean),
             color = "red", linetype = "dashed", size = 1) +
  labs(title = "Distribution of Null Expectations",
       subtitle = sprintf("Median = %.3f", median(results_clean$null_mean)),
       x = "Null Mean",
       y = "Count") +
  theme_minimal()

ggsave(file.path(output_dir, "histogram_null_mean.png"), 
       p3, width = 10, height = 6, dpi = 300)
cat("  ✓ Saved null distribution plot\n")

# 4. Heatmap (if not too large)
if (n_tissues <= 30) {
  
  # Create concordance matrix
  concordance_matrix <- matrix(1, nrow = n_tissues, ncol = n_tissues,
                               dimnames = list(all_tissues, all_tissues))
  
  for (k in 1:nrow(results_clean)) {
    t1 <- results_clean$tissue1[k]
    t2 <- results_clean$tissue2[k]
    val <- results_clean$adjusted_concordance[k]
    concordance_matrix[t1, t2] <- val
    concordance_matrix[t2, t1] <- val
  }
  
  library(reshape2)
  melted <- melt(concordance_matrix)
  
  p4 <- ggplot(melted, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = median(results_clean$adjusted_concordance),
                         name = "Adjusted\nConcordance") +
    labs(title = "Pairwise Tissue Concordance Heatmap (human)",
         x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 8))
  
  ggsave(file.path(output_dir, "heatmap_concordance.png"), 
         p4, width = 12, height = 11, dpi = 300)
  cat("  ✓ Saved heatmap\n")
}

################################################################################
# FINAL SUMMARY
################################################################################

cat("\n")
cat("="  %>% rep(80) %>% paste0(collapse=""), "\n")
cat("ANALYSIS COMPLETE!\n")
cat("="  %>% rep(80) %>% paste0(collapse=""), "\n")
cat("\n")
cat("All results saved to:", output_dir, "\n")
cat("\nFiles created:\n")
cat("  1. human_all_pairwise_results.csv - Full results table\n")
cat("  2. human_all_pairwise_results.rds - R data file\n")
cat("  3. summary_statistics.csv - Summary stats\n")
cat("  4. top_20_concordant_pairs.csv - Highest concordance\n")
cat("  5. bottom_20_concordant_pairs.csv - Lowest concordance\n")
cat("  6. histogram_adjusted_concordance.png - Distribution plot\n")
cat("  7. raw_vs_adjusted.png - Scatter plot\n")
cat("  8. histogram_null_mean.png - Null distribution\n")
if (n_tissues <= 30) {
  cat("  9. heatmap_concordance.png - Heatmap visualization\n")
}
cat("\n")


results <- read.csv("human_pairwise_concordance/human_all_pairwise_results.csv")
head(results[order(-results$adjusted_concordance, results$p_value), ])
