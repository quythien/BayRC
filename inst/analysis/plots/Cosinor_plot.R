# Scatter plot 
current_wd <- "/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects"
load(file.path(current_wd, "Collaborative/GTEXdata/data/CAMO.bab.hum.RData"))
# Baboon Lung
library(dplyr)
library(tidyr)
library(purrr)
source("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Pipeline/one_cosinor_OLS_new.R")

# Plot for baboon
baboon_LUN = list(
  expr = baboon_withTOD$baboon$LUN,
  tod = baboon_withTOD$tod$LUN
)

# Plot for human 
human_LUN = list(
  expr = gtex$CPM.large.clean$LUN,
  tod = gtex$tod$LUN
)

###
# Load required libraries
library(biomaRt)
library(dplyr)

# ============================================
# 1. MAP HUMAN GENES TO SYMBOLS
# ============================================

# Get row names from human_LUN$expr (assuming these are Ensembl IDs)
human_genes <- rownames(human_LUN$expr)

# Connect to Ensembl
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "hsapiens_gene_ensembl",
                      mirror = "useast")

# Query for gene symbols
human_gene_map <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = human_genes,
  mart = ensembl
)

# Create a complete mapping with fallback to original ID
# Start with a data frame of all original IDs
complete_map <- data.frame(
  ensembl_gene_id = human_genes,
  stringsAsFactors = FALSE
)

# Merge with biomaRt results
complete_map <- complete_map %>%
  left_join(human_gene_map, by = "ensembl_gene_id") %>%
  mutate(
    # Use symbol if available and not empty, otherwise use original ID
    final_symbol = ifelse(is.na(hgnc_symbol) | hgnc_symbol == "", 
                          ensembl_gene_id, 
                          hgnc_symbol)
  )
# ============================================
# 2. STANDARDIZE BMAL1/ARNTL NAMING
# ============================================

# Function to standardize gene names
standardize_gene_names <- function(symbols) {
  symbols <- toupper(symbols)  # Convert to uppercase for consistency
  
  # ARNTL (official symbol) -> BMAL1 (common name)
  symbols[symbols == "ARNTL"] <- "BMAL1"
  
  # You might also want to handle other aliases:
  # symbols[symbols == "MOP3"] <- "BMAL1"  # Another ARNTL alias
  # symbols[symbols == "BHLHE5"] <- "BMAL1"  # Yet another alias
  
  return(symbols)
}
complete_map$final_symbol_std <- standardize_gene_names(complete_map$final_symbol)

symbol_map <- setNames(complete_map$final_symbol_std, 
                       complete_map$ensembl_gene_id)

# Add symbol column to human data
human_LUN$expr$Symbol <- symbol_map[rownames(human_LUN$expr)]

# Reorganize to put Symbol as first column
human_LUN$expr <- human_LUN$expr %>%
  dplyr::select(Symbol, everything())

baboon_LUN$expr$Symbol_std <- standardize_gene_names(baboon_LUN$expr$Symbol)


# Plot these selected genes:
library(ggplot2)
library(gridExtra)

# ============================================
# PREPROCESS BABOON DATA (TPM + LOG only)
# ============================================

bab_raw <- baboon_LUN$expr
bab_mat <- bab_raw[, grep("LUN.ZT", colnames(bab_raw)), drop = FALSE]
tod_bab <- as.numeric(sub("LUN\\.ZT", "", colnames(bab_mat)))
gene_symbols <- bab_raw$Symbol

bab_mat <- as.matrix(bab_mat)
rownames(bab_mat) <- gene_symbols

# TPM normalization + log2 transform
tpm_bab <- sweep(bab_mat, 2, colSums(bab_mat), FUN = "/") * 1e6
logtpm_bab <- log2(tpm_bab + 1)

baboon_LUN_processed <- list(
  expr = logtpm_bab,
  tod = tod_bab
)

# ============================================
# CONVERT TIME SCALE
# ============================================

convert_zt <- function(tod) {
  ifelse(tod > 18, tod - 24, tod)
}

human_tod_converted <- convert_zt(human_LUN$tod)
baboon_tod_converted <- convert_zt(baboon_LUN_processed$tod)

# ============================================
# DEFINE GENES AND CREATE PLOTS
# ============================================

genes_to_plot <- c("CRY1", "DBP", "CRY2", "NR1D1", "PER2", "BMAL1", "NR1D2", "MRPL43", "MRPL1")
plots_list <- list()

for(gene in genes_to_plot) {
  
  # Find gene indices
  human_idx <- which(human_LUN$expr$Symbol == gene)
  baboon_idx <- which(rownames(baboon_LUN_processed$expr) == gene)
  
  if(length(human_idx) == 0 || length(baboon_idx) == 0) next
  if(length(human_idx) > 1) human_idx <- human_idx[1]
  if(length(baboon_idx) > 1) baboon_idx <- baboon_idx[1]
  
  # Extract expression
  human_expr <- as.numeric(human_LUN$expr[human_idx, -1])
  baboon_expr <- as.numeric(baboon_LUN_processed$expr[baboon_idx, ])
  
  # Fit cosinor
  human_fit <- tryCatch(
    one_cosinor_OLS(tod = human_tod_converted, y = human_expr, 
                    alpha = 0.05, period = 24, CI = FALSE),
    error = function(e) NULL
  )
  
  baboon_fit <- tryCatch(
    one_cosinor_OLS(tod = baboon_tod_converted, y = baboon_expr, 
                    alpha = 0.05, period = 24, CI = FALSE),
    error = function(e) NULL
  )
  
  if(is.null(human_fit) || is.null(baboon_fit)) next
  
  # Prepare data
  combined_data <- rbind(
    data.frame(tod = human_tod_converted, expression = human_expr, species = "Human"),
    data.frame(tod = baboon_tod_converted, expression = baboon_expr, species = "Baboon")
  )
  
  # Generate fitted curves
  time_seq <- seq(min(c(human_tod_converted, baboon_tod_converted)), 
                  max(c(human_tod_converted, baboon_tod_converted)), 
                  length.out = 200)
  
  fitted_combined <- rbind(
    data.frame(
      tod = time_seq,
      expression = human_fit$M$est + human_fit$A$est * cos(2*pi*time_seq/24 + human_fit$phase$est),
      species = "Human"
    ),
    data.frame(
      tod = time_seq,
      expression = baboon_fit$M$est + baboon_fit$A$est * cos(2*pi*time_seq/24 + baboon_fit$phase$est),
      species = "Baboon"
    )
  )
  
  # Convert peaks
  human_peak <- ifelse(human_fit$peak > 18, human_fit$peak - 24, human_fit$peak)
  baboon_peak <- ifelse(baboon_fit$peak > 18, baboon_fit$peak - 24, baboon_fit$peak)
  
  # Create plot
  p <- ggplot() +
    geom_point(data = combined_data, 
               aes(x = tod, y = expression, color = species),
               size = 2.5, alpha = 0.6) +
    geom_line(data = fitted_combined,
              aes(x = tod, y = expression, color = species),
              linewidth = 1) +
    scale_color_manual(values = c("Human" = "#E41A1C", "Baboon" = "#377EB8")) +
    labs(title = gene,
         subtitle = sprintf("Human: Peak=%.1fh, R²=%.2f, p=%.2e | Baboon: Peak=%.1fh, R²=%.2f, p=%.2e",
                            human_peak, human_fit$test$R2, human_fit$test$pval,
                            baboon_peak, baboon_fit$test$R2, baboon_fit$test$pval),
         x = "Zeitgeber Time (hours)",
         y = "Expression Level",
         color = "Species") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 8),
          axis.title = element_text(size = 10))
  
  plots_list[[gene]] <- p
}

# ============================================
# SAVE PDF
# ============================================

pdf("gene_cosinor_plots.pdf", width = 12, height = 10)

# Circadian genes
circadian_genes <- c("CRY1", "DBP", "CRY2", "NR1D1", "PER2", "BMAL1", "NR1D2")
circadian_plots <- plots_list[names(plots_list) %in% circadian_genes]
if(length(circadian_plots) > 0) {
  do.call(grid.arrange, c(circadian_plots, ncol = 3))
}

# Mitochondrial genes
mito_genes <- c("MRPL43", "MRPL1")
mito_plots <- plots_list[names(plots_list) %in% mito_genes]
if(length(mito_plots) > 0) {
  do.call(grid.arrange, c(mito_plots, ncol = 2))
}

dev.off()

cat("✓ Saved plots to gene_cosinor_plots.pdf\n")