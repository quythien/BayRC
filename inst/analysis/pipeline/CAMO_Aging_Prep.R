setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Collaborative/Paper/Congruence/PNAS_aging")
options(width = "140")

# Reading in 
library(oligo)
pheno <- read.table("data/GSE71620_Phenotype_GEO.txt", 
                             header = TRUE, 
                             sep = "\t", 
                             stringsAsFactors = FALSE)
cel_files <- list.files("data/GSE71620_extracted", 
                        pattern = "CEL.gz", 
                        full.names = TRUE)

data <- read.celfiles(cel_files)

# Checking data 
class(data)
dim(data)         # 1178100 features × 420 samples
sampleNames(data) # list of CEL files / sample IDs
featureNames(data)[1:10] # first 10 probe IDs
pData(data)         # shows the phenoData table
varLabels(phenoData(data))  # variable names

# Normalization 
#  Background correctionn
#  Quantile normalization
#  Log2 transformation

eset <- rlsma(data)        # ExpressionSet object
print(paste("Normalized expression matrix:", nrow(eset), "genes ×", ncol(eset), "samples"))
#[1] "Normalized expression matrix: 33297 genes × 420 samples" 

library(hugene11sttranscriptcluster.db)
probe_ids <- rownames(eset)  # Your probe IDs
gene_symbols <- mapIds(hugene11sttranscriptcluster.db, 
                       keys = probe_ids,
                       column = "SYMBOL", 
                       keytype = "PROBEID")

annotated_genes <- !is.na(gene_symbols)
eset_annotated <- eset[annotated_genes, ]
gene_names_annotated <- gene_symbols[annotated_genes]

print(paste("Original:", nrow(eset), "probe sets"))  # Original: 33297 probe sets
print(paste("Annotated:", nrow(eset_annotated), "genes with symbols")) # 22025 genes with symbols

# View(exprs(eset_annotated))
library(dplyr)
expr_data <- exprs(eset_annotated)
gene_annotation <- data.frame(
  ProbeID = rownames(expr_data),
  GeneSymbol = gene_names_annotated,
  stringsAsFactors = FALSE
)

# Calculate IQR for each probe
gene_annotation$IQR <- apply(expr_data, 1, IQR)

# Keep only the probe with highest IQR per gene
final_data <- gene_annotation %>%
  group_by(GeneSymbol) %>%
  slice_max(IQR, n = 1, with_ties = FALSE) %>%
  ungroup()

final_expr <- expr_data[final_data$ProbeID, ]
rownames(final_expr) <- final_data$GeneSymbol
nrow(final_expr) # Final set: 19998
# paper: 20237

############################################################################################
# Real data: final_expr (20252 x 420)
# pheno (210 x 8) : 210 subjects 
# PMI: mean postmortem interval 
# RIN: RNA integrity number  
# We have both BA11: Brodmann area 11 and BA47: Brodmann area 47 
# TOD: in range of -6 to 18 
# Split the data:
# Sixty-four subjects were removed based on these criteria: death was not witnessed because their time of death (TOD) cannot be precisely determined
pheno <- (pheno[!is.na(pheno$TOD), ]) # 146 patients 
pheno$AgeGroup <- ifelse(pheno$Age < 40, "younger",
                         ifelse(pheno$Age >= 60, "older", NA))
table(pheno$AgeGroup) # 37 older and 31 younger, match paper 

col_ids <- sub(".*BA(11|47)-([0-9]+)\\.CEL\\.gz", "\\2", colnames(final_expr))
keep_cols <- col_ids %in% pheno$ID
expr_filtered <- final_expr[, keep_cols]
# Filtered data that is also split 
BA11_expr <- expr_filtered[, grepl("BA11", colnames(expr_filtered))]
BA47_expr <- expr_filtered[, grepl("BA47", colnames(expr_filtered))]

saveRDS(list(expr = BA11_expr, pheno = pheno), "data/BA11_data.rds")
saveRDS(list(expr = BA47_expr, pheno = pheno), "data/BA47_data.rds")

############################################################################################
# Cosinor Analysis: All 
BA11 <- readRDS("data/BA11_data.rds")
BA47 <- readRDS("data/BA47_data.rds")

source("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/src/fitSinCurve.R")
source("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/src/circadianDrawing_axis.R")
source("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Pipeline/one_cosinor_OLS_new.R")

fit_cosinor <- function(df_expr, tod, region = "region", species_name = NULL) {
  # Keep complete cases on TOD
  keep <- is.finite(tod)
  df <- df_expr[, keep, drop = FALSE]  # Filter columns (samples)
  tod_use <- tod[keep]
  
  # Get gene names (genes are rows)
  genes <- rownames(df_expr)
  
  cosinor_list <- lapply(genes, function(gene) {
    # Extract gene expression values (genes are rows, samples are columns)
    gene_expression_values <- as.numeric(df[gene, ])
    
    # Fit cosinor model
    cosinor_result <- one_cosinor_OLS(tod = tod_use, y = gene_expression_values)
    
    # Create result data frame
    result_df <- data.frame(
      Gene = gene,
      Offset = cosinor_result$M$est,
      Offset_LL = cosinor_result$M$CI[1],
      Offset_UL = cosinor_result$M$CI[2],
      Amplitude = cosinor_result$A$est,
      Amplitude_SD = cosinor_result$A$sd,
      Phase = cosinor_result$phase$est,
      Phase_SD = cosinor_result$phase$sd,
      R2 = cosinor_result$test$R2,
      Fstat = cosinor_result$test$Fstat,
      sigma2 = cosinor_result$test$sigma2,
      Peak = cosinor_result$peak,
      P_Value = cosinor_result$test$pval,
      stringsAsFactors = FALSE
    )
    # Add optional columns
    if (!is.null(species_name)) {
      result_df$Species <- species_name
    }
    
    if (!is.null(region) && region != "region") {
      result_df$Region <- region
    }
    
    return(result_df)
  })
  
  # Combine all results
  do.call(rbind, cosinor_list)
}

circular_variance <- function(samples, P = 24) {
  samples_rad <- circular((samples / P) * (2 * pi), units = "radians")
  var_circ <- 1 - rho.circular(samples_rad)  # 1 - mean resultant length
  return(var_circ)
}

ba11_cosinor <- fit_cosinor(BA11$expr, BA11$pheno$TOD, region = "BA11")
ba47_cosinor <- fit_cosinor(BA47$expr, BA47$pheno$TOD, region = "BA47")
ba11_cosinor$Q_Value = p.adjust(ba11_cosinor$P_Value, method = "BH")
ba47_cosinor$Q_Value = p.adjust(ba47_cosinor$P_Value, method = "BH")

# Summary 
sum(ba11_cosinor$P_Value < 0.05) # 2486
sum(ba11_cosinor$P_Value < 0.01) # 816
sum(ba11_cosinor$P_Value < 0.001) # 163
sum(ba11_cosinor$P_Value < 0.0001) # 47 
sum(ba11_cosinor$Q_Value < 0.05) # 61
sum(ba11_cosinor$Q_Value < 0.01) # 21 

sum(ba47_cosinor$P_Value < 0.05) # 1610
sum(ba47_cosinor$P_Value < 0.01) # 509
sum(ba47_cosinor$P_Value < 0.001) # 119
sum(ba47_cosinor$P_Value < 0.0001) # 38 
sum(ba47_cosinor$Q_Value < 0.05) # 38 
sum(ba47_cosinor$Q_Value < 0.01) # 16 

# View 
View(ba11_cosinor[ba11_cosinor$Q_Value < 0.05, ])
View(ba47_cosinor[ba47_cosinor$Q_Value < 0.05, ])

# Export 
library(writexl)
write_xlsx(ba11_cosinor, "output/ba11_cosinor_full.xlsx")
write_xlsx(ba47_cosinor, "output/ba47_cosinor_full.xlsx")


############################################################################################
# AW-Fisher Meta Analysis 

common_genes <- intersect(ba11_cosinor$Gene, ba47_cosinor$Gene)

p_mat <- data.frame(
  Gene = common_genes,
  BA11 = ba11_cosinor$P_Value[match(common_genes, ba11_cosinor$Gene)],
  BA47 = ba47_cosinor$P_Value[match(common_genes, ba47_cosinor$Gene)]
)

library(AWFisher)

p_input <- as.matrix(p_mat[, -1])
rownames(p_input) = p_mat$Gene
# Apply AW-Fisher
awf_result <- AWFisher_pvalue(p_input)
qvalue <- p.adjust(awf_result$pvalue, "BH") 

# Extract results
sum(qvalue < 0.1) # 255 vs. paper : 465
sum(qvalue < 0.05) # 145 v. paper: 235
sum(qvalue < 0.01) # 63 vs. paper: 84

aw_result <- data.frame(
  Gene   = rownames(p_input),
  P_AW   = awf_result$pvalue,
  Q_AW   = qvalue,
  stringsAsFactors = FALSE
)

sig_genes <- aw_result %>%
  filter(Q_AW < 0.05) %>%
  arrange(P_AW) 

write_xlsx(aw_result, "output/meta_all.xlsx")
write_xlsx(sig_genes, "output/meta_q005.xlsx")

############################################################################################
# BA11
# Age stratification 
expr_ids <- sub(".*BA11-([0-9]+)\\.CEL\\.gz", "\\1", colnames(BA11$expr))
expr_ids <- as.integer(expr_ids)  # make numeric
names(expr_ids) <- colnames(BA11$expr)

# Younger IDs
younger_ids <- BA11$pheno$ID[BA11$pheno$AgeGroup == "younger"]

# Older IDs
older_ids   <- BA11$pheno$ID[BA11$pheno$AgeGroup == "older"]


# Expression subset for younger
BA11_younger <- BA11$expr[, expr_ids %in% younger_ids]

# Expression subset for older
BA11_older   <- BA11$expr[, expr_ids %in% older_ids]

# TOD
tod_map <- setNames(BA11$pheno$TOD, BA11$pheno$ID)

# TOD for each column in BA11_younger
ids_younger <- sub(".*BA11-([0-9]+)\\.CEL\\.gz", "\\1", colnames(BA11_younger))
tod_younger <- tod_map[ids_younger]

# TOD for each column in BA11_older
ids_older <- sub(".*BA11-([0-9]+)\\.CEL\\.gz", "\\1", colnames(BA11_older))
tod_older   <- tod_map[ids_older]


ba11_younger_cosinor <- fit_cosinor(BA11_younger, tod_younger, region = "BA11-younger")
ba11_older_cosinor <- fit_cosinor(BA11_older, tod_older, region = "BA11-older")
ba11_younger_cosinor$Q_Value = p.adjust(ba11_younger_cosinor$P_Value, method = "BH")
ba11_older_cosinor$Q_Value = p.adjust(ba11_older_cosinor$P_Value, method = "BH")


# Summary 
sum(ba11_younger_cosinor$P_Value < 0.05) # 1141
sum(ba11_younger_cosinor$P_Value < 0.01) # 268
sum(ba11_younger_cosinor$P_Value < 0.001) # 34
sum(ba11_younger_cosinor$P_Value < 0.0001) # 7 
sum(ba11_younger_cosinor$Q_Value < 0.05) # 2
sum(ba11_younger_cosinor$Q_Value < 0.01) # 2 

sum(ba11_older_cosinor$P_Value < 0.05) # 1274
sum(ba11_older_cosinor$P_Value < 0.01) # 265
sum(ba11_older_cosinor$P_Value < 0.001) # 23
sum(ba11_older_cosinor$P_Value < 0.0001) # 5 
sum(ba11_older_cosinor$Q_Value < 0.05) # 1 
sum(ba11_older_cosinor$Q_Value < 0.01) # 0 


write_xlsx(ba11_younger_cosinor, "output/ba11_cosinor_younger.xlsx")
write_xlsx(ba11_older_cosinor, "output/ba11_cosinor_older.xlsx")

############################################################################################
# BA47


expr_ids <- sub(".*BA47-([0-9]+)\\.CEL\\.gz", "\\1", colnames(BA47$expr))
expr_ids <- as.integer(expr_ids)  # make numeric
names(expr_ids) <- colnames(BA47$expr)

# Younger IDs
younger_ids <- BA47$pheno$ID[BA47$pheno$AgeGroup == "younger"]

# Older IDs
older_ids   <- BA47$pheno$ID[BA47$pheno$AgeGroup == "older"]


# Expression subset for younger
BA47_younger <- BA47$expr[, expr_ids %in% younger_ids]

# Expression subset for older
BA47_older   <- BA47$expr[, expr_ids %in% older_ids]

# TOD
tod_map <- setNames(BA47$pheno$TOD, BA47$pheno$ID)

# TOD for each column in BA47_younger
ids_younger <- sub(".*BA47-([0-9]+)\\.CEL\\.gz", "\\1", colnames(BA47_younger))
tod_younger <- tod_map[ids_younger]

# TOD for each column in BA47_older
ids_older <- sub(".*BA47-([0-9]+)\\.CEL\\.gz", "\\1", colnames(BA47_older))
tod_older   <- tod_map[ids_older]


ba47_younger_cosinor <- fit_cosinor(BA47_younger, tod_younger, region = "BA47-younger")
ba47_older_cosinor <- fit_cosinor(BA47_older, tod_older, region = "BA47-older")
ba47_younger_cosinor$Q_Value = p.adjust(ba47_younger_cosinor$P_Value, method = "BH")
ba47_older_cosinor$Q_Value = p.adjust(ba47_older_cosinor$P_Value, method = "BH")


# Summary 
sum(ba47_younger_cosinor$P_Value < 0.05) # 1881
sum(ba47_younger_cosinor$P_Value < 0.01) # 448
sum(ba47_younger_cosinor$P_Value < 0.001) # 79
sum(ba47_younger_cosinor$P_Value < 0.0001) # 15 
sum(ba47_younger_cosinor$Q_Value < 0.05) # 4
sum(ba47_younger_cosinor$Q_Value < 0.01) # 0 

sum(ba47_older_cosinor$P_Value < 0.05) # 1019
sum(ba47_older_cosinor$P_Value < 0.01) # 240
sum(ba47_older_cosinor$P_Value < 0.001) # 26
sum(ba47_older_cosinor$P_Value < 0.0001) # 8 
sum(ba47_older_cosinor$Q_Value < 0.05) # 0
sum(ba47_older_cosinor$Q_Value < 0.01) # 0 


write_xlsx(ba47_younger_cosinor, "output/ba47_cosinor_younger.xlsx")
write_xlsx(ba47_older_cosinor, "output/ba47_cosinor_older.xlsx")

# Check against paper: 
ba47_younger_cosinor[ba47_younger_cosinor$Gene == "FKBP5", ]$Peak    # 3.88 vs. paper 4  
ba47_older_cosinor[ba47_older_cosinor$Gene == "FKBP5", ]$Peak.  # 17.15 vs. paper 17

ba11_younger_cosinor[ba11_younger_cosinor$Gene == "FKBP5", ]$Peak # 4.33 vs paper 4
ba11_older_cosinor[ba11_older_cosinor$Gene == "FKBP5", ]$Peak # 17.35 vs paper 17 

ba47_cosinor[ba47_cosinor$Gene == "FKBP5", ]$Peak    # 4.09 vs. paper 4  
ba11_cosinor[ba11_cosinor$Gene == "FKBP5", ]$Peak    # 2.76 vs. paper 3  



############################################################################################
# Goal: Combining the two brain regions together -- stratified by age group 

# Cosinor Analysis: All 
BA11 <- readRDS("data/BA11_data.rds")
BA47 <- readRDS("data/BA47_data.rds")

source("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/src/fitSinCurve.R")
source("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/src/circadianDrawing_axis.R")
source("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Pipeline/one_cosinor_OLS_new.R")


# Step 1: Combine BA11 and BA47 expression data

# Find common genes between the two regions
common_genes <- intersect(rownames(BA11$expr), rownames(BA47$expr))
print(paste("Common genes between BA11 and BA47:", length(common_genes)))

# Subset both datasets to common genes only
BA11_common <- BA11$expr[common_genes, ]
BA47_common <- BA47$expr[common_genes, ]

# Combine expression matrices (genes x samples)
combined_expr <- cbind(BA11_common, BA47_common)
print(paste("Combined expression matrix:", nrow(combined_expr), "genes ×", ncol(combined_expr), "samples"))


ba11_ids <- sub(".*BA11-([0-9]+)\\.CEL\\.gz", "\\1", colnames(BA11_common))
ba47_ids <- sub(".*BA47-([0-9]+)\\.CEL\\.gz", "\\1", colnames(BA47_common))

combined_sample_info <- data.frame(
  sample_name = colnames(combined_expr),
  sample_id = as.integer(c(ba11_ids, ba47_ids)),
  region = c(rep("BA11", length(ba11_ids)), rep("BA47", length(ba47_ids))),
  stringsAsFactors = FALSE
)
pheno_map <- setNames(BA11$pheno$AgeGroup, BA11$pheno$ID)
tod_map <- setNames(BA11$pheno$TOD, BA11$pheno$ID)

combined_sample_info$age_group <- pheno_map[as.character(combined_sample_info$sample_id)]
combined_sample_info$TOD <- tod_map[as.character(combined_sample_info$sample_id)]

complete_samples <- !is.na(combined_sample_info$age_group) & !is.na(combined_sample_info$TOD)
combined_sample_info_clean <- combined_sample_info[complete_samples, ]
combined_expr_clean <- combined_expr[, complete_samples]

print("Sample distribution after combining regions:")
print(table(combined_sample_info_clean$age_group, combined_sample_info_clean$region))

# Step 2: Create age-stratified datasets combining both regions

# Younger group (combining BA11 + BA47)
younger_samples <- combined_sample_info_clean$age_group == "younger"
combined_younger_expr <- combined_expr_clean[, younger_samples]
combined_younger_tod <- combined_sample_info_clean$TOD[younger_samples]

print(paste("Combined younger group:", ncol(combined_younger_expr), "samples"))
print(paste("- BA11 younger samples:", sum(combined_sample_info_clean$region[younger_samples] == "BA11")))
print(paste("- BA47 younger samples:", sum(combined_sample_info_clean$region[younger_samples] == "BA47")))

# Older group (combining BA11 + BA47)
older_samples <- combined_sample_info_clean$age_group == "older"
combined_older_expr <- combined_expr_clean[, older_samples]
combined_older_tod <- combined_sample_info_clean$TOD[older_samples]

print(paste("Combined older group:", ncol(combined_older_expr), "samples"))
print(paste("- BA11 older samples:", sum(combined_sample_info_clean$region[older_samples] == "BA11")))
print(paste("- BA47 older samples:", sum(combined_sample_info_clean$region[older_samples] == "BA47")))

# Step 4: Perform cosinor analysis on combined age groups
# Younger group analysis (BA11 + BA47 combined)
combined_younger_cosinor <- fit_cosinor(combined_younger_expr, combined_younger_tod, region = "Combined-younger")
combined_younger_cosinor$Q_Value <- p.adjust(combined_younger_cosinor$P_Value, method = "BH")

# Older group analysis (BA11 + BA47 combined)
combined_older_cosinor <- fit_cosinor(combined_older_expr, combined_older_tod, region = "Combined-older")
combined_older_cosinor$Q_Value <- p.adjust(combined_older_cosinor$P_Value, method = "BH")

############################################################################################
# Summary statistics

print("=== COMBINED YOUNGER GROUP (BA11 + BA47) RESULTS ===")
print(paste("P < 0.05:", sum(combined_younger_cosinor$P_Value < 0.05))) # 2455
print(paste("P < 0.01:", sum(combined_younger_cosinor$P_Value < 0.01))) # 832
print(paste("P < 0.001:", sum(combined_younger_cosinor$P_Value < 0.001))) # 191
print(paste("P < 0.0001:", sum(combined_younger_cosinor$P_Value < 0.0001))) # 67
print(paste("Q < 0.05:", sum(combined_younger_cosinor$Q_Value < 0.05))) # 86
print(paste("Q < 0.01:", sum(combined_younger_cosinor$Q_Value < 0.01))) # 31

print("\n=== COMBINED OLDER GROUP (BA11 + BA47) RESULTS ===")
print(paste("P < 0.05:", sum(combined_older_cosinor$P_Value < 0.05))) # 2180
print(paste("P < 0.01:", sum(combined_older_cosinor$P_Value < 0.01))) # 667
print(paste("P < 0.001:", sum(combined_older_cosinor$P_Value < 0.001))) # 119
print(paste("P < 0.0001:", sum(combined_older_cosinor$P_Value < 0.0001))) # 34
print(paste("Q < 0.05:", sum(combined_older_cosinor$Q_Value < 0.05))) # 23 
print(paste("Q < 0.01:", sum(combined_older_cosinor$Q_Value < 0.01))) # 5


combined_pheno <- merge(combined_sample_info_clean, BA11$pheno, 
                        by.x = "sample_id", by.y = "ID", all.x = TRUE, sort = FALSE)

saveRDS(list(expr = combined_expr_clean, pheno = combined_pheno), "data/combined_data.rds")

print("Saved combined dataset: data/combined_data.rds")

write_xlsx(combined_younger_cosinor, "output/combined_cosinor_younger.xlsx")
write_xlsx(combined_older_cosinor, "output/combined_cosinor_older.xlsx")
write_xlsx(combined_sample_info_clean, "output/combined_sample_info.xlsx")


print("\n=== SIGNIFICANT GENES (Q < 0.05) ===")
sig_younger <- combined_younger_cosinor[combined_younger_cosinor$Q_Value < 0.05, ]
sig_older <- combined_older_cosinor[combined_older_cosinor$Q_Value < 0.05, ]
nrow(sig_younger)
nrow(sig_older)