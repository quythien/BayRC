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
base_result_dir <- file.path(current_aging, "results", "brain_regions")
# Expects final_results_extended in the session (final_brain_circadian_results_extended.RDS
# under base_result_dir, written by CAMO_Aging_Congruence_1.R)
options(width = Sys.getenv("COLUMNS"))
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

load(file.path(BAYRC_PATHWAY_DIR, "kegg.pathway.list_hsa.RData"))
load(file.path(BAYRC_PATHWAY_DIR, "kegg.pathway.list_cel_GeneNames.RData"))
load(file.path(BAYRC_PATHWAY_DIR, "hw_orth.RData"))
load(file.path(BAYRC_PATHWAY_DIR, "human.pathway.list.RData"))
load(file.path(BAYRC_PATHWAY_DIR, "go.pathway.list_hsa.RData"))

#── Cosinor benchmark ────────────────────────────────────────────────────────────
source(bayrc_file(BAYRC_PIPELINE_DIR, "one_cosinor_OLS_new.R"))

#── Source R scripts ─────────────────────────────────────────────────────────────
WD <- dirname(BAYRC_PACKAGE_DIR)
setwd(WD)
scripts <- list.files("R", full.names=TRUE)
sapply(scripts, source)

#────────────────────────────────────────────────────────────────
### Aim 1A: Pathway-level concordance 
all_equal_rows <- function(...) {
  matrices <- list(...)
  base_names <- rownames(matrices[[1]])
  all(sapply(matrices, function(m) identical(rownames(m), base_names)))
}

# Both groups must share the same gene order
all_equal_rows(
  final_results_extended$brain_results$COMBINED_younger$rho,
  final_results_extended$brain_results$COMBINED_older$rho
) # TRUE

# Gene symbols from Ensembl IDs, needed for pathway enrichment
library(biomaRt)
Sys.setenv(XDG_CACHE_HOME = "~/biocache")
biomartCacheClear()
ensemble <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
younger_symbol  <- match_symbols(final_results_extended$brain_results$COMBINED_younger, 
                                 BF = 3, 
                                 p_rhythmic = 0.2, 
                                 ensemble)

older_symbol <- match_symbols(final_results_extended$brain_results$COMBINED_older, 
                              BF = 3, 
                              p_rhythmic = 0.2, 
                              ensemble)

#################
# Discordance

discordance(final_results_extended$brain_results$BA11_older$rho, final_results_extended$brain_results$BA47_older$rho)

# Directional discordance, each adjusted against its chance level
discordance <- function(rho1, rho2) {
  A  <- mean(rho1 == 1)                     # P(rho1=1)
  D  <- mean(rho2 == 1)                     # P(rho2=1)
  TP <- mean(rho1 == 1 & rho2 == 1)
  FP <- mean(rho1 == 1 & rho2 == 0)
  FN <- mean(rho1 == 0 & rho2 == 1)

  # A → B: proportion of A's rhythmic genes that are not rhythmic in B
  discord_A2B <- if (A > 0) FP / A else NA_real_
  discord_A2B_null <- 1 - D
  discord_A2B_adj <- if (!is.na(discord_A2B)) {
    (discord_A2B - discord_A2B_null) / (1 - discord_A2B_null)
  } else NA_real_

  # B → A: proportion of B's rhythmic genes that are not rhythmic in A
  discord_B2A <- if (D > 0) FN / D else NA_real_
  discord_B2A_null <- 1 - A
  discord_B2A_adj <- if (!is.na(discord_B2A)) {
    (discord_B2A - discord_B2A_null) / (1 - discord_B2A_null)
  } else NA_real_

  return(list(
    discord_A2B = discord_A2B,
    discord_A2B_adj = discord_A2B_adj,
    discord_B2A = discord_B2A,
    discord_B2A_adj = discord_B2A_adj
  ))
}

#################

data_rho <- list(
  younger = younger_symbol$rho,
  older  = older_symbol$rho
)

data_phi <- list(
  younger = younger_symbol$phi,
  older  = older_symbol$phi
)

# Stage 1: union test for pathways with rhythmically active genes
select.pathway.kegg_union <- pathSelect(
  mcmc.merge.list      = list(younger = younger_symbol, older = older_symbol),
  pathway.list         = kegg.pathway.list_hsa,
  dataset.names        = c("younger", "older"),
  ranking.method       = "union",
  score_type           = "pos",
  pathwaysize.lower.cut = 5,
  pathwaysize.upper.cut = 200,
  qvalue.cut           = 0.20,
  nperm                = 1000
)

# Active pathways (p < 0.05 at Stage 1)
active_pathway_idx   <- select.pathway.kegg_union$results$pval < 0.05
active_pathway_list  <- kegg.pathway.list_hsa[
  select.pathway.kegg_union$results$pathway[active_pathway_idx]
]

# Stage 2: gain/loss enrichment within active pathways
select.pathway.kegg_gain <- pathSelect(
  mcmc.merge.list      = list(younger = younger_symbol, older = older_symbol),
  pathway.list         = active_pathway_list,
  dataset.names        = c("younger", "older"),
  ranking.method       = "gain",
  score_type           = "pos",
  pathwaysize.lower.cut = 5,
  pathwaysize.upper.cut = 200,
  qvalue.cut           = 0.20,
  nperm                = 1000
)

select.pathway.kegg_loss <- pathSelect(
  mcmc.merge.list      = list(younger = younger_symbol, older = older_symbol),
  pathway.list         = active_pathway_list,
  dataset.names        = c("younger", "older"),
  ranking.method       = "loss",
  score_type           = "pos",
  pathwaysize.lower.cut = 5,
  pathwaysize.upper.cut = 200,
  qvalue.cut           = 0.20,
  nperm                = 1000
)

#──────────────────────────────────────────────────────────────────
# Aim 1A: Pathway concordance score 
library(parallel)

pairwise_names <- combn(names(data_rho), 2, simplify = FALSE)

# Pathways passed through the Stage 1 union test
all_pathway_names <- select.pathway.kegg_union$results$pathway

all_acs_results <- list()

output_base <- "acs_output"

for (pair in pairwise_names) {
  pair_label <- paste(pair, collapse = "_vs_")
  cat("=== Processing:", pair_label, "===\n")

  dataset <- list(
    data_rho[[pair[1]]],
    data_rho[[pair[2]]]
  )
  names(dataset) <- pair

  # ACS for each pathway
  acs_results <- mclapply(all_pathway_names, function(pw_name) {
    cat("  Pathway:", pw_name, "\n")
    
    single_pathway <- kegg.pathway.list_hsa[pw_name]
    
    result <- tryCatch({
      out <- multi_ACS_ADS_pathway(
        mcmc.merge.list = dataset,
        dataset.names = names(dataset),
        select.pathway.list = single_pathway,
        measure = "Fmeasure",
        B = 10,
        parallel = TRUE,
        n.cores = 20,
        output.dir = file.path(output_base, pair_label, "by_pathway", pw_name)
      )
      
      acs_val <- out$ACS.mat[1, 1]
      return(c(pathway = pw_name, acs = as.numeric(acs_val)))
      
    }, error = function(e) {
      cat("  Error in pathway", pw_name, ":", conditionMessage(e), "\n")
      return(NULL)
    })
    
  }, mc.cores = min(length(all_pathway_names), 30))

  acs_results_clean <- do.call(rbind, acs_results[!sapply(acs_results, is.null)])
  
  if (!is.null(acs_results_clean)) {
    acs_df <- data.frame(
      pathway = acs_results_clean[, "pathway"],
      acs = as.numeric(acs_results_clean[, "acs"]),
      stringsAsFactors = FALSE
    )
    acs_df_sorted <- acs_df[order(-acs_df$acs), ]
    rownames(acs_df_sorted) <- NULL
    all_acs_results[[pair_label]] <- acs_df_sorted
  } else {
    warning("No valid results for:", pair_label)
  }
}

library(openxlsx)

wb <- createWorkbook()
for (sheet in names(all_acs_results)) {
  addWorksheet(wb, sheetName = sheet)
  writeData(wb, sheet = sheet, all_acs_results[[sheet]])
}

saveWorkbook(wb, file = file.path(current_aging, "output/Pairwise_ACS_Pathway_GSEA.xlsx"), overwrite = TRUE)

#──────────────────────────────────────────────────────────────────
# Phase inference 
shift_phi_to_positive <- function(phi_matrix) {
  phi_matrix + 6  # shift [-6, 18] to [0, 24]
}

data_phi <- list(
  younger = shift_phi_to_positive(data_phi$younger),
  older = shift_phi_to_positive(data_phi$older)
)

data_phi <- list(
  younger = data_phi$younger %% 24,
  older = data_phi$older %% 24
)

P <- 24
comparisons <- combn(names(data_phi), 2, simplify = FALSE)

results_list <- mclapply(comparisons, function(pair) {
  pair_label <- paste(pair, collapse = "_vs_")
  cat("Processing:", pair_label, "\n")

  matrix1 <- list(
    phi = data_phi[[pair[1]]],
    rho = data_rho[[pair[1]]]
  )
  
  matrix2 <- list(
    phi = data_phi[[pair[2]]],
    rho = data_rho[[pair[2]]]
  )
  
  res <- phase_inf(matrix1, matrix2,
                   P = 24, credMass = 0.90,
                   shift = 4, a = -12,
                   fdr_thresh = 0.20)

  list(
    label     = pair_label,
    diff_df   = res$phase_difference,
    conserve_df = res$phase_conservation
  )
}, mc.cores = min(length(comparisons), detectCores() - 1))

names(results_list) <- sapply(results_list, function(x) x$label)

# Phase differences and phase conservation, one workbook each
library(openxlsx)
wb_diff <- createWorkbook()
wb_conserve <- createWorkbook()

for (res in results_list) {
  addWorksheet(wb_diff,     sheetName = res$label)
  addWorksheet(wb_conserve, sheetName = res$label)
  
  writeData(wb_diff,     sheet = res$label, res$diff_df)
  writeData(wb_conserve, sheet = res$label, res$conserve_df)
}

saveWorkbook(wb_diff,     file = file.path(current_aging, "results/brain_regions/phase_shift_diff_FDR_4h.xlsx"),     overwrite = TRUE)
saveWorkbook(wb_conserve, file = file.path(current_aging, "results/brain_regions/phase_shift_conserve_FDR_4h.xlsx"), overwrite = TRUE)

####################################
# AIM 1C

options(bitmapType = "cairo")
library(parallel)
library(pathview)

data_rho <- list(
  younger = younger_symbol$rho,
  older  = older_symbol$rho
)

data_phi <- list(
  younger = younger_symbol$phi,
  older  = older_symbol$phi
)

pathway_ids <- c(
  "Staphylococcus aureus infection" = "hsa05150",
  "Circadian rhythm" = "hsa04710",
  "Cell adhesion molecules (CAMs)" = "hsa04514",
  "Antigen processing and presentation" = "hsa04612",
  "Systemic lupus erythematosus" = "hsa05322",
  "Graft-versus-host disease" = "hsa05332",
  "Phagosome" = "hsa04145",
  "Rheumatoid arthritis" = "hsa05323",
  "Allograft rejection" = "hsa05330",
  "Type I diabetes mellitus" = "hsa04940",
  "Autoimmune thyroid disease" = "hsa05320",
  "Osteoclast differentiation" = "hsa04380",
  "Hematopoietic cell lineage" = "hsa04640",
  "Jak-STAT signaling pathway" = "hsa04630",
  "Asthma" = "hsa05310",
  "p53 signaling pathway" = "hsa04115",
  "Complement and coagulation cascades" = "hsa04610",
  "Viral myocarditis" = "hsa05416",
  "Apoptosis" = "hsa04210"
)

pathway_ids <- c(
  "Circadian rhythm" = "hsa04710",
  "Cell adhesion molecules (CAMs)" = "hsa04514",
  "Jak-STAT signaling pathway" = "hsa04630",
  "p53 signaling pathway" = "hsa04115",
  "Complement and coagulation cascades" = "hsa04610",
  "Apoptosis" = "hsa04210",
  "Hematopoietic cell lineage" = "hsa04640"
)

pairs <- list(
  young_old  = c("younger", "older")
)

tasks <- expand.grid(
  gene_type = "all", 
  pathway   = names(pathway_ids),
  contrast  = names(pairs),
  stringsAsFactors = FALSE
)
tasks$pw_id <- sub("hsa", "", pathway_ids[tasks$pathway])

out_base <- file.path(base_result_dir, "out_new")

plot_base <- file.path(base_result_dir, "module_plot_new")
if (!dir.exists(plot_base)) {
  dir.create(plot_base, recursive = TRUE)
}

if (!dir.exists(out_base)) {
  dir.create(out_base, recursive = TRUE)
}

for (i in seq_len(nrow(tasks))) {
  
  gt <- tasks$gene_type[i]
  pw <- tasks$pw_id[i]
  sel <- pairs[[tasks$contrast[i]]]

  out_dir  <- file.path(out_base, gt, tasks$contrast[i], tasks$pw_id[i])
  plot_dir <- file.path(plot_base, gt, tasks$contrast[i], tasks$pw_id[i])
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  result <- tryCatch({
    KEGG_module(
      mcmc.merge.list = data_rho[sel],
      dataset.names   = sel,
      KEGGspecies     = "hsa",
      KEGGpathwayID   = pw,
      data.pair       = sel,
      gene_type       = gt,
      bayesF_cut      = 2,
      p_rhythmic      = 0.2,
      minM            = 2,
      maxM            = NULL,
      B               = 1000,
      cores           = 20,
      search_method   = "SA",
      Elbow_plot      = TRUE,
      filePath        = out_dir,
      ensemble        = ensemble
    )
  }, error = function(e) {
    message(sprintf("✗ Error in task %d (%s, %s, %s): %s", 
                    i, gt, tasks$contrast[i], tasks$pw_id[i], e$message))
    return(NULL)
  })

  if (!is.null(result)) {
    tryCatch({
      KEGG_module_topology_plot(
        res_KEGG_module = result,
        which_to_draw   = "all",
        filePath        = plot_dir
      )
    }, error = function(e) {
      message(sprintf("✗ Failed to plot topology for task %d: %s", i, e$message))
    })
  }
}

#####
data_rho$older
# Density of per-gene posterior rhythmicity probability, older group
png(file.path(current_aging, "results/brain_regions/rho_rowmeans_density.png"), width = 800, height = 600)

plot(
  density(rowMeans(data_rho$older, na.rm = TRUE)),
  main = "Density of Row Means (Older)",
  xlab = "Row mean of rho",
  ylab = "Density",
  lwd = 2,
  col = "blue"
)

dev.off()

