# SAME SPECIE 


rm(list=ls());

current_gtex <- "/home/qtp1/Projects/Collaborative"
current_wd <- "/home/qtp1/Projects/Circadian"
current_aging <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging"

outdir = "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/all_plots"

load(file = "/home/qtp1/Projects/Collaborative/GTEXdata/result/summary/hb/mcmc_rho_BF3.RData")
load(file.path(current_gtex, "GTEXdata/result/summary/hb/phi/mcmc_phi_BF3.RData"))  




# "mcmc_data_baboon" "mcmc_data_human"  "mcmc_phi_baboon"  "mcmc_phi_human

library(KEGGREST)
library(parallel)


#base_result_dir <- file.path(current_aging, "results", "brain_regions")
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

thien_dir <- file.path(current_wd, "Kyle/Circadian-analysis-main/R/v1/BayRC/Thien")
source(file.path(thien_dir, "pathwaySelect.R"))
source(file.path(thien_dir, "multi_pathway.R"))
source(file.path(thien_dir, "multi_global.R"))
source(file.path(thien_dir, "internal.R"))
if (file.exists(file.path(thien_dir, "congruence.cpp"))) {
  Rcpp::sourceCpp(file.path(thien_dir, "congruence.cpp"))
}
if (file.exists(file.path(thien_dir, "permutation_functions.cpp"))) {
  Rcpp::sourceCpp(file.path(thien_dir, "permutation_functions.cpp"))
source(file.path(thien_dir, "plots/heatmap.R"))
}


#load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/kegg.pathway.list_hsa.RData"))
#load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/kegg.pathway.list_cel_GeneNames.RData"))
load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/hw_orth.RData"))
load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/human.pathway.list.RData"))
load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/go.pathway.list_hsa.RData"))

kegg.pathway.list_hsa <- readRDS("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/pathway_data/kegg_pathway_list_hsa.rds")

#── Source R scripts ─────────────────────────────────────────────────────────────
WD <- "Kyle/Circadian-analysis-main/R/v1"
setwd(file.path(current_wd, WD))
scripts <- list.files("BayRC", pattern="\\.R$", full.names = TRUE)
sapply(scripts, source)

# Reset wd
setwd(current_aging)

#---------------------------------------------------------------------------------
to_zt <- function(t_cos) ifelse(t_cos >= 18, t_cos - 24, t_cos)



baboon_LIV <- list(
  rho = mcmc_data_baboon$LIV,
  phi = mcmc_phi_baboon$LIV
)

baboon_ILE <- list(
  rho = mcmc_data_baboon$ILE,
  phi = mcmc_phi_baboon$ILE
)

#---------------------------------------------------------------------------------
# Global concordance score 
outLIV.dir <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/results/baboon/outLIV_final"
if (!dir.exists(outLIV.dir)) {
  dir.create(outLIV.dir, recursive = TRUE, showWarnings = FALSE)
}
cat("Directory ready:", outLIV.dir, "\n")
thien_dir <- file.path(current_wd, "Kyle/Circadian-analysis-main/R/v1/BayRC/Thien")
if (file.exists(file.path(thien_dir, "pathwaySelect.R"))) {
  source(file.path(thien_dir, "pathwaySelect.R"))
}
if (file.exists(file.path(thien_dir, "multi_pathway.R"))) {
  source(file.path(thien_dir, "multi_pathway.R"))
}
if (file.exists(file.path(thien_dir, "multi_global.R"))) {
  source(file.path(thien_dir, "multi_global.R"))
}
if (file.exists(file.path(thien_dir, "Permutation_Sim.R"))) {
  source(file.path(thien_dir, "Permutation_Sim.R"))
}
if (file.exists(file.path(thien_dir, "internal.R"))) {
  source(file.path(thien_dir, "internal.R"))
}
if (file.exists(file.path(thien_dir, "congruence.cpp"))) {
  Rcpp::sourceCpp(file.path(thien_dir, "congruence.cpp"))
}
if (file.exists(file.path(thien_dir, "permutation_functions.cpp"))) {
  Rcpp::sourceCpp(file.path(thien_dir, "permutation_functions.cpp"))
source(file.path(thien_dir, "plots/heatmap.R"))
}

# Analyze ALL genes together
t_start <- Sys.time()

results_global1 <- multi_conservation(
  mcmc.merge.list = list(baboon_LIV, baboon_ILE),
  dataset.names = c("baboon_LIV", "baboon_ILE"),
  select.pathway.list = "global",  
  n_perm = 500,
  n_boot = 100,
  outLIV.dir = file.path(outLIV.dir, "LIV_ILE_concordance"),
  use_cpp = TRUE,
  comLIVe_pvalue = TRUE, 
  comLIVe_ci = TRUE 
)

t_end <- Sys.time()

elapsed_time <- difftime(t_end, t_start, units = "mins")
print(elapsed_time)

# 2.43 minutes for p value only at 500 perm
#  1.03 min for boot only at 100 boots
#  0.8562239 for both 
rhy   <- detect_rhy(baboon_LIV, baboon_ILE, 0.25)
rhy_summary <- data.frame(
  Category = c("Rhythmic in LIV",
               "Rhythmic in ILE"),
  Count = c(rhy$n_rhythmic_A,
            rhy$n_rhythmic_B)
)

# 
# rhy_summary <- rbind(
#   rhy_summary,
#   data.frame(Category = c("LIV", "ILE"),
#              Count = c(rhy$threshold_A, rhy$threshold_B))
# )
# 
# rhy_summary
# 
#---------------------------------------------------------------------------------
pA    <- rowMeans(baboon_LIV$rho)
pB    <- rowMeans(baboon_ILE$rho)

# BFDR: Rhythmic biomarkers 
rhy   <- detect_rhy(baboon_LIV, baboon_ILE, 0.25)
rhy_summary <- data.frame(
  Category = c("Rhythmic in LIV",
               "Rhythmic in ILE"),
  Count = c(rhy$n_rhythmic_A,
            rhy$n_rhythmic_B)
)
rhy_summary

# BFDR: 
trans_outer <- transition_classify(pA, pB, bfdr_alpha = 0.25)

phase_inner <- phase_infer(
  phi_matrix1 = baboon_LIV$phi,
  phi_matrix2 = baboon_ILE$phi,
  gain_loss_status = trans_outer$gain_loss_status,
  bfdr_alpha = 0.25,
  shift = 2,
  P = 24, comLIVe_hdi = TRUE
)

# phase_inner2 <- phase_infer(
#   phi_matrix1 = baboon_LIV$phi,
#   phi_matrix2 = baboon_ILE$phi,
#   gain_loss_status = trans_outer$gain_loss_status,
#   bfdr_alpha = 0.2,
#   shift = 2,
#   P = 24, comLIVe_hdi = TRUE
# )

# === PHASE INFERENCE SUMMARY ===
#   Maintained genes: 677 
# BFDR α = 0.25  | shift threshold = 2 h
# Significant phase-shifted genes: 39 
# Significant phase-conserved genes: 335 
# Undetermined genes: 303 
# HDI comLIVation: enabled




#---------------------------------------------------------------------------------

# Maintained plot 
#───────────────────────────────────────────────────────────────
# Peak Concordance Plot for Maintained Genes - BABOON ILE vs LIV
#───────────────────────────────────────────────────────────────
library(ggplot2)
library(dplyr)
library(ggrepel)

# Helper
to_zt <- function(t_cos) ifelse(t_cos >= 18, t_cos - 24, t_cos)

#───────────────────────────────────────────────────────────────
# Extract data from phase_inner
#───────────────────────────────────────────────────────────────
gene_names <- names(phase_inner$peak1)

# Create dataframe from phase_inner
maintained_df <- data.frame(
  Gene = gene_names,
  Peak_LIV = phase_inner$peak1,
  Peak_ILE = phase_inner$peak2,
  deltaPhi = phase_inner$deltaPhi.Est,
  prob_shift = phase_inner$prob_shift,
  prob_conserved = phase_inner$prob_conserved,
  stringsAsFactors = FALSE
)

# Filter for maintained genes only
maintained_genes <- names(trans_outer$gain_loss_status[trans_outer$gain_loss_status == "Maintained"])

maintained_df <- maintained_df %>%
  filter(Gene %in% maintained_genes) %>%
  mutate(Peak_ILE_ZT = to_zt(Peak_ILE),
         Peak_LIV_ZT = to_zt(Peak_LIV))

#───────────────────────────────────────────────────────────────
# Assign phase classification
#───────────────────────────────────────────────────────────────
phase_class <- rep("Undetermined", length(gene_names))
names(phase_class) <- gene_names
phase_class[phase_inner$flag_cons]  <- "Phase-conserved"
phase_class[phase_inner$flag_shift] <- "Phase-shifted"
phase_class[phase_inner$flag_undetermined] <- "Undetermined"

# Merge into maintained_df
maintained_df <- maintained_df %>%
  mutate(phase_class = phase_class[Gene]) %>%
  mutate(phase_class = ifelse(is.na(phase_class), "Undetermined", phase_class))

#───────────────────────────────────────────────────────────────
# Circadian gene list
#───────────────────────────────────────────────────────────────
circadian_genes <- c(
  "BHLHE40","BHLHE41","BMAL1","BTRC","CLOCK","CREB1",
  "CRY1","CRY2","CSNK1D","CSNK1E","CUL1","DBP","FBXL3","FBXW11",
  "NFIL3","NPAS2","NR1D1","NR1D2","PER1","PER2","PER3",
  "RORA","RORB","RORC"
)
genes_to_label <- maintained_df %>% filter(Gene %in% circadian_genes)

#───────────────────────────────────────────────────────────────
# Color palette
#───────────────────────────────────────────────────────────────
phase_colors <- c(
  "Phase-conserved" = "#1B9E77",
  "Phase-shifted"   = "#D95F02",
  "Undetermined"    = "#7570B3"
)

#───────────────────────────────────────────────────────────────
# Global summary statistics
#───────────────────────────────────────────────────────────────
rhythmic_genes_outer <- names(trans_outer$gain_loss_status)[
  trans_outer$gain_loss_status %in% c("Gain", "Loss", "Maintained")
]

n_rhythmic_outer <- length(rhythmic_genes_outer)
n_phase_conserved <- sum(phase_inner$flag_cons, na.rm = TRUE)
n_phase_shifted <- sum(phase_inner$flag_shift, na.rm = TRUE)
n_undetermined <- sum(phase_inner$flag_undetermined, na.rm = TRUE)

pct_phase_conserved_global <- 100 * n_phase_conserved / n_rhythmic_outer
pct_phase_shifted_global <- 100 * n_phase_shifted / n_rhythmic_outer
pct_undetermined_global <- 100 * n_undetermined / n_rhythmic_outer

cat("\n[Global Summary - ILE vs LIV]\n")
cat("  Rhythmic (outer τ_c ≤ 0.25):", n_rhythmic_outer, "\n")
cat("  Phase-conserved (inner τ_p ≤ 0.10):", n_phase_conserved, 
    sprintf("(%.2f%%)\n", pct_phase_conserved_global))
cat("  Phase-shifted (inner τ_p ≤ 0.10):", n_phase_shifted,
    sprintf("(%.2f%%)\n", pct_phase_shifted_global))
cat("  Undetermined:", n_undetermined,
    sprintf("(%.2f%%)\n", pct_undetermined_global))

#───────────────────────────────────────────────────────────────
# ComLIVe ±3 h concordance among maintained genes
#───────────────────────────────────────────────────────────────
calculate_peak_difference <- function(a, b) {
  d <- abs(a - b)
  ifelse(d > 12, 24 - d, d)
}

maintained_df <- maintained_df %>%
  mutate(peak_diff = calculate_peak_difference(Peak_ILE_ZT, Peak_LIV_ZT),
         within_concordance = peak_diff <= 2)

n_total <- nrow(maintained_df)
n_within <- sum(maintained_df$within_concordance, na.rm = TRUE)
pct_within <- 100 * n_within / n_total

cat("\n[Within ±2 h summary]\n")
cat("  Maintained genes:", n_total, "\n")
cat("  Within ±2 h:", n_within, "\n")
cat(sprintf("  ⇒ %.1f%% within ±2 h interval\n", pct_within))

#───────────────────────────────────────────────────────────────
# Plot preparation with circular remapping to minimize distance to diagonal
#───────────────────────────────────────────────────────────────

# Function to re-map circular time points to minimize distance to diagonal
remap_to_diagonal <- function(x, y, offset = 24) {
  # Calculate all 9 possible transformations considering circular wrapping
  x_options <- c(x, x - offset, x + offset)
  y_options <- c(y, y - offset, y + offset)

  best_dist <- Inf
  best_x <- x
  best_y <- y

  for (x_opt in x_options) {
    for (y_opt in y_options) {
      # Distance to diagonal y = x is |y - x|
      dist <- abs(y_opt - x_opt)
      if (dist < best_dist - 1e-9) {
        best_dist <- dist
        best_x <- x_opt
        best_y <- y_opt
      }
    }
  }

  return(data.frame(x = best_x, y = best_y, dist = best_dist))
}

plot_df <- maintained_df %>%
  rowwise() %>%
  mutate(
    remapped = remap_to_diagonal(Peak_ILE_ZT, Peak_LIV_ZT),
    Peak_ILE_ZT_plot = remapped$x,
    Peak_LIV_ZT_plot = remapped$y
  ) %>%
  select(-remapped)

# Count categories
n_Rc <- nrow(plot_df)
pct_Rc <- round(100 * n_within / n_Rc, 1)

# Subtitle
subtitle_text <- bquote(
  "Rhythmically Conserved Set " ~ R[c] ~
    "(" * n == .(n_Rc) * ", " * .(pct_Rc) * "% within " * "\u00B1" * "2 h interval)"
)

#───────────────────────────────────────────────────────────────
# Create plot
#───────────────────────────────────────────────────────────────
p <- ggplot(plot_df, aes(
  x = Peak_ILE_ZT_plot,
  y = Peak_LIV_ZT_plot,
  color = phase_class
)) +
  # Identity and ±2 h boundaries (no need for ±22 due to remapping)
  geom_abline(intercept = 0, slope = 1,
              color = "black", linetype = "dashed", linewidth = 1.0) +
  geom_abline(intercept = 2, slope = 1,
              color = "darkgreen", linetype = "dotted", linewidth = 1.1) +
  geom_abline(intercept = -2, slope = 1,
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
    point.padding = 1.5,
    min.segment.length = 0,
    force_pull = 0.3,
    
    max.overlaps = Inf
  ) +
  
  labs(
    title = "Circadian Peak Concordance: Baboon ILE versus LIVpocampus",
    subtitle = subtitle_text,
    x = "Peak Hour – ILE (ZT)",
    y = "Peak Hour – LIVpocampus (ZT)",
    color = "Phase class"
  ) +
  
  scale_x_continuous(
    breaks = seq(-6, 18, 6),
    labels = sprintf("ZT%+d", seq(-6, 18, 6))
  ) +
  scale_y_continuous(
    breaks = seq(-6, 18, 6),
    labels = sprintf("ZT%+d", seq(-6, 18, 6))
  ) +
  coord_cartesian(xlim = c(-8, 20), ylim = c(-8, 20)) +

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
save_dir <- file.path(outLIV.dir, "figure/baboon_brain")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = file.path(save_dir, "Baboon_ILE_LIV_Peak_Concordance_0.25_2h_new.pdf"),
  plot = p,
  width = 9,
  height = 8
)



#───────────────────────────────────────────────────────────────
################################################################################
# COMPLETE WORKFLOW: Baboon ILE vs LIV Pathway Analysis
################################################################################

# Set parameters
dataset_names <- c("ILE", "LIV")
qvalue_cut <- 0.25
nperm <- 10000
pathway_size_min <- 10
pathway_size_max <- 300

# Set outLIV directory
# Set output directory
output.dir <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/results/baboon/output_final"
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
# STAGE 1: Identify Active Pathways (Union Test - P-value Filter)
################################################################################

cat("\n=== STAGE 1: Union Test (Active Pathways) ===\n")

result_union <- pathSelect(
  mcmc.merge.list = list(ILE = baboon_ILE, LIV = baboon_LIV),
  pathway.list = kegg.pathway.list_hsa,
  dataset.names = c("ILE", "LIV"),
  ranking.method = "union",
  score_type = "pos",
  qvalue.cut = 0.2,
  pathwaysize.lower.cut = pathway_size_min,
  pathwaysize.upper.cut = pathway_size_max,
  nperm = nperm,
  nproc = 1
)

# Extract the results dataframe
union_results <- result_union$results

# Filter by P-VALUE (not q-value)
active_pathways <- result_union$results %>%
  filter(pval < 0.20) %>%
  pull(pathway)


################################################################################
# STAGE 2: Test Transitions on Active Pathways (Q-value Filter)
################################################################################

cat("\n=== STAGE 2: Transition Tests (Gain, Loss, Conservation) ===\n")

# Ensure correct order
active_pathway_list <- kegg.pathway.list_hsa[match(active_pathways, names(kegg.pathway.list_hsa))]

# Gain enrichment - Q-VALUE FILTERING
result_gain <- pathSelect(
  mcmc.merge.list = list(ILE = baboon_ILE, LIV = baboon_LIV),
  pathway.list = kegg.pathway.list_hsa,
  dataset.names = c("ILE", "LIV"),
  ranking.method = "gain",
  score_type = "pos",
  qvalue.cut = 0.20,
  pathwaysize.lower.cut = 10,
  pathwaysize.upper.cut = 250,
  nperm = nperm,
  nproc = 1
)

# Loss enrichment - Q-VALUE FILTERING
result_loss <- pathSelect(
  mcmc.merge.list = list(ILE = baboon_ILE, LIV = baboon_LIV),
  pathway.list = active_pathway_list,
  dataset.names = c("ILE", "LIV"),
  ranking.method = "loss",
  score_type = "pos",
  qvalue.cut = 0.20,
  pathwaysize.lower.cut = 5,
  pathwaysize.upper.cut = 300,
  nperm = nperm,
  nproc = 1
)

# Conservation enrichment - Q-VALUE FILTERING
result_cons <- pathSelect(
  mcmc.merge.list = list(ILE = baboon_ILE, LIV = baboon_LIV),
  pathway.list = kegg.pathway.list_hsa,
  dataset.names = c("ILE", "LIV"),
  ranking.method = "conserved",
  score_type = "pos",
  qvalue.cut = 0.2,
  pathwaysize.lower.cut = 10,
  pathwaysize.upper.cut = pathway_size_max,
  nperm = nperm,
  nproc = 1
)


# result_gain$results %>%
#   filter(padj < 0.20) %>%
#   select(pathway, size, Gain_Index, Expected_N_Gain, pval, padj, Gain_Loss_Ratio_Arithmetic)
# 
# 
result_loss$results %>%
  filter(pval < 0.10) %>%
  select(pathway, size, Loss_Index, Expected_N_Loss, pval, padj, Gain_Loss_Ratio_Arithmetic)

# Add kegg phagosome in 

result_cons$results %>%
  filter(padj < 0.20) %>%
  select(pathway, size, Conserved_Index, Expected_N_Conserved, pval, padj, Gain_Loss_Ratio_Arithmetic)


cat("  Gain-enriched (LIV):  ", sum(result_gain$results$Significant), "\n")
cat("  Loss-enriched (LIV):  ", sum(result_loss$results$Significant), "\n")
cat("  Conservation-enriched:", sum(result_cons$results$Significant), "\n")

################################################################################
# STAGE 3: Filter to Significant Pathways
################################################################################

cat("\n=== STAGE 3: Filtering Significant Pathways ===\n")
library(tidyr)
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
  # filter(pathway %in% active_pathways) %>% # Chang this if needed 
  mutate(
    gain_sig = replace_na(gain_sig, FALSE),
    loss_sig = replace_na(loss_sig, FALSE),
    cons_sig = replace_na(cons_sig, FALSE)
  ) %>%
  arrange(pathway)

significant_pathways <- significance_table %>%
  filter(gain_sig | loss_sig | cons_sig) %>%
  pull(pathway) %>%
  union(c("KEGG Spliceosome", "KEGG Long-term depression", "KEGG Serotonergic synapse"))

cat("Pathways significant for at least one transition:", length(significant_pathways), "\n")

filtered_pathway_list <- kegg.pathway.list_hsa[match(significant_pathways, names(kegg.pathway.list_hsa))]

################################################################################
# STAGE 4: Run multi_conservation on FILTERED Pathways
################################################################################

cat("\n=== STAGE 4: Descriptive Metrics (multi_conservation) ===\n")
cat("Running on", length(significant_pathways), "filtered pathways...\n")

result_multiconservation <- multi_conservation(
  mcmc.merge.list = list(ILE = baboon_ILE, LIV = baboon_LIV),
  dataset.names = c("ILE", "LIV"),
  select.pathway.list = filtered_pathway_list,
  n_perm = 1000,
  n_boot = 1000,
  outLIV.dir = file.path(outLIV.dir, "multiconservation_ILE_LIV"),
  use_cpp = TRUE
)

################################################################################
# GO ENRICHMENT WORKFLOW (GO PATHWAY LIST)
################################################################################

cat("\n\n################################################################################\n")
cat("# GO ENRICHMENT WORKFLOW: Baboon Lung vs Human Lung\n")
cat("################################################################################\n")
load("/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/BayRC/Thien/pathway/go.pathway.list_hsa.RData")
output.dir <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/results/baboon/output_final"
go_output.dir <- file.path(output.dir, "go_enrichment_lung")
dir.create(go_output.dir, recursive = TRUE, showWarnings = FALSE)

result_gain_go <- pathSelect(
  mcmc.merge.list = list(Human_LUN = human_LUN, Baboon_LUN = baboon_LUN),
  pathway.list = go.pathway.list_hsa,
  dataset.names = c("Human_LUN", "Baboon_LUN"),
  ranking.method = "gain",
  score_type = "pos",
  qvalue.cut = 0.25,
  pathwaysize.lower.cut = 20,
  pathwaysize.upper.cut = 200,
  nperm = nperm,
  nproc = 1
)

top_gain_go <- result_gain_go$results %>%
  arrange(padj) %>%
  slice_head(n = 10) %>%
  dplyr::select(pathway, pval, padj, Top_Gain_Genes)

result_loss_go <- pathSelect(
  mcmc.merge.list = list(Human_LUN = human_LUN, Baboon_LUN = baboon_LUN),
  pathway.list = go.pathway.list_hsa,
  dataset.names = c("Human_LUN", "Baboon_LUN"),
  ranking.method = "loss",
  score_type = "pos",
  qvalue.cut = 0.20,
  pathwaysize.lower.cut = 20,
  pathwaysize.upper.cut = 300,
  nperm = nperm,
  nproc = 1
)

top_loss_go <- result_loss_go$results %>%
  arrange(padj) %>%
  filter(padj < 0.05) %>%
  dplyr::select(pathway, pval, padj, Top_Loss_Genes)

# 27 

result_cons_go <- pathSelect(
  mcmc.merge.list = list(Human_LUN = human_LUN, Baboon_LUN = baboon_LUN),
  pathway.list = go.pathway.list_hsa,
  dataset.names = c("Human_LUN", "Baboon_LUN"),
  ranking.method = "conserved",
  score_type = "pos",
  qvalue.cut = 0.20,
  pathwaysize.lower.cut = 20,
  pathwaysize.upper.cut = 200,
  nperm = nperm,
  nproc = 1
)

top_cons_go <- result_cons_go$results %>%
  arrange(padj) %>%
  filter(padj < 0.05) %>%
  dplyr::select(pathway, pval, padj, Top_Conserved_Genes)

print_pathway_summary(result_gain_go, filter_by = "q", cutoff = 0.2)
print_pathway_summary(result_loss_go, filter_by = "q", cutoff = 0.05)
print_pathway_summary(result_cons_go, filter_by = "q", cutoff = 0.05)


cat("Descriptive metrics calculated for", nrow(result_multiconservation), "pathways\n")

result_multiconservation_filtered <- result_multiconservation %>%
  select(
    Pathway,
    ILE_vs_LIV_AdjustedConcordance,
    ILE_vs_LIV_PValue,
    ILE_vs_LIV_QValue,
    ILE_vs_LIV_GainLossRatio,
  )

relevant_pathways <- result_multiconservation_filtered$Pathway



####################################

#---------------------------------------------------------------------------------
pA    <- rowMeans(baboon_ILE$rho)
pB    <- rowMeans(baboon_LIV$rho)

trans_outer <- transition_classify(pA, pB, bfdr_alpha = 0.25)
phase_inner <- phase_infer(
  phi_matrix1 = baboon_ILE$phi,
  phi_matrix2 = baboon_LIV$phi,
  gain_loss_status = trans_outer$gain_loss_status,
  bfdr_alpha = 0.25,
  shift = 2,
  P = 24, comLIVe_hdi = TRUE
)

base_path <- file.path(outLIV.dir, "heatmap_baboon_ILE_LIV")
dir.create(base_path, showWarnings = FALSE, recursive = TRUE)

###############################################################
# PLOT EACH PATHWAY
###############################################################
circadian_genes <- kegg.pathway.list_hsa[["KEGG Circadian rhythm"]]

# Replace ARNTL with BMAL1
circadian_genes <- ifelse(circadian_genes == "ARNTL", "BMAL1", circadian_genes)

# Save back into list
kegg.pathway.list_hsa[["KEGG Circadian rhythm"]] <- circadian_genes


trans_outer <- transition_classify(pA, pB, bfdr_alpha = 0.25)



if (!dir.exists(tempdir())) {
  dir.create(tempdir(), recursive = TRUE)
}


cat("\n=== PLOTTING PATHWAYS FOR BABOON ILE vs LIV ===\n")

pathways_to_plot = relevant_pathways

for (pathway in pathways_to_plot) {
  
  cat("\nProcessing:", pathway, "\n")
  
  # -------- Pathway available? --------
  if (!pathway %in% names(kegg.pathway.list_hsa)) {
    cat("   ✗ Pathway not found — skip\n")
    next
  }
  
  pathway_genes <- kegg.pathway.list_hsa[[pathway]]
  overlap <- intersect(rownames(baboon_ILE$rho), pathway_genes)
  
  if (length(overlap) == 0) {
    cat("   ✗ No overlapping genes — skip\n")
    next
  }
  
  cat("   Found", length(overlap), "overlapping genes\n")
  
  # -------- Create outLIV folder for this pathway --------
  safe_name <- gsub("[^A-Za-z0-9_]", "_", pathway)
  out_dir <- file.path(base_path, safe_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # -------- Plotting --------
  tryCatch({
    
    # Close any open graphics devices
    while (dev.cur() > 1) dev.off()
    
    plot_heatmap(
      data1 = baboon_ILE,
      data2 = baboon_LIV,
      pathway_genes = pathway_genes,
      pathway_name = pathway,
      phase_results = phase_inner,
      transition_results = trans_outer,
      group_names = c("ILE", "LIV"),
      versions = "both",
      save_path = out_dir
    )
    
    cat("   ✓ Saved to:", out_dir, "\n")
    
  }, error = function(e) {
    cat("   ✗ ERROR:", e$message, "\n")
  })
}

cat("\n=== COMPLETE: ALL PATHWAYS PLOTTED ===\n")
cat("OutLIV directory:", base_path, "\n")


################################################################################
# EXPECTED COUNTS TABLE - PROBABILISTIC ESTIMATION
# Based on posterior probabilities (not binary classification)
################################################################################

#########################################################################
# BABOON ILE vs BABOON LIV : Expected Probabilistic Counts
#########################################################################

library(dplyr)
library(kableExtra)

# ----------------------------------------------------------------------
# InLIV probabilities
# ----------------------------------------------------------------------
pA <- rowMeans(baboon_ILE$rho)   # ILE
pB <- rowMeans(baboon_LIV$rho)   # LIV

# Transition classification
trans_outer <- transition_classify(pA, pB, bfdr_alpha = 0.25)

cat("\n=== Checking Probabilistic Indices ===\n")
print(names(trans_outer))

# Total genes
n_total <- length(pA)

# ----------------------------------------------------------------------
# Extract probabilistic indices
# ----------------------------------------------------------------------
if ("p_gain" %in% names(trans_outer)) {
  
  expected_gain <- sum(trans_outer$p_gain, na.rm = TRUE)
  expected_loss <- sum(trans_outer$p_loss, na.rm = TRUE)
  expected_cons <- sum(trans_outer$p_cons, na.rm = TRUE)
  
  expected_rhythmic_ILE <- sum(pA, na.rm = TRUE)
  expected_rhythmic_LIV <- sum(pB, na.rm = TRUE)
  
} else {
  
  cat("\nCalculating probabilistic indices manually...\n")
  
  p_rhythmic_ILE <- pA
  p_rhythmic_LIV <- pB
  
  expected_rhythmic_ILE <- sum(p_rhythmic_ILE, na.rm = TRUE)
  expected_rhythmic_LIV <- sum(p_rhythmic_LIV, na.rm = TRUE)
  
  p_gain <- pB * (1 - pA)
  expected_gain <- sum(p_gain, na.rm = TRUE)
  
  p_loss <- pA * (1 - pB)
  expected_loss <- sum(p_loss, na.rm = TRUE)
  
  p_cons <- pA * pB
  expected_cons <- sum(p_cons, na.rm = TRUE)
  
  cat("Expected rhythmic ILE:", round(expected_rhythmic_ILE, 1), "\n")
  cat("Expected rhythmic LIV:", round(expected_rhythmic_LIV, 1), "\n")
  cat("Expected Gain (LIV > ILE):", round(expected_gain, 1), "\n")
  cat("Expected Loss (ILE > LIV):", round(expected_loss, 1), "\n")
  cat("Expected Conserved:", round(expected_cons, 1), "\n")
}

# ----------------------------------------------------------------------
# Helper to safely extract thresholds
# ----------------------------------------------------------------------
safe_get <- function(x) {
  if (is.null(x)) return(NA)
  unname(x)
}

thr_ILE  <- safe_get(trans_outer$tau_rhythmic_A)
thr_LIV  <- safe_get(trans_outer$tau_rhythmic_B)
thr_gain <- safe_get(trans_outer$tau_gain)
thr_loss <- safe_get(trans_outer$tau_loss)
thr_cons <- safe_get(trans_outer$tau_cons)

# ----------------------------------------------------------------------
# Expected counts summary table
# ----------------------------------------------------------------------
genome_wide_expected <- data.frame(
  Category = c(
    "Total Genes",
    "Expected Rhythmic in Baboon ILE",
    "Expected Rhythmic in Baboon LIV",
    "Expected Gain",
    "Expected Loss",
    "Expected Conserved"
  ),
  Expected_Count = c(
    n_total,
    expected_rhythmic_ILE,
    expected_rhythmic_LIV,
    expected_gain,
    expected_loss,
    expected_cons
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Percentage = sprintf("%.2f%%", 100 * Expected_Count / n_total)
  )

# Add thresholds safely
genome_wide_expected$Threshold <- c(
  NA,
  thr_ILE,
  thr_LIV,
  thr_gain,
  thr_loss,
  thr_cons
)

# ----------------------------------------------------------------------
# Print results
# ----------------------------------------------------------------------
cat("\n========================================\n")
cat("EXPECTED COUNTS (PROBABILISTIC)\n")
cat("Baboon ILE vs Baboon LIV\n")
cat("========================================\n\n")

print(genome_wide_expected)


##########
gene_names <- names(phase_inner$peak1)

phase_conserved_genes <- gene_names[phase_inner$flag_cons]
length(phase_conserved_genes)

