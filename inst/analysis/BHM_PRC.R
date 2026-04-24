# SAME SPECIE  (Baboon PRC vs Human PRC)


rm(list = ls())

current_gtex   <- "/home/qtp1/Projects/Collaborative"
current_wd     <- "/home/qtp1/Projects/Circadian"
current_aging  <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging"

outdir <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/all_plots"
# 
# load(file = "/home/qtp1/Projects/Collaborative/GTEXdata/result/summary/hb/mcmc_rho_BF3.RData")
# load(file.path(current_gtex, "GTEXdata/result/summary/hb/phi/mcmc_phi_BF3.RData"))
load("/home/qtp1/Projects/Collaborative/GTEXdata/data/CAMO_PRC_hmb.RData") # gtex, baboon_withTOD, mice

# Objects from these:
# gtex , mice,  baboon_withTOD
library(KEGGREST)
library(parallel)

options(width = 10000)

#── Packages ─────────────────────────────────────────────────────────────────────
library(parallel); library(edgeR); library(Rcpp)
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
}
source(file.path(thien_dir, "plots/heatmap.R"))


load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/hw_orth.RData"))
load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/human.pathway.list.RData"))
load(file.path(current_wd, "Kyle/Circadian-analysis-main/R/pathway_data/go.pathway.list_hsa.RData"))

kegg.pathway.list_hsa <- readRDS(
  "/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/pathway_data/kegg_pathway_list_hsa.rds"
)

#── Source R scripts ─────────────────────────────────────────────────────────────
WD <- "Kyle/Circadian-analysis-main/R/v1"
setwd(file.path(current_wd, WD))
scripts <- list.files("BayRC", pattern="\\.R$", full.names = TRUE)
sapply(scripts, source)

# Reset wd
setwd(current_aging)

#---------------------------------------------------------------------------------
to_zt <- function(t_cos) ifelse(t_cos >= 18, t_cos - 24, t_cos)

#---------------------------------------------------------------------------------
# BRAIN/PRC OBJECTS # gtex , mice,  baboon_withTOD

#---------------------------------------------------------------------------------
human_mcmc1  <- readRDS("/home/qtp1/Projects/Collaborative/GTEXdata/result/PRC/human/PRC1/gtex_PRC1_bay_1.RDS")
human_mcmc2  <- readRDS("/home/qtp1/Projects/Collaborative/GTEXdata/result/PRC/human/PRC2/gtex_PRC2_bay_1.RDS")
baboon_mcmc <- readRDS("/home/qtp1/Projects/Collaborative/GTEXdata/result/PRC/baboon/PRC/baboon_PRC_bay_1.RDS")
mouse_mcmc  <- readRDS("/home/qtp1/Projects/Collaborative/GTEXdata/result/PRC/mice/PRC/mice_PRC_bay_1.RDS")

ensemble <- try_any_mirror(dataset = "hsapiens_gene_ensembl")

human_mcmc1 <- match_symbols(human_mcmc1, BF = 3, p_rhythmic = 0.2, ensemble)
human_mcmc2 <- match_symbols(human_mcmc2, BF = 3, p_rhythmic = 0.2, ensemble)
baboon_mcmc <- match_symbols(baboon_mcmc, BF = 3, p_rhythmic = 0.2, ensemble)
mouse_mcmc  <- match_symbols(mouse_mcmc,  BF = 3, p_rhythmic = 0.2, ensemble)

baboon_PRC <- list(
  rho = baboon_mcmc$rho,
  phi = baboon_mcmc$phi
)

human_PRC1 <- list(
  rho = human_mcmc1$rho,
  phi = human_mcmc1$phi
)

human_PRC2 <- list(
  rho = human_mcmc2$rho,
  phi = human_mcmc2$phi
)

mouse_PRC <- list(
  rho = mouse_mcmc$rho,
  phi = mouse_mcmc$phi
)

#-----


#---------------------------------------------------------------------------------
# Global concordance score 
#---------------------------------------------------------------------------------
output.dir <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/results/hmb/output_final/hb"
if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
}
cat("Directory ready:", output.dir, "\n")

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
results_global1 <- multi_conservation(
  mcmc.merge.list   = list(baboon_PRC, human_PRC2),
  dataset.names     = c("Baboon_PRC", "Human_PRC"),
  select.pathway.list = "global",
  n_perm            = 1000,
  n_boot            = 1000,
  output.dir        = file.path(output.dir, "concordance_PRC"),
  use_cpp           = TRUE
)

rhy <- detect_rhy(baboon_PRC, human_PRC2, 0.25)

#---------------------------------------------------------------------------------
# Outer (rhythmic) and inner (phase) analysis
#---------------------------------------------------------------------------------

# dataset1 = Baboon PRC, dataset2 = Human PRC
pA <- rowMeans(baboon_PRC$rho)
pB <- rowMeans(human_PRC2$rho)

trans_outer <- transition_classify(pA, pB, bfdr_alpha = 0.25)

phase_inner <- phase_infer(
  phi_matrix1      = baboon_PRC$phi,
  phi_matrix2      = human_PRC2$phi,
  gain_loss_status = trans_outer$gain_loss_status,
  bfdr_alpha       = 0.25,
  shift            = 2,
  P                = 24
)

#---------------------------------------------------------------------------------
# Maintained plot: Peak Concordance Plot - Baboon PRC vs Human PRC
#---------------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(ggrepel)

# Helper (already defined above, but keep here local if needed)
to_zt <- function(t_cos) ifelse(t_cos >= 18, t_cos - 24, t_cos)

#───────────────────────────────────────────────────────────────
# Extract data from phase_inner
#───────────────────────────────────────────────────────────────
gene_names <- names(phase_inner$peak1)

maintained_df <- data.frame(
  Gene           = gene_names,
  Peak_BaboonPRC = phase_inner$peak1,
  Peak_HumanPRC  = phase_inner$peak2,
  deltaPhi       = phase_inner$deltaPhi.Est,
  prob_shift     = phase_inner$prob_shift,
  prob_conserved = phase_inner$prob_conserved,
  stringsAsFactors = FALSE
)

# Filter for maintained genes only
maintained_genes <- names(trans_outer$gain_loss_status[
  trans_outer$gain_loss_status == "Maintained"
])

maintained_df <- maintained_df %>%
  filter(Gene %in% maintained_genes) %>%
  mutate(
    Peak_BaboonPRC_ZT = to_zt(Peak_BaboonPRC),
    Peak_HumanPRC_ZT  = to_zt(Peak_HumanPRC)
  )

#───────────────────────────────────────────────────────────────
# Assign phase classification
#───────────────────────────────────────────────────────────────
phase_class <- rep("Undetermined", length(gene_names))
names(phase_class) <- gene_names
phase_class[phase_inner$flag_cons]        <- "Phase-conserved"
phase_class[phase_inner$flag_shift]       <- "Phase-shifted"
phase_class[phase_inner$flag_undetermined] <- "Undetermined"

maintained_df <- maintained_df %>%
  mutate(
    phase_class = phase_class[Gene],
    phase_class = ifelse(is.na(phase_class), "Undetermined", phase_class)
  )

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

n_rhythmic_outer   <- length(rhythmic_genes_outer)
n_phase_conserved  <- sum(phase_inner$flag_cons, na.rm = TRUE)
n_phase_shifted    <- sum(phase_inner$flag_shift, na.rm = TRUE)
n_undetermined     <- sum(phase_inner$flag_undetermined, na.rm = TRUE)

pct_phase_conserved_global <- 100 * n_phase_conserved / n_rhythmic_outer
pct_phase_shifted_global   <- 100 * n_phase_shifted / n_rhythmic_outer
pct_undetermined_global    <- 100 * n_undetermined / n_rhythmic_outer

cat("\n[Global Summary - Baboon PRC vs Human PRC]\n")
cat("  Rhythmic (outer τ_c ≤ 0.25):", n_rhythmic_outer, "\n")
cat("  Phase-conserved (inner τ_p ≤ 0.10):", n_phase_conserved,
    sprintf("(%.2f%%)\n", pct_phase_conserved_global))
cat("  Phase-shifted (inner τ_p ≤ 0.10):", n_phase_shifted,
    sprintf("(%.2f%%)\n", pct_phase_shifted_global))
cat("  Undetermined:", n_undetermined,
    sprintf("(%.2f%%)\n", pct_undetermined_global))

#───────────────────────────────────────────────────────────────
# Compute ±3 h concordance among maintained genes
#───────────────────────────────────────────────────────────────
calculate_peak_difference <- function(a, b) {
  d <- abs(a - b)
  ifelse(d > 12, 24 - d, d)
}

maintained_df <- maintained_df %>%
  mutate(
    peak_diff          = calculate_peak_difference(Peak_BaboonPRC_ZT, Peak_HumanPRC_ZT),
    within_concordance = peak_diff <= 3
  )

n_total  <- nrow(maintained_df)
n_within <- sum(maintained_df$within_concordance, na.rm = TRUE)
pct_within <- 100 * n_within / n_total

cat("\n[Within ±3 h summary]\n")
cat("  Maintained genes:", n_total, "\n")
cat("  Within ±3 h:", n_within, "\n")
cat(sprintf("  ⇒ %.1f%% within ±3 h interval\n", pct_within))

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
    remapped = remap_to_diagonal(Peak_BaboonPRC_ZT, Peak_HumanPRC_ZT),
    Peak_BaboonPRC_ZT_plot = remapped$x,
    Peak_HumanPRC_ZT_plot  = remapped$y
  ) %>%
  select(-remapped)

n_Rc   <- nrow(plot_df)
pct_Rc <- round(100 * n_within / n_Rc, 1)

subtitle_text <- bquote(
  "Rhythmically Conserved Set " ~ R[c] ~
    "(" * n == .(n_Rc) * ", " * .(pct_Rc) * "% within " * "\u00B1" * "3 h interval)"
)

#───────────────────────────────────────────────────────────────
# Create plot
#───────────────────────────────────────────────────────────────
p <- ggplot(plot_df, aes(
  x = Peak_BaboonPRC_ZT_plot,
  y = Peak_HumanPRC_ZT_plot,
  color = phase_class
)) +
  geom_abline(intercept = 0,  slope = 1,
              color = "black",    linetype = "dashed", linewidth = 1.0) +
  geom_abline(intercept = 3,  slope = 1,
              color = "darkgreen", linetype = "dotted", linewidth = 1.1) +
  geom_abline(intercept = -3, slope = 1,
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
    title    = "Circadian Peak Concordance: Baboon PRC versus Human PRC",
    subtitle = subtitle_text,
    x        = "Peak Hour – Baboon PRC (ZT)",
    y        = "Peak Hour – Human PRC (ZT)",
    color    = "Phase class"
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
    plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, margin = margin(b = 10)),
    axis.title    = element_text(face = "bold", size = 13),
    axis.text     = element_text(size = 12),
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.box        = "horizontal",
    legend.title      = element_text(face = "bold"),
    legend.background = element_rect(color = "gray70", fill = "white"),
    legend.key        = element_rect(fill = "white", color = NA),
    legend.margin     = margin(t = 5, b = 5),
    plot.margin       = margin(15, 15, 15, 15)
  )

p

#───────────────────────────────────────────────────────────────
# Save
#───────────────────────────────────────────────────────────────
save_dir <- file.path(output.dir, "figure/baboon_human_PRC")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(save_dir, "Baboon_Human_PRC_Peak_Concordance_0.15.pdf"),
  plot = p, width = 9, height = 8
)

#───────────────────────────────────────────────────────────────
#───────────────────────────────────────────────────────────────
################################################################################
# COMPLETE WORKFLOW: Baboon PRC vs Human PRC Pathway Analysis
################################################################################

# Set parameters
dataset_names <- c("Baboon_PRC", "Human_PRC")
qvalue_cut <- 0.25
nperm <- 10000
pathway_size_min <- 10
pathway_size_max <- 300

# Set output directory
output.dir <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/results/baboon/output_final"
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

################################################################################
# STAGE 1: Identify Active Pathways (Union Test - P-value Filter)
################################################################################

cat("\n=== STAGE 1: Union Test (Active Pathways) ===\n")

result_union <- pathSelect(
  mcmc.merge.list = list(Baboon_PRC = baboon_PRC, Human_PRC = human_PRC),
  pathway.list = kegg.pathway.list_hsa,
  dataset.names = c("Baboon_PRC", "Human_PRC"),
  ranking.method = "union",
  score_type = "pos",
  qvalue.cut = 0.9,
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

cat("Active pathways (p < 0.05):", length(active_pathways), "\n")

################################################################################
# STAGE 2: Test Transitions on Active Pathways (Q-value Filter)
################################################################################

cat("\n=== STAGE 2: Transition Tests (Gain, Loss, Conservation) ===\n")

# Ensure correct order
active_pathway_list <- kegg.pathway.list_hsa[match(active_pathways, names(kegg.pathway.list_hsa))]

# Gain enrichment - Q-VALUE FILTERING
result_gain <- pathSelect(
  mcmc.merge.list = list(Baboon_PRC = baboon_PRC, Human_PRC = human_PRC),
  pathway.list = kegg.pathway.list_hsa,
  dataset.names = c("Baboon_PRC", "Human_PRC"),
  ranking.method = "gain",
  score_type = "pos",
  qvalue.cut = 0.05,
  pathwaysize.lower.cut = 10,
  pathwaysize.upper.cut = 500,
  nperm = nperm,
  nproc = 1
)

# Loss enrichment - Q-VALUE FILTERING
result_loss <- pathSelect(
  mcmc.merge.list = list(Baboon_PRC = baboon_PRC, Human_PRC = human_PRC),
  pathway.list = kegg.pathway.list_hsa,
  dataset.names = c("Baboon_PRC", "Human_PRC"),
  ranking.method = "loss",
  score_type = "pos",
  qvalue.cut = 0.05,
  pathwaysize.lower.cut = pathway_size_min,
  pathwaysize.upper.cut = 300,
  nperm = nperm,
  nproc = 1
)


result_gain$results %>%
  filter(pval < 0.05) %>%
  select(pathway, size, Gain_Index, Expected_N_Gain, pval, padj, Gain_Loss_Ratio_Arithmetic)


result_loss$results %>%
  filter(pval < 0.05) %>%
  select(pathway, size, Loss_Index, Expected_N_Loss, pval, padj, Gain_Loss_Ratio_Arithmetic)

result_cons$results %>%
  filter(pval < 0.05) %>%
  select(pathway, size, Conserved_Index, Expected_N_Conserved, pval, padj, Gain_Loss_Ratio_Arithmetic)


result_cons$results %>%
  filter(pathway == "KEGG Ribosome") %>%
  select(pathway, size, Conserved_Index, Expected_N_Conserved, Top_Conserved_Genes, pval, padj, Gain_Loss_Ratio_Arithmetic)


# Conservation enrichment - Q-VALUE FILTERING
result_cons <- pathSelect(
  mcmc.merge.list = list(Baboon_PRC = baboon_PRC, Human_PRC = human_PRC),
  pathway.list = kegg.pathway.list_hsa,
  dataset.names = c("Baboon_PRC", "Human_PRC"),
  ranking.method = "conserved",
  score_type = "pos",
  qvalue.cut = 0.20,
  pathwaysize.lower.cut = pathway_size_min,
  pathwaysize.upper.cut = pathway_size_max,
  nperm = nperm,
  nproc = 1
)

cat("\nEnrichment Results (q < 0.20):\n")
cat("  Gain-enriched (Human_PRC):  ", sum(result_gain$results$Significant), "\n")
cat("  Loss-enriched (Human_PRC):  ", sum(result_loss$results$Significant), "\n")
cat("  Conservation-enriched:", sum(result_cons$results$Significant), "\n")

print_pathway_summary <- function(result_obj, 
                                  filter_by = c("q", "p"),
                                  cutoff = 0.2) {
  
  filter_by <- match.arg(filter_by)   # ensure "p" or "q"
  
  # Extract table
  df <- result_obj$results
  
  # Compute expected union and ratio
  df$Expected_Union <- df$Expected_N_Gain + df$Expected_N_Loss + df$Expected_N_Conserved
  df$Union_Size_Ratio <- df$Expected_Union / df$size
  
  # Decide filtering method
  if (filter_by == "q") {
    df_sig <- df[df$padj < cutoff, ]
  } else {
    df_sig <- df[df$pval < cutoff, ]
  }
  
  # Print formatted output
  apply(df_sig, 1, function(x) {
    cat(sprintf(
      "%s (Expected gain = %.1f; Expected loss = %.1f; Expected conserved = %.1f; Expected union = %.1f; Expected union / size = %.1f / %d = %.3f; Gain/Loss ratio = %.3f; p value = %.4g; q value = %.4g)\n\n",
      x["pathway"],
      as.numeric(x["Expected_N_Gain"]),
      as.numeric(x["Expected_N_Loss"]),
      as.numeric(x["Expected_N_Conserved"]),
      as.numeric(x["Expected_Union"]),
      as.numeric(x["Expected_Union"]),
      as.numeric(x["size"]),
      as.numeric(x["Union_Size_Ratio"]),
      as.numeric(x["Gain_Loss_Ratio_Arithmetic"]),   # <-- used directly
      as.numeric(x["pval"]),
      as.numeric(x["padj"])
    ))
  })
}

print_pathway_summary(result_gain, filter_by = "q", cutoff = 0.2)
print_pathway_summary(result_loss, filter_by = "q", cutoff = 0.2)
print_pathway_summary(result_cons, filter_by = "q", cutoff = 0.2)

################################################################################
# STAGE 3: Filter to Significant Pathways
################################################################################

cat("\n=== STAGE 3: Filtering Significant Pathways ===\n")

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
  mcmc.merge.list = list(Baboon_PRC = baboon_PRC, Human_PRC = human_PRC),
  dataset.names = c("Baboon_PRC", "Human_PRC"),
  select.pathway.list = filtered_pathway_list,
  n_perm = 1000,
  n_boot = 1000,
  output.dir = file.path(output.dir, "multiconservation_filtered_PRC"),
  use_cpp = TRUE
)

cat("Descriptive metrics calculated for", nrow(result_multiconservation), "pathways\n")

result_multiconservation_filtered <- result_multiconservation %>%
  select(
    Pathway,
    Baboon_PRC_vs_Human_PRC_AdjustedCongruence,
    Baboon_PRC_vs_Human_PRC_PValue,
    Baboon_PRC_vs_Human_PRC_QValue,
    Baboon_PRC_vs_Human_PRC_GainLossRatio
  )

relevant_pathways <- names(filtered_pathway_list)

result_multiconservation_filtered <- result_multiconservation_filtered %>%
  filter(Pathway %in% relevant_pathways)

###############################################
# PATHWAY LIST TO PLOT - BABOON vs HUMAN PRC
###############################################


###############################################
# SET UP OUTPUT DIRECTORY
###############################################

base_path <- file.path(output.dir, "heatmap_baboon_human_PRC")
dir.create(base_path, showWarnings = FALSE, recursive = TRUE)

###############################################################
# PLOT EACH PATHWAY
###############################################################
pA <- rowMeans(human_PRC2$rho)
pB <- rowMeans(mouse_PRC$rho)

circadian_genes <- kegg.pathway.list_hsa[["KEGG Circadian rhythm"]]

# Replace ARNTL with BMAL1
circadian_genes <- ifelse(circadian_genes == "ARNTL", "BMAL1", circadian_genes)

# Save back into list
kegg.pathway.list_hsa[["KEGG Circadian rhythm"]] <- circadian_genes

trans_outer <- transition_classify(pA, pB, bfdr_alpha = 0.20)

# phase_inner <- phase_infer(
#   phi_matrix1      = human_PRC$phi,
#   phi_matrix2      = baboon_PRC$phi,
#   gain_loss_status = trans_outer$gain_loss_status,
#   bfdr_alpha       = 0.15,
#   shift            = 3,
#   P                = 24,
#   compute_hdi      = TRUE
# )

pathways_to_plot = relevant_pathways
source("/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/code/heatmap.R")

cat("\n=== PLOTTING PATHWAYS FOR BABOON vs HUMAN PRC ===\n")
if (!dir.exists(tempdir())) {
  dir.create(tempdir(), recursive = TRUE)
}
for (pathway in pathways_to_plot) {
  
  cat("\nProcessing:", pathway, "\n")
  
  # -------- Pathway available? --------
  if (!pathway %in% names(kegg.pathway.list_hsa)) {
    cat("   ✗ Pathway not found — skip\n")
    next
  }
  
  pathway_genes <- kegg.pathway.list_hsa[[pathway]]
  overlap <- intersect(rownames(baboon_PRC$rho), pathway_genes)
  
  if (length(overlap) == 0) {
    cat("   ✗ No overlapping genes — skip\n")
    next
  }
  
  cat("   Found", length(overlap), "overlapping genes\n")
  
  # -------- Create output folder for this pathway --------
  safe_name <- gsub("[^A-Za-z0-9_]", "_", pathway)
  out_dir <- file.path(base_path, safe_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # -------- Plotting --------
  tryCatch({
    
    # Close any open graphics devices
    while (dev.cur() > 1) dev.off()
    
    plot_pathway_integrated(
      data1 = human_PRC,
      data2 = baboon_PRC,
      pathway_genes = pathway_genes,
      pathway_name = pathway,
      phase_results = phase_inner,
      transition_results = trans_outer,
      group_names = c("Human PRC", "Baboon PRC"),
      versions = "both",
      save_path = out_dir
    )
    
    cat("   ✓ Saved to:", out_dir, "\n")
    
  }, error = function(e) {
    cat("   ✗ ERROR:", e$message, "\n")
  })
}

cat("\n=== COMPLETE: ALL PATHWAYS PLOTTED ===\n")
cat("Output directory:", base_path, "\n")


##########

################################################################################
# GENOME-WIDE EXPECTED COUNTS TABLE
################################################################################
################################################################################
# EXPECTED COUNTS TABLE - PROBABILISTIC ESTIMATION
# Based on posterior probabilities (not binary classification)
################################################################################

library(dplyr)
library(kableExtra)

# For Baboon PRC vs Human PRC analysis
pA <- rowMeans(baboon_PRC$rho)
pB <- rowMeans(human_PRC2$rho)

# Get transition classification with probabilistic indices
trans_outer <- transition_classify(pA, pB, bfdr_alpha = 0.25)

cat("\n=== Checking Probabilistic Indices ===\n")
cat("Names in trans_outer:\n")
print(names(trans_outer))

# Extract the probabilistic indices (these should be vectors for each gene)
# Expected counts are the SUM of these probabilities across all genes
n_total <- length(pA)

# Check if trans_outer has the probability vectors
if ("p_gain" %in% names(trans_outer)) {
  expected_gain <- sum(trans_outer$p_gain, na.rm = TRUE)
  expected_loss <- sum(trans_outer$p_loss, na.rm = TRUE)
  expected_cons <- sum(trans_outer$p_cons, na.rm = TRUE)
} else {
  # Calculate manually if not in trans_outer
  cat("\nCalculating probabilistic indices manually...\n")
  
  # Probability of being rhythmic in each condition
  p_rhythmic_A <- pA
  p_rhythmic_B <- pB
  
  # Expected number rhythmic in each condition
  expected_rhythmic_baboon <- sum(p_rhythmic_A, na.rm = TRUE)
  expected_rhythmic_human <- sum(p_rhythmic_B, na.rm = TRUE)
  
  # Gain: rhythmic in B (human) but not in A (baboon)
  # p_gain = pB * (1 - pA)
  p_gain <- pB * (1 - pA)
  expected_gain <- sum(p_gain, na.rm = TRUE)
  
  # Loss: rhythmic in A (baboon) but not in B (human)
  # p_loss = pA * (1 - pB)
  p_loss <- pA * (1 - pB)
  expected_loss <- sum(p_loss, na.rm = TRUE)
  
  # Conserved: rhythmic in both A and B
  # p_cons = pA * pB
  p_cons <- pA * pB
  expected_cons <- sum(p_cons, na.rm = TRUE)
  
  cat("Expected rhythmic Baboon:", round(expected_rhythmic_baboon, 1), "\n")
  cat("Expected rhythmic Human:", round(expected_rhythmic_human, 1), "\n")
  cat("Expected Gain:", round(expected_gain, 1), "\n")
  cat("Expected Loss:", round(expected_loss, 1), "\n")
  cat("Expected Conserved:", round(expected_cons, 1), "\n")
}

# Create expected counts summary
genome_wide_expected <- data.frame(
  Category = c(
    "Total Genes", 
    "Expected Rhythmic in Baboon PRC", 
    "Expected Rhythmic in Human PRC",
    "Expected Gain (Human > Baboon)", 
    "Expected Loss (Human < Baboon)", 
    "Expected Maintained (R_c)"
  ),
  Expected_Count = c(
    n_total,
    expected_rhythmic_baboon,
    expected_rhythmic_human,
    expected_gain,
    expected_loss,
    expected_cons
  ),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Percentage = sprintf("%.2f%%", 100 * Expected_Count / n_total)
  )

# Add thresholds
genome_wide_expected$Threshold <- c(
  NA,
  trans_outer$tau_rhythmic_A,
  trans_outer$tau_rhythmic_B,
  trans_outer$tau_gain,
  trans_outer$tau_loss,
  trans_outer$tau_cons
)

cat("\n========================================\n")
cat("EXPECTED COUNTS (PROBABILISTIC)\n")
cat("Baboon PRC vs Human PRC\n")
cat("========================================\n\n")

print(genome_wide_expected)

# Verification
cat("\n=== VERIFICATION ===\n")
cat("Expected Gain + Loss + Conserved =", 
    round(expected_gain + expected_loss + expected_cons, 1), "\n")




####
library(ggplot2)
library(dplyr)
library(tidyr)
library(circular)

circular_mean_24h <- function(samples, P = 24) {
  rad <- circular((samples / P) * 2 * pi, units = "radians")
  m <- mean.circular(rad)
  t <- (as.numeric(m) / (2 * pi)) * P
  t %% P
}

circadian_genes <- c(
  "BHLHE40","BHLHE41","BMAL1","ARNTL","BTRC","CLOCK","CREB1",
  "CRY1","CRY2","CSNK1D","CSNK1E","CUL1","DBP","FBXL3","FBXW11",
  "NFIL3","NPAS2","NR1D1","NR1D2","PER1","PER2","PER3",
  "RORA","RORB","RORC","CIART","TEF","HLF"
)

maintained_genes <- names(trans_outer$gain_loss_status)[trans_outer$gain_loss_status == "Maintained"]
conserved_genes <- intersect(circadian_genes, maintained_genes)

results <- data.frame(
  Gene = conserved_genes,
  Mean_Baboon_hr = sapply(conserved_genes, function(g) circular_mean_24h(baboon_PRC$phi[g,])),
  Mean_Human_hr = sapply(conserved_genes, function(g) circular_mean_24h(human_PRC$phi[g,])),
  SD_Baboon_hr = sapply(conserved_genes, function(g) {
    r <- circular((baboon_PRC$phi[g,] / 24) * 2 * pi, units = "radians")
    (sd.circular(r) / (2*pi)) * 24
  }),
  SD_Human_hr = sapply(conserved_genes, function(g) {
    r <- circular((human_PRC$phi[g,] / 24) * 2 * pi, units = "radians")
    (sd.circular(r) / (2*pi)) * 24
  })
)

plot_data <- results %>%
  pivot_longer(
    cols = -Gene,
    names_to = c(".value", "Species"),
    names_pattern = "(.+)_(Baboon|Human)_hr"
  ) %>%
  mutate(
    Species = factor(Species, levels = c("Baboon", "Human")),
    gene_id = as.numeric(factor(Gene)) + 2,  # SKIP CENTER - START AT 3
    
    rad = (Mean / 24) * 2 * pi,
    rad_sd = (SD / 24) * 2 * pi,
    rad_lower = rad - rad_sd,
    rad_upper = rad + rad_sd
  )

n_genes <- length(unique(plot_data$Gene))
gene_colors <- setNames(
  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
    "#A6D854", "#A65628", "#F781BF", "#999999", "#66C2A5",
    "#FC8D62", "#8DA0CB", "#E78AC3")[1:n_genes],
  unique(plot_data$Gene)
)

p <- ggplot(plot_data) +
  # CI bands - single layer, low opacity
  geom_segment(
    aes(x = rad_lower, xend = rad_upper, y = gene_id, yend = gene_id, 
        color = Gene, linetype = Species),
    linewidth = 8, alpha = 0.3, lineend = "butt"
  ) +
  
  # Centroids
  geom_point(
    aes(x = rad, y = gene_id, fill = Gene, shape = Species),
    color = "white", size = 5, stroke = 2, alpha = 0.9
  ) +
  
  coord_polar(start = -pi/2, direction = 1) +
  scale_x_continuous(
    limits = c(0, 2*pi),
    breaks = seq(0, 22, 2) * 2*pi/24,
    labels = paste0("ZT", seq(0, 22, 2))
  ) +
  scale_y_continuous(limits = c(0, n_genes + 3)) +
  scale_color_manual(values = gene_colors) +
  scale_fill_manual(values = gene_colors) +
  scale_shape_manual(values = c("Baboon" = 24, "Human" = 21)) +
  scale_linetype_manual(values = c("Baboon" = "solid", "Human" = "dashed")) +
  labs(
    title = "Circadian Peak Timing of Phase-Conserved Genes",
    subtitle = "Posterior phase estimates (ZT) ± SD"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40"),
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid.major.x = element_line(color = "gray88", linewidth = 0.5),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.background = element_rect(fill = "white", color = "gray70", linewidth = 0.5),
    legend.margin = margin(10, 10, 10, 10),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 4, shape = 21, stroke = 1.5, alpha = 1)),
    color = "none",
    shape = guide_legend(order = 2, override.aes = list(size = 4, stroke = 1.5, fill = "gray50", color = "white", alpha = 1)),
    linetype = "none"
  )
save_dir <- file.path(output.dir, "figure/baboon_human_PRC")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(save_dir, "Baboon_Human_PRC_circular_phase_plot.pdf"), p, width = 12, height = 10, dpi = 800)
ggsave(file.path(save_dir, "Baboon_Human_PRC_circular_phase_plot.png"), p, width = 12, height = 10, dpi = 800)

##### Test new function

################################################################################
# TEST: Iteration-Level Jaccard Index
################################################################################

# Compute Jaccard for each iteration
result <- compute_adjusted_jaccard_analytical_pvalue(
  rho_A = baboon_PRC$rho,
  rho_B = human_PRC$rho
)

# Summary statistics
jaccard_dist <- result$jaccard_obs_mean
jaccard_mean <- mean(jaccard_dist, na.rm = TRUE)
jaccard_sd <- sd(jaccard_dist, na.rm = TRUE)
jaccard_ci <- quantile(jaccard_dist, probs = c(0.025, 0.975), na.rm = TRUE)

cat("\nPosterior Distribution of Jaccard Index:\n")
cat(sprintf("  Mean: %.4f\n", jaccard_mean))
cat(sprintf("  SD:   %.4f\n", jaccard_sd))
cat(sprintf("  95%% CI: [%.4f, %.4f]\n", jaccard_ci[1], jaccard_ci[2]))

# Summary of confusion matrix across iterations
confusion_summary <- colMeans(result$confusion, na.rm = TRUE)
cat("\nAverage Confusion Matrix Elements:\n")
cat(sprintf("  Both rhythmic (a):     %.1f genes\n", confusion_summary["a_both"]))
cat(sprintf("  Loss (b):              %.1f genes\n", confusion_summary["b_loss"]))
cat(sprintf("  Gain (c):              %.1f genes\n", confusion_summary["c_gain"]))
cat(sprintf("  Neither rhythmic (d):  %.1f genes\n", confusion_summary["d_neither"]))

# Gain/Loss ratio distribution
gain_loss_ratio <- result$confusion[, "c_gain"] / result$confusion[, "b_loss"]
cat("\nGain/Loss Ratio:\n")
cat(sprintf("  Mean: %.4f\n", mean(gain_loss_ratio, na.rm = TRUE)))
cat(sprintf("  95%% CI: [%.4f, %.4f]\n", 
            quantile(gain_loss_ratio, 0.025, na.rm = TRUE),
            quantile(gain_loss_ratio, 0.975, na.rm = TRUE)))


#####
################################################################################
# CONCISE PAIRWISE TISSUE CONCORDANCE ANALYSIS
################################################################################

# Core function for concordance (minimal output)

################################################################################
# 1. WITHIN-HUMAN PAIRWISE
################################################################################

tissues_human <- names(mcmc_data_human)
n_human <- length(tissues_human)
pairs_human <- combn(tissues_human, 2, simplify = FALSE)

results_within_human <- data.frame(
  Tissue1 = character(length(pairs_human)),
  Tissue2 = character(length(pairs_human)),
  Jaccard_Obs = numeric(length(pairs_human)),
  Jaccard_Adj = numeric(length(pairs_human)),
  Jaccard_Null = numeric(length(pairs_human)),
  stringsAsFactors = FALSE
)

cat(sprintf("Computing %d within-human pairs...\n", length(pairs_human)))
pb <- txtProgressBar(min = 0, max = length(pairs_human), style = 3)

for (i in seq_along(pairs_human)) {
  setTxtProgressBar(pb, i)
  t1 <- pairs_human[[i]][1]
  t2 <- pairs_human[[i]][2]
  
  result <- compute_concordance_minimal(
    mcmc_data_human[[t1]], 
    mcmc_data_human[[t2]]
  )
  
  results_within_human$Tissue1[i] <- t1
  results_within_human$Tissue2[i] <- t2
  results_within_human$Jaccard_Obs[i] <- result$jaccard_obs
  results_within_human$Jaccard_Adj[i] <- result$jaccard_adj
  results_within_human$Jaccard_Null[i] <- result$jaccard_null
}
close(pb)

################################################################################
# 2. WITHIN-BABOON PAIRWISE
################################################################################

tissues_baboon <- names(mcmc_data_baboon)
n_baboon <- length(tissues_baboon)
pairs_baboon <- combn(tissues_baboon, 2, simplify = FALSE)

results_within_baboon <- data.frame(
  Tissue1 = character(length(pairs_baboon)),
  Tissue2 = character(length(pairs_baboon)),
  Jaccard_Obs = numeric(length(pairs_baboon)),
  Jaccard_Adj = numeric(length(pairs_baboon)),
  Jaccard_Null = numeric(length(pairs_baboon)),
  stringsAsFactors = FALSE
)

cat(sprintf("\nComputing %d within-baboon pairs...\n", length(pairs_baboon)))
pb <- txtProgressBar(min = 0, max = length(pairs_baboon), style = 3)

for (i in seq_along(pairs_baboon)) {
  setTxtProgressBar(pb, i)
  t1 <- pairs_baboon[[i]][1]
  t2 <- pairs_baboon[[i]][2]
  
  result <- compute_concordance_minimal(
    mcmc_data_baboon[[t1]], 
    mcmc_data_baboon[[t2]]
  )
  
  results_within_baboon$Tissue1[i] <- t1
  results_within_baboon$Tissue2[i] <- t2
  results_within_baboon$Jaccard_Obs[i] <- result$jaccard_obs
  results_within_baboon$Jaccard_Adj[i] <- result$jaccard_adj
  results_within_baboon$Jaccard_Null[i] <- result$jaccard_null
}
close(pb)

################################################################################
# 3. CROSS-SPECIES (MATCHING TISSUES)
################################################################################

matching_tissues <- intersect(tissues_baboon, tissues_human)

results_cross_species <- data.frame(
  Tissue = character(length(matching_tissues)),
  Jaccard_Obs = numeric(length(matching_tissues)),
  Jaccard_Adj = numeric(length(matching_tissues)),
  Jaccard_Null = numeric(length(matching_tissues)),
  stringsAsFactors = FALSE
)

cat(sprintf("\nComputing %d cross-species pairs...\n", length(matching_tissues)))
pb <- txtProgressBar(min = 0, max = length(matching_tissues), style = 3)

for (i in seq_along(matching_tissues)) {
  setTxtProgressBar(pb, i)
  tissue <- matching_tissues[i]
  
  result <- compute_concordance_minimal(
    mcmc_data_baboon[[tissue]], 
    mcmc_data_human[[tissue]]
  )
  
  results_cross_species$Tissue[i] <- tissue
  results_cross_species$Jaccard_Obs[i] <- result$jaccard_obs
  results_cross_species$Jaccard_Adj[i] <- result$jaccard_adj
  results_cross_species$Jaccard_Null[i] <- result$jaccard_null
}
close(pb)

################################################################################
# SAVE RESULTS
################################################################################

write.csv(results_within_human, 
          file.path(output.dir, "concordance_within_human.csv"), 
          row.names = FALSE)

write.csv(results_within_baboon, 
          file.path(output.dir, "concordance_within_baboon.csv"), 
          row.names = FALSE)

write.csv(results_cross_species, 
          file.path(output.dir, "concordance_cross_species.csv"), 
          row.names = FALSE)

cat("\n\nDone! Results saved.\n")

# Sort and display top 10 for each
cat("\n=== TOP 10 WITHIN-HUMAN ===\n")
print(head(results_within_human[order(-results_within_human$Jaccard_Adj), ], 20))

cat("\n=== TOP 10 WITHIN-BABOON ===\n")
print(head(results_within_baboon[order(-results_within_baboon$Jaccard_Adj), ], 20))

cat("\n=== TOP 10 CROSS-SPECIES ===\n")
print(head(results_cross_species[order(-results_cross_species$Jaccard_Adj), ], 20))
# Humans: small but universal circadian program; Baboons: large but tissue-diversified circadian program


########


trans_outer <- transition_classify_marginal(pA, pB, bfdr_alpha = 0.20)
ribosome_status <- trans_outer$gain_loss_status[ribosome_genes]
ribosome_status <- ribosome_status[!is.na(ribosome_status)]
table(ribosome_status)



