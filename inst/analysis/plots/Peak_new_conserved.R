#---------------------------------------------------------------------------------

# Maintained plot 
#───────────────────────────────────────────────────────────────
# Peak Concordance Plot for Maintained Genes - BABOON SCN vs HIP
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
  Peak_Young = phase_inner$peak1,
  Peak_Old = phase_inner$peak2,
  deltaPhi = phase_inner$deltaPhi.Est,
  prob_shift = phase_inner$prob_shift,
  prob_conserved = phase_inner$prob_conserved,
  stringsAsFactors = FALSE
)

# Filter for maintained genes only
maintained_genes <- names(trans_outer$gain_loss_status[trans_outer$gain_loss_status == "Maintained"])

maintained_df <- maintained_df %>%
  filter(Gene %in% maintained_genes) %>%
  mutate(Peak_SCN_ZT = to_zt(Peak_SCN),
         Peak_HIP_ZT = to_zt(Peak_HIP))

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

cat("\n[Global Summary - SCN vs HIP]\n")
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
  mutate(peak_diff = calculate_peak_difference(Peak_SCN_ZT, Peak_HIP_ZT),
         within_concordance = peak_diff <= 3)

n_total <- nrow(maintained_df)
n_within <- sum(maintained_df$within_concordance, na.rm = TRUE)
pct_within <- 100 * n_within / n_total

cat("\n[Within ±3 h summary]\n")
cat("  Maintained genes:", n_total, "\n")
cat("  Within ±3 h:", n_within, "\n")
cat(sprintf("  ⇒ %.1f%% within ±3 h interval\n", pct_within))

#───────────────────────────────────────────────────────────────
# Plot preparation
#───────────────────────────────────────────────────────────────
plot_df <- maintained_df %>%
  mutate(
    Peak_SCN_ZT_plot = Peak_SCN_ZT,
    Peak_HIP_ZT_plot = Peak_HIP_ZT
  )

# Count categories
n_Rc <- nrow(plot_df)
pct_Rc <- round(100 * n_within / n_Rc, 1)

# Subtitle
subtitle_text <- bquote(
  "Rhythmically Conserved Set " ~ R[c] ~
    "(" * n == .(n_Rc) * ", " * .(pct_Rc) * "% within " * "\u00B1" * "3 h interval)"
)

#───────────────────────────────────────────────────────────────
# Create plot
#───────────────────────────────────────────────────────────────
p <- ggplot(plot_df, aes(
  x = Peak_SCN_ZT_plot,
  y = Peak_HIP_ZT_plot,
  color = phase_class
)) +
  # Identity and ±3 h boundaries
  geom_abline(intercept = 0, slope = 1,
              color = "black", linetype = "dashed", linewidth = 1.0) +
  geom_abline(intercept = 3, slope = 1,
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
    title = "Circadian Peak Concordance: Baboon SCN vs Hippocampus",
    subtitle = subtitle_text,
    x = "Peak Hour – SCN (ZT)",
    y = "Peak Hour – Hippocampus (ZT)",
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
save_dir <- file.path(output.dir, "figure/baboon_brain")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path("Baboon_SCN_HIP_Peak_Concordance_0.15.pdf"),
       plot = p, width = 9, height = 8)
