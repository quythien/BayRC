# Concordance plotting helpers (from former run_concordance_*.R)
# Source this after the corresponding analysis has created phase_inner, trans_outer, and output.dir.
# Regenerate concordance figure for Baboon_Human_LUN
# Source in the LUN screen session where all objects are loaded

library(ggplot2); library(dplyr); library(ggrepel)
if (!dir.exists(tempdir())) dir.create(tempdir(), recursive = TRUE)

to_zt <- function(t_cos) ifelse(t_cos >= 18, t_cos - 24, t_cos)

remap_to_diagonal <- function(x, y, offset = 24) {
  x_options <- c(x, x - offset, x + offset)
  y_options <- c(y, y - offset, y + offset)
  best_dist <- Inf; best_x <- x; best_y <- y
  for (x_opt in x_options) {
    for (y_opt in y_options) {
      dist <- abs(y_opt - x_opt)
      if (dist < best_dist - 1e-9) { best_dist <- dist; best_x <- x_opt; best_y <- y_opt }
    }
  }
  return(data.frame(x = best_x, y = best_y, dist = best_dist))
}

gene_names <- names(phase_inner$peak1)

maintained_df <- data.frame(
  Gene           = gene_names,
  Peak_BaboonLUN = phase_inner$peak1,
  Peak_HumanLUN  = phase_inner$peak2,
  deltaPhi       = phase_inner$deltaPhi.Est,
  prob_shift     = phase_inner$prob_shift,
  prob_conserved = phase_inner$prob_conserved,
  stringsAsFactors = FALSE
)

maintained_genes <- names(trans_outer$gain_loss_status[trans_outer$gain_loss_status == "Maintained"])

maintained_df <- maintained_df %>%
  filter(Gene %in% maintained_genes) %>%
  mutate(Peak_BaboonLUN_ZT = to_zt(Peak_BaboonLUN),
         Peak_HumanLUN_ZT  = to_zt(Peak_HumanLUN))

phase_class <- rep("Undetermined", length(gene_names))
names(phase_class) <- gene_names
phase_class[phase_inner$flag_cons]  <- "Phase-conserved"
phase_class[phase_inner$flag_shift] <- "Phase-shifted"
phase_class[phase_inner$flag_undetermined] <- "Undetermined"

maintained_df <- maintained_df %>%
  mutate(phase_class = phase_class[Gene]) %>%
  mutate(phase_class = ifelse(is.na(phase_class), "Undetermined", phase_class))

circadian_genes <- c(
  "BHLHE40","BHLHE41","BMAL1","BTRC","CLOCK","CREB1",
  "CRY1","CRY2","CSNK1D","CSNK1E","CUL1","DBP","FBXL3","FBXW11",
  "NFIL3","NPAS2","NR1D1","NR1D2","PER1","PER2","PER3",
  "RORA","RORB","RORC"
)

phase_colors <- c(
  "Phase-conserved" = "#1B9E77",
  "Phase-shifted"   = "#D95F02",
  "Undetermined"    = "#7570B3"
)

calculate_peak_difference <- function(a, b) {
  d <- abs(a - b)
  ifelse(d > 12, 24 - d, d)
}

maintained_df <- maintained_df %>%
  mutate(peak_diff = calculate_peak_difference(Peak_BaboonLUN_ZT, Peak_HumanLUN_ZT),
         within_concordance = peak_diff <= 2)

n_total  <- nrow(maintained_df)
n_within <- sum(maintained_df$within_concordance, na.rm = TRUE)
pct_within <- 100 * n_within / n_total

cat("\n[Baboon_Human_LUN: Within +/-2 h summary]\n")
cat("  Maintained genes:", n_total, "\n")
cat("  Within +/-2 h:", n_within, "\n")
cat(sprintf("  => %.1f%% within +/-2 h interval\n", pct_within))

# Diagnostic: which circadian genes are in the maintained set?
circ_in <- circadian_genes[circadian_genes %in% maintained_df$Gene]
circ_out <- circadian_genes[!circadian_genes %in% maintained_df$Gene]
cat("\n[Circadian genes in Maintained set]:", length(circ_in), "/", length(circadian_genes), "\n")
cat("  Present:", paste(circ_in, collapse = ", "), "\n")
if (length(circ_out) > 0) cat("  Missing:", paste(circ_out, collapse = ", "), "\n")

plot_df <- maintained_df %>%
  rowwise() %>%
  mutate(
    remapped = list(remap_to_diagonal(Peak_BaboonLUN_ZT, Peak_HumanLUN_ZT)),
    Peak_BaboonLUN_ZT_plot = remapped$x,
    Peak_HumanLUN_ZT_plot  = remapped$y
  ) %>%
  dplyr::select(-remapped) %>%
  ungroup()

n_Rc   <- nrow(plot_df)
pct_Rc <- round(100 * n_within / n_Rc, 1)

subtitle_text <- bquote(
  "Rhythmically Conserved Set " ~ R[c] ~
    "(" * n == .(n_Rc) * ", " * .(pct_Rc) * "% within " * "\u00B1" * "2 h interval)"
)

p <- ggplot(plot_df, aes(
  x = Peak_BaboonLUN_ZT_plot,
  y = Peak_HumanLUN_ZT_plot,
  color = phase_class
)) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1.0) +
  geom_abline(intercept = 2, slope = 1, color = "darkgreen", linetype = "dotted", linewidth = 1.1) +
  geom_abline(intercept = -2, slope = 1, color = "darkgreen", linetype = "dotted", linewidth = 1.1) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = phase_colors) +
  geom_text_repel(
    data = plot_df %>% filter(Gene %in% circadian_genes),
    aes(label = Gene),
    color = "black", fontface = "bold.italic", size = 3.5,
    segment.color = "gray50", box.padding = 0.2, point.padding = 0.1,
    min.segment.length = 0, force_pull = 1, force = 3,
    max.overlaps = Inf, max.time = 10, max.iter = 100000,
    seed = 42
  ) +
  labs(
    title    = "Circadian Peak Concordance: Baboon Lung versus Human Lung",
    subtitle = subtitle_text,
    x        = "Peak Hour - Baboon Lung (ZT)",
    y        = "Peak Hour - Human Lung (ZT)",
    color    = "Phase class"
  ) +
  scale_x_continuous(breaks = seq(-6, 18, 6), labels = sprintf("ZT%+d", seq(-6, 18, 6))) +
  scale_y_continuous(breaks = seq(-6, 18, 6), labels = sprintf("ZT%+d", seq(-6, 18, 6))) +
  coord_cartesian(xlim = c(-8, 20), ylim = c(-8, 20)) +
  theme_bw(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, margin = margin(b = 10)),
    axis.title    = element_text(face = "bold", size = 13),
    axis.text     = element_text(size = 12),
    legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(color = "gray70", fill = "white"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.margin = margin(t = 5, b = 5),
    plot.margin = margin(15, 15, 15, 15)
  )

save_dir <- file.path(output.dir, "figure/baboon_human_lung")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(save_dir, "Baboon_Human_LUN_Peak_Concordance_0.25_2h.pdf"),
       plot = p, width = 9, height = 8)
cat("\nFigure saved to:", file.path(save_dir, "Baboon_Human_LUN_Peak_Concordance_0.25_2h.pdf"), "\n")


# Regenerate concordance figure for Baboon_LUN_LIV
# Source in the lun_liv screen session where all objects are loaded

library(ggplot2); library(dplyr); library(ggrepel)
if (!dir.exists(tempdir())) dir.create(tempdir(), recursive = TRUE)

to_zt <- function(t_cos) ifelse(t_cos >= 18, t_cos - 24, t_cos)

remap_to_diagonal <- function(x, y, offset = 24) {
  x_options <- c(x, x - offset, x + offset)
  y_options <- c(y, y - offset, y + offset)
  best_dist <- Inf; best_x <- x; best_y <- y
  for (x_opt in x_options) {
    for (y_opt in y_options) {
      dist <- abs(y_opt - x_opt)
      if (dist < best_dist - 1e-9) { best_dist <- dist; best_x <- x_opt; best_y <- y_opt }
    }
  }
  return(data.frame(x = best_x, y = best_y, dist = best_dist))
}

gene_names <- names(phase_inner$peak1)

maintained_df <- data.frame(
  Gene = gene_names,
  Peak_LIV = phase_inner$peak1,
  Peak_LUN = phase_inner$peak2,
  deltaPhi = phase_inner$deltaPhi.Est,
  prob_shift = phase_inner$prob_shift,
  prob_conserved = phase_inner$prob_conserved,
  stringsAsFactors = FALSE
)

maintained_genes <- names(trans_outer$gain_loss_status[trans_outer$gain_loss_status == "Maintained"])

maintained_df <- maintained_df %>%
  filter(Gene %in% maintained_genes) %>%
  mutate(Peak_LUN_ZT = to_zt(Peak_LUN),
         Peak_LIV_ZT = to_zt(Peak_LIV))

phase_class <- rep("Undetermined", length(gene_names))
names(phase_class) <- gene_names
phase_class[phase_inner$flag_cons]  <- "Phase-conserved"
phase_class[phase_inner$flag_shift] <- "Phase-shifted"
phase_class[phase_inner$flag_undetermined] <- "Undetermined"

maintained_df <- maintained_df %>%
  mutate(phase_class = phase_class[Gene]) %>%
  mutate(phase_class = ifelse(is.na(phase_class), "Undetermined", phase_class))

circadian_genes <- c(
  "BHLHE40","BHLHE41","BMAL1","BTRC","CLOCK","CREB1",
  "CRY1","CRY2","CSNK1D","CSNK1E","CUL1","DBP","FBXL3","FBXW11",
  "NFIL3","NPAS2","NR1D1","NR1D2","PER1","PER2","PER3",
  "RORA","RORB","RORC"
)

phase_colors <- c(
  "Phase-conserved" = "#1B9E77",
  "Phase-shifted"   = "#D95F02",
  "Undetermined"    = "#7570B3"
)

calculate_peak_difference <- function(a, b) {
  d <- abs(a - b)
  ifelse(d > 12, 24 - d, d)
}

maintained_df <- maintained_df %>%
  mutate(peak_diff = calculate_peak_difference(Peak_LUN_ZT, Peak_LIV_ZT),
         within_concordance = peak_diff <= 2)

n_total  <- nrow(maintained_df)
n_within <- sum(maintained_df$within_concordance, na.rm = TRUE)
pct_within <- 100 * n_within / n_total

cat("\n[Baboon_LUN_LIV: Within +/-2 h summary]\n")
cat("  Maintained genes:", n_total, "\n")
cat("  Within +/-2 h:", n_within, "\n")
cat(sprintf("  => %.1f%% within +/-2 h interval\n", pct_within))

plot_df <- maintained_df %>%
  rowwise() %>%
  mutate(
    remapped = list(remap_to_diagonal(Peak_LUN_ZT, Peak_LIV_ZT)),
    Peak_LUN_ZT_plot = remapped$x,
    Peak_LIV_ZT_plot = remapped$y
  ) %>%
  dplyr::select(-remapped) %>%
  ungroup()

n_Rc   <- nrow(plot_df)
pct_Rc <- round(100 * n_within / n_Rc, 1)

subtitle_text <- bquote(
  "Rhythmically Conserved Set " ~ R[c] ~
    "(" * n == .(n_Rc) * ", " * .(pct_Rc) * "% within " * "\u00B1" * "2 h interval)"
)

p <- ggplot(plot_df, aes(
  x = Peak_LUN_ZT_plot,
  y = Peak_LIV_ZT_plot,
  color = phase_class
)) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1.0) +
  geom_abline(intercept = 2, slope = 1, color = "darkgreen", linetype = "dotted", linewidth = 1.1) +
  geom_abline(intercept = -2, slope = 1, color = "darkgreen", linetype = "dotted", linewidth = 1.1) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = phase_colors) +
  geom_text_repel(
    data = plot_df %>% filter(Gene %in% circadian_genes),
    aes(label = Gene),
    color = "black", fontface = "bold.italic", size = 4,
    segment.color = "gray50", box.padding = 1.2, point.padding = 1.5,
    min.segment.length = 0, force_pull = 0.3, max.overlaps = Inf
  ) +
  labs(
    title = "Circadian Peak Concordance: Baboon Lung versus Liver",
    subtitle = subtitle_text,
    x = "Peak Hour - Lung (ZT)",
    y = "Peak Hour - Liver (ZT)",
    color = "Phase class"
  ) +
  scale_x_continuous(breaks = seq(-6, 18, 6), labels = sprintf("ZT%+d", seq(-6, 18, 6))) +
  scale_y_continuous(breaks = seq(-6, 18, 6), labels = sprintf("ZT%+d", seq(-6, 18, 6))) +
  coord_cartesian(xlim = c(-8, 20), ylim = c(-8, 20)) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 12),
    legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(color = "gray70", fill = "white"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.margin = margin(t = 5, b = 5),
    plot.margin = margin(15, 15, 15, 15)
  )

save_dir <- file.path(output.dir, "figure/baboon_brain")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(filename = file.path(save_dir, "Baboon_LUN_LIV_Peak_Concordance_0.25_2h_new.pdf"),
       plot = p, width = 9, height = 8)
cat("\nFigure saved to:", file.path(save_dir, "Baboon_LUN_LIV_Peak_Concordance_0.25_2h_new.pdf"), "\n")


# Regenerate concordance figure for Baboon_SCN_HIP
# Source in the SCN screen session where all objects are loaded

library(ggplot2); library(dplyr); library(ggrepel)
if (!dir.exists(tempdir())) dir.create(tempdir(), recursive = TRUE)

to_zt <- function(t_cos) ifelse(t_cos >= 18, t_cos - 24, t_cos)

remap_to_diagonal <- function(x, y, offset = 24) {
  x_options <- c(x, x - offset, x + offset)
  y_options <- c(y, y - offset, y + offset)
  best_dist <- Inf; best_x <- x; best_y <- y
  for (x_opt in x_options) {
    for (y_opt in y_options) {
      dist <- abs(y_opt - x_opt)
      if (dist < best_dist - 1e-9) { best_dist <- dist; best_x <- x_opt; best_y <- y_opt }
    }
  }
  return(data.frame(x = best_x, y = best_y, dist = best_dist))
}

gene_names <- names(phase_inner$peak1)

maintained_df <- data.frame(
  Gene = gene_names,
  Peak_HIP = phase_inner$peak1,
  Peak_SCN = phase_inner$peak2,
  deltaPhi = phase_inner$deltaPhi.Est,
  prob_shift = phase_inner$prob_shift,
  prob_conserved = phase_inner$prob_conserved,
  stringsAsFactors = FALSE
)

maintained_genes <- names(trans_outer$gain_loss_status[trans_outer$gain_loss_status == "Maintained"])

maintained_df <- maintained_df %>%
  filter(Gene %in% maintained_genes) %>%
  mutate(Peak_SCN_ZT = to_zt(Peak_SCN),
         Peak_HIP_ZT = to_zt(Peak_HIP))

phase_class <- rep("Undetermined", length(gene_names))
names(phase_class) <- gene_names
phase_class[phase_inner$flag_cons]  <- "Phase-conserved"
phase_class[phase_inner$flag_shift] <- "Phase-shifted"
phase_class[phase_inner$flag_undetermined] <- "Undetermined"

maintained_df <- maintained_df %>%
  mutate(phase_class = phase_class[Gene]) %>%
  mutate(phase_class = ifelse(is.na(phase_class), "Undetermined", phase_class))

circadian_genes <- c(
  "BHLHE40","BHLHE41","BMAL1","BTRC","CLOCK","CREB1",
  "CRY1","CRY2","CSNK1D","CSNK1E","CUL1","DBP","FBXL3","FBXW11",
  "NFIL3","NPAS2","NR1D1","NR1D2","PER1","PER2","PER3",
  "RORA","RORB","RORC"
)

phase_colors <- c(
  "Phase-conserved" = "#1B9E77",
  "Phase-shifted"   = "#D95F02",
  "Undetermined"    = "#7570B3"
)

calculate_peak_difference <- function(a, b) {
  d <- abs(a - b)
  ifelse(d > 12, 24 - d, d)
}

maintained_df <- maintained_df %>%
  mutate(peak_diff = calculate_peak_difference(Peak_SCN_ZT, Peak_HIP_ZT),
         within_concordance = peak_diff <= 2)

n_total  <- nrow(maintained_df)
n_within <- sum(maintained_df$within_concordance, na.rm = TRUE)
pct_within <- 100 * n_within / n_total

cat("\n[Baboon_SCN_HIP: Within +/-2 h summary]\n")
cat("  Maintained genes:", n_total, "\n")
cat("  Within +/-2 h:", n_within, "\n")
cat(sprintf("  => %.1f%% within +/-2 h interval\n", pct_within))

plot_df <- maintained_df %>%
  rowwise() %>%
  mutate(
    remapped = list(remap_to_diagonal(Peak_SCN_ZT, Peak_HIP_ZT)),
    Peak_SCN_ZT_plot = remapped$x,
    Peak_HIP_ZT_plot = remapped$y
  ) %>%
  dplyr::select(-remapped) %>%
  ungroup()

n_Rc   <- nrow(plot_df)
pct_Rc <- round(100 * n_within / n_Rc, 1)

subtitle_text <- bquote(
  "Rhythmically Conserved Set " ~ R[c] ~
    "(" * n == .(n_Rc) * ", " * .(pct_Rc) * "% within " * "\u00B1" * "2 h interval)"
)

p <- ggplot(plot_df, aes(
  x = Peak_SCN_ZT_plot,
  y = Peak_HIP_ZT_plot,
  color = phase_class
)) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1.0) +
  geom_abline(intercept = 2, slope = 1, color = "darkgreen", linetype = "dotted", linewidth = 1.1) +
  geom_abline(intercept = -2, slope = 1, color = "darkgreen", linetype = "dotted", linewidth = 1.1) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = phase_colors) +
  geom_text_repel(
    data = plot_df %>% filter(Gene %in% circadian_genes),
    aes(label = Gene),
    color = "black", fontface = "bold.italic", size = 4,
    segment.color = "gray50", box.padding = 1.2, point.padding = 1.5,
    min.segment.length = 0, force_pull = 0.3, max.overlaps = Inf
  ) +
  labs(
    title = "Circadian Peak Concordance: Baboon SCN versus Hippocampus",
    subtitle = subtitle_text,
    x = "Peak Hour - SCN (ZT)",
    y = "Peak Hour - Hippocampus (ZT)",
    color = "Phase class"
  ) +
  scale_x_continuous(breaks = seq(-6, 18, 6), labels = sprintf("ZT%+d", seq(-6, 18, 6))) +
  scale_y_continuous(breaks = seq(-6, 18, 6), labels = sprintf("ZT%+d", seq(-6, 18, 6))) +
  coord_cartesian(xlim = c(-8, 20), ylim = c(-8, 20)) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5, margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 12),
    legend.position = "bottom", legend.direction = "horizontal", legend.box = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(color = "gray70", fill = "white"),
    legend.key = element_rect(fill = "white", color = NA),
    legend.margin = margin(t = 5, b = 5),
    plot.margin = margin(15, 15, 15, 15)
  )

save_dir <- file.path(output.dir, "figure/baboon_brain")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(filename = file.path(save_dir, "Baboon_SCN_HIP_Peak_Concordance_0.25_2h_new.pdf"),
       plot = p, width = 9, height = 8)
cat("\nFigure saved to:", file.path(save_dir, "Baboon_SCN_HIP_Peak_Concordance_0.25_2h_new.pdf"), "\n")


# Fixed concordance figure for Baboon SUN vs PUT
# Run this in the sun_put screen session where all objects are loaded

library(ggplot2)
library(dplyr)
library(ggrepel)

if (!dir.exists(tempdir())) dir.create(tempdir(), recursive = TRUE)

# Simple remap: minimize distance to diagonal via circular wrapping
remap_to_diagonal <- function(x, y, offset = 24) {
  x_options <- c(x, x - offset, x + offset)
  y_options <- c(y, y - offset, y + offset)

  best_dist <- Inf
  best_x <- x
  best_y <- y

  for (x_opt in x_options) {
    for (y_opt in y_options) {
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

# --- Rebuild maintained_df from phase_inner and trans_outer ---
to_zt <- function(t_cos) ifelse(t_cos >= 18, t_cos - 24, t_cos)

gene_names <- names(phase_inner$peak1)

maintained_df <- data.frame(
  Gene = gene_names,
  Peak_PUT = phase_inner$peak1,
  Peak_SUN = phase_inner$peak2,
  deltaPhi = phase_inner$deltaPhi.Est,
  prob_shift = phase_inner$prob_shift,
  prob_conserved = phase_inner$prob_conserved,
  stringsAsFactors = FALSE
)

maintained_genes <- names(trans_outer$gain_loss_status[trans_outer$gain_loss_status == "Maintained"])

maintained_df <- maintained_df %>%
  filter(Gene %in% maintained_genes) %>%
  mutate(Peak_SUN_ZT = to_zt(Peak_SUN),
         Peak_PUT_ZT = to_zt(Peak_PUT))

# Assign phase classification
phase_class <- rep("Undetermined", length(gene_names))
names(phase_class) <- gene_names
phase_class[phase_inner$flag_cons]  <- "Phase-conserved"
phase_class[phase_inner$flag_shift] <- "Phase-shifted"
phase_class[phase_inner$flag_undetermined] <- "Undetermined"

maintained_df <- maintained_df %>%
  mutate(phase_class = phase_class[Gene]) %>%
  mutate(phase_class = ifelse(is.na(phase_class), "Undetermined", phase_class))

# Circadian gene labels
circadian_genes <- c(
  "BHLHE40","BHLHE41","BMAL1","BTRC","CLOCK","CREB1",
  "CRY1","CRY2","CSNK1D","CSNK1E","CUL1","DBP","FBXL3","FBXW11",
  "NFIL3","NPAS2","NR1D1","NR1D2","PER1","PER2","PER3",
  "RORA","RORB","RORC"
)
genes_to_label <- maintained_df %>% filter(Gene %in% circadian_genes)

phase_colors <- c(
  "Phase-conserved" = "#1B9E77",
  "Phase-shifted"   = "#D95F02",
  "Undetermined"    = "#7570B3"
)

# Compute +/-2 h concordance
calculate_peak_difference <- function(a, b) {
  d <- abs(a - b)
  ifelse(d > 12, 24 - d, d)
}

maintained_df <- maintained_df %>%
  mutate(peak_diff = calculate_peak_difference(Peak_SUN_ZT, Peak_PUT_ZT),
         within_concordance = peak_diff <= 2)

n_total <- nrow(maintained_df)
n_within <- sum(maintained_df$within_concordance, na.rm = TRUE)
pct_within <- 100 * n_within / n_total

cat("\n[Within +/-2 h summary]\n")
cat("  Maintained genes:", n_total, "\n")
cat("  Within +/-2 h:", n_within, "\n")
cat(sprintf("  => %.1f%% within +/-2 h interval\n", pct_within))

# Remap for plotting
plot_df <- maintained_df %>%
  rowwise() %>%
  mutate(
    remapped = list(remap_to_diagonal(Peak_SUN_ZT, Peak_PUT_ZT)),
    Peak_SUN_ZT_plot = remapped$x,
    Peak_PUT_ZT_plot = remapped$y
  ) %>%
  select(-remapped) %>%
  ungroup()

# Diagnostic: check range of remapped values
cat("\n[Remapped ranges]\n")
cat("  SUN (x):", round(min(plot_df$Peak_SUN_ZT_plot), 2), "to", round(max(plot_df$Peak_SUN_ZT_plot), 2), "\n")
cat("  PUT (y):", round(min(plot_df$Peak_PUT_ZT_plot), 2), "to", round(max(plot_df$Peak_PUT_ZT_plot), 2), "\n")

n_Rc <- nrow(plot_df)
pct_Rc <- round(100 * n_within / n_Rc, 1)

subtitle_text <- bquote(
  "Rhythmically Conserved Set " ~ R[c] ~
    "(" * n == .(n_Rc) * ", " * .(pct_Rc) * "% within " * "\u00B1" * "2 h interval)"
)

p <- ggplot(plot_df, aes(
  x = Peak_SUN_ZT_plot,
  y = Peak_PUT_ZT_plot,
  color = phase_class
)) +
  geom_abline(intercept = 0, slope = 1,
              color = "black", linetype = "dashed", linewidth = 1.0) +
  geom_abline(intercept = 2, slope = 1,
              color = "darkgreen", linetype = "dotted", linewidth = 1.1) +
  geom_abline(intercept = -2, slope = 1,
              color = "darkgreen", linetype = "dotted", linewidth = 1.1) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = phase_colors) +
  geom_text_repel(
    data = plot_df %>% filter(Gene %in% circadian_genes),
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
    title = "Circadian Peak Concordance: Baboon Substantia Nigra versus Putamen",
    subtitle = subtitle_text,
    x = "Peak Hour - Substantia Nigra (ZT)",
    y = "Peak Hour - Putamen (ZT)",
    color = "Phase class"
  ) +
  # Use scale for breaks/labels only (no limits — let coord_cartesian handle clipping)
  scale_x_continuous(
    breaks = seq(-6, 18, 6),
    labels = sprintf("ZT%+d", seq(-6, 18, 6))
  ) +
  scale_y_continuous(
    breaks = seq(-6, 18, 6),
    labels = sprintf("ZT%+d", seq(-6, 18, 6))
  ) +
  # coord_cartesian zooms without removing data points
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

save_dir <- file.path(output.dir, "figure/baboon_brain")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = file.path(save_dir, "Baboon_SUN_PUT_Peak_Concordance_0.25_2h_new.pdf"),
  plot = p,
  width = 9,
  height = 8
)
cat("\nFigure saved to:", file.path(save_dir, "Baboon_SUN_PUT_Peak_Concordance_0.25_2h_new.pdf"), "\n")

