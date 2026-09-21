# Peak-phase concordance scatter shared by the application scripts.

library(ggplot2)
library(ggrepel)

phase_colors <- c("Phase-conserved" = "#1B9E77",
                  "Phase-shifted"   = "#D95F02",
                  "Undetermined"    = "#7570B3")

to_zt <- function(t) ifelse(t >= 18, t - 24, t)

# circular difference in hours, taking the shorter way round
circ_diff <- function(a, b, P = 24) {
  d <- abs(a - b)
  ifelse(d > P / 2, P - d, d)
}

# shift y by whole periods so each point sits nearest the diagonal
remap_to_diagonal <- function(x, y, P = 24) y + P * round((x - y) / P)

peak_concordance_plot <- function(peak_x, peak_y, phase_class, label_genes,
                                  title, xlab, ylab, window = 2) {
  d <- data.frame(Gene = names(peak_x),
                  x = to_zt(peak_x),
                  y = to_zt(peak_y),
                  phase_class = phase_class[names(peak_x)],
                  stringsAsFactors = FALSE)
  d$y <- remap_to_diagonal(d$x, d$y)
  d$phase_class[is.na(d$phase_class)] <- "Undetermined"

  n_within <- sum(circ_diff(peak_x, peak_y) <= window)
  subtitle <- bquote("Rhythmically Conserved Set " ~ R[c] ~
                     "(" * n == .(nrow(d)) * ", " *
                     .(round(100 * n_within / nrow(d), 1)) *
                     "% within " * "±" * .(window) * " h interval)")

  ggplot(d, aes(x = x, y = y, color = phase_class)) +
    geom_abline(intercept = 0, slope = 1, color = "black",
                linetype = "dashed", linewidth = 1.0) +
    geom_abline(intercept = window, slope = 1, color = "darkgreen",
                linetype = "dotted", linewidth = 1.1) +
    geom_abline(intercept = -window, slope = 1, color = "darkgreen",
                linetype = "dotted", linewidth = 1.1) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(values = phase_colors) +
    # the repel search starts from a random layout, so it is seeded to make
    # every rebuild of a panel place its labels the same way
    geom_text_repel(data = d[d$Gene %in% label_genes, ], aes(label = Gene),
                    color = "black", fontface = "bold.italic", size = 4,
                    segment.color = "gray50", box.padding = 1.2,
                    point.padding = 1.5, min.segment.length = 0,
                    force_pull = 0.3, max.overlaps = Inf, seed = 1) +
    labs(title = title, subtitle = subtitle, x = xlab, y = ylab,
         color = "Phase class") +
    scale_x_continuous(breaks = seq(-6, 18, 6),
                       labels = sprintf("ZT%+d", seq(-6, 18, 6))) +
    scale_y_continuous(breaks = seq(-6, 18, 6),
                       labels = sprintf("ZT%+d", seq(-6, 18, 6))) +
    coord_cartesian(xlim = c(-8, 20), ylim = c(-8, 20)) +
    theme_bayrc(base_size = 14) +
    theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          plot.subtitle = element_text(size = 13, hjust = 0.5,
                                       margin = margin(b = 10)),
          axis.title = element_text(face = "bold", size = 13),
          axis.text = element_text(size = 12),
          legend.position = "bottom",
          legend.title = element_text(face = "bold"),
          legend.background = element_rect(color = "gray70", fill = "white"),
          legend.key = element_rect(fill = "white", color = NA),
          plot.margin = margin(15, 15, 15, 15))
}
