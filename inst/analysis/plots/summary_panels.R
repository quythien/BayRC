## Figure 2 panels A and B: rhythmicity transitions and phase classification,
## read from the plot caches the application scripts write.
##
## Usage: Rscript summary_panels.R [outdir]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
bayrc.needs.summary <- FALSE
source(file.path(analysis.dir, "config.R"))
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

args   <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else
  file.path(dirname(BAYRC_OUTPUT_DIR), "paper", "demos", "figure2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bfdr_alpha <- 0.25
shift      <- 2

## each pair, and the directory its application script writes
sources <- list(
  list(label = "PUT → SUN", dir = "baboon_PUT_SUN"),
  list(label = "PUT → VIC", dir = "baboon_PUT_VIC"),
  list(label = "SCN → HIP", dir = "baboon_SCN_HIP"))

read_pair <- function(s) {
  f <- file.path(BAYRC_FIGURE_DIR, s$dir, "plot_data.rds")
  if (!file.exists(f))
    stop("no plot cache for ", s$label, " at ", f,
         "\n  run applications/Baboon_", sub(" → ", "_", s$label), ".R first")
  cache <- readRDS(f)
  status <- cache$trans$gain_loss_status
  cls    <- cache$phase_class[cache$maintained]
  d <- (cache$phase$peak2 - cache$phase$peak1) %% 24
  d <- ifelse(d > 12, d - 24, d)
  # direction of a shift is the sign of peak2 - peak1; ahead when negative
  dm <- d[cache$maintained]
  sh <- cls == "Phase-shifted"
  list(label = s$label,
       trans = c(Gain = sum(status == "Gain"),
                 Conserved = sum(status == "Maintained"),
                 Loss = sum(status == "Loss")),
       phase = c(Aligned = sum(cls == "Phase-conserved"),
                 Ahead = sum(sh & dm < 0, na.rm = TRUE),
                 Behind = sum(sh & dm > 0, na.rm = TRUE),
                 Undetermined = sum(cls == "Undetermined")),
       offset = mean(dm, na.rm = TRUE))
}

pairs_data <- lapply(sources, read_pair)
pairs <- vapply(pairs_data, `[[`, character(1), "label")

trans <- do.call(rbind, lapply(pairs_data, function(p)
  data.frame(pair = p$label, status = names(p$trans), n = as.integer(p$trans))))
phase <- do.call(rbind, lapply(pairs_data, function(p)
  data.frame(pair = p$label, status = names(p$phase), n = as.integer(p$phase),
             total = sum(p$phase))))
phase$percent <- 100 * phase$n / phase$total

trans$pair   <- factor(trans$pair, levels = rev(pairs))
phase$pair   <- factor(phase$pair, levels = rev(pairs))
trans$status <- factor(trans$status, levels = c("Gain", "Conserved", "Loss"))
phase$status <- factor(phase$status,
                       levels = c("Aligned", "Ahead", "Behind", "Undetermined"))

common <- theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = 19),
        plot.subtitle = element_text(size = 14.5, colour = "black",
                                     margin = margin(b = 14)),
        legend.position = "bottom", legend.title = element_blank(),
        axis.title.x = element_text(face = "bold", size = 15),
        axis.text = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 14),
        axis.title.y = element_blank(), plot.margin = margin(10, 18, 10, 10),
        plot.tag = element_text(face = "bold", size = 18))

d <- ggplot(trans, aes(n, pair, fill = status)) +
  geom_col(position = position_dodge(width = .8), width = .7) +
  geom_text(aes(label = n), position = position_dodge(width = .8),
            hjust = -.2, size = 4.6) +
  scale_fill_manual(values = c(Gain = "#4393C3", Conserved = "#35978F",
                               Loss = "#D8A65D")) +
  scale_x_continuous(expand = expansion(mult = c(0, .13))) +
  labs(title = "Rhythmicity transitions",
       subtitle = sprintf("Genome-wide counts | BFDR = %.2f", bfdr_alpha),
       x = "Number of genes", tag = "A") + common

labels <- vapply(pairs_data, function(p)
  sprintf("%s\nAverage peak time difference:\n%+.2f h", p$label, p$offset),
  character(1))

e <- ggplot(phase, aes(percent, pair, fill = status)) +
  geom_col(width = .58, position = position_stack(reverse = TRUE)) +
  # segments under 4% are left unlabelled
  geom_text(aes(label = ifelse(percent >= 4, sprintf("%.0f%%", percent), "")),
            position = position_stack(vjust = .5, reverse = TRUE),
            size = 4.4, colour = "white") +
  scale_y_discrete(labels = setNames(labels, pairs)) +
  scale_x_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 100)) +
  # the two shifted classes share a hue
  scale_fill_manual(values = c(Aligned = "#1B9E77", Ahead = "#FDB863",
                               Behind = "#D95F02", Undetermined = "#8274B5"),
                    labels = c("Phase-conserved", "Second region peaks earlier",
                               "Second region peaks later", "Undetermined")) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = "Timing among conserved genes",
       subtitle = sprintf("Posterior phase classification | ±%g h window, BFDR = %.2f",
                          shift, bfdr_alpha),
       x = "Percentage of conserved genes", tag = "B") + common +
  # key anchored left and set smaller so it stays on the page
  theme(legend.justification = "left", legend.margin = margin(l = 0, r = 0),
        legend.text = element_text(size = 12.5),
        legend.key.width = unit(9, "pt"))

fig <- d + e + plot_layout(widths = c(1, 1.12))
# cairo carries the arrow in the pair labels, which the base pdf device drops
ggsave(file.path(outdir, "summary_panels.pdf"), fig, width = 15, height = 5.2,
       device = cairo_pdf)
write.csv(trans, file.path(outdir, "transition_counts.csv"), row.names = FALSE)
write.csv(phase, file.path(outdir, "phase_counts.csv"), row.names = FALSE)

cat("\ncounts read from the plot caches\n")
print(reshape(trans[, c("pair", "status", "n")], idvar = "pair",
              timevar = "status", direction = "wide"), row.names = FALSE)
print(reshape(phase[, c("pair", "status", "n")], idvar = "pair",
              timevar = "status", direction = "wide"), row.names = FALSE)
cat("\nmean peak difference over maintained genes:\n")
for (p in pairs_data) cat(sprintf("  %-10s %+.2f h\n", p$label, p$offset))
cat("\nwrote", outdir, "\n")
