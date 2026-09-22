## Figure 4: the enrichment on shared pathway rows, and what each enriched
## pathway is made of.
##
## Everything is read from stage2_significant.csv, which the application
## scripts write.
##
## Writes the whole of Figure 4; assemble_figures.R does not touch it.
##
## Usage: Rscript pathway_transition_panels.R [outdir]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
bayrc.needs.summary <- FALSE
source(file.path(analysis.dir, "config.R"))
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})

args   <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else
  file.path(dirname(BAYRC_OUTPUT_DIR), "paper", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

stage2_q <- 0.05
pairs <- list(
  list(label = "Putamen versus substantia nigra", dir = "baboon_PUT_SUN"),
  list(label = "Putamen versus visual cortex",    dir = "baboon_PUT_VIC"))

read_pair <- function(p) {
  sig <- read.csv(file.path(BAYRC_FIGURE_DIR, p$dir, "stage2_significant.csv"))
  sig <- sig[sig$q < stage2_q, ]
  sig$comparison <- p$label
  sig
}

sig <- do.call(rbind, lapply(pairs, read_pair))
sig$pathway <- sub("^KEGG ", "", sig$pathway)
sig$direction <- factor(sig$direction, levels = c("gain", "loss", "conserved"),
                        labels = c("Gain", "Loss", "Conserved"))
sig$comparison <- factor(sig$comparison,
                         levels = vapply(pairs, `[[`, character(1), "label"))

## one row order for both panels, strongest first
best <- tapply(sig$q, sig$pathway, min)
lev  <- names(sort(best, decreasing = TRUE))
sig$pathway <- factor(sig$pathway, levels = lev)

ramp <- c("#DFF2F4", "#AFDEC2", "#7AC6B9", "#CF9166", "#C06C84")
common <- theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(plot.subtitle = element_text(size = 14.5, colour = "black",
                                     margin = margin(b = 10)),
        axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(face = "bold", size = 15),
        strip.text = element_text(face = "bold", size = 14.5),
        strip.background = element_blank(),
        legend.text = element_text(size = 13), legend.title = element_text(size = 14),
        plot.tag = element_text(face = "bold", size = 18))

## A. both comparisons on the same pathway rows
sig$expected <- with(sig, ifelse(direction == "Gain", Expected_N_Gain,
                          ifelse(direction == "Loss", Expected_N_Loss,
                                 Expected_N_Conserved)))
a <- ggplot(sig, aes(direction, pathway)) +
  geom_point(aes(size = expected, colour = -log10(q))) +
  facet_wrap(~comparison) +
  scale_x_discrete(drop = FALSE) +
  scale_colour_gradientn(colours = ramp,
                         name = expression(-log[10](q))) +
  scale_size_area(max_size = 11, name = "expected gene count") +
  labs(x = "Transition", y = NULL, tag = "A") + common

## B. what the enrichment is made of, with a star on the tested transition
comp <- do.call(rbind, lapply(split(sig, list(sig$comparison, sig$pathway),
                                    drop = TRUE), function(s) {
  r <- s[1, ]
  tot <- r$Expected_N_Gain + r$Expected_N_Loss + r$Expected_N_Conserved
  data.frame(comparison = r$comparison, pathway = r$pathway,
             status = c("Gain", "Loss", "Conserved"),
             frac = c(r$Expected_N_Gain, r$Expected_N_Loss,
                      r$Expected_N_Conserved) / tot,
             starred = c("Gain", "Loss", "Conserved") %in% s$direction)
}))
comp$pathway    <- factor(comp$pathway, levels = lev)
comp$status     <- factor(comp$status, levels = c("Gain", "Loss", "Conserved"))
comp$comparison <- factor(comp$comparison, levels = levels(sig$comparison))

## each star sits at the middle of its segment, in stacking order
comp <- comp[order(comp$comparison, comp$pathway, comp$status), ]
comp$mid <- unlist(lapply(split(comp$frac, list(comp$comparison, comp$pathway),
                                drop = TRUE),
                          function(f) 100 * (cumsum(f) - f / 2)))

b <- ggplot(comp, aes(100 * frac, pathway, fill = status)) +
  geom_col(width = .68, position = position_stack(reverse = TRUE)) +
  geom_text(data = comp[comp$starred, ], aes(x = mid, label = "*"),
            colour = "white", size = 9, vjust = .72, show.legend = FALSE) +
  facet_wrap(~comparison) +
  scale_fill_manual(values = c(Gain = "#4393C3", Loss = "#D8A65D",
                               Conserved = "#35978F"), name = NULL) +
  scale_x_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
  labs(subtitle = sprintf("* enriched for that transition at q < %.2f", stage2_q),
       x = "Percentage of expected rhythmic genes", y = NULL, tag = "B") + common

fig <- a / b + plot_layout(heights = c(1, 1))
ggsave(file.path(outdir, "Figure_4.pdf"), fig, width = 15, height = 12)
ggsave(file.path(outdir, "Figure_4.png"), fig, width = 15, height = 12,
       dpi = 200, bg = "white")

cat("\nenriched rows:", nrow(sig), "over", nlevels(droplevels(sig$pathway)),
    "pathways\n")
print(sig[order(sig$comparison, sig$q),
          c("comparison", "pathway", "direction", "q")], row.names = FALSE)
cat("wrote", outdir, "\n")
