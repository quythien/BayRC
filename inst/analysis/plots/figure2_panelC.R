## Figure 2C: which circadian pathway genes are rhythmic in which tissue.
##
## The adjusted c-score in panel B compares the identity of the rhythmic genes
## in a pair, not their timing and not how many there are. This panel shows that
## identity gene by gene, with the tissues in the order panel B's own dendrogram
## gives, so the block the cluster forms is visible rather than asserted.
##
## Usage: Rscript figure2_panelC.R [outdir]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
source(file.path(analysis.dir, "plots", "theme_bayrc.R"))
source(file.path(analysis.dir, "plots", "palette_concordance.R"))
suppressPackageStartupMessages({library(BayRC); library(ggplot2)})

args   <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else file.path(BAYRC_FIGURE_DIR, "figure2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bfdr_alpha <- 0.25
pathway    <- "KEGG Circadian rhythm"

conc <- file.path(BAYRC_OUTPUT_DIR, "heatmap_circadian_pairs", "within_baboon",
                  "pairwise_concordance_baboon_circadian_Matrix.csv")
if (!file.exists(conc))
  stop("no circadian concordance matrix at ", conc,
       "\n  run plots/heatmap_circadian_pairs.R within_baboon first")

m <- as.matrix(read.csv(conc, row.names = 1, check.names = FALSE))
diag(m) <- 1
hc     <- hclust(as.dist(1 - m), method = "ward.D2")
groups <- cutree(hc, k = 2)
within <- vapply(sort(unique(groups)), function(g) {
  s <- m[groups == g, groups == g, drop = FALSE]
  mean(s[row(s) != col(s)])
}, numeric(1))
tight <- names(groups)[groups == which.max(within)]

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
rho <- get(grep("^mcmc_data", ls(), value = TRUE)[1])

kegg <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))
kegg <- lapply(kegg, function(g) replace(g, g == "ARNTL", "BMAL1"))

tissues <- intersect(hc$labels[hc$order], names(rho))
genes   <- intersect(kegg[[pathway]], rownames(rho[[tissues[1]]]))

cells <- do.call(rbind, lapply(tissues, function(t) {
  p   <- rowMeans(rho[[t]])
  tau <- bfdr_from_posterior(p, alpha = bfdr_alpha)$threshold
  data.frame(tissue = t, gene = genes, posterior = as.numeric(p[genes]),
             called = as.numeric(p[genes]) >= tau)
}))

# genes down the panel in order of how widely they are rhythmic, tissues across
# it in the order panel B clusters them
order_gene <- names(sort(tapply(cells$called, cells$gene, sum), decreasing = TRUE))
cells$gene   <- factor(cells$gene, levels = rev(order_gene))
cells$tissue <- factor(cells$tissue, levels = tissues)
cells$block  <- ifelse(cells$tissue %in% tight, "High-concordance cluster",
                       "Remaining tissues")

p <- ggplot(cells, aes(tissue, gene, fill = posterior)) +
  geom_tile(colour = "white", linewidth = .4) +
  # a ring on the cells that clear the threshold, so the panel carries the call
  # as well as the posterior behind it
  geom_point(data = cells[cells$called, ],
             aes(shape = "Called rhythmic at BFDR"), size = 1.15,
             fill = NA, colour = "white", stroke = .55) +
  scale_shape_manual(name = NULL, values = c(`Called rhythmic at BFDR` = 21),
                     labels = sprintf("Called rhythmic at BFDR = %.2f", bfdr_alpha)) +
  facet_grid(~ block, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(colours = concordance_colors, limits = c(0, 1),
                       name = expression(Pr(rho == 1))) +
  guides(shape = "none") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "Rhythmic membership of the circadian pathway",
       subtitle = sprintf("Posterior probability of rhythmicity; ○ marks genes called rhythmic at BFDR = %.2f",
                          bfdr_alpha),
       x = NULL, y = NULL) +
  theme_bayrc(base_size = 14) +
  # this panel spans the figure while A and B take half of it each, so its type
  # is set smaller here to print at the same size as theirs
  theme(plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(size = 14, margin = margin(b = 7)),
        axis.text.x = element_text(size = 15, angle = 90, vjust = .5, hjust = 1),
        axis.text.y = element_text(size = 15, face = "italic"),
        strip.text = element_text(face = "bold", size = 15),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(6, "pt"),
        legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 14),
        legend.key.height = unit(13, "mm"))

cairo_pdf(file.path(outdir, "Fig2C_circadian_membership.pdf"),
          width = 8.5, height = 5.4)
# the column layout gives this panel the whole width of a column that holds two
# panels side by side above it, so it is scaled up more than they are and its
# type is set down by the same ratio
print(p +
      theme(plot.title = element_text(size = 12.5),
            plot.subtitle = element_text(size = 9.5),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            strip.text = element_text(size = 10),
            # a block name is wider than the block it sits over
            strip.clip = "off",
            legend.title = element_text(size = 10.5),
            legend.text = element_text(size = 9.5),
            legend.key.width = unit(4, "mm"),
            legend.key.height = unit(22, "mm")))
dev.off()

# The same panel for the two-by-two layout, where it sits beside the phase
# heatmap. It is drawn on that panel's canvas so the two carry type at one size,
# and the extra height gives each gene name a row taller than the name itself.
cairo_pdf(file.path(outdir, "Fig2C_circadian_membership_nature.pdf"),
          width = 9.375, height = 8.6)
print(p +
      # a blank name gives the ring key the title row the bar has, so the two
      # keys stand on one line
      # plotmath sets its own face, so the name is written out and left to the
      # theme, which is what the panel beside it does
      labs(subtitle = NULL) +
      # plotmath sets its own face, so the name is written out here and left
      # to the theme, which is what the panel beside it does
      guides(fill = guide_colourbar(position = "right", direction = "vertical",
               title = "Pr(\u03c1 = 1)",
               title.hjust = .5,
               barwidth = unit(6, "mm"), barheight = unit(62, "mm")),
             shape = guide_legend(position = "bottom", direction = "horizontal",
               title = NULL,
               override.aes = list(size = 4.2, stroke = 1.2, colour = "white",
                                   fill = NA))) +
      theme(# the block names and the body are set to the heights the phase
            # panel puts them at, so the two read across
            plot.title = element_text(face = "bold", size = 24,
                                      margin = margin(b = 14)),
            strip.text = element_text(face = "bold", size = 17,
                                      margin = margin(t = 0, b = 5.4)),
            axis.text.x = element_text(size = 15),
            axis.text.y = element_text(size = 15),
            legend.position = "bottom",
            legend.box.just = "top",
            legend.key.height = unit(5, "mm"),
            legend.key.width = unit(26, "mm"),
            legend.title = element_text(size = 16, face = "bold",
                                        margin = margin(b = 26)),
            legend.text = element_text(size = 15),
            # the ring is white, so its key needs a dark tile and an edge
            legend.key = element_rect(fill = concordance_colors[185],
                                      colour = "black", linewidth = .4),
            legend.title.position = "top",
            axis.ticks = element_blank(),
            panel.spacing = unit(8.5, "pt"),
            legend.margin = margin(t = 6, b = 34),
            # a block name is wider than the block it sits over
            strip.clip = "off",
            # the phase panel keeps a right-hand column for its group names, so
            # the same width is held back here and the two bodies stand over
            # each other
            plot.margin = margin(t = 7, r = 4, b = 5, l = 10)))
dev.off()
cat("Saving:", file.path(outdir, "Fig2C_circadian_membership.pdf"), "\n")

write.csv(cells, file.path(outdir, "Fig2C_circadian_membership.csv"),
          row.names = FALSE)

shared <- tapply(cells$called, cells$block, function(x) x)
cat("\ngenes called rhythmic, mean per tissue\n")
for (b in unique(cells$block))
  cat(sprintf("  %-26s %.1f of %d\n", b,
              mean(tapply(cells$called[cells$block == b],
                          droplevels(cells$tissue[cells$block == b]), sum)),
              length(genes)))

# how often a pair of tissues in the same block agrees on a gene's call
agreement <- function(ts) {
  k <- vapply(ts, function(t) cells$called[cells$tissue == t][
    order(cells$gene[cells$tissue == t])], logical(length(genes)))
  pr <- combn(ncol(k), 2)
  mean(vapply(seq_len(ncol(pr)), function(i)
    mean(k[, pr[1, i]] == k[, pr[2, i]]), numeric(1)))
}
cat("\npairs in the same block agreeing on a gene's call\n")
cat(sprintf("  cluster %.1f%%   rest %.1f%%\n",
            100 * agreement(tight), 100 * agreement(setdiff(tissues, tight))))
