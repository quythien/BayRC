## Figure 6C: which circadian pathway genes are rhythmic in which tissue, with
## the tissues in the order of panel B's dendrogram.
##
## Usage: Rscript figure6_panelC.R [outdir]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
source(file.path(analysis.dir, "plots", "theme_bayrc.R"))
source(file.path(analysis.dir, "plots", "palette_concordance.R"))
suppressPackageStartupMessages({library(BayRC); library(ggplot2)})

args   <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else file.path(BAYRC_FIGURE_DIR, "figure6")
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

# genes ordered by the number of tissues calling them rhythmic
order_gene <- names(sort(tapply(cells$called, cells$gene, sum), decreasing = TRUE))
cells$gene   <- factor(cells$gene, levels = rev(order_gene))
cells$tissue <- factor(cells$tissue, levels = tissues)
cells$block  <- ifelse(cells$tissue %in% tight, "High-concordance cluster",
                       "Remaining tissues")

p <- ggplot(cells, aes(tissue, gene, fill = posterior)) +
  geom_tile(colour = "white", linewidth = .4) +
  # a ring marks the cells called rhythmic
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
  # full-width panel, so type is set to print at the size of the half-width ones
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

cairo_pdf(file.path(outdir, "Fig6C_circadian_membership.pdf"),
          width = 8.5, height = 5.4)
# column layout: this panel is scaled up more than the pair above it, so its
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

# the same panel for the two-by-two layout, on the canvas of the phase heatmap
# beside it so the two carry type at one size
cairo_pdf(file.path(outdir, "Fig6C_circadian_membership_nature.pdf"),
          width = 9.375, height = 8.6)
print(p +
      labs(subtitle = NULL) +
      # the key title is plain text rather than plotmath, so it takes the theme face
      guides(fill = guide_colourbar(position = "right", direction = "vertical",
               title = "Pr(\u03c1 = 1)",
               title.hjust = .5,
               barwidth = unit(6, "mm"), barheight = unit(62, "mm")),
             shape = guide_legend(position = "bottom", direction = "horizontal",
               title = NULL,
               override.aes = list(size = 4.2, stroke = 1.2, colour = "white",
                                   fill = NA))) +
      theme(# block names and body at the heights the phase panel uses
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
            # margins match the phase panel so the two bodies align
            plot.margin = margin(t = 7, r = 4, b = 5, l = 10)))
dev.off()
cat("Saving:", file.path(outdir, "Fig6C_circadian_membership.pdf"), "\n")

write.csv(cells, file.path(outdir, "Fig6C_circadian_membership.csv"),
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
