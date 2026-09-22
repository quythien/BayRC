## Figure 2C: when the core clock peaks in each tissue, with the tissues split
## by the two-cluster cut of panel B's dendrogram. Each ring is one clock gene
## and each point one tissue at its posterior peak time, with the 95% circular
## HDI drawn as an arc.
##
## Usage: Rscript clock_phase_dial.R [outdir]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
source(file.path(analysis.dir, "plots", "theme_bayrc.R"))
suppressPackageStartupMessages({library(BayRC); library(ggplot2)})

args   <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else file.path(BAYRC_FIGURE_DIR, "figure2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

bfdr_alpha <- 0.25
# loop order: activators, repressors, nuclear receptors, output
clock <- c("BMAL1", "CLOCK", "NPAS2", "PER1", "PER2", "CRY1", "CRY2",
           "NR1D1", "NR1D2", "DBP")

# circular mean and sd in hours, from the resultant of the unit vectors
circ_mean <- function(x, P = 24) {
  a <- 2 * pi * x / P
  (atan2(mean(sin(a)), mean(cos(a))) %% (2 * pi)) * P / (2 * pi)
}
circ_sd <- function(x, P = 24) {
  a <- 2 * pi * x / P
  R <- sqrt(mean(cos(a))^2 + mean(sin(a))^2)
  sqrt(-2 * log(max(R, 1e-12))) * P / (2 * pi)
}

conc <- file.path(BAYRC_OUTPUT_DIR, "heatmap_circadian_pairs", "within_baboon",
                  "pairwise_concordance_baboon_circadian_Matrix.csv")
if (!file.exists(conc))
  stop("no circadian concordance matrix at ", conc,
       "\n  run plots/heatmap_circadian_pairs.R within_baboon first")

m <- as.matrix(read.csv(conc, row.names = 1, check.names = FALSE))
diag(m) <- 1
groups <- cutree(hclust(as.dist(1 - m), method = "ward.D2"), k = 2)

# the tight cluster is the one with the higher mean within-cluster concordance
within <- vapply(sort(unique(groups)), function(g) {
  s <- m[groups == g, groups == g, drop = FALSE]
  mean(s[row(s) != col(s)])
}, numeric(1))
tight <- names(groups)[groups == which.max(within)]
cat("tight cluster (mean within-cluster c-score", sprintf("%.3f", max(within)),
    "):\n  ", paste(sort(tight), collapse = " "), "\n")
cat("remaining tissues (", sprintf("%.3f", min(within)), "):\n  ",
    paste(sort(setdiff(names(groups), tight)), collapse = " "), "\n")

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))
rho <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
phi <- get(grep("^mcmc_phi",  ls(), value = TRUE)[1])

tissues <- intersect(names(groups), names(rho))
clock   <- clock[clock %in% rownames(rho[[tissues[1]]])]
cat("clock genes measured:", paste(clock, collapse = " "), "\n")

# a gene enters the dial for a tissue only where it is called rhythmic there
peaks <- do.call(rbind, lapply(tissues, function(t) {
  p  <- rowMeans(rho[[t]])
  tau <- bfdr_from_posterior(p, alpha = bfdr_alpha)$threshold
  do.call(rbind, lapply(clock, function(g) {
    if (p[[g]] < tau) return(NULL)
    d <- phi[[t]][g, ]
    h <- circular_HDI(d, credMass = 0.95, P = 24)
    data.frame(tissue = t, gene = g,
               peak  = circ_mean(d),
               lower = h$lower, upper = h$upper,
               cluster = if (t %in% tight) "High-concordance cluster" else
                         "Remaining tissues")
  }))
}))

peaks$gene    <- factor(peaks$gene, levels = rev(clock))
peaks$cluster <- factor(peaks$cluster,
                        levels = c("High-concordance cluster", "Remaining tissues"))
peaks$ring    <- as.integer(peaks$gene)

# each HDI as a run of points from its lower to its upper bound, across the seam
arcs <- do.call(rbind, lapply(seq_len(nrow(peaks)), function(i) {
  r <- peaks[i, ]
  span <- (r$upper - r$lower) %% 24
  data.frame(id = i, gene = r$gene, ring = r$ring, cluster = r$cluster,
             peak = (r$lower + seq(0, span, length.out = 24)) %% 24)
}))

# the cluster is drawn just outside each ring and the rest just inside it
peaks$r <- peaks$ring + ifelse(peaks$cluster == "High-concordance cluster", .18, -.18)
arcs$r  <- arcs$ring  + ifelse(arcs$cluster  == "High-concordance cluster", .18, -.18)

cols <- c("High-concordance cluster" = "#1372B1", "Remaining tissues" = "#BEBEBE")

p <- ggplot() +
  geom_path(data = arcs, aes(peak, r, group = id, colour = cluster),
            linewidth = .7, alpha = .5, lineend = "round") +
  geom_point(data = peaks, aes(peak, r, colour = cluster), size = 2.1) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_continuous(limits = c(0, 24), breaks = seq(0, 18, 6),
                     labels = sprintf("ZT%d", seq(0, 18, 6))) +
  scale_y_continuous(breaks = seq_along(levels(peaks$gene)),
                     labels = levels(peaks$gene),
                     limits = c(-1.2, length(levels(peaks$gene)) + .8)) +
  coord_polar(theta = "x", start = 0) +
  labs(title = "Core clock peak times by tissue",
       subtitle = sprintf("Posterior peak and 95%% circular HDI, BFDR = %.2f",
                          bfdr_alpha),
       x = NULL, y = NULL) +
  theme_bayrc(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = .5),
        plot.subtitle = element_text(size = 12.5, hjust = .5,
                                     margin = margin(b = 6)),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10.5, face = "italic"),
        panel.grid.major.y = element_line(colour = "grey88", linewidth = .3),
        panel.grid.major.x = element_line(colour = "grey88", linewidth = .3),
        legend.position = "bottom")

cairo_pdf(file.path(outdir, "Fig2C_clock_phase_dial.pdf"), width = 7, height = 7)
print(p)
dev.off()
cat("Saving:", file.path(outdir, "Fig2C_clock_phase_dial.pdf"), "\n")

# circular sd of peak times across tissues, per gene and group
spread <- do.call(rbind, lapply(levels(peaks$cluster), function(g) {
  do.call(rbind, lapply(clock, function(k) {
    d <- peaks$peak[peaks$cluster == g & peaks$gene == k]
    if (!length(d)) return(NULL)
    data.frame(cluster = g, gene = k, n = length(d),
               sd_hours = round(circ_sd(d), 2))
  }))
}))
write.csv(spread, file.path(outdir, "Fig2C_clock_phase_spread.csv"),
          row.names = FALSE)
cat("\nspread of peak times across tissues, hours\n")
print(reshape(spread[, c("gene", "cluster", "sd_hours")], idvar = "gene",
              timevar = "cluster", direction = "wide"), row.names = FALSE)
