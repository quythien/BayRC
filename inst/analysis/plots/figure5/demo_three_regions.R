# Three baboon brain regions through the three-condition heatmap, written under
# FIG5_DIR/regions3. SCN is the reference; hippocampus and putamen are the
# comparators.

library(BayRC)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Paths come from config.R; override any of them with the matching env var.
this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

here <- Sys.getenv("FIG5_DIR", unset = file.path(BAYRC_FIGURE_DIR, "figure5"))
dir.create(here, recursive = TRUE, showWarnings = FALSE)
res  <- BAYRC_RESULT_DIR
bfdr_alpha <- 0.25
shift      <- 2

r <- new.env(); load(file.path(res, "summary", "hb", "mcmc_rho_BF3.RData"), envir = r)
p <- new.env(); load(file.path(res, "summary", "hb", "phi", "mcmc_phi_BF3.RData"), envir = p)
grab <- function(tis) list(rho = r$mcmc_data_baboon[[tis]],
                           phi = p$mcmc_phi_baboon[[tis]])
scn <- grab("SCN"); hip <- grab("HIP"); put <- grab("PUT")

kegg <- readRDS(system.file("extdata", "kegg_pathway_list_hsa.rds", package = "BayRC"))
kegg <- lapply(kegg, function(g) replace(g, g == "ARNTL", "BMAL1"))
genes <- intersect(kegg[["KEGG Circadian rhythm"]], rownames(scn$rho))
cat("pathway genes:", length(genes), "\n")

# the reference goes first, so gain and loss are the comparator's
pair <- function(ref, cmp) {
  tr <- transition_classify(rowMeans(ref$rho), rowMeans(cmp$rho),
                            bfdr_alpha = bfdr_alpha)
  ph <- phase_infer(phi_matrix1 = ref$phi, phi_matrix2 = cmp$phi,
                    gain_loss_status = tr$gain_loss_status,
                    bfdr_alpha = bfdr_alpha, shift = shift, P = 24,
                    compute_hdi = TRUE)
  list(trans = tr, phase = ph)
}
a <- pair(scn, hip); b <- pair(scn, put)

# rows run by which regions carry the rhythm, all three first and none last
called <- function(x) {
  q <- rowMeans(x$rho)
  bfdr_from_posterior(q, alpha = bfdr_alpha)$rhythmic_genes[match(genes, names(q))]
}
r1 <- called(scn); r2 <- called(hip); r3 <- called(put)
block <- ifelse(r1 & r2 & r3, 1, ifelse(r1 & r2 & !r3, 2,
         ifelse(r1 & !r2 & r3, 3, ifelse(!r1 & r2 & r3, 4,
         ifelse(r1 & !r2 & !r3, 5, ifelse(!r1 & r2 & !r3, 6,
         ifelse(!r1 & !r2 & r3, 7, 8)))))))
row_order <- genes[order(block, -rowMeans(scn$rho[genes, , drop = FALSE]))]
cat("genes per block:", paste(table(factor(block, levels = 1:8)), collapse = " "), "\n")

plot_heatmap(
  data1 = scn, data2 = hip, data3 = put,
  pathway_genes = genes, pathway_name = "KEGG Circadian rhythm",
  phase_results = a$phase, phase_results3 = b$phase,
  transition_results = a$trans, transition_results3 = b$trans,
  group_names  = c("Baboon SCN\n(Reference)", "Baboon HIP", "Baboon PUT"),
  legend_names = c("SCN", "HIP", "PUT"),
  row_order = row_order, canvas_width = 15.5,
  show_legend = TRUE, legend_side = "right",
  versions = "full", save_path = file.path(here, "regions3"))
