# Runnable illustration of the paired pathway profiles; values are synthetic.
# Usage from the package root:
# Rscript inst/analysis/readme_figures/plot_pathway_profiles_demo.R [output.pdf]
library(BayRC)
args <- commandArgs(trailingOnly = TRUE)
file <- if (length(args)) args[1] else "pathway_profiles_demo.pdf"
genes <- paste0("Gene", 1:8)
tissues <- paste0("Tissue", 1:8)
d <- expand.grid(gene = genes, tissue = tissues, stringsAsFactors = FALSE)
g <- match(d$gene, genes)
t <- match(d$tissue, tissues)
d$posterior <- ifelse(g == 8 & t != 1, 0.2, 0.98)
# Illustrative calls; use tissue-specific BFDR calls from all tested genes in real data.
d$called <- d$posterior > 0.9
d$peak <- (c(23, 0, 1, 7, 8, 14, 15, 20)[g] + ifelse(t == 4, 3, 0)) %% 24
d$resultant <- 0.9
d$interval_width <- 4
block <- rep(1:2, each = 4)
concordance <- outer(block, block, function(a, b) ifelse(a == b, 0.7, 0.1))
concordance[1:4, 1:4] <- 0.8
diag(concordance) <- 1
dimnames(concordance) <- list(tissues, tissues)
profiles <- plot_pathway_profiles(
  data = d, pathway_genes = genes, tissues = tissues,
  concordance = concordance, pathway_name = "Illustrative pathway profiles",
  tissue_k = 2, gene_k = 3, min_gene_fraction = 0.5, file = file
)
print(profiles$tissue_groups)
print(profiles$selected_tissues)
print(profiles$hierarchy$gene_groups)
