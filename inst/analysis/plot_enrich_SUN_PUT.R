# Enrichment for phase-shifted genes: Baboon SUN vs PUT
# Requires phase_inner and trans_outer in the current session

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

if (!exists("phase_inner") || !exists("trans_outer")) {
  stop("phase_inner and trans_outer must exist in the session before running this script.")
}

output.dir <- "/home/qtp1/Projects/Collaborative/Paper/Congruence/PNAS_aging/results/baboon/output_final"
plot_dir <- file.path(output.dir, "figure/baboon_brain")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Shifted genes among maintained set
flag_shift <- phase_inner$flag_shift
gene_names <- names(flag_shift)
if (is.null(gene_names)) {
  gene_names <- names(phase_inner$peak1)
}
shifted_genes <- gene_names[flag_shift]
shifted_genes <- shifted_genes[!is.na(shifted_genes)]

if (length(shifted_genes) == 0) {
  stop("No shifted genes found.")
}

# Map symbols to Entrez IDs
entrez <- bitr(shifted_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
if (is.null(entrez) || nrow(entrez) == 0) {
  stop("No Entrez IDs mapped from shifted genes.")
}

gene_ids <- unique(entrez$ENTREZID)

# GO enrichment (BP)
ego_bp <- enrichGO(
  gene = gene_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.2,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# GO enrichment (MF)
ego_mf <- enrichGO(
  gene = gene_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.2,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# KEGG enrichment
# Use human KEGG (hsa) since symbols are human orthologs
kegg <- enrichKEGG(
  gene = gene_ids,
  organism = "hsa",
  pvalueCutoff = 0.2,
  qvalueCutoff = 0.2
)

# Plot helpers
wrap_terms <- function(x, width = 45) {
  ifelse(nchar(x) > width, sapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n")), x)
}

save_dot_bar <- function(enr, prefix, title_text, show_n = 15) {
  if (is.null(enr) || nrow(as.data.frame(enr)) == 0) {
    cat("No enrichment for", prefix, "\n")
    return()
  }
  if (inherits(enr, "enrichResult") && !is.null(enr@result$Description)) {
    enr@result$Description <- wrap_terms(enr@result$Description, width = 45)
  }

  g1 <- dotplot(enr, showCategory = show_n) +
    ggtitle(title_text) +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5)
    )
  g2 <- barplot(enr, showCategory = show_n) +
    ggtitle(title_text) +
    theme(
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5)
    )

  ggsave(file.path(plot_dir, paste0(prefix, "_dotplot.pdf")), g1, width = 9, height = 8)
  ggsave(file.path(plot_dir, paste0(prefix, "_barplot.pdf")), g2, width = 9, height = 8)

  # Save table
  write.csv(as.data.frame(enr), file.path(plot_dir, paste0(prefix, ".csv")), row.names = FALSE)
}

save_dot_bar(ego_bp, "SUN_PUT_shifted_GO_BP", "GO Biological Process: Phase-shifted genes (SUN vs PUT)", show_n = 15)
save_dot_bar(ego_mf, "SUN_PUT_shifted_GO_MF", "GO Molecular Function: Phase-shifted genes (SUN vs PUT)", show_n = 15)
save_dot_bar(kegg,  "SUN_PUT_shifted_KEGG",  "KEGG: Phase-shifted genes (SUN vs PUT)", show_n = 15)

cat("\nAll plots saved to:", plot_dir, "\n")
