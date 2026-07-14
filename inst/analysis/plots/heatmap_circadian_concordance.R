# Heatmap for KEGG Circadian rhythm with adjusted concordance score
# Expect baboon_LUN, human_LUN, kegg.pathway.list_hsa, phase_inner, trans_outer, output.dir in the environment

if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Package 'dplyr' is required.")
}

if (!exists("baboon_LUN") || !exists("human_LUN")) {
  stop("baboon_LUN and human_LUN must exist in the environment.")
}
if (!exists("kegg.pathway.list_hsa")) {
  stop("kegg.pathway.list_hsa must exist in the environment.")
}
if (!exists("output.dir")) {
  stop("output.dir must exist in the environment.")
}

pathway_name <- "KEGG Circadian rhythm"
if (!pathway_name %in% names(kegg.pathway.list_hsa)) {
  stop("KEGG Circadian rhythm not found in kegg.pathway.list_hsa.")
}

# Compute adjusted concordance score for the circadian pathway
result_circadian <- multi_conservation(
  mcmc.merge.list = list(Baboon_LUN = baboon_LUN, Human_LUN = human_LUN),
  dataset.names = c("Baboon_LUN", "Human_LUN"),
  select.pathway.list = kegg.pathway.list_hsa[pathway_name],
  n_perm = 1000,
  n_boot = 1000,
  output.dir = file.path(output.dir, "circadian_concordance"),
  use_cpp = TRUE
)

# Extract adjusted concordance (fallback to NA if not found)
adj_conc <- NA_real_
adj_col <- "Baboon_LUN_vs_Human_LUN_AdjustedConcordance"
if (adj_col %in% colnames(result_circadian)) {
  adj_conc <- result_circadian[[adj_col]][1]
}

subtitle_text <- if (is.na(adj_conc)) {
  "Adjusted Concordance: NA"
} else {
  sprintf("Adjusted Concordance: %.4f", adj_conc)
}

# Build heatmap using existing integrated plotting
plot_heatmap(
  data1 = baboon_LUN,
  data2 = human_LUN,
  pathway_genes = kegg.pathway.list_hsa[[pathway_name]],
  pathway_name = pathway_name,
  phase_results = phase_inner,
  transition_results = trans_outer,
  group_names = c("Baboon Lung", "Human Lung"),
  versions = "both",
  save_path = file.path(output.dir, "heatmap_circadian_concordance"),
  subtitle_override = subtitle_text
)
