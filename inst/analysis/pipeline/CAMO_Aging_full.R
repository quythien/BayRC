# Pathway conservation test (B = 1e6) for younger against older, then per-pathway
# rhythmicity heatmaps. data/mcmc_young_old.rds holds match_symbols() output for
# the COMBINED_younger and COMBINED_older MCMC runs.

rm(list = ls())
gc()

# Paths come from config.R; override any of them with the matching env var.
bayrc.needs.summary <- FALSE
this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
current_aging <- BAYRC_AGING_DIR

mcmc_age = readRDS(file.path(current_aging, "data/mcmc_young_old.rds"))


library(openxlsx)

output.dir = file.path(current_aging, "results/brain_regions/Conservation_Pathway_Test_1e6")  
if (!dir.exists(output.dir)) {
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
}

timing <- system.time({
result <- multi_conservation_pathway(
  mcmc.merge.list = mcmc_age,
  dataset.names = c("younger", "older"),
  select.pathway.list = kegg.pathway.list_hsa,
  B = 1000000,
  rounds = 5000,                    
  save_intermediate = FALSE,        
  intermediate_dir = "Be6_intermediate",
  ncores = 30,                     
  parallel = TRUE,# <= 125 GB
  output.dir = output.dir
)
})
# 2000 rounds take about 200 GB
print(timing)
print(output.dir)

cat("Total elapsed time:", timing["elapsed"], "seconds\n")
cat("That's", round(timing["elapsed"]/60, 1), "minutes\n") # 3.68 hours 

# Top 5 pathways by each conservation measure
as.data.frame(lapply(conservation_results$Conservation_results, function(x) {
  values <- as.vector(x)
  names(values) <- rownames(x)

  head(sort(values, decreasing = TRUE), 5)
})$Cp)

conservation_results$PValue_results$congruence_index["KEGG Jak-STAT signaling pathway", 1]
conservation_results$Conservation_results$congruence_index["KEGG Jak-STAT signaling pathway", 1]
###########
# Heatmap of log Bayes factor per pathway gene, annotated by rhythmicity concordance
plot_pathway_concordance <- function(younger_data, older_data, pathway_genes,
                                     pathway_name, BF_threshold = 3, p_rhythmic = 0.20,
                                     save_path = NULL) {

  library(pheatmap)

  all_genes <- rownames(younger_data$rho)
  overlap_genes <- intersect(all_genes, pathway_genes)

  # Posterior odds over prior odds of rhythmicity
  calculate_BF <- function(rho_matrix, p_rhythmic = 0.20) {
    rho_mean <- rowMeans(rho_matrix, na.rm = TRUE)
    BF <- (rho_mean * (1 - p_rhythmic)) / ((1 - rho_mean + 1e-10) * p_rhythmic)
    return(BF)
  }

  younger_BF <- calculate_BF(younger_data$rho[overlap_genes, , drop = FALSE], p_rhythmic)
  older_BF <- calculate_BF(older_data$rho[overlap_genes, , drop = FALSE], p_rhythmic)

  younger_rhythmic <- younger_BF > BF_threshold
  older_rhythmic <- older_BF > BF_threshold

  concordance <- ifelse(younger_rhythmic & older_rhythmic, "Concordant Rhythmic",
                        ifelse(!younger_rhythmic & !older_rhythmic, "Concordant Arrhythmic", 
                               "Discordant"))

  heatmap_matrix <- cbind(
    Younger = log(younger_BF),
    Older = log(older_BF)
  )
  rownames(heatmap_matrix) <- overlap_genes

  annotation_row <- data.frame(
    Concordance = factor(concordance, levels = c("Concordant Rhythmic", "Concordant Arrhythmic", "Discordant")),
    row.names = overlap_genes
  )

  annotation_colors <- list(
    Concordance = c("Concordant Rhythmic" = "#FFA500",
                    "Concordant Arrhythmic" = "white",
                    "Discordant" = "#4169E1")
  )

  # White at log BF <= 0 to red at the top of the scale
  bf_colors <- colorRampPalette(c("white", "#d73027"))(100)

  # Scale runs from 0 to the largest positive log BF, capped at 3
  positive_values <- heatmap_matrix[heatmap_matrix > 0]
  if(length(positive_values) > 0) {
    max_bf <- max(positive_values, na.rm = TRUE)
    max_bf <- min(max_bf, 3)
  } else {
    max_bf <- 1
  }

  # Negative log BF drawn as 0
  heatmap_matrix_adj <- heatmap_matrix
  heatmap_matrix_adj[heatmap_matrix_adj < 0] <- 0

  bf_breaks <- seq(0, max_bf, length.out = 101)

  p <- pheatmap(
    heatmap_matrix_adj,
    annotation_row = annotation_row,
    annotation_colors = annotation_colors,
    color = bf_colors,
    breaks = bf_breaks,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    main = paste(pathway_name),
    fontsize_row = 8,
    fontsize_col = 12,
    show_rownames = TRUE,
    show_colnames = TRUE,
    cellwidth = 50,
    cellheight = 15,
    width = 8,
    height = 40,
    filename = if(!is.null(save_path)) {
      if(!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)
      paste0(save_path, "/", gsub("[^A-Za-z0-9]", "_", pathway_name), "_heatmap.pdf")
    } else NA
  )
  
  return(p)
}

pathways_to_analyze <- c(
  "KEGG Adipocytokine signaling pathway",
  "KEGG Aldosterone-regulated sodium reabsorption",
  "KEGG Allograft rejection",
  "KEGG Amoebiasis",
  "KEGG Antigen processing and presentation",
  "KEGG Apoptosis",
  "KEGG Asthma",
  "KEGG Autoimmune thyroid disease",
  "KEGG Axon guidance",
  "KEGG Basal transcription factors",
  "KEGG Bile secretion",
  "KEGG Biosynthesis of unsaturated fatty acids",
  "KEGG Bladder cancer",
  "KEGG Cardiac muscle contraction",
  "KEGG Cell adhesion molecules (CAMs)",
  "KEGG Circadian rhythm",
  "KEGG Complement and coagulation cascades",
  "KEGG Cyanoamino acid metabolism",
  "KEGG Cytokine-cytokine receptor interaction",
  "KEGG Cytosolic DNA-sensing pathway",
  "KEGG Dilated cardiomyopathy",
  "KEGG Dorso-ventral axis formation",
  "KEGG ECM-receptor interaction",
  "KEGG Endocytosis",
  "KEGG Fatty acid biosynthesis",
  "KEGG Fatty acid degradation",
  "KEGG Fatty acid elongation",
  "KEGG Fc gamma R-mediated phagocytosis",
  "KEGG Focal adhesion",
  "KEGG Glioma",
  "KEGG Glycolysis / Gluconeogenesis",
  "KEGG Glycosaminoglycan biosynthesis - ganglio series",
  "KEGG Glycosaminoglycan biosynthesis - heparan sulfate / heparin",
  "KEGG Glycosaminoglycan biosynthesis - keratan sulfate",
  "KEGG Glycosphingolipid biosynthesis - ganglio series",
  "KEGG Hepatitis C",
  "KEGG Homologous recombination",
  "KEGG Hypertrophic cardiomyopathy (HCM)",
  "KEGG Intestinal immune network for IgA production",
  "KEGG Jak-STAT signaling pathway",
  "KEGG Leishmaniasis",
  "KEGG Long-term depression",
  "KEGG Long-term potentiation",
  "KEGG Lysosome",
  "KEGG Lysine biosynthesis",
  "KEGG Lysine degradation",
  "KEGG MAPK signaling pathway",
  "KEGG mRNA surveillance pathway",
  "KEGG mTOR signaling pathway",
  "KEGG Maturity onset diabetes of the young",
  "KEGG Mismatch repair",
  "KEGG Mucin type O-Glycan biosynthesis",
  "KEGG Natural killer cell mediated cytotoxicity",
  "KEGG Neuroactive ligand-receptor interaction",
  "KEGG Nitrogen metabolism",
  "KEGG Nucleotide excision repair",
  "KEGG Olfactory transduction",
  "KEGG Oocyte meiosis",
  "KEGG Osteoclast differentiation",
  "KEGG PPAR signaling pathway",
  "KEGG Pancreatic secretion",
  "KEGG Peroxisome",
  "KEGG Phagosome",
  "KEGG Phosphatidylinositol signaling system",
  "KEGG Prion diseases",
  "KEGG Proteasome",
  "KEGG Protein processing in endoplasmic reticulum",
  "KEGG Pyruvate metabolism",
  "KEGG Regulation of actin cytoskeleton",
  "KEGG Rheumatoid arthritis",
  "KEGG Sphingolipid metabolism",
  "KEGG Spliceosome",
  "KEGG Staphylococcus aureus infection",
  "KEGG Systemic lupus erythematosus",
  "KEGG Taurine and hypotaurine metabolism",
  "KEGG Type I diabetes mellitus",
  "KEGG Valine, leucine and isoleucine biosynthesis",
  "KEGG Viral myocarditis",
  "KEGG Vitamin digestion and absorption",
  "KEGG Wnt signaling pathway"
)

for(test_pathway in pathways_to_analyze) {

  cat("Processing:", test_pathway, "\n")

  if(!test_pathway %in% names(kegg.pathway.list_hsa)) {
    cat("  Pathway not found, skipping...\n")
    next
  }

  # Circadian rhythm: ARNTL renamed BMAL1, NR1D1 and NR1D2 added
  if(test_pathway == "KEGG Circadian rhythm") {
    pathway_genes <- gsub("ARNTL", "BMAL1", kegg.pathway.list_hsa[[test_pathway]])
    pathway_genes <- unique(c(pathway_genes, "NR1D1", "NR1D2"))
  } else {
    pathway_genes <- kegg.pathway.list_hsa[[test_pathway]]
  }

  overlap_check <- intersect(rownames(mcmc_age[[1]]$rho), pathway_genes)
  if(length(overlap_check) == 0) {
    cat("  No overlapping genes found, skipping...\n")
    next
  }
  
  cat("  Found", length(overlap_check), "overlapping genes\n")

  tryCatch({
    # Close open devices so the loop does not exhaust them
    while(dev.cur() > 1) dev.off()
    
    p <- plot_pathway_concordance(
      younger_data = mcmc_age[[1]],
      older_data = mcmc_age[[2]], 
      pathway_genes = pathway_genes,
      pathway_name = test_pathway,
      BF_threshold = 2,
      p_rhythmic = 0.20,
      save_path = file.path(current_aging, "results/brain_regions/plot/heatmap")
    )
    cat("  Successfully created plot for", test_pathway, "\n")
  }, error = function(e) {
    while(dev.cur() > 1) dev.off()
    cat("  Error creating plot for", test_pathway, ":", e$message, "\n")
  })
  
  cat("\n")
}