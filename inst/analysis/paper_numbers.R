## Numbers the paper quotes that no other analysis script prints, recomputed
## from the stored summaries, the case-study caches and the tables the other
## scripts wrote, and printed in the order of the paper's sections.
## README_numbers.md maps every quoted number to the line here, or to the
## script that prints it instead.
##
## Reads summary/hb and summary/hm, the caches under BAYRC_OUTPUT_DIR/figures,
## the tables under BAYRC_OUTPUT_DIR and BAYRC_FIGURE_DIR, the two atlas objects
## under BAYRC_GTEX_DIR/data, and the pooled and per-cycle mouse chains for
## Supplementary S9. Writes nothing.

suppressPackageStartupMessages(library(BayRC))

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
source(file.path(analysis.dir, "config.R"))

alpha <- 0.25
shift <- 2
## the application scripts leave one cache directory per comparison here
case.dir <- file.path(BAYRC_OUTPUT_DIR, "figures")

heading <- function(x) cat("\n==", x, "==\n")
line <- function(label, value) cat(sprintf("  %-50s %s\n", label, value))
pct <- function(k, n) sprintf("%d of %d (%.1f%%)", k, n, 100 * k / n)
wrap <- function(d) { d <- d %% 24; ifelse(d > 12, d - 24, d) }
circ_mean <- function(h) {
  a <- 2 * pi * h / 24
  (atan2(mean(sin(a)), mean(cos(a))) %% (2 * pi)) * 24 / (2 * pi)
}
n_called <- function(p) sum(bfdr_from_posterior(p, alpha = alpha)$rhythmic_genes)
## genes clearing more than one transition threshold, which take the last label
cleared_twice <- function(tr) {
  hit <- cbind(tr$p_gain >= tr$tau_gain, tr$p_loss >= tr$tau_loss,
               tr$p_cons >= tr$tau_cons)
  names(tr$gain_loss_status)[rowSums(hit) > 1]
}

kegg <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))
kegg <- lapply(kegg, function(g) replace(g, g == "ARNTL", "BMAL1"))
clock <- c("BMAL1", "CLOCK", "NPAS2", "PER1", "PER2", "PER3", "CRY1", "CRY2",
           "NR1D1", "NR1D2", "RORA", "RORB", "RORC", "DBP", "TEF", "HLF",
           "NFIL3", "BHLHE40", "BHLHE41", "CIART", "CSNK1D", "CSNK1E")

hb <- new.env(); load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"), envir = hb)
hm <- new.env(); load(file.path(BAYRC_RESULT_DIR, "summary", "hm", "mcmc_rho_BF3.RData"), envir = hm)
p_baboon <- sapply(hb$mcmc_data_baboon, rowMeans)
p_human  <- sapply(hb$mcmc_data_human, rowMeans)
genes_hb <- rownames(p_baboon)
genes_hm <- rownames(hm$mcmc_data_mouse[[1]])


heading("Applications: data and gene sets")
b <- new.env(); load(bayrc_file(BAYRC_GTEX_DIR, "data", "CAMO.bab.hum.RData"), envir = b)
zt <- as.numeric(sub("LUN[.]ZT", "", grep("^LUN[.]ZT", colnames(b$baboon_withTOD$baboon$LUN), value = TRUE)))
line("baboon tissues matched to GTEx", length(b$baboon_withTOD$baboon))
line("baboon lung samples, spacing", sprintf("%d at ZT %s, %g h apart", length(zt),
     paste(range(zt), collapse = "-"), unique(diff(zt))))
rm(b); invisible(gc())
m <- new.env(); load(bayrc_file(BAYRC_GTEX_DIR, "data", "CAMO.mouse.hum.RData"), envir = m)
ct <- as.numeric(sub(".*_CT", "", colnames(m$mice$count_clean$LUN)))
line("mouse tissues", paste(names(m$mice$count_clean), collapse = " "))
line("mouse samples, spacing, distinct phases", sprintf("%d at CT %s, %g h apart, %d distinct",
     length(ct), paste(range(ct), collapse = "-"), unique(diff(ct)), length(unique(ct %% 24))))
line("mouse-human matched genes", format(nrow(m$mice$count_clean$LUN), big.mark = ","))
rm(m); invisible(gc())
line("baboon-human matched genes", format(length(genes_hb), big.mark = ","))
line("three-species genes", format(length(intersect(genes_hb, genes_hm)), big.mark = ","))
only_hb <- setdiff(genes_hb, genes_hm)
line("baboon-human genes with no mouse counterpart", length(only_hb))
line("  of them clock genes / KEGG circadian genes", sprintf("%d / %d",
     length(intersect(only_hb, clock)), length(intersect(only_hb, kegg[["KEGG Circadian rhythm"]]))))
for (u in list(genes_hb, genes_hm))
  line(sprintf("KEGG pathways with >= 15 of %s genes", format(length(u), big.mark = ",")),
       sprintf("%d of %d; KEGG Circadian rhythm measured %d", sum(lengths(lapply(kegg, intersect, u)) >= 15),
               length(kegg), length(intersect(kegg[["KEGG Circadian rhythm"]], u))))


heading("Posterior against cosinor, baboon lung (Results; Supplementary S5)")
prior_odds <- 0.2 / 0.8
bf <- p_baboon[, "LUN"] / (1 - p_baboon[, "LUN"]) / prior_odds
for (k in c(1, 3, 5))
  line(sprintf("BF > %d: posterior cut, genes", k),
       sprintf("%.2f, %s", k * prior_odds / (1 + k * prior_odds), format(sum(bf > k), big.mark = ",")))
grid <- read.csv(file.path(BAYRC_FIGURE_DIR, "S5_threshold_grid.csv"), check.names = FALSE)
cat("  Table S5, from S5_threshold_grid.csv:\n")
print(grid, row.names = FALSE)


heading("Case Study 1: genome-wide and circadian concordance")
cg <- as.matrix(read.csv(file.path(BAYRC_OUTPUT_DIR, "all_plots", "Baboon_Concordance_Matrix.csv"),
                         row.names = 1, check.names = FALSE))
off <- cg[upper.tri(cg)]
line("tissues, pairs", sprintf("%d, %d", ncol(cg), length(off)))
line("median / max genome-wide score", sprintf("%.4f / %.4f", median(off), max(off)))
line("SCN-HIP score, rank among pairs", sprintf("%.4f, %d of %d", cg["SCN", "HIP"],
     sum(off >= cg["SCN", "HIP"]), length(off)))
diag(cg) <- 1
tree <- hclust(as.dist(1 - cg), method = "ward.D2")
## the smallest clade holding a tissue, read off the dendrogram Figure 2A draws
clade <- function(t) {
  for (k in seq(ncol(cg) - 1, 2)) {
    g <- cutree(tree, k = k)
    if (sum(g == g[t]) >= 4) return(paste(sort(names(g)[g == g[t]]), collapse = " "))
  }
}
line("Figure 2A clade holding PUT", clade("PUT"))
line("Figure 2A clade holding CER", clade("CER"))
pv <- read.csv(file.path(BAYRC_OUTPUT_DIR, "genome_concordance_pvalues.csv"))
for (i in seq_len(nrow(pv)))
  line(sprintf("%s score, p (%d permutations)", pv$pair[i], pv$permutations[i]),
       sprintf("%.4f, %.4f", pv$adjusted_concordance[i], pv$p_value[i]))
cc <- as.matrix(read.csv(file.path(BAYRC_OUTPUT_DIR, "heatmap_circadian_pairs", "within_baboon",
                                   "pairwise_concordance_baboon_circadian_Matrix.csv"),
                         row.names = 1, check.names = FALSE))
line("circadian matrix tissues, SCN included", sprintf("%d, %s", ncol(cc), "SCN" %in% colnames(cc)))


heading("Case Studies 1 and 2: transitions and phase in the three baboon pairs")
cache <- list(SCN_HIP = readRDS(file.path(case.dir, "baboon_SCN_HIP", "plot_data.rds")),
              PUT_SUN = readRDS(file.path(case.dir, "baboon_PUT_SUN", "plot_data.rds")),
              PUT_VIC = readRDS(file.path(case.dir, "baboon_PUT_VIC", "plot_data.rds")))
for (nm in names(cache)) {
  x  <- cache[[nm]]
  st <- x$trans$gain_loss_status
  mt <- x$maintained
  d  <- wrap(x$phase$peak2 - x$phase$peak1)[mt]
  cl <- x$phase_class[mt]
  sh <- cl == "Phase-shifted"
  cat(sprintf("\n  %s, over %s genes\n", nm, format(length(st), big.mark = ",")))
  line("gain / conserved / loss", sprintf("%d / %d / %d", sum(st == "Gain"), length(mt), sum(st == "Loss")))
  line("conserved within 2 h", pct(sum(abs(d) <= shift), length(mt)))
  for (k in c("Phase-conserved", "Phase-shifted", "Undetermined")) line(k, pct(sum(cl == k), length(mt)))
  line("shifted, comparator later / earlier", sprintf("%d / %d", sum(sh & d > 0), sum(sh & d < 0)))
  line("mean comparator minus reference (Figure 3E)", sprintf("%+.3f h", mean(d)))
  line("genes clearing two transition thresholds", length(cleared_twice(x$trans)))
}
cat("\n"); line("conserved PUT-VIC over PUT-SUN", sprintf("%.2f",
     length(cache$PUT_VIC$maintained) / length(cache$PUT_SUN$maintained)))


heading("Case Study 2: two-stage enrichment and Table 1")
pd <- "KEGG Parkinson disease"
for (cmp in c("SUN", "VIC")) {
  dir <- file.path(case.dir, paste0("baboon_PUT_", cmp))
  s1 <- read.csv(file.path(dir, "stage1_union.csv"))
  s2 <- read.csv(file.path(dir, "stage2_significant.csv"))
  pm <- read.csv(file.path(dir, "pathway_metrics.csv"), check.names = FALSE)
  names(pm) <- sub("^PUT_vs_[A-Z]+_", "", names(pm))
  cat(sprintf("\n  PUT against %s\n", cmp))
  line("pathways tested / Stage 1 retained (q < 0.20)", sprintf("%d / %d", nrow(s1), sum(s1$q < 0.20)))
  line("Stage 2 pathways / pathway-transition pairs", sprintf("%d / %d", length(unique(s2$pathway)), nrow(s2)))
  line("enriched pairs by transition", paste(names(table(s2$direction)), table(s2$direction), collapse = ", "))
  e <- unique(s2[, c("pathway", "size", "Expected_N_Gain", "Expected_N_Loss", "Expected_N_Conserved")])
  e$union <- e$Expected_N_Gain + e$Expected_N_Loss + e$Expected_N_Conserved
  share <- 100 * cbind(Gain = e$Expected_N_Gain, Loss = e$Expected_N_Loss,
                       Conserved = e$Expected_N_Conserved) / e$union
  line("largest slice in each enriched pathway", paste(names(table(colnames(share)[max.col(share)])),
       table(colnames(share)[max.col(share)]), collapse = ", "))
  line("loss share across enriched pathways", sprintf("%.1f%% to %.1f%%", min(share[, "Loss"]), max(share[, "Loss"])))
  ## Table 1: one row per enriched pathway, as the table prints it
  e$adj_c <- pm$AdjustedConcordance[match(e$pathway, pm$Pathway)]
  e$c_p   <- pm$PValue[match(e$pathway, pm$Pathway)]
  e$GLR   <- e$Expected_N_Gain / e$Expected_N_Loss
  e$enrichment <- sapply(e$pathway, function(p) {
    r <- s2[s2$pathway == p, ]
    paste(sprintf("%s q=%.2g", r$direction, r$q), collapse = "; ")
  })
  tab <- data.frame(pathway = sub("^KEGG ", "", e$pathway), size = e$size,
                    adj_c = round(e$adj_c, 3), GLR = round(e$GLR, 2), union = round(e$union, 1),
                    gain = round(e$Expected_N_Gain, 1), loss = round(e$Expected_N_Loss, 1),
                    cons = round(e$Expected_N_Conserved, 1), enrichment = e$enrichment)
  print(tab[order(-tab$size), ], row.names = FALSE)
  line("Parkinson disease pathway score, p", sprintf("%.4f, %.6f", e$adj_c[e$pathway == pd], e$c_p[e$pathway == pd]))
}
cat("\n"); line("conserved genes PUT-SUN, share of all", pct(length(cache$PUT_SUN$maintained), length(cache$PUT_SUN$measured)))


heading("Case Study 2: KEGG Parkinson disease gene by gene")
shifted <- list()
for (nm in c("PUT_SUN", "PUT_VIC")) {
  x  <- cache[[nm]]
  g  <- intersect(kegg[[pd]], x$measured)
  st <- x$trans$gain_loss_status[g]
  cl <- x$phase_class[g]
  d  <- wrap(x$phase$peak2 - x$phase$peak1)[g]
  lost <- g[st == "Loss"]
  cat(sprintf("\n  %s\n", nm))
  line("measured / lost / conserved / gained", sprintf("%d / %d / %d / %d", length(g),
       sum(st == "Loss"), sum(st == "Maintained"), sum(st == "Gain")))
  mt <- g[st == "Maintained"]
  for (k in c("Phase-shifted", "Phase-conserved", "Undetermined")) {
    s <- mt[cl[mt] == k]
    line(sprintf("%s (%d)", k, length(s)),
         paste(sprintf("%s %+.2f", s[order(d[s])], sort(d[s])), collapse = ", "))
  }
  s <- mt[cl[mt] == "Phase-shifted"]
  line("shifted range, all positive", sprintf("%+.2f to %+.2f h, %s", min(d[s]), max(d[s]), all(d[s] > 0)))
  shifted[[nm]] <- d[s]
  line("lost proteasome subunits", sprintf("%d: %s", sum(grepl("^PSM", lost)),
       paste(sort(lost[grepl("^PSM", lost)]), collapse = " ")))
  line("lost oxidative phosphorylation genes", sum(lost %in% kegg[["KEGG Oxidative phosphorylation"]]))
}
both <- intersect(names(shifted$PUT_SUN), names(shifted$PUT_VIC))
cat("\n"); line("shifted in both, SUN / VIC", paste(sprintf("%s %+.2f / %+.2f", both,
     shifted$PUT_SUN[both], shifted$PUT_VIC[both]), collapse = ", "))
v <- cache$PUT_VIC$trans$gain_loss_status
line("HSPA5, NDUFB2 against VIC", paste(v[c("HSPA5", "NDUFB2")], collapse = ", "))
vic <- unique(read.csv(file.path(case.dir, "baboon_PUT_VIC", "stage2_significant.csv"))$pathway)
ox  <- intersect(kegg[["KEGG Oxidative phosphorylation"]], genes_hb)
line("VIC-enriched sets sharing >= 40 OxPhos genes",
     sum(sapply(vic, function(p) length(intersect(intersect(kegg[[p]], genes_hb), ox)) >= 40)))


heading("Case Study 3: lung across three species, human as reference")
hbp <- new.env(); load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"), envir = hbp)
hmp <- new.env(); load(file.path(BAYRC_RESULT_DIR, "summary", "hm", "phi", "mcmc_phi_BF3.RData"), envir = hmp)
universe <- intersect(genes_hb, genes_hm)
lung <- function(rho, phi) list(rho = rho$LUN[universe, , drop = FALSE], phi = phi$LUN[universe, , drop = FALSE])
hum <- lung(hb$mcmc_data_human, hbp$mcmc_phi_human)
bab <- lung(hb$mcmc_data_baboon, hbp$mcmc_phi_baboon)
mou <- lung(hm$mcmc_data_mouse, hmp$mcmc_phi_mouse)
rm(hbp, hmp); invisible(gc())
## the reference goes first, as in plots/figure6/make_figure6.R
pair <- function(ref, cmp) {
  tr <- transition_classify(rowMeans(ref$rho), rowMeans(cmp$rho), bfdr_alpha = alpha)
  ph <- phase_infer(phi_matrix1 = ref$phi, phi_matrix2 = cmp$phi,
                    gain_loss_status = tr$gain_loss_status, bfdr_alpha = alpha,
                    shift = shift, P = 24, compute_hdi = TRUE)
  list(trans = tr, phase = ph, d = wrap(ph$peak2 - ph$peak1))
}
res <- list(baboon = pair(hum, bab), mouse = pair(hum, mou))
ph <- rowMeans(hum$rho)
line("genes", format(length(universe), big.mark = ","))
line("rhythmic human / baboon / mouse", sprintf("%d (%.1f%%) / %s / %d", n_called(ph),
     100 * n_called(ph) / length(ph), format(n_called(rowMeans(bab$rho)), big.mark = ","),
     n_called(rowMeans(mou$rho))))
circ <- intersect(kegg[["KEGG Circadian rhythm"]], universe)
line("KEGG Circadian rhythm genes", length(circ))
for (sp in names(res)) {
  r  <- res[[sp]]
  mt <- names(r$trans$gain_loss_status)[r$trans$gain_loss_status == "Maintained"]
  cl <- ifelse(r$phase$flag_cons[mt], "conserved", ifelse(r$phase$flag_shift[mt], "shifted", "undetermined"))
  cat(sprintf("\n  %s against human\n", sp))
  line("tau_cons, rhythm-conserved", sprintf("%.4f, %d", r$trans$tau_cons, length(mt)))
  ## p_cons is p_human x p_comparator, so the conserved set cannot be larger
  ## than the number of human genes that reach tau_cons on their own
  line("human genes reaching tau_cons", sum(ph >= r$trans$tau_cons))
  line("within 2 h of human", pct(sum(abs(r$d[mt]) <= shift), length(mt)))
  line("phase conserved / shifted / undetermined", paste(table(factor(cl,
       levels = c("conserved", "shifted", "undetermined"))), collapse = " / "))
  two <- cleared_twice(r$trans)
  line("genes clearing two thresholds", if (length(two)) paste(two, collapse = " ") else "none")
  for (g in two) line(sprintf("  %s p_gain, tau_gain, p_cons, tau_cons", g),
       sprintf("%.4f %.4f %.4f %.4f", r$trans$p_gain[match(g, names(r$trans$gain_loss_status))],
               r$trans$tau_gain, r$trans$p_cons[match(g, names(r$trans$gain_loss_status))], r$trans$tau_cons))
  pg <- intersect(circ, mt)
  for (k in c("conserved", "shifted", "undetermined")) {
    s <- pg[cl[match(pg, mt)] == k]
    line(sprintf("pathway genes phase-%s", k), paste(sprintf("%s %+.2f", s, r$d[s]), collapse = ", "))
  }
}
held <- setdiff(circ, unlist(lapply(res, function(r)
  names(r$trans$gain_loss_status)[r$trans$gain_loss_status == "Maintained"])))
cat("\n"); line("pathway genes conserved in neither comparator", length(held))
line("human, baboon posterior PER1 CUL1 CSNK1D", paste(sprintf("%.3f/%.3f",
     ph[c("PER1", "CUL1", "CSNK1D")], rowMeans(bab$rho)[c("PER1", "CUL1", "CSNK1D")]), collapse = " "))
line("highest human posterior among those", sprintf("%s %.3f", names(which.max(ph[held])), max(ph[held])))

six <- read.csv(file.path(BAYRC_OUTPUT_DIR, "cross_species_mouse", "cross_species_six_tissues.csv"))
cat("\n  six tissues, rhythm-conserved and phase-conserved (cross_species_six_tissues.csv)\n")
print(six[, c("comparison", "tissue", "maintained", "within_pct", "phase_conserved", "phase_shifted")],
      row.names = FALSE)

cat("\n  Stage 1 of the twelve comparisons, KEGG Circadian rhythm\n")
active <- list()
for (sp in c("baboon", "mouse")) for (t in c("AOR", "CER", "HEA", "KIC", "LIV", "LUN")) {
  dir <- file.path(case.dir, sprintf("%s_human_%s_genome", sp, t))
  if (sp == "baboon" && t == "LUN") dir <- file.path(case.dir, "baboon_human_LUN")
  s <- read.csv(file.path(dir, "stage1_union.csv"))
  s <- s[order(s$pval), ]
  i <- match("KEGG Circadian rhythm", s$pathway)
  active[[sp]][[t]] <- s$pathway[s$q < 0.20]
  line(sprintf("%s %s: genes, tested, retained, rank, q", sp, t),
       sprintf("%d, %d, %d, %d, %.2g", length(readRDS(file.path(dir, "plot_data.rds"))$measured),
               nrow(s), sum(s$q < 0.20), i, s$q[i]))
}
for (t in names(active$mouse))
  line(sprintf("retained by both comparators, %s", t),
       paste(sub("^KEGG ", "", intersect(active$baboon[[t]], active$mouse[[t]])), collapse = "; "))


heading("Discussion: human rhythmic calls across the 26 tissues")
hc <- sort(apply(p_human, 2, n_called), decreasing = TRUE)
line("fewest, most", sprintf("%s %d, %s %d", names(hc)[length(hc)], min(hc), names(hc)[1], max(hc)))
line("lung calls, rank of 26", sprintf("%d, %d", hc[["LUN"]], match("LUN", names(hc))))


heading("Supplementary S4: residual tests over every gene")
for (t in c("LUN", "PUT", "SUN")) {
  r <- read.csv(file.path(BAYRC_FIGURE_DIR, sprintf("Cosinor_residual_diagnostics_%s.csv", t)))
  line(sprintf("%s: genes, Shapiro p >= 0.05, variance p >= 0.05", t),
       sprintf("%s, %.1f%%, %.1f%%", format(nrow(r), big.mark = ","),
               100 * mean(r$shapiro.pval >= 0.05, na.rm = TRUE),
               100 * mean(r$hetero.lm.pval >= 0.05, na.rm = TRUE)))
}


heading("Supplementary S8: calibration and operating characteristics")
cal <- file.path(BAYRC_OUTPUT_DIR, "calibration")
rec <- trimws(grep("rhythmic_prop|kept_draws", readLines(file.path(cal, "run_record.txt")), value = TRUE))
line("rhythmic share realised, kept draws", paste(rec, collapse = "; "))
tb <- read.csv(file.path(cal, "bfdr_calibration_summary.csv"))
line("realised over nominal, strictest to loosest", sprintf("%.2f to %.2f", tb$ratio[1], tb$ratio[nrow(tb)]))
oc  <- read.csv(file.path(cal, "bfdr_operating_characteristics.csv"))
auc <- read.csv(file.path(cal, "bayrc_cosinor_auc.csv"))
line("largest |AUC BayRC - AUC cosinor|", sprintf("%.4f", max(abs(auc$BayRC - auc$Cosinor))))
line("smallest power, n = 48 and A/sigma >= 1", sprintf("%.4f", min(oc$power[oc$n == 48 & oc$A >= 1])))
line("largest type I error minus alpha", sprintf("%+.4f", max(oc$type1 - oc$alpha)))
line("largest FDR minus alpha, n = 48, A/sigma <= 1", sprintf("%+.4f", max((oc$fdr - oc$alpha)[oc$n == 48 & oc$A <= 1])))
for (a in c(0.05, 0.10, 0.25)) {
  f <- oc$fdr[oc$n == 48 & oc$A >= 1.5 & oc$alpha == a]
  line(sprintf("FDR at n = 48, A/sigma >= 1.5, alpha %.2f", a), sprintf("%.3f to %.3f", min(f), max(f)))
}


heading("Supplementary S9: the two mouse cycles")
## the chains are indexed by mouse identifier; the summary maps them to symbols
obj <- hm$mcmc_data_mouse[[1]]
sym <- setNames(rownames(obj), attr(obj, "ensembl_gene_ids"))
m <- new.env(); load(bayrc_file(BAYRC_GTEX_DIR, "data", "CAMO.mouse.hum.RData"), envir = m)
core <- c("BMAL1", "PER1", "PER2", "NR1D1", "NR1D2", "DBP", "CRY1", "CRY2")
arms <- c(pooled = "mice", cycle1 = "mice_cycle/cycle1", cycle2 = "mice_cycle/cycle2")
gap <- c(); pooled_best <- 0
for (t in names(m$mice$count_clean)) {
  i <- match(core, unname(sym[rownames(m$mice$count_clean[[t]])])); i <- i[!is.na(i)]
  arm <- lapply(arms, function(a) {
    x <- readRDS(list.files(file.path(BAYRC_RESULT_DIR, a, t), "\\.RDS$", full.names = TRUE)[1])
    list(post = mean(rowMeans(x$rho)[i]), peak = apply(x$phi[i, , drop = FALSE], 1, circ_mean))
  })
  gap[t] <- median(abs(wrap(arm$cycle1$peak - arm$cycle2$peak)))
  best <- arm$pooled$post > max(arm$cycle1$post, arm$cycle2$post)
  pooled_best <- pooled_best + best
  line(sprintf("%s: clock posterior pooled / c1 / c2, pooled highest", t),
       sprintf("%.3f / %.3f / %.3f, %s", arm$pooled$post, arm$cycle1$post, arm$cycle2$post, best))
}
line("tissues where pooled clock posterior is highest", sprintf("%d of %d", pooled_best, length(gap)))
line("median |cycle 1 - cycle 2| clock peak, range", sprintf("%.2f to %.2f h", min(gap), max(gap)))
cat("\n")
