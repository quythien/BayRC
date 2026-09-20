## Which conservation-enriched SUN-PUT pathways could carry a heatmap panel.
##
## Enrichment q-values come from the 229-pathway release, the gene sets from
## the 354-pathway release, so the two are joined on a normalised name rather
## than on string equality. An unmapped name stops the script.
##
## Usage: Rscript sunput_panel_candidates.R [alpha1] [alpha2] [shift]

suppressPackageStartupMessages(library(BayRC))

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

args  <- commandArgs(trailingOnly = TRUE)
a1    <- if (length(args) >= 1) as.numeric(args[1]) else 0.25
a2    <- if (length(args) >= 2) as.numeric(args[2]) else 0.30
shift <- if (length(args) >= 3) as.numeric(args[3]) else 2

out.dir <- file.path(BAYRC_OUTPUT_DIR, "param_search")
cache <- readRDS(file.path(out.dir, "cache_PUT_SUN.rds"))
datA <- cache$A; datB <- cache$B           # A = PUT, B = SUN
genes <- rownames(datA$rho)

kegg.env <- new.env()
load(file.path(BAYRC_PATHWAY_DIR, "kegg.pathway.list_hsa.RData"), envir = kegg.env)
kegg229 <- kegg.env$kegg.pathway.list_hsa
kegg354 <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))

## join the two releases on a punctuation- and case-insensitive name
norm_name <- function(x) {
  x <- sub("^KEGG ", "", x)
  x <- tolower(x)
  x <- gsub("['’]s\\b", "", x)      # Parkinson's -> Parkinson
  x <- gsub("['’]", "", x)
  x <- gsub("[^a-z0-9]+", " ", x)
  trimws(gsub(" +", " ", x))
}
map354 <- setNames(names(kegg354), norm_name(names(kegg354)))
to354 <- function(nm229) {
  hit <- map354[norm_name(nm229)]
  unname(hit)
}

## stage 1 and stage 2 on the 229 release, size filter on measured genes
raw <- readRDS(file.path(out.dir, "pathselect_raw_nperm10000.rds"))
u <- raw$kegg229_union
u <- u[u$measured >= 15, ]
active <- u$pathway[u$pval < 0.05]
cons <- raw$kegg229_conserved
cons <- cons[cons$pathway %in% active, ]
cons$q <- p.adjust(cons$pval, "BH")
cons <- cons[order(cons$pval), ]
cand <- cons[cons$q < 0.20, ]
cat("stage 1:", length(active), "active of", nrow(u), "tested (229 release)\n")
cat("conservation-enriched at q < 0.20:", nrow(cand), "\n\n")

unmapped <- cand$pathway[is.na(to354(cand$pathway))]
if (length(unmapped)) {
  for (nm in unmapped) {
    near <- names(kegg354)[agrepl(sub("^KEGG ", "", nm), names(kegg354),
                                  max.distance = 0.25, ignore.case = TRUE)]
    cat("unmapped:", nm, "| closest in 354:",
        if (length(near)) paste(head(near, 3), collapse = " / ") else "none", "\n")
  }
  stop(length(unmapped), " pathway name(s) have no 354-release match")
}

## signed phase difference, SUN minus PUT
dphi <- ((datB$phi - datA$phi + 12) %% 24) - 12
med  <- apply(dphi, 1, median)
rm(dphi)

pA <- rowMeans(datA$rho); pB <- rowMeans(datB$rho)
fits <- lapply(c(a1, a2), function(a) {
  invisible(capture.output(tr <- transition_classify(pA, pB, bfdr_alpha = a)))
  ph <- phase_infer(phi_matrix1 = datA$phi, phi_matrix2 = datB$phi,
                    gain_loss_status = tr$gain_loss_status, bfdr_alpha = a,
                    shift = shift, P = 24, compute_hdi = FALSE)
  list(alpha = a, status = tr$gain_loss_status, phase = ph)
})

rows <- list()
for (i in seq_len(nrow(cand))) {
  nm229 <- cand$pathway[i]; nm354 <- to354(nm229)
  gset <- intersect(kegg354[[nm354]], genes)
  for (f in fits) {
    st <- f$status[gset]
    maint <- names(st)[st == "Maintained"]
    sh <- intersect(maint, names(f$phase$flag_shift)[f$phase$flag_shift])
    co <- intersect(maint, names(f$phase$flag_cons)[f$phase$flag_cons])
    d <- med[sh]
    rows[[length(rows) + 1]] <- data.frame(
      pathway = nm354, q_cons = signif(cand$q[i], 3), alpha = f$alpha,
      measured = length(gset), maintained = length(maint),
      shifted = length(sh), conserved = length(co),
      undetermined = length(maint) - length(sh) - length(co),
      mean_d = if (length(d)) round(mean(d), 2) else NA,
      sd_d = if (length(d) > 1) round(sd(d), 2) else NA,
      min_d = if (length(d)) round(min(d), 1) else NA,
      max_d = if (length(d)) round(max(d), 1) else NA,
      pct_dominant = if (length(d))
        round(100 * max(mean(d > 0), mean(d < 0))) else NA)
  }
}
d <- do.call(rbind, rows)

## a panel needs enough maintained genes to fill it and one clear direction
score <- with(d[d$alpha == a2, ],
              ifelse(maintained >= 8 & pct_dominant >= 80, maintained, 0))
ord <- d[d$alpha == a2, ][order(-score, -d$maintained[d$alpha == a2]), "pathway"]
d <- d[order(match(d$pathway, ord), d$alpha), ]

write.csv(d, file.path(out.dir, "sunput_panel_candidates.csv"), row.names = FALSE)
cat("signed phase difference is SUN minus PUT; window =", shift, "h\n\n")
print(d, row.names = FALSE)
