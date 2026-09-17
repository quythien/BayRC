################################################################################
# How far do the within-baboon tissue pairs move between the two samplers?
#
# The genome-wide concordance heatmap (Figure 2) is built from the adjusted
# concordance of all choose(26, 2) = 325 within-baboon pairs. Under the fixed
# sampler the most-concordant pairs reorder a good deal more than they do in
# human, so this writes out the reordering in full: the top 20 under each arm
# with its rank under the other, the rank-shift distribution over all 325
# pairs, and every pair that crosses the top-10 boundary in either direction.
#
# Tissue codes are expanded from the mapping in
# GTEXdata/R/25_tissues_preprocessing.R so the pairs read without a lookup.
#
# Usage:
#   Rscript baboon_pair_rank_shift.R [baseline.csv] [fixed.rds] [out.csv]
################################################################################

# Paths come from inst/analysis/config.R; override with the matching env var.
this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(if (is.na(this.file)) dirname(getwd()) else
                   dirname(dirname(normalizePath(this.file))), "config.R"))

args <- commandArgs(trailingOnly = TRUE)

# The 2025 arm's table is the one plots/heatmap_baboon.R wrote; it is still on
# disk under the aging results tree, so nothing has to be recomputed.
baseline.file <- if (length(args) >= 1) args[1] else
  file.path(BAYRC_AGING_DIR, "results", "baboon", "output_final",
            "concordance_within_baboon.csv")
fixed.file <- if (length(args) >= 2) args[2] else
  file.path(BAYRC_OUTPUT_DIR, "within_baboon_pairwise_concordance.rds")
out.file <- if (length(args) >= 3) args[3] else
  file.path(BAYRC_OUTPUT_DIR, "baboon_pair_rank_shift.csv")

# Tissue names --------------------------------------------------------------
# Same codes 25_tissues_preprocessing.R assigns, as the GTEx tissue they were
# matched to. The baboon regions carry the same labels.

tissue.name <- c(
  AMY = "Brain - Amygdala",
  AOR = "Artery - Aorta",
  ASC = "Colon - Transverse",
  CER = "Brain - Cerebellum",
  HEA = "Heart - Atrial Appendage",
  HIP = "Brain - Hippocampus",
  ILE = "Small Intestine - Terminal Ileum",
  KIC = "Kidney - Cortex",
  LIV = "Liver",
  LUN = "Lung",
  MUA = "Muscle - Skeletal",
  OES = "Esophagus - Muscularis",
  OMF = "Adipose - Visceral (Omentum)",
  PAN = "Pancreas",
  PIT = "Pituitary",
  PRC = "Brain - Frontal Cortex (BA9)",
  PUT = "Brain - Putamen (basal ganglia)",
  SCN = "Brain - Hypothalamus",
  SKI = "Skin - Not Sun Exposed (Suprapubic)",
  SPL = "Spleen",
  STF = "Stomach",
  SUN = "Brain - Substantia nigra",
  TES = "Testis",
  THR = "Thyroid",
  VIC = "Brain - Cortex",
  WAS = "Adipose - Subcutaneous")

expand <- function(code) ifelse(code %in% names(tissue.name),
                                tissue.name[code], code)

# Load and align ------------------------------------------------------------

base <- utils::read.csv(baseline.file, stringsAsFactors = FALSE)
fix  <- readRDS(fixed.file)

# a pair is unordered, so key on the sorted tissue codes
pair.key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "_")
base$key <- pair.key(base$Tissue1, base$Tissue2)
fix$key  <- pair.key(fix$tissue1,  fix$tissue2)

if (anyDuplicated(base$key) || anyDuplicated(fix$key))
  stop("a tissue pair appears twice")

d <- merge(base[, c("key", "Tissue1", "Tissue2", "Jaccard_Obs", "Jaccard_Adj")],
           fix[,  c("key", "raw", "adjusted")], by = "key")
if (nrow(d) != nrow(base) || nrow(d) != nrow(fix))
  stop("the two arms do not cover the same pairs: ",
       nrow(base), " baseline, ", nrow(fix), " fixed, ", nrow(d), " matched")

d$rank_base  <- rank(-d$Jaccard_Adj, ties.method = "min")
d$rank_fixed <- rank(-d$adjusted,    ties.method = "min")
d$rank_shift <- d$rank_fixed - d$rank_base
d$pair <- paste(d$Tissue1, d$Tissue2, sep = "-")
d$pair_long <- paste(expand(d$Tissue1), expand(d$Tissue2), sep = "  vs  ")

d <- d[, c("pair", "pair_long", "Tissue1", "Tissue2",
           "Jaccard_Obs", "Jaccard_Adj", "raw", "adjusted",
           "rank_base", "rank_fixed", "rank_shift")]
names(d)[names(d) == "Jaccard_Obs"] <- "raw_base"
names(d)[names(d) == "Jaccard_Adj"] <- "adj_base"
names(d)[names(d) == "raw"]         <- "raw_fixed"
names(d)[names(d) == "adjusted"]    <- "adj_fixed"

utils::write.csv(d[order(d$rank_base), ], out.file, row.names = FALSE)

# Report --------------------------------------------------------------------

options(width = 200)
cat("within-baboon adjusted concordance, 2025 baseline vs fixed sampler\n")
cat(nrow(d), "pairs\n\n")

cat("agreement over all", nrow(d), "pairs\n")
cat(sprintf("  Pearson  %.3f\n", stats::cor(d$adj_base, d$adj_fixed)))
cat(sprintf("  Spearman %.3f\n", stats::cor(d$adj_base, d$adj_fixed,
                                            method = "spearman")))
cat(sprintf("  median   %.4f -> %.4f\n", stats::median(d$adj_base),
            stats::median(d$adj_fixed)))
cat(sprintf("  range    %.4f to %.4f -> %.4f to %.4f\n",
            min(d$adj_base), max(d$adj_base),
            min(d$adj_fixed), max(d$adj_fixed)))

cat("\nrank shift (fixed rank minus baseline rank)\n")
cat(sprintf("  median |shift| %.0f    mean |shift| %.1f    max |shift| %d\n",
            stats::median(abs(d$rank_shift)), mean(abs(d$rank_shift)),
            max(abs(d$rank_shift))))
cat(sprintf("  within 10 ranks: %d of %d (%.0f%%)\n",
            sum(abs(d$rank_shift) <= 10), nrow(d),
            100 * mean(abs(d$rank_shift) <= 10)))
cat(sprintf("  moved more than 50 ranks: %d\n", sum(abs(d$rank_shift) > 50)))

show <- function(x, label) {
  cat("\n", label, "\n", sep = "")
  print(data.frame(
    pair  = x$pair,
    tissues = x$pair_long,
    adj_base  = round(x$adj_base, 4),
    adj_fixed = round(x$adj_fixed, 4),
    rank_base  = x$rank_base,
    rank_fixed = x$rank_fixed,
    shift = sprintf("%+d", x$rank_shift),
    stringsAsFactors = FALSE), row.names = FALSE)
}

show(head(d[order(d$rank_base), ], 20),
     "top 20 in the 2025 baseline, with where they land in the fixed arm")
show(head(d[order(d$rank_fixed), ], 20),
     "top 20 in the fixed arm, with where they came from in the baseline")

top10.base  <- d$pair[d$rank_base  <= 10]
top10.fixed <- d$pair[d$rank_fixed <= 10]

cat("\ntop-10 membership: ", length(intersect(top10.base, top10.fixed)),
    " of 10 shared\n", sep = "")

dropped <- d[d$pair %in% setdiff(top10.base, top10.fixed), ]
entered <- d[d$pair %in% setdiff(top10.fixed, top10.base), ]
if (nrow(dropped)) show(dropped[order(dropped$rank_base), ],
                        "LEAVES the top 10 under the fixed sampler")
if (nrow(entered)) show(entered[order(entered$rank_fixed), ],
                        "ENTERS the top 10 under the fixed sampler")

# Why the pairs move ---------------------------------------------------------
# The obvious suspicion is that a pair moves because the fix changed one of its
# two tissues more than the others. It does not: the shift tracks how many
# rhythmic genes a tissue has, not how much its count changed. Sparse tissues
# rise and dense ones fall, because the adjusted concordance keeps a residual
# dependence on the size of the two gene sets and that dependence changes sign
# between the arms.

counts.file <- Sys.getenv("BAYRC_PAPER_TABLE",
  unset = file.path(dirname(BAYRC_OUTPUT_DIR), "paper_table.rds"))
if (file.exists(counts.file)) {
  tab <- readRDS(counts.file)$tissue
  bab <- tab[tab$species == "baboon", ]
  cnt.base  <- stats::setNames(bab$bf3_base,  bab$tissue)
  cnt.fixed <- stats::setNames(bab$bf3_fixed, bab$tissue)

  mean.shift <- vapply(names(cnt.fixed), function(x)
    mean(d$rank_shift[d$Tissue1 == x | d$Tissue2 == x]), numeric(1))

  cat("\nwhat predicts a tissue's mean rank shift\n")
  cat(sprintf("  its rhythmic count, fixed arm     Pearson %+.3f  Spearman %+.3f\n",
              stats::cor(cnt.fixed, mean.shift),
              stats::cor(cnt.fixed, mean.shift, method = "spearman")))
  cat(sprintf("  how much the fix changed it       Pearson %+.3f  Spearman %+.3f\n",
              stats::cor(cnt.fixed/cnt.base, mean.shift),
              stats::cor(cnt.fixed/cnt.base, mean.shift, method = "spearman")))

  gm.base  <- sqrt(cnt.base[d$Tissue1]  * cnt.base[d$Tissue2])
  gm.fixed <- sqrt(cnt.fixed[d$Tissue1] * cnt.fixed[d$Tissue2])
  cat("\nadjusted concordance against the geometric mean of the pair's two counts\n")
  cat(sprintf("  baseline arm  Pearson %+.3f\n",
              stats::cor(gm.base,  d$adj_base)))
  cat(sprintf("  fixed arm     Pearson %+.3f\n",
              stats::cor(gm.fixed, d$adj_fixed)))
}

cat("\nwrote", out.file, "\n")
