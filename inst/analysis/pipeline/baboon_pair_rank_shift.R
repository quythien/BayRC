################################################################################
# How far do the within-baboon tissue pairs move between two runs?
#
# The genome-wide concordance heatmap (Figure 2) is built from the adjusted
# concordance of all choose(26, 2) = 325 within-baboon pairs. The most
# concordant pairs reorder more between runs in baboon than they do in human,
# so this writes out the reordering in full: the top 20 under each run with its
# rank under the other, the rank-shift distribution over all 325 pairs, and
# every pair that crosses the top-10 boundary in either direction.
#
# Tissue codes are expanded from the mapping in
# GTEXdata/R/25_tissues_preprocessing.R so the pairs read without a lookup.
#
# Usage:
#   Rscript baboon_pair_rank_shift.R [previous.csv] [current.rds] [out.csv]
################################################################################

# Paths come from inst/analysis/config.R; override with the matching env var.
this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
source(file.path(if (is.na(this.file)) dirname(getwd()) else
                   dirname(dirname(normalizePath(this.file))), "config.R"))

args <- commandArgs(trailingOnly = TRUE)

# previous run: the 2025 table plots/heatmap_baboon.R wrote under the aging tree
prev.file <- if (length(args) >= 1) args[1] else
  file.path(BAYRC_AGING_DIR, "results", "baboon", "output_final",
            "concordance_within_baboon.csv")
curr.file <- if (length(args) >= 2) args[2] else
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

prev <- utils::read.csv(prev.file, stringsAsFactors = FALSE)
curr <- readRDS(curr.file)

# a pair is unordered, so key on the sorted tissue codes
pair.key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "_")
prev$key <- pair.key(prev$Tissue1, prev$Tissue2)
curr$key <- pair.key(curr$tissue1, curr$tissue2)

if (anyDuplicated(prev$key) || anyDuplicated(curr$key))
  stop("a tissue pair appears twice")

d <- merge(prev[, c("key", "Tissue1", "Tissue2", "Jaccard_Obs", "Jaccard_Adj")],
           curr[, c("key", "raw", "adjusted")], by = "key")
if (nrow(d) != nrow(prev) || nrow(d) != nrow(curr))
  stop("the two runs do not cover the same pairs: ",
       nrow(prev), " previous, ", nrow(curr), " current, ", nrow(d), " matched")

d$rank_prev <- rank(-d$Jaccard_Adj, ties.method = "min")
d$rank_curr <- rank(-d$adjusted,    ties.method = "min")
d$rank_shift <- d$rank_curr - d$rank_prev
d$pair <- paste(d$Tissue1, d$Tissue2, sep = "-")
d$pair_long <- paste(expand(d$Tissue1), expand(d$Tissue2), sep = "  vs  ")

d <- d[, c("pair", "pair_long", "Tissue1", "Tissue2",
           "Jaccard_Obs", "Jaccard_Adj", "raw", "adjusted",
           "rank_prev", "rank_curr", "rank_shift")]
names(d)[names(d) == "Jaccard_Obs"] <- "raw_prev"
names(d)[names(d) == "Jaccard_Adj"] <- "adj_prev"
names(d)[names(d) == "raw"]         <- "raw_curr"
names(d)[names(d) == "adjusted"]    <- "adj_curr"

utils::write.csv(d[order(d$rank_prev), ], out.file, row.names = FALSE)

# Report --------------------------------------------------------------------

options(width = 200)
cat("within-baboon adjusted concordance, 2025 run vs current run\n")
cat(nrow(d), "pairs\n\n")

cat("agreement over all", nrow(d), "pairs\n")
cat(sprintf("  Pearson  %.3f\n", stats::cor(d$adj_prev, d$adj_curr)))
cat(sprintf("  Spearman %.3f\n", stats::cor(d$adj_prev, d$adj_curr,
                                            method = "spearman")))
cat(sprintf("  median   %.4f -> %.4f\n", stats::median(d$adj_prev),
            stats::median(d$adj_curr)))
cat(sprintf("  range    %.4f to %.4f -> %.4f to %.4f\n",
            min(d$adj_prev), max(d$adj_prev),
            min(d$adj_curr), max(d$adj_curr)))

cat("\nrank shift (current rank minus previous rank)\n")
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
    adj_prev = round(x$adj_prev, 4),
    adj_curr = round(x$adj_curr, 4),
    rank_prev = x$rank_prev,
    rank_curr = x$rank_curr,
    shift = sprintf("%+d", x$rank_shift),
    stringsAsFactors = FALSE), row.names = FALSE)
}

show(head(d[order(d$rank_prev), ], 20),
     "top 20 in the 2025 run, with where they land in the current run")
show(head(d[order(d$rank_curr), ], 20),
     "top 20 in the current run, with where they came from")

top10.prev <- d$pair[d$rank_prev <= 10]
top10.curr <- d$pair[d$rank_curr <= 10]

cat("\ntop-10 membership: ", length(intersect(top10.prev, top10.curr)),
    " of 10 shared\n", sep = "")

dropped <- d[d$pair %in% setdiff(top10.prev, top10.curr), ]
entered <- d[d$pair %in% setdiff(top10.curr, top10.prev), ]
if (nrow(dropped)) show(dropped[order(dropped$rank_prev), ],
                        "leaves the top 10 in the current run")
if (nrow(entered)) show(entered[order(entered$rank_curr), ],
                        "enters the top 10 in the current run")

# What moves a pair ----------------------------------------------------------
# Rank shift against each tissue's rhythmic count and its ratio between runs.

counts.file <- Sys.getenv("BAYRC_PAPER_TABLE",
  unset = file.path(dirname(BAYRC_OUTPUT_DIR), "paper_table.rds"))
if (file.exists(counts.file)) {
  tab <- readRDS(counts.file)$tissue
  bab <- tab[tab$species == "baboon", ]
  cnt.prev <- stats::setNames(bab$bf3_base,  bab$tissue)
  cnt.curr <- stats::setNames(bab$bf3_fixed, bab$tissue)

  mean.shift <- vapply(names(cnt.curr), function(x)
    mean(d$rank_shift[d$Tissue1 == x | d$Tissue2 == x]), numeric(1))

  cat("\nwhat predicts a tissue's mean rank shift\n")
  cat(sprintf("  its rhythmic count, current run   Pearson %+.3f  Spearman %+.3f\n",
              stats::cor(cnt.curr, mean.shift),
              stats::cor(cnt.curr, mean.shift, method = "spearman")))
  cat(sprintf("  the ratio between the two runs    Pearson %+.3f  Spearman %+.3f\n",
              stats::cor(cnt.curr/cnt.prev, mean.shift),
              stats::cor(cnt.curr/cnt.prev, mean.shift, method = "spearman")))

  gm.prev <- sqrt(cnt.prev[d$Tissue1] * cnt.prev[d$Tissue2])
  gm.curr <- sqrt(cnt.curr[d$Tissue1] * cnt.curr[d$Tissue2])
  cat("\nadjusted concordance against the geometric mean of the pair's two counts\n")
  cat(sprintf("  2025 run     Pearson %+.3f\n",
              stats::cor(gm.prev, d$adj_prev)))
  cat(sprintf("  current run  Pearson %+.3f\n",
              stats::cor(gm.curr, d$adj_curr)))
}

cat("\nwrote", out.file, "\n")
