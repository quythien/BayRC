################################################################################
# Build the per-species rho / phi summary objects the downstream scripts load.
#
# Every figure and concordance script starts with
#   load(".../summary/hb/mcmc_rho_BF3.RData")
#   load(".../summary/hb/phi/mcmc_phi_BF3.RData")
# but nothing in the repository ever wrote those two files. This script does.
#
# It walks the per-tissue MCMC output written by CAMO_h_b.R (or by the
# fixed-arm runner), keeps the rho and phi chains, puts the rows in Ensembl-ID
# order, renames them to HGNC symbols, and attaches the three attributes the
# downstream code reads off these matrices: "symbols", "RHYindex" (Bayes
# factor > BF at the p_rhythmic prior) and "ensembl_gene_ids". That is the
# layout of the 2025 artifacts, and --validate checks the output against them.
#
# Usage:
#   Rscript summarize_rho_phi.R [result.dir] [summary.dir] [--validate]
#
# result.dir   directory holding human/<TIS>/gtex_<TIS>_bay_1.RDS and
#              baboon/<TIS>/baboon_<TIS>_bay_1.RDS  (default: the 2025 run)
# summary.dir  where to write mcmc_rho_BF3.RData and phi/mcmc_phi_BF3.RData
# --validate   compare the rebuilt objects against the ones already in
#              summary.dir instead of overwriting them
#
# Both directories may also be set through BAYRC_RESULT_DIR and
# BAYRC_SUMMARY_DIR.
################################################################################

# config.R lives one directory up, in inst/analysis/
this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else
  dirname(dirname(normalizePath(this.file)))
source(file.path(analysis.dir, "config.R"))
suppressPackageStartupMessages(library(BayRC))   # summarize_bay()

args     <- commandArgs(trailingOnly = TRUE)
validate <- "--validate" %in% args
args     <- setdiff(args, "--validate")

result.dir <- if (length(args) >= 1) args[1] else
  Sys.getenv("BAYRC_RESULT_DIR", unset = file.path(BAYRC_GTEX_DIR, "result"))
summary.dir <- if (length(args) >= 2) args[2] else
  Sys.getenv("BAYRC_SUMMARY_DIR",
             unset = file.path(BAYRC_GTEX_DIR, "result", "summary", "hb"))

BF <- 3
P_RHYTHMIC <- 0.2

# Gene naming ------------------------------------------------------------------
# The sampler writes Ensembl IDs (human) or leaves the rows unnamed; the
# downstream scripts want HGNC symbols. The map is the one the 2025 run used,
# shipped as a CSV so nothing here has to call biomaRt and drift from it.

map.file <- Sys.getenv("BAYRC_SYMBOL_MAP",
  unset = file.path(BAYRC_PACKAGE_DIR, "inst", "extdata",
                    "ensembl_symbol_map.csv"))
gene.map <- utils::read.csv(map.file, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(gene.map$ensembl_gene_id),
          !anyDuplicated(gene.map$symbol))

# Rows come out in Ensembl-ID order, which is how the 2025 artifacts are laid
# out and how biomaRt returned them.
gene.map <- gene.map[order(gene.map$ensembl_gene_id), ]

load(file.path(BAYRC_GTEX_DIR, "data", "CAMO.bab.hum.RData"))
tissues <- sort(intersect(names(gtex$CPM.large.clean),
                          names(baboon_withTOD$baboon)))

# BAYRC_TISSUES restricts the run to a few tissues, which is how --validate is
# usually exercised: the layout rule is the same for every tissue, so three of
# them settle whether the builder reproduces the stored artifacts, without
# holding 8 GB of chains in memory.
if (nzchar(Sys.getenv("BAYRC_TISSUES"))) {
  want <- strsplit(Sys.getenv("BAYRC_TISSUES"), "[ ,]+")[[1]]
  if (!all(want %in% tissues)) stop("unknown tissue in BAYRC_TISSUES")
  tissues <- want
}

ens.of <- list(
  human  = lapply(gtex$CPM.large.clean[tissues], rownames),
  baboon = lapply(baboon_withTOD$baboon[tissues], rownames)
)

rds.path <- function(species, tissue) {
  prefix <- if (species == "human") "gtex" else "baboon"
  file.path(result.dir, species, tissue,
            paste0(prefix, "_", tissue, "_bay_1.RDS"))
}

# One chain, every tissue of one species ---------------------------------------

collect <- function(species, what) {
  out <- lapply(tissues, function(tis) {
    f <- rds.path(species, tis)
    if (!file.exists(f)) stop("missing MCMC output: ", f)
    cat(species, what, tis, "\n")
    res <- readRDS(f)
    m   <- res[[what]]
    stopifnot(ncol(m) == 2001, sum(is.na(m)) == 0)

    ens <- ens.of[[species]][[tis]]
    if (!is.null(rownames(m)) && !identical(rownames(m), ens))
      stop("stored rownames disagree with the expression matrix for ",
           species, " ", tis)
    if (length(ens) != nrow(m))
      stop("gene count mismatch for ", species, " ", tis)

    idx <- match(gene.map$ensembl_gene_id, ens)
    if (anyNA(idx))
      stop("genes in the symbol map are missing from ", species, " ", tis)
    m <- m[idx, , drop = FALSE]

    # RHYindex is a property of rho, so it is computed once on rho and reused
    # for phi.
    rhy <- summarize_bay(res$rho[idx, , drop = FALSE], BF, P_RHYTHMIC)$Rhythmicity

    rownames(m) <- gene.map$symbol
    attr(m, "symbols")          <- gene.map$symbol
    attr(m, "RHYindex")         <- rhy
    attr(m, "ensembl_gene_ids") <- gene.map$ensembl_gene_id
    m
  })
  names(out) <- tissues
  out
}

# Compare a rebuilt list against whatever is already on disk -------------------

compare <- function(built, stored, label) {
  ok <- TRUE
  for (tis in names(built)) {
    a <- built[[tis]]; b <- stored[[tis]]
    checks <- c(
      dims    = identical(dim(a), dim(b)),
      rows    = identical(rownames(a), rownames(b)),
      ens     = identical(attr(a, "ensembl_gene_ids"), attr(b, "ensembl_gene_ids")),
      rhy     = identical(attr(a, "RHYindex"), attr(b, "RHYindex")),
      numbers = isTRUE(all.equal(unname(as.matrix(a)), unname(as.matrix(b))))
    )
    if (!all(checks)) {
      ok <- FALSE
      cat(label, tis, "MISMATCH:",
          paste(names(checks)[!checks], collapse = ", "), "\n")
    }
  }
  cat(label, if (ok) "matches the stored artifact\n" else "DIFFERS\n")
  ok
}

# rho --------------------------------------------------------------------------

mcmc_data_human  <- collect("human",  "rho")
mcmc_data_baboon <- collect("baboon", "rho")

if (validate) {
  stored <- new.env()
  load(file.path(summary.dir, "mcmc_rho_BF3.RData"), stored)
  ok.rho <- compare(mcmc_data_human,  stored$mcmc_data_human,  "rho human") &&
            compare(mcmc_data_baboon, stored$mcmc_data_baboon, "rho baboon")
  rm(stored)
} else {
  save(mcmc_data_human, mcmc_data_baboon,
       file = file.path(summary.dir, "mcmc_rho_BF3.RData"))
}
rm(mcmc_data_human, mcmc_data_baboon)
invisible(gc())

# phi --------------------------------------------------------------------------

mcmc_phi_human  <- collect("human",  "phi")
mcmc_phi_baboon <- collect("baboon", "phi")

if (validate) {
  stored <- new.env()
  load(file.path(summary.dir, "phi", "mcmc_phi_BF3.RData"), stored)
  ok.phi <- compare(mcmc_phi_human,  stored$mcmc_phi_human,  "phi human") &&
            compare(mcmc_phi_baboon, stored$mcmc_phi_baboon, "phi baboon")
  rm(stored)
  if (!(ok.rho && ok.phi)) stop("rebuilt summaries do not match the stored ones")
  cat("validation passed\n")
} else {
  dir.create(file.path(summary.dir, "phi"), recursive = TRUE,
             showWarnings = FALSE)
  save(mcmc_phi_human, mcmc_phi_baboon,
       file = file.path(summary.dir, "phi", "mcmc_phi_BF3.RData"))
  cat("wrote", file.path(summary.dir, "mcmc_rho_BF3.RData"), "and",
      file.path(summary.dir, "phi", "mcmc_phi_BF3.RData"), "\n")
}
