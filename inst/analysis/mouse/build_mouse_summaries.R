## Assemble the human/mouse rho and phi summaries from the six mouse tissue
## chains and the matching human chains.
##
## Same layout rule as pipeline/summarize_rho_phi.R: rows in Ensembl-ID order,
## named by HGNC symbol, carrying the symbols / RHYindex / ensembl_gene_ids
## attributes. The gene set is the three-way match, the genes the symbol map
## shares with the mouse-to-human object, so human, baboon and mouse can be put
## on one set of rows.
##
## The human side is taken from the human/baboon summary rather than sampled
## again. The sampler treats genes independently: every prior is a fixed input,
## p_rhythmic is held at 0.2 rather than estimated across genes, and each
## update is elementwise over rows, so a gene's posterior does not depend on
## which other genes share the matrix.

args    <- commandArgs(trailingOnly = TRUE)
resdir  <- args[1]
hbdir   <- args[2]
outdir  <- args[3]
mapfile <- args[4]

BF <- 3
P_RHYTHMIC <- 0.2

dir.create(file.path(outdir, "phi"), recursive = TRUE, showWarnings = FALSE)

gene.map <- utils::read.csv(mapfile, stringsAsFactors = FALSE)

# CAMO.mouse.hum.RData sits under BAYRC_GTEX_DIR/data, beside CAMO.bab.hum.RData
datafile <- file.path(Sys.getenv("BAYRC_GTEX_DIR"), "data", "CAMO.mouse.hum.RData")
if (!nzchar(Sys.getenv("BAYRC_GTEX_DIR")) || !file.exists(datafile))
  stop("set BAYRC_GTEX_DIR to the directory holding data/CAMO.mouse.hum.RData")
load(datafile)
tissues <- sort(names(mice$count_clean))

## The three-way set: the symbol map cut to the genes the mouse object carries.
mouse.ens <- lapply(mice$count_clean[tissues], rownames)
gene.map <- gene.map[gene.map$ensembl_gene_id %in% mouse.ens[[1]], ]
gene.map <- gene.map[order(gene.map$ensembl_gene_id), ]
cat("three-way matched genes:", nrow(gene.map), "\n")

rm(mice, gtex); invisible(gc())

summarize_bay_rhy <- function(rho) {
  p <- rowMeans(rho)
  bf <- (p/(1 - p + 1e-20)) / (P_RHYTHMIC/(1 - P_RHYTHMIC))
  as.integer(bf > BF)
}

lay_out <- function(m, rhy) {
  rownames(m) <- gene.map$symbol
  attr(m, "symbols")          <- gene.map$symbol
  attr(m, "RHYindex")         <- rhy
  attr(m, "ensembl_gene_ids") <- gene.map$ensembl_gene_id
  m
}

## Mouse, read from the per-tissue chains ---------------------------------------

collect_mouse <- function(what) {
  out <- lapply(tissues, function(tis) {
    f <- file.path(resdir, "mice", tis, paste0("mice_", tis, "_bay_1.RDS"))
    if (!file.exists(f)) stop("missing MCMC output: ", f)
    cat("mouse", what, tis, "\n")
    res <- readRDS(f)
    ens <- mouse.ens[[tis]]
    if (length(ens) != nrow(res$rho))
      stop("gene count mismatch for mouse ", tis)
    idx <- match(gene.map$ensembl_gene_id, ens)
    if (anyNA(idx)) stop("symbol map does not cover mouse ", tis)
    stopifnot(ncol(res$rho) == 2001, sum(is.na(res$rho)) == 0)
    rhy <- summarize_bay_rhy(res$rho[idx, , drop = FALSE])
    m <- lay_out(res[[what]][idx, , drop = FALSE], rhy)
    rm(res); invisible(gc())
    m
  })
  names(out) <- tissues
  out
}

## Human, cut from the human/baboon summary ------------------------------------

collect_human <- function(file, object) {
  e <- new.env()
  load(file, envir = e)
  stored <- get(object, envir = e)
  out <- lapply(tissues, function(tis) {
    cat("human", object, tis, "\n")
    m <- stored[[tis]]
    idx <- match(gene.map$ensembl_gene_id, attr(m, "ensembl_gene_ids"))
    if (anyNA(idx)) stop("symbol map does not cover human ", tis)
    stopifnot(ncol(m) == 2001, sum(is.na(m)) == 0)
    lay_out(m[idx, , drop = FALSE], attr(m, "RHYindex")[idx])
  })
  names(out) <- tissues
  rm(stored, e); invisible(gc())
  out
}

mcmc_data_mouse <- collect_mouse("rho")
mcmc_data_human <- collect_human(file.path(hbdir, "mcmc_rho_BF3.RData"),
                                 "mcmc_data_human")
save(mcmc_data_human, mcmc_data_mouse,
     file = file.path(outdir, "mcmc_rho_BF3.RData"))
rm(mcmc_data_human, mcmc_data_mouse); invisible(gc())

mcmc_phi_mouse <- collect_mouse("phi")
mcmc_phi_human <- collect_human(file.path(hbdir, "phi", "mcmc_phi_BF3.RData"),
                                "mcmc_phi_human")
save(mcmc_phi_human, mcmc_phi_mouse,
     file = file.path(outdir, "phi", "mcmc_phi_BF3.RData"))

record <- c(
  paste("script       ", "src/build_mouse_summaries.R"),
  paste("mode         ", "full run"),
  paste("BayRC version", as.character(utils::packageVersion("BayRC"))),
  paste("mouse chains ", file.path(resdir, "mice")),
  paste("human source ", hbdir),
  paste("symbol map   ", mapfile),
  paste("summaries    ", outdir),
  paste("run at       ", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "parameters",
  sprintf("  %-14s %s", "tissues", paste(tissues, collapse = ", ")),
  sprintf("  %-14s %s", "genes", nrow(gene.map)),
  sprintf("  %-14s %s", "draws", 2001),
  sprintf("  %-14s %s", "BF", BF),
  sprintf("  %-14s %s", "p_rhythmic", P_RHYTHMIC))
writeLines(record, file.path(outdir, "run_record.txt"))

cat("wrote human/mouse summaries to", outdir, "\n")
