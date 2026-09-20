## Frozen two-stage enrichment run over several baboon tissue pairs, for
## comparing candidate case studies. Same setting as the paper: 354-pathway
## release, gene sets cut to the measured genes and kept at 15 or more, union
## stage 1 at BH q < 0.20, transition stage 2 at q < 0.05, 10000 permutations.
##
## Usage: Rscript pair_enrichment_scan.R [A-B A-B ...]

suppressPackageStartupMessages({library(BayRC); library(dplyr)})

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

STAGE1_Q <- 0.20
STAGE2_Q <- 0.05
NPERM    <- 10000
MIN_MEASURED <- 15
SHIFT <- 2
BFDR_ALPHA <- 0.25

args <- commandArgs(trailingOnly = TRUE)
if (!length(args))
  args <- c("HEA-THR", "STF-THR", "ILE-THR", "HEA-STF",
            "ILE-STF", "OMF-STF", "PUT-STF", "PUT-THR")
pairs <- strsplit(args, "-")

out.dir <- file.path(BAYRC_OUTPUT_DIR, "pair_scan")
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))
measured <- rownames(mcmc_data_baboon[[1]])

kegg <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))
kegg <- lapply(kegg, function(g) replace(g, g == "ARNTL", "BMAL1"))
kegg <- lapply(kegg, intersect, measured)
kegg <- kegg[lengths(kegg) >= MIN_MEASURED]

summary.rows <- list(); enriched.rows <- list()
for (pr in pairs) {
  tA <- pr[1]; tB <- pr[2]
  cat("\n########", tA, "vs", tB, "########\n")
  datA <- list(rho = mcmc_data_baboon[[tA]], phi = mcmc_phi_baboon[[tA]])
  datB <- list(rho = mcmc_data_baboon[[tB]], phi = mcmc_phi_baboon[[tB]])

  pA <- rowMeans(datA$rho); pB <- rowMeans(datB$rho)
  invisible(capture.output(tr <- transition_classify(pA, pB, bfdr_alpha = BFDR_ALPHA)))
  st <- tr$gain_loss_status
  maint <- names(st)[st == "Maintained"]

  ph <- phase_infer(phi_matrix1 = datA$phi, phi_matrix2 = datB$phi,
                    gain_loss_status = st, bfdr_alpha = BFDR_ALPHA,
                    shift = SHIFT, P = 24, compute_hdi = FALSE)
  delta <- ((ph$peak2 - ph$peak1 + 12) %% 24) - 12
  d.maint <- delta[maint]

  run <- function(plist, method)
    pathSelect(mcmc.merge.list = setNames(list(datA, datB), c(tA, tB)),
               pathway.list = plist, dataset.names = c(tA, tB),
               ranking.method = method, score_type = "pos", qvalue.cut = STAGE2_Q,
               pathwaysize.lower.cut = MIN_MEASURED,
               pathwaysize.upper.cut = length(measured),
               nperm = NPERM, nproc = 1)$results

  u <- run(kegg, "union"); u$q <- p.adjust(u$pval, "BH")
  active <- u$pathway[u$q < STAGE1_Q]

  sig <- NULL
  if (length(active)) {
    s2 <- do.call(rbind, lapply(c("gain", "loss", "conserved"), function(m) {
      r <- run(kegg[active], m)[, c("pathway", "size", "pval")]
      r$q <- p.adjust(r$pval, "BH"); r$direction <- m; r
    }))
    sig <- s2[s2$q < STAGE2_Q, ]
  }
  sel <- unique(sig$pathway)

  ## redundancy of the enriched set, distinct genes against summed sizes
  slots <- if (length(sel)) sum(lengths(kegg[sel])) else 0
  distinct <- if (length(sel)) length(unique(unlist(kegg[sel]))) else 0

  summary.rows[[length(summary.rows) + 1]] <- data.frame(
    pair = paste(tA, tB, sep = "-"),
    maintained = length(maint), gain = sum(st == "Gain"), loss = sum(st == "Loss"),
    shifted = sum(ph$flag_shift), conserved = sum(ph$flag_cons),
    offset = round(mean(d.maint), 2), pct_sign = round(100 * max(mean(d.maint > 0),
                                                                 mean(d.maint < 0)), 1),
    stage1 = length(active), enriched = length(sel),
    n_rows = if (is.null(sig)) 0 else nrow(sig),
    slots = slots, distinct = distinct,
    redundancy = if (slots) round(1 - distinct / slots, 2) else NA)
  if (length(sel)) {
    sig$pair <- paste(tA, tB, sep = "-")
    enriched.rows[[length(enriched.rows) + 1]] <- sig
  }
  print(summary.rows[[length(summary.rows)]], row.names = FALSE)
}

sm <- do.call(rbind, summary.rows)
en <- if (length(enriched.rows)) do.call(rbind, enriched.rows) else NULL
write.csv(sm, file.path(out.dir, "pair_scan_summary.csv"), row.names = FALSE)
if (!is.null(en)) write.csv(en[order(en$pair, en$direction, en$q), ],
                            file.path(out.dir, "pair_scan_enriched.csv"), row.names = FALSE)

cat("\n=== summary ===\n"); print(sm, row.names = FALSE)
if (!is.null(en)) { cat("\n=== enriched pathways ===\n")
  print(en[order(en$pair, en$direction, en$q),
           c("pair", "pathway", "direction", "size", "q")], row.names = FALSE) }
