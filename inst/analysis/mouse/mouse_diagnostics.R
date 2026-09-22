# Every mouse number the paper and its supplement quote, from one script.
#
# Reads the chains under BAYRC_RESULT_DIR and the summaries beside them, and
# prints each block in the order the manuscript uses it. Run after
# run_mouse.R, run_mouse_cycle.R and build_mouse_summaries.R.

suppressPackageStartupMessages(library(BayRC))

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

tissues <- c("AOR", "CER", "HEA", "KIC", "LIV", "LUN")
clock   <- c("BMAL1", "NR1D1", "NR1D2", "DBP", "PER1", "PER2", "CRY1", "CRY2")
alpha   <- 0.25
arms <- list(pooled = file.path(BAYRC_RESULT_DIR, "mice"),
             cycle1 = file.path(BAYRC_RESULT_DIR, "mice_cycle", "cycle1"),
             cycle2 = file.path(BAYRC_RESULT_DIR, "mice_cycle", "cycle2"))

circ_mean <- function(h) {
  a <- 2 * pi * h / 24
  (atan2(mean(sin(a)), mean(cos(a))) %% (2 * pi)) * 24 / (2 * pi)
}
wrap <- function(d) { d <- d %% 24; ifelse(d > 12, d - 24, d) }

hm <- new.env()
load(file.path(BAYRC_RESULT_DIR, "summary", "hm", "mcmc_rho_BF3.RData"), envir = hm)
hmp <- new.env()
load(file.path(BAYRC_RESULT_DIR, "summary", "hm", "phi", "mcmc_phi_BF3.RData"),
     envir = hmp)

## 1. rhythmic content of the reference and the mouse comparator -------------
cat("\n== rhythmic genes of", nrow(hm$mcmc_data_mouse[[1]]), "at BFDR", alpha, "==\n")
for (t in tissues) {
  h <- rowMeans(hm$mcmc_data_human[[t]]); m <- rowMeans(hm$mcmc_data_mouse[[t]])
  cat(sprintf("  %s  human %4d (%.1f%%)   mouse %4d (%.1f%%)\n", t,
              sum(bfdr_from_posterior(h, alpha = alpha)$rhythmic_genes),
              100 * mean(bfdr_from_posterior(h, alpha = alpha)$rhythmic_genes),
              sum(bfdr_from_posterior(m, alpha = alpha)$rhythmic_genes),
              100 * mean(bfdr_from_posterior(m, alpha = alpha)$rhythmic_genes)))
}

## 2. the two cycles against the pooled series ------------------------------
cat("\n== inclusion and core-clock amplitude, by arm ==\n")
sym <- setNames(rownames(hm$mcmc_data_mouse[[1]]),
                attr(hm$mcmc_data_mouse[[1]], "ensembl_gene_ids"))
e <- new.env(); load(file.path(BAYRC_GTEX_DIR, "data", "CAMO.mouse.hum.RData"), envir = e)
for (t in tissues) {
  for (arm in names(arms)) {
    f <- list.files(file.path(arms[[arm]], t), pattern = "\\.RDS$", full.names = TRUE)
    if (!length(f)) next
    x <- readRDS(f[1])
    ids <- unname(sym[rownames(e$mice$count_clean[[t]])])
    i <- match(clock, ids); i <- i[!is.na(i)]
    rho <- rowMeans(x$rho)
    cat(sprintf("  %s %-7s rhythmic %4d (%.1f%%)  clock A %.3f  clock rho %.3f\n",
                t, arm, sum(bfdr_from_posterior(rho, alpha = alpha)$rhythmic_genes),
                100 * mean(bfdr_from_posterior(rho, alpha = alpha)$rhythmic_genes),
                mean(rowMeans(x$A)[i]), mean(rho[i])))
    rm(x); gc(verbose = FALSE)
  }
}

## 3. peak time of each species, and the offsets the figure draws -----------
cat("\n== peak time and offset against the reference, lung ==\n")
hb <- new.env()
load(file.path(BAYRC_RESULT_DIR, "summary", "hb", "phi", "mcmc_phi_BF3.RData"),
     envir = hb)
cat(sprintf("  %-7s %8s %8s %8s %10s %10s\n", "gene", "human", "baboon",
            "mouse", "bab-hum", "mou-hum"))
for (g in clock) {
  h <- if (g %in% rownames(hb$mcmc_phi_human$LUN)) circ_mean(hb$mcmc_phi_human$LUN[g, ]) else NA
  b <- if (g %in% rownames(hb$mcmc_phi_baboon$LUN)) circ_mean(hb$mcmc_phi_baboon$LUN[g, ]) else NA
  m <- if (g %in% rownames(hmp$mcmc_phi_mouse$LUN)) circ_mean(hmp$mcmc_phi_mouse$LUN[g, ]) else NA
  cat(sprintf("  %-7s %8.2f %8.2f %8.2f %+10.2f %+10.2f\n", g, h, b, m,
              wrap(b - h), wrap(m - h)))
}
cat("\n  The figure and the classification take these from phase_infer's own\n",
    " peak1 and peak2, which the case study reports; the circular means above\n",
    " are a separate summary of the same posterior.\n", sep = "")

## 4. gene-to-gene intervals, which carry no origin -------------------------
cat("\n== peak measured from BMAL1 within each species ==\n")
for (t in c("LUN", "LIV")) {
  ok <- Reduce(intersect, list(clock, rownames(hb$mcmc_phi_human[[t]]),
                               rownames(hmp$mcmc_phi_mouse[[t]])))
  ph <- sapply(ok, function(g) circ_mean(hb$mcmc_phi_human[[t]][g, ]))
  pm <- sapply(ok, function(g) circ_mean(hmp$mcmc_phi_mouse[[t]][g, ]))
  d <- wrap(wrap(ph - ph["BMAL1"]) - wrap(pm - pm["BMAL1"]))[setdiff(ok, "BMAL1")]
  cat(sprintf("  %s  mean |interval difference| %.2f h, within 2 h %d of %d\n",
              t, mean(abs(d)), sum(abs(d) <= 2), length(d)))
}

## 5. the gene universes the two pairs are measured on ----------------------
# The mouse summary is cut to the genes the baboon pair also carries, so both
# of its sides are the three-species set. run_mouse.R reports the wider set the
# mouse chains themselves were run on.
hbr <- new.env()
load(file.path(BAYRC_RESULT_DIR, "summary", "hb", "mcmc_rho_BF3.RData"), envir = hbr)
cat("\n== gene universes ==\n")
cat(sprintf("  human and baboon      %5d\n", nrow(hbr$mcmc_data_human[[1]])))
cat(sprintf("  human, baboon, mouse  %5d   what a three-species panel can draw\n",
            nrow(hm$mcmc_data_human[[1]])))
