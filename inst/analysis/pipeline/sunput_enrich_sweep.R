## Sweep the pathway-list choice, size filter, stage-1 cut and stage-2 q cut
## over the stored pathSelect tables and report what each cell selects.
##
## Usage: Rscript sunput_enrich_sweep.R [nperm]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
nperm <- if (length(args) >= 1) as.integer(args[1]) else 10000
out.dir <- file.path(BAYRC_OUTPUT_DIR, "param_search")
raw <- readRDS(file.path(out.dir, sprintf("pathselect_raw_nperm%d.rds", nperm)))

## pathway names carrying the mitochondrial and proteostasis story
key <- c("KEGG Oxidative phosphorylation", "KEGG Proteasome",
         "KEGG Parkinson's disease", "KEGG Parkinson disease",
         "KEGG Huntington's disease", "KEGG Huntington disease",
         "KEGG Alzheimer's disease", "KEGG Alzheimer disease")

size.rules <- list("raw5-500"  = function(d) d$raw_size >= 5   & d$raw_size <= 500,
                   "raw10-200" = function(d) d$raw_size >= 10  & d$raw_size <= 200,
                   "raw15-500" = function(d) d$raw_size >= 15  & d$raw_size <= 500,
                   "meas10+"   = function(d) d$measured >= 10,
                   "meas15+"   = function(d) d$measured >= 15,
                   "raw5-500+meas15" = function(d) d$raw_size >= 5 &
                     d$raw_size <= 500 & d$measured >= 15,
                   "none"      = function(d) rep(TRUE, nrow(d)))
stage1.rules <- list("p<0.01" = function(u) u$pval < 0.01,
                     "p<0.05" = function(u) u$pval < 0.05,
                     "p<0.10" = function(u) u$pval < 0.10,
                     "p<0.20" = function(u) u$pval < 0.20,
                     "q<0.05" = function(u) p.adjust(u$pval, "BH") < 0.05,
                     "q<0.20" = function(u) p.adjust(u$pval, "BH") < 0.20)
qcuts <- c(0.05, 0.10, 0.20, 0.25)

rows <- list(); hits <- list()
for (ln in c("kegg229", "kegg354")) {
  for (sn in names(size.rules)) {
    u <- raw[[paste0(ln, "_union")]]
    u <- u[size.rules[[sn]](u), ]
    for (s1 in names(stage1.rules)) {
      active <- u$pathway[stage1.rules[[s1]](u)]
      for (q in qcuts) {
        sig <- list()
        for (rk in c("gain", "loss", "conserved")) {
          t2 <- raw[[paste0(ln, "_", rk)]]
          t2 <- t2[t2$pathway %in% active, ]
          if (!nrow(t2)) { sig[[rk]] <- character(0); next }
          t2$q <- p.adjust(t2$pval, "BH")
          sig[[rk]] <- t2$pathway[t2$q < q]
        }
        allsig <- unique(unlist(sig))
        rows[[length(rows) + 1]] <- data.frame(
          list = ln, size = sn, stage1 = s1, qcut = q,
          tested = nrow(u), active = length(active),
          gain = length(sig$gain), loss = length(sig$loss),
          cons = length(sig$conserved), any = length(allsig),
          oxphos_cons = "KEGG Oxidative phosphorylation" %in% sig$conserved,
          proteasome_loss = "KEGG Proteasome" %in% sig$loss,
          parkinson_cons = any(grepl("Parkinson", sig$conserved)),
          neurodeg_cons = sum(grepl("Parkinson|Huntington|Alzheimer",
                                    sig$conserved)))
        hits[[sprintf("%s|%s|%s|%g", ln, sn, s1, q)]] <- sig
      }
    }
  }
}

d <- do.call(rbind, rows)
write.csv(d, file.path(out.dir, "enrich_sweep.csv"), row.names = FALSE)
saveRDS(hits, file.path(out.dir, "enrich_sweep_hits.rds"))

cat("\n=== stage-1 union p-values for the key pathways ===\n")
for (ln in c("kegg229", "kegg354")) {
  u <- raw[[paste0(ln, "_union")]]
  k <- u[u$pathway %in% key, c("pathway", "raw_size", "measured", "pval", "padj")]
  if (nrow(k)) { cat("\n--", ln, "--\n"); print(k, row.names = FALSE) }
}
cat("\n=== stage-2 raw p-values for the key pathways (full list) ===\n")
for (ln in c("kegg229", "kegg354")) for (rk in c("gain", "loss", "conserved")) {
  t2 <- raw[[paste0(ln, "_", rk)]]
  k <- t2[t2$pathway %in% key, c("pathway", "measured", "pval")]
  if (nrow(k)) { cat("\n--", ln, rk, "--\n"); print(k, row.names = FALSE) }
}
cat("\n=== sweep ===\n")
print(d, row.names = FALSE)
