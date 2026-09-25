# Summarize comparator-minus-human clock phases from saved case-study caches.
# Usage: Rscript lung_clock_phase_summary.R <result_fixed> <output-directory>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
base <- normalizePath(args[1]); out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
clock <- c("BMAL1", "NR1D1", "NR1D2", "DBP", "PER1", "PER2", "CRY1", "CRY2")
wrap <- function(x) (x + 12) %% 24 - 12
cm <- function(x) wrap(Arg(mean(exp(2i * pi * x / 24))) * 24 / (2*pi))
paths <- c(baboon = "baboon_human_LUN", mouse = "mouse_human_LUN_genome")
rows <- lapply(names(paths), function(sp) {
  f <- file.path(base, "analysis", "figures", paths[[sp]], "plot_data.rds")
  x <- readRDS(f)
  genes <- rownames(x$phase)
  if (is.null(genes)) genes <- names(x$pA)
  stopifnot(length(genes) == length(x$phase$peak1))
  # Both source application caches put comparator in condition A and human in B.
  keep <- genes %in% clock
  data.frame(comparator = sp, gene = genes[keep],
    comparator_peak = x$phase$peak1[keep], human_peak = x$phase$peak2[keep],
    difference = wrap(x$phase$peak1[keep] - x$phase$peak2[keep]),
    rhythm_conserved = as.character(x$trans$gain_loss_status)[keep] == "Maintained",
    measured_genes = length(genes), source = f)
})
d <- do.call(rbind, rows)
write.csv(d, file.path(out, "lung_clock_phase_genes.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(d, d$comparator), function(z) {
  x <- z$difference[z$rhythm_conserved & is.finite(z$difference)]
  data.frame(comparator = z$comparator[1], n = length(x),
    circular_mean_h = cm(x), median_h = median(x),
    mean_absolute_h = mean(abs(x)), minimum_h = min(x), maximum_h = max(x))
}))
write.csv(summary, file.path(out, "lung_clock_phase_summary.csv"), row.names = FALSE)
print(summary)
writeLines(c(paste("Generated:", Sys.time()), paste("Result directory:", base),
  "Clock genes: BMAL1 NR1D1 NR1D2 DBP PER1 PER2 CRY1 CRY2",
  "Selection: cached Maintained transition classification (BFDR 0.25 in manuscript runs)",
  "Peak estimates: cached phase_infer peak1 and peak2, comparator minus human",
  "Differences wrapped to [-12,12); circular mean also wrapped to [-12,12)",
  capture.output(sessionInfo())), file.path(out, "lung_clock_phase_run_record.txt"))
