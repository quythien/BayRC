# Everything a figure needs that sits downstream of the MCMC and upstream of the
# palette, so --replot draws without loading the draws or re-running inference.
# Posterior draws are kept only for the genes a panel draws.

plot_cache_stamp <- function() {
  rho <- file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData")
  list(summaries = BAYRC_SUMMARY_DIR,
       rho_written = format(file.info(rho)$mtime, "%Y-%m-%d %H:%M"))
}

# panel_genes restricts the draws; everything else is per-gene and already small
write_plot_cache <- function(path, dataA, dataB, panel_genes, ...) {
  keep <- intersect(rownames(dataA$rho), panel_genes)
  cache <- c(plot_cache_stamp(),
             list(panelA = list(rho = dataA$rho[keep, , drop = FALSE],
                                phi = dataA$phi[keep, , drop = FALSE]),
                  panelB = list(rho = dataB$rho[keep, , drop = FALSE],
                                phi = dataB$phi[keep, , drop = FALSE])),
             list(...))
  saveRDS(cache, path, compress = "xz")
  cat("plot cache:", path, "\n")
}

# refuses a cache written against another summary directory or write date
read_plot_cache <- function(path) {
  if (!file.exists(path))
    stop("--replot: no plot cache at ", path,
         "\n  run the script once without --replot first.")
  cache <- readRDS(path)
  now <- plot_cache_stamp()
  if (!identical(cache$summaries, now$summaries))
    stop("--replot: the cache was written against ", cache$summaries,
         "\n  but config.R resolves ", now$summaries, ".")
  if (!identical(cache$rho_written, now$rho_written))
    stop("--replot: the cache was written against summaries dated ",
         cache$rho_written, "\n  but the current ones are dated ",
         now$rho_written, ".")
  cache
}
