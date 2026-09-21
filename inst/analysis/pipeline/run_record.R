# Record which MCMC run and which parameters produced a set of figures.

# Fields of a record written earlier, for a replot to check and to cite.
read_run_record <- function(path) {
  if (!file.exists(path)) return(NULL)
  lines <- readLines(path)
  field <- function(key) {
    hit <- grep(paste0("^", key), lines, value = TRUE)[1]
    if (is.na(hit)) NA_character_ else trimws(sub(paste0("^", key), "", hit))
  }
  list(summaries = field("summaries"), rho = field("rho written"),
       run_at = field("run at"))
}

# A replot draws from tables a full run left behind, so it refuses to run
# against summaries other than the ones those tables came from.
require_run_record <- function(fig.dir, tables) {
  path <- file.path(fig.dir, "run_record.txt")
  record <- read_run_record(path)
  if (is.null(record))
    stop("--replot: no run_record.txt under ", fig.dir,
         "\n  run the script once without --replot first.")
  missing <- tables[!file.exists(file.path(fig.dir, tables))]
  if (length(missing))
    stop("--replot: missing under ", fig.dir, ": ", paste(missing, collapse = ", "),
         "\n  run the script once without --replot first.")
  if (!identical(record$summaries, BAYRC_SUMMARY_DIR))
    stop("--replot: those tables came from ", record$summaries,
         "\n  but config.R resolves ", BAYRC_SUMMARY_DIR, ".")
  record
}

write_run_record <- function(path, script, params, repo = NULL, replot_of = NULL) {
  rho <- file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData")
  commit <- NA_character_
  if (!is.null(repo))
    commit <- tryCatch(system2("git", c("-C", repo, "rev-parse", "--short", "HEAD"),
                               stdout = TRUE, stderr = NULL)[1],
                       warning = function(w) NA_character_,
                       error = function(e) NA_character_)
  lines <- c(
    paste("script       ", script),
    paste("mode         ", if (is.null(replot_of)) "full run" else
                           paste("replot of the run at", replot_of)),
    paste("commit       ", commit),
    paste("BayRC version", as.character(utils::packageVersion("BayRC"))),
    paste("summaries    ", BAYRC_SUMMARY_DIR),
    paste("rho written  ", format(file.info(rho)$mtime, "%Y-%m-%d %H:%M")),
    paste("run at       ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    "parameters",
    sprintf("  %-14s %s", names(params),
            vapply(params, function(x) paste(x, collapse = ", "), character(1))))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path)
  cat("run record:", path, "\n")
}
