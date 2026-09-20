# Record which MCMC run and which parameters produced a set of figures.

write_run_record <- function(path, script, params, repo = NULL) {
  rho <- file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData")
  commit <- NA_character_
  if (!is.null(repo))
    commit <- tryCatch(system2("git", c("-C", repo, "rev-parse", "--short", "HEAD"),
                               stdout = TRUE, stderr = NULL)[1],
                       warning = function(w) NA_character_,
                       error = function(e) NA_character_)
  lines <- c(
    paste("script       ", script),
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
