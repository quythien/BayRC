# tests/testthat/helper-source-package.R
# Sources all R package files so unexported functions are available in tests.
# This helper is automatically loaded by testthat before any test file.

pkg_root <- system.file(package = "BayRC")

# If running via devtools::test() or R CMD check, functions are already loaded.
# Source manually only when running interactively or the package is not yet
# installed (e.g., during development with devtools::load_all()).

if (!isNamespaceLoaded("BayRC")) {
  r_dir <- file.path(dirname(dirname(dirname(getwd()))), "R")
  # Fallback: use the known package root
  pkg_r_dir <- "/home/qtp1/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/BayRC/R"
  if (dir.exists(pkg_r_dir)) {
    scripts <- list.files(pkg_r_dir, pattern = "\\.R$", full.names = TRUE)
    # Source in dependency order (utils first, then cosinor, then heavier files)
    priority <- c("utils.R", "cosinor.R", "estimates.R", "clustering.R", "dpm.R",
                  "simulate.R", "init.R", "init_time.R", "permutation_sim.R",
                  "internal.R", "global.R", "pathway.R", "merge.R",
                  "path_select.R", "mcmc.R", "mcmc_time.R")
    ordered <- c(
      file.path(pkg_r_dir, priority[priority %in% basename(scripts)]),
      scripts[!basename(scripts) %in% priority]
    )
    ordered <- ordered[file.exists(ordered)]
    for (s in ordered) {
      tryCatch(source(s, local = FALSE), error = function(e) {
        message("Could not source ", basename(s), ": ", conditionMessage(e))
      })
    }
  }
}
