################################################################################
# Drive the case-study scripts that produce the Figure 3-6 panels.
#
# Four of the plotting scripts are not standalone: they read phase_inner,
# trans_outer, the baboon_*/human_* objects and output.dir straight out of the
# workspace their case study left behind. Running them on their own stops with
# an error, and running them after the wrong case study silently plots the
# wrong comparison. So each case study and its companions are sourced into one
# environment here, in order, and each group gets a fresh R process so nothing
# leaks between comparisons.
#
# Usage:
#   Rscript run_figure_chain.R [group]
#
#   scn_hip    Baboon_SCN_HIP.R      -> Figure 3 panel A
#              plots/Peak_new_conserved.R
#   sun_put    Baboon_SUN_PUT.R      -> Figure 3 panel B, Figure 5
#              plot_enrich_SUN_PUT.R -> Figure 4
#   hb_lun     Baboon_Human_LUN.R    -> Figure 6
#              plots/heatmap_circadian_concordance.R
#
# With no argument every group runs, one process each.
################################################################################

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))

groups <- list(
  scn_hip = c("Baboon_SCN_HIP.R", "plots/Peak_new_conserved.R"),
  sun_put = c("Baboon_SUN_PUT.R", "plot_enrich_SUN_PUT.R"),
  hb_lun  = c("Baboon_Human_LUN.R", "plots/heatmap_circadian_concordance.R")
)

args <- commandArgs(trailingOnly = TRUE)
want <- if (length(args) && nzchar(args[1])) args[1] else NA_character_

if (is.na(want)) {
  # one process per group, so a failure in one does not take the others with it
  for (g in names(groups)) {
    cat("\n================ group", g, "================\n")
    st <- system2("Rscript", c(shQuote(normalizePath(this.file)), g))
    cat("group", g, "exited", st, "\n")
  }
  quit(save = "no")
}

if (!want %in% names(groups))
  stop("group must be one of: ", paste(names(groups), collapse = ", "))

# The case-study scripts setwd() partway through and never fully restore it,
# so the starting directory is captured and put back between scripts.
here <- getwd()
on.exit(setwd(here), add = TRUE)

for (s in groups[[want]]) {
  f <- file.path(analysis.dir, s)
  if (!file.exists(f)) {
    cat("MISSING", s, "\n")
    next
  }
  cat("\n---- sourcing", s, "at", format(Sys.time(), "%H:%M:%S"), "----\n")
  setwd(here)
  t0 <- Sys.time()
  # the companions need the case study's objects, so everything shares one
  # environment: the global one
  source(f, echo = FALSE, local = FALSE)
  cat("---- finished", s, "in",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min ----\n")
}

setwd(here)
cat("\ngroup", want, "complete\n")
