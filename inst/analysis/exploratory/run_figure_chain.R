################################################################################
# Drive the case-study scripts that produce the Figure 3-6 panels.
#
# The companion plotting scripts read phase_inner, trans_outer, the
# baboon_*/human_* objects and output.dir from their case study's workspace, so
# each group is sourced in order into one session, one R process per group.
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

# Driver state uses dotted names, which rm(list = ls()) in the case studies skips
.bayrc_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
.bayrc_dir  <- if (is.na(.bayrc_file)) getwd() else
  dirname(normalizePath(.bayrc_file))

.bayrc_groups <- list(
  scn_hip = c("Baboon_SCN_HIP.R", "plots/Peak_new_conserved.R"),
  sun_put = c("Baboon_SUN_PUT.R", "plot_enrich_SUN_PUT.R"),
  hb_lun  = c("Baboon_Human_LUN.R", "plots/heatmap_circadian_concordance.R")
)

.bayrc_args <- commandArgs(trailingOnly = TRUE)
.bayrc_want <- if (length(.bayrc_args) && nzchar(.bayrc_args[1]))
  .bayrc_args[1] else NA_character_

if (is.na(.bayrc_want)) {
  # one process per group, so a failure in one does not take the others with it
  for (.g in names(.bayrc_groups)) {
    cat("\n================ group", .g, "================\n")
    .st <- system2("Rscript", c(shQuote(normalizePath(.bayrc_file)), .g))
    cat("group", .g, "exited", .st, "\n")
  }
  quit(save = "no")
}

if (!.bayrc_want %in% names(.bayrc_groups))
  stop("group must be one of: ", paste(names(.bayrc_groups), collapse = ", "))

# The case studies change directory, so the start directory is restored per script
.bayrc_home <- getwd()

for (.s in .bayrc_groups[[.bayrc_want]]) {
  .f <- file.path(.bayrc_dir, .s)
  if (!file.exists(.f)) {
    cat("MISSING", .s, "\n")
    next
  }
  cat("\n---- sourcing", .s, "at", format(Sys.time(), "%H:%M:%S"), "----\n")
  setwd(.bayrc_home)
  .t0 <- Sys.time()
  # companions read the case study's objects from the global environment
  source(.f, echo = FALSE, local = FALSE)
  cat("---- finished", .s, "in",
      round(as.numeric(difftime(Sys.time(), .t0, units = "mins")), 1),
      "min ----\n")
}

setwd(.bayrc_home)
cat("\ngroup", .bayrc_want, "complete\n")
