################################################################################
# Drive the concordance heatmap scripts behind Figure 2.
#
# The three scripts call multi_conservation() and pheatmap() but load neither,
# expecting a session that already has them. This attaches the package first
# and then sources each one in its own process, so a failure in one does not
# take the others with it.
#
# Usage:
#   Rscript run_figure2_chain.R [script]
#
#   baboon     plots/heatmap_baboon.R           within-baboon, 325 pairs
#   human      plots/heatmap_human.R            within-human, 325 pairs
#   circadian  plots/heatmap_circadian_pairs.R  KEGG Circadian subset
#
# With no argument all three run in turn.
################################################################################

.bayrc_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
.bayrc_dir  <- if (is.na(.bayrc_file)) getwd() else
  dirname(normalizePath(.bayrc_file))

.bayrc_scripts <- c(
  baboon    = "plots/heatmap_baboon.R",
  human     = "plots/heatmap_human.R",
  circadian = "plots/heatmap_circadian_pairs.R")

# heatmap_circadian_pairs.R reads its mode from commandArgs() unless `mode`
# already exists, and under a driver those arguments are the driver's own.
.bayrc_modes <- c(circadian = "within_baboon")

.bayrc_args <- commandArgs(trailingOnly = TRUE)
.bayrc_want <- if (length(.bayrc_args) && nzchar(.bayrc_args[1]))
  .bayrc_args[1] else NA_character_

if (is.na(.bayrc_want)) {
  for (.s in names(.bayrc_scripts)) {
    cat("\n================", .s, "================\n")
    cat("exited", system2("Rscript", c(shQuote(normalizePath(.bayrc_file)), .s)), "\n")
  }
  quit(save = "no")
}

if (!.bayrc_want %in% names(.bayrc_scripts))
  stop("script must be one of: ", paste(names(.bayrc_scripts), collapse = ", "))

suppressPackageStartupMessages(library(BayRC))

.bayrc_home <- getwd()
if (.bayrc_want %in% names(.bayrc_modes)) {
  mode <- unname(.bayrc_modes[.bayrc_want])
  cat("mode:", mode, "\n")
}
.f <- file.path(.bayrc_dir, .bayrc_scripts[[.bayrc_want]])
cat("---- sourcing", .bayrc_scripts[[.bayrc_want]], "at",
    format(Sys.time(), "%H:%M:%S"), "----\n")
.t0 <- Sys.time()
source(.f, echo = FALSE, local = FALSE)
setwd(.bayrc_home)
cat("---- finished in",
    round(as.numeric(difftime(Sys.time(), .t0, units = "mins")), 1), "min ----\n")
