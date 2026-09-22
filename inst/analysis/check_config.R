# What every analysis script resolves its paths to, and whether they exist.
#
# Run this first on a new machine. It prints each BAYRC_* path config.R
# settles on, reports the ones that are missing, and checks that the two
# guards behave: bayrc_file() names the environment variable to set, and a
# script that arranges existing files can opt out of the summary requirement.

bayrc.needs.summary <- FALSE
this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
source(file.path(analysis.dir, "config.R"))

cat("paths\n")
vars <- sort(grep("^BAYRC_", ls(), value = TRUE))
missing <- character()
for (v in vars) {
  p <- get(v)
  ok <- file.exists(p)
  if (!ok) missing <- c(missing, v)
  cat(sprintf("  %-22s %-4s %s\n", v, if (ok) "ok" else "--", p))
}

cat("\nguards\n")
cat("  bayrc_file on a missing file: ",
    tryCatch({ bayrc_file(BAYRC_PATHWAY_DIR, "no_such_file.rds"); "returned a path" },
             error = function(e) sub("\n.*", "", conditionMessage(e))), "\n", sep = "")
cat("  summary requirement waived:   ", !isTRUE(bayrc.needs.summary), "\n", sep = "")

cat("\n", length(vars) - length(missing), " of ", length(vars),
    " paths present\n", sep = "")
if (length(missing))
  cat("set these before running a case study: ", paste(missing, collapse = ", "),
      "\n", sep = "")
