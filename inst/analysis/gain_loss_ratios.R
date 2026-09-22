## Gain-loss ratios as the Methods define them, from each tissue's marginal
## posterior probability of rhythmicity: E[Gain] / E[Loss], with
## E[Gain] = sum (1 - pA) pB and E[Loss] = sum pA (1 - pB).
##
## Prints the genome-wide ratio for putamen against substantia nigra and against
## visual cortex, the ratio for every pathway in Table 1 from the expected counts
## in stage2_significant.csv, and the Spearman correlation between the ratio and
## the ratio of rhythmic content, sum pB / sum pA, over the 325 baboon tissue
## pairs. Writes gain_loss_ratios.csv under BAYRC_OUTPUT_DIR.

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
source(file.path(analysis.dir, "config.R"))

hb <- new.env()
load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"), envir = hb)
p <- sapply(hb$mcmc_data_baboon, rowMeans)

glr <- function(a, b) sum((1 - a) * b) / sum(a * (1 - b))

rows <- list()
for (cmp in c("SUN", "VIC")) {
  rows[[length(rows) + 1]] <- data.frame(pair = paste0("PUT_", cmp), pathway = "genome-wide",
                                         size = nrow(p), GLR = glr(p[, "PUT"], p[, cmp]))
  s <- read.csv(file.path(BAYRC_OUTPUT_DIR, "figures", paste0("baboon_PUT_", cmp),
                          "stage2_significant.csv"))
  s <- unique(s[, c("pathway", "size", "Expected_N_Gain", "Expected_N_Loss")])
  rows[[length(rows) + 1]] <- data.frame(pair = paste0("PUT_", cmp), pathway = s$pathway,
                                         size = s$size,
                                         GLR = s$Expected_N_Gain / s$Expected_N_Loss)
}
out <- do.call(rbind, rows)
out$GLR <- round(out$GLR, 3)
print(out, row.names = FALSE)
write.csv(out, file.path(BAYRC_OUTPUT_DIR, "gain_loss_ratios.csv"), row.names = FALSE)

pairs   <- combn(colnames(p), 2)
ratio   <- apply(pairs, 2, function(k) glr(p[, k[1]], p[, k[2]]))
content <- apply(pairs, 2, function(k) sum(p[, k[2]]) / sum(p[, k[1]]))
cat(sprintf("\n%d tissue pairs: Spearman correlation of GLR with rhythmic content ratio = %.3f\n",
            ncol(pairs), cor(ratio, content, method = "spearman")))
