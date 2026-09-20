## Test whether any pathway's phase difference departs from the genome-wide
## offset, using the maintained genes outside the pathway as the null.
##
## Usage: Rscript sunput_module_offset_test.R [shift]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
suppressPackageStartupMessages(library(BayRC))

args  <- commandArgs(trailingOnly = TRUE)
shift <- if (length(args) >= 1) as.numeric(args[1]) else 2

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))
rho <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
phi <- get(grep("^mcmc_phi",  ls(), value = TRUE)[1])
kegg <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))

tr <- transition_classify(rowMeans(rho[["PUT"]]), rowMeans(rho[["SUN"]]),
                          bfdr_alpha = 0.25)
ph <- phase_infer(phi_matrix1 = phi[["PUT"]], phi_matrix2 = phi[["SUN"]],
                  gain_loss_status = tr$gain_loss_status,
                  bfdr_alpha = 0.25, shift = shift, P = 24, compute_hdi = FALSE)

d <- (ph$peak2 - ph$peak1) %% 24
d <- ifelse(d > 12, d - 24, d)
genes <- rownames(rho[["SUN"]])
maint <- tr$gain_loss_status == "Maintained" & is.finite(d)
cat(sprintf("maintained genes: %d   global offset %+.2f h\n\n", sum(maint), mean(d[maint])))

tested <- names(kegg)[sapply(kegg, function(p) sum(genes %in% p) >= 15)]
res <- do.call(rbind, lapply(tested, function(p) {
  inp <- maint & genes %in% kegg[[p]]
  if (sum(inp) < 5) return(NULL)
  out <- maint & !(genes %in% kegg[[p]])
  tt  <- t.test(d[inp], d[out])
  data.frame(pathway = sub("^KEGG ", "", p), n = sum(inp),
             mean_in = mean(d[inp]), mean_out = mean(d[out]),
             diff = mean(d[inp]) - mean(d[out]), p = tt$p.value)
}))
res$q <- p.adjust(res$p, "BH")
res <- res[order(res$p), ]
write.csv(res, file.path(BAYRC_OUTPUT_DIR, "module_offset_test.csv"), row.names = FALSE)

cat(sprintf("pathways tested: %d   any q < 0.05: %d   any q < 0.20: %d\n\n",
            nrow(res), sum(res$q < 0.05), sum(res$q < 0.20)))
cat("ten smallest p-values:\n")
print(head(res, 10), row.names = FALSE, digits = 3)
