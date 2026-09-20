## Per-pathway phase composition for one tissue pair: how many maintained genes
## a pathway carries and how they split between phase-shifted and conserved.
##
## Usage: Rscript pair_phase_composition.R tissueA tissueB [bfdr_alpha] [shift]

this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))
suppressPackageStartupMessages(library(BayRC))

args  <- commandArgs(trailingOnly = TRUE)
tA    <- args[1]; tB <- args[2]
alpha <- if (length(args) >= 3) as.numeric(args[3]) else 0.25
shf   <- if (length(args) >= 4) as.numeric(args[4]) else 2

load(file.path(BAYRC_SUMMARY_DIR, "mcmc_rho_BF3.RData"))
load(file.path(BAYRC_SUMMARY_DIR, "phi", "mcmc_phi_BF3.RData"))
rho  <- get(grep("^mcmc_data", ls(), value = TRUE)[1])
phi  <- get(grep("^mcmc_phi",  ls(), value = TRUE)[1])
kegg <- readRDS(file.path(BAYRC_PATHWAY_DIR, "kegg_pathway_list_hsa.rds"))

tr <- transition_classify(rowMeans(rho[[tA]]), rowMeans(rho[[tB]]), bfdr_alpha = alpha)
ph <- phase_infer(phi_matrix1 = phi[[tA]], phi_matrix2 = phi[[tB]],
                  gain_loss_status = tr$gain_loss_status,
                  bfdr_alpha = alpha, shift = shf, P = 24, compute_hdi = FALSE)

genes <- rownames(rho[[tA]])
st    <- tr$gain_loss_status
d     <- (ph$peak2 - ph$peak1) %% 24; d <- ifelse(d > 12, d - 24, d)
maint <- st == "Maintained" & is.finite(d)

cat(sprintf("\n%s vs %s   alpha %.2f  shift %g\n", tA, tB, alpha, shf))
cat(sprintf("gain %d  loss %d  maintained %d   dominant: %s\n",
            sum(st == "Gain"), sum(st == "Loss"), sum(maint),
            c("Gain", "Loss", "Maintained")[which.max(c(sum(st == "Gain"),
              sum(st == "Loss"), sum(maint)))]))
cat(sprintf("global offset %+.2f h  SD %.2f  %.1f%% share sign\n",
            mean(d[maint]), sd(d[maint]), 100 * max(mean(d[maint] > 0), mean(d[maint] < 0))))

res <- do.call(rbind, lapply(names(kegg), function(p) {
  i <- which(genes %in% kegg[[p]]); if (length(i) < 15) return(NULL)
  m <- intersect(i, which(maint)); if (!length(m)) return(NULL)
  data.frame(pathway = substr(sub("^KEGG ", "", p), 1, 36),
             meas = length(i), maint = length(m),
             shift = sum(ph$flag_shift[i], na.rm = TRUE),
             cons  = sum(ph$flag_cons[i],  na.rm = TRUE),
             loss  = sum(st[i] == "Loss"), gain = sum(st[i] == "Gain"),
             mean_d = round(mean(d[m]), 2))
}))
res$balance <- pmin(res$shift, res$cons)
cat(sprintf("\nspread of pathway mean phase difference: SD %.2f, range %+.2f to %+.2f\n",
            sd(res$mean_d), min(res$mean_d), max(res$mean_d)))
cat("\ntop 12 pathways by min(shifted, conserved):\n")
print(head(res[order(-res$balance, -res$maint), ], 12), row.names = FALSE)
