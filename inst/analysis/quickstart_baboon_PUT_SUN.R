# End-to-end BayRC Quick Start script: baboon putamen (PUT) vs. substantia
# nigra (SUN), mirroring the paper's Case Study 2 on the bundled real
# GSE98965 data (Mure et al. 2018). Produces the exact numbers documented
# in README.md's Quick Start section. MCMC takes about 30-40 minutes per
# condition on the full 5,066-gene set.

library(BayRC)

# ── STEP 1: Load the bundled data and run MCMC ────────────────────────────
baboon <- readRDS(system.file("extdata", "baboon_PUT_SUN_GSE98965.rds", package = "BayRC"))
n_genes <- nrow(baboon$expr_PUT)

data_list_PUT <- list(data = as.data.frame(log2(baboon$expr_PUT + 1)),
                      time = baboon$zt, gname = baboon$gene_symbol)
data_list_SUN <- list(data = as.data.frame(log2(baboon$expr_SUN + 1)),
                      time = baboon$zt, gname = baboon$gene_symbol)

run_mcmc <- function(dat, seed) {
  init <- CBt_init_single(Data.list = dat, P = 24, FitCosinor = TRUE,
                          mu_M = 0, sigma_M = 10, mu_A = 1, sigma_A = 10, seed = seed)
  CB_MCMC_single_rj_slice(
    Data.list = dat, Init.value = init, P = 24,
    iteration = 2500, thin = 1, n.burn = 500, seed = seed,
    p_rhythmic = rep(0.2, n_genes), rj.p.stay = 0.5,
    A_prior = "trunc_Normal_OLS_condi", mu_A = 1, sigma_A = 10^2, A.min = 0,
    A_wb_beta2 = 2, A_gm_shape = 1.99, A_gm_rate = 0.5,
    rj.phi = TRUE, rj.A = TRUE, mu_M = 0, sigma_M = 10^2,
    sigma_prior_v = 2, sigma_prior_s = 0
  )
}
mcmc_PUT <- run_mcmc(data_list_PUT, seed = 1)
mcmc_SUN <- run_mcmc(data_list_SUN, seed = 1)

# ── STEP 2: Annotate (must run before anything downstream) ───────────────
mcmc_PUT <- match_symbols(mcmc_PUT, BF = 3, p_rhythmic = 0.2)
mcmc_SUN <- match_symbols(mcmc_SUN, BF = 3, p_rhythmic = 0.2)

est_PUT <- CB_getAllEst(mcmc_PUT, burn = 50)
names(est_PUT) <- c("A", "phi", "M", "sigma")

# ═══ PART 1: single-group biomarker detection ════════════════════════════
bf_PUT <- summarize_bay(mcmc_PUT$rho, BF = 3, p_rhythmic = 0.2)
print(head(bf_PUT[order(-bf_PUT$BayesF), c("RowAverage", "BayesF")], 5))

detected <- detect_rhy(mcmc_PUT, mcmc_SUN, bfdr_alpha = 0.25)
cat("PUT rhythmic:", detected$n_rhythmic_A, "/", detected$n_total, "\n")
cat("SUN rhythmic:", detected$n_rhythmic_B, "/", detected$n_total, "\n")

# ═══ PART 2: two-group comparison ═════════════════════════════════════════
pA <- rowMeans(mcmc_PUT$rho)
pB <- rowMeans(mcmc_SUN$rho)

trans <- transition_classify(pA, pB, bfdr_alpha = 0.25)

phase <- phase_infer(phi_matrix1 = mcmc_PUT$phi, phi_matrix2 = mcmc_SUN$phi,
                     gain_loss_status = trans$gain_loss_status,
                     shift = 2, P = 24, bfdr_alpha = 0.25, compute_hdi = TRUE)
cat("Phase-shifted:", sum(phase$flag_shift, na.rm = TRUE),
    " Phase-conserved:", sum(phase$flag_cons, na.rm = TRUE),
    " Undetermined:", sum(phase$flag_undetermined, na.rm = TRUE), "\n")

# ── STEP 3: Two-stage pathway enrichment ──────────────────────────────────
kegg <- readRDS(system.file("extdata", "kegg_pathway_list_hsa.rds", package = "BayRC"))

result_union <- pathSelect(mcmc.merge.list = list(A = mcmc_PUT, B = mcmc_SUN),
                           pathway.list = kegg, dataset.names = c("A", "B"),
                           ranking.method = "union", score_type = "pos",
                           qvalue.cut = 0.20, nperm = 500)
active <- dplyr::filter(result_union$results, pval < 0.05)$pathway
cat("Active pathways (Stage 1, pval<0.05):", length(active), "of",
    nrow(result_union$results), "\n")

# ── STEP 4: Genome-wide concordance ───────────────────────────────────────
global <- multi_conservation(mcmc.merge.list = list(A = mcmc_PUT, B = mcmc_SUN),
                             dataset.names = c("A", "B"),
                             select.pathway.list = "global",
                             n_perm = 200, n_boot = 200, use_cpp = TRUE,
                             save_output = FALSE)
print(global)
