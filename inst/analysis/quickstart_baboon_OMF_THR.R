# End-to-end BayRC Quick Start script: baboon omental fat (OMF) vs. thyroid
# (THR) on the bundled real GSE98965 data (Mure et al. 2018).

library(BayRC)

baboon <- readRDS(system.file("extdata", "baboon_OMF_THR_GSE98965.rds", package = "BayRC"))
n_genes <- nrow(baboon$expr_OMF)

data_list_OMF <- list(data = as.data.frame(log2(baboon$expr_OMF + 1)),
                      time = baboon$zt, gname = baboon$gene_symbol)
data_list_THR <- list(data = as.data.frame(log2(baboon$expr_THR + 1)),
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
mcmc_OMF <- run_mcmc(data_list_OMF, seed = 1)
mcmc_THR <- run_mcmc(data_list_THR, seed = 1)

mcmc_OMF <- match_symbols(mcmc_OMF, BF = 3, p_rhythmic = 0.2)
mcmc_THR <- match_symbols(mcmc_THR, BF = 3, p_rhythmic = 0.2)

bf_OMF <- summarize_bay(mcmc_OMF$rho, BF = 3, p_rhythmic = 0.2)
head(bf_OMF[order(-bf_OMF$BayesF), c("RowAverage", "BayesF")], 5)

detected <- detect_rhy(mcmc_OMF, mcmc_THR, bfdr_alpha = 0.20)

pA <- rowMeans(mcmc_OMF$rho)
pB <- rowMeans(mcmc_THR$rho)
trans <- transition_classify(pA, pB, bfdr_alpha = 0.20)
cat("tau_gain =", trans$tau_gain, " n_gain =", trans$n_gain,
    " | tau_loss =", trans$tau_loss, " n_loss =", trans$n_loss,
    " | n_cons =", trans$n_cons, "\n")

phase <- phase_infer(phi_matrix1 = mcmc_OMF$phi, phi_matrix2 = mcmc_THR$phi,
                     gain_loss_status = trans$gain_loss_status,
                     shift = 2, P = 24, bfdr_alpha = 0.20, compute_hdi = TRUE)

kegg <- readRDS(system.file("extdata", "kegg_pathway_list_hsa.rds", package = "BayRC"))

result_union <- pathSelect(mcmc.merge.list = list(A = mcmc_OMF, B = mcmc_THR),
                           pathway.list = kegg, dataset.names = c("A", "B"),
                           ranking.method = "union", score_type = "pos",
                           qvalue.cut = 0.20, nperm = 500)
active <- dplyr::filter(result_union$results, pval < 0.05)$pathway
cat("Active pathways (Stage 1, pval<0.05):", length(active), "of",
    nrow(result_union$results), "\n")

global <- multi_conservation(mcmc.merge.list = list(A = mcmc_OMF, B = mcmc_THR),
                             dataset.names = c("A", "B"),
                             select.pathway.list = "global",
                             n_perm = 200, n_boot = 200, use_cpp = TRUE,
                             save_output = FALSE)
print(global)
