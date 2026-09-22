#' Run MCMC for a single circadian dataset with Zeitgeber-time uncertainty
#'
#' @title CBt MCMC single chain with time-error model
#'
#' @noRd
#'
#' @description
#' Extends \code{CB_MCMC_single_rj_slice} to jointly infer the gene
#' expression parameters and the latent Zeitgeber-time measurement error
#' (t_p) for each sample.  The time-error t_p is modelled with a
#' von Mises or spike-and-slab prior and updated via Metropolis-Hastings
#' within each iteration.
#'
#' @param Data.list A named list with three elements: \code{data} (G x N
#'   data.frame), \code{time} (length-N Zeitgeber time in hours), and
#'   \code{gname} (length-G gene identifiers).
#' @param Init.value List returned by \code{CBt_init_single} containing
#'   starting values for \code{rho}, \code{M}, \code{A}, \code{phi},
#'   \code{sigma}, and \code{t_p}.
#' @param P Numeric; period in hours (default 24).
#' @param iteration Integer; total iterations including burn-in (default 1000).
#' @param thin Integer; thinning interval (default 20).
#' @param n.burn Integer; burn-in iterations (default 1000).
#' @param seed Integer; random seed (default 15213).
#' @param p_rhythmic Numeric vector length G; prior Pr(rho = 1) per gene.
#' @param rj.p.stay Numeric; probability of skipping between-model move.
#' @param t_p_gene Character; which genes inform the t_p update: one of
#'   \code{"rhythmic"}, \code{"fixed"}, or \code{"p_rhythmic"}.
#' @param fix_rhythm Length-G binary vector used when
#'   \code{t_p_gene = "fixed"}.
#' @param p_rhythm_cut Numeric threshold for \code{t_p_gene = "p_rhythmic"}.
#' @param t_p_prior Character; prior on t_p: \code{"VM"} (von Mises),
#'   \code{"mass0_spike"}, \code{"cont_spike"}, or \code{"cnot0_spike"}.
#' @param MH_kappa_t_p Numeric; Metropolis-Hastings step size for t_p.
#' @param theta_t_p,kappa_t_p,kappa_t0_p,alpha_t1,alpha_t2 Numeric;
#'   hyperparameters for the t_p prior.
#' @param A_prior Character; amplitude prior name (see
#'   \code{CB_MCMC_single_rj_slice}).
#' @param mu_A,sigma_A,A.min Numeric; truncated-Normal amplitude prior.
#' @param A_wb_beta2,A_gm_shape,A_gm_rate Numeric; alternative amplitude
#'   prior parameters.
#' @param rj.phi,rj.A Logical; jointly propose phi/A in RJMCMC jump.
#' @param mu_M,sigma_M Numeric; MESOR prior mean and variance.
#' @param sigma_prior_v,sigma_prior_s Numeric; inverse-gamma prior on
#'   residual variance.
#' @param save.file,save.file2,save.file3 Character; paths for intermediate saves.
#'
#' @return Same structure as \code{CB_MCMC_single_rj_slice} plus a
#'   \code{t_p} matrix (N x K) of posterior samples of time-error.
#'
#'
#' @examples
#' \dontrun{
#' sim  <- CBt_sim_data()
#' dat  <- list(data  = as.data.frame(sim$data[[1]]$dat),
#'              time  = sim$data[[1]]$x$time,
#'              gname = paste0("G", seq_len(nrow(sim$data[[1]]$dat))))
#' init <- CBt_init_single(dat)
#' mcmc_out <- CBt_MCMC_single(dat, init,
#'               iteration = 200, thin = 10, n.burn = 100,
#'               theta_t_p = 0, kappa_t_p = 0.5, kappa_t0_p = 0.1,
#'               alpha_t1 = 1, alpha_t2 = 1)
#' }
CBt_MCMC_single = function(Data.list, Init.value, P = 24,
                   iteration = 1000, thin = 20, n.burn=1000, 
                   seed = 15213, 
                   p_rhythmic=rep(0.2, 1000), rj.p.stay = 0.5,
                   t_p_gene = "rhythmic",
                   # "rhythmic": genes with rho = 1 in the current iteration
                   # "fixed": the 0/1 vector fix_rhythm
                   # "p_rhythmic": genes with posterior Pr(rho = 1) above p_rhythm_cut
                   fix_rhythm = c(rep(1, 10), rep(0, 990)),
                   p_rhythm_cut = 0.7,
                   t_p_prior = "VM",
                   # "VM": von Mises prior
                   # "mass0_spike": spike and slab with a point mass at 0
                   # "cont_spike": continuous spike
                   # "cnot0_spike": continuous spike with y dependent on zeta
                   MH_kappa_t_p = 0.5, theta_t_p, kappa_t_p, kappa_t0_p,
                   alpha_t1, alpha_t2,
                   A_prior= "Jeffreys_OLS_condi",
                   mu_A = 1, sigma_A = 10^2, A.min = 0, # truncated-Normal prior on A
                   A_wb_beta2=2, # A ~ sq_expo(2, beta)
                   A_gm_shape=1.99, # gamma shape <= 2; flatter near 2
                   A_gm_rate = 0.5, # smaller is flatter
                   rj.phi = TRUE, rj.A = TRUE,
                   mu_M = 0, sigma_M = 10^2, 
                   sigma_prior_v = 4, sigma_prior_s = 1,
                   save.file = "MCMC_save.rds", 
                   save.file2 = "MCMC_save.rds",
                   save.file3 = "MCMC_save.rds"
){
  omega = 2*pi/P
  Y = as.matrix(Data.list[[1]])
  A.max = apply(Y, 1, function(x){(max(x)-min(x))/2})
  t = Data.list[[2]]
  G = nrow(Y)
  N = ncol(Y)
  t.c = cos(omega*t)
  t.s = sin(omega*t)
  t.c.sum = sum(t.c)
  t.s.sum = sum(t.s)
  c.t.phi = sum(t.c^2)
  s.t.phi = sum(t.s^2)
  cs.t.phi = sum(t.c*t.s)
  y.t.c.sum = apply(t(Y)*t.c, 2, sum)
  y.t.s.sum = apply(t(Y)*t.s, 2, sum)
  # design and OLS projection for the slice sampler
  X = cbind(t.c, t.s)
  XtX = t(X)%*%X
  beta_hat0 =  t(solve(XtX)%*%t(X))
  beta_cov0 = solve(XtX) # scaled by sigma in each iteration
  beta_cov0_11 = beta_cov0[1, 1]
  beta_cov0_22 = beta_cov0[2, 2]
  beta_cov0_12 = beta_cov0[1, 2]
  beta_cov0_rho = beta_cov0_12/sqrt(beta_cov0_11*beta_cov0_22)
  # conditional-mean coefficients of one coordinate on the other
  beta_mean0_1 = beta_cov0_12/beta_cov0_22
  beta_mean0_2 = beta_cov0_12/beta_cov0_11

  # Initialize the chain  ---------------------------------------------------
  rho = Init.value$rho;
  M = Init.value$M;
  A = Init.value$A;
  phi = Init.value$phi;
  sigma = Init.value$sigma; 
  AcosPhi = A*cos(omega*phi);
  AsinPhi = A*sin(omega*phi);
  t_p = Init.value$t_p
  zeta = rep(0, length(t_p))
  t_t_p = t+t_p
  t_t_p.c = cos(omega*(t_t_p))
  t_t_p.s = sin(omega*(t_t_p))
  t_t_p.c.sum = sum(t_t_p.c)
  t_t_p.s.sum = sum(t_t_p.s)
  c.t_t_p.phi = sum(t_t_p.c^2)
  s.t_t_p.phi = sum(t_t_p.s^2)
  cs.t_t_p.phi = sum(t_t_p.c*t_t_p.s)
  y.t_t_p.c.sum = apply(t(Y)*t_t_p.c, 2, sum)
  y.t_t_p.s.sum = apply(t(Y)*t_t_p.s, 2, sum)
  # the same quantities on the shifted times t + t_p
  X_t_p = cbind(t_t_p.c, t_t_p.s)
  XtX_t_p = t(X_t_p)%*%X_t_p
  beta_hat0_t_p =  t(solve(XtX_t_p)%*%t(X_t_p))
  beta_cov0_t_p = solve(XtX_t_p)
  beta_cov0_11_t_p = beta_cov0_t_p[1, 1]
  beta_cov0_22_t_p = beta_cov0_t_p[2, 2]
  beta_cov0_12_t_p = beta_cov0_t_p[1, 2]
  beta_cov0_rho_t_p = beta_cov0_12_t_p/sqrt(beta_cov0_11_t_p*beta_cov0_22_t_p)
  beta_mean0_1_t_p = beta_cov0_12_t_p/beta_cov0_22_t_p
  beta_mean0_2_t_p = beta_cov0_12_t_p/beta_cov0_11_t_p
  rho.store = rep(1, G)

  # Burn-in -----------------------------------------------------------------
  for(iter in 1:n.burn){
    
    set.seed(seed+iter)
    print(paste0("Burn iter = ", iter))
    
    # update btw model --------------------------------------------------------
    stay = stats::rbinom(1, 1, rj.p.stay)
    if(!stay){
      rho.res = try_save(RJMCMC_single_slice(Y, t_t_p.c, t_t_p.s, N,
                                             t_t_p.c.sum, t_t_p.s.sum,
                                             c.t_t_p.phi, s.t_t_p.phi, cs.t_t_p.phi,
                                             y.t_t_p.c.sum, y.t_t_p.s.sum,
                                             AcosPhi, AsinPhi, A, phi,
                                             M, sigma, p_rhythmic, rho,
                                             rj.phi, rj.A,
                                             A_prior,
                                             mu_A, sigma_A, A.min, A.max,
                                             A_wb_beta2,
                                             A_gm_shape, A_gm_rate,
                                             omega, G, P, save.file3),
                         save, save.file)
      if.accept.rj = rho.res$rho!=rho
      rho = rho.res$rho
      rho.store = cbind(rho.store, rho)
    }
    time1 = Sys.time()

    # within model move -------------------------------------------------------

    M = try_save(update_M_single(Y, t_t_p.c, t_t_p.s, N,
                                 AcosPhi, AsinPhi, sigma, rho,
                                 sigma_M, mu_M, omega, G),
                 save, save.file)

    set.seed(seed+iter+1001)
    csAphi = try_save(update_A_phi_slice(Y, X_t_p, XtX_t_p, beta_hat0_t_p,
                                         beta_cov0_11_t_p, beta_cov0_22_t_p, beta_cov0_rho_t_p,
                                         beta_mean0_1_t_p, beta_mean0_2_t_p,
                                         M, sigma, omega,
                                         AcosPhi, AsinPhi, G, P, A_prior,
                                         A.min, A.max, 
                                         mu_A, sigma_A,
                                         A_wb_beta2,
                                         A_gm_shape, A_gm_rate),
                      save, save.file)
    AcosPhi = csAphi$AcosPhi
    AsinPhi = csAphi$AsinPhi
    A = sqrt(AcosPhi^2+AsinPhi^2)
    phi = sapply(atan2(AsinPhi, AcosPhi), adjust.to.2pi)/omega

    sigma = try_save(update_sigma_single(Y, t_t_p.c, t_t_p.s, N, AcosPhi, AsinPhi, M, rho,
                                         sigma_prior_v, sigma_prior_s, omega, G),
                     save, save.file)

    set.seed(seed+iter+1002)
      t_p_res = try_save(
        update_t_p_single_TestPriors(Y, t, t.c, t.s,
                                     t_p, t_t_p.c, t_t_p.s, zeta, pt,    
                                     rho, A, phi, AcosPhi, AsinPhi, M, sigma, 
                                     theta_t_p, kappa_t_p, kappa_t0_p, 
                                     MH_kappa_t_p, alpha_t1, alpha_t2, 
                                     prior = t_p_prior, t_p_gene,
                                     fix_rhythm, p_rhythm_cut, rho.store, 
                                     N, G, omega, P, save.file2),
                     save, save.file)
    if.accept.tp = t_p_res$if.accept
    zeta = t_p_res$zeta
    pt = t_p_res$pt
    t_p = t_p_res$t_p
    R1 = t_p_res$R1
    R2 = t_p_res$R2
    alpha1 = t_p_res$alpha1
    alpha2 = t_p_res$alpha2
    t_t_p = t+t_p
    t_t_p.c = cos(omega*(t_t_p))
    t_t_p.s = sin(omega*(t_t_p))
    t_t_p.c.sum = sum(t_t_p.c)
    t_t_p.s.sum = sum(t_t_p.s)
    c.t_t_p.phi = sum(t_t_p.c^2)
    s.t_t_p.phi = sum(t_t_p.s^2)
    cs.t_t_p.phi = sum(t_t_p.c*t_t_p.s)
    y.t_t_p.c.sum = apply(t(Y)*t_t_p.c, 2, sum)
    y.t_t_p.s.sum = apply(t(Y)*t_t_p.s, 2, sum)
    # recompute the design on the updated t + t_p
    X_t_p = cbind(t_t_p.c, t_t_p.s)
    XtX_t_p = t(X_t_p)%*%X_t_p
    beta_hat0_t_p =  t(solve(XtX_t_p)%*%t(X_t_p))
    beta_cov0_t_p = solve(XtX_t_p)
    beta_cov0_11_t_p = beta_cov0_t_p[1, 1]
    beta_cov0_22_t_p = beta_cov0_t_p[2, 2]
    beta_cov0_12_t_p = beta_cov0_t_p[1, 2]
    beta_cov0_rho_t_p = beta_cov0_12_t_p/sqrt(beta_cov0_11_t_p*beta_cov0_22_t_p)
    beta_mean0_1_t_p = beta_cov0_12_t_p/beta_cov0_22_t_p
    beta_mean0_2_t_p = beta_cov0_12_t_p/beta_cov0_11_t_p

  }

  rho.store = rho
  M.store = M;
  sigma.store = sigma;
  AcosPhi.store = AcosPhi; 
  AsinPhi.store = AsinPhi; 
  A.store = A
  phi.store = phi
  log.r1.store = rho.res$log.r1
  log.r1.SS.store = rho.res$log.r1.SS
  log.r3_A.store = rho.res$log.r3_A
  log.r3_phi.store = rho.res$log.r3_phi
  if.accept.rj.store = if.accept.rj
  if.accept.tp.store = if.accept.tp
  t_p.store = t_p
  zeta.store = zeta
  pt.store = pt
  alpha1.store = alpha1
  alpha2.store = alpha2
  R1.store = R1
  R2.store = R2

  # Start sampling ----------------------------------------------------------

  for(iter in 1:(iteration-n.burn)){

    set.seed(seed+iter)
    print(paste0("iter = ", iter))

    # update btw model --------------------------------------------------------
    stay = stats::rbinom(1, 1, rj.p.stay)
    if(stay){
      rho.store = cbind(rho.store, rho)
      if.accept.rj = cbind(if.accept.rj, -1) 
      log.r1.store = cbind(log.r1.store, -99)
      log.r1.SS.store = cbind(log.r1.SS.store, -99) 
    }else{
      rho.res = try_save(RJMCMC_single_slice(Y, t_t_p.c, t_t_p.s, N,
                                             t_t_p.c.sum, t_t_p.s.sum,
                                             c.t_t_p.phi, s.t_t_p.phi, cs.t_t_p.phi,
                                             y.t_t_p.c.sum, y.t_t_p.s.sum,
                                             AcosPhi, AsinPhi, A, phi,
                                             M, sigma, p_rhythmic, rho,
                                             rj.phi, rj.A,
                                             A_prior,
                                             mu_A, sigma_A, A.min, A.max, 
                                             A_wb_beta2,
                                             A_gm_shape, A_gm_rate,
                                             omega, G, P, save.file3),
                         save, save.file)
      if.accept.rj.store = cbind(if.accept.rj.store, rho.res$rho!=rho) # a change in rho is an accepted jump
      log.r1.store = cbind(log.r1.store, rho.res$log.r1)
      log.r1.SS.store = cbind(log.r1.SS.store, rho.res$log.r1.SS) 
      log.r3_A.store = cbind(log.r3_A.store, rho.res$log.r3_A)
      log.r3_phi.store = cbind(log.r3_phi.store, rho.res$log.r3_phi)
      rho = rho.res$rho
      rho.store = cbind(rho.store, rho)
    }

    # within model move -------------------------------------------------------

    M = try_save(update_M_single(Y, t_t_p.c, t_t_p.s, N,
                                 AcosPhi, AsinPhi, sigma, rho,
                                 sigma_M, mu_M, omega, G),
                 save, save.file)
    M.store = cbind(M.store, M)

    set.seed(seed+iter+1001)
    
    csAphi = try_save(update_A_phi_slice(Y, X_t_p, XtX_t_p, beta_hat0_t_p,
                                         beta_cov0_11_t_p, beta_cov0_22_t_p, beta_cov0_rho_t_p,
                                         beta_mean0_1_t_p, beta_mean0_2_t_p,
                                         M, sigma, omega,
                                         AcosPhi, AsinPhi, G, P, A_prior,
                                         A.min, A.max, 
                                         mu_A, sigma_A,
                                         A_wb_beta2,
                                         A_gm_shape, A_gm_rate),
                      save, save.file)

    if(sum(is.na(csAphi$AcosPhi))>0|sum(is.na(csAphi$AsinPhi))>0){
      save(Y, X_t_p, XtX_t_p, beta_hat0_t_p,
           beta_cov0_11_t_p, beta_cov0_22_t_p, beta_cov0_rho_t_p,
           beta_mean0_1_t_p, beta_mean0_2_t_p,
           M, sigma, omega,
           AcosPhi, AsinPhi, G, P, A_prior,
           mu_A, sigma_A,
           A_wb_beta2,
           A_gm_shape, A_gm_rate, file = paste0(save.file3, "1"))
      stop("Aphi error")
    }
    AcosPhi = csAphi$AcosPhi
    AsinPhi = csAphi$AsinPhi
    A = sqrt(AcosPhi^2+AsinPhi^2)
    phi = sapply(atan2(AsinPhi, AcosPhi), adjust.to.2pi)/omega

    AcosPhi.store = cbind(AcosPhi.store, AcosPhi); 
    AsinPhi.store = cbind(AsinPhi.store, AsinPhi); 
    A.store = cbind(A.store, A)
    phi.store = cbind(phi.store, phi)
    
    sigma = try_save(update_sigma_single(Y, t_t_p.c, t_t_p.s, N, AcosPhi, AsinPhi, M, rho,
                                         sigma_prior_v, sigma_prior_s, omega, G),
                     save, save.file)
    sigma.store = cbind(sigma.store, sigma)

    set.seed(seed+iter+1002)
    
    t_p_res = try_save(
      update_t_p_single_TestPriors(Y, t, t.c, t.s,
                                   t_p, t_t_p.c, t_t_p.s, zeta, pt,    
                                   rho, A, phi, AcosPhi, AsinPhi, M, sigma, 
                                   theta_t_p, kappa_t_p, kappa_t0_p, 
                                   MH_kappa_t_p, alpha_t1, alpha_t2, 
                                   prior = t_p_prior, t_p_gene,
                                   fix_rhythm, p_rhythm_cut, rho.store, 
                                   N, G, omega, P, save.file2),
      save, save.file)

    if.accept.tp = t_p_res$if.accept
    zeta = t_p_res$zeta
    pt = t_p_res$pt
    t_p = t_p_res$t_p
    R1 = t_p_res$R1
    R2 = t_p_res$R2
    alpha1 = t_p_res$alpha1
    alpha2 = t_p_res$alpha2
    zeta.store = cbind(zeta.store, zeta)
    pt.store = cbind(pt.store, pt)
    t_p.store = cbind(t_p.store, t_p)
    R1.store = cbind(R1.store, R1)
    R2.store = cbind(R2.store, R2)
    alpha1.store = cbind(alpha1.store, alpha1)
    alpha2.store = cbind(alpha2.store, alpha2)
    t_t_p.c = cos(omega*(t+t_p))
    t_t_p.s = sin(omega*(t+t_p))
    t_t_p.c.sum = sum(t_t_p.c)
    t_t_p.s.sum = sum(t_t_p.s)
    c.t_t_p.phi = sum(t_t_p.c^2)
    s.t_t_p.phi = sum(t_t_p.s^2)
    cs.t_t_p.phi = sum(t_t_p.c*t_t_p.s)
    y.t_t_p.c.sum = apply(t(Y)*t_t_p.c, 2, sum)
    y.t_t_p.s.sum = apply(t(Y)*t_t_p.s, 2, sum)
    # recompute the design on the updated t + t_p
    X_t_p = cbind(t_t_p.c, t_t_p.s)
    XtX_t_p = t(X_t_p)%*%X_t_p
    beta_hat0_t_p =  t(solve(XtX_t_p)%*%t(X_t_p))
    beta_cov0_t_p = solve(XtX_t_p)
    beta_cov0_11_t_p = beta_cov0_t_p[1, 1]
    beta_cov0_22_t_p = beta_cov0_t_p[2, 2]
    beta_cov0_12_t_p = beta_cov0_t_p[1, 2]
    beta_cov0_rho_t_p = beta_cov0_12_t_p/sqrt(beta_cov0_11_t_p*beta_cov0_22_t_p)
    beta_mean0_1_t_p = beta_cov0_12_t_p/beta_cov0_22_t_p
    beta_mean0_2_t_p = beta_cov0_12_t_p/beta_cov0_11_t_p

    save = list(rho = rho.store,
                M = M.store,
                AcosPhi = AcosPhi.store,
                AsinPhi = AsinPhi.store, 
                A = A.store,
                phi = phi.store,
                sigma = sigma.store,
                t_p = t_p.store,
                log.r1 = log.r1.store,
                log.r3_A = log.r3_A.store,
                log.r3_phi = log.r3_phi.store,
                if.accept.rj = if.accept.rj.store, 
                if.accept.tp = if.accept.tp.store, 
                pt.store = pt.store, 
                zeta.store = zeta.store,
                R1 = R1.store,
                R2 = R2.store,
                alpha1= alpha1.store,
                alpha2= alpha2.store)
  }
  gname <- Data.list[[3]]
  if (!is.null(gname) && length(gname) == G) {
    for (mat_name in names(save)) {
      if (is.matrix(save[[mat_name]]) && nrow(save[[mat_name]]) == G) {
        rownames(save[[mat_name]]) <- gname
      }
    }
  }
  return(save)
}

update_t_p_single_TestPriors = function(Y, t, t.c, t.s,
                                    t_p, t_t_p.c, t_t_p.s, zeta, pt,    
                                    rho, A, phi, AcosPhi, AsinPhi, M, sigma, 
                                    theta_t_p, kappa_t_p, kappa_t0_p, 
                                    MH_kappa_t_p, alpha_t1, alpha_t2, 
                                    prior, t_p_gene,
                                    fix_rhythm, p_rhythm_cut, rho.store, 
                                    N, G, omega, P, save.file2){
  # dispatches the t_p update to the sampler for the chosen prior

  # coefficients of the omega and 2*omega terms in the t_p full conditional
  x1 = cbind(AcosPhi, AsinPhi) %*% rbind(t.c, t.s)
  x2 = cbind(AcosPhi, AsinPhi) %*% rbind(-1*t.s, t.c)
  y = Y-M
  # component 1, frequency omega
  c1.cos = diag(t(x1*rho/sigma)%*%y)+kappa_t_p*cos(omega*theta_t_p)
  c1.sin = diag(t(x2*rho/sigma)%*%y)+kappa_t_p*sin(omega*theta_t_p)
  # component 2, frequency 2*omega
  c2.1 = apply(x1^2*rho/sigma-x2^2*rho/sigma, 2, sum)*(-1/4)
  c2.2 = apply(x1*x2*rho/sigma, 2, sum)/(-2)
  
  R_alpha1 = get_t_p_R_alpha(c1.cos, c1.sin)
  R_alpha2 = get_t_p_R_alpha(c2.1, c2.2)
  R1 = R_alpha1[, 1]; alpha1 = R_alpha1[, 2]
  R2 = R_alpha2[, 1]; alpha2 = R_alpha2[, 2]
  t_p.save = t_p

  if(t_p_gene=="fixed"){
    rho = fix_rhythm
  }else if(t_p_gene=="p_rhythmic"){
    rho = apply(rho.store, 1, mean)>p_rhythm_cut
  } # otherwise the current rho is used

  if(prior=="VM"){
    # t_p ~ VM(theta_t_p, kappa_t_p), updated by Metropolis-Hastings
    a.res = update_t_p_single_MH(Y, t, t_p, t_t_p.c, t_t_p.s, 
                                 rho, AcosPhi, AsinPhi, M, sigma, 
                                 theta_t_p, kappa_t_p, 
                                 MH_kappa_t_p, 
                                 N, G, omega, P)
    zeta = NULL
    pt = NULL
    t_p = a.res$t_p
    if.accept = a.res$if.accept
    
  }else if(prior=="VM_slice"){
    a.res = update_t_p_single_slice(Y, t.c, t.s, t_p, 
                                    rho, AcosPhi, AsinPhi, M, sigma, 
                                    theta_t_p, kappa_t_p, 
                                    N, G, omega, P, save.file2)
    zeta = NULL
    pt = NULL
    t_p = a.res$t_p
    if.accept = NULL
    
  }else if(prior=="mass0_spike"){
    a.res = update_t_p_mass0_spike(Y, t, t_p, t_t_p.c, t_t_p.s, zeta, 
                                   rho, A, phi, AcosPhi, AsinPhi, M, sigma, 
                                   theta_t_p, kappa_t_p, 
                                   MH_kappa_t_p, alpha_t1, alpha_t2, 
                                   N, G, omega, P)
    zeta = a.res$zeta
    pt = a.res$pt
    t_p = a.res$t_p
    if.accept = a.res$if.accept
    
  }else if(prior=="cont_spike_t_zeta"){
    a.res = update_t_p_cont_spike_t_zeta(Y, t, t_p, t_t_p.c, t_t_p.s, zeta, pt,
                                         rho, AcosPhi, AsinPhi, M, sigma, 
                                         theta_t_p, kappa_t_p, kappa_t0_p,
                                         MH_kappa_t_p, alpha_t1, alpha_t2, 
                                         N, G, omega, P)
    zeta = a.res$zeta
    pt = a.res$pt
    t_p = a.res$t_p
    if.accept = a.res$if.accept
    
  }else if(prior=="cont_spike_zeta_t"){
    a.res = update_t_p_cont_spike_zeta_t(Y, t, t_p, t_t_p.c, t_t_p.s, zeta, pt,
                                 rho, AcosPhi, AsinPhi, M, sigma, 
                                 theta_t_p, kappa_t_p, kappa_t0_p,
                                 MH_kappa_t_p, alpha_t1, alpha_t2, 
                                 N, G, omega, P)
    zeta = a.res$zeta
    pt = a.res$pt
    t_p = a.res$t_p
    if.accept = a.res$if.accept
    
  }else if(prior=="cnot0_spike_t_zeta"){
    a.res = update_t_p_cont0_spike_t_zeta(Y, t, t_p, t_t_p.c, t_t_p.s, zeta, pt,
                                          rho, AcosPhi, AsinPhi, M, sigma, 
                                          theta_t_p, kappa_t_p, kappa_t0_p,
                                          MH_kappa_t_p, alpha_t1, alpha_t2, 
                                          N, G, omega, P)
    zeta = a.res$zeta
    pt = a.res$pt
    t_p = a.res$t_p
    if.accept = a.res$if.accept
    
  }else if(prior=="cont0_spike_zeta_t"){
    a.res = update_t_p_cont0_spike_zeta_t(Y, t, t_p, t_t_p.c, t_t_p.s, zeta, pt,
                                          rho, AcosPhi, AsinPhi, M, sigma, 
                                          theta_t_p, kappa_t_p, kappa_t0_p,
                                          MH_kappa_t_p, alpha_t1, alpha_t2, 
                                          N, G, omega, P)
    zeta = a.res$zeta
    pt = a.res$pt
    t_p = a.res$t_p
    if.accept = a.res$if.accept
  }
  
  if(sum(is.na(t_p))>0|length(t_p)<N){
    case="sum(is.na(t_p))>0|length(t_p)<N"
    save(Y, t, t.c, t.s,
         t_p, t_t_p.c, t_t_p.s, 
         rho, A, phi, AcosPhi, AsinPhi, M, sigma, 
         theta_t_p, kappa_t_p, kappa_t0_p, 
         MH_kappa_t_p, alpha_t1, alpha_t2, 
         prior, t_p_gene,
         fix_rhythm, p_rhythm_cut, rho.store, 
         N, G, omega, P, case, t_p.save, file = save.file2)
    stop("NA_tp")
  }
  
  t_p = adjust.t_p(t_p, omega)
  return(list(t_p = t_p, 
              zeta = zeta,
              pt = pt,
              if.accept = if.accept, 
              R1 = R1, 
              R2 = R2, 
              alpha1 = alpha1, 
              alpha2 = alpha2))
  
}

update_t_p_mass0_spike = function(Y, t, t_p, t_t_p.c, t_t_p.s, zeta, 
                                  rho, A, phi, AcosPhi, AsinPhi, M, sigma, 
                                  theta_t_p, kappa_t_p, 
                                  MH_kappa_t_p, alpha_t1, alpha_t2, 
                                  N, G, omega, P){
  
  t_p_zeta1_0 = update_t_p_single_MH(Y, t, t_p, t_t_p.c, t_t_p.s, 
                                     rho, AcosPhi, AsinPhi, M, sigma, 
                                     theta_t_p, kappa_t_p, 
                                     MH_kappa_t_p, 
                                     N, G, omega, P)
  t_p_zeta1 = t_p_zeta1_0$t_p
  if.accept = t_p_zeta1_0$if.accept
  
  pt = rbeta(1, alpha_t1 + sum(zeta), alpha_t2 + sum(1 - zeta))
  p_zeta1= sapply(1:N, function(i){
    intee = integrate_t_p_i(i, Y, t_p, M, sigma,
                          AcosPhi, AsinPhi, t.c, t.s, rho, 
                          theta_t_p, kappa_t_p, omega, P)
    intee.scale = intee$scale
    intee = intee$intee
    
    l_p_zeta0 = get_t_p_loglik(Y, t_t_p.c, t_t_p.s,
                               rho, AcosPhi, AsinPhi, M, sigma)-intee.scale
    p_zeta1 = intee
    
    odd_zeta1 = pt/2/pi/besselI(kappa_t_p, 0)/(1-pt)*p_zeta1/exp(l_p_zeta0)
    p_zeta1 = odd_zeta1/(odd_zeta1+1)
  })
  new.zeta = rbinom(N, 1, p_zeta1)
  a.t_p = ifelse(new.zeta, t_p_zeta1, 0)
  
  return(list(t_p = a.t_p,
              zeta = new.zeta,
              pt = pt,
              if.accept = if.accept&(new.zeta!=zeta)))
}

update_t_p_cont_spike_t_zeta = function(Y, t, t_p, t_t_p.c, t_t_p.s, zeta, pt,
                                  rho, AcosPhi, AsinPhi, M, sigma, 
                                  theta_t_p, kappa_t_p, kappa_t0_p,
                                  MH_kappa_t_p, alpha_t1, alpha_t2, 
                                  N, G, omega, P){
  
  kappa_t_p_use = ifelse(zeta, kappa_t_p, kappa_t0_p)
  a.res = update_t_p_single_MH(Y, t, t_p, t_t_p.c, t_t_p.s, 
                               rho, AcosPhi, AsinPhi, M, sigma, 
                               theta_t_p, kappa_t_p_use, 
                               MH_kappa_t_p, 
                               N, G, omega, P)
  a.t_p = a.res$t_p
  if.accept = a.res$if.accept
  
  p_zeta1 = get_p_zeta1(a.t_p, theta_t_p, kappa_t_p, kappa_t0_p, pt, omega)
  new.zeta = rbinom(N, 1, p_zeta1)
  
  pt = rbeta(1, alpha_t1 + sum(new.zeta), alpha_t2 + sum(1 - new.zeta))
  

  return(list(t_p = a.t_p,
              zeta = new.zeta,
              pt = pt, 
              if.accept = if.accept))
}

update_t_p_cont_spike_zeta_t = function(Y, t, t_p, t_t_p.c, t_t_p.s, zeta, pt,
                                        rho, AcosPhi, AsinPhi, M, sigma, 
                                        theta_t_p, kappa_t_p, kappa_t0_p,
                                        MH_kappa_t_p, alpha_t1, alpha_t2, 
                                        N, G, omega, P){
  
  
  a.res = update_t_p_single_MH_zeta(Y, t, t.c, t.s, 
                                    t_p, t_t_p.c, t_t_p.s, zeta, pt,
                                    rho, AcosPhi, AsinPhi, M, sigma, 
                                    theta_t_p, kappa_t_p, kappa_t0_p, 
                                    MH_kappa_t_p, 
                                    N, G, omega, P)
  a.t_p = a.res$t_p
  a.zeta = a.res$zeta
  if.accept = a.res$if.accept
  
  pt = rbeta(1, alpha_t1 + sum(a.zeta), alpha_t2 + sum(1 - a.zeta))
  
  return(list(t_p = a.t_p,
              zeta = a.zeta,
              pt = pt, 
              if.accept = if.accept))
}

update_t_p_cont0_spike_t_zeta = function(Y, t, t_p, t_t_p.c, t_t_p.s, zeta, pt,
                                        rho, AcosPhi, AsinPhi, M, sigma, 
                                        theta_t_p, kappa_t_p, kappa_t0_p,
                                        MH_kappa_t_p, alpha_t1, alpha_t2, 
                                        N, G, omega, P){
  
  kappa_t_p_use = ifelse(zeta, kappa_t_p, kappa_t0_p)
  a.res = update_t_p_single_MH(Y, t, t_p*zeta, t_t_p.c, t_t_p.s, 
                               rho, AcosPhi, AsinPhi, M, sigma, 
                               theta_t_p, kappa_t_p_use, 
                               MH_kappa_t_p, 
                               N, G, omega, P)
  a.t_p = a.res$t_p
  if.accept = a.res$if.accept
  
  p_zeta1 = get_p_zeta1_Y(a.t_p, theta_t_p, kappa_t_p, kappa_t0_p, pt, omega)
  new.zeta = rbinom(N, 1, p_zeta1)
  
  pt = rbeta(1, alpha_t1 + sum(new.zeta), alpha_t2 + sum(1 - new.zeta))
  
  return(list(t_p = a.t_p,
              zeta = new.zeta,
              pt = pt, 
              if.accept = if.accept))
}

update_t_p_cont0_spike_zeta_t = function(Y, t, t_p, t_t_p.c, t_t_p.s, zeta, pt,
                                        rho, AcosPhi, AsinPhi, M, sigma, 
                                        theta_t_p, kappa_t_p, kappa_t0_p,
                                        MH_kappa_t_p, alpha_t1, alpha_t2, 
                                        N, G, omega, P){
  
  
  a.res = update_t_p_single_MH_zeta_cont0(Y, t, t.c, t.s, 
                                    t_p, t_t_p.c, t_t_p.s, zeta, pt,
                                    rho, AcosPhi, AsinPhi, M, sigma, 
                                    theta_t_p, kappa_t_p, kappa_t0_p, 
                                    MH_kappa_t_p, 
                                    N, G, omega, P)
  a.t_p = a.res$t_p
  a.zeta = a.res$zeta
  if.accept = a.res$if.accept
  
  pt = rbeta(1, alpha_t1 + sum(a.zeta), alpha_t2 + sum(1 - a.zeta))
  
  return(list(t_p = a.t_p,
              zeta = a.zeta,
              pt = pt, 
              if.accept = if.accept))
}

get_p_zeta1 = function(t_p, theta_t_p, kappa_t_p, kappa_t0_p, pt, omega){
  odds_zeta1 = exp((kappa_t_p-kappa_t0_p)*cos(omega*t_p-theta_t_p))*
    pt/(1-pt)*besselI(kappa_t0_p, 0)/besselI(kappa_t_p, 0)
  p_zeta1 = odd_zeta1/(odd_zeta1+1)
}

get_p_zeta1_Y = function(Y, t.c, t.s, 
                         t_p, t_t_p.c, t_t_p.s, 
                         rho, AcosPhi, AsinPhi, M, sigma, 
                         theta_t_p, kappa_t_p, kappa_t0_p, pt, omega){
  
  loglik_term = get_t_p_loglik(Y, t_t_p.c, t_t_p.s,
                               rho, AcosPhi, AsinPhi, M, sigma)-
    get_t_p_loglik(Y, t.c, t.s,
                   rho, AcosPhi, AsinPhi, M, sigma)
  
  odds_zeta1 = exp(loglik_term+(kappa_t_p-kappa_t0_p)*cos(omega*t_p-theta_t_p))*
    pt/(1-pt)*besselI(kappa_t0_p, 0)/besselI(kappa_t_p, 0)
  p_zeta1 = odd_zeta1/(odd_zeta1+1)
}

get_t_p_loglik = function(Y, t_t_p.c, t_t_p.s,
                          rho, AcosPhi, AsinPhi, M, sigma){
  y = Y-M 
  yhat = cbind(AcosPhi, AsinPhi) %*% rbind(t_t_p.c, t_t_p.s)
  lik = apply((y-yhat)^2/(-2*sigma)*rho, 2, sum)
}

integrate_t_p_i = function(i, Y, t_p, M, sigma,
                         AcosPhi, AsinPhi, t.c, t.s, rho, 
                         theta_t_p, kappa_t_p, omega, P){
  
  x1 = cbind(AcosPhi, AsinPhi) %*% rbind(t.c[i], t.s[i])
  x2 = cbind(AcosPhi, AsinPhi) %*% rbind(-1*t.s[i], t.c[i])
  y = Y[, i]-M
  # component 1, frequency omega
  c1.cos = x1*rho/sigma*y+kappa_t_p*cos(omega*theta_t_p)
  c1.sin = x2*rho/sigma*y+kappa_t_p*sin(omega*theta_t_p)
  # component 2, frequency 2*omega
  c2.1 = sum(x1^2*rho/sigma-x2^2*rho/sigma)*(-1/4)
  c2.2 = sum(x1*x2*rho/sigma)/(-2)
  
  R_alpha1 = get_t_p_R_alpha(c1.cos, c1.sin)
  R_alpha2 = get_t_p_R_alpha(c2.1, c2.2)
  R1 = R_alpha1[, 1]; alpha1 = R_alpha1[, 2]
  R2 = R_alpha2[, 1]; alpha2 = R_alpha2[, 2]
  
  l_ff = function(t){
    R1*cos(omega*t-alpha1)+R2*cos(2*omega*t-alpha2)
  }
  l_ff_out = sapply(seq(from = -P/2, to = P/2, by = 0.1), l_ff)
  # shift by the 80th percentile of the log integrand
  l_ff_out.0.8 = quantile(l_ff_out, probs = 0.8)
  ff = function(t){
    exp(R1*cos(omega*t-alpha1)+R2*cos(2*omega*t-alpha2))-l_ff_out.0.8
  }
  
  intee = integrate(ff, lower = -1*P/2, upper = P/2)
  
  return(list(intee = intee$value, 
              scale = l_ff_out.0.8))
}

update_t_p_single_MH = function(Y, t, t_p, t_t_p.c, t_t_p.s, 
                                rho, AcosPhi, AsinPhi, M, sigma, 
                                theta_t_p, kappa_t_p, 
                                MH_kappa_t_p, 
                                N, G, omega, P){
  
  old.log.t_p.post = get_t_p_loglik(Y, t_t_p.c, t_t_p.s,
                                    rho, AcosPhi, AsinPhi, M, sigma)+
    kappa_t_p*cos(omega*(t_p-theta_t_p))

  new.t_p = sapply(1:N, function(i){
    CircStats::rvm(1, t_p[i]/P*2*pi, MH_kappa_t_p)/(2*pi)*P})
  
  new.t_t_p.c = cos(omega*(t+new.t_p))
  new.t_t_p.s = sin(omega*(t+new.t_p))
  new.log.t_p.post = get_t_p_loglik(Y, new.t_t_p.c, new.t_t_p.s,
                                    rho, AcosPhi, AsinPhi, M, sigma)+
    kappa_t_p*cos(omega*(t_p-theta_t_p))
  
  log.r = new.log.t_p.post-old.log.t_p.post
  a.MH.log.u = log(runif(1, 0, 1))
  if.accept = ifelse(log.r>a.MH.log.u, 1, 0)
  a.t_p = ifelse(if.accept, new.t_p, t_p)
  return(list(t_p = a.t_p, 
              if.accept = if.accept))
}

update_t_p_single_MH_zeta = function(Y, t, t.c, t.s, 
                                     t_p, t_t_p.c, t_t_p.s, zeta, pt,
                                     rho, AcosPhi, AsinPhi, M, sigma, 
                                     theta_t_p, kappa_t_p, kappa_t0_p, 
                                     MH_kappa_t_p, 
                                     N, G, omega, P){
  
  # propose zeta
  p_zeta1_new = get_p_zeta1(t_p, theta_t_p, kappa_t_p, kappa_t0_p, pt, omega)
  new.zeta = rbinom(N, 1, p_zeta1_new)

  kappa_old = ifelse(zeta, kappa_t_p, kappa_t0_p)
  kappa_new = ifelse(new.zeta, kappa_t_p, kappa_t0_p)

  # propose t_p from the von Mises approximation to its full conditional
  x1 = cbind(AcosPhi, AsinPhi) %*% rbind(t.c, t.s)
  x2 = cbind(AcosPhi, AsinPhi) %*% rbind(-1*t.s, t.c)
  y = Y-M
  c1.cos0 = diag(t(x1*rho/sigma)%*%y)
  c1.sin0 = diag(t(x2*rho/sigma)%*%y)
  c1.cos.new = c1.cos0+kappa_new*cos(omega*theta_t_p)
  c1.sin.new = c1.sin0+kappa_new*sin(omega*theta_t_p)
  R_alpha1.new = get_t_p_R_alpha(c1.cos, c1.sin)
  R1.new  = R_alpha1[, 1]; alpha1.new  = R_alpha1[, 2]
  new.t_p = sapply(1:N, function(i){
    a.samp = suppressWarnings(circular::rvonmises(1, alpha1.new[i], R1.new[i]))
    a.samp = adjust.t_p(a.samp/omega, omega)
  })
  new.t_t_p.c = cos(omega*(t+new.t_p))
  new.t_t_p.s = sin(omega*(t+new.t_p))

  # log posterior ratio
  term.log.lik = get_t_p_loglik(Y, new.t_t_p.c, new.t_t_p.s, 
                                rho, AcosPhi, AsinPhi, M, sigma)-
    get_t_p_loglik(Y, t_t_p.c, t_t_p.s, rho, AcosPhi, AsinPhi, M, sigma)
  
  term.t_p.prior = kappa_new*cos(omega*(new.t_p-theta_t_p))-
    kappa_old*cos(omega*(t_p-theta_t_p))+
    log(besselI(kappa_old, 0)/besselI(kappa_new, 0))
  
  term.zeta.prior = log(ifelse(new.zeta, pt, 1-pt)/ifelse(zeta, pt, 1-pt))
  
  # proposal probability ratio
  p_zeta1_new_to_old = get_p_zeta1(new.t_p, theta_t_p, 
                                   kappa_t_p, kappa_t0_p, pt, omega)
  proposal.p.zeta_new = ifelse(new.zeta, p_zeta1_new, 1-p_zeta1_new)
  proposal.p.zeta_old = ifelse(zeta, p_zeta1_new_to_old, 1-p_zeta1_new_to_old)
  log.prop.zeta = log(proposal.p.zeta_old/proposal.p.zeta_new)
  
  c1.cos.old = c1.cos0+kappa_old*cos(omega*theta_t_p)
  c1.sin.old = c1.sin0+kappa_old*sin(omega*theta_t_p)
  R_alpha1.old = get_t_p_R_alpha(c1.cos.old, c1.sin.old)
  R1.old = R_alpha1.old[, 1]; alpha1.old = R_alpha1.old[, 2]
  proposal.p.new_t_p = sapply(1:N, function(i){
    suppressWarnings(circular::dvonmises(new.t_p[i]*omega, 
                                         alpha1.new[i], R1.new[i], log = TRUE))
  })
  proposal.p.old_t_p = sapply(1:N, function(i){
    suppressWarnings(circular::dvonmises(t_p[i]*omega, 
                                         alpha1.old[i], R1.old[i], log = TRUE))
  })
  log.prop.t_p = proposal.p.old_t_p-proposal.p.new_t_p
  
  log.r = term.log.lik+term.t_p.prior+term.zeta.prior+log.prop.zeta+log.prop.t_p
  
  a.MH.log.u = log(runif(1, 0, 1))
  if.accept = ifelse(log.r>a.MH.log.u, 1, 0)
  a.t_p = ifelse(if.accept, new.t_p, t_p)
  a.zeta = ifelse(if.accept, new.zeta, zeta)
  return(list(t_p = a.t_p, 
              zeta = a.zeta,
              if.accept = if.accept))
}

update_t_p_single_MH_zeta_cont0 = function(Y, t, t.c, t.s, 
                                     t_p, t_t_p.c, t_t_p.s, zeta, pt,
                                     rho, AcosPhi, AsinPhi, M, sigma, 
                                     theta_t_p, kappa_t_p, kappa_t0_p, 
                                     MH_kappa_t_p, 
                                     N, G, omega, P){
  
  # propose zeta
  p_zeta1_new = get_p_zeta1_Y(Y, t.c, t.s,
                              t_p, t_t_p.c, t_t_p.s,
                              rho, AcosPhi, AsinPhi, M, sigma,
                              theta_t_p, kappa_t_p, kappa_t0_p, pt, omega)
  new.zeta = rbinom(N, 1, p_zeta1_new)

  kappa_old = ifelse(zeta, kappa_t_p, kappa_t0_p)
  kappa_new = ifelse(new.zeta, kappa_t_p, kappa_t0_p)

  # propose t_p from the von Mises approximation to its full conditional
  x1 = cbind(AcosPhi, AsinPhi) %*% rbind(t.c, t.s)
  x2 = cbind(AcosPhi, AsinPhi) %*% rbind(-1*t.s, t.c)
  y = Y-M
  c1.cos0 = diag(t(x1*rho/sigma)%*%y)
  c1.sin0 = diag(t(x2*rho/sigma)%*%y)
  c1.cos.new = c1.cos0+kappa_new*cos(omega*theta_t_p)
  c1.sin.new = c1.sin0+kappa_new*sin(omega*theta_t_p)
  R_alpha1.new = get_t_p_R_alpha(c1.cos, c1.sin)
  R1.new  = R_alpha1[, 1]; alpha1.new  = R_alpha1[, 2]
  new.t_p = sapply(1:N, function(i){
    a.samp = suppressWarnings(circular::rvonmises(1, alpha1.new[i], R1.new[i]))
    a.samp = adjust.t_p(a.samp/omega, omega)
  })
  new.t_t_p.c = cos(omega*(t+new.t_p*new.zeta))
  new.t_t_p.s = sin(omega*(t+new.t_p*new.zeta))

  # log posterior ratio
  term.log.lik = get_t_p_loglik(Y, new.t_t_p.c, new.t_t_p.s, 
                                rho, AcosPhi, AsinPhi, M, sigma)-
    get_t_p_loglik(Y, t_t_p.c, t_t_p.s, rho, AcosPhi, AsinPhi, M, sigma)
  
  term.t_p.prior = kappa_new*cos(omega*(new.t_p-theta_t_p))-
    kappa_old*cos(omega*(t_p-theta_t_p))+
    log(besselI(kappa_old, 0)/besselI(kappa_new, 0))
  
  term.zeta.prior = log(ifelse(new.zeta, pt, 1-pt)/ifelse(zeta, pt, 1-pt))
  
  # proposal probability ratio
  p_zeta1_new_to_old = get_p_zeta1_Y(Y, t.c, t.s, 
                                     new.t_p, new.t_t_p.c, new.t_t_p.s, 
                                     rho, AcosPhi, AsinPhi, M, sigma, 
                                     theta_t_p, kappa_t_p, kappa_t0_p, pt, omega)
  proposal.p.zeta_new = ifelse(new.zeta, p_zeta1_new, 1-p_zeta1_new)
  proposal.p.zeta_old = ifelse(zeta, p_zeta1_new_to_old, 1-p_zeta1_new_to_old)
  log.prop.zeta = log(proposal.p.zeta_old/proposal.p.zeta_new)
  
  c1.cos.old = c1.cos0+kappa_old*cos(omega*theta_t_p)
  c1.sin.old = c1.sin0+kappa_old*sin(omega*theta_t_p)
  R_alpha1.old = get_t_p_R_alpha(c1.cos.old, c1.sin.old)
  R1.old = R_alpha1.old[, 1]; alpha1.old = R_alpha1.old[, 2]
  proposal.p.new_t_p = sapply(1:N, function(i){
    suppressWarnings(circular::dvonmises(new.t_p[i]*omega, 
                                         alpha1.new[i], R1.new[i], log = TRUE))
  })
  proposal.p.old_t_p = sapply(1:N, function(i){
    suppressWarnings(circular::dvonmises(t_p[i]*omega, 
                                         alpha1.old[i], R1.old[i], log = TRUE))
  })
  log.prop.t_p = proposal.p.old_t_p-proposal.p.new_t_p
  
  log.r = term.log.lik+term.t_p.prior+term.zeta.prior+log.prop.zeta+log.prop.t_p
  
  a.MH.log.u = log(runif(1, 0, 1))
  if.accept = ifelse(log.r>a.MH.log.u, 1, 0)
  a.t_p = ifelse(if.accept, new.t_p, t_p)
  a.zeta = ifelse(if.accept, new.zeta, zeta)
  return(list(t_p = a.t_p, 
              zeta = a.zeta,
              if.accept = if.accept))
}

update_t_p_single_slice = function(Y,  t.c, t.s, t_p, 
                                   rho, AcosPhi, AsinPhi, M, sigma, 
                                   theta_t_p, kappa_t_p, 
                                   N, G, omega, P, save.file2){

  x1 = cbind(AcosPhi, AsinPhi) %*% rbind(t.c, t.s)
  x2 = cbind(AcosPhi, AsinPhi) %*% rbind(-1*t.s, t.c)
  y = Y-M
  # component 1, frequency omega
  c1.cos = diag(t(x1*rho/sigma)%*%y)+kappa_t_p*cos(omega*theta_t_p)
  c1.sin = diag(t(x2*rho/sigma)%*%y)+kappa_t_p*sin(omega*theta_t_p)
  # component 2, frequency 2*omega
  c2.1 = apply(x1^2*rho/sigma-x2^2*rho/sigma, 2, sum)*(-1/4)
  c2.2 = apply(x1*x2*rho/sigma, 2, sum)/(-2)
  R_alpha1 = get_t_p_R_alpha(c1.cos, c1.sin)
  R_alpha2 = get_t_p_R_alpha(c2.1, c2.2)
  R1 = R_alpha1[, 1]; alpha1 = R_alpha1[, 2]
  R2 = R_alpha2[, 1]; alpha2 = R_alpha2[, 2]

  new.t_p = sapply(1:N, function(i){
    a.t_p=NA
      fx = function(t_p){
        R1[i]*cos(omega*t_p-alpha1[i])+R2[i]*cos(2*omega*t_p-alpha2[i])
      }
      fxdx = function(t_p){
        # derivative up to the positive factor omega; only its sign is used
        -1*R1[i]*sin(omega*t_p-alpha1[i])-2*R2[i]*sin(2*omega*t_p-alpha2[i])
      }

      # locate the global minimum and maximum of fx on [-P/2, P/2]
      initial_guesses = seq(-P/2, P/2, by=1)
      find_root = function(guess) {
        tryCatch(
          rootSolve::multiroot(fxdx, start = guess, maxiter = 1000, useFortran = TRUE)$root,
          error = function(e) NA
        )
      }

      roots = sapply(initial_guesses, find_root)
      roots = unique(na.omit(roots))
      roots = sort(cleanRoot(roots, range = c(-1*P/2, P/2)))
      
      roots.type = rootsMinMax(fxdx, roots, range = c(-1*P/2, P/2))
      global.min.t = roots[which.min(fx(roots))]
      global.min = fx(global.min.t)
      global.max.t = roots[which.max(fx(roots))]
      global.max = fx(global.max.t)
      global.range = global.max-global.min

      # slice sampling: log u = fx(t_p) - Exp(1)
      z = rexp(1)
      log.u = fx(t_p[i])-z
      # resolution below which differences in t_p are ignored (0.24 h at P = 24)
      ef.abs = P/100

      # slice covers the whole circle: draw uniformly
      if(log.u-global.min< 1e-3|
         fx(global.min.t+ef.abs)>log.u|
         fx(global.min.t-ef.abs)>log.u){
        a.t_p = runif(1, -1*P/2, P/2)
      # slice is within ef.abs of the global maximum: take the maximum
      }else if(fx(global.max.t-ef.abs)<log.u|
               fx(global.max.t+ef.abs)<log.u|
               ((fxdx(global.max.t-ef.abs)*fxdx(global.max.t-ef.abs)>0)&
                abs(global.max-log.u)<global.range*1e-2)){
        a.t_p = global.max.t
      }else{
        fx1 = function(t){
          fx(t)-log.u
        }
        
        root.out = list(R1 = R1[i], R2 = R2[i], R22 = R2[i],
                        alpha1 = alpha1[i], alpha2 = alpha2[i],
                        tp = t_p[i],
                        i=i,
                        z = z)
        save(root.out, file = save.file2)
        # endpoints of the slice are the roots of fx - log.u
        initial_guesses = seq(-P/2, P/2, by=1)
        find_root = function(guess) {
          tryCatch(
            rootSolve::multiroot(fx1, start = guess, maxiter = 1000, useFortran = TRUE)$root,
            error = function(e) NA
          )
        }
        roots = sapply(initial_guesses, find_root)
        roots = unique(na.omit(roots))
        roots = sort(cleanRoot(roots, range = c(-1*P/2, P/2)))

        roots = roots[check.roots(roots, fx1, eps = P/100)]
        dxroots = fxdx(roots)

        # one interval; wrap it if it crosses the seam
        if(length(roots)==2){
          if(dxroots[1]<0&dxroots[2]>0){
            roots[1] = roots[1]+P
            roots = sort(roots)
          }
          a.t_p = stats::runif(1, roots[1], roots[2])

        # two intervals, chosen in proportion to their length
        }else if(length(roots)==4){
          if(dxroots[1]<0&dxroots[2]>0&dxroots[3]<0&dxroots[4]>0){
            roots[1] = roots[1]+P
            roots = sort(roots)
            dxroots = fxdx(roots)
          }
          if(dxroots[1]>0&dxroots[2]<0&dxroots[3]>0&dxroots[4]<0){
            int1 = roots[2] - roots[1]
            int2 = roots[4] - roots[3]
            a.t_p1 = stats::runif(1, roots[1], roots[2])
            a.t_p2 = stats::runif(1, roots[3], roots[4])
            a.t_p = ifelse(stats::rbinom(1, 1, int1/(int1+int2)), 
                           a.t_p1, a.t_p2)
          }else{
            root.out = list(R1 = R1[i], R2 = R2[i],
                            alpha1 = alpha1[i], alpha2 = alpha2[i],
                            tp = t_p[i],
                            case=4, 
                            z = z)
            save(root.out, file = save.file2)
            stop("to adjust: ", paste0("roots=", roots,  ", dx=", dxroots))
          }
        }else if(length(roots)==3){
          # the interval around the maximum crosses the seam
          if(all(roots>global.max.t)){
            roots[3]=roots[3]-P
            roots = sort(roots)
          }else if(all(roots<global.max.t)){
            roots[1]=roots[1]+P
            roots = sort(roots)
          }
          
          if((roots[1]-global.max.t)*(roots[2]-global.max.t)<0){
            # the global maximum lies between roots[1] and roots[2]
            a.t_p = stats::runif(1, roots[1], roots[2])
          }else{
            a.t_p = stats::runif(1, roots[2], roots[3])
          }
        }

      }
      
      if(is.na(a.t_p)){
        root.out = list(R1 = R1[i], R2 = R2[i], R22 = R2[i],
                        alpha1 = alpha1[i], alpha2 = alpha2[i],
                        tp = t_p[i],
                        i=i,
                        case="stillNA",
                        z = z)
        save(root.out, file = save.file2)
        stop("stillNA inside")
      }
    return(a.t_p)
  })

  return(list(t_p = new.t_p))
}

cleanRoot = function(roots, range = c(-12, 12)){
  lower = range[1]
  upper = range[2]
  roots = roots[roots>lower&roots<upper]
  if(length(roots)==0){
    roots.keep = "No Roots in range"
  }else{
    roots.keep = roots[1]
    if(length(roots)>1){
      for(i in 2:length(roots)){
        if(min(abs(roots[i]-roots.keep))>0.01){
          roots.keep = c(roots.keep, roots[i])
        }
      }
    }
  }
  return(roots.keep)
}

rootsMinMax = function(f, roots, range){
  range.min = range[1]
  range.max = range[2]
  n = length(roots)
  consecutive.mean = sapply(1:(n-1), function(i){
    mean(roots[i:(i+1)])
  })
  check.points = c(range.min, consecutive.mean, range.max)
  check.val = f(check.points)
  
  MinMax=sapply(1:n, function(i){
    if(check.val[i]*check.val[i+1]<0){
      if(check.val[i]<0){
        type="min"
      }else{
        type="max"
      }
    }else{
      type="saddle"
    }
    return(type)
  })
  return(MinMax)
}

acos_largeRange = function(x){
  if(abs(x)<1){
    acos(x)
  }else if(x<= -1){
    99
  }else{
    NA
  }
}

acos_largeRange2 = Vectorize(acos_largeRange)

getRange = function(Y, alpha){
  
  lower = (-acos_largeRange2(Y)+alpha)
  upper = (acos_largeRange2(Y)+alpha)
  
  correct.ind = which(abs(lower)>50)
  lower[correct.ind] = 0
  upper[correct.ind] = 2*pi+1 # a range wider than 2*pi means no truncation
  out = cbind(lower, upper)
  
  return(out)
}

rtruncVM = function(a, b, alpha, kappa){
  if(a>=b){
    stop( "argument a is greater than or equal to b" )
  }

  if(b-a>=2*pi){ #no truncation
    new.k = circular::rvonmises(1, alpha, kappa)
  }else{
    # shift alpha by 2*pi so that a and b lie within pi of it
    nomove = abs(a-alpha)<pi&abs(b-alpha)<pi
    left = abs(a-alpha+2*pi)<pi&abs(b-alpha+2*pi)<pi
    right = abs(a-alpha-2*pi)<pi&abs(b-alpha-2*pi)<pi     

    if(nomove|left|right){ #the truncation is in one period (alpha-pi, alpha+pi)
      alpha = ifelse(left, alpha-2*pi, 
                     ifelse(right, alpha+2*pi, alpha))
      
    }else{#the truncation breaks into two periods
      
      p.interval1 = as.numeric(circular::pvonmises(b, alpha, kappa, tol=1e-25))
      p.interval2 = 1- as.numeric(circular::pvonmises(a, alpha, kappa, tol=1e-25))
      p.scale = 1/(p.interval1+p.interval2)
      p.interval1 = p.interval1*p.scale
      # choose one interval to sample from
      choose1 = rbinom(1, 1, p.interval1)
      if(choose1){
        a = alpha-pi+1e-5
      }else{
        b = alpha+pi-1e-5
      }
    }
    G.min = as.numeric(circular::pvonmises(a, alpha, kappa, tol=1e-25))
    G.max = as.numeric(circular::pvonmises(b, alpha, kappa, tol=1e-25)) 
    u = runif(1)
    a_new_p = G.min+u*(G.max-G.min)
    if(G.max-G.min<0){stop("G.max-G.min<0")}
    if(a_new_p>G.max){
      a_new_p=G.max # guards against rounding above G.max
    }

    new.k = as.numeric(circular::qvonmises(a_new_p, alpha, kappa))
  }
  return(new.k)
}

adjust.t_p = function(t, omega){
  period2 = pi/omega
  t = ifelse(t<period2*-1, t+2*period2, 
             ifelse(t>period2, t-2*period2, t))
}

get_t_p_R_alpha = function(x1, x2){
  # x1 and x2 are the cos and sin coefficients
  R = sqrt(x1^2+x2^2)
  alpha = atan2(x2, x1)
  return(cbind(R, alpha))
}

check.roots = function(x, f, eps=P/100){
  f(x-eps)*f(x+eps)<0
}
