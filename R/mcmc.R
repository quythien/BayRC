#' Run reversible-jump MCMC for a single circadian dataset
#'
#' @title CB MCMC single reversible-jump slice sampler
#'
#' @description
#' Runs the BayRC RJMCMC sampler for a single-species/condition circadian
#' expression dataset.  Within each iteration the sampler proposes
#' model-space jumps (rho: 0/1 rhythmicity indicator) via an RJMCMC
#' acceptance ratio and updates amplitude, phase, MESOR, and residual
#' variance via Gibbs or slice steps.
#'
#' @param Data.list A named list with three elements: \code{data} (G x N
#'   data.frame of expression values), \code{time} (length-N numeric
#'   Zeitgeber time vector in hours), and \code{gname} (length-G character
#'   vector of gene identifiers).
#' @param Init.value List returned by \code{CB_init_single} containing
#'   starting values for \code{rho}, \code{M}, \code{A}, \code{phi}, and
#'   \code{sigma}.
#' @param P Numeric; circadian period in hours (default 24).
#' @param iteration Integer; total number of MCMC iterations including
#'   burn-in (default 3000).
#' @param NP_Z1 Integer; quadrature nodes for the phi marginal (default 64).
#' @param A.max Numeric; upper truncation bound for amplitude, either one
#'   value or one per gene. \code{NULL} (the default) uses half the observed
#'   range of each gene; \code{Inf} leaves the amplitude prior one-sided on
#'   \code{[A.min, Inf)}.
#' @param thin Integer; thinning interval; every \code{thin}-th
#'   post-burn-in sample is stored (default 1).
#' @param n.burn Integer; number of burn-in iterations discarded before
#'   storing samples (default 1000).
#' @param seed Integer; random seed for reproducibility (default 15213).
#' @param diagnostics Logical; if \code{TRUE} (default), run
#'   \code{mcmc_diagnostics()} after sampling and attach the result as
#'   \code{$diagnostics}. The same diagnostics can be computed later from the
#'   returned chains.
#' @param p_rhythmic Numeric vector of length G; prior probability
#'   Pr(rho_g = 1) for each gene.  Enters the RJMCMC log-prior odds as
#'   log(p / (1 - p)).  Default is \code{rep(0.2, 100)}.
#' @param rj.p.stay Numeric in (0, 1); probability of skipping the
#'   between-model move (staying in the current model) at each iteration
#'   (default 0.5).
#' @param A_prior Character; prior on amplitude.  One of
#'   \code{"Jeffreys_OLS_condi"}, \code{"trunc_Normal"},
#'   \code{"Jeffreys"}, \code{"Jeffreys_ridge"}, \code{"sq_expo"}, or
#'   \code{"gamma"} (default \code{"Jeffreys_OLS_condi"}).
#' @param mu_A Numeric; mean of the truncated-Normal amplitude prior when
#'   \code{A_prior = "trunc_Normal"} (default 1).
#' @param sigma_A Numeric; variance of the truncated-Normal amplitude prior
#'   (default \code{10^2}).
#' @param A.min Numeric; lower truncation bound for amplitude (default 0).
#' @param rj.phi Logical; if \code{TRUE}, phase phi is jointly proposed in
#'   the RJMCMC birth/death step (default \code{TRUE}).
#' @param rj.A Logical; if \code{TRUE}, amplitude A is jointly proposed in
#'   the RJMCMC birth/death step (default \code{TRUE}).
#' @param mu_M Numeric; prior mean for MESOR (default 0).
#' @param sigma_M Numeric; prior variance for MESOR (default \code{10^2}).
#' @param sigma_prior_v Numeric; degrees-of-freedom hyperparameter of the
#'   inverse-gamma prior on residual variance (default 4).
#' @param sigma_prior_s Numeric; scale hyperparameter of the inverse-gamma
#'   prior on residual variance (default 1).
#'
#' @return A named list of G x K matrices (K = floor(iteration / thin)):
#'   \describe{
#'     \item{rho}{Posterior samples of the rhythmicity indicator in \{0, 1\}.}
#'     \item{phi}{Posterior samples of acrophase in hours (Zeitgeber scale).}
#'     \item{A}{Posterior samples of amplitude.}
#'     \item{M}{Posterior samples of MESOR.}
#'     \item{sigma}{Posterior samples of residual variance.}
#'     \item{log.r1, log.r1.SS, log.r3_A, log.r3_phi}{RJMCMC log-ratio
#'       components (for diagnostics).}
#'     \item{if.accept.rj}{Binary matrix indicating accepted RJMCMC moves.}
#'     \item{Z1_gap}{3 x K matrix: median, 99th percentile and maximum across
#'       genes of \code{abs(Z1_check - logZ1)} at each stored iteration.}
#'   }
#'   All G x K matrices have \code{rownames} equal to \code{Data.list$gname}.
#'   After calling \code{match_symbols()}, the attributes
#'   \code{attr(rho, "symbols")} and \code{attr(rho, "RHYindex")} are set.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sim <- CBt_sim_data()
#' dat <- list(data  = as.data.frame(sim$data[[1]]$dat),
#'             time  = sim$data[[1]]$x$time,
#'             gname = paste0("G", seq_len(nrow(sim$data[[1]]$dat))))
#' init <- CB_init_single(dat)
#' mcmc_out <- CB_MCMC_single_rj_slice(dat, init,
#'               iteration = 200, thin = 10, n.burn = 100)
#' }
CB_MCMC_single_rj_slice = function(Data.list, Init.value, P = 24,
                                iteration = 3000, thin = 1, n.burn = 1000,
                                seed = 15213, diagnostics = TRUE,

                                p_rhythmic=rep(0.2, 100), rj.p.stay = 0.5,
                                A_prior= "Jeffreys_OLS_condi",
                                mu_A = 1, sigma_A = 10^2, A.min = 0, # truncated-Normal prior on A
                                A.max = NULL,
                                rj.phi = TRUE, rj.A = TRUE,
                                mu_M = 0, sigma_M = 10^2,
                                sigma_prior_v = 4, sigma_prior_s = 1,
                                NP_Z1 = 64   # phi-quadrature nodes for Z1
){

  omega = 2*pi/P
  Y = as.matrix(Data.list[[1]])
  ## default bound is half the observed range of each gene
  if(is.null(A.max))
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
  
  for(iter in 1:n.burn){
    
    set.seed(seed+iter)
    print(paste0(A_prior, " Burn iter = ", iter))
    
    # update btw model --------------------------------------------------------
    stay = stats::rbinom(1, 1, rj.p.stay)
    if(!stay){
      rho.res = RJMCMC_single_slice(Y, t.c, t.s, N,
                                    t.c.sum, t.s.sum,
                                    c.t.phi, s.t.phi, cs.t.phi,
                                    y.t.c.sum, y.t.s.sum,
                                    AcosPhi, AsinPhi, A, phi,
                                    M, sigma, p_rhythmic, rho,
                                    rj.phi, rj.A,
                                    A_prior,
                                    mu_A, sigma_A, A.min, A.max, 
                                    omega, G, P, NP = NP_Z1)

      if.accept.rj = rho.res$rho!=rho
      rho = rho.res$rho
    }
    
    
    # within model move -------------------------------------------------------
    M = update_M_single(Y, t.c, t.s, N,
                        AcosPhi, AsinPhi, sigma, rho,
                        sigma_M, mu_M, omega, G)

    csAphi = update_A_phi_slice(Y, X, XtX, beta_hat0,
                                beta_cov0_11, beta_cov0_22, beta_cov0_rho,
                                beta_mean0_1, beta_mean0_2,
                                M, sigma, omega,
                                AcosPhi, AsinPhi, G, P, A_prior,
                                A.min, A.max, 
                                mu_A, sigma_A)
    AcosPhi = csAphi$AcosPhi
    AsinPhi = csAphi$AsinPhi
    A = sqrt(AcosPhi^2+AsinPhi^2)
    phi = sapply(atan2(AsinPhi, AcosPhi), adjust.to.2pi)/omega
    stop_on_na(rho, AcosPhi, AsinPhi, iter = iter, stage = "burn-in")

    sigma = update_sigma_single(Y, t.c, t.s, N, AcosPhi, AsinPhi, M, rho,
                                sigma_prior_v, sigma_prior_s, omega, G)
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
  ## keep the spread of the discrepancy, not the full G x K matrices
  Z1_gap = function(res) {
    d = abs(res$Z1_check - res$logZ1)
    d = d[is.finite(d)]
    if (!length(d)) c(median = NA_real_, q99 = NA_real_, max = NA_real_)
    else c(median = stats::median(d),
           q99 = unname(stats::quantile(d, 0.99)),
           max = max(d))
  }
  Z1_gap.store = Z1_gap(rho.res)
  if.accept.rj.store = if.accept.rj
  time0 = Sys.time()

  # Start sampling ----------------------------------------------------------
  for(iter in 1:(iteration-n.burn)){
    
    time0 = Sys.time()
    set.seed(seed+iter)
    print(paste0("iter = ", iter))
    
    # update btw model --------------------------------------------------------
    stay = stats::rbinom(1, 1, rj.p.stay)
    if(stay){
      rho.store = cbind(rho.store, rho)
      if.accept.rj = cbind(if.accept.rj, -1) 
      log.r1.store = cbind(log.r1.store, -99)
      log.r1.SS.store = cbind(log.r1.SS.store, -99) 
      Z1_gap.store = cbind(Z1_gap.store, c(NA_real_, NA_real_, NA_real_))
    }else{
      rho.res = RJMCMC_single_slice(Y, t.c, t.s, N,
                                    t.c.sum, t.s.sum,
                                    c.t.phi, s.t.phi, cs.t.phi,
                                    y.t.c.sum, y.t.s.sum,
                                    AcosPhi, AsinPhi, A, phi,
                                    M, sigma, p_rhythmic, rho,
                                    rj.phi, rj.A,
                                    A_prior,
                                    mu_A, sigma_A, A.min, A.max,
                                    omega, G, P, NP = NP_Z1)
      if.accept.rj = cbind(if.accept.rj, rho.res$rho!=rho) # a change in rho is an accepted jump
      log.r1.store = cbind(log.r1.store, rho.res$log.r1)
      log.r1.SS.store = cbind(log.r1.SS.store, rho.res$log.r1.SS) 
      log.r3_A.store = cbind(log.r3_A.store, rho.res$log.r3_A)
      log.r3_phi.store = cbind(log.r3_phi.store, rho.res$log.r3_phi)
      Z1_gap.store = cbind(Z1_gap.store, Z1_gap(rho.res))
      rho = rho.res$rho
      rho.store = cbind(rho.store, rho)
    }    
    time1 = Sys.time()

    # within model move -------------------------------------------------------

    M = update_M_single(Y, t.c, t.s, N,
                        AcosPhi, AsinPhi, sigma, rho,
                        sigma_M, mu_M, omega, G)
    M.store = cbind(M.store, M)
    time2 = Sys.time()

    csAphi = update_A_phi_slice(Y, X, XtX, beta_hat0,
                                beta_cov0_11, beta_cov0_22, beta_cov0_rho,
                                beta_mean0_1, beta_mean0_2,
                                M, sigma, omega,
                                AcosPhi, AsinPhi, G, P, A_prior,
                                A.min, A.max, 
                                mu_A, sigma_A)
    
    AcosPhi = csAphi$AcosPhi
    AsinPhi = csAphi$AsinPhi
    A = sqrt(AcosPhi^2+AsinPhi^2)
    phi = sapply(atan2(AsinPhi, AcosPhi), adjust.to.2pi)/omega
    stop_on_na(rho, AcosPhi, AsinPhi, iter = iter, stage = "sampling")

    AcosPhi.store = cbind(AcosPhi.store, AcosPhi);
    AsinPhi.store = cbind(AsinPhi.store, AsinPhi); 
    A.store = cbind(A.store, A)
    phi.store = cbind(phi.store, phi)
    
    time4 = Sys.time()

    sigma = update_sigma_single(Y, t.c, t.s, N, AcosPhi, AsinPhi, M, rho,
                                sigma_prior_v, sigma_prior_s, omega, G)
    sigma.store = cbind(sigma.store, sigma)
    time5 = Sys.time()

    save = list(rho = rho.store,
                M = M.store,
                AcosPhi = AcosPhi.store,
                AsinPhi = AsinPhi.store, 
                A = A.store, 
                phi = phi.store,
                sigma = sigma.store,
                Z1_gap = Z1_gap.store,
                if.accept.rj = if.accept.rj)
  }
  if (thin > 1) {
    K <- ncol(save$rho)
    keep <- seq(1, K, by = thin)
    save <- lapply(save, function(x)
      if (is.matrix(x) && ncol(x) == K) x[, keep, drop = FALSE] else x)
  }

  gname <- Data.list[[3]]
  if (!is.null(gname) && length(gname) == G) {
    for (mat_name in names(save)) {
      if (is.matrix(save[[mat_name]]) && nrow(save[[mat_name]]) == G) {
        rownames(save[[mat_name]]) <- gname
      }
    }
  }
  if (diagnostics) save$diagnostics <- mcmc_diagnostics(save, P = P)

  return(save)
}

# Updating functions ------------------------------------------------------

#' Gibbs update for MESOR
#'
#' @description
#' Samples the MESOR (M) vector from its Gaussian full conditional given the
#' current amplitude, phase, residual variance, and rhythmicity indicators.
#'
#' @param Y G x N expression matrix.
#' @param t.c Length-N vector of cos(omega * t).
#' @param t.s Length-N vector of sin(omega * t).
#' @param N Integer; number of samples.
#' @param AcosPhi Length-G vector A * cos(omega * phi).
#' @param AsinPhi Length-G vector A * sin(omega * phi).
#' @param sigma Length-G residual variance vector.
#' @param rho Length-G rhythmicity indicator vector in \{0, 1\}.
#' @param sigma_M Numeric; prior variance on M.
#' @param mu_M Numeric; prior mean on M.
#' @param omega Numeric; angular frequency 2*pi/P.
#' @param G Integer; number of genes.
#'
#' @return Length-G numeric vector of updated MESOR samples.
#'
#' @keywords internal
update_M_single = function(Y, t.c, t.s, N,
                           AcosPhi, AsinPhi, sigma, rho,
                           sigma_M, mu_M, omega, G){
  a.var = 1/(N/sigma+1/sigma_M)
  b.c = AcosPhi*rho #G*1 mat
  b.s = AsinPhi*rho
  a.cosine = cbind(b.c, b.s) %*% rbind(t.c, t.s)
  a.mean = (apply((Y-a.cosine), 1, sum)/sigma+mu_M/sigma_M)*a.var
  a.M = stats::rnorm(G, a.mean, sqrt(a.var))
}

#' Gibbs update for residual variance
#'
#' @description
#' Samples the gene-wise residual variance vector (sigma) from its
#' inverse-gamma full conditional.
#'
#' @param Y G x N expression matrix.
#' @param t.c Length-N vector of cos(omega * t).
#' @param t.s Length-N vector of sin(omega * t).
#' @param a.N Integer; number of samples.
#' @param AcosPhi Length-G vector A * cos(omega * phi).
#' @param AsinPhi Length-G vector A * sin(omega * phi).
#' @param M.vec Length-G MESOR vector.
#' @param rho Length-G rhythmicity indicator vector.
#' @param sigma_prior_v Numeric; degrees-of-freedom hyperparameter.
#' @param sigma_prior_s Numeric; scale hyperparameter.
#' @param omega Numeric; angular frequency.
#' @param G Integer; number of genes.
#'
#' @return Length-G numeric vector of updated residual variance samples.
#'
#' @keywords internal
update_sigma_single = function(Y, t.c, t.s, a.N,
                               AcosPhi, AsinPhi, M.vec, rho,
                               sigma_prior_v, sigma_prior_s, omega, G){
  b.c = AcosPhi*rho #G*1 mat
  b.s = AsinPhi*rho
  a.cosine = cbind(M.vec, b.c, b.s) %*% rbind(1, t.c, t.s)
  vs2 = sigma_prior_v*sigma_prior_s
  sigma_post_vs2 = (apply((Y-a.cosine)^2, 1, sum)+vs2)/2
         a.sigma = invgamma::rinvgamma(G, shape = (sigma_prior_v+a.N)/2, rate = sigma_post_vs2)
         a.sigma = ifelse(a.sigma<1e-5, 1e-5, a.sigma)
  return(a.sigma)
}

#' Slice-sampling update for amplitude and phase
#'
#' @description
#' Updates the pair (A, phi) within the current model by slice sampling the
#' Cartesian coordinates \code{AcosPhi} and \code{AsinPhi}, one conditional on
#' the other. The slice height sets a per-gene radius \code{A_max}, and each
#' coordinate is then drawn from a truncated Normal on the two arcs the radius
#' allows.
#'
#' @param Y G x N expression matrix.
#' @param X N x 2 design matrix of cos and sin terms.
#' @param XtX Crossproduct of \code{X}.
#' @param beta_hat0 Precomputed OLS projection matrix.
#' @param beta_cov0_11,beta_cov0_22,beta_cov0_rho Prior covariance terms for
#'   the two coordinates.
#' @param beta_mean0_1,beta_mean0_2 Conditional-mean coefficients.
#' @param M,sigma Length-G MESOR and residual-variance vectors.
#' @param omega Numeric; angular frequency 2*pi/P.
#' @param AcosPhi,AsinPhi Length-G current coordinates.
#' @param G Integer; number of genes.
#' @param P Numeric; period in hours.
#' @param A_prior Character; amplitude prior.
#' @param A.min,A.max Numeric; amplitude truncation bounds.
#' @param mu_A,sigma_A Numeric; truncated-Normal amplitude prior parameters.
#'
#' @return A list with the updated \code{AcosPhi} and \code{AsinPhi}.
#'
#' @keywords internal
update_A_phi_slice = function(Y, X, XtX, beta_hat0,
                              beta_cov0_11, beta_cov0_22, beta_cov0_rho, 
                              beta_mean0_1, beta_mean0_2, 
                              M, sigma, omega, 
                              AcosPhi, AsinPhi, G, P, A_prior,
                              A.min, A.max, 
                              mu_A, sigma_A  # A ~ truncated Normal(mu_A, sigma_A)
){
  beta1 = AcosPhi
  beta2 = AsinPhi
  Y_new = Y-M
  A = sqrt(beta1^2+beta2^2)
  # slice height function for the truncated-Normal prior
  f_upper_trunc = function(A0, u0){
    1/A0*exp(-1/sigma_A*(A0-mu_A)^2)-u0
  }
  
  if(A_prior=="Jeffreys_OLS_condi"){
    
# A_prior=="Jeffreys_OLS_condi" -------------------------------------------
    
    beta_hat = Y_new%*%beta_hat0
    u = runif(G, 1e-5, 1/A)
    u = ifelse(u>1e3, 1e3, u) # cap keeps the radius 1/u away from zero
    beta_cov_11 = beta_cov0_11*sigma
    beta_cov_22 = beta_cov0_22*sigma
    max_beta1 = sqrt(1/u^2-beta2^2)
    beta1_mean = beta_hat[, 1]+beta_mean0_1*(beta2-beta_hat[, 2])
    beta1_var = beta_cov_11*(1-beta_cov0_rho^2)
    beta1_new = truncnorm::rtruncnorm(G, -1*max_beta1, max_beta1, 
                                      beta1_mean, sqrt(beta1_var))
    max_beta2 = sqrt(1/u^2-beta1_new^2)
    beta2_mean = beta_hat[, 2]+beta_mean0_2*(beta1_new-beta_hat[, 1])
    beta2_var = beta_cov_22*(1-beta_cov0_rho^2)
    beta2_new = truncnorm::rtruncnorm(G, -1*max_beta2, max_beta2, 
                                      beta2_mean, sqrt(beta2_var))

  }else if(A_prior=="Jeffreys_ridge_condi"){
    
    # A_prior=="Jeffreys_ridge_condi" -----------------------------------------
    
    # ridge solution, one gene at a time
    u = runif(G, 1e-5, 1/A)
    u = ifelse(u>1e3, 1e3, u)
    get_beta_hat0 = sapply(1:G, function(g){
      XtXinv = solve(XtX+1/u[g]^2)
      beta_hat = Y_new[g, ]%*%t(XtXinv%*%t(X))
      beta_cov0 = XtXinv%*%XtX%*%XtXinv
      beta_cov_11 = beta_cov0[1, 1]*sigma[g]
      beta_cov_22 = beta_cov0[2, 2]*sigma[g]
      beta_cov_rho = beta_cov0[1, 2]/sqrt(beta_cov_11*beta_cov_22)
      return(c(beta_cov_11, beta_cov_22, beta_cov_rho, beta_hat))
    })
    get_beta_hat = t(get_beta_hat0)
    beta_cov_11 = get_beta_hat[, 1]
    beta_cov_22 = get_beta_hat[, 2]
    beta_cov_rho = get_beta_hat[, 3]
    beta_hat = get_beta_hat[, 4:5]

    max_beta1 = sqrt(1/u^2-beta2^2)
    beta1_mean = beta_hat[, 1]+beta_cov_rho*sqrt(beta_cov_11/beta_cov_22)*(beta2-beta_hat[, 2])
    beta1_var = beta_cov_11*(1-beta_cov_rho^2)
    beta1_new = truncnorm::rtruncnorm(G, -1*max_beta1, max_beta1, 
                                      beta1_mean, sqrt(beta1_var)) 
    max_beta2 = sqrt(1/u^2-beta1_new^2)
    beta2_mean = beta_hat[, 2]+beta_cov_rho*sqrt(beta_cov_22/beta_cov_11)*(beta1_new-beta_hat[, 1])
    beta2_var = beta_cov_22*(1-beta_cov_rho^2)
    beta2_new = truncnorm::rtruncnorm(G, -1*max_beta2, max_beta2, 
                                      beta2_mean, sqrt(beta2_var)) 
  }else if(A_prior=="trunc_Normal_OLS_condi"){
    
    # A_prior=="trunc_Normal_OLS_condi" ---------------------------------------
    
    beta_hat = Y_new%*%beta_hat0
    u = runif(G, 1e-3, f_upper_trunc(A, 0))
    u = ifelse(u>1e3, 1e3, u)
    beta_cov_11 = beta_cov0_11*sigma
    beta_cov_22 = beta_cov0_22*sigma

    # radius of the slice at height u
    A_max = sapply(1:G, function(g){
      uniroot(f_upper_trunc, c(0, 1e3), u0=u[g], tol=1e-6)$root})
    max_beta1 = sqrt(A_max^2-beta2^2)
    max_beta1[is.na(max_beta1)]=1e-3
    beta1_mean = beta_hat[, 1]+beta_mean0_1*(beta2-beta_hat[, 2])
    beta1_var = beta_cov_11*(1-beta_cov0_rho^2)

    beta1_new_1 = truncnorm::rtruncnorm(G, -1*max_beta1, -1*A.min, 
                                        beta1_mean, sqrt(beta1_var)) 
    beta1_new_2 = truncnorm::rtruncnorm(G, A.min, max_beta1, 
                                        beta1_mean, sqrt(beta1_var)) 
    p1=pnorm(-1*A.min, beta1_mean, sqrt(beta1_var))-
      pnorm(-1*max_beta1, beta1_mean, sqrt(beta1_var))
    p2=pnorm(max_beta1, beta1_mean, sqrt(beta1_var))-
      pnorm(A.min, beta1_mean, sqrt(beta1_var))
    ## both masses can underflow to 0 on a narrow slice; split evenly then
    den = p1+p2
    p1 = ifelse(den>0, p1/den, 0.5)
    beta1_new = ifelse(rbinom(rep(1, G), rep(1, G), p1),
                       beta1_new_1, beta1_new_2)
    
    max_beta2 = sqrt(A_max^2-beta1_new^2)
    max_beta2[is.na(max_beta2)]=1e-3
    beta2_mean = beta_hat[, 2]+beta_mean0_2*(beta1_new-beta_hat[, 1])
    beta2_var = beta_cov_22*(1-beta_cov0_rho^2)

    beta2_new_1 = truncnorm::rtruncnorm(G, -1*max_beta2, -1*A.min, 
                                        beta2_mean, sqrt(beta2_var)) 
    beta2_new_2 = truncnorm::rtruncnorm(G, A.min, max_beta2, 
                                        beta2_mean, sqrt(beta2_var)) 
    p1=pnorm(-1*A.min, beta2_mean, sqrt(beta2_var))-
      pnorm(-1*max_beta2, beta2_mean, sqrt(beta2_var))
    p2=pnorm(max_beta2, beta2_mean, sqrt(beta2_var))-
      pnorm(A.min, beta2_mean, sqrt(beta2_var))
    den = p1+p2
    p1 = ifelse(den>0, p1/den, 0.5)
    beta2_new = ifelse(rbinom(rep(1, G), rep(1, G), p1),
                       beta2_new_1, beta2_new_2)
    
  }else if(A_prior=="trunc_Normal_ridge_condi"){
    
    # A_prior=="trunc_Normal_ridge_condi" -------------------------------------
    
    u = runif(G, 1e-5, f_upper_trunc(A, 0)) 
    u = ifelse(u>1e3, 1e3, u)
    A_max = sapply(1:G, function(g){
      uniroot(f_upper_trunc, c(0, 1000), u0=u[g])$root})
    get_beta_hat0 = sapply(1:G, function(g){
      a.A_max = A_max[g]
      XtXinv = solve(XtX+a.A_max^2)
      beta_hat = Y_new[g, ]%*%t(XtXinv%*%t(X))
      beta_cov0 = XtXinv%*%XtX%*%XtXinv
      beta_cov_11 = beta_cov0[1, 1]*sigma[g]
      beta_cov_22 = beta_cov0[2, 2]*sigma[g]
      beta_cov_rho = beta_cov0[1, 2]/sqrt(beta_cov_11*beta_cov_22)
      return(c(beta_cov_11, beta_cov_22, beta_cov_rho, beta_hat))
    })
    get_beta_hat = t(get_beta_hat0)
    beta_cov_11 = get_beta_hat[, 1]
    beta_cov_22 = get_beta_hat[, 2]
    beta_cov_rho = get_beta_hat[, 3]
    beta_hat = get_beta_hat[, 4:5]
    max_beta1 = sqrt(A_max^2-beta2^2)
    beta1_mean = beta_hat[, 1]+beta_cov_rho*sqrt(beta_cov_11/beta_cov_22)*(beta2-beta_hat[, 2])
    beta1_var = beta_cov_11*(1-beta_cov_rho^2)
    beta1_new = truncnorm::rtruncnorm(G, -1*max_beta1, max_beta1, 
                                      beta1_mean, sqrt(beta1_var)) 
    max_beta2 = sqrt(A_max^2-beta1_new^2)
    beta2_mean = beta_hat[, 2]+beta_cov_rho*sqrt(beta_cov_22/beta_cov_11)*(beta1_new-beta_hat[, 1])
    beta2_var = beta_cov_22*(1-beta_cov_rho^2)
    beta2_new = truncnorm::rtruncnorm(G, -1*max_beta2, max_beta2, 
                                      beta2_mean, sqrt(beta2_var)) 
  }
  
  stopifnot(sum(is.na(beta1_new))==0)
  stopifnot(sum(is.na(beta2_new))==0)
  
  return(list(AcosPhi = beta1_new, 
              AsinPhi = beta2_new))
}

get_post_A_params = function(Y, t.c, t.s, N, 
                             AcosPhi, AsinPhi, A, phi,  
                             M.vec, sigma.vec, 
                             A_prior, 
                             mu_A, sigma_A, omega){
  if(A_prior=="trunc_Normal_OLS_condi"){
    b.c = cos(omega*phi) #G*1 mat
    b.s = sin(omega*phi)
    t.phi = cbind(b.c, b.s) %*% rbind(t.c, t.s)
    T.x = apply(t.phi^2, 1, sum)
    Y.minu.M = Y - as.vector(M.vec)
    U.x = apply(t.phi*Y.minu.M, 1, sum)
    a.var = 1/(T.x/sigma.vec+1/sigma_A)
    a.mean = (U.x/sigma.vec+mu_A/sigma_A)*a.var
  }
  return(list(mean = a.mean, var = a.var))
}

get_post_A = function(Y, t.c, t.s, N, 
                      AcosPhi, AsinPhi, A, phi, 
                      M, sigma, 
                      A_prior, 
                      mu_A, sigma_A, A.min, A.max, # A ~ truncated Normal(mu_A, sigma_A) on (A.min, A.max)
                      omega, G){
  b.c = cos(omega*phi) #G*1 mat
  b.s = sin(omega*phi)
  t.phi = cbind(b.c, b.s) %*% rbind(t.c, t.s)
  T.x = apply(t.phi^2, 1, sum)
  Y.minu.M = Y - as.vector(M)
  U.x = apply(t.phi*Y.minu.M, 1, sum)
  
  if(A_prior=="trunc_Normal_OLS_condi"){

    a.var = 1/(T.x/sigma+1/sigma_A)
    a.mean = (U.x/sigma+mu_A/sigma_A)*a.var
    A.post = truncnorm::dtruncnorm(A, a = A.min, b = A.max, a.mean, sqrt(a.var))
    
  }else if(A_prior == "Jeffreys_OLS_condi"){

    a.var = sigma/T.x 
    a.mean = U.x/T.x 
    A.post = truncnorm::dtruncnorm(A, a = A.min, b = A.max, a.mean, sqrt(a.var))
    
  }
  return(A.post)
}

get_prior_A = function(Y, t.c, t.s, N, 
                      AcosPhi, AsinPhi, A, phi, 
                      M, sigma, 
                      A_prior, 
                      mu_A, sigma_A, A.min, A.max, # A ~ truncated Normal(mu_A, sigma_A) on (A.min, A.max)
                      omega){
  
  if(A_prior=="trunc_Normal_OLS_condi"){
    A.prior = truncnorm::dtruncnorm(A, a = A.min, b = A.max, mu_A, sqrt(sigma_A))
    
  }else if(A_prior == "Jeffreys_OLS_condi"){
    A.prior = 1
    
  }
  return(A.prior)
}

get_log_phi_post_single = function(a.Y.list, a.t.c, a.t.s, a.n,
                                   a.M.vec, AcosPhi, AsinPhi, a.sigma.vec,
                                   omega){
  #when A is present, phi is always present, so use rho = 1 for all
  b.c = AcosPhi #G*1 mat
  b.s = AsinPhi
  a.cosine.1 = cbind(a.M.vec, b.c, b.s) %*% rbind(1, a.t.c, a.t.s)
  diff1 = apply((a.Y.list- a.cosine.1)^2, 1, sum)
  
  -1/2*diff1/a.sigma.vec
}

update_phi_single = function(a.Y, a.t.c, a.t.s, a.N,
                             A.vec, M.vec, phi.vec, sigma.vec, rho.vec,
                             delta_MH_phi, omega, G, P){
  # Metropolis step with a von Mises proposal
  old.log.phi.post = get_log_phi_post_single(a.Y, a.t.c, a.t.s, a.N,
                                             M.vec, A.vec, phi.vec, sigma.vec,
                                             omega)
  new.phi = sapply(1:G, function(g){
    CircStats::rvm(1, phi.vec[g]/P*2*pi, delta_MH_phi)/(2*pi)*P})
  
  new.log.phi.post = get_log_phi_post_single(a.Y, a.t.c, a.t.s, a.N,
                                             M.vec, A.vec, new.phi, sigma.vec,
                                             omega)
  log.r = new.log.phi.post-old.log.phi.post
  a.MH.log.u = log(runif(1, 0, 1))
  if.accept = ifelse(log.r>a.MH.log.u, 1, 0)
  a.phi = ifelse(if.accept, new.phi, phi.vec)
  return(list(phi = a.phi, 
              if.accept = if.accept))
}

update_rho = function(Y.list, t.c.list, t.s.list, N.vec, 
                      A.mat, phi.mat, M.mat, sigma.mat, p.mat, 
                      omega, G){
  rho.per.group = lapply(1:length(Y.list), function(a){
    a.Y.list = Y.list[[a]]
    a.t.c = t.c.list[[a]] #n_j*1 vec
    a.t.s = t.s.list[[a]]
    b.c = A.mat[, a]*cos(omega*phi.mat[, a]) #G*1 mat
    b.s = A.mat[, a]*sin(omega*phi.mat[, a])
    a.cosine.1= cbind(M.mat[, a], b.c, b.s) %*% rbind(1, a.t.c, a.t.s)
    lik1 = apply((a.Y.list- a.cosine.1)^2, 1, sum)
    lik0 = apply((a.Y.list-
                    M.mat[, a] %*% matrix(1, ncol = N.vec[a], nrow = 1))^2,
                 1, sum)
    lik.ratio = -(1/2*sigma.mat[, a])*(lik1-lik0)
    odds_rho1 = (p.mat[, a]+1e-5)/(1-p.mat[, a]+1e-5)*
      exp(ifelse(lik.ratio>100, 100, lik.ratio))
    p.rho1 = odds_rho1/(odds_rho1+1)
    a.rho = stats::rbinom(G, 1, p.rho1)

    return(a.rho)
  })
  do.call(cbind, rho.per.group)
}

## Marginal posterior of phi under M1 and its normalizing constant Z1 ----
##
## The joint M1 kernel for one gene, conditional on M_g and sigma^2_g, is
##     Psi(A, phi) = L(A, phi) * pi(A) * pi(phi),      pi(phi) = 1/P
##     L(A, phi)   = exp{ -[A^2 T(phi) - 2 A U(phi)] / (2 sigma^2_g) }
##     T(phi) = sum_i cos^2(omega(t_i - phi)),  U(phi) = sum_i (Y_gi - M_g) c_i
##
## Collecting the A-terms and completing the square,
##     1/v_A(phi) = T(phi)/sigma^2_g + 1/sigma_A^2
##     m_A(phi)   = v_A(phi) [ U(phi)/sigma^2_g + mu_A/sigma_A^2 ]
## so that, with Z_tr the truncation constant of the prior,
##     Psi(A, phi) = (1/P) (1 / (sqrt(2 pi) sigma_A Z_tr))
##                   exp{ m_A^2/(2 v_A) - mu_A^2/(2 sigma_A^2) }
##                   exp{ -(A - m_A)^2 / (2 v_A) }.
##
## Integrating A over (A.min, A.max) leaves a Gaussian integral between two
## bounds, hence the closed form
##
##   h(phi) = (1/P) (sqrt(v_A) / (sigma_A Z_tr))
##            [ Phi((A.max - m_A)/sqrt(v_A)) - Phi((A.min - m_A)/sqrt(v_A)) ]
##            exp{ m_A^2/(2 v_A) - mu_A^2/(2 sigma_A^2) }
##
## and the marginal posterior of phi is  J(phi) = h(phi) / Z1  with
## Z1 = int_0^P h(phi) dphi.
##
## As sigma_A^2 -> Inf, m_A^2/(2 v_A) tends to U(phi)^2 / (2 sigma^2_g T(phi)),
## the Lomb-Scargle periodogram.
##
## h(phi) is smooth and periodic on [0, P], so the uniform-grid trapezoid rule
## converges geometrically (NP = 32 reaches machine precision). Only Z1 uses
## the grid; log h is evaluated exactly at the current phi.

## log h(phi) at one phi per gene (exact, vectorized over genes)
log_h_single = function(Y, t.c, t.s, M, sigma, phi, A_prior,
                        mu_A, sigma_A, A.min, A.max, omega, P){
  ## cmat[g, i] = cos(omega*(t_i - phi_g))
  cmat = outer(cos(omega*phi), t.c) + outer(sin(omega*phi), t.s)
  Tg = rowSums(cmat^2)
  Ug = rowSums(cmat * (as.matrix(Y) - as.vector(M)))
  log_h_core(Tg, Ug, sigma, A_prior, mu_A, sigma_A, A.min, A.max, P)
}

## shared core: works elementwise on T and U of matching shape
log_h_core = function(Tx, Ux, sigma, A_prior,
                      mu_A, sigma_A, A.min, A.max, P){
  if(A_prior == "trunc_Normal_OLS_condi"){
    vA   = 1/(Tx/sigma + 1/sigma_A)
    mA   = vA*(Ux/sigma + mu_A/sigma_A)
    sdA  = sqrt(sigma_A)
    lPri = -log(sdA) - log(pnorm((A.max - mu_A)/sdA) -
                           pnorm((A.min - mu_A)/sdA))
    lQ   = mA^2/(2*vA) - mu_A^2/(2*sigma_A)
  }else if(A_prior == "Jeffreys_OLS_condi"){
    vA   = sigma/Tx
    mA   = Ux/Tx
    lPri = -log(A.max - A.min) + 0.5*log(2*pi)
    lQ   = mA^2/(2*vA)
  }else{
    stop("log_h_core: closed-form A-marginal implemented for A_prior = ",
         "'trunc_Normal_OLS_condi' and 'Jeffreys_OLS_condi' only.")
  }
  sA    = sqrt(vA)
  lspan = log(pmax(pnorm((A.max - mA)/sA) - pnorm((A.min - mA)/sA), 1e-300))
  -log(P) + 0.5*log(vA) + lPri + lspan + lQ
}

## log Z1 = log int_0^P h(phi) dphi, per gene, by periodic quadrature
get_logZ1_single = function(Y, t.c, t.s, M, sigma, A_prior,
                            mu_A, sigma_A, A.min, A.max,
                            omega, P, G, NP = 64, chunk = 2000){
  ## A.max may be a length-G vector, hence any()
  if(A_prior == "Jeffreys_OLS_condi" && any(!is.finite(A.max)))
    stop("A_prior = 'Jeffreys_OLS_condi' with A.max = Inf is an IMPROPER ",
         "prior on A, so Z1 -- and hence the Bayes factor between M0 and M1 ",
         "-- is defined only up to an arbitrary constant. Use A_prior = ",
         "'trunc_Normal_OLS_condi' (the prior in the paper, eq. 2) or supply ",
         "a finite A.max.")

  Y = as.matrix(Y)
  if(nrow(Y) != G)
    stop("get_logZ1_single: nrow(Y) = ", nrow(Y), " but G = ", G)
  if(length(t.c) != ncol(Y) || length(t.s) != ncol(Y))
    stop("get_logZ1_single: length(t.c) = ", length(t.c),
         " but ncol(Y) = ", ncol(Y))
  if(length(M) != G)     stop("get_logZ1_single: length(M) = ", length(M))
  if(length(sigma) != G) stop("get_logZ1_single: length(sigma) = ", length(sigma))
  if(any(!is.finite(sigma)) || any(sigma <= 0))
    stop("get_logZ1_single: sigma must be finite and positive; ",
         sum(!is.finite(sigma) | sigma <= 0), " gene(s) violate this ",
         "(min = ", min(sigma), ").")
  if(any(!is.finite(M)))
    stop("get_logZ1_single: M has ", sum(!is.finite(M)), " non-finite value(s).")

  phi.grid = seq(0, P, length.out = NP + 1)[-(NP + 1)]
  dphi     = P / NP
  Cmat = outer(cos(omega*phi.grid), t.c) + outer(sin(omega*phi.grid), t.s)
  tC   = t(Cmat)                                           # n x NP
  Tvec = colSums(tC^2)                                     # NP

  ## genes are processed in chunks to bound the size of the G x NP temporaries
  out = numeric(G)
  starts = seq(1, G, by = chunk)
  for(st in starts){
    id   = st:min(st + chunk - 1L, G)
    Yc   = Y[id, , drop = FALSE] - M[id]
    Umat = Yc %*% tC                                       # |id| x NP
    Tmat = matrix(Tvec, nrow = length(id), ncol = NP, byrow = TRUE)
    ## a per-gene A.max is subset to the chunk
    amx  = if(length(A.max) == 1L) A.max else A.max[id]
    lh   = log_h_core(Tmat, Umat, sigma[id], A_prior,
                      mu_A, sigma_A, A.min, amx, P)
    if(any(!is.finite(lh)))
      stop("get_logZ1_single: log h(phi) is non-finite for ",
           sum(rowSums(!is.finite(lh)) > 0), " gene(s) in rows ",
           min(id), "-", max(id), ". First offending gene: ",
           id[which(rowSums(!is.finite(lh)) > 0)[1]],
           " (sigma = ", sigma[id[which(rowSums(!is.finite(lh)) > 0)[1]]], ").")
    mx  = apply(lh, 1, max)
    out[id] = mx + log(rowSums(exp(lh - mx))) + log(dphi)
    rm(Yc, Umat, Tmat, lh, mx)
  }
  out
}

#' One reversible-jump birth/death sweep over the rhythmicity indicators
#'
#' @description
#' Proposes a model switch for every gene at once (birth for genes at
#' rho = 0, death for genes at rho = 1) and accepts or rejects each
#' independently. The acceptance ratio is assembled from the likelihood term
#' (\code{log.r1}), the prior odds on rhythmicity (\code{log.r2}) and the
#' proposal/prior terms for amplitude and phase (\code{log.r3}).
#'
#' The proposal density factorises as J(A, phi) = J(phi) J(A | phi). The
#' phase term uses the marginal posterior of phi,
#' \code{log h(phi) - log Z1}, with \code{Z1} from \code{get_logZ1_single()}.
#'
#' @param Y G x N expression matrix.
#' @param t.c,t.s Length-N vectors of cos and sin of omega * t.
#' @param N Integer; number of samples.
#' @param t.c.sum,t.s.sum Numeric; column sums of \code{t.c} and \code{t.s}.
#' @param c.t,s.t,cs.t Precomputed cross-products of \code{t.c} and \code{t.s}.
#' @param y.t.c.sum.num,y.t.s.sum.num Length-G vectors of Y-weighted sums.
#' @param AcosPhi,AsinPhi Length-G vectors A * cos and A * sin of omega * phi.
#' @param A,phi Length-G amplitude and acrophase vectors.
#' @param M,sigma Length-G MESOR and residual-variance vectors.
#' @param p.vec Length-G prior probability of rhythmicity.
#' @param rho Length-G current rhythmicity indicator in \{0, 1\}.
#' @param propose.phi,propose.A Logical; include phase and amplitude in the
#'   jump proposal.
#' @param A_prior Character; amplitude prior, as in
#'   \code{CB_MCMC_single_rj_slice}.
#' @param mu_A,sigma_A Numeric; truncated-Normal amplitude prior parameters.
#' @param A.min,A.max Numeric; amplitude truncation bounds. \code{A.max} may
#'   be one value or one per gene.
#' @param omega Numeric; angular frequency 2*pi/P.
#' @param G Integer; number of genes.
#' @param P Numeric; period in hours.
#' @param NP Integer; quadrature nodes for the phi marginal.
#'
#' @return A list with the updated \code{rho}, the normalising constant
#'   \code{logZ1}, the diagnostic \code{Z1_check}, and the components of the
#'   log acceptance ratio.
#'
#' @keywords internal
RJMCMC_single_slice = function(Y, t.c, t.s, N,
                               t.c.sum, t.s.sum, 
                               c.t, s.t, cs.t, 
                               y.t.c.sum.num, y.t.s.sum.num,
                               AcosPhi, AsinPhi, A, phi,  
                               M, sigma, p.vec = rep(0.5, G), rho, 
                               propose.phi = FALSE, propose.A = TRUE,
                               A_prior, 
                               mu_A, sigma_A, A.min, A.max, # A ~ truncated Normal(mu_A, sigma_A) on (A.min, A.max)
                               omega, G, P,
                               NP = 64){   # phi-quadrature nodes for Z1
  ## omega, G and P are formals 27-29; check them here for positional callers
  if(missing(omega) || missing(G) || missing(P))
    stop("RJMCMC_single_slice: omega/G/P are missing. The caller supplied ",
         nargs(), " positional arguments; 29 are required before NP.")

  # current state
  b.c.cur = AcosPhi*rho #G*1 mat
  b.s.cur = AsinPhi*rho 
  a.cosine.cur = cbind(M, b.c.cur, b.s.cur) %*% rbind(1, t.c, t.s)
  log.lik_cur0 =  apply((Y-a.cosine.cur)^2, 1, sum) #without -(2*sigma^2)-1
  b.c.jump = AcosPhi*(1-rho) #G*1 mat
  b.s.jump = AsinPhi*(1-rho)
  a.cosine.jump = cbind(M, b.c.jump, b.s.jump) %*% rbind(1, t.c, t.s)
  log.lik_jump0 =  apply((Y-a.cosine.jump)^2, 1, sum) #without -(2*sigma^2)-1
  
  if(propose.A){
    A.post = get_post_A(Y, t.c, t.s, N,
                        AcosPhi, AsinPhi, A, phi,
                        M, sigma,
                        A_prior,
                        mu_A, sigma_A, A.min, A.max,
                        omega, G)
    ## floor at 1e-300 so a zero density gives a very negative log, not -Inf
    log.A.post = log(pmax(A.post, 1e-300))
    A.prior = get_prior_A(Y, t.c, t.s, N,
                         AcosPhi, AsinPhi, A, phi,
                         M, sigma,
                         A_prior,
                         mu_A, sigma_A, A.min, A.max,
                         omega)
    log.A.prior = log(pmax(A.prior, 1e-300))
  }else{
    log.A.prior = 1; log.A.post = 1
  }
  
  if(propose.phi){
    ## log J(phi) = log h(phi) - log Z1; get_post_A supplies J(A | phi)
    logZ1 = get_logZ1_single(Y, t.c, t.s, M, sigma, A_prior,
                             mu_A, sigma_A, A.min, A.max,
                             omega, P, G, NP)
    log.phi.post = log_h_single(Y, t.c, t.s, M, sigma, phi, A_prior,
                                mu_A, sigma_A, A.min, A.max, omega, P) - logZ1
    ## phi ~ Uniform(0, P), so the log prior density is -log(P)
    log.phi.prior = log(1/P)
  }else{
    log.phi.prior = 1; log.phi.post = 1; logZ1 = rep(NA_real_, G)
  }
  log.r1 = -1/(2*sigma)*(log.lik_jump0-log.lik_cur0)
  log.r2 = log((p.vec+1e-5)/(1-p.vec+1e-5))
  log.r3 = log.A.prior+log.phi.prior-log.A.post-log.phi.post
  log.r = log.r1+(log.r2+log.r3)*((-1)^rho) # sign flips for a death move
  ## a non-finite ratio is a rejection
  log.r[!is.finite(log.r)] = -Inf

  a.jump.log.u = log(runif(G, 0, 1))
  a.rho = ifelse(log.r>a.jump.log.u, 1-rho, rho)

  return(list(rho = a.rho,
              logZ1 = logZ1,
                ## diagnostic only; the sampler does not read it
                Z1_check = log.r1 + log.r3,
              log.r1 = log.r1*((-1)^rho),
              log.r1.SS = (log.lik_jump0-log.lik_cur0)*((-1)^rho),
              log.r3_A = (log.A.prior-log.A.post)*((-1)^rho),
              log.r3_phi = (log.phi.prior-log.phi.post)*((-1)^rho)
  ))
}

# Other functions ---------------------------------------------------------

#' Evaluate expression with fallback save on error
#'
#' @description
#' Wraps an expression in a \code{try()} call.  If the expression throws an
#' error, the current \code{out} object is serialised to \code{save.file}
#' and the error is re-thrown, preserving progress on disk.
#'
#' Used by the time-course sampler in \code{mcmc_time.R}. The single-tissue
#' sampler stops at the failing iteration instead of saving a partial chain.
#'
#' @param expr Expression to evaluate.
#' @param out Object to save on failure.
#' @param save.file Character; file path for the emergency \code{saveRDS}.
#'
#' @return Result of \code{expr} on success.
#'
#' @keywords internal
try_save = function(expr, out, save.file){
  tryCatch(expr, error = function(e) {
    saveRDS(out, save.file)
    stop(e)
  })
}

## Stop at the iteration that first produces an NA in rho, AcosPhi or AsinPhi.
stop_on_na = function(rho, AcosPhi, AsinPhi, iter, stage){
  bad = c(rho = anyNA(rho),
          AcosPhi = anyNA(AcosPhi),
          AsinPhi = anyNA(AsinPhi))
  if(any(bad))
    stop("NA in ", paste(names(bad)[bad], collapse = ", "),
         " at ", stage, " iteration ", iter,
         " (", sum(is.na(rho)), " genes in rho); stopping the chain.",
         call. = FALSE)
  invisible(NULL)
}

adjust.to.2pi = function(x){
  d = x/(2*pi)
  d.abs = floor(abs(d))
  if(d>=0&d<1){
    return(x)
  }else if(d<0){
    return(x+(d.abs+1)*2*pi)
  }else if(d>1){
    return(x-(d.abs)*2*pi)
  }
}

#' Sample from a truncated gamma distribution
#'
#' @description
#' Draws \code{n} samples from a Gamma(shape, rate) distribution truncated
#' to the interval (a, b) using the inverse-CDF method on the log scale.
#'
#' @param n Integer; number of samples.
#' @param a Numeric; lower truncation bound (must be < b).
#' @param b Numeric; upper truncation bound.
#' @param shape Numeric; gamma shape parameter.
#' @param rate Numeric; gamma rate parameter.
#'
#' @return Length-n numeric vector of truncated-gamma samples.
#'
#' @keywords internal
rtruncgamma = function(n, a, b, shape, rate){
  if(a>=b){
    stop( "argument a is greater than or equal to b" )
  }
  G.min = stats::pgamma(a, shape, rate, lower.tail = TRUE, log.p = TRUE);
  G.max = stats::pgamma(b, shape, rate, lower.tail = TRUE, log.p = TRUE);
  u = runif(n)
  
  log_a_new_p = G.min+log(1-u+u*exp(G.max-G.min))
  if(log_a_new_p>G.max){
    log_a_new_p=G.max # guards against rounding above G.max
  }

  new.k = stats::qgamma(log_a_new_p, shape, rate, lower.tail = TRUE, log.p = TRUE)
  if(is.infinite(new.k)|is.na(new.k)){
    print(a)
    print(b)
    print(rate)
    stop("rtruncgamma!!!")
  }
  return(new.k)
}

rtrunclgamma = function(n, a, b, shape, rate){
  if(a>=b){
    stop( "argument a is greater than or equal to b" )
  }
  G.min = VGAM::plgamma(a, shape, scale = 1/rate, log.p = TRUE); 
  G.max = VGAM::plgamma(b, shape, scale = 1/rate, log.p = TRUE);
  u = runif(n)
  
  a.new.p = exp(Rmpfr::mpfr(as.character(G.min), 128))+
    u*(exp(Rmpfr::mpfr(as.character(G.max), 128))-exp(Rmpfr::mpfr(as.character(G.min), 128)))
  
  new.k = stats::qgamma(as.numeric(log(a.new.p)), shape, scale = 1/rate, log.p = TRUE)
  # an infinite quantile falls back to a uniform draw on the lower half
  if(is.infinite(new.k)){
    new.k = runif(1, a, (a+b)/2)
  }
  return(new.k)
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

# Diagnostics -------------------------------------------------------------

#' Post-hoc MCMC convergence diagnostics
#'
#' @description
#' Computes the same convergence diagnostics that
#' \code{CB_MCMC_single_rj_slice()} attaches when called with
#' \code{diagnostics = TRUE} (its default): per-gene RJMCMC acceptance rate,
#' per-gene effective sample size (ESS) for the rhythmicity indicator
#' \code{rho}, and per-gene ESS for phase \code{phi}. It applies to any
#' result, since \code{rho}, \code{phi} and \code{if.accept.rj} are stored
#' whatever the \code{diagnostics} setting.
#'
#' @param mcmc_result List returned by \code{CB_MCMC_single_rj_slice()} (or
#'   compatible), containing \code{$rho}, \code{$phi}, and
#'   \code{$if.accept.rj}.
#' @param P Numeric; circadian period in hours, matching the value used for
#'   the MCMC run (default 24). Needed to convert \code{phi} to radians for
#'   its circular autocorrelation.
#'
#' @details
#' \code{phi} is periodic, so its lag-1 autocorrelation is computed on the
#' circle as the mean cosine of successive phase differences in radians,
#' \code{mean(cos(phi_rad[-1] - phi_rad[-K]))}. The ESS for the binary
#' \code{rho} uses the linear lag-1 autocorrelation.
#'
#' @return A list: \code{acceptance_rate}, \code{ess_rho}, \code{ess_phi}
#'   (all per-gene, named by gene), \code{mean_acceptance_rate},
#'   \code{mean_ess_rho}, \code{mean_ess_phi}, \code{n_samples},
#'   \code{p_rhythmic_posterior} and \code{phi_posterior_median}. Also printed as a summary; warns if mean
#'   ESS for rho falls below 100.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mcmc_out <- CB_MCMC_single_rj_slice(dat, init, diagnostics = FALSE)
#' diag <- mcmc_diagnostics(mcmc_out)
#' }
mcmc_diagnostics = function(mcmc_result, P = 24) {
  if (is.null(mcmc_result$rho) || is.null(mcmc_result$if.accept.rj))
    stop("mcmc_result must contain $rho and $if.accept.rj (as returned by CB_MCMC_single_rj_slice()).")

  # acceptance rate over proposed jumps; -1 marks an iteration with no proposal
  proposed <- mcmc_result$if.accept.rj != -1
  n_proposed <- rowSums(proposed)
  n_accepted <- rowSums(mcmc_result$if.accept.rj == 1, na.rm = TRUE)
  accept_rate <- ifelse(n_proposed > 0, n_accepted / n_proposed, NA_real_)

  # ESS for rho from the lag-1 autocorrelation
  ess_rho <- apply(mcmc_result$rho, 1, function(x) {
    K <- length(x)
    if (K < 4 || var(x) < 1e-10) return(NA_real_)
    r1 <- tryCatch(cor(x[-K], x[-1]), error = function(e) NA_real_)
    if (is.na(r1) || abs(r1) >= 1) return(as.numeric(K))
    K * (1 - r1) / (1 + r1)
  })

  # ESS for phi from the circular lag-1 autocorrelation
  ess_phi <- if (is.null(mcmc_result$phi)) {
    NA_real_
  } else {
    apply(mcmc_result$phi, 1, function(x) {
      K <- length(x)
      if (K < 4) return(NA_real_)
      x_rad <- x * (2 * pi / P)
      r1 <- mean(cos(x_rad[-1] - x_rad[-K]))
      if (is.na(r1) || abs(r1) >= 1) return(as.numeric(K))
      K * (1 - r1) / (1 + r1)
    })
  }

  # circular median of phi over the draws with rho = 1
  phi_posterior_median <- if (is.null(mcmc_result$phi)) {
    NA_real_
  } else {
    vapply(seq_len(nrow(mcmc_result$rho)), function(g) {
      phi_rhy <- mcmc_result$phi[g, mcmc_result$rho[g, ] == 1]
      if (length(phi_rhy) < 4) return(NA_real_)
      circular_median(phi_rhy, P = P)
    }, numeric(1))
  }
  if (is.numeric(phi_posterior_median)) names(phi_posterior_median) <- rownames(mcmc_result$rho)

  diagnostics <- list(
    acceptance_rate      = accept_rate,
    ess_rho              = ess_rho,
    ess_phi              = ess_phi,
    mean_acceptance_rate = mean(accept_rate, na.rm = TRUE),
    mean_ess_rho         = mean(ess_rho, na.rm = TRUE),
    mean_ess_phi         = mean(ess_phi, na.rm = TRUE),
    n_samples            = ncol(mcmc_result$rho),
    p_rhythmic_posterior = rowMeans(mcmc_result$rho),
    phi_posterior_median = phi_posterior_median
  )

  cat("\n=== MCMC Diagnostics ===\n")
  cat("Samples stored:        ", diagnostics$n_samples, "\n")
  cat("Mean acceptance rate:  ", round(diagnostics$mean_acceptance_rate, 3), "\n")
  cat("Mean ESS (rho):        ", round(diagnostics$mean_ess_rho, 1), "\n")
  cat("Mean ESS (phi):        ", round(diagnostics$mean_ess_phi, 1), "\n")
  if (!is.na(diagnostics$mean_ess_rho) && diagnostics$mean_ess_rho < 100)
    warning("Low mean ESS for rho (", round(diagnostics$mean_ess_rho, 1),
            " < 100). Consider increasing `iteration` or decreasing `thin`.")

  diagnostics
}
