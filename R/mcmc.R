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
#'   post-burn-in sample is stored (default 20).
#' @param n.burn Integer; number of burn-in iterations discarded before
#'   storing samples (default 1000).
#' @param seed Integer; random seed for reproducibility (default 15213).
#' @param diagnostics Logical; if \code{TRUE}, compute and print basic
#'   post-run convergence diagnostics (RJMCMC acceptance rate, ESS for rho
#'   and phi) and attach them as \code{$diagnostics} on the returned list,
#'   via \code{mcmc_diagnostics()}. Warns if mean ESS for rho falls below
#'   100 (default \code{TRUE}). Can also be computed later on any MCMC
#'   result by calling \code{mcmc_diagnostics()} directly, even if this was
#'   \code{FALSE} at run time.
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
#'     \item{Z1_gap}{3 x K matrix, not G x K: the median, 99th percentile and
#'       maximum across genes of \code{abs(Z1_check - logZ1)} at each stored
#'       iteration. A birth proposal equal to pi_1 would make this zero; it
#'       is not, and the gap does not shrink with \code{NP_Z1}, so it is not
#'       quadrature error. Diagnostic only.}
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
                                #options: "trunc_Normal",
                                # Jeffreys: p(log(A))\propto 1
                                # Jeffreys_ridge
                                # sq_expo: squared exponential p(A) \propto A*exp(-A^2/2)    
                                # gamma: gamma prior p(A) \ propto A^(alpha-1)exp(-rate*A)
                                mu_A = 1, sigma_A = 10^2, A.min = 0, # this applies to A~truncated normal
                                A.max = NULL,
                                rj.phi = TRUE, rj.A = TRUE,
                                mu_M = 0, sigma_M = 10^2, 
                                sigma_prior_v = 4, sigma_prior_s = 1,
                                NP_Z1 = 64   # phi-quadrature nodes for Z1
){
  
  # Data.list=dat.input; Init.value=a.init;
  # P = 24; A.min = 0;
  # iteration = n.iter; thin = n.thin; n.burn=n.burn;
  # seed = 15213;
  # p_rhythmic=0.5; rj.p.stay = 0.5;
  # A_prior= "Jeffreys";
  # mu_A = 1; sigma_A = 10^2; # this applies to A~truncated normal
  # rj.phi = TRUE; rj.A = TRUE;
  # mu_M = 0; sigma_M = 10^2;
  # sigma_prior_v = 4; sigma_prior_s = 1;
  # delta_MH_phi = pi/2;
  # save.file = paste0(out.dir, "/", out.name, "/",
  #                    "/tempSaveXXX", a.ind, ".rds")

  #observations
  omega = 2*pi/P
  Y = as.matrix(Data.list[[1]])
  ## Default bound: half the observed range of each gene. Pass A.max = Inf for
  ## the plain one-sided truncated-Normal amplitude prior on [A.min, Inf),
  ## which is proper for A_prior = "trunc_Normal_OLS_condi" but improper for
  ## "Jeffreys_OLS_condi" (get_logZ1_single refuses that combination).
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
  #Below are for slice sampling
  X = cbind(t.c, t.s)
  #OLS solution for the Jeffreys prior
  # beta_hat =  as.matrix(Y)%*%t(solve(t(X)%*%X)%*%t(X))
  XtX = t(X)%*%X
  beta_hat0 =  t(solve(XtX)%*%t(X))
  beta_cov0 = solve(XtX) #time sigma in each iteration
  beta_cov0_11 = beta_cov0[1, 1]
  beta_cov0_22 = beta_cov0[2, 2]
  beta_cov0_12 = beta_cov0[1, 2]
  beta_cov0_rho = beta_cov0_12/sqrt(beta_cov0_11*beta_cov0_22)
  #calculate the coefficients before(beta-mu1), need to time sigma 
  beta_mean0_1 = beta_cov0_12/beta_cov0_22
  beta_mean0_2 = beta_cov0_12/beta_cov0_11
  
  #test use
  # Y.list = Y; t.c.list = t.c; t.s.list = t.s; N.vec = N;
  # A.mat = A; phi.mat = phi; sigma.mat = sigma; M.mat = M; rho.mat = rho; p.mat = p;
  
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
    #update M    
    M = update_M_single(Y, t.c, t.s, N,
                        AcosPhi, AsinPhi, sigma, rho,
                        sigma_M, mu_M, omega, G)
    
    #udpate A
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
  ## The two G x K diagnostic matrices cost 162 MB per tissue at paper scale
  ## and nothing downstream reads them gene by gene, so only the spread of
  ## the discrepancy is kept: three numbers per iteration instead of G.
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
    # iter = iter+1
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
      if.accept.rj = cbind(if.accept.rj, rho.res$rho!=rho) 
      #if the new rho is not the same as the old rho, it is accepted
      log.r1.store = cbind(log.r1.store, rho.res$log.r1)
      log.r1.SS.store = cbind(log.r1.SS.store, rho.res$log.r1.SS) 
      log.r3_A.store = cbind(log.r3_A.store, rho.res$log.r3_A)
      log.r3_phi.store = cbind(log.r3_phi.store, rho.res$log.r3_phi)
      Z1_gap.store = cbind(Z1_gap.store, Z1_gap(rho.res))
      rho = rho.res$rho
      rho.store = cbind(rho.store, rho)
    }    
    time1 = Sys.time()
    # samp.time$rho[iter]=time1-time0
    
    # within model move -------------------------------------------------------
    
    #update M    
    M = update_M_single(Y, t.c, t.s, N,
                        AcosPhi, AsinPhi, sigma, rho,
                        sigma_M, mu_M, omega, G)
    M.store = cbind(M.store, M)
    # print("Finished M")
    time2 = Sys.time()
    # samp.time$M[iter]=time2-time1
    
    #udpate A
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
    # samp.time$phi[iter]=time4-time2
    
    sigma = update_sigma_single(Y, t.c, t.s, N, AcosPhi, AsinPhi, M, rho,
                                sigma_prior_v, sigma_prior_s, omega, G)
    sigma.store = cbind(sigma.store, sigma)
    # print("Finished sigma")
    time5 = Sys.time()
    # samp.time$sigma[iter]=time5-time4
    
    save = list(rho = rho.store,
                M = M.store,
                AcosPhi = AcosPhi.store,
                AsinPhi = AsinPhi.store, 
                A = A.store, 
                phi = phi.store,
                sigma = sigma.store,
                Z1_gap = Z1_gap.store,
                #log.r1 = log.r1.store,
                #log.r3_A = log.r3_A.store,
                #log.r3_phi = log.r3_phi.store,
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
#Gibbs sampling

#should update rho, A, and phi together with a reversible jump procedure
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
                           #Y.list, t.c.list, t.s.list, N.vec, 
                           AcosPhi, AsinPhi, sigma, rho, 
                           # A.mat, phi.mat, sigma.mat, rho.mat, 
                           sigma_M, mu_M, omega, G){#output a G*J matrix for M samples
  # Y.list = Y; t.c.list = t.c; t.s.list = t.s; N.vec = N; A.mat = A; phi.mat = phi; sigma.mat = sigma; rho.mat = rho
  # t.c=t_t_p.c; t.s = t_t_p.s
  a.var = 1/(N/sigma+1/sigma_M)
  b.c = AcosPhi*rho #G*1 mat
  b.s = AsinPhi*rho
  a.cosine = cbind(b.c, b.s) %*% rbind(t.c, t.s)
  a.mean = (apply((Y-a.cosine), 1, sum)/sigma+mu_M/sigma_M)*a.var
  # xM <<- list(a.mean = a.mean, a.var = a.var)
  # stopifnot("M problem" = all(is.finite(a.mean))&all(a.var<10000))
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
                               # Y.list, t.c.list, t.s.list, N.vec,
                               AcosPhi, AsinPhi, M.vec, rho, 
                               # A.mat, phi.mat, M.mat, rho.mat, 
                               sigma_prior_v, sigma_prior_s, omega, G){
  # Y.list = Y; t.c.list = t.c; t.s.list = t.s; N.vec = N; A.mat = A; phi.mat = phi; M.mat = M; rho.mat = rho
  b.c = AcosPhi*rho #G*1 mat
  b.s = AsinPhi*rho
  a.cosine = cbind(M.vec, b.c, b.s) %*% rbind(1, t.c, t.s)
  vs2 = sigma_prior_v*sigma_prior_s
  sigma_post_vs2 = (apply((Y-a.cosine)^2, 1, sum)+vs2)/2
  # print(paste0("SS 1:", round(mean(sigma_post_vs2[1:200]), 2), ", 0:", round(mean(sigma_post_vs2[201:400]), 2)))
  # xsigma <<- list(shape = sigma_prior_v/2+a.n/2, rate = sigma_post_vs2)
  # stopifnot("sigma problem" = all(xsigma$shape<10000)&all(is.finite(xsigma$rate)))
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
                              mu_A, sigma_A  # this applies to A~truncated normal(mu_A, sigma_A)
){
  beta1 = AcosPhi
  beta2 = AsinPhi
  Y_new = Y-M
  A = sqrt(beta1^2+beta2^2)
  #for the truncation model only
  f_upper_trunc = function(A0, u0){
    1/A0*exp(-1/sigma_A*(A0-mu_A)^2)-u0
  }
  
  if(A_prior=="Jeffreys_OLS_condi"){
    
# A_prior=="Jeffreys_OLS_condi" -------------------------------------------
    
    #OLS solution for the Jeffreys prior
    beta_hat = Y_new%*%beta_hat0
    u = runif(G, 1e-5, 1/A) #sqrt(beta1^2+beta2^2) is just A
    u = ifelse(u>1e3, 1e3, u) #avoid numerical issue for later
    beta_cov_11 = beta_cov0_11*sigma
    beta_cov_22 = beta_cov0_22*sigma
    max_beta1 = sqrt(1/u^2-beta2^2)
    beta1_mean = beta_hat[, 1]+beta_mean0_1*(beta2-beta_hat[, 2])
    beta1_var = beta_cov_11*(1-beta_cov0_rho^2)
    beta1_new = truncnorm::rtruncnorm(G, -1*max_beta1, max_beta1, 
                                      beta1_mean, sqrt(beta1_var)) 
    # print(paste0("beta1_new 1:", round(mean(beta1_new[1:200]), 2), ", 0:", round(mean(beta1_new[201:400]), 2)))
    # print(paste0("beta1_mean sd 1:", round(sd(beta1_mean[1:200]), 2), ", 0:", round(sd(beta1_mean[201:400]), 2)))
    # print(paste0("beta1_new sd 1:", round(sd(beta1_new[1:200]), 2), ", 0:", round(sd(beta1_new[201:400]), 2)))
    max_beta2 = sqrt(1/u^2-beta1_new^2)
    beta2_mean = beta_hat[, 2]+beta_mean0_2*(beta1_new-beta_hat[, 1])
    beta2_var = beta_cov_22*(1-beta_cov0_rho^2)
    beta2_new = truncnorm::rtruncnorm(G, -1*max_beta2, max_beta2, 
                                      beta2_mean, sqrt(beta2_var)) 
    # print(paste0("beta2_new 1:", round(mean(beta2_new[1:200]), 2), ", 0:", round(mean(beta2_new[201:400]), 2)))
    # print(paste0("beta2_mean sd 1:", round(sd(beta2_mean[1:200]), 2), ", 0:", round(sd(beta2_mean[201:400]), 2)))
    # print(paste0("beta2_new sd 1:", round(sd(beta2_new[1:200]), 2), ", 0:", round(sd(beta2_new[201:400]), 2)))

  }else if(A_prior=="Jeffreys_ridge_condi"){
    
    # A_prior=="Jeffreys_ridge_condi" -----------------------------------------
    
    #ridge regression solution 
    # u = runif(G, 0, 1/A)
    u = runif(G, 1e-5, 1/A) #avoid extreme small vals 
    u = ifelse(u>1e3, 1e3, u)
    get_beta_hat0 = sapply(1:G, function(g){
      # for(g in 1:G){
      XtXinv = solve(XtX+1/u[g]^2)
      beta_hat = Y_new[g, ]%*%t(XtXinv%*%t(X))
      beta_cov0 = XtXinv%*%XtX%*%XtXinv
      beta_cov_11 = beta_cov0[1, 1]*sigma[g]
      beta_cov_22 = beta_cov0[2, 2]*sigma[g]
      beta_cov_rho = beta_cov0[1, 2]/sqrt(beta_cov_11*beta_cov_22)
      # }
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
    
    #OLS solution 
    beta_hat = Y_new%*%beta_hat0
    u = runif(G, 1e-3, f_upper_trunc(A, 0))
    u = ifelse(u>1e3, 1e3, u)
    beta_cov_11 = beta_cov0_11*sigma
    beta_cov_22 = beta_cov0_22*sigma

    #A_max = sapply(1:G, function(g){
      # print(paste0(a, "_", g))
     # x=tryCatch({
      #  uniroot(f_upper_trunc, c(0, 1e3), u0=u[g], tol=1e-6)$root
      #}, error = function(e) {
      #  cat("An error occurred: ", conditionMessage(e), "\n")
      #  return(conditionMessage(e))
      #})
      #if(grepl("not of opposite sign", x)){
      #  if(f_upper_trunc(1e4)>0){
      #    return(A.max[g])
      #  }else{
      #    return(min(1e-3, A.max[g]))
      #  }
      #}else{
      #  return(min(x, A.max[g]))
      #}
    #})
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
    ## When the slice is narrow both tail probabilities underflow to 0 and
    ## p1/(p1+p2) is 0/0; rbinom then draws NA and the NA spreads through the
    ## whole sweep. Fall back to an even split, which is what the two equal
    ## (zero) masses imply.
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
    #below are the same as A_prior=="Jeffreys"
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
                      mu_A, sigma_A, A.min, A.max, # this applies to A~truncated normal(mu_A, sigma_A)I(A>A_min)
                      omega, G){
  b.c = cos(omega*phi) #G*1 mat
  b.s = sin(omega*phi)
  t.phi = cbind(b.c, b.s) %*% rbind(t.c, t.s)
  T.x = apply(t.phi^2, 1, sum)
  Y.minu.M = Y - as.vector(M)
  U.x = apply(t.phi*Y.minu.M, 1, sum)
  
  if(A_prior=="trunc_Normal_OLS_condi"){
    #Here we are able to simplify the calculation by using AcosPhi/A and AsinPhi/A. Work later

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
                      mu_A, sigma_A, A.min, A.max, # this applies to A~truncated normal(mu_A, sigma_A)I(A>A_min)
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
  
  -1/2*diff1/a.sigma.vec#+kappa.vec*cos(a.phi.vec*omega-theta.vec)
}

update_phi_single = function(a.Y, a.t.c, a.t.s, a.N,
                             A.vec, M.vec, phi.vec, sigma.vec, rho.vec,
                             delta_MH_phi, omega, G, P){
  #using metropolis algorithm
  old.log.phi.post = get_log_phi_post_single(a.Y, a.t.c, a.t.s, a.N,
                                             M.vec, A.vec, phi.vec, sigma.vec, #rho.mat[, a],
                                             omega)
  # new.phi = Rfast::rvonmises(G, phi.mat[, a]/P*2*pi, delta_MH_phi)/(2*pi)*P
  # new.phi = CircStats::rvm(G, phi.mat[, a]/P*2*pi, delta_MH_phi)/(2*pi)*P
  new.phi = sapply(1:G, function(g){
    CircStats::rvm(1, phi.vec[g]/P*2*pi, delta_MH_phi)/(2*pi)*P})
  
  new.log.phi.post = get_log_phi_post_single(a.Y, a.t.c, a.t.s, a.N,
                                             M.vec, A.vec, new.phi, sigma.vec, #rho.mat[, a], 
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

## NOTE: no longer called by RJMCMC_single_slice. It returns the CONDITIONAL
## p(phi | A); the RJ ratio needs the MARGINAL p(phi), now supplied by
## log_h_single() / get_logZ1_single(). Retained for reference.
get_log_phi_post_full_single = function(a.Y.vec, a.t.c.vec, a.t.s.vec, 
                                        a.t.c.sum.num, a.t.s.sum.num,
                                        c.t, s.t, cs.t,
                                        a.y.t.c.sum.num, a.y.t.s.sum.num,
                                        a.sigma.num, a.A.num, a.phi.num, a.M.num, 
                                        omega, P){
  
  # g = 1
  # get_log_phi_post_full_single(Y[g, ], t.c, t.s, t.c.sum, t.s.sum, c.t, s.t, cs.t, y.t.c.sum.num[g], y.t.s.sum.num[g], sigma[g], A[g], M[g], phi[g], omega, P)
  # a.Y.vec = Y[g, ]; a.t.c.vec = t.c; a.t.s.vec = t.s;a.t.c.sum.num=t.c.sum; a.t.s.sum.num =t.s.sum ;a.y.t.c.sum.num =y.t.c.sum.num[g];a.y.t.s.sum.num =y.t.s.sum.num[g];a.sigma.num = sigma[g]; a.A.num = A[g]; a.M.num = M[g]; a.phi.num = phi[g];
  # a.t.c.vec = t.c.vec; a.t.s.vec = t.s.vec; a.Y.vec = Y.mat[g, ]; a.sigma.num = sigma.vec[g]; a.A.num = A.vec[g]; a.M.num = M.vec[g]; a.phi.num = phi.vec[g]; a.kappa.num = kappa; a.theta.num = theta
  # a.t.c.vec = a.t.c; a.t.s.vec = a.t.s; a.Y.vec = a.Y.list[g, ]; a.sigma.num = sigma.mat[g, a]; a.A.num = A.mat[g, a]; a.M.num = M.mat[g, a]; a.phi.num = phi.mat[g, a]; a.kappa.num = kappa[g]; a.theta.num = theta[g]
  
  c.t_Y_M = a.y.t.c.sum.num-a.M.num*a.t.c.sum.num
  s.t_Y_M = a.y.t.s.sum.num-a.M.num*a.t.s.sum.num
  a.r = a.A.num/a.sigma.num
  r2 = a.A.num^2/a.sigma.num
  phi.c = cos(omega*a.phi.num)
  phi.s = sin(omega*a.phi.num)
  
  exp.scale=0 #will only be changed under very extreme circumstances
  
  log.a.phi = -r2/2*c.t*phi.c^2-r2/2*s.t*phi.s^2-
    r2*cs.t*phi.s*phi.c+
    a.r*c.t_Y_M*phi.c+
    a.r*s.t_Y_M*phi.s
  a1 = log.a.phi
  b1.scale=0

  phi.fun = function(phi, scale.val){
    exp((-r2/2*c.t*cos(omega*phi)^2-r2/2*s.t*sin(omega*phi)^2-
          r2*cs.t*sin(omega*phi)*cos(omega*phi)+
          a.r*c.t_Y_M*cos(omega*phi)+
          a.r*s.t_Y_M*sin(omega*phi))+scale.val
        )
  }
  # tryCatch(integrate(phi.fun, lower = 0, upper = P), error=function(e) e)
  
  phi.fun.log = function(phi){
    -r2/2*c.t*cos(omega*phi)^2-r2/2*s.t*sin(omega*phi)^2-
      r2*cs.t*sin(omega*phi)*cos(omega*phi)+
      a.r*c.t_Y_M*cos(omega*phi)+
      a.r*s.t_Y_M*sin(omega*phi)
  }
  
  b1.exp = tryCatch(integrate(phi.fun,
                              lower = 0,
                              upper = P, b1.scale),
                    error=function(e) e)
  
  if(("error" %in% class(b1.exp))|b1.exp[[1]]==0){
    
    phi.fun.log.dx = function(phi){
      r2/2*(c.t-s.t)*sin(2*omega*phi)-
        r2*cs.t*cos(2*omega*phi)-
        a.r*c.t_Y_M*sin(omega*phi)+
        a.r*s.t_Y_M*cos(omega*phi)
    }

    initial_guesses = seq(a.phi.num, a.phi.num+P, by=2)
    initial_guesses[initial_guesses>P]=initial_guesses[initial_guesses>P]-P
    initial_guesses[initial_guesses<0]=initial_guesses[initial_guesses<0]+P
    find_root = function(guess) {
      tryCatch(
        rootSolve::multiroot(phi.fun.log.dx,
                             start = guess, maxiter = 100)$root,
        error = function(e) NA
      )
    }
    
    
    # Apply find_root to each initial guess
    roots = sapply(initial_guesses, find_root)
    # Remove NA values and duplicates
    roots = unique(na.omit(roots))
    roots = sort(cleanRoot(roots, range = c(0, P)))
    
    if(length(roots)<2){
      initial_guesses = seq(a.phi.num, a.phi.num+P, by=1)
      
      initial_guesses[initial_guesses>P]=initial_guesses[initial_guesses>P]-P
      initial_guesses[initial_guesses<0]=initial_guesses[initial_guesses<0]+P
      find_root = function(guess) {
        tryCatch(
          rootSolve::multiroot(phi.fun.log.dx,
                               start = guess, maxiter = 1000)$root,
          error = function(e) NA
        )
      }
      # Apply find_root to each initial guess
      roots = sapply(initial_guesses, find_root)
      # Remove NA values and duplicates
      roots = unique(na.omit(roots))
      roots = sort(cleanRoot(roots, range = c(0, P)))
    }
    
    roots.type = rootsMinMax(phi.fun.log.dx, roots, range = c(0, P))
    roots.max = roots[roots.type=="max"]
    roots.min = roots[roots.type=="min"]
    # max.vals = phi.fun(roots.max)
    # min.vals = phi.fun(roots.min)
    max.log.vals = phi.fun.log(roots.max)
    min.log.vals = phi.fun.log(roots.min)
    max.global.log = max(max.log.vals)
    min.global.log = min(min.log.vals)
    range.global.log = max.global.log-min.global.log
    
    scale.target = c(seq(400, 2300, by = 100), seq(0, -2300, by=-100))
    flag=0
    i=1
    
    while(flag==0&i<=length(scale.target)){
      b1.scale = scale.target[i]-max.global.log
      b1.exp = tryCatch(integrate(phi.fun,
                                  lower = 0,
                                  upper = P, b1.scale),
                        error=function(e) e)
      
      if(("error" %in% class(b1.exp))|b1.exp[[1]]==0){
        a.range.size = P/2
        b1.exp = tryCatch(integrate(phi.fun,
                                    lower = roots.max[1]-a.range.size,
                                    upper = roots.max[1]+a.range.size, b1.scale),
                          error=function(e) e)
      }else{flag=1}
      
      if(("error" %in% class(b1.exp))|b1.exp[[1]]==0){
        if(length(roots.min)==2){
          #two interval separation
          b11.exp = tryCatch(integrate(phi.fun, lower = roots.min[1], upper = roots.min[2], b1.scale), error=function(e) e)
          b12.exp = tryCatch(integrate(phi.fun, lower = roots.min[2], upper = roots.min[1]+P, b1.scale), error=function(e) e)
          b1.exp = tryCatch(list(exp=b11.exp[[1]]+b12.exp[[1]]) , error = function(e) e)
        }
      }else{flag=1}
      
      for(eps in 0.1^(c(0:50, seq(-P/2.5/10, -0.001, by=0.02)))){
        if(("error" %in% class(b1.exp))|b1.exp[[1]]==0){
          if(length(roots.max)==2){
            #two interval separation
            b11.exp = tryCatch(integrate(phi.fun, lower = roots.max[1]-eps, upper = roots.max[1]+eps, b1.scale), error=function(e) e)
            b12.exp = tryCatch(integrate(phi.fun, lower = roots.max[2]-eps, upper = roots.max[2]+eps, b1.scale), error=function(e) e)
            b1.exp = tryCatch(list(exp=b11.exp[[1]]+b12.exp[[1]]) , error = function(e) e)
          }else if(length(roots.max)==1){
            b1.exp = tryCatch(integrate(phi.fun, lower = roots.max[1]-eps, upper = roots.max[1]+eps, b1.scale), error=function(e) e)
          }
        }else{flag=1}
      }

      i=i+1
    }
  }
  
  #this will be a very crude approximation
  if(("error" %in% class(b1.exp))|b1.exp[[1]]==0){
    target.range = 1e4
    exp.scale = round(log10(range.global.log/target.range))
    log.a.phi = -r2/2*c.t*phi.c^2-r2/2*s.t*phi.s^2-
      r2*cs.t*phi.s*phi.c+
      a.r*c.t_Y_M*phi.c+
      a.r*s.t_Y_M*phi.s
    a1 = log.a.phi/10^exp.scale
    
    b1.scale=0
    phi.fun = function(phi, scale.val){
      exp((-r2/2*c.t*cos(omega*phi)^2-r2/2*s.t*sin(omega*phi)^2-
             r2*cs.t*sin(omega*phi)*cos(omega*phi)+
             a.r*c.t_Y_M*cos(omega*phi)+
             a.r*s.t_Y_M*sin(omega*phi))/10^exp.scale+scale.val
      )
    }

    phi.fun.log = function(phi){
      (-r2/2*c.t*cos(omega*phi)^2-r2/2*s.t*sin(omega*phi)^2-
        r2*cs.t*sin(omega*phi)*cos(omega*phi)+
        a.r*c.t_Y_M*cos(omega*phi)+
        a.r*s.t_Y_M*sin(omega*phi))/10^exp.scale
    }
    
    b1.exp = tryCatch(integrate(phi.fun,
                                lower = 0,
                                upper = P, b1.scale),
                      error=function(e) e)
    
    if(("error" %in% class(b1.exp))|b1.exp[[1]]==0){
      
      max.log.vals = phi.fun.log(roots.max)
      min.log.vals = phi.fun.log(roots.min)
      max.global.log = max(max.log.vals)
      min.global.log = min(min.log.vals)
      range.global.log = max.global.log-min.global.log
      
      scale.target = c(seq(400, 2300, by = 100), seq(0, -2300, by=-100))
      flag=0
      i=1
      
      while(flag==0&i<=length(scale.target)){
        b1.scale = scale.target[i]-max.global.log
        b1.exp = tryCatch(integrate(phi.fun,
                                    lower = 0,
                                    upper = P, b1.scale),
                          error=function(e) e)
        
        if(("error" %in% class(b1.exp))|b1.exp[[1]]==0){
          a.range.size = P/2
          b1.exp = tryCatch(integrate(phi.fun,
                                      lower = roots.max[1]-a.range.size,
                                      upper = roots.max[1]+a.range.size, b1.scale),
                            error=function(e) e)
        }else{flag=1}
        
        if(("error" %in% class(b1.exp))|b1.exp[[1]]==0){
          if(length(roots.min)==2){
            #two interval separation
            b11.exp = tryCatch(integrate(phi.fun, lower = roots.min[1], upper = roots.min[2], b1.scale), error=function(e) e)
            b12.exp = tryCatch(integrate(phi.fun, lower = roots.min[2], upper = roots.min[1]+P, b1.scale), error=function(e) e)
            b1.exp = tryCatch(list(exp=b11.exp[[1]]+b12.exp[[1]]) , error = function(e) e)
          }
        }else{flag=1}
        
        for(eps in 0.1^(c(0:50, seq(-P/2.5/10, -0.001, by=0.02)))){
          if(("error" %in% class(b1.exp))|b1.exp[[1]]==0){
            if(length(roots.max)==2){
              #two interval separation
              b11.exp = tryCatch(integrate(phi.fun, lower = roots.max[1]-eps, upper = roots.max[1]+eps, b1.scale), error=function(e) e)
              b12.exp = tryCatch(integrate(phi.fun, lower = roots.max[2]-eps, upper = roots.max[2]+eps, b1.scale), error=function(e) e)
              b1.exp = tryCatch(list(exp=b11.exp[[1]]+b12.exp[[1]]) , error = function(e) e)
            }else if(length(roots.max)==1){
              b1.exp = tryCatch(integrate(phi.fun, lower = roots.max[1]-eps, upper = roots.max[1]+eps, b1.scale), error=function(e) e)
            }
          }else{flag=1}
        }
        
        i=i+1
      }
    }
        
  }else{
    b1 = log(b1.exp[[1]])
    out = a1+b1.scale-b1
  }
  
  if(("error" %in% class(b1.exp))|b1.exp[[1]]==0){
    # out = 0 #when there is still no return (i reaches 20), it is probabily because of a very large integral, so we just give out = 0.
    stop("out is NULL")
  }else{
    b1 = log(b1.exp[[1]])
    out = a1+b1.scale-b1+exp.scale
  }
  return(out)
}

## ===========================================================================
## MARGINAL posterior of phi under M1, and its normalizing constant Z1.
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
## and the MARGINAL posterior of phi is  J(phi) = h(phi) / Z1  with
## Z1 = int_0^P h(phi) dphi.
##
## Note m_A^2/(2 v_A) is the profiled quadratic form: as sigma_A^2 -> Inf it
## becomes U(phi)^2 / (2 sigma^2_g T(phi)), i.e. the Lomb-Scargle periodogram.
## The marginal therefore favours phases whose best-fitting amplitude explains
## the data well -- unlike the CONDITIONAL p(phi | A), which holds A fixed.
##
## h(phi) is smooth and PERIODIC on [0, P], so the uniform-grid trapezoid rule
## converges geometrically: NP = 32 already reaches machine precision. Only Z1
## uses the grid; log h is evaluated exactly at the current phi.
## ===========================================================================

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
  ## any(): A.max may be a length-G vector (per-gene bound). `&&` with a
  ## length > 1 right-hand side is an ERROR in R >= 4.3.
  if(A_prior == "Jeffreys_OLS_condi" && any(!is.finite(A.max)))
    stop("A_prior = 'Jeffreys_OLS_condi' with A.max = Inf is an IMPROPER ",
         "prior on A, so Z1 -- and hence the Bayes factor between M0 and M1 ",
         "-- is defined only up to an arbitrary constant. Use A_prior = ",
         "'trunc_Normal_OLS_condi' (the prior in the paper, eq. 2) or supply ",
         "a finite A.max.")

  Y = as.matrix(Y)
  ## informative failures rather than a silent NaN cascade
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

  ## Chunk over genes. The previous version built ~20 temporary G x NP
  ## matrices at once; at G ~ 15000, NP = 64 that is >100 MB of transients per
  ## fork, and with several forks under mclapply the children can be killed,
  ## which surfaces only as "all scheduled cores encountered errors".
  out = numeric(G)
  starts = seq(1, G, by = chunk)
  for(st in starts){
    id   = st:min(st + chunk - 1L, G)
    Yc   = Y[id, , drop = FALSE] - M[id]
    Umat = Yc %*% tC                                       # |id| x NP
    Tmat = matrix(Tvec, nrow = length(id), ncol = NP, byrow = TRUE)
    ## A.max may be per-gene; it MUST be subset with the chunk or it recycles
    ## column-major against a smaller matrix and silently gives the wrong
    ## prior span for all but the first phi column.
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
#' Proposes a model switch for every gene at once -- birth for genes currently
#' at rho = 0, death for genes at rho = 1 -- and accepts or rejects each
#' independently. The acceptance ratio is assembled from the likelihood term
#' (\code{log.r1}), the prior odds on rhythmicity (\code{log.r2}) and the
#' proposal/prior terms for amplitude and phase (\code{log.r3}).
#'
#' The phase half of \code{log.r3} uses the MARGINAL posterior of phi,
#' \code{log h(phi) - log Z1}, with \code{Z1} from
#' \code{get_logZ1_single()}. Pairing a conditional in A with a conditional in
#' phi would not give a joint density in (A, phi); the factorization has to be
#' J(A, phi) = J(phi) J(A | phi).
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
                               mu_A, sigma_A, A.min, A.max, # this applies to A~truncated normal(mu_A, sigma_A)I(A>A_min)
                               omega, G, P,
                               NP = 64){   # phi-quadrature nodes for Z1
  ## Positional-call guard: omega/G/P are formals 27-29 and BayCT.R supplies
  ## them positionally, so a caller with the wrong arity fails here with a
  ## message naming the problem rather than three frames deep in get_post_A.
  if(missing(omega) || missing(G) || missing(P))
    stop("RJMCMC_single_slice: omega/G/P are missing. The caller supplied ",
         nargs(), " positional arguments; 29 are required before NP. ",
         "Check the call site in BayCT.R.")
                      # c.t=c.t.phi; s.t=s.t.phi; cs.t= cs.t.phi;
                      # y.t.c.sum.num=y.t.c.sum; y.t.s.sum.num=y.t.s.sum;
                      # p.vec = p_rhythmic

# t.c=t_t_p.c; t.s=t_t_p.s;
#                       t.c.sum=t_t_p.c.sum; t.s.sum=t_t_p.s.sum;
#                       c.t=c.t_t_p.phi; s.t=s.t_t_p.phi; cs.t=cs.t_t_p.phi;
#                       y.t.c.sum.num=y.t_t_p.c.sum; y.t.s.sum.num=y.t_t_p.s.sum;

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
    # A.vec = sqrt(AcosPhi^2+AsinPhi^2)
    A.post = get_post_A(Y, t.c, t.s, N, 
                        AcosPhi, AsinPhi, A, phi, 
                        M, sigma, 
                        A_prior, 
                        mu_A, sigma_A, A.min, A.max, # this applies to A~truncated normal(mu_A, sigma_A)I(A>A_min)
                        omega, G)
    ## ifelse(x == 0, 1e-7, .) substituted a density of ~1 where -Inf was meant
    log.A.post = log(pmax(A.post, 1e-300))
    A.prior = get_prior_A(Y, t.c, t.s, N, 
                         AcosPhi, AsinPhi, A, phi, 
                         M, sigma, 
                         A_prior, 
                         mu_A, sigma_A, A.min, A.max, # this applies to A~truncated normal(mu_A, sigma_A)I(A>A_min)
                         omega)
    log.A.prior = log(pmax(A.prior, 1e-300))
  }else{
    log.A.prior = 1; log.A.post = 1
  }
  
  if(propose.phi){
    ## MARGINAL posterior of phi, log J(phi) = log h(phi) - log Z1.
    ## Previously this slot held get_log_phi_post_full_single, the CONDITIONAL
    ## p(phi | A). Pairing a conditional in A with a conditional in phi does
    ## not give a joint density in (A, phi): the product integrates to roughly
    ## 1.2-1.4, not 1. A valid factorization needs one MARGINAL and one
    ## conditional, J(A, phi) = J(phi) J(A | phi), which is what this is.
    logZ1 = get_logZ1_single(Y, t.c, t.s, M, sigma, A_prior,
                             mu_A, sigma_A, A.min, A.max,
                             omega, P, G, NP)
    log.phi.post = log_h_single(Y, t.c, t.s, M, sigma, phi, A_prior,
                                mu_A, sigma_A, A.min, A.max, omega, P) - logZ1
    ## phi ~ Uniform(0, P) since kappa_0 = 0, so the density is 1/P and its
    ## log is -log(P) = -3.178. The previous line supplied 1/P = +0.042, the
    ## density itself rather than its log, inflating r by exp(3.22) = 25x in
    ## the M0 -> M1 direction and deflating it 25x in reverse.
    log.phi.prior = log(1/P)
  }else{
    log.phi.prior = 1; log.phi.post = 1; logZ1 = rep(NA_real_, G)
  }
  log.r1 = -1/(2*sigma)*(log.lik_jump0-log.lik_cur0)
  log.r2 = log((p.vec+1e-5)/(1-p.vec+1e-5))#*((-1)^rho) #already taken care below
  log.r3 = log.A.prior+log.phi.prior-log.A.post-log.phi.post
  log.r = log.r1+(log.r2+log.r3)*((-1)^rho) #if rho.mat = 1, then reverse
  ## Without this, a single non-finite log.r makes ifelse() return NA, rho
  ## becomes NA, and the NA propagates into update_M_single and kills the sweep
  ## at iteration 1. Treat a non-finite ratio as "reject".
  log.r[!is.finite(log.r)] = -Inf

  a.jump.log.u = log(runif(G, 0, 1))
  a.rho = ifelse(log.r>a.jump.log.u, 1-rho, rho)

  return(list(rho = a.rho,
              logZ1 = logZ1,
                ## log.r1 + log.r3 would equal log Z1 identically if the
                ## birth proposal were exactly pi_1. It does not: the gap is
                ## unchanged (to four decimals) at NP = 32, 64, 128 and 256,
                ## so it is not phi-grid error. The single-tissue sampler
                ## draws a proposal where the multi-tissue one draws none,
                ## and that is the likeliest source. Kept as a diagnostic;
                ## the sampler's decisions do not read it.
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
#' The single-tissue sampler no longer uses this -- it stops at the iteration
#' that goes wrong instead of saving a partial chain -- but the time-course
#' sampler in \code{mcmc_time.R} still does.
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

## An NA anywhere in rho or in the (AcosPhi, AsinPhi) pair spreads through the
## rest of the sweep and comes back as a chain of NAs hours later. Stop at the
## iteration that produced it instead, and say which state went bad.
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
    log_a_new_p=G.max #By calculation, log_a_new_p will never be larger than G.max, but there is numerical issue...
  }
  
  # a.new.p = exp(Rmpfr::mpfr(as.character(G.min), 128))+
  #   u*(exp(Rmpfr::mpfr(as.character(G.max), 128))-exp(Rmpfr::mpfr(as.character(G.min), 128)))
  # u = runif(n)
  # (log_a_new_p = G.min+log(1-u+u*exp(G.max-G.min)))
  # (log_a_new_p<G.max)
  # stats::qgamma(log_a_new_p, shape, rate, lower.tail = TRUE, log.p = TRUE)
  new.k = stats::qgamma(log_a_new_p, shape, rate, lower.tail = TRUE, log.p = TRUE) #still gets inf
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
  
  new.k = stats::qgamma(as.numeric(log(a.new.p)), shape, scale = 1/rate, log.p = TRUE) #still gets inf
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

####################################################
####################################################

####################################################
####################################################

#' Post-hoc MCMC convergence diagnostics
#'
#' @description
#' Computes the same convergence diagnostics that
#' \code{CB_MCMC_single_rj_slice()} attaches when called with
#' \code{diagnostics = TRUE} (its default): per-gene RJMCMC acceptance rate,
#' per-gene effective sample size (ESS) for the rhythmicity indicator
#' \code{rho}, and per-gene ESS for phase \code{phi}. Use this on an MCMC
#' result that was run with \code{diagnostics = FALSE}; the underlying
#' \code{rho}, \code{phi}, and \code{if.accept.rj} matrices are always
#' stored in the returned object regardless of that flag, so diagnostics can
#' always be computed after the fact.
#'
#' @param mcmc_result List returned by \code{CB_MCMC_single_rj_slice()} (or
#'   compatible), containing \code{$rho}, \code{$phi}, and
#'   \code{$if.accept.rj}.
#' @param P Numeric; circadian period in hours, matching the value used for
#'   the MCMC run (default 24). Needed to convert \code{phi} to radians for
#'   its circular autocorrelation.
#'
#' @details
#' \code{phi} is periodic (an hour near \code{P} is one step from an hour
#' near 0), so its lag-1 autocorrelation is computed on the circle: phase is
#' converted to radians and the autocorrelation is the mean cosine of
#' successive phase differences, \code{mean(cos(phi_rad[-1] -
#' phi_rad[-K]))}, rather than a linear Pearson correlation on the raw
#' hour values (which would be distorted by wraparound near \code{0/P}).
#' \code{rho} is binary and unaffected by this, so its ESS still uses linear
#' lag-1 autocorrelation.
#'
#' @return A list: \code{acceptance_rate}, \code{ess_rho}, \code{ess_phi}
#'   (all per-gene, named by gene), \code{mean_acceptance_rate},
#'   \code{mean_ess_rho}, \code{mean_ess_phi}, \code{n_samples}, and
#'   \code{p_rhythmic_posterior}. Also printed as a summary; warns if mean
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

  # Per-gene RJMCMC acceptance rate: proportion of proposed jumps accepted.
  # if.accept.rj == -1 means stay was drawn (no jump proposed).
  proposed <- mcmc_result$if.accept.rj != -1
  n_proposed <- rowSums(proposed)
  n_accepted <- rowSums(mcmc_result$if.accept.rj == 1, na.rm = TRUE)
  accept_rate <- ifelse(n_proposed > 0, n_accepted / n_proposed, NA_real_)

  # Per-gene ESS for rho using lag-1 autocorrelation (no external dependencies).
  ess_rho <- apply(mcmc_result$rho, 1, function(x) {
    K <- length(x)
    if (K < 4 || var(x) < 1e-10) return(NA_real_)
    r1 <- tryCatch(cor(x[-K], x[-1]), error = function(e) NA_real_)
    if (is.na(r1) || abs(r1) >= 1) return(as.numeric(K))
    K * (1 - r1) / (1 + r1)
  })

  # Per-gene ESS for phi using circular lag-1 autocorrelation (phi wraps at
  # 0/P, so a linear correlation on the raw hour values would be wrong).
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

  # Per-gene posterior phase estimate: circular median over the rhythmic
  # draws only (phi is only biologically meaningful when rho = 1 for that
  # draw), paralleling p_rhythmic_posterior's per-gene summary for rho.
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
