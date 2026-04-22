CB_MCMC_single_rj_slice = function(Data.list, Init.value, P = 24,
                                iteration = 3000, thin = 20, n.burn=1000,
                                seed = 15213,
                                # p_rhythmic: prior Pr(rho=1) per gene (length G).
                                # Used in the RJMCMC log prior odds: log(p/(1-p)).
                                # Uniform prior: rep(0.2, G) means 20% prior rhythmicity.
                                # Gene-specific priors (e.g. from cosinor p-values) are supported.
                                p_rhythmic=rep(0.2, 100), rj.p.stay = 0.5,
                                A_prior= "Jeffreys_OLS_condi",
                                #options: "trunc_Normal",
                                # Jeffreys: p(log(A))\propto 1
                                # Jeffreys_ridge
                                # sq_expo: squared exponential p(A) \propto A*exp(-A^2/2)    
                                # gamma: gamma prior p(A) \ propto A^(alpha-1)exp(-rate*A)
                                mu_A = 1, sigma_A = 10^2, A.min = 0, # this applies to A~truncated normal
                                A_wb_beta2=2, #this applies to A~sq_expo(2, beta)
                                A_gm_shape=1.99, #A_shape<=2, flatter when closer to 2, 
                                A_gm_rate = 0.5, #the smaller the flatter
                                rj.phi = TRUE, rj.A = TRUE,
                                mu_M = 0, sigma_M = 10^2, 
                                sigma_prior_v = 4, sigma_prior_s = 1,
                                save.file = "MCMC_save.rds", 
                                save.file2 = "MCMC_save.rds"
){
  
  # Data.list=dat.input; Init.value=a.init;
  # P = 24; A.min = 0;
  # iteration = n.iter; thin = n.thin; n.burn=n.burn;
  # seed = 15213;
  # p_rhythmic=0.5; rj.p.stay = 0.5;
  # A_prior= "Jeffreys";
  # mu_A = 1; sigma_A = 10^2; # this applies to A~truncated normal
  # A_wb_beta2=2; #this applies to A~sq_expo(2; beta)
  # A_gm_shape=1.99; #A_shape<=2; flatter when closer to 2;
  # A_gm_rate = 0.5; #the smaller the flatter
  # rj.phi = TRUE; rj.A = TRUE;
  # mu_M = 0; sigma_M = 10^2;
  # sigma_prior_v = 4; sigma_prior_s = 1;
  # delta_MH_phi = pi/2;
  # save.file = paste0(out.dir, "/", out.name, "/",
  #                    "/tempSaveXXX", a.ind, ".rds")

  #observations
  omega = 2*pi/P
  Y = as.matrix(Data.list[[1]])
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
      rho.res = try_save(RJMCMC_single_slice(Y, t.c, t.s, N,
                                             t.c.sum, t.s.sum,
                                             c.t.phi, s.t.phi, cs.t.phi,
                                             y.t.c.sum, y.t.s.sum,
                                             AcosPhi, AsinPhi, A, phi,
                                             M, sigma, p_rhythmic, rho,
                                             rj.phi, rj.A,
                                             A_prior,
                                             mu_A, sigma_A, A.min,
                                             A_wb_beta2,
                                             A_gm_shape, A_gm_rate,
                                             omega, G, P, save.file2),
                         save, save.file)
      # rho.res=RJMCMC_single_slice(Y, t.c, t.s, N, 
      #                             t.c.sum, t.s.sum, 
      #                             c.t.phi, s.t.phi, cs.t.phi, 
      #                             y.t.c.sum, y.t.s.sum,
      #                             AcosPhi, AsinPhi, A, phi, 
      #                             M, sigma, rep(p_rhythmic, G), rho, 
      #                             rj.phi, rj.A,
      #                             A_prior, 
      #                             mu_A, sigma_A, A.min, 
      #                             A_wb_beta2, 
      #                             A_gm_shape, A_gm_rate,  
      #                             omega, G, P)
      if.accept.rj = rho.res$rho!=rho
      # log.r1 = rho.res$log.r1
      # log.r1.SS = rho.res$log.r1.SS
      # log.r3_A = rho.res$log.r3_A
      # log.r3_phi = rho.res$log.r3_phi
      rho = rho.res$rho
    }
    
    
    # within model move -------------------------------------------------------
    #update M    
    M = try_save(update_M_single(Y, t.c, t.s, N,
                                 AcosPhi, AsinPhi, sigma, rho,
                                 sigma_M, mu_M, omega, G),
                 save, save.file)
    # M = update_M_single(Y, t.c, t.s, N,
    #                     AcosPhi, AsinPhi, sigma, rho,
    #                     sigma_M, mu_M, omega, G)
    
    #udpate A
    csAphi = try_save(update_A_phi_slice(Y, X, XtX, beta_hat0,
                                         beta_cov0_11, beta_cov0_22, beta_cov0_rho,
                                         beta_mean0_1, beta_mean0_2,
                                         M, sigma, omega,
                                         AcosPhi, AsinPhi, G, P, A_prior,
                                         mu_A, sigma_A,
                                         A_wb_beta2,
                                         A_gm_shape, A_gm_rate),
                      save, save.file)
    # csAphi = update_A_phi_slice(Y, X, XtX, beta_hat0,
    #                             beta_cov0_11, beta_cov0_22, beta_cov0_rho,
    #                             beta_mean0_1, beta_mean0_2,
    #                             M, sigma, omega,
    #                             AcosPhi, AsinPhi, G, P, A_prior,
    #                             mu_A, sigma_A,
    #                             A_wb_beta2,
    #                             A_gm_shape, A_gm_rate)
    AcosPhi = csAphi$AcosPhi
    AsinPhi = csAphi$AsinPhi
    A = sqrt(AcosPhi^2+AsinPhi^2)
    phi = sapply(atan2(AsinPhi, AcosPhi), adjust.to.2pi)/omega
    
    sigma = try_save(update_sigma_single(Y, t.c, t.s, N, AcosPhi, AsinPhi, M, rho,
                                         sigma_prior_v, sigma_prior_s, omega, G),
                     save, save.file)
    # sigma = update_sigma_single(Y, t.c, t.s, N, AcosPhi, AsinPhi, M, rho,
    #                             sigma_prior_v, sigma_prior_s, omega, G)
    
    
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
    }else{
      rho.res = try_save(RJMCMC_single_slice(Y, t.c, t.s, N,
                                             t.c.sum, t.s.sum,
                                             c.t.phi, s.t.phi, cs.t.phi,
                                             y.t.c.sum, y.t.s.sum,
                                             AcosPhi, AsinPhi, A, phi,
                                             M, sigma, p_rhythmic, rho,
                                             rj.phi, rj.A,
                                             A_prior,
                                             mu_A, sigma_A, A.min,
                                             A_wb_beta2,
                                             A_gm_shape, A_gm_rate,
                                             omega, G, P, save.file2),
                         save, save.file)
      # rho.res=RJMCMC_single_slice(Y, t.c, t.s, N,
      #                             t.c.sum, t.s.sum,
      #                             c.t.phi, s.t.phi, cs.t.phi,
      #                             y.t.c.sum, y.t.s.sum,
      #                             AcosPhi, AsinPhi, A, phi,
      #                             M, sigma, rep(p_rhythmic, G), rho,
      #                             rj.phi, rj.A,
      #                             A_prior,
      #                             mu_A, sigma_A, A.min,
      #                             A_wb_beta2,
      #                             A_gm_shape, A_gm_rate,
      #                             omega, G, P)
      if.accept.rj = cbind(if.accept.rj, rho.res$rho!=rho) 
      #if the new rho is not the same as the old rho, it is accepted
      log.r1.store = cbind(log.r1.store, rho.res$log.r1)
      log.r1.SS.store = cbind(log.r1.SS.store, rho.res$log.r1.SS) 
      log.r3_A.store = cbind(log.r3_A.store, rho.res$log.r3_A)
      log.r3_phi.store = cbind(log.r3_phi.store, rho.res$log.r3_phi)
      rho = rho.res$rho
      rho.store = cbind(rho.store, rho)
    }    
    time1 = Sys.time()
    # samp.time$rho[iter]=time1-time0
    
    # within model move -------------------------------------------------------
    
    #update M    
    M = try_save(update_M_single(Y, t.c, t.s, N,
                                 AcosPhi, AsinPhi, sigma, rho,
                                 sigma_M, mu_M, omega, G),
                 save, save.file)
    M.store = cbind(M.store, M)
    # print("Finished M")
    time2 = Sys.time()
    # samp.time$M[iter]=time2-time1
    
    #udpate A
    csAphi = try_save(update_A_phi_slice(Y, X, XtX, beta_hat0,
                                         beta_cov0_11, beta_cov0_22, beta_cov0_rho,
                                         beta_mean0_1, beta_mean0_2,
                                         M, sigma, omega,
                                         AcosPhi, AsinPhi, G, P, A_prior,
                                         mu_A, sigma_A,
                                         A_wb_beta2,
                                         A_gm_shape, A_gm_rate),
                      save, save.file)
    # csAphi = update_A_phi_slice(Y, X, XtX, beta_hat0,
    #                             beta_cov0_11, beta_cov0_22, beta_cov0_rho,
    #                             beta_mean0_1, beta_mean0_2,
    #                             M, sigma, omega,
    #                             AcosPhi, AsinPhi, G, P, A_prior,
    #                             mu_A, sigma_A,
    #                             A_wb_beta2,
    #                             A_gm_shape, A_gm_rate)
    # stopifnot(sum(is.na(csAphi$AcosPhi))==0)
    # stopifnot(sum(is.na(csAphi$AsinPhi))==0)
    
    AcosPhi = csAphi$AcosPhi
    AsinPhi = csAphi$AsinPhi
    A = sqrt(AcosPhi^2+AsinPhi^2)
    phi = sapply(atan2(AsinPhi, AcosPhi), adjust.to.2pi)/omega
    
    AcosPhi.store = cbind(AcosPhi.store, AcosPhi); 
    AsinPhi.store = cbind(AsinPhi.store, AsinPhi); 
    A.store = cbind(A.store, A)
    phi.store = cbind(phi.store, phi)
    
    time4 = Sys.time()
    # samp.time$phi[iter]=time4-time2
    
    sigma = try_save(update_sigma_single(Y, t.c, t.s, N, AcosPhi, AsinPhi, M, rho,
                                         sigma_prior_v, sigma_prior_s, omega, G),
                     save, save.file)
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
                log.r1 = log.r1.store,
                log.r3_A = log.r3_A.store,
                log.r3_phi = log.r3_phi.store,
                if.accept.rj = if.accept.rj)
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

# Updating functions ------------------------------------------------------
#Gibbs sampling

#should update rho, A, and phi together with a reversible jump procedure
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

update_A_phi_slice = function(Y, X, XtX, beta_hat0, 
                              beta_cov0_11, beta_cov0_22, beta_cov0_rho, 
                              beta_mean0_1, beta_mean0_2, 
                              M, sigma, omega, 
                              AcosPhi, AsinPhi, G, P, A_prior,
                              mu_A, sigma_A, # this applies to A~truncated normal(mu_A, sigma_A)
                              A_wb_beta2, #this applies to A~sq_expo(2, beta)
                              A_gm_shape, A_gm_rate 
){
  beta1 = AcosPhi
  beta2 = AsinPhi
  Y_new = Y-M
  A = sqrt(beta1^2+beta2^2)
  #for the truncation model only
  f_upper_trunc = function(A0, u0){
    1/A0*exp(-1/sigma_A*(A0-mu_A)^2)-u0
  }
  f_upper_gamma = function(A0, u0){
    (A0)^(A_gm_shape-2)*exp(-1*A_gm_rate*A0)-u0
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
    # for(i in 1:500){
    #   print(i)
      # u = runif(G, 1e-3, f_upper_trunc(A, 0))
      # TEMP.SAVE.u <<-u
      # u = ifelse(u>1e3, 1e3, u)
      # beta_cov_11 = beta_cov0_11*sigma
      # beta_cov_22 = beta_cov0_22*sigma
      # A_max = sapply(1:G, function(g){
      #   uniroot(f_upper_trunc, c(0, 1000), u0=u[g], tol=1e-6)$root})
      # TEMP.SAVE.A_max <<-A_max
      # TEMP.SAVE.beta2 <<-beta2
      # max_beta1 = sqrt(A_max^2-beta2^2)
      # TEMP.SAVE.max_beta1 <<-max_beta1
      # max_beta1[is.na(max_beta1)]=0
      # beta1_mean = beta_hat[, 1]+beta_mean0_1*(beta2-beta_hat[, 2])
      # beta1_var = beta_cov_11*(1-beta_cov0_rho^2)
      # beta1_new = truncnorm::rtruncnorm(G, -1*max_beta1, max_beta1,
      #                                   beta1_mean, sqrt(beta1_var))
      # TEMP.SAVE.beta1_new <<-beta1_new
      # max_beta2 = sqrt(A_max^2-beta1_new^2)
      # TEMP.SAVE.max_beta2 <<-max_beta2
      # max_beta2[is.na(max_beta2)]=0
      # beta2_mean = beta_hat[, 2]+beta_mean0_2*(beta1_new-beta_hat[, 1])
      # beta2_var = beta_cov_22*(1-beta_cov0_rho^2)
      # beta2_new = truncnorm::rtruncnorm(G, -1*max_beta2, max_beta2,
      #                                   beta2_mean, sqrt(beta2_var))
      # TEMP.SAVE.beta2_new <<-beta2_new
      # stopifnot(sum(is.na(beta1_new))==0)
      # stopifnot(sum(is.na(beta2_new))==0)
    # }
    u = runif(G, 1e-3, f_upper_trunc(A, 0))
    u = ifelse(u>1e3, 1e3, u)
    beta_cov_11 = beta_cov0_11*sigma
    beta_cov_22 = beta_cov0_22*sigma
    A_max = sapply(1:G, function(g){
      uniroot(f_upper_trunc, c(0, 1000), u0=u[g], tol=1e-6)$root})
    max_beta1 = sqrt(A_max^2-beta2^2)
    max_beta1[is.na(max_beta1)]=1e-3
    beta1_mean = beta_hat[, 1]+beta_mean0_1*(beta2-beta_hat[, 2])
    beta1_var = beta_cov_11*(1-beta_cov0_rho^2)
    beta1_new = truncnorm::rtruncnorm(G, -1*max_beta1, max_beta1,
                                      beta1_mean, sqrt(beta1_var))
    max_beta2 = sqrt(A_max^2-beta1_new^2)
    max_beta2[is.na(max_beta2)]=1e-3
    beta2_mean = beta_hat[, 2]+beta_mean0_2*(beta1_new-beta_hat[, 1])
    beta2_var = beta_cov_22*(1-beta_cov0_rho^2)
    beta2_new = truncnorm::rtruncnorm(G, -1*max_beta2, max_beta2,
                                      beta2_mean, sqrt(beta2_var))
    
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
  }else if(A_prior=="sq_expo"){
    
    # A_prior=="sq_expo" ------------------------------------------------------
    
    get_beta_hat0 = sapply(1:G, function(g){
      # for(g in 1:G){
      XtXinv = solve(XtX+2*sigma[g]/A_wb_beta2) #woodbury matrix identity 
      beta_hat = Y_new[g, ]%*%t(XtXinv%*%t(X))
      beta_cov0 = XtXinv%*%XtX%*%XtXinv
      beta_cov = beta_cov0*sigma[g]
      new_beta = MASS::mvrnorm(n = 1, beta_hat, beta_cov, tol = 1e-5, empirical = FALSE, EISPACK = FALSE)
      return(new_beta)
    })
    beta1_new = get_beta_hat0[1, ]
    beta2_new = get_beta_hat0[2, ]
  }else if(A_prior=="gamma_OLS"){
    
    # A_prior=="gamma_OLS" ---------------------------------------
    
    #OLS solution 
    beta_hat = Y_new%*%beta_hat0
    u = runif(G, 1e-5, f_upper_gamma(A, 0)) 
    u = ifelse(u>1e3, 1e3, u)
    beta_cov_11 = beta_cov0_11*sigma
    beta_cov_22 = beta_cov0_22*sigma
    A_max = sapply(1:G, function(g){
      uniroot(f_upper_gamma, c(0, 1000), u0=u[g])$root})
    max_beta1 = sqrt(A_max^2-beta2^2)
    beta1_mean = beta_hat[, 1]+beta_mean0_1*(beta2-beta_hat[, 2])
    beta1_var = beta_cov_11*(1-beta_cov0_rho^2)
    beta1_new = truncnorm::rtruncnorm(G, -1*max_beta1, max_beta1, 
                                      beta1_mean, sqrt(beta1_var)) 
    max_beta2 = sqrt(A_max^2-beta1_new^2)
    beta2_mean = beta_hat[, 2]+beta_mean0_2*(beta1_new-beta_hat[, 1])
    beta2_var = beta_cov_22*(1-beta_cov0_rho^2)
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
                      mu_A, sigma_A, A.min, # this applies to A~truncated normal(mu_A, sigma_A)I(A>A_min)
                      A_wb_beta2, #this applies to A~sq_expo(2, beta)
                      A_gm_shape, A_gm_rate,  
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
    A.post = truncnorm::dtruncnorm(A, a = A.min, b = Inf, a.mean, sqrt(a.var))
    
  }else if(A_prior == "Jeffreys_OLS_condi"){

    a.var = sigma/T.x 
    a.mean = U.x/T.x 
    A.post = truncnorm::dtruncnorm(A, a = A.min, b = Inf, a.mean, sqrt(a.var))
    
  }else if(A_prior=="sq_expo"){

    alpha = 2
    beta = T.x/sigma/2+1/A_wb_beta2
    gamma = U.x/sigma
    # fox=sapply(1:G, function(g){
    #   FoxWright(alpha, beta[g], gamma[g])
    # })
    # fox1=sapply(1:G, function(g){
    #   FoxWright1(alpha, beta[g], gamma[g])
    # })
    # fox2=sapply(1:G, function(g){
    #   FoxWright2(alpha, beta[g], gamma[g])
    # })
    # library(Rmpfr)
    # fox3=sapply(1:G, function(g){
    #   FoxWright3(alpha, beta[g], gamma[g])
    # })
    # A.post = 2*beta^(alpha/2)*A^(alpha-1)*exp(-1*beta*A^2+gamma*A)/fox
    # pdf("/home/xix66/CircaBayes/program_V1_test/sim1_Aim1_slice/test/fox023.pdf")
    # a.min = min(c(fox, fox2, Rmpfr::asNumeric(fox3)))
    # par(mfrow=c(2, 2))
    # plot(fox, fox1, xlim = c(a.min, 2))
    # abline(a=0, b=1)
    # plot(fox, fox2, xlim = c(a.min, 2))
    # abline(a=0, b=1)
    # plot(fox, Rmpfr::asNumeric(fox3), xlim = c(a.min, 2))
    # abline(a=0, b=1)
    # plot(fox2, Rmpfr::asNumeric(fox3), xlim = c(a.min, 2))
    # abline(a=0, b=1)
    # dev.off()
    
    A.post = A_post_sq_expo(A, alpha, beta, gamma)
    
  }else if(A_prior=="gamma_OLS"){

    alpha = A_gm_shape
    beta = T.x/sigma/2+A_gm_rate
    gamma = U.x/sigma
    fox=sapply(1:G, function(g){
      FoxWright(alpha, beta[g], gamma[g])
    })
    A.post = 2*beta^(alpha/2)*A^(alpha-1)*exp(-1*beta*A^2+gamma*A)/fox
  }
  return(A.post)
}

get_prior_A = function(Y, t.c, t.s, N, 
                      AcosPhi, AsinPhi, A, phi, 
                      M, sigma, 
                      A_prior, 
                      mu_A, sigma_A, A.min, # this applies to A~truncated normal(mu_A, sigma_A)I(A>A_min)
                      A_wb_beta2, #this applies to A~sq_expo(2, beta)
                      A_gm_shape, A_gm_rate,  
                      omega){
  
  if(A_prior=="trunc_Normal_OLS_condi"){
    A.prior = truncnorm::dtruncnorm(A, a = A.min, b = Inf, mu_A, sqrt(sigma_A))
    
  }else if(A_prior == "Jeffreys_OLS_condi"){
    A.prior = 1 #???
    
  }else if(A_prior=="sq_expo"){
    A.prior = stats::dweibull(A, 2, A_wb_beta2)
    
  }else if(A_prior=="gamma_OLS"){
    A.prior = stats::dgamma(A, shape = A_gm_shape, rate = A_gm_rate)
    
  }
  return(A.prior)
}



update_A_single = function(a.Y, a.t.c, a.t.s, a.N,
                           # Y.list, t.c.list, t.s.list, 
                           phi.vec, M.vec, sigma.vec, 
                           # N.vec, phi.mat, M.mat, sigma.mat, 
                           mu_A, sigma_A, omega, G, A.min){
  post.A.params = get_post_A_params(a.Y, a.t.c, a.t.s, a.N,
                                    phi.vec, M.vec, sigma.vec, 
                                    mu_A, sigma_A, omega)
  # xA <<- list(a.var = post.A.params$var, a.mean = post.A.params$mean)
  # stopifnot("A problem" = (all(xA$a.var<10000))&all(is.finite(xA$a.mean)))
  a.A = truncnorm::rtruncnorm(G, a = A.min, b = Inf, 
                              post.A.params$mean, sqrt(post.A.params$var))
}

# get_log_phi_post = function(a.Y.list, a.t.c, a.t.s, a.n,
#                             a.M.vec, a.A.vec, a.phi.vec, a.sigma.vec, a.rho.vec,
#                             kappa.vec, theta.vec, 
#                             omega){
#   # b.c = a.A.vec*a.rho.vec*cos(omega*a.phi.vec) #G*1 mat
#   # b.s = a.A.vec*a.rho.vec*sin(omega*a.phi.vec)
#   #when A is present, phi is always present, so use rho = 1 for all
#   b.c = a.A.vec*cos(omega*a.phi.vec) #G*1 mat
#   b.s = a.A.vec*sin(omega*a.phi.vec)
#   a.cosine.1 = cbind(a.M.vec, b.c, b.s) %*% rbind(1, a.t.c, a.t.s)
#   diff1 = apply((a.Y.list- a.cosine.1)^2*a.rho.vec, 1, sum)#/sum(rho.mat[, a])*a.n
# 
#   -1/2*diff1/a.sigma.vec+kappa.vec*cos(a.phi.vec*omega-theta.vec)
# }

#fix: delete rho in phi calc
get_log_phi_post_single = function(a.Y.list, a.t.c, a.t.s, a.n,
                                   a.M.vec, AcosPhi, AsinPhi, a.sigma.vec,
                                   # kappa.vec, theta.vec, 
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
    # xrho <<- list(a.cosine.1 = a.cosine.1,
    #               lik1 = lik1,
    #               lik0 = lik0,
    #               odds_rho1 = odds_rho1,
    #               p.rho1 = p.rho1,
    #               a.rho = a.rho,
    #               p.mat = p.mat,
    #               sigma.mat = sigma.mat)
    # stopifnot("rho problem" = all(a.rho==1|a.rho==0))
    return(a.rho)
  })
  do.call(cbind, rho.per.group)
}

# update_p = function(rho.mat){
#   #we update p_gj based on rho_gj
#   rho.mean = apply(rho.mat, 1, mean) #mean of all groups
#   p.per.group = lapply(1:ncol(rho.mat), function(a){
#     #since rho_gj is only one observation, try a weighted version of global mean and single-group mean
#     a.rho = stats::rbeta(G, rho.mat[, a]+rho.mean+1, 1-rho.mat[, a]+1-rho.mean+1)
#   })
#   do.call(cbind, p.per.group)
# }

# get_log_phi_post_full_single(a.Y.vec, a.t.c.vec, a.t.s.vec,
#                                         a.t.c.sum.num, a.t.s.sum.num,
#                                         c.t, s.t, cs.t,
#                                         a.y.t.c.sum.num, a.y.t.s.sum.num,
#                                         a.sigma.num, a.A.num, a.phi.num, a.M.num,
#                                         omega, P, save.file2=paste0(out.dir, "/", out.name, "/CBt_temp_debug3_", 700, ".rds"))
get_log_phi_post_full_single = function(a.Y.vec, a.t.c.vec, a.t.s.vec, 
                                        a.t.c.sum.num, a.t.s.sum.num,
                                        c.t, s.t, cs.t,
                                        a.y.t.c.sum.num, a.y.t.s.sum.num,
                                        a.sigma.num, a.A.num, a.phi.num, a.M.num, 
                                        omega, P, save.file2){
  
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
  # phi.fun = function(phi){
  #   exp(-r2/2*c.t*cos(omega*phi)^2-r2/2*s.t*sin(omega*phi)^2-r2*cs.t*sin(omega*phi)*cos(omega*phi)
  #       +a.r/a.sigma.num*c.t_Y_M*cos(omega*phi)+a.r*s.t_Y_M/a.sigma.num*sin(omega*phi)+c.kt*cos(omega*phi)+s.kt*sin(omega*phi))
  # }
  #
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
    
    # a.phi=seq(5.5, 6.5, by=0.001)
    # a.phi=seq(17.5, 18.5, by=0.001)
    # a.phi=seq(13, 16, by=0.0000001)
    # a.phi=seq(0.0001, P, by=0.001)
    # par(mfrow=c(1, 2))
    # plot(a.phi, phi.fun.log(a.phi), type="l")
    # # plot(a.phi, phi.fun.log.dx(a.phi), type="l")
    # plot(a.phi, phi.fun(a.phi, b1.scale), type="l")

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
    #out
    #2.729303
  }
  
  if(("error" %in% class(b1.exp))|b1.exp[[1]]==0){
    # out = 0 #when there is still no return (i reaches 20), it is probabily because of a very large integral, so we just give out = 0. #check later!
    save(a.Y.vec, a.t.c.vec, a.t.s.vec, 
         a.t.c.sum.num, a.t.s.sum.num,
         c.t, s.t, cs.t,
         a.y.t.c.sum.num, a.y.t.s.sum.num,
         a.sigma.num, a.A.num, a.phi.num, a.M.num, 
         omega, P,
         file=save.file2)
    stop("out is NULL")
  }else{
    b1 = log(b1.exp[[1]])
    out = a1+b1.scale-b1+exp.scale
    #out
    #2.729303
  }
  return(out)
}

RJMCMC_single_slice = function(Y, t.c, t.s, N, 
                               t.c.sum, t.s.sum, 
                               c.t, s.t, cs.t, 
                               y.t.c.sum.num, y.t.s.sum.num,
                               AcosPhi, AsinPhi, A, phi,  
                               M, sigma, p.vec = rep(0.5, G), rho, 
                               propose.phi = FALSE, propose.A = TRUE,
                               A_prior, 
                               mu_A, sigma_A, A.min, # this applies to A~truncated normal(mu_A, sigma_A)I(A>A_min)
                               A_wb_beta2, #this applies to A~sq_expo(2, beta)
                               A_gm_shape, A_gm_rate,  
                               omega, G, P, save.file2){
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
                        mu_A, sigma_A, A.min, # this applies to A~truncated normal(mu_A, sigma_A)I(A>A_min)
                        A_wb_beta2, #this applies to A~sq_expo(2, beta)
                        A_gm_shape, A_gm_rate,  
                        omega, G)
    log.A.post = ifelse(A.post==0, 1e-7, log(A.post))
    A.prior = get_prior_A(Y, t.c, t.s, N, 
                         AcosPhi, AsinPhi, A, phi, 
                         M, sigma, 
                         A_prior, 
                         mu_A, sigma_A, A.min, # this applies to A~truncated normal(mu_A, sigma_A)I(A>A_min)
                         A_wb_beta2, #this applies to A~sq_expo(2, beta)
                         A_gm_shape, A_gm_rate,  
                         omega)
    log.A.prior = ifelse(A.prior==0, 1e-7, log(A.prior))
  }else{
    log.A.prior = 1; log.A.post = 1
  }
  
  if(propose.phi){
    log.phi.post = sapply(1:G, function(g){
      # print(g)
      get_log_phi_post_full_single(Y[g, ], t.c, t.s, 
                                   t.c.sum, t.s.sum, 
                                   c.t, s.t, cs.t, 
                                   y.t.c.sum.num[g], y.t.s.sum.num[g],       
                                   sigma[g], A[g], phi[g], M[g],  
                                   omega, P, save.file2)
    })
    # for(g in 1:G){
    #   print(g)
    #   get_log_phi_post_full_single(Y[g, ], t.c, t.s, t.c.sum, t.s.sum, c.t, s.t, cs.t, y.t.c.sum.num[g], y.t.s.sum.num[g], sigma[g], A[g], M[g], phi[g], omega, P)
    # }
    log.phi.prior = 1/P #the prior of phi is U(0, 1), so log(1)=0
  }else{
    log.phi.prior = 1; log.phi.post = 1
  }
  log.r1 = -1/(2*sigma)*(log.lik_jump0-log.lik_cur0)
  log.r2 = log((p.vec+1e-5)/(1-p.vec+1e-5))#*((-1)^rho) #already taken care below
  log.r3 = log.A.prior+log.phi.prior-log.A.post-log.phi.post
  log.r = log.r1+(log.r2+log.r3)*((-1)^rho) #if rho.mat = 1, then reverse
  # gene.idx = 5; log.r1[gene.idx]; log.r2[gene.idx]; log.r3[gene.idx]; log.r[gene.idx]; rho.mat[gene.idx, a]
  # hist(log.r3[1:39]*((-1)^rho.mat[1:39, a]))
  # hist((log.r3[40:72]*((-1)^rho.mat[40:72, a]))[rho.mat[40:72, a]==0])
  # hist((log.r1[40:72])[rho.mat[40:72, a]==0]) #should not jump
  # hist(log.r[40:72][rho.mat[40:72, a]==0])
  # hist(log.r[40:72][rho.mat[40:72, a]==1])
  a.jump.log.u = log(runif(G, 0, 1))
  a.rho = ifelse(log.r>a.jump.log.u, 1-rho, rho)
  # xrho <<- list(a.rho = a.rho,
  #               log.r1 = log.r1,
  #               log.r3_A = log.A.prior-log.A.post,
  #               log.r3_phi = log.phi.prior-log.phi.post)
  #   return(list(a.rho = a.rho,
  #               log.r1 = log.r1,
  #               log.r3_A = (log.A.prior-log.A.post)*((-1)^rho),
  #               log.r3_phi = (log.phi.prior-log.phi.post)*((-1)^rho)
  #   ))
  # })
  # do.call(cbind, rho.per.group)
  return(list(rho = a.rho,
              log.r1 = log.r1*((-1)^rho),
              log.r1.SS = (log.lik_jump0-log.lik_cur0)*((-1)^rho),
              log.r3_A = (log.A.prior-log.A.post)*((-1)^rho),
              log.r3_phi = (log.phi.prior-log.phi.post)*((-1)^rho)
  ))
}


update_theta_kappa = function(rho.mat, phi.mat, theta.vec, kappa.vec, VM_w.vec,
                              # VM_x.vec, VM_v.vec, VM_u.mat, 
                              VM_theta0, VM_R0, VM_c, VM_Bessel_k, omega, G, 
                              t_p = FALSE){
  # rho.mat=rho; phi.mat=phi; theta.vec=theta; kappa.vec=kappa;
  # VM_x.vec=VM_x; VM_v.vec=VM_v; VM_u.mat=VM_u; VM_w.vec=VM_w;
  # save(rho, phi, theta, kappa, VM_x, VM_v, VM_u, VM_w, a.seed, paste0("RJMCMC_tests/theta_kappa_error.rda"))
  # a.seed = 3
  # set.seed(a.seed)
  #Von mises latent variables initiation
  if(t_p){
    #when t_p = TRUE, phi.mat is the t_p matrix. 
    VM_Rn.cos_theta_n = VM_R0*cos(VM_theta0)+sum(cos(phi.mat*omega))
    VM_Rn.sin_theta_n = VM_R0*sin(VM_theta0)+sum(sin(phi.mat*omega))
    VM_m = VM_c+length(phi.mat) 
  }else{
    VM_Rn.cos_theta_n = VM_R0*cos(VM_theta0)+apply(cos(phi.mat*omega)*rho.mat, 1, sum)
    VM_Rn.sin_theta_n = VM_R0*sin(VM_theta0)+apply(sin(phi.mat*omega)*rho.mat, 1, sum)
    VM_m = VM_c+apply(rho.mat, 1, sum)
  }
  VM_Rn = sqrt(VM_Rn.cos_theta_n^2+VM_Rn.sin_theta_n^2)
  cos_theta_n = VM_Rn.cos_theta_n/VM_Rn
  sin_theta_n = VM_Rn.sin_theta_n/VM_Rn
  tan_theta_n = sin_theta_n/cos_theta_n
  VM_theta_n0 = atan(tan_theta_n)
  VM_theta_n = ifelse(cos_theta_n<0, VM_theta_n0+pi, VM_theta_n0+2*pi*(tan_theta_n<0)) #adjust the atan output to [0, 2*pi)
  # #update
  # a.VM_x = runif(G, 0, VM_w.vec^(VM_m-1)) #can be 0.
  a.log.VM_x = (VM_m-1)*log(VM_w.vec) - stats::rexp(G, 1) #can be 0.
  # a.VM_v = runif(G, 0, exp(VM_Rn*kappa.vec*(1+cos(theta.vec-VM_theta_n))))
  a.log.VM_v =  VM_Rn*kappa.vec*(1+cos(theta.vec-VM_theta_n)) - stats::rexp(G, 1)
  a.VM_vn = a.log.VM_v/(VM_Rn*(1+cos(theta.vec-VM_theta_n)))
  stopifnot("a.VM_vn NA" = sum(is.na(a.VM_vn))==0)
  # stopifnot("a.VM_vn should be smaller than kappa" = all(a.VM_vn < kappa.vec))
  #update theta.vec
  rhs = a.log.VM_v/(VM_Rn*kappa.vec)-1 #the right hand side of Gibbs sampling step (4) in Damien&Walker
  # my.save <<- list(VM_Rn = VM_Rn, VM_theta_n = VM_theta_n, kappa.vec = kappa.vec, a.VM_v = a.VM_v, rhs = rhs,
  #                  phi.mat = phi.mat, rho.mat = rho.mat)
  stopifnot("rhs NA" = sum(is.na(rhs))==0)
  stopifnot("cos(x) will always be smaller than 1" = all(rhs<=1)) #later: if never triggered will be deleted
  
  a.theta = sapply(1:length(rhs), function(a){
    a.rhs = rhs[a]
    if(a.rhs< -1){
      a.theta = runif(1, 0, 2*pi)
    }else{
      a.rhs.radis = acos(a.rhs)
      a.rhs.radis2 = 2*pi-a.rhs.radis
      a.theta1 = runif(1, VM_theta_n[a], a.rhs.radis+VM_theta_n[a])
      a.theta2 = runif(1, VM_theta_n[a]+a.rhs.radis2, 2*pi+VM_theta_n[a])
      a.temp.samp = sample(c(1, 2), 1, replace = TRUE)
      a.theta = ifelse(a.temp.samp==1, a.theta1, a.theta2)
      a.theta = adjust.to.2pi(a.theta)
    }
    return(a.theta)
  })
  a.theta = unlist(a.theta)
  stopifnot("theta problem"= sum(is.na(a.theta))==0&all(is.finite(a.theta)))
  #update VM_E and VM_N instead of VM_u
  a.VM_E = VM_w.vec+stats::rexp(G, rate = (besselI(kappa.vec, 0)-1)) #mean is (besselI(kappa.vec)-1)^(-1)
  a.VM_N = do.call(cbind, lapply(1:VM_Bessel_k, function(a.k){
    lambda.a.k = factorial(a.k)^(-2)*((1/2)^(2*a.k))
    F_k = stats::rexp(G, rate = VM_w.vec*lambda.a.k*kappa.vec^(2*a.k))
    stopifnot("F_k NA" = sum(is.na(F_k))==0)
    stopifnot("N_k should be larger than kappa" = all(kappa.vec*(1+F_k)^(1/2/a.k) > kappa.vec))
    kappa.vec*(1+F_k)^(1/2/a.k)
  }))
  a.VM_N.min = apply(a.VM_N, 1, min)
  #update VM_w
  # a.VM_w = sapply(1:G, function(g){truncdist::rtrunc(1, "gamma", a = a.VM_x[g]^(1/(VM_m[g]-1)), b = a.VM_E[g], shape=1, scale=1)})
  a.val = ifelse(exp(a.log.VM_x/(VM_m-1))==0, min(0.00001, runif(1, min(a.VM_E), 0.00001)), exp(a.log.VM_x/(VM_m-1)))
  a.VM_w = sapply(1:G, function(g){rtruncgamma(1, a = a.val[g],# a.VM_x[g]^(1/(VM_m[g]-1)),
                                               b = a.VM_E[g], shape=1, rate = 1)})
  # print(paste0(which.min(a.VM_w), ": ", min(a.VM_w)))
  # print("Finished w")
  #update kappa
  # a.kappa = sapply(1:G, function(g){truncdist::rtrunc(1, "gamma", a = max(0, a.VM_vn[g]), b = a.VM_N.min[g], shape=1, rate = VM_Rn[g])})
  a.kappa = sapply(1:G, function(g){rtruncgamma(1, a = max(0.00001, a.VM_vn[g]), b = a.VM_N.min[g], shape=1, rate = VM_Rn[g])})
  # print(paste0(which.max(a.kappa), ": ", max(a.kappa)))
  # print("Finished kappa")
  stopifnot("kappa problem" = all(a.VM_vn<a.VM_N.min))
  stopifnot("Kappa problem"= sum(is.na(a.kappa))==0&all(is.finite(a.kappa)))
  
  return(list(#a.VM_x = a.VM_x,
    #a.VM_v = a.log.VM_v,
    #a.VM_E = a.VM_E,
    #a.VM_N = a.VM_N,
    a.VM_w = a.VM_w,
    a.theta = a.theta,
    a.kappa = a.kappa))
}

# update_theta_kappa = function(rho.mat, phi.mat, theta.vec, kappa.vec,
#                               VM_x.vec, VM_v.vec, VM_u.mat, VM_w.vec, 
#                               VM_theta0, VM_R0, VM_c, VM_Bessel_k, omega, G){
#   #Von mises latent variables initiation
#   VM_Rn.cos_theta_n = VM_R0*cos(VM_theta0)+apply(cos(phi.mat*omega)*rho.mat, 1, sum)
#   VM_Rn.sin_theta_n = VM_R0*sin(VM_theta0)+apply(sin(phi.mat*omega)*rho.mat, 1, sum)
#   VM_Rn = sqrt(VM_Rn.cos_theta_n^2+VM_Rn.sin_theta_n^2)
#   cos_theta_n = VM_Rn.cos_theta_n/VM_Rn
#   sin_theta_n = VM_Rn.sin_theta_n/VM_Rn
#   tan_theta_n = sin_theta_n/cos_theta_n
#   VM_theta_n0 = atan(tan_theta_n)
#   VM_theta_n = ifelse(cos_theta_n<0, VM_theta_n0+pi, VM_theta_n0+2*pi*tan_theta_n<0) #adjust the atan output to [0, 2*pi)
#   VM_m = VM_c+apply(rho.mat, 1, sum)
#   # a.VM_x = runif(G, 0, VM_w.vec^(VM_m-1))
#   a.log.VM_x = (VM_m-1)*log(VM_w.vec)- stats::rexp(G, 1)
#   # a.VM_v = runif(G, 0, exp(VM_Rn*kappa.vec*(1+cos(theta.vec-VM_theta_n))))
#   a.log.VM_v =  VM_Rn*kappa.vec*(1+cos(theta.vec-VM_theta_n)) - stats::rexp(G, 1)
#   a.VM_vn = a.log.VM_v/(VM_Rn*(1+cos(theta.vec-VM_theta_n)))
#   stopifnot("a.VM_vn should be smaller than kappa" = all(a.VM_vn < kappa.vec))
#   #update theta.vec
#   rhs = a.log.VM_v/(VM_Rn*kappa.vec)-1 #the right hand side of Gibbs sampling step (4) in Damien&Walker
#   stopifnot("cos(x) will always be smaller than 1" = all(rhs<=1)) #later: if never triggered will be deleted
#   
#   a.theta = sapply(1:length(rhs), function(a){
#     a.rhs = rhs[a]
#     if(a.rhs< -1){
#       a.theta = runif(1, 0, 2*pi)
#     }else{
#       a.rhs.radis = acos(a.rhs)
#       a.rhs.radis2 = 2*pi-a.rhs.radis
#       a.theta1 = runif(1, VM_theta_n[a], a.rhs.radis+VM_theta_n[a])
#       a.theta2 = runif(1, VM_theta_n[a]+a.rhs.radis2, 2*pi+VM_theta_n[a])
#       a.temp.samp = sample(c(1, 2), 1, replace = TRUE)
#       a.theta = ifelse(a.temp.samp==1, a.theta1, a.theta2)
#       a.theta = adjust.to.2pi(a.theta)
#     }
#     return(a.theta)
#   })
#   #update VM_E and VM_N instead of VM_u
#   a.VM_E = VM_w.vec+stats::rexp(G, rate = (besselI(kappa.vec, 0)-1)) #mean is (besselI(kappa.vec)-1)^(-1)
#   a.VM_N = do.call(cbind, lapply(1:VM_Bessel_k, function(a.k){
#     lambda.a.k = factorial(a.k)^(-2)*((1/2)^(2*a.k))
#     F_k = stats::rexp(G, rate = VM_w.vec*lambda.a.k*kappa.vec^(2*a.k))
#     stopifnot("N_k should be larger than kappa" = all(kappa.vec*(1+F_k)^(1/2/a.k) > kappa.vec))
#     kappa.vec*(1+F_k)^(1/2/a.k)
#   }))
#   a.VM_N.min = apply(a.VM_N, 1, min)
#   #update VM_w
#   # a.VM_w = sapply(1:G, function(g){truncdist::rtrunc(1, "gamma", a = a.VM_x[g]^(1/(VM_m[g]-1)), b = a.VM_E[g], shape=1, scale=1)})
#   # a.VM_w = sapply(1:G, function(g){rtruncgamma(1, a = a.VM_x[g]^(1/(VM_m[g]-1)), b = a.VM_E[g], shape=1, rate = 1)})
#   a.log.VM_w = sapply(1:G, function(g){rtrunclgamma(1, a = (1/(VM_m[g]-1))*a.log.VM_x[g], b = log(a.VM_E[g]), shape=1, rate = 1)})
#   
#   # print("Finished w")
#   #update kappa
#   # a.kappa = sapply(1:G, function(g){truncdist::rtrunc(1, "gamma", a = max(0, a.VM_vn[g]), b = a.VM_N.min[g], shape=1, rate = VM_Rn[g])})
#   # a.kappa = sapply(1:G, function(g){rtruncgamma(1, a = max(0, a.VM_vn[g]), b = a.VM_N.min[g], shape=1, rate = VM_Rn[g])})
#   a.kappa = sapply(1:G, function(g){rtrunclgamma(1, a = max(-10^10, log(a.VM_vn[g])), b = log(a.VM_N.min[g]), shape=1, rate = VM_Rn[g])})
#   # print("Finished kappa")
#   # xthekp <<- list(a.VM_N = a.VM_N, a.kappa = a.kappa, VM_Rn = VM_Rn, a.log.VM_v = a.log.VM_v, a.VM_vn = a.VM_vn, VM_theta_n = VM_theta_n,
#   #                 VM_w.vec = VM_w.vec, kappa.vec = kappa.vec, a.theta = a.theta)
#   stopifnot("VM_w problem" = all(a.VM_x^(1/(VM_m-1))<a.VM_E))
#   stopifnot("kappa problem" = all(a.VM_vn<a.VM_N.min))
#   
#   stopifnot("Kappa problem"= sum(is.na(a.kappa))==0&all(is.finite(a.kappa)))
#   
#   return(list(a.VM_x = a.VM_x,
#               a.VM_v = a.log.VM_v,
#               a.VM_E = a.VM_E,
#               a.VM_N = a.VM_N,
#               a.VM_w = a.VM_w,
#               a.theta = a.theta,
#               a.kappa = a.kappa))
# }


# Other functions ---------------------------------------------------------
try_save = function(expr, out, save.file){
  a.res = try(expr)
  if("try-error" %in% class(a.res)){
    saveRDS(out, save.file)
  }else{
    return(a.res)
  }
}
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
#test lgamma
# install.packages("VGAM")
# x.l = VGAM::rlgamma(1000, location = 0, shape = 1, scale = 1)
# x = rgamma(1000, shape = 1, scale = 2)
# # x = rgamma(1000, shape = 1, rate = 1/10)
# par(mfrow = c(2, 1))
# hist(x.l)
# hist(log(x))
# hist(exp(x.l))
# hist(x)
# pgamma(2, 1, 1)
# VGAM::plgamma(log(2), 1, 1)

####################################################
# FoxWright Method 1
####################################################

#' @export
FoxWright1<-function(alpha,beta,gamma,eps=.00001,log=FALSE){
  u=alpha;v=gamma;w=beta;
  j=0
  x0=0
  x1=lgamma(u/2+j/2)+j*log(abs(v)/sqrt(w))-lgamma(1)-lfactorial(j)
  while(abs(x1-x0)>eps && exp(x1)!=0){
    j=j+1
    x0=x1
    x1=lgamma(u/2+j/2)+j*log(abs(v)/sqrt(w))-lgamma(1)-lfactorial(j)
  }
  
  SumTimes=j
  ######################
  
  i=seq(0,SumTimes,by=1)
  if(v>0){
    logitem=lgamma(u/2+i/2)+i*log(v/sqrt(w))-lgamma(1)-lfactorial(i)
    Max_of_terms=max(logitem)
    logitem_Minus_Maxlog=logitem-Max_of_terms
    item_star=exp(logitem_Minus_Maxlog)
    fox=exp(Max_of_terms)*sum(item_star)
    logfox=log(fox)}
  else if (v<0) {
    logitem=lgamma(u/2+i/2)+i*log(abs(v)/sqrt(w))-lgamma(1)-lfactorial(i)
    Max_of_terms=max(logitem)
    logitem_Minus_Maxlog=logitem-Max_of_terms
    item_star=exp(logitem_Minus_Maxlog)
    fox=exp(Max_of_terms)*sum(item_star*(-1)^i)
    logfox=log(fox)
  } else {
    fox=0
    logfox=-Inf
  }
  
  result=ifelse(log==FALSE,fox,logfox)
  return(result)
  
}

####################################################
# FoxWright Method 2
####################################################

#' @export
FoxWright2<-function(alpha,beta,gamma,eps=.00001,log=FALSE){
  u=alpha;v=gamma;w=beta;
  
  RobustFox<-function(u,v,w,log,SumTimes){
    i=seq(0,SumTimes,by=1)
    if(v>0){
      logitem=lgamma(u/2+i/2)+i*log(v/sqrt(w))-lgamma(1)-lfactorial(i)
      item=exp(logitem)
      fox=sum(item)
      logfox=log(fox)}
    else if (v<0) {
      logitem=lgamma(u/2+i/2)+i*log(abs(v)/sqrt(w))-lgamma(1)-lfactorial(i)
      item=exp(logitem)*(-1)^i
      fox=sum(item)
      logfox=log(fox)
    } else {
      fox=0
      logfox=-Inf
    }
    
    result=ifelse(log==F,fox,logfox)
    return(result)
  }
  
  
  # eps: converge level (difference of two numbers)
  ## find the SumTimes which makes RobustFox function converge
  i=1
  x0=0
  x1=RobustFox(u,v,w,SumTimes=i,log)
  while(abs(x1-x0)>eps){
    i=i+1
    x0=x1
    x1=RobustFox(u,v,w,SumTimes=i,log)
  }
  
  return(FoxWright=x1)
}

####################################################
# FoxWright Method 3
####################################################

#install.packages('Rmpfr')
#' @export
FoxWright3<-function(alpha,beta,gamma,eps=.00001,log=FALSE){
  u=alpha;v=gamma;w=beta;
  j=1
  x0=0
  x1=lgamma(u/2+j/2)+j*log(abs(v)/sqrt(w))-lgamma(1)-lfactorial(j)
  while(abs(x1-x0)>eps && round(exp(x1),2)!=0){
    j=j+1
    x0=x1
    x1=lgamma(u/2+j/2)+j*log(abs(v)/sqrt(w))-lgamma(1)-lfactorial(j)
  }
  
  SumTimes=j
  ########################################
  i=seq(0,SumTimes,by=1)
  if(v>0){
    logitem=lgamma(u/2+i/2)+i*log(v/sqrt(w))-lgamma(1)-lfactorial(i)
    
    logitem <- mpfr(logitem, precBits = 106)
    item=exp(logitem)
    fox=sum(item)
    logfox=log(fox)}
  else if (v<0) {
    logitem=lgamma(u/2+i/2)+i*log(abs(v)/sqrt(w))-lgamma(1)-lfactorial(i)
    logitem <- mpfr(logitem, precBits = 106)
    item=exp(logitem)*(-1)^i
    fox=sum(item)   # this is a list type
    logfox=log(fox)
  } else {
    fox=0
    logfox=-Inf
  }
  
  
  fox <- capture.output(fox)[2]
  fox <- substr(fox,5,nchar(fox))
  
  logfox <- capture.output(logfox)[2]
  logfox <- substr(logfox,5,nchar(logfox))
  
  result=ifelse(log==FALSE,fox,logfox)
  return(noquote(result))
}

####################################################
# FoxWright Method 4
####################################################

#' @export
FoxWright<-function(alpha,beta,gamma,eps=.00001,log=FALSE){
  u=alpha;v=gamma;w=beta;
  ###### part 1 #######
  k=0
  x1=lgamma(u/2+k)-lgamma(u/2)+k*log(v^2)-k*log(w)-lfactorial(2*k)
  x0=x1+1
  while(exp(lgamma(u/2))*exp(x1)!=0 |
        abs(exp(x1)-exp(x0))>0 |
        k<20){
    k=k+1
    x0=x1
    x1=lgamma(u/2+k)-lgamma(u/2)+k*log(v^2)-k*log(w)-lfactorial(2*k)
  }
  T1=k
  ###### part 2 #######
  k=0
  x1=lgamma((u+1)/2+k)-lgamma((u+1)/2)+k*log(v^2)-k*log(w)-lfactorial(2*k+1)
  x0=x1+1
  while(exp(lgamma((u+1)/2))*(v/sqrt(w))*exp(x1)!=0 |
        abs(exp(x1)-exp(x0))>0 |
        k<20){
    k=k+1
    x0=x1
    x1=lgamma((u+1)/2+k)-lgamma((u+1)/2)+k*log(v^2)-k*log(w)-lfactorial(2*k+1)
  }
  T2=k
  ###### SUM #######
  SumTimes=max(T1,T2)
  #  print(paste("SumTimes=",SumTimes))
  
  #k=seq(0,SumTimes,by=1)
  if(SumTimes<180){k=seq(0,SumTimes,by=1)} else {k=seq(0,SumTimes+300,by=1)}
  #  print(paste("k=",k[length(k)]))
  
  logM1=lgamma(u/2+k)-lgamma(u/2)+k*log(v^2)-k*log(w)-lfactorial(2*k)
  
  logM2=lgamma((u+1)/2+k)-lgamma((u+1)/2)+k*log(v^2)-k*log(w)-lfactorial(2*k+1)
  
  fox=exp(lgamma(u/2))*sum(exp(logM1))+exp(lgamma((u+1)/2))*(v/sqrt(w))*sum(exp(logM2))
  
  #fox=sum(c(exp(lgamma(u/2))*exp(logM1),exp(lgamma((u+1)/2))*(v/sqrt(w))*exp(logM2)))
  result=ifelse(log==F,fox,log(fox))
  return(result)
}




A_post_sq_expo = function(x, alpha=2, beta, gamma){
  #alpha=2 for weibull #this formula only holds for alpha=2!!!
  bb = 1/beta/2
  aa = bb*gamma 
  
  k0 = aa*sqrt(2*pi*bb)*(1-pnorm(-1*aa, 0, sqrt(bb)))+bb*(1-pweibull(aa, alpha, sqrt(bb*2)))
  k=1/k0
  k*(aa*sqrt(2*pi*bb)*dnorm(x-aa, 0, sqrt(bb))+bb*sign(x-aa)*dweibull(abs(x-aa), alpha, scale=sqrt(bb*2)))
}

