# my.TOJR.grid = expand.grid(Rhy1 = c(0, 1),
#                            Rhy2 = c(0, 1),
#                            Rhy3 = c(0, 1),
#                            dA2 = c(0), 
#                            dPhase2 = c(0, 2), 
#                            dM2 = c(0), 
#                            dA3 = c(0), 
#                            dPhase3 = c(0, 2), 
#                            dM3 = 0)
# my.TOJR.grid = my.TOJR.grid[c(1, 2, 4, 12, 24, 32), ]
# #1(1): arrhy; 
# #2(2): rhyI; 
# #3(4): rhyI&II, same phase, 
# #4(12): rhyI&II, different phase
# #5(24): rhyI&II&III, phase 1 and 2 the same, phase 3 different
# #6(32): rhyI&II&III, phase 1, 2, 3 all different. 
# my.TOJR.grid$n = rep(1, 6)

# history 
# 2023-06-08: change the function to allow heterogeneous time measurement error, also make the samples have paired sigma.t. 
# 2023-11-05: 
  # A new version of function to allow proportions of samples to have t_p and fixed t_p
  # phase of the data should be randomized across P

CBt_sim_data = function(P = 24, 
                        TOJR.grid = my.TOJR.grid[, 1:3],
                        dParam = list(g2 = my.TOJR.grid[, 4:6], 
                                      g3 = my.TOJR.grid[, 4:6]),
                        G = my.TOJR.grid$n,
                        n = 30, 
                        fixed_tp = TRUE, #if TRUE then the sigma.t will be the exact tp, however the signs will be evenly distributed above and below 0. 
                        sigma.t = 2, 
                        n.sigma.t = 16,
                        p.sigma.t = NULL, #percentage of samples with measurement error 
                        # one of n.sigma.t or p.sigma.t must to NULL
                        Params = list(A = 0.5, phi = 0, M = 5, sigma = 1),
                        TimePoints = NULL, #used to specify time and evenly distribute sample
                        PhaseClust = NULL){
  # 
  #   P = 24;
  #   TOJR.grid = my.TOJR.grid[, 1:3];
  #   dParam = list(g2 = my.TOJR.grid[, 4:6],
  #                 g3 = my.TOJR.grid[, 7:9])
  #   G = my.TOJR.grid$n;
  #   n = 30; Params = c(A = 1.2, phi = 6, M = 5, sigma = 1);
  #   sigma.t = 2;
  #   TimePoints = NULL; #used to specify time and evenly distribute sample
  #   PhaseClust = NULL
  if(is.null(n.sigma.t)&is.null(p.sigma.t)){
    stop("n.sigma.t and  p.sigma.t can not be both NULL. ")
  }else if((!is.null(n.sigma.t))&(!is.null(p.sigma.t))){
    stop("One of n.sigma.t and p.sigma.t must be NULL. ")
  }
  
  n.group = ncol(TOJR.grid)
  n.genes = sum(G)
  # stopifnot("The length of sample size vector n is not the same as the number of columns in TOJR.grid. " =
  #             length(n)==n.group)
  
  TOJR.grid.long = do.call(rbind.data.frame, lapply(1:nrow(TOJR.grid), function(a){
    matrix(replicate(G[a], TOJR.grid[a, ]), nrow =G[a], byrow = TRUE)
  }))
  colnames(TOJR.grid.long) = colnames(TOJR.grid)
  TOJR.dParam.long = lapply(seq_along(dParam), function(j){
    a.TOJR.dParam.long = do.call(rbind.data.frame, 
                                 lapply(1:nrow(dParam[[j]]), function(a){
                                   matrix(replicate(G[a], dParam[[j]][a, ]), nrow =G[a], byrow = TRUE)
                                 }))
    colnames(a.TOJR.dParam.long) = colnames(dParam[[j]])
    a.TOJR.dParam.long = apply(a.TOJR.dParam.long, c(1, 2), as.numeric)
    a.TOJR.dParam.long = as.data.frame(a.TOJR.dParam.long)
    return(a.TOJR.dParam.long)
  })
  g1 = as.data.frame(apply(TOJR.dParam.long[[1]], c(1, 2), function(x){0}))
  TOJR.dParam.long = append(TOJR.dParam.long, list(g1), after=0)
  
  sigma = lapply(1:n.group, function(j){
    a.sigma = rep(Params["sigma"], n.genes)
  })
  
  # phi1 = rep(Params["phi"], n.genes)
  phi1 = runif(n.genes, 0, P)
  A1 = rep(Params["A"], n.genes)
  M1 = rep(Params["M"], n.genes)
  
  A = lapply(1:n.group, function(j){
    a.A = A1+TOJR.dParam.long[[j]]$dA
    a.A = as.numeric(TOJR.grid.long[, j])*a.A
  })
  phase = lapply(1:n.group, function(j){
    # a.phi = rnorm(n.genes, phi1, 0.3)+TOJR.dParam.long[[j]]$dPhase
    # # add a little vairance to avoid numerical issue
    a.phi = phi1 +TOJR.dParam.long[[j]]$dPhase
    #do not add variance ...
  })
  M = lapply(1:n.group, function(j){
    a.M = M1+TOJR.dParam.long[[j]]$dM
  })
  
  # simulate t for matched non-missing samples  
  if(is.null(TimePoints)){
    t = runif(n, 0, P) #the observed time
    #the miss measurement of time: MCT-CT
    ## decide which samples should have pt \neq 0
    if(is.null(n.sigma.t)){
      t.p.ind = rbinom(n, 1, p.sigma.t)
    }else{
      t.p.ind1 = sample.int(n, n.sigma.t)
      t.p.ind = rep(0, n)
      t.p.ind[t.p.ind1] = 1
    }
    ## assign a t_p value for samples with pt \neq 0
    if(fixed_tp){
      #further randomize half of t.p.ind1 to be smaller than 0
      n.t.p = sum(t.p.ind)
      n.t.p.minus = floor(n.t.p/2)
      t.p.minus.ind = sample.int(n.t.p, n.t.p.minus)
      t.p.ind[t.p.ind==1][t.p.minus.ind] = -1
      t.p = sigma.t * t.p.ind
    }else{
      t.p = rnorm(n, 0, sigma.t)*t.p.ind
    }
  }
  #change history: 
  #20231001: sigma.t should be the same in one simulation. 
  
  dat = lapply(1:n.group, function(j){
    t1 = cos((t+t.p)*2*pi/P)
    t2 = sin((t+t.p)*2*pi/P)
    beta1 = A[[j]]*cos(phase[[j]]*2*pi/P)
    beta2 = A[[j]]*sin(phase[[j]]*2*pi/P)
    x = data.frame(x1 = 1, t1, t2)
    beta = data.frame(M = M[[j]], beta1, beta2)
    epsilon = matrix(rnorm(n*n.genes, sd = rep(sigma[[j]], each = n)), ncol = n, byrow = TRUE)
    Y = as.matrix(beta) %*% t(as.matrix(x))+epsilon
    x$time = t #observed time
    x$time.p = t.p #MCT
    x$P = P
    beta$A = A[[j]]
    beta$phase = phase[[j]]
    beta$sigma = sigma[[j]]
    return(list(x = x,
                beta = beta,
                dat = Y))
  })
  return(list(data = dat,
              TOJR = cbind(TOJR.grid.long, TOJR.dParam.long)))
}

# previous: without heterogeneous TME ----------
# CBt_sim_data = function(P = 24, 
#                         TOJR.grid = my.TOJR.grid[, 1:3],
#                         dParam = list(g2 = my.TOJR.grid[, 4:6], 
#                                       g3 = my.TOJR.grid[, 4:6]),
#                         sigma.t = 2,
#                         G = my.TOJR.grid$n,
#                         n = rep(20, 3), Params = list(A = 0.5, phi = 0, M = 5, sigma = 1),
#                         TimePoints = NULL, #used to specify time and evenly distribute sample
#                         PhaseClust = NULL){
#   # 
#   #   P = 24;
#   #   TOJR.grid = my.TOJR.grid[, 1:3];
#   #   dParam = list(g2 = my.TOJR.grid[, 4:6],
#   #                 g3 = my.TOJR.grid[, 7:9])
#   #   G = my.TOJR.grid$n;
#   #   n = 30; Params = c(A = 1.2, phi = 6, M = 5, sigma = 1);
#   #   sigma.t = 2;
#   #   TimePoints = NULL; #used to specify time and evenly distribute sample
#   #   PhaseClust = NULL
#   
#   n.group = ncol(TOJR.grid)
#   n.genes = sum(G)
#   # stopifnot("The length of sample size vector n is not the same as the number of columns in TOJR.grid. " =
#   #             length(n)==n.group)
#   
#   TOJR.grid.long = do.call(rbind.data.frame, lapply(1:nrow(TOJR.grid), function(a){
#     matrix(replicate(G[a], TOJR.grid[a, ]), nrow =G[a], byrow = TRUE)
#   }))
#   colnames(TOJR.grid.long) = colnames(TOJR.grid)
#   TOJR.dParam.long = lapply(seq_along(dParam), function(j){
#     a.TOJR.dParam.long = do.call(rbind.data.frame, 
#                                  lapply(1:nrow(dParam[[j]]), function(a){
#                                    matrix(replicate(G[a], dParam[[j]][a, ]), nrow =G[a], byrow = TRUE)
#                                  }))
#     colnames(a.TOJR.dParam.long) = colnames(dParam[[j]])
#     a.TOJR.dParam.long = apply(a.TOJR.dParam.long, c(1, 2), as.numeric)
#     a.TOJR.dParam.long = as.data.frame(a.TOJR.dParam.long)
#     return(a.TOJR.dParam.long)
#   })
#   g1 = as.data.frame(apply(TOJR.dParam.long[[1]], c(1, 2), function(x){0}))
#   TOJR.dParam.long = append(TOJR.dParam.long, list(g1), after=0)
#   
#   sigma = lapply(1:n.group, function(j){
#     a.sigma = rep(Params$sigma, n.genes)
#   })
#   phi1 = rep_len(Params$phi, n.genes)
#   A1 = rep(Params$A, n.genes)
#   M1 = rep(Params$M, n.genes)
#   
#   A = lapply(1:n.group, function(j){
#     a.A = A1+TOJR.dParam.long[[j]]$dA
#     a.A = as.numeric(TOJR.grid.long[, j])*a.A
#   })
#   phase = lapply(1:n.group, function(j){
#     # a.phi = rnorm(n.genes, phi1, 0.3)+TOJR.dParam.long[[j]]$dPhase
#     # # add a little vairance to avoid numerical issue
#     a.phi = phi1 +TOJR.dParam.long[[j]]$dPhase
#     #do not add variance ...
#   })
#   M = lapply(1:n.group, function(j){
#     a.M = M1+TOJR.dParam.long[[j]]$dM
#   })
#   
#   # simulate t for matched non-missing samples  
#   if(is.null(TimePoints)){
#     t = runif(n, 0, P)
#     #the miss measurement of time/MCT-CT
#     t.p = rnorm(n, 0, sigma.t)
#   }
#   
#   dat = lapply(1:n.group, function(j){
#     t1 = cos((t+t.p)*2*pi/P)
#     t2 = sin((t+t.p)*2*pi/P)
#     beta1 = A[[j]]*cos(phase[[j]]*2*pi/P)
#     beta2 = A[[j]]*sin(phase[[j]]*2*pi/P)
#     x = data.frame(x1 = 1, t1, t2)
#     beta = data.frame(M = M[[j]], beta1, beta2)
#     epsilon = matrix(rnorm(n*n.genes, sd = rep(sigma[[j]], each = n)), ncol = n, byrow = TRUE)
#     Y = as.matrix(beta) %*% t(as.matrix(x))+epsilon
#     x$time = t 
#     x$time.p = t.p
#     x$P = P
#     beta$A = A[[j]]
#     beta$phase = phase[[j]]
#     beta$sigma = sigma[[j]]
#     return(list(x = x,
#                 beta = beta,
#                 dat = Y))
#   })
#   return(list(data = dat,
#               TOJR = cbind(TOJR.grid.long, TOJR.dParam.long)))
# }