CB_DPM = function(MCMC, init, DPM.alpha = 1, type = c("rho", "rho&phi", "rho&phi_no_trans"), thin = 40, burn = 300){
  
  # thin the MCMC samples
  G = dim(MCMC$rho)[1]
  J = dim(MCMC$rho)[2]
  n.iter = dim(MCMC$rho)[3]
  iter.thin.seq = seq(burn+1, n.iter, by = thin)
  iter.thin.seq = iter.thin.seq[-1]
  
  #initialize the chain
  l = init; l.store = l;
  p.store = matrix(nrow = G, ncol = J)
  bc.store = matrix(nrow = G, ncol = J)
  bs.store = matrix(nrow = G, ncol = J) 
  z.mat.store = matrix(nrow = G, ncol = ifelse(type=="rho", J, J*3))
  
  #start updating 
  for(iter in iter.thin.seq){
    if(type == "rho"){
      input = MCMC$rho[ , , (iter-thin+1):iter]
      input = apply(input, c(1, 2), mean)
      # p.store = abind::abind(p.store, input, along = 3)
      z.mat = qnorm(input)
    }else if(type == "rho&phi"){
      input.phi = MCMC$phi[ , , (iter-thin+1):iter]
      input.p = MCMC$rho[ , , (iter-thin+1):iter]
      
      input.phi = input.p*input.phi
      input.bc = apply(input.phi, c(1, 2, 3), function(x){
        if(x==0){NA}else{(cos(x)+1)/2}})
      input.bs = apply(input.phi, c(1, 2, 3), function(x){
        if(x==0){NA}else{(sin(x)+1)/2}})
      
      input.p = apply(input.p, c(1, 2), mean)
      # p.store = abind::abind(p.store, input.p, along = 3)
      input.bc = apply(input.bc, c(1, 2), mean, na.rm = TRUE)
      input.bc[is.na(input.bc)] = 0.5
      # bc.store = abind::abind(bc.store, input.bc, along = 3)
      input.bs = apply(input.bs, c(1, 2), mean, na.rm = TRUE)
      input.bs[is.na(input.bs)] = 0.5
      # bs.store = abind::abind(bs.store, input.bs, along = 3)
      input = cbind(input.p, input.bc, input.bs)
      z.mat = qnorm(input)
    }else if(type == "rho&phi_no_trans"){
      input.phi = MCMC$phi[ , , (iter-thin+1):iter]
      input.p = MCMC$rho[ , , (iter-thin+1):iter]
      
      input.phi = input.p*input.phi
      input.bc = apply(input.phi, c(1, 2, 3), function(x){
        if(x==0){NA}else{cos(x)}})
      input.bs = apply(input.phi, c(1, 2, 3), function(x){
        if(x==0){NA}else{sin(x)}})
      
      input.p = apply(input.p, c(1, 2), mean)
      # p.store = abind::abind(p.store, input.p, along = 3)
      input.bc = apply(input.bc, c(1, 2), mean, na.rm = TRUE)
      input.bc[is.na(input.bc)] = 0
      # bc.store = abind::abind(bc.store, input.bc, along = 3)
      input.bs = apply(input.bs, c(1, 2), mean, na.rm = TRUE)
      input.bs[is.na(input.bs)] = 0
      # bs.store = abind::abind(bs.store, input.bs, along = 3)
      z.mat = cbind(qnorm(input.p), input.bc, input.bs)
    }
    #convert p to z, which is the value used for clustering.
    #the following two line are from the BayesMetaSeq package
    z.mat[, seq_len(J)][which(z.mat[, seq_len(J)]> 3)] <- 4 #XX: when z=3, pi^* is 0.9986, then pi is 0.9973 with positive beta. z=4 -> p=0.9999. This steps makes extreme cases more extrme?
    z.mat[, seq_len(J)][which(z.mat[, seq_len(J)]< -3)] <- -4
    z.mat.store = abind::abind(z.mat.store, z.mat, along = 3)
    
    new.l = sapply(seq_len(G), function(g){
      update_l_probs(z.mat, l, i=g, DPM.alpha)
    })
    l.store = cbind(l.store, new.l)
    l = new.l
  }
  
  return(list(l = l.store, 
              z.mat.store = z.mat.store,
              type = type))
}

update_l_probs = function(z.mat, l.vec, i, DPM.alpha){
  l.vec.no.i = l.vec[-i]
  if(sum(l.vec==l.vec[i])==1){
    l.pool = c(unique(l.vec.no.i), l.vec[i]) #when i forms a unique cluster, put l_i to the end as a possible new cluster (no need to generate new cluster label)
  }else{
    l.pool = c(unique(l.vec),max(unique(l.vec))+1) #append a new empty cluster at the end
  }
  
  z.i <- z.mat[i,]; J = ncol(z.mat); G = nrow(z.mat)
  prob.to.l = sapply(unique(l.vec.no.i), function(a.l){
    cluster.size = length(which(l.vec.no.i==a.l))
    cluster.ind = which(l.vec.no.i==a.l)
    if(cluster.size > 1){
      z.center <- apply(z.mat[cluster.ind,], 2, mean)
    }else{
      z.center <- z.mat[cluster.ind,]
    } #calculate cluster center
    
    post.mean = cluster.size/(cluster.size+1)*z.center
    post.sigma = diag((cluster.size+2)/(cluster.size+1),J)
    
    ifelse(cluster.size!=0,
           (cluster.size/(G-1+DPM.alpha))*mvtnorm::dmvnorm(x=z.i,mean=post.mean,sigma=post.sigma,log=FALSE),
           0)
  }, simplify=TRUE)
  
  prob.to.l.plus.1 <- (DPM.alpha/(G-1+DPM.alpha))*mvtnorm::dmvnorm(x=z.i,mean=rep(0,J),sigma=diag(2,J),log=FALSE) 
  #XX: for gene i, the probability of being from a new cluster.
  
  std.prob <- c(prob.to.l,prob.to.l.plus.1)/(sum(prob.to.l)+prob.to.l.plus.1)  #XX: normalize the probability
  sample(l.pool,size=1,prob=std.prob) #XX: draw cluster membership using updated probability
}



