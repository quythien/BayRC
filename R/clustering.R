CB_clustering = function(MCMC, k, main = "Rhythmicity clusters"){
  
  cluster.top = MCMC$l
  rownames(cluster.top) = paste('gene',seq_len(nrow(cluster.top)),sep=' ')
  
  # dist.mat <- mydist(y=cluster.top)
  an.out = ConsensusClusterPlus::ConsensusClusterPlus(d=t(cluster.top),
                                                      maxK=k,
                                                      reps=50,
                                                      pItem=0.8,
                                                      pFeature=1,
                                                      title= main,
                                                      clusterAlg="hc",
                                                      innerLinkage="ward.D2",
                                                      finalLinkage="ward.D2", 
                                                      distance="mydist",
                                                      seed=1262118388.71279,
                                                      plot=NULL)
  
}

copy.upper.lower <- function(m)
{
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}

mydist <- function(y) {
  diss.mat <- matrix(,nrow=nrow(y),ncol=nrow(y))
  rownames(diss.mat) <- colnames(diss.mat) <- rownames(y)
  diag(diss.mat) <- 1
  for (i in 1:nrow(diss.mat)) {
    for (j in i:ncol(diss.mat)){
      diss.mat[i,j] <- sum(y[i,]==y[j,])/ncol(y)
    }
  }
  diss.mat <- copy.upper.lower(diss.mat)
  return(1-diss.mat)
}