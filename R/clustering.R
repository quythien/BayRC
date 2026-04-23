#' Consensus clustering on MCMC rhythmicity indicator samples
#'
#' @title BayRC consensus clustering
#'
#' @description
#' Applies \code{ConsensusClusterPlus} to the posterior rhythmicity
#' indicator matrix (\code{MCMC$l}) to identify groups of genes with
#' similar circadian behaviour.  A custom dissimilarity based on the
#' proportion of MCMC samples where two genes share the same indicator
#' value is used (\code{mydist}).
#'
#' @param MCMC Named list; MCMC output containing an element \code{l}
#'   (G x K binary matrix of rhythmicity indicator samples).
#' @param k Integer; maximum number of clusters to evaluate.
#' @param main Character; title for the consensus clustering plot
#'   (default \code{"Rhythmicity clusters"}).
#'
#' @return The output of \code{ConsensusClusterPlus::ConsensusClusterPlus},
#'   a list of consensus matrices and cluster assignments for each k from
#'   2 to \code{k}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' clust <- CB_clustering(mcmc_out, k = 5)
#' }
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

#' Copy upper triangle to lower triangle of a matrix
#'
#' @description
#' Fills the lower triangle of matrix \code{m} with the transpose of the
#' upper triangle, making the matrix symmetric.
#'
#' @param m Square numeric matrix with upper triangle filled.
#'
#' @return Symmetric numeric matrix.
#'
#' @keywords internal
copy.upper.lower <- function(m)
{
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  m
}

#' Compute MCMC-based dissimilarity matrix for gene clustering
#'
#' @description
#' Computes a G x G dissimilarity matrix where entry (i, j) is one minus
#' the proportion of MCMC iterations in which genes i and j share the same
#' rhythmicity indicator value.  Used as the distance function inside
#' \code{CB_clustering}.
#'
#' @param y G x K binary matrix; posterior rhythmicity indicator samples.
#'
#' @return G x G dissimilarity matrix (values in \[0, 1\]) with row and
#'   column names from \code{rownames(y)}.
#'
#' @keywords internal
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