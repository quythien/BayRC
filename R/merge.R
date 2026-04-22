##' The \code{merge} is function to merge multiple MCMCout matrices by
##' matching genes/orthologs.
##' @title Merge multiple MCMCout matrices by orthologs
##' @param mcmc.list: a list of MCMC output matrices.
##' @param species: a vector specie names of same length as mcmc.list indicating
##' the species of each study.
##' @param ortholog.db: a data.frame/matrix match orthologs between species.
##' Column names should be consistent with the species input. If not provided,
##' datasets will be merges without orthologs matching.
##' @param reference: the index of the reference MCMC matrix. Merged
##' list will be named using the rownames of this matrix.
##' @param uniqG: TRUE: only keep the gene with greatest posterior DE signal
##' when multiple matches provided in the orhthologs file. FALSE: keep duplicated
##' genes when multiple matches exist. default is TRUE.

##' @return an merged list of multiple MCMCout datasets (with same number of
##' rows and rownames)
##' @export
##' @examples
##' \dontrun{
##' data(hb)
##' summaryDE = indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' hb_pData = summaryDE[,c(3,1)]
##' hb_MCMCout = bayes(hb_pData, seed=12345)
##' data(hs)
##' summaryDE = indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' hs_pData = summaryDE[,c(3,1)]
##' hs_MCMCout = bayes(hs_pData, seed=12345)
##' data(ht)
##' summaryDE <- indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' ht_pData <- summaryDE[,c(3,1)]
##' ht_MCMCout <- bayes(ht_pData, seed=12345)
##' data(ha)
##' summaryDE = indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' ha_pData = summaryDE[,c(3,1)]
##' ha_MCMCout = bayes(ha_pData, seed=12345)
##' data(hi)
##' summaryDE = indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' hi_pData = summaryDE[,c(3,1)]
##' hi_MCMCout = bayes(hi_pData, seed=12345)
##' data(hl)
##' summaryDE = indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' hl_pData = summaryDE[,c(3,1)]
##' hl_MCMCout = bayes(hl_pData, seed=12345)
##' data(hb)
##' summaryDE = indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' mb_pData = summaryDE[,c(3,1)]
##' mb_MCMCout = bayes(mb_pData, seed=12345)
##' data(hs)
##' summaryDE = indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' ms_pData = summaryDE[,c(3,1)]
##' ms_MCMCout = bayes(ms_pData, seed=12345)
##' data(ht)
##' summaryDE <- indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' mt_pData <- summaryDE[,c(3,1)]
##' mt_MCMCout <- bayes(mt_pData, seed=12345)
##' data(ma)
##' summaryDE = indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' ma_pData = summaryDE[,c(3,1)]
##' ma_MCMCout = bayes(ma_pData, seed=12345)
##' data(hi)
##' summaryDE = indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' mi_pData = summaryDE[,c(3,1)]
##' mi_MCMCout = bayes(mi_pData, seed=12345)
##' data(ml)
##' summaryDE = indDE(data=data,group=as.factor(group),data.type="microarray",
##'                   case.label="2", ctrl.label="1")
##' ml_pData = summaryDE[,c(3,1)]
##' ml_MCMCout = bayes(ml_pData, seed=12345)
##'
##' mcmc.list <- list(hb_MCMCout,hs_MCMCout,ht_MCMCout,
##' ha_MCMCout,hi_MCMCout,hl_MCMCout,
##' mb_MCMCout,ms_MCMCout,mt_MCMCout,
##' ma_MCMCout,mi_MCMCout,ml_MCMCout)
##' data(hm_orth)
##' species <- c(rep("human",6), rep("mouse",6))
##' mcmc.merge.list <- merge(mcmc.list,species = species,
##' ortholog.db = hm_orth, reference=1,uniqG=T)
##' }

merge <- function(mcmc.list,species,ortholog.db = NULL,
                  reference=1,uniqG=T){
  if(is.null(ortholog.db) | length(unique(species)) == 1){
    print("Only one species, merge datasets without ortholog match")
    M <- length(mcmc.list)
    mcmc.merge.list <- vector("list",M)

    gene.list <- lapply(mcmc.list,rownames)
    DEgene.list <- lapply(mcmc.list,function(x) rownames(x)[attr(x,"DEindex")])

    int_genes <- Reduce(intersect, gene.list)

    for(m in 1:M){
      mcmc.merge.list[[m]] <- mcmc.list[[m]][int_genes,]
      DEgenes <- int_genes[which(int_genes %in% DEgene.list[[m]])]
      DEindex <- which(rownames(mcmc.merge.list[[m]])%in% DEgenes)
      attr(mcmc.merge.list[[m]],"DEindex") <- DEindex
    }
    return(mcmc.merge.list)
  }else{
    if(!all(species %in% colnames(ortholog.db))){
      stop('Found species without corresponding orthologs in the ortholog.db provides.')
    }

    if(length(species) != length(mcmc.list)){
      stop('The number of species does not match with the number of mcmc matrices.')
    }
    M <- length(mcmc.list)
    mcmc.merge.list <- DEgene.merge.list <- vector("list",M)

    gene.list <- lapply(mcmc.list,rownames)
    DEgene.list <- lapply(mcmc.list,function(x) rownames(x)[attr(x,"DEindex")])

    match_gene <- orthMatch(gene.list,species,ortholog.db)
    ref_gene <- match_gene[[reference]]
    ## match_gene and ref_gene of same dimension
    for(m in 1:M){
      mcmc.merge.list[[m]] <- mcmc.list[[m]][match_gene[[m]],]
      rownames(mcmc.merge.list[[m]]) <- ref_gene
      DEgenes <- ref_gene[which(match_gene[[m]] %in% DEgene.list[[m]])]
      DEindex <- which(rownames(mcmc.merge.list[[m]])%in% DEgenes)
      attr(mcmc.merge.list[[m]],"DEindex") <- DEindex
    }
    if(uniqG==T){
      mcmc.merge.list.uniq = lapply(mcmc.merge.list, function(mat){
        geneSplit = split(1:nrow(mat),row.names(mat))
        uniq.index = sapply(geneSplit,function(g){
          if(length(g)==1){
            return(g)
          }else{
            DEmean = abs(apply(mat[g,],1,mean))
            selected.row = which.max(DEmean)
            return(g[selected.row])
          }
        })
        uniq.mat = mat[uniq.index,]
        attr(uniq.mat,"DEindex") = which(!is.na(match(uniq.index,attr(mat,"DEindex"))))
        return(uniq.mat)
      })
      return(mcmc.merge.list.uniq)
    }else{
      return(mcmc.merge.list)
    }
  }
}


orthMatch <- function(gene.list,species,ortholog.db){
  M <- length(gene.list)
  index.out <- gene.out <- vector("list",M)
  for(m in 1:M){
    spec <- species[m]
    gene <- gene.list[[m]]
    index.out[[m]] <- which(ortholog.db[,spec] %in% gene)
  }
    common.index <- Reduce(intersect,index.out)#genes in othtolog that appeared in all mcmclists with corresponding species

  for(m in 1:M){
    spec <- species[m]
    gene.out[[m]] <- as.character(ortholog.db[common.index,spec])#otherwise factor used
  }

  return(gene.out)
}
