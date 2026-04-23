#' Merge multiple MCMC outputs by matching orthologous genes
#'
#' @title Merge MCMC output lists by orthologs
#'
#' @description
#' Intersects gene sets across a list of MCMC outputs from different
#' species or conditions.  When \code{ortholog.db} is supplied and multiple
#' species are present, the ortholog database maps each species' genes to a
#' common reference gene space (defined by the \code{reference} dataset).
#' When datasets share a single species (or no ortholog database is
#' provided), the function simply takes the intersection of gene names.
#' If \code{uniqG = TRUE}, duplicate ortholog matches are resolved by
#' keeping the entry with the highest mean absolute posterior signal.
#'
#' @param mcmc.list A list of MCMC output matrices or lists, one per
#'   condition or study.  Row names must be gene identifiers.
#' @param species Character vector of the same length as \code{mcmc.list};
#'   species label for each dataset.
#' @param ortholog.db Data.frame or matrix with one column per species
#'   (column names matching \code{species}); each row maps orthologous
#'   genes across species.  If \code{NULL}, datasets are merged by gene
#'   name intersection only (default \code{NULL}).
#' @param reference Integer; index of the reference dataset whose gene
#'   names define the merged row names (default 1).
#' @param uniqG Logical; if \code{TRUE} (default), keep only the
#'   highest-signal match when multiple orthologs map to the same reference
#'   gene.
#'
#' @return A list of the same length as \code{mcmc.list}, each element
#'   row-subsetted and reordered to the common gene set, with row names
#'   from the reference dataset.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' merged <- merge_mcmc(list(human_res, mouse_res),
#'                      species      = c("human", "mouse"),
#'                      ortholog.db  = hm_orth,
#'                      reference    = 1)
#' }
merge_mcmc <- function(mcmc.list,species,ortholog.db = NULL,
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


#' Match gene lists to a common ortholog set
#'
#' @description
#' For each species in \code{gene.list}, finds the rows in
#' \code{ortholog.db} that contain at least one gene from that species,
#' then takes the intersection across all species to return a common set
#' of orthologous genes.
#'
#' @param gene.list List of character vectors; one per species, each
#'   containing gene identifiers present in that species' dataset.
#' @param species Character vector; species label for each element of
#'   \code{gene.list}.
#' @param ortholog.db Data.frame; ortholog database (see \code{merge_mcmc}).
#'
#' @return A list of character vectors, one per species, giving the gene
#'   identifiers (in the original species' namespace) that correspond to
#'   the common ortholog set.
#'
#' @keywords internal
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
