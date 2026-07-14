#### Thien
##########################
###### Check functions ###
##########################

# check.rawData <- function(data){
#   G <- nrow(data)
#   N <- ncol(data)
#   if(!is.numeric(data)) {
#     stop("expression data not numeric")
#   }
#   if(is.null(row.names(data))) {
#     stop("gene symbol missing")
#   }
#   if(N <= 3) {
#     stop("too few samples")
#   }
#   if(sum(duplicated(row.names(data)))>0){
#     stop("duplicate gene symbols")
#   }
# }

# check.groupData <- function(group){
#   l <- nlevels(group)
#   if( l != 2 ) {
#     stop("not a two-class comparison")
#   }
# }
# 
# check.pData <- function(pData){
#   G <- nrow(pData)
#   if(!is.numeric(pData)) {
#     stop("pvalue data not numeric")
#   }
#   if(is.null(row.names(pData))) {
#     stop("gene symbol missing")
#   }
#   if(ncol(pData) <2) {
#     stop("missing either p-value or effect size")
#   }
# }


# check.compatibility <- function(data, group, case.label, ctrl.label){
#   G <- nrow(data)
#   N1 <- ncol(data)
#   N2 <- length(group)
#   if(N1 != N2) {
#     stop("expression data and class label have unmatched sample size")
#   }
#   if(!all(group %in% c(case.label,ctrl.label))){
#     stop("including class labels other than the case and control")
#   }
#   if(sum(group==case.label) <= 1 ||sum(group==ctrl.label) <= 1){
#     stop("not enough samples in either case or control group")
#   }
# }

##########################
###### BayesP part ###
##########################

# PtoZ <- function(p2, lfc) {
#   sgn <- sign(lfc)
#   z <- ifelse(sgn>0, qnorm(p2/2), qnorm(p2/2,lower.tail = F))
#   return(z)
# }
# 
# SelectGamma <- function(p){
#   ## Gamma is DE proportion  = 1-pi0
#   m <- length(p)
#   lambda <- seq(0,0.95,by=0.01)
#   pi0 <- sapply(lambda, function(x) sum(p>x)/(m*(1-x))  )
#   
#   # fit a natural cubic spline
#   library(splines)
#   dat <- data.frame(pi0=pi0, lambda=lambda)
#   lfit <- lm(pi0 ~ ns(lambda, df = 3), data=dat)
#   pi0hat <- predict(lfit, data.frame(lambda=1))
#   gamma <- 1 - pi0hat
#   return(gamma)
# }

##################################
###### ACS scores ################
##################################
# 
# CSY <- function(d1,d2,index1){
#   Sens <- calcSensC(d1[index1,],d2[index1,])
#   Spec <- calcSpecC(d1[-index1,],d2[-index1,])
#   
#   CS <- Sens + Spec - 1
#   if(is.nan(CS)) {
#     CS <- 0
#   }
#   return(CS)
# }
# 
# CSF <- function(d1,d2,index1,index2){
#   Sens <- calcSensC(d1[index1,],d2[index1,])
#   Prec <- calcPrecC(d1[index2,],d2[index2,])
#   
#   CS <- (2*Sens*Prec)/(Sens+Prec)
#   if(is.nan(CS)) {
#     CS <- 0
#   }
#   return(CS)
# }
# 
# CSG <- function(d1,d2,index1,index2){
#   Sens <- calcSensC(d1[index1,],d2[index1,])
#   Spec <- calcSpecC(d1[-index1,],d2[-index1,])
#   
#   CS <- sqrt(Sens*Spec)
#   if(is.nan(CS)) {
#     CS <- 0
#   }
#   return(CS)
# }
# 
# 
# ECSY <- function(d1,d2){
#   ESens <- calcESensC(d1,d2)
#   ESpec <- calcESpecC(d1,d2)
#   
#   ECS <- ESens + ESpec - 1
#   if(is.nan(ECS)) {
#     ECS <- 0
#   }
#   return(ECS)
# }
# 
# ECSF <- function(d1,d2){
#   ESens <- calcESensC(d1,d2)
#   EPrec <- calcEPrecC(d1,d2)
#   
#   ECS <- (2*ESens*EPrec)/(ESens+EPrec)
#   if(is.nan(ECS)) {
#     ECS <- 0
#   }
#   return(ECS)
# }
# 
# ECSG <- function(d1,d2){
#   ESens <- calcESensC(d1,d2)
#   ESpec <- calcESpecC(d1,d2)
#   
#   ECS <- sqrt(ESens*ESpec)
#   if(is.nan(ECS)) {
#     ECS <- 0
#   }
#   return(ECS)
# }
# 
# permCSY <- function(d1,d2){
#   Sens <- calcSensC(d1,d2)
#   Spec <- calcSpecC(d1,d2)
#   
#   permCS <- Sens + Spec - 1
#   if(is.nan(permCS)) {
#     permCS <- 0
#   }
#   return(permCS)
# }
# 
# permCSF <- function(d1,d2){
#   Sens <- calcSensC(d1,d2)
#   Prec <- calcPrecC(d1,d2)
#   
#   permCS <- (2*Sens*Prec)/(Sens+Prec)
#   if(is.nan(permCS)) {
#     permCS <- 0
#   }
#   return(permCS)
# }
# 
# permCSG <- function(d1,d2){
#   Sens <- calcSensC(d1,d2)
#   Spec <- calcSpecC(d1,d2)
#   
#   permCS <- sqrt(Sens*Spec)
#   if(is.nan(permCS)) {
#     permCS <- 0
#   }
#   return(permCS)
# }


# CS <- function(dat1,dat2,deIndex1,deIndex2,measure="Fmeasure"){
#   ## same for both global and pathway
#   if(measure=="youden"){
#     CS <- CSY(dat1,dat2,deIndex1)
#   } else if(measure=="Fmeasure"){
#     CS <- CSF(dat1,dat2,deIndex1,deIndex2)
#   } else if(measure=="geo.mean"){
#     CS <- CSG(dat1,dat2,deIndex1)
#   }
#   return(CS)
# }
# 
# 
# ECS <- function(dat1,dat2,measure="Fmeasure"){
#   ## same for both global and pathway
#   if(measure=="youden"){
#     ECS <- ECSY(dat1,dat2)
#   } else if(measure=="Fmeasure"){
#     ECS <- ECSF(dat1,dat2)
#   } else if(measure=="geo.mean"){
#     ECS <- ECSG(dat1,dat2)
#   }
#   return(ECS)
# }
# 
# permCS <- function(dat1,dat2,measure="Fmeasure"){
#   ## same for both global and pathway
#   if(measure=="youden"){
#     pemrCS <- permCSY(dat1,dat2)
#   } else if(measure=="Fmeasure"){
#     permCS <- permCSF(dat1,dat2)
#   } else if(measure=="geo.mean"){
#     permCS <- permCSG(dat1,dat2)
#   }
#   return(permCS)
# }

##################################
#### Global (expected value from marginal,
#### permutate genes to get p-value)
##################################

# perm_global <- function(dat1,dat2,measure="Fmeasure",B){
#   G <- nrow(dat1) 
#   
#   #rawCS <- CS(dat1,dat2,deIndex1,deIndex2,measure)
#   #expCS <- ECS(dat1,dat2,measure)
#   #rawDS <- DS(dat1,dat2,deIndex1,deIndex2,measure)
#   #expDS <- EDS(dat1,dat2,measure)
#   
#   out <- matrix(NA,B,2)
#   colnames(out) <- c("permCS","permECS")
#   
#   for(b in 1:B){
#     dat1perm <- dat1[sample(1:G,G,replace = F),]
#     dat2perm <- dat2[sample(1:G,G,replace = F),]
#     out[b,"permCS"] <- permCS(dat1perm,dat2perm,measure)
#     out[b,"permECS"] <- ECS(dat1perm,dat2perm,measure)
#   }
#   return(out)
# }
# add absolute value here :
# ACS_global <- function(dat1,dat2,deIndex1,deIndex2,
#                        measure="Fmeasure"){
#   cs <- CS(dat1,dat2,deIndex1,deIndex2,measure)
#   ecs <- ECS(dat1,dat2,measure)
#   acs <- (cs - ecs)/(1-ecs)
#   return(acs)
# }
# # add absolute value here :
# pACS_global <- function(dat1,dat2,deIndex1,deIndex2,
#                         measure="Fmeasure",acs,permOut){
#   
#   permcs <- permOut[,"permCS"]
#   permecs <- permOut[,"permECS"]
#   permacs <- (permcs - permecs)/(1-permecs)
#   
#   p_acs <- (sum(permacs>=acs) + 1)/(length(permacs)+1)
#   return(p_acs)
# }

##################################
#### Pathway (expected value from global,
#### permute genes to get p-value)
##################################

# margin_pathway <-  function(dat1,dat2,
#                             select.pathway.list,
#                             measure="Fmeasure"){
#   select.pathways <- names(select.pathway.list)
#   
#   dat1 = as.matrix(dat1)
#   dat2 = as.matrix(dat2)
#   # Add
#   data_genes <- rownames(dat1)
#   pathway.size <- sapply(select.pathway.list,function(x) {
#     length(intersect(data_genes,x))})
#   K <- length(select.pathways)
#   G <- nrow(dat1)
#   out <- matrix(NA,nrow=K,ncol=1)
#   rownames(out) <- select.pathways
#   colnames(out) <- c("ECS")
#   
#   R <- 20 ##fairly enough
#   
#   for(k in 1:K){
#     #print(k)
#     ecsk <- edsk <- rep(NA,R)
#     pathsizek <- pathway.size[k]
#     for(j in 1:R){
#       index <- sample(1:G,pathsizek,replace=F)
#       dat1.select <- dat1[index,]
#       dat2.select <- dat2[index,]
#       ecsk[j] <- ECS(dat1.select,dat2.select,measure)
#     }
#     out[k,"ECS"] <- mean(ecsk)
#   }
#   
#   return(out)
#   
# }

# perm_pathway <- function(dat1,dat2,
#                          select.pathway.list,
#                          measure="Fmeasure",B,parallel=F,n.cores=4){
#   
#   select.pathways <- names(select.pathway.list)
#   data_genes <- attr(dat1, "symbols")
#   pathway.size <- sapply(select.pathway.list,function(x) {
#     length(intersect(data_genes,x))})
#   K <- length(select.pathways)
#   G <- nrow(dat1)
#   
#   #out <- array(1,dim=c(B,K,4),dimnames=
#   #list(1:B,select.pathways,
#   #c("permCS","permECS","permDS","permEDS")))
#   
#   #rawCS <- CS(dat1,dat2,deIndex1,deIndex2,measure)
#   #expCS <- ECS(dat1,dat2,measure)
#   #rawDS <- DS(dat1,dat2,deIndex1,deIndex2,measure)
#   #expDS <- EDS(dat1,dat2,measure)
#   
#   out <- array(1,dim=c(B,K,1),dimnames=
#                  list(1:B,select.pathways,c("permCS")))
#   
#   for(k in 1:K){
#     #print(k)
#     pathsizek <- pathway.size[k]
#     if(parallel == T){
#       require(parallel)
#       permFunc = function(b){
#         dat1perm <- dat1[sample(1:G,G,replace = F),]
#         dat2perm <- dat2[sample(1:G,G,replace = F),]
#         index <- sample(1:G,pathsizek,replace=F)
#         dat1perm.select <- dat1perm[index,]
#         dat2perm.select <- dat2perm[index,]
#         permCS_res <- permCS(dat1perm.select,dat2perm.select,measure)
#         return(list(permCS_res = permCS_res))
#       }
#       out.ls = mclapply(1:B, permFunc, mc.cores = n.cores)
#       for(b in 1:B){
#         out[b,k,"permCS"] <- out.ls[[b]]$permCS_res
#       }
#     }else{
#       for(b in 1:B){
#         dat1perm <- dat1[sample(1:G,G,replace = F),]
#         dat2perm <- dat2[sample(1:G,G,replace = F),]
#         index <- sample(1:G,pathsizek,replace=F)
#         dat1perm.select <- dat1perm[index,]
#         dat2perm.select <- dat2perm[index,]
#         
#         out[b,k,"permCS"] <- permCS(as.matrix(dat1perm.select),as.matrix(dat2perm.select,measure))
#         #out[b,k,"permECS"] <- ECS(dat1perm.select,dat2perm.select,measure)
#         #out[b,k,"permEDS"] <- EDS(dat1perm.select,dat2perm.select,measure)
#       }
#     }
#   }
#   return(out)
# }
# Add absolute value here 
# ACS_pathway <- function(dat1,dat2,deIndex1,deIndex2,
#                         select.pathway.list,
#                         measure="Fmeasure",marginOut){
#   
#   select.pathways <- names(select.pathway.list)
#   data_genes <- attr(dat1, "symbols")
#   pathway.size <- sapply(select.pathway.list,function(x) {
#     length(intersect(data_genes,x))})
#   K <- length(select.pathways)
#   G <- nrow(dat1)
#   
#   acs <- rep(NA,K)
#   names(acs) <- select.pathways
#   # The issue is here 
#   for(k in 1:K){
#     # path_genek <- select.pathway.list[[k]]
#     # genek <- intersect(path_genek,data_genes)
#     # 
#     # # Subset without including NA values
#     # valid_indices <- match(genek, data_genes)
#     # valid_indices <- valid_indices[!is.na(valid_indices)]  # Remove NA values
#     # 
#     # dat1_k <- dat1[valid_indices,]
#     # dat2_k <- dat2[valid_indices,]
#     # 
#     # # Correctly assign the RHYindex and symbols attributes
#     # attr(dat1_k, "RHYindex") <- attr(dat1, "RHYindex")[valid_indices]
#     # attr(dat2_k, "RHYindex") <- attr(dat2, "RHYindex")[valid_indices]
#     # 
#     # attr(dat1_k, "symbols") <- attr(dat1, "symbols")[valid_indices]
#     # attr(dat2_k, "symbols") <- attr(dat2, "symbols")[valid_indices]
#     # Pathway genes and intersection with data genes
#     path_genek <- select.pathway.list[[k]]
#     genek <- intersect(path_genek, data_genes)
#     
#     # Find valid indices for both datasets, aligning by gene symbols
#     valid_indices1 <- match(genek, attr(dat1, "symbols"))
#     valid_indices2 <- match(genek, attr(dat2, "symbols"))
#     
#     # Remove NA values (i.e., genes not found in either dataset)
#     valid_indices1 <- valid_indices1[!is.na(valid_indices1)]
#     valid_indices2 <- valid_indices2[!is.na(valid_indices2)]
#     
#     # Now align dat1 and dat2 based on the same genes
#     aligned_genes <- intersect(attr(dat1, "symbols")[valid_indices1], attr(dat2, "symbols")[valid_indices2])
#     
#     # Get indices of the aligned genes in each dataset
#     aligned_indices1 <- match(aligned_genes, attr(dat1, "symbols"))
#     aligned_indices2 <- match(aligned_genes, attr(dat2, "symbols"))
#     
#     # Subset the data with the aligned indices
#     dat1_k <- dat1[aligned_indices1, ]
#     dat2_k <- dat2[aligned_indices2, ]
#     
#     # Correctly assign the RHYindex and symbols attributes, using aligned indices
#     attr(dat1_k, "RHYindex") <- attr(dat1, "RHYindex")[aligned_indices1]
#     attr(dat2_k, "RHYindex") <- attr(dat2, "RHYindex")[aligned_indices2]
#     
#     attr(dat1_k, "symbols") <- attr(dat1, "symbols")[aligned_indices1]
#     attr(dat2_k, "symbols") <- attr(dat2, "symbols")[aligned_indices2]
#     
#     # Now dat1_k and dat2_k should have matching rows for the same genes
#     
#     
#     if(length(intersect(names(deIndex1), genek))<=3 ){
#       deIndex1_k <- 1:nrow(dat1_k)
#     } else {
#       deIndex1_k <- which(attr(dat1_k, "symbols")%in%intersect(names(deIndex1), genek))
#     }
#     
#     if(length(intersect(names(deIndex2), genek))<=3 ){
#       deIndex2_k <- 1:nrow(dat2_k)
#     } else {
#       deIndex2_k <- which(attr(dat2_k, "symbols")%in%intersect(names(deIndex2), genek))
#     }
#     
#     cs <- CS(as.matrix(dat1_k),as.matrix(dat2_k),as.matrix(deIndex1_k),as.matrix(deIndex2_k),measure)
#     ecs <- marginOut[k,"ECS"]
#     acs[k] <- (cs - ecs)/(1-ecs)
#   }
#   return(acs)
# }
# # Add absolute value here 
# pACS_pathway <- function(dat1,dat2,deIndex1,deIndex2,
#                          select.pathway.list,
#                          measure="Fmeasure",acs,permOut,marginOut){
#   
#   select.pathways <- names(select.pathway.list)
#   K <- length(select.pathways)
#   
#   p_acs <- rep(NA,K)
#   names(p_acs) <- select.pathways
#   
#   for(k in 1:K){
#     
#     permcs <- permOut[,k,"permCS"]
#     ecs <- marginOut[k,"ECS"]
#     #permecs <- permOut[,k,"permECS"]
#     permacs <- (permcs - ecs)/(1-ecs)
#     
#     p_acs[k] <- (sum(permacs>=acs[k]) + 1)/(length(permacs)+1)
#   }
#   
#   return(p_acs)
#   
# }


##########################
##Pathway enrich analysis#
##########################

gsa.fisher <- function(x, background, pathway) {
  ####x is the list of query genes
  ####backgroud is a list of background genes that query genes from
  ####pathway is a list of different pathway genes
  count_table<-matrix(0,2,2)
  x<-toupper(x)
  background<-toupper(background)
  index<-which(toupper(background) %in% toupper(x)==FALSE)
  background_non_gene_list<-background[index]
  x<-toupper(x)
  pathway<-lapply(pathway,function(x) intersect(toupper(background),toupper(x)))
  get.fisher <- function(path) {
    res <- NA
    ####in the gene list and in the pathway
    count_table[1,1]<-sum(x %in% path)
    #count_table[1,1]<-sum(is.na(charmatch(x,path))==0)
    ####in the gene list but not in the pathway
    count_table[1,2]<-length(x)-count_table[1,1]
    ####not in the gene list but in the pathway
    count_table[2,1]<-sum(background_non_gene_list%in% path)
    ####not in the gene list and not in the pathway
    count_table[2,2]<-length(background_non_gene_list)-count_table[2,1]
    matched_gene<-x[x %in% path]
    match_num<-length(matched_gene)
    overlap_info<-array(0,dim=4)
    names(overlap_info)<-c("DE in Geneset","DE not in Genese","NonDE in Geneset","NonDE out of Geneset")
    overlap_info[1]=count_table[1,1]
    overlap_info[2]=count_table[1,2]
    overlap_info[3]=count_table[2,1]
    overlap_info[4]=count_table[2,2]
    if(length(count_table)==4){
      res <- fisher.test(count_table, alternative="greater")$p}
    return(list(p_value=res,match_gene=matched_gene,match_num=match_num,
                fisher_table=overlap_info))
  }
  p_val<-rep(0,length(pathway))
  
  match_gene_list<-list(length(pathway))
  match_gene <- array(0,dim=length(pathway))
  num1<-array(0,dim=length(pathway))
  num2<-matrix(0,nrow=length(pathway),ncol=4)
  colnames(num2)<-c("DE in Geneset","DE not in Genese","NonDE in Geneset","NonDE out of Geneset")
  for(i in 1:length(pathway)){
    result<-get.fisher(pathway[[i]])
    p_val[i]<-result$p_value
    match_gene_list[[i]]<-result$match_gene
    match_gene[i]<-paste(match_gene_list[[i]],collapse="/")
    num1[i]<-result$match_num
    num2[i,]<-result$fisher_table
  }
  names(p_val) <- names(pathway)
  q_val <- p.adjust(p_val, "BH")
  
  summary<-data.frame(pvalue=p_val,
                      qvalue=q_val,
                      DE_in_Set=num2[,1],
                      DE_not_in_Set=num2[,2],
                      NonDE_in_Set=num2[,3],
                      NonDE_not_in_Set=num2[,4])
  
  a<-format(summary,digits=3)
  
  return(a)
}

# Update 10/17/24
# gsa.fisher.circadian <- function(rhythmic_genes, arrhythmic_genes, background, pathway) {
#   #### rhythmic_genes and arrhythmic_genes: binary (1 = rhythmic, 0 = arrhythmic)
#   #### background: a list of all genes in the background
#   #### pathway: a list of pathway genes
#   
#   count_table <- matrix(0, 2, 2)  # Fisher's Exact Test contingency table
#   
#   # Pathway genes adjusted to intersect with the background
#   pathway <- lapply(pathway, function(path) intersect(toupper(background), toupper(path)))
#   
#   # Function to perform Fisher's Exact Test
#   get.fisher <- function(path) {
#     res <- NA
#     # Rhythmic genes in the pathway
#     count_table[1, 1] <- sum(rhythmic_genes %in% path)
#     
#     # Rhythmic genes but not in the pathway
#     count_table[1, 2] <- sum(rhythmic_genes %in% background) - count_table[1, 1]
#     
#     # Arrhythmic genes in the pathway
#     count_table[2, 1] <- sum(arrhythmic_genes %in% path)
#     
#     # Arrhythmic genes not in the pathway
#     count_table[2, 2] <- sum(arrhythmic_genes %in% background) - count_table[2, 1]
#     
#     matched_genes <- rhythmic_genes[rhythmic_genes %in% path]
#     match_num <- length(matched_genes)
#     
#     overlap_info <- c(
#       "Rhythmic in Pathway" = count_table[1, 1],
#       "Rhythmic not in Pathway" = count_table[1, 2],
#       "Arrhythmic in Pathway" = count_table[2, 1],
#       "Arrhythmic not in Pathway" = count_table[2, 2]
#     )
#     if (length(count_table) == 4) {
#       res <- fisher.test(count_table, alternative = "greater")$p
#     }
#     
#     return(list(p_value = res, match_genes = matched_genes, match_num = match_num, fisher_table = overlap_info))
#   }
#   
#   # Prepare to collect results for each pathway
#   p_vals <- numeric(length(pathway))
#   fisher_results <- list()
#   
#   for (i in seq_along(pathway)) {
#     result <- get.fisher(pathway[[i]])
#     p_vals[i] <- result$p_value
#     fisher_results[[i]] <- result
#   }
#   
#   # Adjust p-values using Benjamini-Hochberg correction
#   q_vals <- p.adjust(p_vals, method = "BH")
#   
#   # Return summary table with the results
#   summary <- data.frame(
#     pvalue = p_vals,
#     qvalue = q_vals,
#     Rhythmic_in_Set = sapply(fisher_results, function(res) res$fisher_table["Rhythmic in Pathway"]),
#     Rhythmic_not_in_Set = sapply(fisher_results, function(res) res$fisher_table["Rhythmic not in Pathway"]),
#     Arrhythmic_in_Set = sapply(fisher_results, function(res) res$fisher_table["Arrhythmic in Pathway"]),
#     Arrhythmic_not_in_Set = sapply(fisher_results, function(res) res$fisher_table["Arrhythmic not in Pathway"])
#   )
#   
#   return(summary)
# }
# 
# 
# fisher <- function(x){
#   n <- length(x)
#   y <- -2*log(x)
#   Tf <- sum(y)
#   return(1-pchisq(Tf,2*n))
# }

# Update 03/21


# Fisher’s method for combining p-values
fisher <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  n <- length(x)
  y <- -2 * log(x)
  Tf <- sum(y)
  return(1 - pchisq(Tf, 2 * n))
}

# Fisher’s Exact Test for circadian genes
gsa.fisher.circadian <- function(rhythmic_genes, arrhythmic_genes, background, pathway) {
  count_table <- matrix(0, 2, 2)
  pathway <- lapply(pathway, function(path) intersect(toupper(background), toupper(path)))
  
  get.fisher <- function(path) {
    res <- NA
    count_table[1, 1] <- sum(rhythmic_genes %in% path)
    count_table[1, 2] <- sum(rhythmic_genes %in% background) - count_table[1, 1]
    count_table[2, 1] <- sum(arrhythmic_genes %in% path)
    count_table[2, 2] <- sum(arrhythmic_genes %in% background) - count_table[2, 1]
    
    overlap_info <- c(
      "Rhythmic in Pathway" = count_table[1, 1],
      "Rhythmic not in Pathway" = count_table[1, 2],
      "Arrhythmic in Pathway" = count_table[2, 1],
      "Arrhythmic not in Pathway" = count_table[2, 2]
    )
    
    if (all(count_table >= 0) && length(count_table) == 4) {
      res <- fisher.test(count_table, alternative = "greater")$p
    }
    return(list(p_value = res, match_genes = rhythmic_genes[rhythmic_genes %in% path], 
                match_num = length(rhythmic_genes[rhythmic_genes %in% path]), 
                fisher_table = overlap_info))
  }
  
  p_vals <- numeric(length(pathway))
  fisher_results <- list()
  for (i in seq_along(pathway)) {
    result <- get.fisher(pathway[[i]])
    p_vals[i] <- result$p_value
    fisher_results[[i]] <- result
  }
  
  q_vals <- p.adjust(p_vals, method = "BH")
  summary <- data.frame(
    pvalue = p_vals,
    qvalue = q_vals,
    Rhythmic_in_Set = sapply(fisher_results, function(res) res$fisher_table["Rhythmic in Pathway"]),
    Rhythmic_not_in_Set = sapply(fisher_results, function(res) res$fisher_table["Rhythmic not in Pathway"]),
    Arrhythmic_in_Set = sapply(fisher_results, function(res) res$fisher_table["Arrhythmic in Pathway"]),
    Arrhythmic_not_in_Set = sapply(fisher_results, function(res) res$fisher_table["Arrhythmic not in Pathway"])
  )
  return(summary)
}

##########################
###### SA functions ######
##########################

## energy function

E_tot <- function(delta.mat,a,delta.est){
  ## vector "a" of length n
  ## vector "delta_est" of length K+1: start from theta_0, then ordered from k=1 to K
  n <- nrow(delta.mat)
  K <- length(unique(a))
  #theta_0 <- delta.est[1]
  theta_0 <- 0
  E <- sum(sapply(1:n, function(x) {
    sum(sapply(1:n, function(y){
      if(a[x]==a[y]){
        (delta.mat[x,y] - delta.est[as.character(a[x])])^2
      } else{
        (delta.mat[x,y] - theta_0)^2
      }
    },simplify=T))
  }, simplify=T))
  
  return(E)
}

## Estimate of delta.est (the means)

Est_mean <- function(delta.mat,a){
  # the first element is always the off-diagonal parts
  n <- nrow(delta.mat)
  K <- length(unique(a))
  total <- rep(0,K+1)
  size <- rep(0,K+1)
  deltamean <- rep(0,K+1)
  names(deltamean) <- c(0,sort(unique(a)))
  for(k in 1:K){
    a_k <- sort(unique(a))[k]
    total[1+k] <- sum(delta.mat[a==a_k,a==a_k])
    size[1+k] <- sum(a==a_k)^2
    deltamean[1+k] <- total[1+k]/size[1+k]
  }
  deltamean[1] <- (sum(delta.mat) - sum(total[-1]))/(n*n - sum(size[-1]))
  return(deltamean)
}

## Trial = split or relocate

Split <- function(a) {
  n <- length(a)
  ua <- unique(a)
  if(length(ua)==n) {
    return(a)
  } else{
    #ua.pick <- sample(x=ua,size=1)
    #a[names(sample(x=which(a==ua.pick),size=1))] <- max(ua)+1
    a.pick <- sample(x=a,size=1)
    pick.ind <- sample(x=which(a==a.pick),size=1)
    a[pick.ind] <- max(a)+1
    return(a)
  }
}

Relocate <- function(a){
  n <- length(a)
  ua <- unique(a)
  if(length(ua)==1) {
    return(a)
  } else{
    pick.ind <- sample(x=1:n,size=1)
    #a.pick <- a[pick.ind]
    #a[pick.ind] <- sample(x=a[-which(a==a.pick)],size=1)
    a[pick.ind] <- sample(x=a[-pick.ind],size=1)
    return(a)
  }
}


##########################
###### Scatterness #######
##########################


scatter = function(dat,cluster.assign,sil_cut=0.1){
  sd_check = apply(dat, 1, sd)
  
  if(ncol(dat) == 1|!all(sd_check != 0) ){
    dist.dat = as.matrix(dist(dat))
  }else{
    dist.dat = 1 - cor(t(dat))
  }
  sil = silhouette(cluster.assign, dist=dist.dat,diss=T)
  
  new.dist = dist.dat
  cluster.assign2 = cluster.assign
  
  while(min(sil[,3]) < sil_cut){
    cluster.assign2<-cluster.assign2[-which(sil[,3] == min(sil[,3]))]
    for (i in unique(cluster.assign2)){
      if(length(which(cluster.assign2==i))==1){
        cluster.assign2<-cluster.assign2[-which(cluster.assign2==i)]
      }
    }
    temp<-cluster.assign2
    for(d in 1:length(cluster.assign2)){
      cluster.assign2[d]<-rank(unique(temp))[which(unique(temp)==temp[d])]
    }#rename cluster index, so it is integer from 1 to k
    
    new.dist<-new.dist[rownames(new.dist)%in%names(cluster.assign2),
                       colnames(new.dist)%in%names(cluster.assign2)]
    sil <- silhouette(cluster.assign2, dist=new.dist, diss=T)#recalculate silhoutte
  }
  
  scatter.index = which(!names(cluster.assign)%in%names(cluster.assign2))
  return(scatter.index)
}

##########################
###### Text mining #######
##########################
TextMine <- function(hashtb, pathways, pathway, result, scatter.index=NULL,permutation="all"){
  cat("Performing Text Mining Analysis...\n")
  hashtbAll = hashtb
  k <- length(unique(result))
  if(!is.null(scatter.index)){
    cluster = result
    cluster[scatter.index] = k+1
    hashtb = hashtb[hashtb [,3]%in%which(pathways %in% pathway),]
    tmk = list()
    nperm = 1000
    if (nrow(hashtb) == 0){
      for (i in 1:(k-1)){
        tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)
      }
    }
    else{
      for (i in 1:k){
        e = cluster[cluster == i]
        e = which(pathways %in% names(e))
        hashcl = hashtb[hashtb [,3]%in%e,]
        hashcl = hashcl[duplicated(hashcl[,2]) | duplicated(hashcl[,2], fromLast=TRUE),]
        if (nrow(hashcl) != 0){
          hashf = hashcl
          hashf[,2] = 1
          hashf = aggregate(hashf[,c("row","value")],by = hashf["phrase"],FUN=sum)
          
          rownames(hashf) = hashf[,"phrase"]
          hashf = hashf[,-1]
          colnames(hashf) = c("count","sum")
          hashap = hashcl
          hashap[,c(2,3,4)] = 0
          mperm = matrix(nrow = nrow(hashf),ncol = nperm)
          for (j in 1:nperm){
            if (permutation=="all"){
              subtb = hashtbAll[hashtbAll[,3]%in%sample(1:length(pathways),length(e)),]
            }else if(permutation=="enriched"){
              subtb = hashtb[hashtb[,3]%in%sample(unique(hashtb[,3]),length(e)),]
            }else{
              stop("Permutation should be 'all' or 'enriched' ")
            }
            subtb = rbind(subtb,hashap)
            subtb = subtb[subtb$phrase %in% hashap$phrase,]
            subtb[,2] = 1
            subtb = aggregate(subtb[,c("row","value")],by = subtb["phrase"],FUN=sum)
            rownames(subtb) = subtb[,"phrase"]
            subtb = subtb[,-1]
            colnames(subtb) = c("count","sum")
            mperm[,j] = subtb[,2]
          }
          hashf[,"p-value"] = apply(cbind(hashf[,2],mperm),1,
                                    function(x)((nperm + 2)-rank(x)[1])/(nperm + 1))
          hashf[,"q-vlaue"] = p.adjust(hashf[,"p-value"],method = "BH")
          tmk[[i]] = hashf[order(hashf[,3],-hashf[,2]),]
        }
        else {tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)}
      }
      tmk[[k+1]] = matrix(NA,nrow = 1,ncol = 4)
    }
  }else{
    hashtb = hashtb[hashtb [,2]%in%which(pathways %in% pathway),]
    tmk = list()
    nperm = 1000
    if (nrow(hashtb) == 0){
      for (i in 1:(k-1)){
        tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)
      }
    }
    else{
      for (i in 1:k){
        e = cluster[cluster == i]
        e = which(pathways %in% names(e))
        hashcl = hashtb[hashtb [,3]%in%e,]
        hashcl = hashcl[duplicated(hashcl[,2]) | duplicated(hashcl[,2], fromLast=TRUE),]
        if (nrow(hashcl) != 0){
          hashf = hashcl
          hashf[,2] = 1
          hashf = aggregate(hashf[,c("row","value")],by = hashf["phrase"],FUN=sum)
          
          rownames(hashf) = hashf[,"phrase"]
          hashf = hashf[,-1]
          colnames(hashf) = c("count","sum")
          hashap = hashcl
          hashap[,c(2,3,4)] = 0
          mperm = matrix(nrow = nrow(hashf),ncol = nperm)
          for (j in 1:nperm){
            if (permutation=="all"){
              subtb = hashtbAll[hashtbAll[,3]%in%sample(1:length(pathways),length(e)),]
            }else if(permutation=="enriched"){
              subtb = hashtb[hashtb[,3]%in%sample(unique(hashtb[,3]),length(e)),]
            }else{
              stop("Permutation should be 'all' or 'enriched' ")
            }
            subtb = rbind(subtb,hashap)
            subtb = subtb[subtb$phrase %in% hashap$phrase,]
            subtb[,2] = 1
            subtb = aggregate(subtb[,c("row","value")],by = subtb["phrase"],FUN=sum)
            rownames(subtb) = subtb[,"phrase"]
            subtb = subtb[,-1]
            colnames(subtb) = c("count","sum")
            
            mperm[,j] = subtb[,2]
          }
          hashf[,"p-value"] = apply(cbind(hashf[,2],mperm),1,
                                    function(x)((nperm + 2)-rank(x)[1])/(nperm + 1))
          hashf[,"q-vlaue"] = p.adjust(hashf[,"p-value"],method = "BH")
          tmk[[i]] = hashf[order(hashf[,3],-hashf[,2]),]
        }
        else {tmk[[i]] = matrix(NA,nrow = 1,ncol = 4)}
      }
    }
  }
  return(tmk)
} # End of Text Mining


writeTextOut <- function(tm_filtered,k,pathway.summary,scatter.index=NULL) {
  if(is.null(dim(tm_filtered[[1]]))==TRUE|dim(tm_filtered[[1]])[1] == 0){
    print(paste("No phrase pass q-value threshold in cluster 1"))
    cat("Cluster 1\n", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(pathway.summary[[1]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T, append = T, row.names=F,col.names=F)
  } else {
    cat("Cluster 1\n", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    cat("Key words,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(as.character(rownames(tm_filtered[[1]])[1:15])), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    cat("q_value,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(tm_filtered[[1]][1:15,4]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    cat("count,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
    write.table(t(tm_filtered[[1]][1:15,1]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
    write.table(pathway.summary[[1]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T, append = T, row.names=F,col.names=F)
    
  }
  
  for (i in 2:k){
    if(is.null(dim(tm_filtered[[i]]))==TRUE|dim(tm_filtered[[i]])[1] == 0){
      print(paste("No phrase pass q-value threshold in cluster ",k,sep = ""))
      cat(paste("\nCluster ", i, "\n", sep = ""), file = paste("Clustering_Summary_K",k,".csv",sep=""), append = T)
      write.table(pathway.summary[[i]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T, append = T, row.names=F,col.names=F)
    } else {
      cat(paste("\nCluster ", i, "\n", sep = ""), file = paste("Clustering_Summary_K",k,".csv",sep=""), append = T)
      cat("Key words,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
      write.table(t(as.character(rownames(tm_filtered[[i]])[1:15])), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("q_value,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
      write.table(t(tm_filtered[[i]][1:15,4]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      cat("count,", file = paste("Clustering_Summary_K",k,".csv",sep=""),append=T)
      write.table(t(tm_filtered[[i]][1:15,1]), paste("Clustering_Summary_K",k,".csv",sep=""), sep=',',quote=F, append = T, row.names=F,col.names=F,na="")
      write.table(pathway.summary[[i]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T, append = T, row.names=F,col.names=F)
    }
    
  }
  if(!is.null(scatter.index)){
    cat(paste("\nSingleton Term", "\n", sep = ""), file = paste("Clustering_Summary_K",k,".csv",sep=""), append = T)
    write.table(pathway.summary[[k+1]], paste("Clustering_Summary_K",k,".csv",sep=""), sep=",",quote=T,
                append = T, row.names=F,col.names=F)
  }
}

textMine <- function(hashtb,pathways,cluster.assign,scatter.index=NULL,thres=0.05,permutation="all"){
  tmk <- TextMine(hashtb=hashtb, pathways= pathways,
                  pathway=names(cluster.assign), result=cluster.assign, scatter.index,permutation=permutation)
  C <- length(unique(cluster.assign))
  tm_filtered <- list()
  for (i in 1:C){
    tm_filtered[[i]] <- tmk[[i]][which((as.numeric(tmk[[i]][,4]) < thres)), ]
  }
  if(!is.null(scatter.index)){
    tm_filtered[[C+1]] = tmk[[C+1]]
    cluster = cluster.assign
    cluster[scatter.index] = C+1
    pathway.summary <- lapply(1:(C+1), function(x) names(which(cluster==x)))
    writeTextOut(tm_filtered,C,pathway.summary,scatter.index=scatter.index)
  }else{
    pathway.summary <- lapply(1:C, function(x) names(which(cluster.assign==x)))
    writeTextOut(tm_filtered,C,pathway.summary)
  }
  return(tm_filtered)
}

##########################
###### ACS/ADS_DE plot ###
##########################
# ARS_to_size <- function(ARSp,factor=2){
#   if(ARSp > 0.05) {
#     return(1)
#   } else {
#     return(-log10(ARSp)*factor)
#   }
# }
# 
# ACS_ADS_DE <- function(ds1,ds2,DEevid1,DEevid2,ACSp,ADSp,cluster=NULL,
#                        highlight.pathways=NULL,lb=0,ub=1,size.scale=4){
#   
#   P <- length(ACSp)
#   ACS_size=sapply(ACSp,ARS_to_size)
#   ADS_size=sapply(ADSp,ARS_to_size)
#   
#   
#   if(!is.null(cluster)){
#     RB = rainbow(length(unique(cluster)))
#     color = c()
#     for (i in 1:length(cluster)) {
#       if(cluster[i] != "scatter"){
#         color[i] = RB[as.numeric(cluster[i])]
#       }else{
#         color[i] = "grey50"
#       }
#     }
#   }else if(!is.null(highlight.pathways)){
#     color = rep("grey50",P)
#     color[highlight.pathways] = "red"
#   }else{
#     color = rep("black",P)
#   }
#   
#   
#   if(!is.null(highlight.pathways)){
#     index = 1:P
#     index[-highlight.pathways] = ""
#   }else{
#     index = rep("",P)
#   }
#   data_ACS <- data.frame(ds1_score=DEevid1,ds2_score=DEevid2,
#                          ACS_size=ACS_size,index=index,
#                          color_pos=color)
#   data_ADS <- data.frame(ds1_score=DEevid1,ds2_score=DEevid2,
#                          ADS_size=ADS_size,index=index,
#                          color_neg=color)
#   
#   p_pos <-ggplot(data_ACS, aes(x=ds2_score, y=ds1_score,label=index)) +
#     geom_point(size = data_ACS$ACS_size*size.scale, shape=16, color = data_ACS$color_pos)+
#     geom_text(size=8*size.scale,parse=TRUE,color="black",hjust = -0.05,vjust=-0.05) +
#     theme_bw() +
#     coord_fixed(ylim=c(lb,ub),xlim=c(lb,ub)) +
#     labs(x="",y="") +
#     scale_x_continuous(name="",breaks=seq(0,1,by=0.5),limits=c(0,1)) +
#     scale_y_continuous(name="",breaks=seq(0,1,by=0.5),limits=c(0,1)) +
#     theme(#legend.title = element_blank(),
#       axis.line = element_line(colour = "black"),
#       #axis.line=element_blank(),
#       axis.text.x = element_text(size = 40,face = "bold"),
#       axis.text.y = element_text(size = 40,face = "bold"),
#       panel.border = element_blank(),
#       panel.grid.major = element_line(linetype = 'solid',#size = 2,
#                                       colour = "white"),
#       panel.grid.minor = element_line(linetype = 'solid',
#                                       colour = "white"),
#       panel.background = element_rect(fill = "#FFF1E1")) +
#     annotate("text", x = (ub-0.15), y = lb, fontface=2,
#              label = paste(ds2,sep=""),
#              size=8*size.scale,colour="blue",hjust=0.6,vjust=0.1) +
#     annotate("text", x = lb, y = (ub-0.15), fontface=2,
#              label=paste(ds1,sep=""),
#              size=8*size.scale,colour="blue",vjust=0,hjust=0.2)
#   
#   p_neg <-ggplot(data_ADS, aes(x=ds1_score, y=ds2_score,label=index)) + #label=index
#     geom_point(size = data_ADS$ADS_size*size.scale, shape=16, color = data_ADS$color_neg)+
#     geom_text(size=8*size.scale,parse=TRUE,color="black",hjust = -0.05,vjust=-0.05) +
#     theme_bw() +
#     coord_fixed(ylim=c(lb,ub),xlim=c(lb,ub)) +
#     labs(x="",y="") +
#     scale_x_continuous(name="",breaks=seq(0,1,by=0.5),limits=c(0,1)) +
#     scale_y_continuous(name="",breaks=seq(0,1,by=0.5),limits=c(0,1)) +
#     theme(#legend.title = element_blank(),
#       axis.line = element_line(colour = "black"),
#       #axis.line=element_blank(),
#       axis.text.x = element_text(size = 40,face = "bold"),
#       axis.text.y = element_text(size = 40,face = "bold"),
#       panel.border = element_blank(),
#       panel.grid.major = element_line(linetype = 'solid',#size = 2,
#                                       colour = "white"),
#       panel.grid.minor = element_line(linetype = 'solid',
#                                       colour = "white"),
#       panel.background = element_rect(fill = "#EBF5FF")) +
#     annotate("text", x = (ub-0.15), y = lb, fontface=2,
#              label = paste(ds1,sep=""),
#              size=8*size.scale,colour="blue",hjust=0.6,vjust=0.1) +
#     annotate("text", x = lb, y = (ub-0.15), fontface=2,
#              label=paste(ds2,sep=""),
#              size=8*size.scale,colour="blue",vjust=0,hjust=0.2)
#   
#   #ggsave(filename=paste(ds2,"_",ds1,"_ACS_figure",".pdf",sep=""),p_pos,
#   #width = 10, height = 10)
#   
#   #ggsave(filename=paste(ds1,"_",ds2,"_ADS_figure",".pdf",sep=""),p_neg,
#   #width = 10, height = 10)
#   
#   plist <- list(p_pos,p_neg)
#   return(plist)
# }

##########################
######    parseXML     ###
##########################
parseRelation <- function(pathwayID, keggSpecies="hsa", binary = T, sep = "-") {
  # download xml file
  pathview::download.kegg(pathway.id = pathwayID, keggSpecies, kegg.dir = ".", file.type="xml")
  # generate relation matrix
  KEGG.pathID2name = lapply(KEGGREST::keggList("pathway",keggSpecies),function(x) strsplit(x," - ")[[1]][-length(strsplit(x," - ")[[1]])])
  names(KEGG.pathID2name) = gsub(paste0("path:",keggSpecies),"",names(KEGG.pathID2name))
  
  pathName = unlist(KEGG.pathID2name[pathwayID])
  
  xmlFile = paste0(getwd(), "/",keggSpecies, pathwayID,".xml")
  pathway = KEGGgraph::parseKGML(xmlFile)
  pathway = KEGGgraph::splitKEGGgroup(pathway)
  
  entries = KEGGgraph::nodes(pathway)
  types = sapply(entries, KEGGgraph::getType)
  relations = unique(KEGGgraph::edges(pathway)) ## to avoid duplicated edges
  relationNum = length(relations)
  entryNames = as.list(sapply(entries, KEGGgraph::getName))
  if(any(types == "group") || any(types=="map")){
    entryNames = entryNames[!(types %in% c("group","map"))]
  }
  entryIds = names(entryNames)
  entryNames = lapply(1:length(entryNames), function(i) paste(entryNames[[i]],collapse=sep))
  names(entryNames) = entryIds
  
  entryNames.unique = unique(entryNames)
  entryNum = length(entryNames.unique)
  
  relation.mat = matrix(0, entryNum, entryNum)
  rownames(relation.mat) = colnames(relation.mat) = entryNames.unique
  
  ## if no relation edge, just return
  if(relationNum == 0){
    print(paste0("There is no topological connected gene nodes in ", pathName))
    return(relation.mat)
  }
  
  entry1 = KEGGgraph::getEntryID(relations)[,1]
  entry2 = KEGGgraph::getEntryID(relations)[,2]
  for(i in 1:length(relations)){
    if(entry1[i] %in% names(entryNames) && entry2[i] %in% names(entryNames)){
      relation.mat[entryNames[[entry1[i]]],entryNames[[entry2[i]]]]=1
    }
    else{
      print(paste("connections not included:",entry1[i], entry2[i], sep=" "))
    }
  }
  
  file.remove(xmlFile)
  return(relation.mat)
}
##########################
###### KEGG module SA ####
##########################
SA_module_M = function(sp.mat, xmlG, M, nodes, B = 1000,
                       G.ini.list=NULL, reps_eachM = 100,topG_from_previous=10,
                       Tm0=10,mu=0.95,epsilon=1e-5,
                       N=1000,run=10000,seed=12345,sub.num=1){
  #Null distribution for M
  set.seed(seed)
  null.sp.dist = rep(NA,B)
  for(b in 1:B){
    permute.set <- sample(xmlG,M)
    permute.mat <- sp.mat[match(permute.set,row.names(sp.mat)),
                          match(permute.set,colnames(sp.mat))]
    null.sp.dist[b] <- mean(c(permute.mat[lower.tri(permute.mat)]))
  }
  null.sp.mean = mean(null.sp.dist)
  null.sp.median = median(null.sp.dist)
  if(is.null(G.ini.list)){
    p.mean.ls = c()
    G.module.ls = list()
    SP.ls = c()
    for (l in 1:reps_eachM) {
      ##Initialize
      if(length(nodes) == M){
        G.module = nodes
      }else{
        G.module = sample(nodes, M)
      }
      SPc = avgSP(G.module, sp.mat)
      r = 0
      count = 0
      Tm = Tm0
      while((length(nodes)>M) & (r < run) & (count < N) & (Tm >= epsilon)) {
        #pi = exp(-GPc/Tm) ## Boltzmann dist #may need a different Tm or -logP to be comparable?
        #print(SPc)
        ##New trial
        r = r+1
        a.node = sample(setdiff(nodes,G.module),sub.num)
        G.module_new = G.module
        G.module_new[sample(M,sub.num)] = a.node
        
        SPn = avgSP(G.module_new, sp.mat)
        
        if(SPn < SPc | SPc == Inf) {
          ##accept
          SPc = SPn;
          G.module = G.module_new;
        }else{
          count = count + 1;
          p = exp((SPc-SPn)/Tm)
          #print(p)
          r = min(1,p); ## acceptance prob.
          u <- runif(1);
          if(u>r) {
            ##not accept
            Tm <- Tm*mu
          } else {
            SPc = SPn;
            G.module = G.module_new;
          }
        }
      }
      G.module.ls[[l]] = G.module
      p.mean.ls[l] = (sum(null.sp.dist<= SPc) + 1)/(B+1)
      SP.ls[l] = SPc
    }
  }else{
    each.times = round(reps_eachM/topG_from_previous)
    case.index = expand.grid(1:length(G.ini.list),1:each.times)
    p.mean.ls = c()
    G.module.ls = list()
    SP.ls = c()
    for (l in 1:nrow(case.index)) {
      G.ini = G.ini.list[[case.index[l,1]]]
      ##Initialize
      if(length(nodes) == M){
        G.module = nodes
      }else{
        x = M-length(G.ini)
        G.module = c(G.ini,sample(setdiff(nodes,G.ini),x))
      }
      SPc = avgSP(G.module, sp.mat)
      r = 0
      count = 0
      Tm = Tm0
      while((length(nodes)>M) & (r < run) & (count < N) & (Tm >= epsilon)) {
        #pi = exp(-GPc/Tm) ## Boltzmann dist #may need a different Tm or -logP to be comparable?
        #print(SPc)
        ##New trial
        r = r+1
        a.node = sample(setdiff(nodes,G.module),sub.num)
        G.module_new = G.module
        G.module_new[sample(M,sub.num)] = a.node
        
        SPn = avgSP(G.module_new, sp.mat)
        
        if(SPn < SPc | SPc == Inf) {
          ##accept
          SPc = SPn;
          G.module = G.module_new;
        }else{
          count = count + 1;
          p = exp((SPc-SPn)/Tm)
          #print(p)
          r = min(1,p); ## acceptance prob.
          u <- runif(1);
          if(u>r) {
            ##not accept
            Tm <- Tm*mu
          } else {
            SPc = SPn;
            G.module = G.module_new;
          }
        }
      }
      G.module.ls[[l]] = G.module
      p.mean.ls[l] = (sum(null.sp.dist<= SPc) + 1)/(B+1)
      SP.ls[l] = SPc
    }
    
  }
  p.sd.ls = sqrt(p.mean.ls*(1-p.mean.ls)/B)
  r.p = rank(p.mean.ls,ties.method = "random")
  index = match(1:topG_from_previous,r.p)
  best.index = which(r.p == 1)
  
  minG = G.module.ls[[best.index]]
  sp = SP.ls[best.index]
  p.mean = p.mean.ls[best.index]
  p.sd = p.sd.ls[best.index]
  
  top.G = G.module.ls[index]
  top.pmean = p.mean.ls[index]
  top.psd = p.sd.ls[index]
  top.sp = SP.ls[index]
  
  return(list(minG = minG,sp = sp,p.mean = p.mean,p.sd = p.sd,
              top.G = top.G,top.sp = top.sp,top.pmean = top.pmean,top.psd = top.psd,
              null.sp.mean = null.sp.mean,null.sp.median = null.sp.median))
}

avgSP = function(G.module, sp.mat){
  m = length(G.module)
  set.mat = sp.mat[match(G.module,row.names(sp.mat)),
                   match(G.module,colnames(sp.mat))]
  G.sp = mean(c(set.mat[lower.tri(set.mat)]))
  return(G.sp)
}

#' Annotate MCMC output with gene symbols and rhythmicity index
#'
#' @title Map Ensembl IDs to gene symbols and set rhythmicity attributes
#'
#' @description
#' Takes the raw list output from \code{CB_MCMC_single_rj_slice} and
#' (optionally) converts Ensembl row names to HGNC gene symbols via
#' biomaRt.  If row names are already symbols they are used directly.
#' Computes a Bayes Factor for each gene using \code{summarize_bay} and
#' sets \code{attr(rho, "symbols")} and \code{attr(rho, "RHYindex")} so
#' that downstream functions can access gene identifiers and rhythmicity
#' classifications.  Duplicate symbols are resolved by keeping the row
#' with the highest Bayes Factor.
#'
#' @param input_df Named list; raw MCMC output (e.g. from
#'   \code{CB_MCMC_single_rj_slice}).  Must contain a \code{rho} matrix
#'   with gene identifiers as row names.
#' @param BF Numeric; Bayes Factor threshold above which a gene is called
#'   rhythmic.
#' @param p_rhythmic Numeric in (0, 1); prior probability of rhythmicity
#'   used to compute the Bayes Factor (default 0.5).
#' @param ensemble A biomaRt Mart object; required when row names are
#'   Ensembl IDs (default \code{NULL}).
#'
#' @return The input list filtered to unique symbols, with all matrix
#'   elements row-subsetted and \code{attr(*, "symbols")},
#'   \code{attr(*, "RHYindex")}, and (if Ensembl IDs were present)
#'   \code{attr(*, "ensembl_gene_ids")} set on each matrix.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' res  <- CB_MCMC_single_rj_slice(dat, init)
#' res2 <- match_symbols(res, BF = 3, p_rhythmic = 0.2)
#' }
match_symbols <- function(input_df, BF, p_rhythmic = 0.5, ensemble = NULL) {
  if (!requireNamespace("biomaRt", quietly = TRUE))
    stop("Package 'biomaRt' is required. Install with: BiocManager::install('biomaRt')")
  
  # Extract the rho matrix from the input dataframe.
  summary_data <- input_df$rho
  
  # Get gene IDs from row names
  gene_ids <- row.names(summary_data)
  
  # Check if row names are already gene symbols (not Ensembl IDs)
  # More robust check: majority should be gene symbols (not ENS IDs)
  ens_pattern <- grepl("^ENS[MG]", gene_ids)  # ENSEMBL or ENSMUSG patterns
  pct_ensembl <- mean(ens_pattern)
  avg_length <- mean(nchar(gene_ids))
  pct_letters <- mean(grepl("^[A-Za-z]", gene_ids))
  
  print(paste("Percentage of Ensembl IDs:", round(pct_ensembl * 100, 1), "%"))
  print(paste("Average gene ID length:", round(avg_length, 1)))
  print(paste("Percentage starting with letters:", round(pct_letters * 100, 1), "%"))
  
  # If less than 10% are Ensembl IDs, assume they're already symbols
  already_symbols <- pct_ensembl < 0.1 && 
    avg_length < 20 &&  # Gene symbols are typically shorter
    pct_letters > 0.9   # Most start with letters
  
  if (already_symbols) {
    print("Row names appear to be gene symbols already. Adding attributes directly...")
    symbols <- gene_ids
  } else {
    print("Row names appear to be Ensembl IDs. Converting to gene symbols...")
    
    if (is.null(ensemble)) {
      stop("Ensemble mart object required for Ensembl ID conversion")
    }
    
    ensembl = ensemble
    Sys.setenv(XDG_CACHE_HOME = "~/biocache")
    
    get_gene_symbols <- function(gene_ids, mart) {
      if (all(grepl("^ENSG", gene_ids))) {
        filt <- "ensembl_gene_id"
      } else {
        filt <- "hgnc_symbol"
      }
      
      getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol"),
        filters = filt,
        values  = gene_ids,
        mart    = mart
      )
    }
    
    gene_symbols <- get_gene_symbols(gene_ids, ensembl)
    
    # Replace missing or empty hgnc_symbol entries with the ensembl_gene_id.
    missing_idx <- which(gene_symbols$hgnc_symbol == "" | is.na(gene_symbols$hgnc_symbol))
    gene_symbols$hgnc_symbol[missing_idx] <- gene_symbols$ensembl_gene_id[missing_idx]
    
    # Map the original gene_ids to their corresponding symbols.
    symbols <- if (all(grepl("^ENS", gene_ids))) {
      gene_symbols$hgnc_symbol[match(gene_ids, gene_symbols$ensembl_gene_id)]
    } else {
      gene_ids
    }
  }
  
  # Set the "symbols" attribute for the rho matrix.
  attr(summary_data, "symbols") <- symbols
  
  # Compute summary statistics (e.g., Bayes Factor, rhythmicity)
  summary_df <- summarize_bay(input_df$rho, BF, p_rhythmic)
  
  # Set the rhythmicity index attribute ("RHYindex") for the rho matrix.
  attr(summary_data, "RHYindex") <- summary_df$Rhythmicity
  
  # Resolve duplicate gene symbols by keeping the row with the highest Bayes Factor.
  symbols_unique <- unique(symbols)
  data_filtered <- sapply(symbols_unique, function(sym) {
    indices <- which(symbols == sym)
    if (length(indices) > 1) {
      max_idx <- indices[which.max(summary_df$BayesF[indices])]
      return(max_idx)
    } else {
      return(indices)
    }
  })
  data_filtered <- unlist(data_filtered)
  
  # Filter every element of input_df (if matrix or data.frame) using the selected indices.
  input_df_filtered <- lapply(input_df, function(element) {
    if (is.matrix(element) || is.data.frame(element)) {
      return(element[data_filtered, , drop = FALSE])
    } else {
      return(element)
    }
  })
  
  final_symbols  <- symbols[data_filtered]
  final_RHYindex <- summary_df$Rhythmicity[data_filtered]
  final_ens_ids  <- if (!already_symbols) gene_ids[data_filtered] else NULL
  G_final        <- length(final_symbols)

  # Propagate symbols, RHYindex (and Ensembl IDs if applicable) to ALL G×K matrices.
  # This ensures rownames and attributes survive downstream subsetting operations.
  for (mat_name in names(input_df_filtered)) {
    el <- input_df_filtered[[mat_name]]
    if (is.matrix(el) && nrow(el) == G_final) {
      rownames(el)              <- final_symbols
      attr(el, "symbols")      <- final_symbols
      attr(el, "RHYindex")     <- final_RHYindex
      if (!is.null(final_ens_ids))
        attr(el, "ensembl_gene_ids") <- final_ens_ids
      input_df_filtered[[mat_name]] <- el
    }
  }

  # Legacy duplicate-check kept for user transparency
  tryCatch({
    if (length(unique(final_symbols)) < G_final)
      warning("Duplicate gene symbols found after deduplication — check BF threshold.")
  })
  
  return(input_df_filtered)
}
#' Align multiple MCMC outputs to a common set of homologous genes
#'
#' @title Match homologs across species MCMC outputs
#'
#' @description
#' Uses biomaRt to retrieve one-to-one orthologs between each non-reference
#' species and the reference species (default human), then intersects the
#' mapped gene sets across all datasets.  Each input list is row-filtered
#' and reordered so that all datasets share the same genes in the same order,
#' enabling direct cross-species comparisons.
#'
#' @param input_dfs A list of MCMC output lists (one per dataset), each
#'   with matrices whose row names are Ensembl gene IDs and whose
#'   \code{attr(rho, "ensembl_gene_ids")} attribute is set.
#' @param species_from Character vector of the same length as
#'   \code{input_dfs}; species label for each dataset (currently supports
#'   \code{"human"} and \code{"mouse"}).
#' @param ref Character; reference species for the common gene space
#'   (default \code{"human"}).
#'
#' @return A list of the same length as \code{input_dfs}, each element
#'   row-filtered to the intersection of homologous reference-space gene
#'   IDs and reordered consistently.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' aligned <- match_homologs(list(human_res, mouse_res),
#'                           species_from = c("human", "mouse"))
#' }
# Matching homologs
match_homologs <- function(input_dfs, species_from ,ref = "human") {
  if (length(input_dfs) != length(species_from)) {
    stop("The number of datasets must match the number of species provided.")
  }
  
  # Identify unique species and map to Ensembl names
  unique_species <- unique(species_from)
  species_map <- c("human" = "homo_sapiens", 
                   "mouse" = "mus_musculus")
  
  # Choose a reference species, e.g., human
  if(missing(ref) || ref == "human") {
    ref_species_name <- "human"
  } else if (ref != "human") {
    ref_species_name <- ref
    } else {
    ref_species_name <- species_from[1]
    }
  
  ref_species_ensembl <- species_map[ref_species_name]
  if (! ref_species_ensembl %in% species_map[unique_species]) {
    stop("Reference species (human) not found in your species_from.")
  }
  message("Using ", ref_species_ensembl, " as the reference.")
  
  # ============= 1) Build Mappings to Reference Space =============
  #
  # For each species S != reference, build a mapping:
  #   M_S: S_gene_ID -> human_gene_ID
  # For the reference species, the mapping is identity.
  mapping_list <- list()
  
  # (A) Gather all gene IDs from the reference species' datasets
  ref_dfs_indices <- which(species_map[species_from] == ref_species_ensembl)
  if (length(ref_dfs_indices) == 0) {
    stop("No datasets correspond to the reference species (human).")
  }
  
  ref_gene_ids <- unique(unlist(
    lapply(input_dfs[ref_dfs_indices], function(df) {
      attr(df$rho, "ensembl_gene_ids")
    })
  ))
  message("Reference species (human) gene IDs found: ", length(ref_gene_ids))
  
  # (B) Identity mapping for the reference species
  mapping_list[[ref_species_ensembl]] <- setNames(ref_gene_ids, ref_gene_ids)
  
  # (C) For each other species, retrieve homologs (using reference as source)
  other_species <- setdiff(unique_species, ref_species_name)
  for (sp in other_species) {
    sp_ensembl <- species_map[sp]
    message("Retrieving homologs from ", ref_species_ensembl, " to ", sp_ensembl)
    
    homologs <- getHomologs(
      unique(ref_gene_ids),
      species_from = ref_species_ensembl,
      species_to   = sp_ensembl
    )
    if (nrow(homologs) == 0) {
      stop("No homolog mappings found from ", ref_species_ensembl, " to ", sp_ensembl)
    }
    colnames(homologs) <- c("source_ensembl_id", "homolog_ensembl_id")
    
    # Create mapping: target species gene ID -> human gene ID
    mapping_list[[sp_ensembl]] <- setNames(
      homologs$source_ensembl_id,
      homologs$homolog_ensembl_id
    )
  }
  
  # ============= 2) Convert Each Dataset & Compute Intersection =============
  #
  # Map each dataset’s gene IDs into reference space and collect unique reference IDs.
  mapped_ids_list <- vector("list", length(input_dfs))
  
  for (i in seq_along(input_dfs)) {
    df <- input_dfs[[i]]
    sp_ensembl <- species_map[species_from[i]]
    
    original_ids <- attr(df$rho, "ensembl_gene_ids")
    mapped_ref_ids <- mapping_list[[sp_ensembl]][original_ids]
    mapped_ref_ids <- mapped_ref_ids[!is.na(mapped_ref_ids)]
    
    mapped_ids_list[[i]] <- unique(mapped_ref_ids)
  }
  
  # Compute the intersection of mapped (reference) IDs across all datasets.
  common_ref_ids <- Reduce(intersect, mapped_ids_list)
  message("Intersection of reference gene IDs across all datasets: ", length(common_ref_ids))
  
  if (length(common_ref_ids) == 0) {
    stop("No common homologs found across all datasets.")
  }
  
  # ============= 3) Final Filtering & Reordering of Each Dataset =============
  #
  # Filter each dataset to keep only rows corresponding to the common reference IDs,
  # and set the rownames to these stable Ensembl IDs.
  input_dfs_filtered <- lapply(seq_along(input_dfs), function(i) {
    df <- input_dfs[[i]]
    sp_ensembl <- species_map[species_from[i]]
    
    original_ids <- attr(df$rho, "ensembl_gene_ids")
    mapped_ref_ids <- mapping_list[[sp_ensembl]][original_ids]
    
    # Identify valid rows: those with a valid mapping and in the common set.
    valid_idx <- !is.na(mapped_ref_ids) & (mapped_ref_ids %in% common_ref_ids)
    if (sum(valid_idx) == 0) {
      message("Warning: No homologs matched for dataset ", i)
    }
    
    # Initial filtering based on valid_idx
    df_filtered <- lapply(df, function(element) {
      if (is.matrix(element) || is.data.frame(element)) {
        return(element[valid_idx, , drop = FALSE])
      } else {
        return(element)
      }
    })
    
    # Get the filtered mapped reference IDs and original symbols
    new_ids <- mapped_ref_ids[valid_idx]
    new_symbols <- attr(df$rho, "symbols")[valid_idx]
    
    # Reorder rows so that they exactly follow the common_ref_ids order.
    common_order <- common_ref_ids
    row_idx <- match(common_order, new_ids)
    row_idx <- row_idx[!is.na(row_idx)]
    
    df_filtered <- lapply(df_filtered, function(element) {
      if (is.matrix(element) || is.data.frame(element)) {
        element <- element[row_idx, , drop = FALSE]
      }
      element
    })
    
    final_symbols  <- new_symbols[row_idx]
    final_ens_ids  <- new_ids[row_idx]

    # Propagate gene symbols as rownames and set attributes on ALL G×K matrices.
    # Downstream functions (pathSelect, detect_rhy, phase_infer) all use rownames
    # for gene-level matching; symbols must be consistent across rho, phi, A, M, etc.
    G_final <- length(final_symbols)
    for (mat_name in names(df_filtered)) {
      el <- df_filtered[[mat_name]]
      if (is.matrix(el) && nrow(el) == G_final) {
        rownames(el) <- final_symbols
        attr(el, "symbols")          <- final_symbols
        attr(el, "ensembl_gene_ids") <- final_ens_ids
        df_filtered[[mat_name]] <- el
      }
    }

    df_filtered
  })
  
  return(input_dfs_filtered)
}

#' Compute pairwise rhythmicity concordance between two conditions
#'
#' @title Pairwise concordance (Jaccard and conditional) for two MCMC outputs
#'
#' @description
#' Compares binary posterior rhythmicity classifications (\code{rho} rows
#' thresholded at 0.5) between two conditions.  Returns the Jaccard index
#' (intersection / union), conditional proportions (how many rhythmic genes
#' in species 1 are also rhythmic in species 2, and vice versa), and raw
#' proportion counts.
#'
#' @param matrix1 Named list with element \code{rho} (G x K binary matrix);
#'   MCMC output for condition 1.
#' @param matrix2 Named list with element \code{rho} (G x K binary matrix);
#'   MCMC output for condition 2.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{conditional_species1_to_species2}{Proportion of species-1
#'       rhythmic genes also rhythmic in species 2.}
#'     \item{conditional_species2_to_species1}{Proportion of species-2
#'       rhythmic genes also rhythmic in species 1.}
#'     \item{jaccard_concordance}{Jaccard index (intersection / union).}
#'     \item{shared_rhythmic, species1_total_rhythmic,
#'       species2_total_rhythmic}{Raw proportions.}
#'     \item{interpretation}{Human-readable strings for the conditional
#'       proportions.}
#'   }
#'
concordance2 <- function(matrix1, matrix2) {
  rho1 <- matrix1$rho
  rho2 <- matrix2$rho
  
  # Standard 2x2 contingency table
  TP <- mean(rho1 == 1 & rho2 == 1)  # Both rhythmic
  FP <- mean(rho1 == 1 & rho2 == 0)  # Only matrix1 rhythmic
  FN <- mean(rho1 == 0 & rho2 == 1)  # Only matrix2 rhythmic
  TN <- mean(rho1 == 0 & rho2 == 0)  # Neither rhythmic
  
  # Published approach: Conditional probabilities
  
  # "Of genes rhythmic in species 1, what % are also rhythmic in species 2?"
  conditional_1to2 <- ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
  
  # "Of genes rhythmic in species 2, what % are also rhythmic in species 1?"  
  conditional_2to1 <- ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_)
  
  # Our bidirectional approach for comparison
  jaccard_index <- TP / (TP + FP + FN)  # Shared / Union
  
  # Total rhythmic genes in each species
  total_rhythmic_species1 <- TP + FP
  total_rhythmic_species2 <- TP + FN
  total_shared <- TP
  
  return(list(
    # Published approach
    conditional_species1_to_species2 = conditional_1to2,
    conditional_species2_to_species1 = conditional_2to1,
    
    # Our approach for comparison
    jaccard_concordance = jaccard_index,
    
    # Raw counts (as proportions)
    shared_rhythmic = total_shared,
    species1_total_rhythmic = total_rhythmic_species1,
    species2_total_rhythmic = total_rhythmic_species2,
    
    # Summary
    interpretation = list(
      conditional_1to2_meaning = paste0(round(conditional_1to2*100, 1), "% of species1 rhythmic genes are also rhythmic in species2"),
      conditional_2to1_meaning = paste0(round(conditional_2to1*100, 1), "% of species2 rhythmic genes are also rhythmic in species1")
    )
  ))
}

# concordance <- function(matrix1, matrix2) {
#   rho1 = matrix1$rho
#   rho2 = matrix2$rho
#   # Pairwise contingency table
#   TP <- mean(rho1 == 1 & rho2 == 1)
#   FP <- mean(rho1 == 1 & rho2 == 0)
#   FN <- mean(rho1 == 0 & rho2 == 1)
#   
#   # Jaccard index = intersection / union
#   congruence_index <- TP / (TP + FN + FP)
#   
#   return(congruence_index)
# }

# discordance <- function(matrix1, matrix2, eps = .Machine$double.eps) {
#   rho1 = matrix1$rho
#   rho2 = matrix2$rho
#   stopifnot(length(rho1) == length(rho2))
#   
#   TP <- mean(rho1 == 1 & rho2 == 1)
#   FP <- mean(rho1 == 1 & rho2 == 0)
#   FN <- mean(rho1 == 0 & rho2 == 1)
#   
#   # Total denominator (union of rhythmic genes)
#   total_denom <- TP + FN + FP
#   
#   # Gain and Loss as proportions of the union
#   gain <- if (total_denom > eps) FP / total_denom else NA_real_
#   loss <- if (total_denom > eps) FN / total_denom else NA_real_
#   
#   list(
#     gain = gain,
#     loss = loss
#   )
# }

#' Compute circular phase difference between two phase vectors
#'
#' @description
#' Computes phi2 - phi1 and wraps the result to \[-P/2, P/2\] (for hours)
#' or \[-pi, pi\] (for radians) to give the signed shortest-arc difference.
#'
#' @param phi1 Numeric vector; acrophase of condition 1 in hours (or
#'   radians if \code{units = "radians"}).
#' @param phi2 Numeric vector; acrophase of condition 2.
#' @param units Character; \code{"hours"} (default) or \code{"radians"}.
#'
#' @return Numeric vector of signed circular phase differences, wrapped to
#'   \[-12, 12\] hours or \[-pi, pi\] radians.
#'
# Helper function for phase differences
phase_difference <- function(phi1, phi2, units = "hours") {
  if (units == "hours") {
    diff <- phi2 - phi1
    # Wrap to [-12, +12] range
    diff <- ((diff + 12) %% 24) - 12
  } else if (units == "radians") {
    diff <- phi2 - phi1
    # Wrap to [-π, +π] range  
    diff <- ((diff + pi) %% (2*pi)) - pi
  }
  return(diff)
}

# congruence <- function(matrix1, matrix2, delta = 3, units = "hours") {
#   require(circular)
#   
#   rho1 <- matrix1$rho
#   rho2 <- matrix2$rho
#   phi1 <- matrix1$phi
#   phi2 <- matrix2$phi
#   
#   # Genes rhythmic in either 
#   rhythmic_union <- (rho1 == 1) | (rho2 == 1)
#   num_union <- sum(rhythmic_union)
#   
#   if (num_union == 0) {
#     return(list(
#       congruence_index = NA_real_,
#       gain_index = NA_real_,
#       loss_index = NA_real_,
#       Cp = NA_real_, 
#       Dp = NA_real_
#     ))
#   }
#   
#   # Conserved rhythmic genes (TP: rhythmic in both)
#   tp <- (rho1 == 1) & (rho2 == 1)
#   num_tp <- sum(tp)
#   
#   num_concordant <- 0
#   num_discordant <- 0
#   
#   if (num_tp > 0) {
#     phase_diffs <- phase_difference(phi1[tp], phi2[tp], units = units)
#     phase_distances <- abs(phase_diffs)
#     
#     # Concordant: small phase differences (<= delta)
#     concordant <- phase_distances <= delta
#     num_concordant <- sum(concordant, na.rm = TRUE)
#     
#     # Discordant: large phase differences (> delta)
#     num_discordant <- num_tp - num_concordant
#   }
#   
#   # Gain: rhythmic in matrix1 but not matrix2
#   num_gain <- sum(rho1 == 1 & rho2 == 0)
#   
#   # Loss: rhythmic in matrix2 but not matrix1
#   num_loss <- sum(rho1 == 0 & rho2 == 1)
#   
#   # Proportions over the union of rhythmic genes
#   conserved_concordant <- num_concordant / num_union
#   conserved_discordant <- num_discordant / num_union
#   loss <- num_loss / num_union
#   gain <- num_gain / num_union
#   
#   # Proportion of union that is conserved with similar phase
#   phase_jaccard <- conserved_concordant
#   
#   # Also compute standard rhythmicity Jaccard for reference
#   rhythm_jaccard <- num_tp / num_union
#   
#   return(list(
#     congruence_index = rhythm_jaccard,
#     gain_index = gain,
#     loss_index = loss,
#     Cp = conserved_concordant, 
#     Dp = conserved_discordant
#   ))
# }

#' Compute probabilistic congruence (c-score) between two conditions
#'
#' @title Probabilistic congruence index
#'
#' @description
#' Computes the expected Jaccard-based congruence index using posterior
#' rhythmicity probabilities p_A = rowMeans(rho_A) and
#' p_B = rowMeans(rho_B).  The c-score is defined as
#' E[intersection] / E[union] where
#' E[intersection] = sum(p_A * p_B) and
#' E[union] = sum(p_A) + sum(p_B) - E[intersection].
#' Gain and loss indices are the fractions of the expected union gained or
#' lost relative to condition A.
#'
#' @param matrix1 Named list with element \code{rho} (G x K matrix);
#'   MCMC posterior samples for condition 1 (reference / earlier time point).
#' @param matrix2 Named list with element \code{rho} (G x K matrix);
#'   MCMC posterior samples for condition 2.
#' @param delta Numeric; phase-shift threshold in \code{units} (currently
#'   unused in the probabilistic c-score; reserved for future use,
#'   default 3).
#' @param units Character; \code{"hours"} (default) or \code{"radians"}.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{congruence_index}{Numeric; E[intersection] / E[union].}
#'     \item{gain_index}{Numeric; expected gain / E[union].}
#'     \item{loss_index}{Numeric; expected loss / E[union].}
#'     \item{gain_loss_ratio}{Numeric; gain_index / loss_index.}
#'     \item{n_genes}{Integer; number of genes.}
#'   }
#'
congruence <- function(matrix1, matrix2, delta = 3, units = "hours") {
  if (!missing(delta) && delta != 3)
    warning("`delta` is not used by the probabilistic c-score (paper Eq. 3); it is reserved for a future phase-threshold variant.")
  if (!missing(units) && units != "hours")
    warning("`units` is not used by the probabilistic c-score; it is reserved for a future phase-threshold variant.")

  # Convert MCMC matrices to posterior probabilities
  p_A <- rowMeans(matrix1$rho)  # Condition 1 (e.g., younger)
  p_B <- rowMeans(matrix2$rho)  # Condition 2 (e.g., older)
  
  intersection_probs <- p_A * p_B
  expected_intersection <- sum(intersection_probs)
  expected_union <- sum(p_A) + sum(p_B) - expected_intersection
  
  if (expected_union == 0) {
    return(list(
      congruence_index = 0,
      gain_index = 0,
      loss_index = 0,
      gain_loss_ratio = NA_real_,
      n_genes = length(p_A)
    ))
  }
  
  # Gain: genes rhythmic in B but not A (gained rhythmicity)
  # Loss: genes rhythmic in A but not B (lost rhythmicity)
  expected_loss <- sum(p_A) - expected_intersection
  expected_gain <- sum(p_B) - expected_intersection
  
  # Normalize by expected union
  union_inv <- 1 / expected_union
  congruence_index <- expected_intersection * union_inv
  loss_index <- expected_loss * union_inv
  gain_index <- expected_gain * union_inv
  
  # Gain/Loss ratio — paper Eq. 4: GLR = Gain / Loss.
  # Inf when Loss=0 and Gain>0 is the mathematically correct value.
  gain_loss_ratio <- if (loss_index == 0) {
    if (gain_index == 0) {
      NA_real_  # Both zero — ratio undefined
    } else {
      Inf       # Loss=0, Gain>0 — paper Eq. 4 gives +Inf
    }
  } else {
    gain_index / loss_index
  }
  
  return(list(
    congruence_index = congruence_index,
    gain_index = gain_index,
    loss_index = loss_index,
    gain_loss_ratio = gain_loss_ratio,
    n_genes = length(p_A)
  ))
}

#' Compute pairwise concordance matrix across multiple tissues
#'
#' @description
#' For each pair of tissues in \code{human_data}, computes the Jaccard
#' concordance index between their posterior rhythmicity classifications.
#' Optionally restricts the comparison to the top \code{n_gene} genes
#' ranked by mean posterior rhythmicity across all tissues.
#'
#' @param human_data Named list of MCMC output lists; each element must
#'   have a \code{rho} matrix with matching gene row names.
#' @param n_gene Integer or \code{NULL}; if specified, restrict to the
#'   top \code{n_gene} genes by mean rhythmicity across tissues.
#'
#' @return Square symmetric numeric matrix of Jaccard concordance values,
#'   with row and column names equal to \code{names(human_data)}.
#'
#' @export
pairwise_concordance <- function(human_data, n_gene = NULL) {
  tissue_names <- names(human_data)
  
  # If n_gene is specified, compute overall mean of rho for each gene
  if (!is.null(n_gene)) {
    # Find common genes across all tissues
    common_genes <- Reduce(intersect, lapply(human_data, function(x) rownames(x$rho)))
    
    # Compute mean rho for each gene:
    gene_means <- sapply(common_genes, function(gene) {
      mean(sapply(human_data, function(x) as.numeric(x$rho[gene, 1])))
    })
    
    # Select the top n_gene genes based on \bar{\rho}(g)
    top_genes <- names(sort(gene_means, decreasing = TRUE))[1:n_gene]
  }
  
  n <- length(tissue_names)
  concordance_matrix <- matrix(NA, nrow = n, ncol = n, 
                               dimnames = list(tissue_names, tissue_names))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i < j) {
        rho1 <- human_data[[tissue_names[i]]]$rho
        rho2 <- human_data[[tissue_names[j]]]$rho
        
        # Restrict to common genes for the pair
        common_genes <- intersect(rownames(rho1), rownames(rho2))
        if (!is.null(n_gene)) {
          # Further restrict to the top genes
          common_genes <- intersect(common_genes, top_genes)
        }
        rho1 <- rho1[common_genes, , drop = FALSE]
        rho2 <- rho2[common_genes, , drop = FALSE]
        
        # Now concordance() returns just the Jaccard index
        jaccard <- concordance(rho1, rho2)
        concordance_matrix[i, j] <- jaccard
        concordance_matrix[j, i] <- jaccard
      } else if (i == j) {
        concordance_matrix[i, j] <- 1  # perfect concordance with self
      }
    }
  }
  
  return(concordance_matrix)
}



try_any_mirror <- function(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl") {
  mirrors <- c("useast", "www", "asia")  # List of mirrors to try
  ensembl <- NULL
  
  for (mirror in mirrors) {
    message(paste("Trying", mirror, "mirror..."))
    try({
      ensembl <- useEnsembl(biomart = biomart, dataset = dataset, mirror = mirror)
    }, silent = TRUE)
    
    if (!is.null(ensembl)) {
      message(paste("Connected to", mirror, "mirror."))
      return(ensembl)
    }
  }
  
  stop("Unable to connect to any Ensembl mirror.")
}


#' Summarise posterior rhythmicity and Bayes Factor per gene
#'
#' @description
#' Computes the posterior mean rhythmicity (\code{RowAverage}), the number
#' of zero samples (\code{Zeroes}), and the Bayes Factor for each gene.
#' The Bayes Factor is computed as posterior odds / prior odds:
#' BF = (rho_bar / (1 - rho_bar)) / (p_rhythmic / (1 - p_rhythmic)).
#' Genes with BF > the threshold are called rhythmic (Rhythmicity = 1).
#'
#' @param input_df G x K binary matrix of posterior rho samples.
#' @param BF Numeric; Bayes Factor threshold for calling a gene rhythmic.
#' @param p_rhythmic Numeric; prior probability of rhythmicity (default 0.5).
#'
#' @return Data.frame with columns \code{Zeroes}, \code{RowAverage},
#'   \code{Rhythmicity} (0/1), and \code{BayesF}.
#'
#' @export
summarize_bay <- function(input_df, BF, p_rhythmic = 0.2) {
  # Default 0.2 matches CB_MCMC_single_rj_slice default (paper §2.1).
  # BF = posterior_odds / prior_odds; using a different p_rhythmic here than
  # was passed to the MCMC will produce miscalibrated Bayes Factors.
  if (!isTRUE(all.equal(p_rhythmic, 0.2)))
    message("summarize_bay: p_rhythmic = ", p_rhythmic,
            ". Ensure this matches the value used in CB_MCMC_single_rj_slice ",
            "(default 0.2) or BF values will be miscalibrated.")

  row_counts  <- rowSums(input_df == 0)
  row_average <- rowMeans(input_df)

  prior_odds <- p_rhythmic / (1 - p_rhythmic + 1e-20)
  
  # Calculate the posterior odds
  posterior_odds <- row_average / (1 - row_average + 1e-20)
  
  # Calculate the Bayes Factor correctly: BF = posterior odds / prior odds
  BayesF <- posterior_odds / prior_odds
  
  # Decide if the gene is rhythmic or non-rhythmic based on the BF threshold
  rhythmicity <- ifelse(BayesF > BF, 1, 0)
  
  # Returning new data frame
  summary_df <- data.frame(
    Zeroes = row_counts,
    RowAverage = row_average,
    Rhythmicity = rhythmicity,
    BayesF = BayesF
  )
  
  return(summary_df)
}


#' Intersect phase matrices to common genes
#'
#' @description
#' Finds the common gene symbols between two MCMC output lists using the
#' \code{attr(rho, "symbols")} attribute and returns subsetted phase
#' matrices aligned to that common gene set.
#'
#' @param X_summary Named list; MCMC output for condition X with
#'   \code{attr(rho, "symbols")} set.
#' @param Y_summary Named list; MCMC output for condition Y.
#'
#' @return List with elements \code{X_c} and \code{Y_c} (phase matrices
#'   restricted to common genes, in matching row order).
#'
# Filter common genes between datasets
phi_filter <- function(X_summary, Y_summary) {
  common_genes <- intersect(attr(X_summary$rho, "symbols"), attr(Y_summary$rho, "symbols"))

  if (length(common_genes) == 0) {
    stop("No common genes found between datasets X_summary and Y_summary.")
  }

  # Subset the phase data to only include common genes
  X_c <- X_summary$phi[match(common_genes, attr(X_summary$rho, "symbols")), ]
  Y_c <- Y_summary$phi[match(common_genes, attr(Y_summary$rho, "symbols")), ]
  
  return(list(
    X_c = X_c,
    Y_c = Y_c
  )) 
}  

# Phase Inference 
phase_inf <- function(matrix1, matrix2,
                      P           = 24,
                      credMass    = 0.95,
                      shift       = 4,
                      a           = -P/2,
                      fdr_thresh  = 0.05,
                      p_rhythmic  = 0.20,
                      BF_cut      = 3) {
  
  phi_matrix1 = matrix1$phi
  phi_matrix2 = matrix2$phi
  rho_matrix1 = matrix1$rho
  rho_matrix2 = matrix2$rho
  
  
  # 1. Helper: circular median in [−P/2,P/2)
  circular_median <- function(samples, P) {
    rad <- circular((samples / P) * 2 * pi, units = "radians")
    med <- median.circular(rad, type = "median")
    t   <- (med / (2 * pi)) * P
    ((t + P/2) %% P) - P/2
  }
  
  circular_median_24h <- function(samples, P = 24) {
    rad <- circular((samples / P) * 2 * pi, units = "radians")
    med <- median.circular(rad, type = "median")
    t   <- (med / (2 * pi)) * P
    t %% P
  }
  
  
  # Harmonize genes
  common_genes <- Reduce(intersect, list(
    rownames(phi_matrix1), rownames(phi_matrix2),
    rownames(rho_matrix1), rownames(rho_matrix2)
  ))
  phi1 <- phi_matrix1
  phi2 <- phi_matrix2
  rho1 <- rho_matrix1
  rho2 <- rho_matrix2
  
  
  # Compute phase differences
  phase_diff_matrix <- ((phi1 - phi2 + P/2) %% P) - P/2
  peak1 <- apply(phi1, 1, function(x) circular_median_24h(x, P))
  peak2 <- apply(phi2, 1, function(x) circular_median_24h(x, P))
  
  # 2. Estimate Bayes Factor per gene from rho_matrix
  rho_bar1 <- rowMeans(rho1)
  rho_bar2 <- rowMeans(rho2)
  BF1 <- (rho_bar1 * (1 - p_rhythmic)) / ((1 - rho_bar1 + 1e-20) * p_rhythmic)
  BF2 <- (rho_bar2 * (1 - p_rhythmic)) / ((1 - rho_bar2 + 1e-20) * p_rhythmic)
  names(BF1) <- names(BF2) <- rownames(phase_diff_matrix)
  
  # 3. Analyze one test type ("difference" or "conservation")
  analyze_test <- function(test_type) {
    gene_names <- rownames(phase_diff_matrix)
    
    results <- lapply(seq_len(nrow(phase_diff_matrix)), function(i) {
      vec     <- phase_diff_matrix[i, ]
      hdi     <- circular_HDI(vec, credMass = credMass, P = P, a = a)
      phi.Est <- circular_median(vec, P)
      
      p_shift <- mean(abs(vec) >= shift)
      p_cons  <- 1 - p_shift
      
      # Calculate in_HDI once - same for both test types
      in_HDI <- in_circular_interval(hdi$lower, hdi$upper, 0, P, a)
      
      # Determine which probability to use for this test
      if (test_type == "difference") {
        test_prob <- p_shift
      } else {
        test_prob <- p_cons
      }
      
      hdi_width <- if (hdi$upper >= hdi$lower) {
        hdi$upper - hdi$lower
      } else {
        (P/2 - hdi$lower) + (hdi$upper - (-P/2))
      }
      
      # Calculate log difference
      log_BF_diff <- log(BF2[gene_names[i]]) - log(BF1[gene_names[i]])
      
      
      # Determine gain/loss status
      gain_loss <- ifelse(log_BF_diff < 0, "Loss",
                          ifelse(log_BF_diff > 0, "Gain", "Neutral"))
      
      data.frame(
        gene            = gene_names[i],
        phi1.Est           = peak1[i],
        phi2.Est          = peak2[i],
        phi.Est         = phi.Est,
        phi.Lower       = hdi$lower,
        phi.Upper       = hdi$upper,
        hdi.width       = hdi_width,
        test_prob       = test_prob,
        prob_shift      = p_shift,
        prob_conserved  = p_cons,
        in_HDI          = in_HDI,  
        BF1             = BF1[i],
        BF2             = BF2[i],
        log_BF_diff     = log_BF_diff,
        Gain_Loss       = gain_loss,
        stringsAsFactors = FALSE
      )
    })
    
    df <- do.call(rbind, results)
    
    # Order by the relevant probability for this test
    ord <- order(-df$test_prob)
    df  <- df[ord, ]
    
    # Calculate BFDR using the relevant probability
    df$BFDR <- cumsum(1 - df$test_prob) / seq_len(nrow(df))
    
    # Set flags based on test type 
    if (test_type == "difference") {
      df$flag <- (!df$in_HDI &                    # Zero NOT in HDI = phase difference
                    df$test_prob > 0.5 &
                    df$BFDR <= fdr_thresh &
                    df$BF1 > BF_cut &
                    df$BF2 > BF_cut)
    } else {
      df$flag <- (df$in_HDI &                     # Zero IN HDI = phase conservation
                    df$test_prob > 0.5 &
                    df$BFDR <= fdr_thresh &
                    df$BF1 > BF_cut &
                    df$BF2 > BF_cut)
    }
    
    return(df)
  }
  
  list(
    phase_difference   = analyze_test("difference"),
    phase_conservation = analyze_test("conservation")
  )
}


pathShiftEnrich <- function(phi1, phi2, pathway.list,
                                 P = 24, shift = 2, credMass = 0.95,
                                 bfdr_thresh = 0.1, prob_shift_thresh = 0.9,
                                 min_overlap = 5, minRecurrence = 1) {
  
  # Step 1: Compute circular phase difference matrix
  common_genes <- intersect(rownames(phi1), rownames(phi2))
  phi1 <- phi1[common_genes, , drop = FALSE]
  phi2 <- phi2[common_genes, , drop = FALSE]
  
  phase_diff_mat <- (phi1 - phi2 + P/2) %% P - P/2  # wrap into (-P/2, P/2]
  rownames(phase_diff_mat) <- common_genes
  
  # Step 2: Get significance for each gene using phase_inf
  phase_df <- phase_inf(phase_diff_matrix = phase_diff_mat, P = P, shift = shift, credMass = credMass)
  sig_genes <- with(phase_df, gene[prob_shift >= prob_shift_thresh & BFDR <= bfdr_thresh & !in_HDI])
  
  # Step 3: Pathway-level summarization
  result <- lapply(names(pathway.list), function(pw) {
    genes_in_path <- intersect(rownames(phase_diff_mat), pathway.list[[pw]])
    sig_in_path <- intersect(sig_genes, genes_in_path)
    n_overlap <- length(genes_in_path)
    n_sig <- length(sig_in_path)
    
    data.frame(
      Pathway = pw,
      PathwaySize = length(pathway.list[[pw]]),
      Overlap = n_overlap,
      ShiftedCount = n_sig,
      ShiftedProp = ifelse(n_overlap > 0, n_sig / n_overlap, NA)
    )
  })
  
  df <- do.call(rbind, result)
  
  # Step 4: Filter
  df <- subset(df, Overlap >= min_overlap & ShiftedCount >= minRecurrence)
  df <- df[order(-df$ShiftedCount, -df$ShiftedProp), ]
  rownames(df) <- NULL
  
  return(df)
}

############### PLOTTING

#' Plot posterior phase distributions for selected genes
#'
#' @description
#' Produces a faceted density plot of posterior acrophase samples for a
#' set of genes across two species/tissues.  One panel per gene x species
#' combination; saves a PNG to \code{output_dir}.
#'
#' @param gene_aliases Character vector; gene names to plot (case-insensitive
#'   match against row names of \code{phi_mat1} and \code{phi_mat2}).
#' @param phi_mat1 G x K matrix; posterior phase samples for condition 1.
#' @param phi_mat2 G x K matrix; posterior phase samples for condition 2.
#' @param species_names Character vector of length 2; labels for conditions
#'   (default \code{c("Species1", "Species2")}).
#' @param tissue_names Character vector of length 2; tissue labels used in
#'   the plot subtitle and filename.
#' @param output_dir Character; directory in which to save the PNG file.
#'
#' @return Called for side effects; invisibly returns \code{NULL}.  Saves
#'   a PNG named \code{<genes>_<tissue1>_<tissue2>_posterior_phase.png}.
#'
plotGenePosteriorPhase <- function(gene_aliases,
                                   phi_mat1, phi_mat2,
                                   species_names = c("Species1", "Species2"),
                                   tissue_names = c("Tissue1", "Tissue2"),
                                   output_dir) {
  # Helper: match gene row index (case-insensitive)
  match_gene <- function(mat, gene) {
    idx <- which(tolower(rownames(mat)) == tolower(gene))
    if (length(idx) == 0) return(NA_integer_) else return(idx[1])
  }
  
  # Initialize list to store data
  df_list <- list()
  
  for (gene in gene_aliases) {
    idx1 <- match_gene(phi_mat1, gene)
    idx2 <- match_gene(phi_mat2, gene)
    
    if (is.na(idx1) || is.na(idx2)) {
      message(paste(gene, "not found in one or both datasets. Skipping."))
      next
    }
    
    # Extract samples and filter NAs
    df1 <- tibble(
      phi = as.numeric(phi_mat1[idx1, ]),
      species = species_names[1],
      gene = gene,
      tissue = tissue_names[1]
    ) %>% filter(is.finite(phi))
    
    df2 <- tibble(
      phi = as.numeric(phi_mat2[idx2, ]),
      species = species_names[2],
      gene = gene,
      tissue = tissue_names[2]
    ) %>% filter(is.finite(phi))
    
    df_gene <- bind_rows(df1, df2) %>%
      mutate(phi_hours = ((12 / pi) * phi) %% 24)
    
    df_list[[gene]] <- df_gene
  }
  
  if (length(df_list) == 0) {
    message("None of the specified genes were found in both datasets.")
    return(NULL)
  }
  
  # Combine all gene data
  df_combined <- bind_rows(df_list)
  df_combined$species <- factor(df_combined$species, levels = species_names)
  df_combined$gene <- factor(df_combined$gene, levels = gene_aliases)
  
  # Generate color palette for genes
  n_colors <- length(gene_aliases)
  palette <- brewer.pal(min(n_colors, 8), "Set2")
  if (n_colors > 8) {
    palette <- colorRampPalette(palette)(n_colors)
  }
  
  # Create density plot
  p <- ggplot(df_combined, aes(x = phi_hours, fill = gene)) +
    geom_density(alpha = 0.4) +
    facet_grid(gene ~ species) +
    scale_fill_manual(values = palette) +
    scale_x_continuous(breaks = seq(0, 24, 6), limits = c(0, 24)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = "Posterior Phase Distribution",
         subtitle = paste(tissue_names[1], "&", tissue_names[2]),
         x = "Phase (hours)", y = "Density", fill = "Gene") +
    theme_minimal(base_size = 14) +
    theme_minimal(base_size = 14) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_line(color = "gray80", linewidth = 0.5),
      panel.grid.minor = element_line(color = "gray90", linewidth = 0.25),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank()
    )
  
  
  # Create filename and save
  gene_part <- paste(gene_aliases, collapse = "_")
  filename <- paste0(gene_part, "_", tissue_names[1], "_", tissue_names[2], "_posterior_phase.png")
  output_path <- file.path(output_dir, filename)
  
  ggsave(output_path, plot = p, width = 15, height = 10, dpi = 300)
  message("Plot saved to: ", output_path)
}

############### PLOTTING

#' Plot core clock gene median phases on a polar axis
#'
#' @description
#' Computes the circular median phase for each core clock gene and draws a
#' polar bar chart using ggplot2, saving to \code{output_path}.
#'
#' @param phi_mat G x K matrix; posterior phase samples.
#' @param gene_aliases Character vector; core clock gene names to plot.
#' @param tissue_label Character; plot title / tissue label.
#' @param output_path Character; file path for the saved PNG.
#' @param P Numeric; period in hours (default 24).
#'
#' @return Called for side effects; saves a PNG and returns \code{NULL}.
#'
plotCoreClockMedians <- function(phi_mat, gene_aliases, tissue_label, output_path, P = 24) {
  # Step 1: Compute median phase for core clock genes
  get_median_phase <- function(phi_mat, gene) {
    i <- which(tolower(rownames(phi_mat)) == tolower(gene))
    if (length(i) == 0) return(NA)
    samples <- ((12 / pi) * as.numeric(phi_mat[i, ])) %% P
    median_phase <- circular_median(samples, P = P)
    return(median_phase)
  }
  
  # Step 2: Get median for each core clock gene
  df <- purrr::map_dfr(gene_aliases, function(g) {
    median_phase <- get_median_phase(phi_mat, g)
    tibble(
      gene = toupper(g),
      median = median_phase
    )
  })
  
  # Step 3: Prepare plot dataframe
  df$gene <- factor(df$gene, levels = toupper(gene_aliases))
  df <- df %>% mutate(
    ymin = 0.5,
    ymax = 1,
    theta = 2 * pi * (1 - as.numeric(median) / P)
  )
  
  # Step 4: Plot with ggplot2
  p <- ggplot(df, aes(x = theta, y = ymin, fill = gene)) +
    geom_bar(stat = "identity", width = 1, color = "black", alpha = 0.7) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(length(gene_aliases), "Set1")) +
    coord_polar(start = -pi/2) +
    scale_x_continuous(breaks = seq(0, P, by = 3), labels = paste0(seq(0, P, by = 3), "h")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.2)) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(title = tissue_label, x = "", y = "", fill = NULL)
  
  # Step 5: Save plot
  ggsave(output_path, plot = p, width = 8, height = 8, dpi = 300)
  message("Core clock gene median plot saved to: ", output_path)
}


#' Plot circular HDI arcs for core clock genes on a polar plot
#'
#' @description
#' Computes the circular highest-density interval for each core clock gene
#' using \code{circular_HDI} and renders coloured arcs on a 24-hour polar
#' axis, saving to \code{output_path}.
#'
#' @param phi_mat G x K matrix; posterior phase samples.
#' @param gene_aliases Character vector; gene names to plot.
#' @param tissue_label Character; plot title / tissue label.
#' @param output_path Character; file path for the saved PNG.
#' @param credMass Numeric; HDI coverage (default 0.95).
#' @param P Numeric; period in hours (default 24).
#'
#' @return Called for side effects; saves a PNG and returns \code{NULL}.
#'
plotHDIClockPolar <- function(phi_mat, gene_aliases, tissue_label, output_path,
                              credMass = 0.95, P = 24) {
  # Build HDI data frame using your own circular_HDI
  hdi_df <- purrr::map_dfr(gene_aliases, function(g) {
    i <- which(tolower(rownames(phi_mat)) == tolower(g))
    if (length(i) == 0) return(NULL)
    samples <- ((12 / pi) * as.numeric(phi_mat[i, ])) %% P
    hdi <- circular_HDI(samples, credMass = credMass, P = P)
    tibble(
      gene = toupper(g),
      phi.Lower = hdi$lower,
      phi.Upper = hdi$upper
    )
  })
  
  if (nrow(hdi_df) == 0) {
    warning("No valid genes found.")
    return(NULL)
  }
  
  hdi_df$gene <- factor(hdi_df$gene, levels = toupper(gene_aliases))
  
  # Plot
  p <- ggplot(hdi_df, aes(xmin = phi.Lower, xmax = phi.Upper, ymin = 0, ymax = 1, fill = gene)) +
    geom_rect(color = "grey", alpha = 0.1) +
    scale_x_continuous("", limits = c(0, P), breaks = 0:P, labels = paste0(0:P, "h")) +
    coord_polar(start = 0) +
    theme_minimal() +
    labs(title = paste("     ", tissue_label), x = "", y = "") +
    theme(
      axis.text.y = element_blank(),
      legend.title = element_blank(),
      plot.margin = unit(c(1, 0, 0, 0), "cm"),
      plot.title = element_text(hjust = 0.5)
    )
  
  ggsave(output_path, plot = p, width = 5.5, height = 5.5, dpi = 300)
  message("Polar HDI plot saved to: ", output_path)
}
################################################################################
# BFDR (Bayesian False Discovery Rate) Functions
################################################################################
# Main BFDR function - calculates threshold from posterior probabilities
################################################################################
# Bayesian False Discovery Rate (BFDR) threshold estimation
# Implements Eq. (25): τ_c = max{ τ : BFDR(τ) ≤ α }
################################################################################
#' Estimate the Bayesian False Discovery Rate threshold from posterior probabilities
#'
#' @title BFDR threshold from posterior rhythmicity probabilities
#'
#' @description
#' Implements the BFDR rule described in the BayRC paper (Eq. 25):
#' sort genes by decreasing posterior probability p_g, compute the running
#' average of (1 - p_g), and return the largest threshold tau_c such that
#' BFDR(tau_c) = mean(1 - p_g | p_g >= tau_c) <= alpha.
#' This controls the expected proportion of false rhythmic calls at level alpha.
#'
#' @param posterior_probs Numeric vector; marginal posterior probability of
#'   rhythmicity for each gene (values in \[0, 1\]).
#' @param alpha Numeric; desired BFDR level (default 0.05).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{threshold}{Numeric; the BFDR-optimal cutoff tau_c.}
#'     \item{rhythmic_genes}{Logical vector; \code{TRUE} for called rhythmic.}
#'     \item{n_rhythmic}{Integer; number of rhythmic calls.}
#'     \item{bfdr_values}{Numeric vector; BFDR at each ranked position.}
#'     \item{sorted_probs}{Sorted (descending) posterior probabilities.}
#'   }
#'
#' @export
bfdr_from_posterior <- function(posterior_probs, alpha = 0.05) {
  # Sort posterior probabilities in descending order
  ord <- order(posterior_probs, decreasing = TRUE)
  sorted_probs <- posterior_probs[ord]
  n <- length(sorted_probs)
  
  # Expected false discoveries and BFDR curve
  cumsum_false <- cumsum(1 - sorted_probs)
  discoveries <- seq_len(n)
  bfdr_values <- cumsum_false / discoveries
  
  # Identify largest index satisfying BFDR ≤ α
  valid_idx <- which(bfdr_values <= alpha)
  
  if (length(valid_idx) == 0) {
    return(list(
      threshold = 1.0,
      rhythmic_genes = rep(FALSE, n),
      n_rhythmic = 0,
      bfdr_values = bfdr_values,
      sorted_probs = sorted_probs
    ))
  }
  
  max_idx <- max(valid_idx)
  threshold <- sorted_probs[max_idx]
  
  rhythmic_genes <- posterior_probs >= threshold
  
  list(
    threshold = threshold,
    rhythmic_genes = rhythmic_genes,
    n_rhythmic = sum(rhythmic_genes),
    bfdr_values = bfdr_values,
    sorted_probs = sorted_probs
  )
}


################################################################################
# Bayesian rhythmic-gene detection with BFDR control
################################################################################
# Step 1: Bayesian rhythmic-gene detection with BFDR control
################################################################################

#' Detect rhythmic genes in two conditions using BFDR control
#'
#' @description
#' Applies \code{bfdr_from_posterior} separately to the marginal posterior
#' rhythmicity probabilities of two conditions and returns the BFDR-optimal
#' rhythmic gene sets for each, together with data.frames of annotated gene
#' information.
#'
#' @param dat1 Named list; MCMC output for condition A with a \code{rho}
#'   matrix.
#' @param dat2 Named list; MCMC output for condition B with a \code{rho}
#'   matrix.
#' @param bfdr_alpha Numeric; BFDR level (default 0.05).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{rhythmic_A_logical, rhythmic_B_logical}{Logical vectors for
#'       conditions A and B.}
#'     \item{rhythmic_A, rhythmic_B}{Data.frames of rhythmic genes with
#'       posterior probability columns.}
#'     \item{threshold_A, threshold_B}{BFDR thresholds for each condition.}
#'     \item{n_rhythmic_A, n_rhythmic_B, n_total, bfdr_alpha}{Counts and
#'       settings.}
#'   }
#'
#' @export
detect_rhy <- function(dat1, dat2, bfdr_alpha = 0.05) {
  # 1. Posterior rhythmicity probabilities
  p_A_all <- rowMeans(dat1$rho)  # p_g^(A)
  p_B_all <- rowMeans(dat2$rho)  # p_g^(B)
  
  # 2. Gene identifiers
  gene_names <- attr(dat1$rho, "symbols")
  if (is.null(gene_names))
    gene_names <- rownames(dat1$rho)
  if (is.null(gene_names))
    gene_names <- paste0("Gene_", seq_along(p_A_all))
  
  # 3. Determine thresholds τ_A, τ_B via BFDR rule (Eq. 25)
  bfdr_A <- bfdr_from_posterior(p_A_all, alpha = bfdr_alpha)
  bfdr_B <- bfdr_from_posterior(p_B_all, alpha = bfdr_alpha)
  tau_A <- bfdr_A$threshold
  tau_B <- bfdr_B$threshold
  
  # 4. Bayes-optimal classifications (Eq. 26)
  rhythmic_A <- p_A_all >= tau_A
  rhythmic_B <- p_B_all >= tau_B
  rhythmic_either <- rhythmic_A | rhythmic_B
  rhythmic_both   <- rhythmic_A & rhythmic_B
  
  # 5. Annotated data frames
  df_rhythmic_A <- data.frame(
    gene = gene_names[rhythmic_A],
    posterior_A = p_A_all[rhythmic_A],
    posterior_B = p_B_all[rhythmic_A],
    stringsAsFactors = FALSE
  )
  df_rhythmic_B <- data.frame(
    gene = gene_names[rhythmic_B],
    posterior_A = p_A_all[rhythmic_B],
    posterior_B = p_B_all[rhythmic_B],
    stringsAsFactors = FALSE
  )

  # 6. Structured output
  list(
    rhythmic_A_logical = rhythmic_A,
    rhythmic_B_logical = rhythmic_B,
    rhythmic_A = df_rhythmic_A,
    rhythmic_B = df_rhythmic_B,
    threshold_A = tau_A,
    threshold_B = tau_B,
    n_rhythmic_A = sum(rhythmic_A),
    n_rhythmic_B = sum(rhythmic_B),
    n_total = length(p_A_all),
    bfdr_alpha = bfdr_alpha
  )
}

################################################################################
# Phase analysis with Bayesian rhythmicity & transition classification
################################################################################
################################################################################
# Step 2: Bayesian transition classification (gain / loss / maintained)
################################################################################
#' Classify rhythmicity transitions using joint BFDR on transition posteriors
#'
#' @title Bayesian transition classification (gain / loss / maintained)
#'
#' @description
#' Computes joint transition posterior probabilities:
#' p_gain = (1 - p_A) * p_B (gained rhythmicity),
#' p_loss = p_A * (1 - p_B) (lost rhythmicity),
#' p_cons = p_A * p_B (maintained rhythmicity),
#' then applies \code{bfdr_from_posterior} to each to identify gain, loss,
#' and conserved gene sets at the specified BFDR level.
#'
#' @param pA Numeric vector; marginal posterior rhythmicity probability for
#'   condition A (reference; length G).
#' @param pB Numeric vector; marginal posterior rhythmicity probability for
#'   condition B (length G).
#' @param bfdr_alpha Numeric; BFDR control level (default 0.05).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{gain_loss_status}{Named character vector of length G with
#'       values \code{"Gain"}, \code{"Loss"}, \code{"Maintained"}, or
#'       \code{"Non-rhythmic"}.}
#'     \item{gain_genes, loss_genes, cons_genes}{Integer indices of genes
#'       in each class.}
#'     \item{tau_gain, tau_loss, tau_cons}{BFDR thresholds for each
#'       transition type.}
#'     \item{n_gain, n_loss, n_cons}{Counts.}
#'     \item{results}{Data.frame with per-gene classification details.}
#'   }
#'
#' @export
transition_classify <- function(pA, pB, bfdr_alpha = 0.05) {
  
  n_genes <- length(pA)
  gene_names <- names(pA)
  if (is.null(gene_names)) gene_names <- paste0("Gene_", seq_len(n_genes))
  
  # 1. Transition posteriors
  p_gain <- (1 - pA) * pB   # gained rhythmicity
  p_loss <- pA * (1 - pB)   # lost rhythmicity
  p_cons <- pA * pB         # maintained rhythmicity
  
  # 2. Condition-specific BFDR thresholds (paper §2.2, τ_A and τ_B)
  bfdr_A <- bfdr_from_posterior(pA, alpha = bfdr_alpha)
  bfdr_B <- bfdr_from_posterior(pB, alpha = bfdr_alpha)
  tau_rhythmic_A <- bfdr_A$threshold
  tau_rhythmic_B <- bfdr_B$threshold

  # 3. Apply BFDR to each transition type (paper §2.2, τ_gain, τ_loss, τ_cons)
  bfdr_gain <- bfdr_from_posterior(p_gain, alpha = bfdr_alpha)
  bfdr_loss <- bfdr_from_posterior(p_loss, alpha = bfdr_alpha)
  bfdr_cons <- bfdr_from_posterior(p_cons, alpha = bfdr_alpha)

  tau_gain <- bfdr_gain$threshold
  tau_loss <- bfdr_loss$threshold
  tau_cons <- bfdr_cons$threshold
  
  # 3. Classify genes
  gain_genes <- which(p_gain >= tau_gain)
  loss_genes <- which(p_loss >= tau_loss)
  cons_genes <- which(p_cons >= tau_cons)
  
  gain_loss_status <- rep("Non-rhythmic", n_genes)
  names(gain_loss_status) <- gene_names
  gain_loss_status[gain_genes] <- "Gain"
  gain_loss_status[loss_genes] <- "Loss"
  gain_loss_status[cons_genes] <- "Maintained"
  
  # Summary counts
  n_gain <- length(gain_genes)
  n_loss <- length(loss_genes)
  n_cons <- length(cons_genes)
  
  
  cat("=== Transition-level BFDR results ===\n")
  cat("τ_gain =", round(tau_gain, 3), "| n_gain =", n_gain, "\n")
  cat("τ_loss =", round(tau_loss, 3), "| n_loss =", n_loss, "\n")
  cat("τ_cons =", round(tau_cons, 3), "| n_cons =", n_cons, "\n\n")
  
  
  results_df <- data.frame(
    gene = gene_names,
    pA = pA,
    pB = pB,
    p_gain = p_gain,
    p_loss = p_loss,
    p_cons = p_cons,
    classification = gain_loss_status,
    stringsAsFactors = FALSE
  )
  
  
  
  list(
    # Condition-specific thresholds (paper §2.2 τ_A, τ_B)
    tau_rhythmic_A = tau_rhythmic_A,
    tau_rhythmic_B = tau_rhythmic_B,
    # Transition thresholds (paper §2.2 τ_gain, τ_loss, τ_cons)
    tau_gain = tau_gain,
    tau_loss = tau_loss,
    tau_cons = tau_cons,
    # Transition posterior vectors (paper §2.2 p_gain, p_loss, p_cons)
    p_gain = p_gain,
    p_loss = p_loss,
    p_cons = p_cons,
    gain_genes = gain_genes,
    loss_genes = loss_genes,
    cons_genes = cons_genes,
    n_gain = length(gain_genes),
    n_loss = length(loss_genes),
    n_cons = length(cons_genes),
    gain_loss_status = gain_loss_status,
    results = results_df
  )
}

# Marginal variant: applies BFDR separately to each condition's marginal
# rhythmicity probability (pA, pB), then classifies based on individual
# thresholds. Contrast with transition_classify() which uses joint
# transition probabilities (pA*pB, etc.).
#' Classify rhythmicity transitions using marginal BFDR on each condition
#'
#' @description
#' Applies \code{bfdr_from_posterior} separately to p_A and p_B (marginal
#' rhythmicity probabilities) to define tau_A and tau_B, then classifies
#' genes as Gain, Loss, Maintained, or Non-rhythmic based on whether each
#' gene exceeds its condition's threshold.  Contrast with
#' \code{transition_classify}, which uses joint transition probabilities.
#'
#' @param pA Numeric vector; marginal posterior rhythmicity probability for
#'   condition A (length G).
#' @param pB Numeric vector; marginal posterior rhythmicity probability for
#'   condition B (length G).
#' @param bfdr_alpha Numeric; BFDR control level (default 0.05).
#'
#' @return A list with the same structure as \code{transition_classify}
#'   but with \code{tau_A} and \code{tau_B} instead of \code{tau_gain/
#'   loss/cons}, and a \code{results} data.frame with columns
#'   \code{rhythmic_A} and \code{rhythmic_B}.
#'
#' @export
transition_classify_marginal <- function(pA, pB, bfdr_alpha = 0.05) {

  n_genes <- length(pA)
  gene_names <- names(pA)
  if (is.null(gene_names)) gene_names <- paste0("Gene_", seq_len(n_genes))

  bfdr_A <- bfdr_from_posterior(pA, alpha = bfdr_alpha)
  bfdr_B <- bfdr_from_posterior(pB, alpha = bfdr_alpha)

  tau_A <- bfdr_A$threshold
  tau_B <- bfdr_B$threshold

  rhythmic_A <- pA >= tau_A
  rhythmic_B <- pB >= tau_B

  gain_genes  <- which(!rhythmic_A &  rhythmic_B)
  loss_genes  <- which( rhythmic_A & !rhythmic_B)
  cons_genes  <- which( rhythmic_A &  rhythmic_B)

  gain_loss_status <- rep("Non-rhythmic", n_genes)
  names(gain_loss_status) <- gene_names
  gain_loss_status[gain_genes] <- "Gain"
  gain_loss_status[loss_genes] <- "Loss"
  gain_loss_status[cons_genes] <- "Maintained"

  cat("=== Marginal BFDR results ===\n")
  cat("τ_A =", round(tau_A, 3), "| n_rhythmic_A =", sum(rhythmic_A), "\n")
  cat("τ_B =", round(tau_B, 3), "| n_rhythmic_B =", sum(rhythmic_B), "\n")
  cat("n_gain =", length(gain_genes),
      "| n_loss =", length(loss_genes),
      "| n_cons =", length(cons_genes), "\n\n")

  results_df <- data.frame(
    gene           = gene_names,
    pA             = pA,
    pB             = pB,
    rhythmic_A     = rhythmic_A,
    rhythmic_B     = rhythmic_B,
    classification = gain_loss_status,
    stringsAsFactors = FALSE
  )

  list(
    tau_A            = tau_A,
    tau_B            = tau_B,
    gain_genes       = gain_genes,
    loss_genes       = loss_genes,
    cons_genes       = cons_genes,
    n_gain           = length(gain_genes),
    n_loss           = length(loss_genes),
    n_cons           = length(cons_genes),
    gain_loss_status = gain_loss_status,
    results          = results_df
  )
}

################################################################################
# Step 3: Phase inference (Δφ estimation, shift/conservation BFDR)
################################################################################
# phase_infer <- function(phi_matrix1, phi_matrix2, gain_loss_status,
#                         P = 24, credMass = 0.95, shift = 4,
#                         a = -P/2, bfdr_alpha = 0.05) {
#   require(circular)
#   n_genes <- nrow(phi_matrix1)
#   gene_names <- rownames(phi_matrix1)
#   
#   # helper functions -----------------------------------------------------------
#   circ_med <- function(x, P) {
#     r <- circular((x / P) * 2 * pi, units = "radians")
#     t <- (median.circular(r, type = "median") / (2 * pi)) * P
#     ((t + P / 2) %% P) - P / 2
#   }
#   circ_med_24 <- function(x, P = 24) {
#     r <- circular((x / P) * 2 * pi, units = "radians")
#     t <- (median.circular(r, type = "median") / (2 * pi)) * P
#     t %% P
#   }
#   
#   # phase medians --------------------------------------------------------------
#   peak1 <- apply(phi_matrix1, 1, circ_med_24, P)
#   peak2 <- apply(phi_matrix2, 1, circ_med_24, P)
#   phase_diff_matrix <- ((phi_matrix1 - phi_matrix2 + P / 2) %% P) - P / 2
#   
#   maintained_idx <- which(gain_loss_status == "Maintained")
#   deltaPhi.Est <- deltaPhi.Lower <- deltaPhi.Upper <-
#     hdi.width <- prob_shift <- prob_conserved <- in_HDI <- rep(NA_real_, n_genes)
#   
#   # progress bar ---------------------------------------------------------------
#   if (length(maintained_idx) > 0) {
#     message("Computing HDI and phase metrics for maintained genes...")
#     pb <- txtProgressBar(min = 0, max = length(maintained_idx), style = 3)
#     
#     for (k in seq_along(maintained_idx)) {
#       i <- maintained_idx[k]
#       vec <- phase_diff_matrix[i, ]
#       
#       hdi <- circular_HDI(vec, credMass, P, a)
#       phi_est <- circ_med(vec, P)
#       
#       prob_shift[i] <- mean(abs(vec) >= shift)
#       prob_conserved[i] <- 1 - prob_shift[i]
#       in_HDI[i] <- in_circular_interval(hdi$lower, hdi$upper, 0, P, a)
#       
#       deltaPhi.Est[i] <- phi_est
#       deltaPhi.Lower[i] <- hdi$lower
#       deltaPhi.Upper[i] <- hdi$upper
#       
#       hdi.width[i] <- if (hdi$upper >= hdi$lower)
#         hdi$upper - hdi$lower else (P / 2 - hdi$lower) + (hdi$upper + P / 2)
#       
#       setTxtProgressBar(pb, k)
#     }
#     close(pb)
#   }
#   
#   # BFDR control for phase-shift / conserved -----------------------------------
#   flag_shift <- rep(FALSE, n_genes)
#   flag_cons  <- rep(FALSE, n_genes)
#   BFDR_shift_vec <- rep(NA_real_, n_genes)
#   BFDR_cons_vec  <- rep(NA_real_, n_genes)
#   
#   if (length(maintained_idx) > 0) {
#     ps  <- prob_shift[maintained_idx]
#     pcs <- prob_conserved[maintained_idx]
#     
#     # ----- SHIFT -----
#     ord_shift <- order(-ps)
#     BFDR_shift_sorted <- cumsum(1 - ps[ord_shift]) / seq_along(ps[ord_shift])
#     sel_shift <- (ps[ord_shift] > 0.5 & BFDR_shift_sorted <= bfdr_alpha)
#     
#     # remap back to original order
#     BFDR_shift_vec[maintained_idx[ord_shift]] <- BFDR_shift_sorted
#     flag_shift[maintained_idx[ord_shift[sel_shift]]] <- TRUE
#     
#     # ----- CONSERVED -----
#     ord_cons <- order(-pcs)
#     BFDR_cons_sorted <- cumsum(1 - pcs[ord_cons]) / seq_along(pcs[ord_cons])
#     sel_cons <- (pcs[ord_cons] > 0.5 & BFDR_cons_sorted <= bfdr_alpha)
#     
#     # remap back to original order
#     BFDR_cons_vec[maintained_idx[ord_cons]] <- BFDR_cons_sorted
#     flag_cons[maintained_idx[ord_cons[sel_cons]]] <- TRUE
#   }
#   
#   # summary print --------------------------------------------------------------
#   cat("\n=== PHASE INFERENCE SUMMARY ===\n")
#   cat("Maintained genes:", length(maintained_idx), "\n")
#   cat("BFDR α =", bfdr_alpha, " | shift threshold =", shift, "h\n")
#   cat("Significant phase-shifted genes:", sum(flag_shift, na.rm = TRUE), "\n")
#   cat("Significant phase-conserved genes:", sum(flag_cons, na.rm = TRUE), "\n\n")
#   
#   
#   
#   # output ---------------------------------------------------------------------
#   list(
#     peak1 = peak1,
#     peak2 = peak2,
#     deltaPhi.Est = deltaPhi.Est,
#     deltaPhi.Lower = deltaPhi.Lower,
#     deltaPhi.Upper = deltaPhi.Upper,
#     prob_shift = prob_shift,
#     flag_shift = flag_shift,
#     BFDR_shift = BFDR_shift_vec,
#     prob_conserved = prob_conserved,
#     flag_cons = flag_cons,
#     BFDR_cons  = BFDR_cons_vec,
#     hdi.width = hdi.width,
#     in_HDI = in_HDI
#   )
# }
# 
# 
# phase_analysis <- function(matrix1, matrix2,
#                            P = 24, credMass = 0.95,
#                            shift = 4, a = -P/2,
#                            bfdr_alpha = 0.05) {
#   rhy <- detect_rhy(matrix1, matrix2, bfdr_alpha)
#   pA <- rowMeans(matrix1$rho)
#   pB <- rowMeans(matrix2$rho)
#   trans <- transition_classify(pA, pB, bfdr_alpha)
#   
#   gene_names <- rownames(matrix1$phi)
#   gain_loss_status <- rep("Non-rhythmic", length(gene_names))
#   names(gain_loss_status) <- gene_names
#   gain_loss_status[trans$gain_genes] <- "Gain"
#   gain_loss_status[trans$loss_genes] <- "Loss"
#   gain_loss_status[trans$cons_genes] <- "Maintained"
#   
#   # updated call — rhy removed
#   phase <- phase_infer(matrix1$phi, matrix2$phi, gain_loss_status,
#                        P, credMass, shift, a, bfdr_alpha)
#   
#   list(
#     rhythmic_summary   = rhy,
#     transition_summary = trans,
#     phase_summary      = phase
#   )
# }
#

#' Infer phase shifts and conservation among maintained rhythmic genes
#'
#' @title Bayesian phase inference with BFDR-controlled shift/conservation flags
#'
#' @description
#' For genes classified as "Maintained" by \code{transition_classify} or
#' \code{transition_classify_marginal}, computes the posterior distribution
#' of the phase difference delta_phi = phi1 - phi2 (wrapped to [-P/2, P/2)),
#' estimates the probability of a significant shift (|delta_phi| >= shift),
#' and applies BFDR control to flag significantly shifted or conserved genes.
#' Optionally computes the full circular HDI of delta_phi when
#' \code{compute_hdi = TRUE}.
#'
#' @param phi_matrix1 G x K matrix; posterior phase samples for condition 1.
#' @param phi_matrix2 G x K matrix; posterior phase samples for condition 2.
#' @param gain_loss_status Named character vector of length G with values
#'   \code{"Gain"}, \code{"Loss"}, \code{"Maintained"}, or
#'   \code{"Non-rhythmic"} (from \code{transition_classify}).
#' @param P Numeric; period in hours (default 24).
#' @param credMass Numeric; HDI coverage (default 0.95); only used when
#'   \code{compute_hdi = TRUE}.
#' @param shift Numeric; phase-shift threshold in hours (default 4).
#' @param a Numeric; left endpoint for circular normalisation (default
#'   \code{-P/2}).
#' @param bfdr_alpha Numeric; BFDR control level (default 0.05).
#' @param compute_hdi Logical; if \code{TRUE}, compute the circular HDI
#'   for each maintained gene (slower; default \code{FALSE}).
#'
#' @return A list with per-gene vectors:
#'   \describe{
#'     \item{peak1, peak2}{Circular median phases for each condition.}
#'     \item{deltaPhi.Est, deltaPhi.Lower, deltaPhi.Upper}{Phase-difference
#'       median and HDI (filled only when \code{compute_hdi = TRUE}).}
#'     \item{prob_shift, prob_conserved}{Posterior probabilities of shift
#'       and conservation.}
#'     \item{flag_shift, flag_cons, flag_undetermined}{BFDR-significant flags.}
#'     \item{BFDR_shift, BFDR_cons}{Running BFDR values.}
#'   }
#'
#' @export
phase_infer <- function(phi_matrix1, phi_matrix2, gain_loss_status,
                        P = 24, credMass = 0.95, shift = 2,
                        a = -P/2, bfdr_alpha = 0.05,
                        compute_hdi = FALSE) {
  n_genes <- nrow(phi_matrix1)
  gene_names <- rownames(phi_matrix1)
  
  # helper functions -----------------------------------------------------------
  circ_med <- function(x, P) {
    r <- circular((x / P) * 2 * pi, units = "radians")
    t <- (median.circular(r, type = "median") / (2 * pi)) * P
    ((t + P / 2) %% P) - P / 2
  }
  circ_med_24 <- function(x, P = 24) {
    r <- circular((x / P) * 2 * pi, units = "radians")
    t <- (median.circular(r, type = "median") / (2 * pi)) * P
    t %% P
  }
  
  # phase medians --------------------------------------------------------------
  peak1 <- apply(phi_matrix1, 1, circ_med_24, P)
  peak2 <- apply(phi_matrix2, 1, circ_med_24, P)
  phase_diff_matrix <- ((phi_matrix1 - phi_matrix2 + P / 2) %% P) - P / 2
  
  maintained_idx <- which(gain_loss_status == "Maintained")
  deltaPhi.Est <- deltaPhi.Lower <- deltaPhi.Upper <-
    hdi.width <- prob_shift <- prob_conserved <- in_HDI <- rep(NA_real_, n_genes)
  
  # progress bar ---------------------------------------------------------------
  if (length(maintained_idx) > 0) {
    message("Computing phase metrics for maintained genes...")

    for (k in seq_along(maintained_idx)) {
      i <- maintained_idx[k]
      vec <- phase_diff_matrix[i, ]
      
      # Always compute probabilities (needed for BFDR)
      prob_shift[i] <- mean(abs(vec) >= shift)
      prob_conserved[i] <- 1 - prob_shift[i]
      
      # Only compute HDI if requested
      if (compute_hdi) {
        hdi <- circular_HDI(vec, credMass, P, a)
        phi_est <- circ_med(vec, P)
        
        in_HDI[i] <- in_circular_interval(hdi$lower, hdi$upper, 0, P, a)
        deltaPhi.Est[i] <- phi_est
        deltaPhi.Lower[i] <- hdi$lower
        deltaPhi.Upper[i] <- hdi$upper
        
        hdi.width[i] <- if (hdi$upper >= hdi$lower)
          hdi$upper - hdi$lower else (P / 2 - hdi$lower) + (hdi$upper + P / 2)
      }
      
          }
  }
  
  # BFDR control for phase-shift / conserved -----------------------------------
  flag_shift <- rep(FALSE, n_genes)
  flag_cons  <- rep(FALSE, n_genes)
  flag_undetermined <- rep(FALSE, n_genes)
  BFDR_shift_vec <- rep(NA_real_, n_genes)
  BFDR_cons_vec  <- rep(NA_real_, n_genes)
  
  if (length(maintained_idx) > 0) {
    ps  <- prob_shift[maintained_idx]
    pcs <- prob_conserved[maintained_idx]
    
    # ----- SHIFT -----
    ord_shift <- order(-ps)
    BFDR_shift_sorted <- cumsum(1 - ps[ord_shift]) / seq_along(ps[ord_shift])
    sel_shift <- (BFDR_shift_sorted <= bfdr_alpha)
    
    BFDR_shift_vec[maintained_idx[ord_shift]] <- BFDR_shift_sorted
    idx_shift_local <- ord_shift[sel_shift]
    
    # ----- CONSERVED -----
    ord_cons <- order(-pcs)
    BFDR_cons_sorted <- cumsum(1 - pcs[ord_cons]) / seq_along(pcs[ord_cons])
    sel_cons <- (BFDR_cons_sorted <= bfdr_alpha)
    
    BFDR_cons_vec[maintained_idx[ord_cons]] <- BFDR_cons_sorted
    idx_cons_local <- ord_cons[sel_cons]
    
    # ----- EXCLUSIVE CATEGORIZATION -----
    overlap_local <- intersect(idx_shift_local, idx_cons_local)
    
    if (length(overlap_local) > 0) {
      idx_shift_local <- setdiff(idx_shift_local, overlap_local)
      idx_cons_local <- setdiff(idx_cons_local, overlap_local)
    }
    
    idx_union_local <- union(idx_shift_local, idx_cons_local)
    idx_undetermined_local <- setdiff(seq_along(maintained_idx), idx_union_local)
    
    flag_shift[maintained_idx[idx_shift_local]] <- TRUE
    flag_cons[maintained_idx[idx_cons_local]] <- TRUE
    flag_undetermined[maintained_idx[idx_undetermined_local]] <- TRUE
  }
  
  # summary print --------------------------------------------------------------
  cat("\n=== PHASE INFERENCE SUMMARY ===\n")
  cat("Maintained genes:", length(maintained_idx), "\n")
  cat("BFDR α =", bfdr_alpha, " | shift threshold =", shift, "h\n")
  cat("Significant phase-shifted genes:", sum(flag_shift, na.rm = TRUE), "\n")
  cat("Significant phase-conserved genes:", sum(flag_cons, na.rm = TRUE), "\n")
  cat("Undetermined genes:", sum(flag_undetermined, na.rm = TRUE), "\n")
  if (compute_hdi) {
    cat("HDI computation: enabled\n\n")
  } else {
    cat("HDI computation: skipped (use compute_hdi=TRUE to enable)\n\n")
  }
  
  # output ---------------------------------------------------------------------
  # Name every per-gene vector by gene_names so downstream code can index by
  # gene symbol (matching transition_classify()'s gain_loss_status, which is
  # already named) instead of silently returning NA on character indexing.
  names(peak1) <- names(peak2) <- gene_names
  names(deltaPhi.Est) <- names(deltaPhi.Lower) <- names(deltaPhi.Upper) <- gene_names
  names(prob_shift) <- names(flag_shift) <- names(BFDR_shift_vec) <- gene_names
  names(prob_conserved) <- names(flag_cons) <- names(BFDR_cons_vec) <- gene_names
  names(flag_undetermined) <- gene_names
  names(hdi.width) <- names(in_HDI) <- gene_names

  list(
    peak1 = peak1,
    peak2 = peak2,
    deltaPhi.Est = deltaPhi.Est,
    deltaPhi.Lower = deltaPhi.Lower,
    deltaPhi.Upper = deltaPhi.Upper,
    prob_shift = prob_shift,
    flag_shift = flag_shift,
    BFDR_shift = BFDR_shift_vec,
    prob_conserved = prob_conserved,
    flag_cons = flag_cons,
    BFDR_cons  = BFDR_cons_vec,
    flag_undetermined = flag_undetermined,
    hdi.width = hdi.width,
    in_HDI = in_HDI
  )
}



#' Full BayRC phase analysis pipeline
#'
#' @title Integrated rhythmicity detection, transition classification, and phase inference
#'
#' @description
#' Convenience wrapper that calls \code{detect_rhy}, \code{transition_classify},
#' and \code{phase_infer} in sequence.  Returns a consolidated list covering
#' rhythmic gene detection, gain/loss/maintained classification, and
#' phase-shift / conservation flags.
#'
#' @param matrix1 Named list; MCMC output for condition 1 (with \code{rho}
#'   and \code{phi} matrices).
#' @param matrix2 Named list; MCMC output for condition 2.
#' @param P Numeric; period in hours (default 24).
#' @param credMass Numeric; HDI coverage (default 0.95).
#' @param shift Numeric; phase-shift threshold in hours (default 4).
#' @param a Numeric; left endpoint of the circular interval (default
#'   \code{-P/2}).
#' @param bfdr_alpha Numeric; BFDR control level (default 0.05).
#' @param compute_hdi Logical; passed to \code{phase_infer} (default
#'   \code{FALSE}).
#'
#' @return A list with elements \code{rhythmic_summary} (from
#'   \code{detect_rhy}), \code{transition_summary} (from
#'   \code{transition_classify}), and \code{phase_summary} (from
#'   \code{phase_infer}).
#'
phase_analysis <- function(matrix1, matrix2,
                           P = 24, credMass = 0.95,
                           shift = 4, a = -P/2,
                           bfdr_alpha = 0.05,
                           compute_hdi = FALSE) {  # NEW
  rhy <- detect_rhy(matrix1, matrix2, bfdr_alpha)
  pA <- rowMeans(matrix1$rho)
  pB <- rowMeans(matrix2$rho)
  trans <- transition_classify(pA, pB, bfdr_alpha)
  
  gene_names <- rownames(matrix1$phi)
  gain_loss_status <- rep("Non-rhythmic", length(gene_names))
  names(gain_loss_status) <- gene_names
  gain_loss_status[trans$gain_genes] <- "Gain"
  gain_loss_status[trans$loss_genes] <- "Loss"
  gain_loss_status[trans$cons_genes] <- "Maintained"
  
  phase <- phase_infer(matrix1$phi, matrix2$phi, gain_loss_status,
                       P, credMass, shift, a, bfdr_alpha, compute_hdi)  # PASS IT
  
  list(
    rhythmic_summary   = rhy,
    transition_summary = trans,
    phase_summary      = phase
  )
}




# phase_analysis <- function(matrix1, matrix2,
#                            P = 24,
#                            credMass = 0.95,
#                            shift = 4,
#                            a = -P/2,
#                            bfdr_alpha = 0.05) {
#   require(circular)
#   
#   # ------------------------------------------------------------------
#   # basic objects
#   # ------------------------------------------------------------------
#   phi_matrix1 <- matrix1$phi
#   phi_matrix2 <- matrix2$phi
#   gene_names  <- rownames(phi_matrix1)
#   n_genes     <- length(gene_names)
#   
#   # helpers for circular medians
#   circ_med <- function(x, P) {
#     r <- circular((x / P) * 2 * pi, units = "radians")
#     t <- (median.circular(r, type = "median") / (2 * pi)) * P
#     ((t + P/2) %% P) - P/2
#   }
#   circ_med_24 <- function(x, P = 24) {
#     r <- circular((x / P) * 2 * pi, units = "radians")
#     t <- (median.circular(r, type = "median") / (2 * pi)) * P
#     t %% P
#   }
#   
#   # ------------------------------------------------------------------
#   # (1) rhythmic BFDR (per condition)
#   # ------------------------------------------------------------------
#   rhy <- detect_rhy(matrix1, matrix2, bfdr_alpha = bfdr_alpha)
#   pA  <- rowMeans(matrix1$rho)
#   pB  <- rowMeans(matrix2$rho)
#   
#   # ------------------------------------------------------------------
#   # (2) transition probabilities + BFDR on gain/loss/cons
#   # ------------------------------------------------------------------
#   p_gain <- (1 - pA) * pB           # A=0, B=1
#   p_loss <- pA * (1 - pB)           # A=1, B=0
#   p_cons <- pA * pB                 # A=1, B=1
#   
#   bfdr_gain <- bfdr_from_posterior(p_gain, alpha = bfdr_alpha)
#   bfdr_loss <- bfdr_from_posterior(p_loss, alpha = bfdr_alpha)
#   bfdr_cons <- bfdr_from_posterior(p_cons, alpha = bfdr_alpha)
#   
#   tau_gain <- bfdr_gain$threshold
#   tau_loss <- bfdr_loss$threshold
#   tau_cons <- bfdr_cons$threshold
#   
#   gain_genes <- which(p_gain >= tau_gain)
#   loss_genes <- which(p_loss >= tau_loss)
#   cons_genes <- which(p_cons >= tau_cons)
#   
#   n_gain <- length(gain_genes)
#   n_loss <- length(loss_genes)
#   n_cons <- length(cons_genes)
#   
#   cat("=== Transition-level BFDR results ===\n")
#   cat("τ_gain =", round(tau_gain, 3), "| n_gain =", n_gain, "\n")
#   cat("τ_loss =", round(tau_loss, 3), "| n_loss =", n_loss, "\n")
#   cat("τ_cons =", round(tau_cons, 3), "| n_cons =", n_cons, "\n\n")
#   
#   # posterior-based classification (replaces Venn-style classification)
#   gain_loss_status <- rep("Non-rhythmic", n_genes)
#   names(gain_loss_status) <- gene_names
#   gain_loss_status[gain_genes] <- "Gain"
#   gain_loss_status[loss_genes] <- "Loss"
#   gain_loss_status[cons_genes] <- "Maintained"
#   
#   # ------------------------------------------------------------------
#   # (3) phase estimates for all genes
#   # ------------------------------------------------------------------
#   peak1 <- apply(phi_matrix1, 1, circ_med_24, P)
#   peak2 <- apply(phi_matrix2, 1, circ_med_24, P)
#   
#   phase_diff_matrix <- ((phi_matrix1 - phi_matrix2 + P/2) %% P) - P/2
#   
#   phi1.Var   <- numeric(n_genes)
#   phi1.Lower <- numeric(n_genes)
#   phi1.Upper <- numeric(n_genes)
#   phi2.Var   <- numeric(n_genes)
#   phi2.Lower <- numeric(n_genes)
#   phi2.Upper <- numeric(n_genes)
#   
#   cat("Computing circular variance and HDI...\n")
#   for (i in seq_len(n_genes)) {
#     phi1 <- phi_matrix1[i, ]
#     phi2 <- phi_matrix2[i, ]
#     
#     r1 <- circular((phi1 / P) * 2 * pi, units = "radians")
#     r2 <- circular((phi2 / P) * 2 * pi, units = "radians")
#     
#     phi1.Var[i] <- var.circular(r1)
#     phi2.Var[i] <- var.circular(r2)
#     
#     h1 <- circular_HDI(phi1, credMass = credMass, P = P, a = 0)
#     h2 <- circular_HDI(phi2, credMass = credMass, P = P, a = 0)
#     
#     phi1.Lower[i] <- h1$lower
#     phi1.Upper[i] <- h1$upper
#     phi2.Lower[i] <- h2$lower
#     phi2.Upper[i] <- h2$upper
#   }
#   
#   # ------------------------------------------------------------------
#   # (4) phase difference analysis (only for genes rhythmic in BOTH)
#   # ------------------------------------------------------------------
#   maintained_idx <- which(rhy$rhythmic_both_logical)
#   
#   deltaPhi.Est     <- rep(NA_real_, n_genes)
#   deltaPhi.Lower   <- rep(NA_real_, n_genes)
#   deltaPhi.Upper   <- rep(NA_real_, n_genes)
#   hdi.width        <- rep(NA_real_, n_genes)
#   prob_shift       <- rep(NA_real_, n_genes)
#   prob_conserved   <- rep(NA_real_, n_genes)
#   in_HDI           <- rep(NA,        n_genes)
#   
#   if (length(maintained_idx) > 0) {
#     for (i in maintained_idx) {
#       vec <- phase_diff_matrix[i, ]
#       
#       # HDI on phase difference
#       hdi <- circular_HDI(vec, credMass = credMass, P = P, a = a)
#       
#       # circular median of phase difference
#       phi_est <- circ_med(vec, P)
#       
#       # posterior prob of shift ≥ user-specified "shift"
#       prob_shift[i]     <- mean(abs(vec) >= shift)
#       prob_conserved[i] <- 1 - prob_shift[i]
#       
#       in_HDI[i] <- in_circular_interval(hdi$lower, hdi$upper, 0, P, a)
#       
#       deltaPhi.Est[i]   <- phi_est
#       deltaPhi.Lower[i] <- hdi$lower
#       deltaPhi.Upper[i] <- hdi$upper
#       
#       # wrap-length of HDI
#       hdi.width[i] <- if (hdi$upper >= hdi$lower) {
#         hdi$upper - hdi$lower
#       } else {
#         (P/2 - hdi$lower) + (hdi$upper + P/2)
#       }
#     }
#   }
#   
#   # ------------------------------------------------------------------
#   # (5) BFDR on phase shift / phase conserved (only among maintained)
#   # ------------------------------------------------------------------
#   flag_shift <- rep(FALSE, n_genes)
#   flag_cons  <- rep(FALSE, n_genes)
#   
#   if (length(maintained_idx) > 0) {
#     # work on maintained only to avoid NA pollution
#     ps  <- prob_shift[maintained_idx]
#     pcs <- prob_conserved[maintained_idx]
#     inH <- in_HDI[maintained_idx]
#     
#     # order by decreasing prob_shift
#     ord_shift <- order(-ps)
#     BFDR_shift <- cumsum(1 - ps[ord_shift]) / seq_along(ps[ord_shift])
#     sel_shift <- (!inH[ord_shift] & ps[ord_shift] > 0.5 &
#                     BFDR_shift <= bfdr_alpha)
#     flag_shift[maintained_idx[ord_shift[sel_shift]]] <- TRUE
#     
#     # order by decreasing prob_conserved
#     ord_cons <- order(-pcs)
#     BFDR_cons <- cumsum(1 - pcs[ord_cons]) / seq_along(pcs[ord_cons])
#     sel_cons <- (inH[ord_cons] & pcs[ord_cons] > 0.5 &
#                    BFDR_cons <= bfdr_alpha)
#     flag_cons[maintained_idx[ord_cons[sel_cons]]] <- TRUE
#   }
#   
#   # ------------------------------------------------------------------
#   # (6) print summary and return
#   # ------------------------------------------------------------------
#   cat("\n=== RHYTHMICITY SUMMARY ===\n")
#   cat("BFDR α =", bfdr_alpha, "\n")
#   cat("τ_A =", round(rhy$threshold_A, 3),
#       "τ_B =", round(rhy$threshold_B, 3), "\n")
#   cat("Rhythmic (A):", rhy$n_rhythmic_A,
#       "Rhythmic (B):", rhy$n_rhythmic_B,
#       "Both:", rhy$n_rhythmic_both, "\n")
#   
#   cat("\n=== PHASE ANALYSIS (maintained genes) ===\n")
#   cat("Shifted genes:", sum(flag_shift, na.rm = TRUE),
#       "Conserved:", sum(flag_cons, na.rm = TRUE), "\n")
#   
#   return(list(
#     all_genes = data.frame(
#       gene         = gene_names,
#       phi1.Est     = peak1,
#       phi1.Var     = phi1.Var,
#       phi1.Lower   = phi1.Lower,
#       phi1.Upper   = phi1.Upper,
#       phi2.Est     = peak2,
#       phi2.Var     = phi2.Var,
#       phi2.Lower   = phi2.Lower,
#       phi2.Upper   = phi2.Upper,
#       deltaPhi.Est = deltaPhi.Est,
#       deltaPhi.Lower = deltaPhi.Lower,
#       deltaPhi.Upper = deltaPhi.Upper,
#       hdi.width    = hdi.width,
#       prob_shift   = prob_shift,
#       prob_conserved = prob_conserved,
#       in_HDI       = in_HDI,
#       flag_shift   = flag_shift,
#       flag_conserved = flag_cons,
#       posterior_A  = pA,
#       posterior_B  = pB,
#       passed_BFDR_A = rhy$rhythmic_A_logical,
#       passed_BFDR_B = rhy$rhythmic_B_logical,
#       gain_loss_status = gain_loss_status,
#       stringsAsFactors = FALSE
#     ),
#     rhythmic_summary = rhy,
#     parameters = list(
#       P          = P,
#       credMass   = credMass,
#       shift      = shift,
#       bfdr_alpha = bfdr_alpha,
#       tau_A      = rhy$threshold_A,
#       tau_B      = rhy$threshold_B,
#       tau_gain   = tau_gain,
#       tau_loss   = tau_loss,
#       tau_cons   = tau_cons,
#       n_gain     = n_gain,
#       n_loss     = n_loss,
#       n_cons     = n_cons
#     ),
#     gain_genes = sort(gene_names[gain_genes]),
#     loss_genes = sort(gene_names[loss_genes]),
#     cons_genes = sort(gene_names[cons_genes])
#   ))
# }
