phase_fixture <- function() {
  d <- expand.grid(gene=paste0('g',1:4),tissue=paste0('t',1:4),stringsAsFactors=FALSE)
  d$peak <- rep(c(23,1,10,12),4)+rep(c(0,.2,.4,3),each=4)
  d$resultant <- .9; d$called <- TRUE; d$posterior <- .95
  d
}
test_that('circular distances are invariant to time origin and gene input order', {
  d <- phase_fixture()
  run <- function(d) pathway_phase_hierarchy(d,paste0('g',1:4),paste0('t',1:4),k=2)
  x <- run(d); d$peak <- (d$peak+7)%%24; y <- run(d)
  expect_equal(x$gene_distance,y$gene_distance,tolerance=1e-12)
  expect_equal(x$tissue_distance,y$tissue_distance,tolerance=1e-12)
  z <- run(d[nrow(d):1,])
  expect_equal(y$gene_distance,z$gene_distance)
  expect_equal(x$gene_distance[1,2],1-.9^2*cos(2*pi*2/24))
})
test_that('diffuse phase posteriors give neutral loss and sparse pairs stop', {
  d <- phase_fixture(); d$resultant <- 0
  x <- pathway_phase_hierarchy(d,paste0('g',1:4),paste0('t',1:4),k=2)
  expect_equal(unname(x$gene_distance[1,2]),1)
  d$called[d$tissue=='t4'] <- FALSE
  expect_error(pathway_phase_hierarchy(d,paste0('g',1:4),paste0('t',1:4),k=2),'Insufficient tissue overlap')
  expect_error(pathway_phase_hierarchy(rbind(d,d[1,]),paste0('g',1:4),paste0('t',1:4),k=2),'Duplicate')
})
test_that('rhythmic tissue membership is learned without reference to phase', {
  d <- phase_fixture()
  m <- matrix(.1,4,4,dimnames=list(paste0('t',1:4),paste0('t',1:4)))
  diag(m)<-1; m[1,2]<-m[2,1]<-.9; m[3,4]<-m[4,3]<-.7
  run <- function(d) plot_pathway_profiles(d,paste0('g',1:4),paste0('t',1:4),m,gene_k=2,draw=FALSE)
  x <- run(d); expect_setequal(x$selected_tissues,c('t1','t2'))
  d$peak <- rev(d$peak)
  y <- run(d); expect_equal(x$tissue_groups,y$tissue_groups)
  expect_equal(x$selected_tissues,y$selected_tissues)
})
test_that('insufficient phase coverage preserves rhythmicity profiles', {
  d <- phase_fixture(); d$called[d$gene!='g1'] <- FALSE
  m <- matrix(.1,4,4,dimnames=list(paste0('t',1:4),paste0('t',1:4)))
  diag(m)<-1;m[1,2]<-m[2,1]<-.9;m[3,4]<-m[4,3]<-.7
  expect_warning(x <- plot_pathway_profiles(d,paste0('g',1:4),paste0('t',1:4),m,draw=FALSE),'Phase clustering unavailable')
  expect_null(x$hierarchy)
  expect_match(x$phase_status,'Too few coverage-eligible')
  expect_equal(unname(x$posterior),matrix(.95,4,4))
  expect_length(x$gene_order,4)
})
test_that('sparse comparisons shrink toward neutral rather than stopping', {
  d <- phase_fixture()
  d$called <- (d$gene %in% c('g1','g2') & d$tissue %in% c('t1','t2')) |
    (d$gene %in% c('g3','g4') & d$tissue %in% c('t3','t4'))
  run <- function(mode) pathway_phase_hierarchy(d,paste0('g',1:4),paste0('t',1:4),
    min_gene_tissues=2,k=2,sparse_method=mode,shrink_strength=2)
  expect_error(run('strict'),'Insufficient gene overlap')
  x <- run('shrink')
  expect_equal(unname(x$gene_distance['g1','g3']),1)
  expect_equal(unname(x$gene_distance['g1','g2']),1+.5*((1-.9^2*cos(2*pi*2/24))-1))
  expect_equal(unname(x$gene_overlap['g1','g3',1]),0)
  expect_true(all(is.finite(x$tissue_distance)))
})
