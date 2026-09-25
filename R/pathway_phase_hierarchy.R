#' Uncertainty-aware gene and tissue clustering within a pathway
#'
#' @description First clusters genes by circular peak-time profiles, then
#' clusters tissues using distances averaged equally across the resulting gene
#' groups. This is exploratory clustering, not a test of tissue outlier status.
#' @param data Data frame with unique gene/tissue rows and columns gene, tissue,
#'   peak (hours), resultant (posterior mean resultant length, between 0 and 1),
#'   and called (logical rhythmicity call). Optional interval_width is the
#'   circular credible-interval arc length in hours, used only for plotting.
#' @param pathway_genes Character vector of pathway gene identifiers.
#' @param tissues Character vector of tissues to cluster and display.
#' @param pathway_name Label for the selected pathway.
#' @param gene_tissues Tissues used to learn gene groups; defaults to tissues.
#' @param tissue_blocks Optional named vector assigning each gene_tissues entry
#'   to a block. Gene distances give each block equal weight.
#' @param min_gene_tissues Minimum called tissues per retained gene, either a
#'   single count across gene_tissues or a named count per tissue block.
#' @param min_gene_overlap Minimum shared tissues for each gene pair within
#'   each block. Insufficient overlap stops clustering rather than imputing.
#' @param min_tissue_overlap Minimum shared genes for each tissue pair.
#' @param min_groups Minimum represented gene groups for each tissue pair.
#' @param k Number of gene groups; fixed explicitly, not optimized for a tissue.
#' @param period Cycle length in hours.
#' @param sparse_method Either strict (original overlap checks) or shrink
#'   (exploratory regularization of sparse distances toward neutral loss 1).
#' @param shrink_strength Positive regularization constant. In shrink mode,
#'   a comparison supported by n observations receives weight n/(n+strength).
#' @details For two phase summaries, the dissimilarity is
#'   \eqn{1-R_1 R_2 \cos(2\pi(\mu_1-\mu_2)/period)}. It is the expected circular
#'   loss under independent marginal posteriors. Diffuse posteriors move the
#'   distance toward the neutral value 1; this is not inverse-variance weighting.
#'   Both trees use average linkage. Tissue distances first average over shared
#'   called genes within each gene group, then equally across represented groups.
#'   Coverage requirements are necessary even when uncertainty is incorporated.
#'   Interval widths alone do not identify resultant lengths without additional
#'   distributional assumptions and are not converted to clustering weights.
#'   In exploratory shrink mode, each within-block mean loss d is replaced by
#'   1 + n/(n+shrink_strength)*(d-1), with neutral loss 1 for zero overlap.
#'   The same rule applies within gene groups when comparing tissues. This is
#'   an explicit regularization assumption, not evidence for missing phases.
#'   The requested number of groups is capped at the eligible gene count;
#'   fewer than two eligible genes or entirely neutral gene distances cannot
#'   support a tree. Strict remains the default.
#' @return List containing peak/resultant matrices, gene and tissue hclust
#'   objects, gene-group assignments, distances, overlap counts, excluded genes,
#'   and settings. Non-called cells remain missing.
#' @export
pathway_phase_hierarchy <- function(data, pathway_genes, tissues,
    pathway_name = "Selected pathway", gene_tissues = tissues,
    tissue_blocks = NULL, min_gene_tissues = 3L,
    min_gene_overlap = 2L, min_tissue_overlap = 3L,
    min_groups = 1L, k = 3L, period = 24, sparse_method = "strict", shrink_strength = 2) {
  sparse_method <- match.arg(sparse_method,c("strict","shrink"))
  if(length(shrink_strength)!=1 || !is.finite(shrink_strength) || shrink_strength<=0) stop("shrink_strength must be positive.")
  required <- c("gene", "tissue", "peak", "resultant", "called")
  if (!all(required %in% names(data))) stop("Missing columns: ", paste(setdiff(required, names(data)), collapse = ", "))
  if (!is.logical(data$called) || anyNA(data$called)) stop("called must be logical without missing values.")
  if (!is.numeric(data$peak) || !is.numeric(data$resultant)) stop("peak and resultant must be numeric.")
  integer_scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x)
  if (!all(vapply(list(min_gene_overlap, min_tissue_overlap, min_groups, k), integer_scalar, logical(1)))) stop("Overlap thresholds and k must be positive integers.")
  if (length(period) != 1 || !is.finite(period) || period <= 0) stop("period must be positive.")
  tissues <- unique(as.character(tissues)); gene_tissues <- unique(as.character(gene_tissues))
  all_tissues <- union(gene_tissues, tissues)
  if (length(tissues) < 2 || length(gene_tissues) < 2 || !all(all_tissues %in% data$tissue)) stop("Supply at least two observed tissues for each stage.")
  genes <- intersect(unique(as.character(pathway_genes)), as.character(data$gene))
  d <- data[data$gene %in% genes & data$tissue %in% all_tissues, , drop = FALSE]
  if (anyDuplicated(d[c("gene", "tissue")])) stop("Duplicate gene/tissue rows.")
  if (any(d$called & (!is.finite(d$peak) | !is.finite(d$resultant) | d$resultant < 0 | d$resultant > 1))) stop("Called cells require finite peaks and resultant lengths in [0,1].")
  P <- R <- W <- matrix(NA_real_, length(genes), length(all_tissues), dimnames = list(genes, all_tissues))
  ix <- cbind(match(d$gene[d$called], genes), match(d$tissue[d$called], all_tissues))
  P[ix] <- d$peak[d$called] %% period; R[ix] <- d$resultant[d$called]
  if ("interval_width" %in% names(d)) {
    w <- d$interval_width[d$called]
    if (!is.numeric(w) || any(!is.na(w) & (!is.finite(w) | w < 0 | w > period))) stop("interval_width must be an arc length in [0,period] or NA.")
    W[ix] <- w
  }
  if (is.null(tissue_blocks)) tissue_blocks <- setNames(rep("all", length(gene_tissues)), gene_tissues)
  if (is.null(names(tissue_blocks)) || anyDuplicated(names(tissue_blocks)) || !all(gene_tissues %in% names(tissue_blocks)) || anyNA(tissue_blocks[gene_tissues])) stop("tissue_blocks must name each gene_tissues entry once.")
  blocks <- split(gene_tissues, tissue_blocks[gene_tissues])
  if (!is.numeric(min_gene_tissues) || !length(min_gene_tissues) || any(!is.finite(min_gene_tissues) | min_gene_tissues < 1 | min_gene_tissues != floor(min_gene_tissues))) stop("min_gene_tissues must contain positive integer counts.")
  if (is.null(names(min_gene_tissues)) && length(min_gene_tissues) == 1) {
    keep <- rowSums(is.finite(P[, gene_tissues, drop = FALSE])) >= min_gene_tissues
  } else {
    if (is.null(names(min_gene_tissues)) || !all(names(blocks) %in% names(min_gene_tissues))) stop("Provide a minimum count for every tissue block.")
    keep <- Reduce(`&`, lapply(names(blocks), function(b) rowSums(is.finite(P[, blocks[[b]], drop = FALSE])) >= min_gene_tissues[[b]]))
  }
  excluded <- genes[!keep]; genes <- genes[keep]
  requested_k <- k
  if (sparse_method=="shrink") k <- min(k,length(genes))
  if (length(genes) < max(2, k)) stop("Too few coverage-eligible genes for k groups.")
  P <- P[genes, , drop = FALSE]; R <- R[genes, , drop = FALSE]; W <- W[genes, , drop = FALSE]
  loss <- function(p, q, r, s) 1 - r * s * cos(2 * pi * (p - q) / period)
  G <- matrix(0, length(genes), length(genes), dimnames = list(genes, genes))
  GO <- array(0L, c(length(genes), length(genes), length(blocks)), dimnames = list(genes, genes, names(blocks)))
  for (i in seq_len(length(genes) - 1L)) for (j in seq.int(i + 1L, length(genes))) {
    v <- vapply(seq_along(blocks), function(b) {
      t <- blocks[[b]]; t <- t[is.finite(P[i,t]) & is.finite(P[j,t])]
      GO[i,j,b] <<- GO[j,i,b] <<- length(t)
      if (sparse_method=="strict" && length(t) < min_gene_overlap) stop("Insufficient gene overlap: ", genes[i], " / ", genes[j], " in block ", names(blocks)[b], ". Revise coverage or subset explicitly.")
      if (!length(t)) return(1)
      raw <- mean(loss(P[i,t], P[j,t], R[i,t], R[j,t]))
      if(sparse_method=="shrink") 1 + length(t)/(length(t)+shrink_strength)*(raw-1) else raw
    }, numeric(1))
    G[i,j] <- G[j,i] <- mean(v)
  }
  if(sparse_method=="shrink" && all(abs(G[upper.tri(G)]-1)<1e-12)) stop("Insufficient phase information: all gene distances are neutral.")
  gh <- stats::hclust(stats::as.dist(G), method = "average")
  groups <- stats::cutree(gh, k = k)
  groups[] <- match(groups,unique(groups[gh$order]))
  T <- O <- C <- matrix(0, length(tissues), length(tissues), dimnames = list(tissues, tissues))
  for (i in seq_len(length(tissues) - 1L)) for (j in seq.int(i + 1L, length(tissues))) {
    a <- tissues[i]; b <- tissues[j]
    shared <- genes[is.finite(P[,a]) & is.finite(P[,b])]
    O[i,j] <- O[j,i] <- length(shared)
    v <- vapply(seq_len(k), function(g) {
      s <- shared[groups[shared] == g]
      if (!length(s)) return(if(sparse_method=="shrink")1 else NA_real_)
      raw <- mean(loss(P[s,a], P[s,b], R[s,a], R[s,b]))
      if(sparse_method=="shrink") 1 + length(s)/(length(s)+shrink_strength)*(raw-1) else raw
    }, numeric(1))
    C[i,j] <- C[j,i] <- length(unique(groups[shared]))
    if (sparse_method=="strict" && (length(shared) < min_tissue_overlap || C[i,j] < min_groups)) stop("Insufficient tissue overlap: ", a, " / ", b, ". Revise coverage or subset explicitly.")
    T[i,j] <- T[j,i] <- mean(v, na.rm = TRUE)
  }
  diag(O) <- colSums(is.finite(P[,tissues,drop=FALSE]))
  list(peak = P, resultant = R, interval_width = W,
       gene_tree = gh, tissue_tree = stats::hclust(stats::as.dist(T), "average"),
       gene_groups = groups, gene_distance = G, tissue_distance = T,
       gene_overlap = GO, tissue_overlap = O, group_overlap = C,
       excluded_genes = excluded, tissues = tissues, gene_tissues = gene_tissues,
       pathway_name = pathway_name, period = period,
       settings = list(k=k,requested_k=requested_k,sparse_method=sparse_method,shrink_strength=shrink_strength,tissue_blocks=tissue_blocks, min_gene_tissues=min_gene_tissues,
                       min_gene_overlap=min_gene_overlap, min_tissue_overlap=min_tissue_overlap,
                       min_groups=min_groups))
}

#' Plot gene and tissue trees in one pathway phase heatmap
#' @param x Result of \code{pathway_phase_hierarchy}.
#' @param file Optional PDF output path. Otherwise draws on the current device.
#' @param width,height PDF dimensions in inches.
#' @return Invisibly returns x. The left tree groups genes using gene_tissues;
#' the top tree groups the displayed tissues using equally weighted gene groups.
#' White-disc area represents credible-interval width, if supplied; gray cells
#' are not called rhythmic.
#' @export
plot_pathway_phase_hierarchy <- function(x, file = NULL, width = 10, height = 8) {
  if (!requireNamespace("ComplexHeatmap", quietly=TRUE) || !requireNamespace("circlize", quietly=TRUE)) stop("Install ComplexHeatmap and circlize to plot the hierarchy.")
  if (!is.null(file)) { grDevices::pdf(file, width=width, height=height); on.exit(grDevices::dev.off(), add=TRUE) }
  col <- circlize::colorRamp2(seq(0, x$period, length.out=13), grDevices::hcl(h=seq(15,375,length.out=13), c=68,l=66))
  p <- x$peak[,x$tissues,drop=FALSE]; w <- x$interval_width[,x$tissues,drop=FALSE]
  groups <- paste0("G",x$gene_groups[rownames(p)])
  group_levels <- unique(groups[x$gene_tree$order])
  group_colors <- setNames(grDevices::hcl.colors(length(group_levels),"Dark 3"),group_levels)
  ann <- ComplexHeatmap::rowAnnotation(`Phase group`=groups,
     col=list(`Phase group`=group_colors),show_annotation_name=FALSE,
     annotation_legend_param=list(`Phase group`=list(title="Gene phase group")))
  ht <- ComplexHeatmap::Heatmap(p, name="Peak (h)", col=col, na_col="#e8e8e8",
      cluster_rows=x$gene_tree, row_dend_reorder=FALSE,
      cluster_columns=x$tissue_tree, column_dend_reorder=FALSE,
      row_dend_width=grid::unit(25,"mm"),column_dend_height=grid::unit(25,"mm"),
      left_annotation=ann,column_title=paste0(x$pathway_name,"\n"),
      column_title_gp=grid::gpar(fontsize=15),
      row_names_gp=grid::gpar(fontface="italic",fontsize=12),
      column_names_gp=grid::gpar(fontsize=12),
      cell_fun=function(j,i,xx,yy,cw,ch,fill) {
        if (is.finite(w[i,j])) grid::grid.circle(xx,yy,r=grid::unit.pmin(cw,ch)*.4*sqrt(w[i,j]/x$period),
          gp=grid::gpar(fill="white",col="#666666",lwd=.5))
      }, heatmap_legend_param=list(at=seq(0,x$period,length.out=5)))
  symbols <- ComplexHeatmap::Legend(title="Phase profiles",
    labels=c("95% interval width","Not called rhythmic"),
    graphics=list(
      function(x,y,w,h) grid::grid.circle(x,y,r=grid::unit(1.5,"mm"),gp=grid::gpar(fill="white",col="#666666")),
      function(x,y,w,h) grid::grid.rect(x,y,width=grid::unit(3,"mm"),height=grid::unit(3,"mm"),gp=grid::gpar(fill="#e8e8e8",col=NA))))
  ComplexHeatmap::draw(ht,heatmap_legend_side="right",annotation_legend_side="right",
     annotation_legend_list=list(symbols),padding=grid::unit(c(8,16,8,8),"mm"))
  ComplexHeatmap::decorate_column_dend("Peak (h)", {
    grid::grid.text("Tissue phase clustering",x=.5,y=grid::unit(1,"npc")+grid::unit(2,"mm"),
                   just="bottom",gp=grid::gpar(fontsize=11))
  })
  ComplexHeatmap::decorate_row_dend("Peak (h)", {
    grid::grid.text("Gene phase clustering",x=grid::unit(-4,"mm"),y=.5,rot=90,
                   gp=grid::gpar(fontsize=11))
  })
  invisible(x)
}

#' Paired rhythmicity and phase profiles for a selected pathway
#'
#' @param data Phase-summary data as in \code{pathway_phase_hierarchy}, with
#'   an additional posterior column containing rhythmicity probabilities.
#' @param pathway_genes,tissues,pathway_name See \code{pathway_phase_hierarchy}.
#' @param concordance Named symmetric pathway rhythmic-concordance matrix on
#'   [-1,1]. This must come from rhythmicity analysis, not phase clustering.
#' @param tissue_k Number of rhythmic tissue groups (specified in advance).
#' @param tissue_linkage Linkage for rhythmic concordance clustering. ward.D2
#'   reproduces the original manuscript convention applied to 1-concordance;
#'   its variance-minimization interpretation requires Euclidean distances.
#' @param gene_k Number of gene phase groups.
#' @param min_gene_fraction Minimum called fraction per rhythmic tissue block.
#' @param min_gene_overlap,min_tissue_overlap,min_groups See
#'   \code{pathway_phase_hierarchy}.
#' @param period Cycle length in hours.
#' @param sparse_method Either strict (original overlap checks) or shrink
#'   (exploratory regularization of sparse distances toward neutral loss 1).
#' @param shrink_strength Positive regularization constant. In shrink mode,
#'   a comparison supported by n observations receives weight n/(n+strength).
#' @param file Optional PDF path; NULL draws on the current device.
#' @param draw Whether to draw the panels.
#' @param show_limited Whether to display genes excluded by the phase-coverage
#'   requirement in a separate limited-coverage block (default TRUE). This block
#'   is not a learned phase group; hiding it does not change clustering.
#' @param width,height PDF dimensions in inches.
#' @details Learns tissue groups from rhythmic concordance, identifies the group
#'   with highest mean off-diagonal concordance, then learns gene phase groups
#'   across all supplied tissues with equal block weights. Within the selected
#'   rhythmic group, phase clustering orders tissues without changing membership.
#'   Both panels retain all observed pathway genes; genes failing phase coverage
#'   are shown after eligible genes. Gray phase cells are not called rhythmic.
#'   The same gene and tissue order is used in both panels. Selection of the
#'   most concordant group is explicit and descriptive, not a significance test.
#'   If coverage or overlap prevents phase clustering, a warning is issued,
#'   hierarchy is NULL, and phase_status explains why. Both data panels remain
#'   available in rhythmic tissue order, without assigned gene phase groups.
#'   Phase-group and coverage labels appear only beside the phase panel.
#' @return Invisibly returns rhythmic tree, learned memberships and group scores,
#'   selected group, phase hierarchy, plotted matrices and display orders.
#' @export
plot_pathway_profiles <- function(data, pathway_genes, tissues, concordance,
    pathway_name = "Selected pathway", tissue_k = 2L, tissue_linkage = "ward.D2",
    gene_k = 3L, min_gene_fraction = .5, min_gene_overlap = 2L,
    min_tissue_overlap = 3L, min_groups = 1L, period = 24,
    file = NULL, draw = TRUE, width = 16, height = 9, show_limited = TRUE,
    sparse_method = "strict", shrink_strength = 2) {
  if (!is.logical(show_limited) || length(show_limited)!=1 || is.na(show_limited)) stop("show_limited must be TRUE or FALSE.")
  tissues <- unique(as.character(tissues))
  if (!is.matrix(concordance) || !is.numeric(concordance) ||
      is.null(rownames(concordance)) || is.null(colnames(concordance)) ||
      anyDuplicated(rownames(concordance)) || anyDuplicated(colnames(concordance)) ||
      !all(tissues %in% rownames(concordance)) || !all(tissues %in% colnames(concordance))) stop("Provide a named concordance matrix covering all tissues.")
  m <- concordance[tissues,tissues,drop=FALSE]
  if (any(!is.finite(m)) || any(m < -1 | m > 1) || max(abs(m-t(m))) > 1e-8) stop("Concordance must be finite, symmetric and in [-1,1].")
  if (length(tissue_k)!=1 || !is.finite(tissue_k) || tissue_k < 2 || tissue_k >= length(tissues) || tissue_k != floor(tissue_k)) stop("tissue_k must be an integer between 2 and number of tissues minus 1.")
  if (length(min_gene_fraction)!=1 || !is.finite(min_gene_fraction) || min_gene_fraction<=0 || min_gene_fraction>1) stop("min_gene_fraction must be in (0,1].")
  if (!"posterior" %in% names(data) || !is.numeric(data$posterior) || any(!is.na(data$posterior) & (!is.finite(data$posterior) | data$posterior < 0 | data$posterior > 1))) stop("Supply posterior rhythmicity probabilities in [0,1].")
  diag(m) <- 1
  rh <- stats::hclust(stats::as.dist(1-m), method=tissue_linkage)
  membership <- stats::cutree(rh,k=tissue_k)
  score <- vapply(seq_len(tissue_k), function(g) {
    s <- m[membership==g,membership==g,drop=FALSE]
    if(nrow(s)<2) return(NA_real_)
    mean(s[lower.tri(s)])
  }, numeric(1))
  if(all(is.na(score))) stop("No non-singleton rhythmic tissue group.")
  selected <- which.max(replace(score,is.na(score),-Inf))
  chosen <- names(membership)[membership==selected]
  blocks <- setNames(as.character(membership),names(membership))
  minimum <- ceiling(table(blocks)*min_gene_fraction)
  phase_status <- "Clustered"
  h <- tryCatch(pathway_phase_hierarchy(data,pathway_genes,chosen,pathway_name,
       gene_tissues=tissues,tissue_blocks=blocks,min_gene_tissues=minimum,
       min_gene_overlap=min_gene_overlap,min_tissue_overlap=min_tissue_overlap,
       min_groups=min_groups,k=gene_k,period=period,sparse_method=sparse_method,shrink_strength=shrink_strength), error=function(e) {
         msg <- conditionMessage(e)
         if (!grepl("^(Too few coverage-eligible genes|Insufficient gene overlap|Insufficient tissue overlap|Insufficient phase information)",msg)) stop(e)
         phase_status <<- msg
         warning("Phase clustering unavailable: ",msg," Rhythmicity and available phase profiles remain displayable.",call.=FALSE)
         NULL
       })
  genes <- intersect(unique(as.character(pathway_genes)),as.character(data$gene))
  d <- data[data$gene %in% genes & data$tissue %in% tissues,,drop=FALSE]
  if(anyDuplicated(d[c("gene","tissue")])) stop("Duplicate gene/tissue rows.")
  P <- Q <- W <- matrix(NA_real_,length(genes),length(tissues),dimnames=list(genes,tissues))
  ix <- cbind(match(d$gene,genes),match(d$tissue,tissues)); Q[ix] <- d$posterior
  P[ix[d$called,,drop=FALSE]] <- d$peak[d$called] %% period
  if("interval_width" %in% names(d)) W[ix[d$called,,drop=FALSE]] <- d$interval_width[d$called]
  if (is.null(h)) {
    coverage <- lapply(names(minimum),function(b) {
      ts <- names(blocks)[blocks==b]
      rowSums(is.finite(P[,ts,drop=FALSE])) >= minimum[[b]]
    })
    eligible <- genes[Reduce(`&`,coverage)]
    go <- eligible
  } else {
    eligible <- h$gene_tree$labels
    go <- h$gene_tree$labels[h$gene_tree$order]
  }
  if (show_limited) go <- c(go,setdiff(genes,eligible))
  if (!length(go)) stop("No coverage-eligible genes to display; use show_limited=TRUE.")
  to <- rh$labels[rh$order]
  if (!is.null(h)) to[to %in% chosen] <- h$tissue_tree$labels[h$tissue_tree$order]
  result <- list(rhythmic_tree=rh,tissue_groups=membership,group_scores=score,
                 selected_group=selected,selected_tissues=chosen,hierarchy=h,phase_status=phase_status,
                 posterior=Q,peak=P,interval_width=W,gene_order=go,tissue_order=to)
  if(draw) {
    if(!requireNamespace("ComplexHeatmap",quietly=TRUE) || !requireNamespace("circlize",quietly=TRUE)) stop("Install ComplexHeatmap and circlize to plot profiles.")
    if(!is.null(file)) { grDevices::pdf(file,width=width,height=height); on.exit(grDevices::dev.off(),add=TRUE) }
    phase_col <- circlize::colorRamp2(seq(0,period,length.out=13),grDevices::hcl(h=seq(15,375,length.out=13),c=68,l=66))
    # Match the manuscript concordance palette, scaled to probabilities [0,1].
    prob_colors <- grDevices::colorRampPalette(c(
      "#DFF2F4", "#C6E9DC", "#B3E2C6", "#82CEC1",
      "#47AACC", "#3394C2", "#2282B9", "#1372B1"))(200)
    prob_col <- circlize::colorRamp2(seq(0,1,length.out=200),prob_colors)
    group_labels <- paste0("Group ",membership[to]); group_labels[membership[to]==selected] <- paste0("Group ",selected," (highest concordance)")
    split <- factor(group_labels,levels=unique(group_labels))
    row_groups <- if(is.null(h)) {
      ifelse(go %in% eligible,"Coverage eligible\nUnclustered","Below clustering\ncoverage threshold")
    } else ifelse(go %in% names(h$gene_groups),paste0("Phase group ",h$gene_groups[go]),"Below clustering\ncoverage threshold")
    row_split <- factor(row_groups,levels=unique(row_groups))
    grid::grid.newpage(); grid::pushViewport(grid::viewport(layout=grid::grid.layout(1,2),y=.51,height=.86))
    for(i in 1:2) {
      grid::pushViewport(grid::viewport(layout.pos.row=1,layout.pos.col=i))
      mat <- if(i==1) Q[go,to,drop=FALSE] else P[go,to,drop=FALSE]
      ww <- W[go,to,drop=FALSE]
      ht <- ComplexHeatmap::Heatmap(mat,name=if(i==1)"Rhythmicity" else "Peak (h)",
        col=if(i==1)prob_col else phase_col,na_col="#e8e8e8",
        heatmap_legend_param=list(at=if(i==1)c(0,.5,1) else seq(0,period,length.out=5),
          title=if(i==1)expression(atop("Rhythmicity",Pr(rho==1~"|"~data))) else "Peak (h)"),
        cluster_rows=FALSE,cluster_columns=FALSE,row_split=row_split,
        column_split=split,cluster_row_slices=FALSE,cluster_column_slices=FALSE,
        row_gap=grid::unit(5,"mm"),row_title=if(i==1)NULL else levels(row_split),row_title_rot=0,
        row_title_gp=grid::gpar(fontsize=9),column_title_gp=grid::gpar(fontsize=9),
        row_names_gp=grid::gpar(fontsize=9,fontface="italic"),column_names_gp=grid::gpar(fontsize=9),
        cell_fun=if(i==1)NULL else function(j,i,x,y,w,h,fill) {
          if(is.finite(ww[i,j])) grid::grid.circle(x,y,r=grid::unit.pmin(w,h)*.4*sqrt(ww[i,j]/period),gp=grid::gpar(fill="white",col="#666666",lwd=.4))
        })
      ComplexHeatmap::draw(ht,newpage=FALSE,column_title=if(i==1)"C  Rhythmicity profiles" else "D  Phase profiles")
      grid::popViewport()
    }
    grid::popViewport(); grid::grid.text(pathway_name,y=.975,gp=grid::gpar(fontsize=15,fontface="bold"))
  }
  invisible(result)
}
