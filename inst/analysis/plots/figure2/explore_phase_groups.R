# Figure 2 panel D: peak timing of the KEGG circadian pathway across the 25
# tissues, with the genes grouped by phase. Reads pathway_phase_summary.csv from
# export_pathway.R and writes option8_<mode>.pdf, option8_<mode>_nature.pdf and
# option8_<mode>_long.pdf for three tissue bases, plus the clustering tables.
# The paper's panel is option8_balanced_global_nature.pdf.
library(ComplexHeatmap)
library(circlize)
library(cluster)
library(grid)

# Paths come from config.R; override any of them with the matching env var.
bayrc.needs.summary <- FALSE
this.file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
analysis.dir <- if (is.na(this.file)) getwd() else dirname(normalizePath(this.file))
while (!file.exists(file.path(analysis.dir, "config.R")) &&
       dirname(analysis.dir) != analysis.dir) analysis.dir <- dirname(analysis.dir)
source(file.path(analysis.dir, "config.R"))

# the panel and its tables sit beside the other Figure 2 panels
out <- Sys.getenv("FIG2_PHASE_DIR", unset = file.path(BAYRC_FIGURE_DIR, "figure2"))
if (!file.exists(file.path(out, "pathway_phase_summary.csv")))
  stop("no pathway_phase_summary.csv under ", out, "; run plots/figure2/export_pathway.R first")
# cairo embeds its fonts and carries the hyphen the block names are written
# with, so the panel that goes into the paper is drawn with it where it is
# available and with the base device otherwise
fig_device <- if (nzchar(Sys.getenv("FIG2_CAIRO"))) grDevices::cairo_pdf else grDevices::pdf
d <- read.csv(file.path(out,"pathway_phase_summary.csv"))
genes<-unique(d$gene);tissues<-unique(d$tissue)
blocks<-c("High-concordance cluster","Remaining tissues")
ts<-lapply(blocks,function(b)unique(d$tissue[d$block==b]));names(ts)<-blocks
P<-R<-H<-matrix(NA_real_,length(genes),length(tissues),dimnames=list(genes,tissues))
for(i in seq_len(nrow(d))) if(d$called[i]) {
 g<-d$gene[i];t<-d$tissue[i];P[g,t]<-d$peak[i];R[g,t]<-d$resultant[i]
 H[g,t]<-(d$upper[i]-d$lower[i])%%24
}
coverage<-data.frame(gene=genes,high=rowSums(is.finite(P[,ts[[1]]])),
                    remaining=rowSums(is.finite(P[,ts[[2]]])))
# A shared eligible set makes the three clustering approaches comparable.
core<-coverage$gene[coverage$high>=6 & coverage$remaining>=4]
stopifnot(length(core)>=4)
write.csv(coverage,file.path(out,"phase_coverage.csv"),row.names=FALSE)
distances<-function(g,sets,uncertainty=TRUE) {
 D<-matrix(0,length(g),length(g),dimnames=list(g,g))
 for(i in seq_along(g)) for(j in seq_along(g)) if(j>i) {
  vals<-vapply(sets,function(t) {
   ok<-is.finite(P[g[i],t])&is.finite(P[g[j],t]);t<-t[ok]
   if(!length(t)) return(NA_real_)
   r<-if(uncertainty) R[g[i],t]*R[g[j],t] else rep(1,length(t))
   mean(1-r*cos(2*pi*(P[g[i],t]-P[g[j],t])/24))
  },numeric(1))
  # No-overlap bootstrap samples contribute neutral similarity, not opposition.
  vals[!is.finite(vals)]<-1
  D[i,j]<-D[j,i]<-mean(vals)
 }
 as.dist(D)
}
overlap<-do.call(rbind,lapply(blocks,function(b) {
 comb<-combn(core,2)
 data.frame(group=b,gene1=comb[1,],gene2=comb[2,],
 n=vapply(seq_len(ncol(comb)),function(i)sum(is.finite(P[comb[1,i],ts[[b]]])&
                                               is.finite(P[comb[2,i],ts[[b]]])),integer(1)))
}))
write.csv(overlap,file.path(out,"shared_tissue_coverage.csv"),row.names=FALSE)
if(any(overlap$n==0)) stop("Eligible genes have no overlapping tissue coverage; inspect shared_tissue_coverage.csv")
sets<-list(high_only=ts[1],remaining_only=ts[2],balanced_global=ts)
titles<-c(high_only="Peak timing of the circadian pathway",
          remaining_only="Peak timing of the circadian pathway",
          balanced_global="Peak timing of the circadian pathway")
subtitles<-c(
 high_only="Posterior peak time where the gene is called rhythmic; disc area is the width of the 95% credible interval. Genes grouped by phase using the high-concordance tissues",
 remaining_only="Posterior peak time where the gene is called rhythmic; disc area is the width of the 95% credible interval. Genes grouped by phase using the remaining tissues",
 balanced_global="Posterior peak time where the gene is called rhythmic; disc area is the width of the 95% credible interval. Genes grouped by phase, both tissue groups weighted equally")
cyc<-colorRamp2(seq(0,24,length.out=13),hcl(h=seq(15,375,length.out=13),c=68,l=66))
all_stats<-list();assignments<-list();set.seed(927)
for(mode in names(sets)) {
 D<-distances(core,sets[[mode]])
 hc<-hclust(D,method="average")
 ks<-2:min(5,length(core)-1)
 scores<-vapply(ks,function(k)mean(silhouette(cutree(hc,k),D)[,"sil_width"]),numeric(1))
 k<-ks[which.max(scores)];grp<-cutree(hc,k)
 co<-matrix(0,length(core),length(core),dimnames=list(core,core))
 for(b in 1:200) {
  sampled<-lapply(sets[[mode]],function(t)sample(t,length(t),replace=TRUE))
  z<-cutree(hclust(distances(core,sampled),"average"),k)
  co<-co+outer(z,z,"==")/200
 }
 write.csv(co,file.path(out,paste0(mode,"_bootstrap_coassignment.csv")))
 # Compare uncertainty-aware grouping to the same rule on point estimates.
 raw<-cutree(hclust(distances(core,sets[[mode]],FALSE),"average"),k)
 agreement<-mean(outer(raw,raw,"==")[lower.tri(co)]==outer(grp,grp,"==")[lower.tri(co)])
 stability<-vapply(sort(unique(grp)),function(cl) {
  ix<-which(grp==cl);if(length(ix)<2)return(NA_real_)
  x<-co[ix,ix,drop=FALSE];mean(x[lower.tri(x)])
 },numeric(1))
 all_stats[[mode]]<-data.frame(mode,k=ks,silhouette=scores,selected=ks==k,
                              point_estimate_pair_agreement=agreement)
 coreorder<-hc$labels[hc$order]
 sparse<-setdiff(genes,core)
 sparse<-sparse[order(-rowSums(is.finite(P[sparse,,drop=FALSE])))]
 ord<-c(coreorder,sparse)
 groups_order<-unique(grp[coreorder])
 lab<-setNames(paste0("Phase group ",seq_along(groups_order)),groups_order)
 rowgroup<-c(unname(lab[as.character(grp[coreorder])]),rep("Limited coverage",length(sparse)))
 rowgroup<-factor(rowgroup,levels=c(unname(lab),"Limited coverage"))
 assignments[[mode]]<-data.frame(mode,gene=ord,phase_group=as.character(rowgroup),
                                stability=c(stability[grp[coreorder]],rep(NA_real_,length(sparse))))
 mat<-P[ord,,drop=FALSE];wid<-H[ord,,drop=FALSE]
 reference_r<-2
 cell<-function(j,i,x,y,w,h,fill) {
  grid.rect(x,y,w,h,gp=gpar(fill=if(is.finite(mat[i,j]))cyc(mat[i,j]) else "#E8E8E8",col="white",lwd=.45))
  if(is.finite(wid[i,j])) {
   maxr<-min(convertWidth(w,"mm",valueOnly=TRUE),convertHeight(h,"mm",valueOnly=TRUE))*.38
   reference_r<<-maxr
   grid.circle(x,y,r=unit(maxr*sqrt(min(wid[i,j],12)/12),"mm"),gp=gpar(fill="white",col="#666666",lwd=.3))
  }
 }
 split<-factor(vapply(tissues,function(t)as.character(d$block[match(t,d$tissue)]),character(1)),levels=blocks)
 ht<-Heatmap(mat,name="Peak time (ZT)",col=cyc,na_col="#E8E8E8",show_heatmap_legend=FALSE,
    rect_gp=gpar(type="none"),cell_fun=cell,cluster_rows=FALSE,cluster_columns=FALSE,
    row_split=rowgroup,cluster_row_slices=FALSE,row_gap=unit(2,"mm"),
    column_split=split,cluster_column_slices=FALSE,column_gap=unit(3,"mm"),
    row_names_side="left",row_names_gp=gpar(fontsize=15,fontface="italic"),
    row_title=c(paste0("G",seq_along(groups_order)),"Limited\ncoverage"),
    row_title_side="right",row_title_rot=0,row_title_gp=gpar(fontsize=15,fontface="bold"),
    column_names_gp=gpar(fontsize=15),
    column_title_gp=gpar(fontsize=17,fontface="bold"),
    heatmap_legend_param=list(at=seq(0,24,6),direction="horizontal",legend_width=unit(50,"mm"),
                             title_position="topcenter",
                             title_gp=gpar(fontsize=11,fontface="bold"),labels_gp=gpar(fontsize=10)))
 # the disk key is built as a legend object rather than drawn by hand, so it
 # sits on the same row as the colour bar and takes the same title styling
 ref <- c(2,4,8,12)
 disk_lgd <- Legend(labels=paste0(ref," h"),title="95% credible interval width",
   type="points",pch=21,size=unit(5.0*sqrt(ref/12),"mm"),nrow=1,
   legend_gp=gpar(fill="white",col="#666666"),background="white",
   title_gp=gpar(fontsize=16,fontface="bold"),labels_gp=gpar(fontsize=15),
   title_position="topcenter")
 col_lgd <- Legend(col_fun=cyc,title="Peak time (ZT)",at=seq(0,24,6),
   direction="horizontal",legend_width=unit(80,"mm"),title_position="topcenter",
   title_gp=gpar(fontsize=16,fontface="bold"),labels_gp=gpar(fontsize=15))
 both <- packLegend(col_lgd,disk_lgd,direction="horizontal",gap=unit(26,"mm"))
 pdf(file.path(out,paste0("option8_",mode,".pdf")),width=9.375,height=7.4)
 draw(ht,heatmap_legend_list=list(both),
      heatmap_legend_side="bottom",padding=unit(c(10,5,12,5),"mm"),
      column_title=titles[[mode]],column_title_gp=gpar(fontsize=15,fontface="bold"))
 dev.off()
 # the same heatmap on a taller canvas, for the two-by-two layout where it
 # has to stand the same height as the membership panel beside it
 fig_device(file.path(out,paste0("option8_",mode,"_nature.pdf")),width=9.375,height=8.6)
 # draw() seats a column title directly above the blocks, so the heading is
 # written here instead, at the heights the membership panel beside it uses
 draw(ht,heatmap_legend_list=list(both),
      heatmap_legend_side="bottom",padding=unit(c(7.5,5,13.4,2),"mm"))
 upViewport(0)
 grid.text(titles[[mode]],y=unit(1,"npc")-unit(7,"bigpts"),just="top",
           gp=gpar(fontsize=24,fontface="bold"))
 dev.off()

 # the same panel transposed, for a tall column beside A, B and C
 matL <- t(mat); widL <- t(wid)
 cellL <- function(j,i,x,y,w,h,f){
  grid.rect(x,y,w,h,gp=gpar(fill=if(is.finite(matL[i,j]))cyc(matL[i,j]) else "#E8E8E8",col="white",lwd=.45))
  if(is.finite(matL[i,j])){
   mr<-min(convertWidth(w,"mm",valueOnly=TRUE),convertHeight(h,"mm",valueOnly=TRUE))*.38
   grid.circle(x,y,r=unit(mr*sqrt(min(widL[i,j],12)/12),"mm"),gp=gpar(fill="white",col="#666666",lwd=.3))}}
 htL <- Heatmap(matL,name="Peak time (ZT)",col=cyc,na_col="#E8E8E8",
   rect_gp=gpar(type="none"),cell_fun=cellL,cluster_rows=FALSE,cluster_columns=FALSE,
   row_split=split,cluster_row_slices=FALSE,row_gap=unit(3,"mm"),
   column_split=factor(ifelse(as.character(rowgroup)=="Limited coverage","Limited coverage",
                              sub("Phase group ","G",as.character(rowgroup))),
                       levels=c(paste0("G",seq_along(groups_order)),"Limited coverage")),
   cluster_column_slices=FALSE,column_gap=unit(2,"mm"),
   row_title_side="right",row_title_rot=-90,row_title_gp=gpar(fontsize=11,fontface="bold"),
   column_title_side="top",column_title_gp=gpar(fontsize=10,fontface="bold"),
   column_names_side="top",column_names_gp=gpar(fontsize=11,fontface="italic"),
   row_names_gp=gpar(fontsize=11),show_heatmap_legend=FALSE,
   right_annotation=rowAnnotation(sp=anno_empty(border=FALSE,width=unit(3.5,"mm")),
                                  show_annotation_name=FALSE),
   row_title=c("High-concordance cluster","Remaining tissues"))
 col_lgdL <- Legend(col_fun=cyc,title="Peak time (ZT)",at=seq(0,24,6),
   direction="horizontal",legend_width=unit(34,"mm"),title_position="topcenter",
   title_gp=gpar(fontsize=11,fontface="bold"),labels_gp=gpar(fontsize=10))
 disk_lgdL <- Legend(labels=paste0(ref," h"),title="95% credible interval width",
   type="points",pch=21,size=unit(4.0*sqrt(ref/12),"mm"),nrow=1,
   legend_gp=gpar(fill="white",col="#666666"),background="white",
   title_gp=gpar(fontsize=11,fontface="bold"),labels_gp=gpar(fontsize=10),
   title_position="topcenter")
 bothL <- packLegend(col_lgdL,disk_lgdL,direction="horizontal",gap=unit(7,"mm"))
 ht_opt$TITLE_PADDING <- unit(7,"mm")
 pdf(file.path(out,paste0("option8_",mode,"_long.pdf")),width=5.4,height=10.4)
 draw(htL,heatmap_legend_list=list(bothL),heatmap_legend_side="bottom",
      padding=unit(c(6,4,8,4),"mm"),
      column_title=titles[[mode]],column_title_gp=gpar(fontsize=12,fontface="bold"))
 dev.off()
 ht_opt$TITLE_PADDING <- unit(2.5,"mm")
 cat(mode,": k=",k," silhouette=",round(max(scores),3),"\n")
 print(split(core,grp));print(stability)
}
write.csv(do.call(rbind,all_stats),file.path(out,"clustering_comparison.csv"),row.names=FALSE)
write.csv(do.call(rbind,assignments),file.path(out,"phase_group_assignments.csv"),row.names=FALSE)
