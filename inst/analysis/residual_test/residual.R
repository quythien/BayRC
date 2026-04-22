# setwd("~/Documents/BayCT")
setwd("/Users/thienpham/Library/CloudStorage/OneDrive-UniversityofPittsburgh/Projects/Circadian/Kyle/Circadian-analysis-main/R/v1/BayRC/Thien/analysis")
observed_para_resid_pval = read.csv("./residual_test/observed_para_resid.csv", row.names = 1)
resid.list = get(load("./residual_test/residual_putamen.RData"))
fitted.list = get(load("./residual_test/fitted_putamen.RData"))

library(org.Hs.eg.db)

# Create a mapping vector
gene_map <- mapIds(org.Hs.eg.db, 
                   keys = rownames(observed_para_resid_pval), 
                   column = "SYMBOL", 
                   keytype = "ENSEMBL", 
                   multiVals = "first")

observed_para_resid_pval$symbol <- gene_map[rownames(observed_para_resid_pval)]

#make visual plots of normality and homoscedasticity
library(gridExtra)
library(ggpubr)
#top result
sorted_observed_para_resid_pval = observed_para_resid_pval[order(observed_para_resid_pval$pvalue), ]
top.genes = rownames(sorted_observed_para_resid_pval)[1:6]

pdf("./residual/top_rhyth_residual_plot.pdf", width = 10, height = 6, onefile = FALSE)
plts.list = list()
for(i in 1:6){
  g.name = top.genes[i]
  g.symbol = sorted_observed_para_resid_pval[g.name, "symbol"]
  display.name = ifelse(is.na(g.symbol), g.name, g.symbol)
  p = ggqqplot(resid.list[[g.name]], font.main = c(14, "bold", "black"),
               title = paste0(display.name,  ": \nR2 = ", round(sorted_observed_para_resid_pval[g.name, "R2"], 2), 
                              ", pval = ", ifelse(round(sorted_observed_para_resid_pval[g.name, "pvalue"], 2)==0, 
                                                  round(sorted_observed_para_resid_pval[g.name, "pvalue"], 8), 
                                                  round(sorted_observed_para_resid_pval[g.name, "pvalue"], 2))))
  plts.list[[i]] = ggplotGrob(p)
}
grid.arrange(grobs = plts.list, ncol = 3)
dev.off()


##residual~TOD
DataListFormatedCRG <- readRDS("~/Documents/BayCT/Code/Multi-Tissue/DataListFormatedCRG.rds")
identical(resid.list[["ENSG00000000003"]] %>% names(), DataListFormatedCRG[["Putamen"]][["data"]] %>% colnames())
tod = DataListFormatedCRG[["Putamen"]][["time"]]
pdf("./residual/top_rhyth_residual_to_TOD.pdf", width = 10, height = 6, onefile = FALSE)
plts.list = list()
for(i in 1:6){
  g.name = top.genes[i]
  g.symbol = sorted_observed_para_resid_pval[g.name, "symbol"]
  display.name = ifelse(is.na(g.symbol), g.name, g.symbol)
  #plot(x = tod, y = resid.list[[g.name]])
  df = as.data.frame(list(Residual = resid.list[[g.name]], 
                          TOD = tod))
  p = ggplot(df, aes(x = Residual))+
    geom_point(aes(y = tod, x = resid.list[[g.name]]))+
    coord_flip() +
    ylab("TOD")+xlab("Residual")+
    ggtitle(display.name)
  plts.list[[i]] = ggplotGrob(p)
}
grid.arrange(grobs = plts.list, ncol = 3)
dev.off()

