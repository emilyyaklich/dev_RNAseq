# Name: wgcna X deg
# Author: EY 
# Date: 09/03/2024 
# Version:4.4.2
# Description: Will cross check DEG results with WGCNA results

library(WGCNA)
library(mltools)
library(data.table)
library(FDRestimation)
#library(car)
library(emmeans)
library(gtools)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
#library(ggpubr)
library(goseq)
library(GO.db)
library(tidyr)
library(VennDiagram)
library(RColorBrewer)
setwd('/home/ely67071/dev_RNAseq/')

input_mat_filt<-readRDS("sunflower/wgcna/input_matrix.RData")
color2gene_list_data<-readRDS("sunflower/wgcna/color2gene_list_data.RData")

module_names <-names(color2gene_list_data)

all_results_stars_sorted <-read.csv(file='sunflower/wgcna/wgcna_anova_results.csv')

# cross check to see what genes in each module are differentially expressed from DESEq2 analysis

DEData_pairwise_cs<-ImportCSVs('sunflower/deseq_results/pairwise',0.05)
# filter out significant results
mydataSig_pairwise_cs<-lapply(DEData_pairwise_cs,SigDEdf,PvaluesCol=7,CritP=0.05)

for (i in 1:length(mydataSig_pairwise_cs)) {
  assign(paste0("mydatasig", i), as.data.frame(mydataSig_pairwise_cs[[i]]))
}



DE_modules_10v20<-lapply(color2gene_list_data,function(x) {intersect(x$unlist.colnames.input_mat_filt.., mydatasig1$Gene)})
DE_modules_10v20<-DE_modules_10v20[(module_names)]
d1 <- data.frame(module_name=names(DE_modules_10v20), genes=matrix(DE_modules_10v20))
d1$genes <- sapply(d1$genes, paste, collapse=",")
write.csv((d1), file="sunflower/wgcna/DE_modules_10_v_20.csv")




DE_modules_20v30<-lapply(color2gene_list_data,function(x) {intersect(x$unlist.colnames.input_mat_filt.., mydatasig2$Gene)})
DE_modules_20v30<-DE_modules_20v30[(module_names)]
d2 <- data.frame(module_name=names(DE_modules_20v30), genes=matrix(DE_modules_20v30))
d2$genes <- sapply(d2$genes, paste, collapse=",")
write.csv(as.data.frame(d2), file="sunflower/wgcna/DE_modules_20_v_30.csv")



DE_modules_30v35<-lapply(color2gene_list_data,function(x) {intersect(x$unlist.colnames.input_mat_filt.., mydatasig3$Gene)})
DE_modules_30v35<-DE_modules_30v35[(module_names)]
d3 <- data.frame(module_name=names(DE_modules_30v35), genes=matrix(DE_modules_30v35))
d3$genes <- sapply(d3$genes, paste, collapse=",")
write.csv(as.data.frame(d3), file="sunflower/wgcna/DE_modules_30_v_35.csv")


# create a list of the DE genes per module/dev stage that will be plotted into the venn diagram
venn_list<-list()
for(module in seq_along(module_names)){
  print(module_names[module])
  venn_list[[module_names[module]]]<-list(x10v20=c(DE_modules_10v20[module_names[module]]),x20v30=c(DE_modules_20v30[module_names[module]]),x30v35=c(DE_modules_30v35[module_names[module]]))
}




myCol<-brewer.pal(3, "Pastel2")


# individual venn diagram
venn_blue <- list(x10v20=c(venn_list$MEblue$x10v20[[1]]),x20v30=c(venn_list$MEblue$x20v30[[1]]),x30v35=c(venn_list$MEblue$x30v35[[1]]))
venn.diagram(venn_blue,filename="wgcna/blue_venn.png",
             lwd = 2,
             lty = 'blank',
             fill = myCol,
             
             # Numbers
             cex = .8,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.8,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             rotation = 1)




# create venn diagrams for overlap in differentially expressed genes
for(module in seq_along(venn_list)){
  module_name <- module_names[module]
  plot_name <- paste("sunflower/plots/wgcna/change_plots/venn_DE_stages/",module_name,"_venn.png", sep="")
  for(dev_stages in venn_list[module]){
    venn_list_plot<-list(x10v20=c(dev_stages$x10v20[[1]]),x20v30=c(dev_stages$x20v30[[1]]),x30v35=c(dev_stages$x30v35[[1]]))
    venn.diagram(venn_list_plot,filename=plot_name,
                 lwd = 2,
                 lty = 'blank',
                 fill = myCol,
                 
                 # Numbers
                 cex = 1.2,
                 fontface = "bold",
                 fontfamily = "sans",
                 
                 # Set names
                 cat.cex = 1.2,
                 cat.fontface = "bold",
                 cat.default.pos = "outer",
                 cat.pos = c(-27, 27, 135),
                 cat.dist = c(0.055, 0.055, 0.085),
                 cat.fontfamily = "sans",
                 rotation = 1)
    
  }
}








plot_list_w_DE<-list()
i<-0
for(module in seq_along(module_names)){
  print(module_names[module])
  i <- i+1
  # to access lengths of DE genes from DFneed to remove ME from the module name
  mod1<-as.name(sub('..', '', module_names[module]))
  module_subset <-subset(all_results_stars_sorted, Module %in% module_names[module])
  plotting_df<-data.frame(dev_stage=c("10D","20D","30D", "35D"), eigenval=c(module_subset$emmean.X10D, module_subset$emmean.X20D, module_subset$emmean.X30D, module_subset$emmean.X35D))
  module_size=subset(module_sizes, mergedColors ==module_names[module])
  plot_title<-paste(module_subset$Module, module_size$Freq, sep=": ")
  x10v20label<-as.character(length(DE_modules_10v20[[module]]))
  
  x20v30label<-as.character(length(DE_modules_20v30[[module]]))
  x30v35label<-as.character(length(DE_modules_30v35[[module]]))
  print(x10v20label)
  print(x20v30label)
  print(x30v35label)
  dev_sig<-paste("pval: ", as.character(module_subset$dev_stage_stars))
  plot_list_w_DE[[i]]<-ggplot(data=plotting_df, aes(x=dev_stage, y=eigenval, group=1)) + geom_line(size=1)+ geom_point(size=2)+ ggtitle(plot_title) + annotate('text',x="35D",y=0.2,label=dev_sig, size =2)+ ylim(-0.5,0.5)+theme_bw() + theme(plot.title = element_text(size=8),axis.line=element_line(colour = "black"),
                                                                                                                                                                                                                                               panel.grid.major=element_blank(),
                                                                                                                                                                                                                                               panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                               panel.border=element_blank(),
                                                                                                                                                                                                                                               panel.background = element_blank()) + 
    geom_text(show.legend = FALSE, size=3, x = 1.5, y = 0, label = x10v20label, color="red") +
    geom_text(show.legend = FALSE,size=3, x = 2.5, y = 0, label = x20v30label, color="red") +
    geom_text(show.legend = FALSE,size=3, x = 3.5, y = 0, label = x30v35label, color="red")
  print(i)
  
  
}


# plot them in subsets (55 plots seems way too large)
subset_1<-plot_list_w_DE[1:12]
png('sunflower/plots/wgcna/change_plots/subset1_w_DE.png', width=2000, height =2200, res=300)
plot_1 <- do.call(grid.arrange, subset_1)
dev.off()


subset_2 <- plot_list_w_DE[13:24]
png('sunflower/plots/wgcna/change_plots/subset2_w_DE.png', width=2000, height =2200, res=300)
plot_2 <- do.call(grid.arrange, subset_2)
dev.off()

subset_3<- plot_list_w_DE[25:36]
png('sunflower/plots/wgcna/change_plots/subset3_w_DE.png', width=2000, height =2200, res=300)
plot_3<-do.call(grid.arrange, subset_3)
dev.off()

subset_4<- plot_list_w_DE[37:48]
png('sunflower/plots/wgcna/change_plots/subset4_w_DE.png', width=2000, height =2200, res=300)
plot_4<-do.call(grid.arrange, subset_4)
dev.off()

subset_5 <- plot_list_w_DE[49:55]
png('sunflower/plots/wgcna/change_plots/subset5_w_DE.png', width=2000, height =2200, res=300)
plot_5<-do.call(grid.arrange, subset_5)
dev.off()



# print individual plots
j<-0
for (module in seq_along(module_names)) {
  j <- j+1
  filename<-paste(module_names[module], "DE",sep='_')
  path_to_plots<-paste('sunflower/plots/wgcna/change_plots/individual_plots/',filename,".png", sep="")
  png(path_to_plots, width=1000, height =1000, res=300)
  print(plot_list_w_DE[[j]])
  dev.off()
}
