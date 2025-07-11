# Name: analyze wgcna
# Author: EY
# Date: 10/04/2023
# Version:4.2.1
# Description: will analyze wgcna network for the developmental expression

#install.packages("abind")
#install.packages("VennDiagram")

library(WGCNA)
library(mltools)
library(data.table)
library(FDRestimation)
library(car)
library(emmeans)
library(gtools)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)
library(goseq)
library(GO.db)
library(tidyr)
library(VennDiagram)
library(RColorBrewer)
source("lettuce/Functions.R")

setwd('/home/ely67071/dev_RNAseq/')


# read in metadata
metadata<-read.csv('lettuce/metadata.csv', row.names=1)

# remove CM
metadata <- metadata[-c(1:4), ]


input_mat_filt<-readRDS("lettuce/wgcna/input_matrix.RData")

#table(dimnames(input_mat_filt)[[1]]==metadata$samples)

# read in the network
netwk<-readRDS("lettuce/wgcna/netwk_bicor_blocksize_default.RData")

# columns and rows of the matrix
nGenes<-ncol(input_mat_filt)
nsamples<-nrow(input_mat_filt)
# turn the labels into colors
mergedColors=labels2colors(netwk$colors)
netwk$colors


# create a table of module sizes
module_sizes <- as.data.frame(table(mergedColors))
module_sizes$mergedColors <-paste("ME", module_sizes$mergedColors,sep="")




MEs0 <- moduleEigengenes(input_mat_filt, mergedColors)$eigengenes
MEs <- orderMEs(MEs0)
help("orderMEs")


df_long <- MEs %>%rownames_to_column(var = "row_id") %>%gather(column_name, value, -row_id)
# Create a new column with modified labels



df_long <- MEs %>% gather(column_name, value, -row_id)
long_df <- pivot_longer(MEs, cols = everything(), names_to = "module", values_to = "value", values_drop_na = TRUE)



# Suppose you have a dataset named 'wide_data' that you want to make long
# 'wide_data' is assumed to be your dataset in wide format

# Use the pivot_longer function to make the dataset long
long_data <- pivot_longer(MEs,
                          cols = -id, # Specify the column(s) to remain unchanged (e.g., ID column)
                          names_to = "variable", # Name of the new column that will contain variable names
                          values_to = "value") # Name of the new column that will contain values



# linear regression 
# process the metadata





metadata$dev_stage <- factor(make.names(metadata$dev_stage),levels = c("VM", "TM", "IM", "IMFM"))
levels(metadata$dev_stage)
MEs$"samples" <-rownames(MEs)
MEs$samples

AllData<-merge(metadata[,c(1,2)], MEs, by="samples")



# think about recoding this
LR_Mod <- function(Yvar,dataset) {
  mod <- lm(Yvar ~dev_stage,
            data=dataset)
  return(mod)
}


LR_mod_results <- lapply(AllData[,c(3:27)], function(x) {LR_Mod(x,AllData)})
str(LR_mod_results)


LR_mod_ANOVA <- lapply(LR_mod_results, function(x) {as.data.frame(Anova(x,test="F",type=2))})
?Anova
LR_mod_ANOVA$MEyellow
# save f and p vals
anova_columns <- lapply(LR_mod_ANOVA, function(x) {x[c(1),c(3,4)]})





anova_wide<-lapply(anova_columns,function(x){as.data.frame(t(x))})
anova_widewlabels<-lapply(names(anova_wide),function(x) {anova_wide[[x]]$Module <- x;return(anova_wide[[x]])} ) 

#f vals
anova_fvals <- lapply(anova_widewlabels, function(x) {x[1,]})

# combine into 1 df
all_fvals <- do.call("rbind", anova_fvals)
colnames(all_fvals) <- c("dev_stage_f","Module")

anova_pvals <- lapply(anova_widewlabels, function(x) {x[2,]})
all_pvals <- do.call("rbind", anova_pvals)
colnames(all_pvals) <- c("dev_stage_p","Module")

anova_results <- merge(all_fvals, all_pvals, by="Module")

# adjust p-value using Bonferroni-Holm (also known as seqeuential bonferroni)
all_pvals_adj<-lapply((anova_results[3]), function(x) {p.adjust(x, method="holm", n=length(x))})
all_pvals_adj<-as.data.frame(do.call(cbind, all_pvals_adj))


# LS means
# for all modules
mod_means<-lapply(LR_mod_results, function(x) {emmeans(x, ~dev_stage, type="response")})
?emmeans

#identical(mod_means,mod_means_test)
mod_means_df<-lapply(mod_means, function(x) {as.data.frame(x)[,c(1:4)]})

mod_means_w_labels<-lapply(names(mod_means_df), function(x) {mod_means_df[[x]]$Module <- x;return(mod_means_df[[x]])})
all_mod_means<-do.call("rbind", mod_means_w_labels)



all_mod_means_wide <-reshape(all_mod_means[,c(1:3,5)], idvar="Module", timevar = "dev_stage", direction="wide")

# combine and save results

all_results<-merge(anova_results, all_mod_means_wide, by="Module")


# adjust p-value using Bonferroni-Holm (also known as seqeuential bonferroni)
all_pvals_adj<-lapply((all_results[3]), function(x) {p.adjust(x, method="holm", n=length(x))})
all_pvals_adj<-as.data.frame(do.call(cbind, all_pvals_adj))

# rename columns
colnames(all_pvals_adj) <- c("dev_stage_p_adj")

# combine the dataframes
all_results_adj<- cbind(all_results, all_pvals_adj)

# turn the p-values into significance stars
all_stars<-lapply((all_results_adj[12]), function(x) {stars.pval(x)})
all_stars<-as.data.frame(do.call(cbind, all_stars))

# name columns
colnames(all_stars) <- c("dev_stage_stars")
all_results_stars <- cbind(all_results_adj, all_stars)

# set the stars as factors
all_results_stars$dev_stage_stars<-factor(all_results_stars$dev_stage_stars, levels=c("***", "**",   "*",  " " ))


# remove the NA
all_results_stars[is.na(all_results_stars)] <- " "

# sort the columns based on signficance (star levels)
all_results_stars_sorted<-all_results_stars %>% arrange(dev_stage_stars)
write.csv(all_results_stars_sorted, file="lettuce/wgcna/wgcna_anova_results.csv", row.names=FALSE)

# subset just the significant results from above
all_results_stars_sorted_sig <- all_results_stars_sorted[ which(all_results_stars_sorted$dev_stage_p_adj < 0.05), ]

# get a characterlist of just the significant modules names
sig_modules<-all_results_stars_sorted_sig$Module


all_data_sig<-AllData[sig_modules]
all_data_sig$dev_stage<-paste0(AllData$dev_stage)

# run t-test on all ANOVA significant results
t_test<-lapply(all_data_sig, function(x) {pairwise.t.test(x, AllData$dev_stage, p.adj="holm")})
print(t_test)


t_test_all<-lapply(AllData[3:27], function(x) {pairwise.t.test(x, AllData$dev_stage, p.adj="holm")})
print(t_test_all)
t_test_all$MEturquoise

# extract module names from the columns
module_names<-(all_results_stars_sorted[,1])
                                                                                                                                                                        
# loop through and create variables to plot each module
plot_list<-list()
i<-0
for(module in module_names){
  i <- i+1

  module_subset <-subset(all_results_stars_sorted, Module %in% module)
  plotting_df<-data.frame(dev_stage=factor(c("VM","TM","IM", "IMFM"),levels = c("VM", "TM", "IM", "IMFM")), eigenval=c(module_subset$emmean.VM, module_subset$emmean.TM, module_subset$emmean.IM, module_subset$emmean.IMFM), se=c(module_subset$SE.VM))
  module_size=subset(module_sizes, mergedColors ==module)
  plot_title<-paste(module_subset$Module, module_size$Freq, sep=": ")
  dev_sig<-paste("pval: ", as.character(module_subset$dev_stage_stars))
  plot_list[[i]]<-ggplot(data=plotting_df, aes(x=dev_stage, y=eigenval, group=1)) + geom_ribbon(aes(ymin = module_subset$SE.VM/2, ymax = module_subset$SE.VM/2), fill = "grey70") +
geom_line(size=1)+ geom_point(size=2)+ geom_errorbar(aes(x= dev_stage,ymin = eigenval-se, ymax = eigenval+se))+ ggtitle(plot_title) + annotate('text',x="VM",y=-0.48,label=dev_sig, size =2, hjust=0)+ ylim(-0.5,0.5)+theme_bw() + theme(plot.title = element_text(size=8),axis.line=element_line(colour = "black"),
                                                                                                                                                                                                                         panel.grid.major=element_blank(),
                                                                                                                                                                                                                         panel.grid.minor = element_blank(),
                                                                                                                                                                                                                         panel.border=element_blank(),
                                                                                                                                                                                                                         panel.background = element_blank())

}


ggplot(data=plotting_df, aes(x=dev_stage, y=eigenval, group=1)) + 
  geom_line(size=1)+geom_point(size=2)+ geom_errorbar(aes(x= dev_stage,ymin = eigenval-se, ymax = eigenval+se))+ggtitle(plot_title) + annotate('text',x="35D",y=0.2,label=dev_sig, size =2)+ ylim(-0.5,0.5)+theme_bw() + theme(plot.title = element_text(size=8),axis.line=element_line(colour = "black"),
                                                                                                                                                        panel.grid.major=element_blank(),
                                                                                                                                                        panel.grid.minor = element_blank(),
                                                                                                                                                        panel.border=element_blank(),
                                                                                                                                                        panel.background = element_blank())








# print individual plots
j<-0
for (module in module_names) {
  j <- j+1
  path_to_plots<-paste('lettuce/plots/wgcna/change_plots/individual_plots/',module,"w_error.png", sep="")
  png(path_to_plots, width=1000, height =1000, res=300)
  print(plot_list[[j]])
  dev.off()
}

# plot them in subsets (55 lots seems way too large)
subset_1<-plot_list[1:12]
png('lettuce/plots/wgcna/change_plots/subset1_w_error.png', width=2000, height =2200, res=300)
plot_1 <- do.call(grid.arrange, subset_1)
dev.off()


subset_2 <- plot_list[13:24]
png('lettuce/plots/wgcna/change_plots/subset2_w_error.png', width=2000, height =2200, res=300)
plot_2 <- do.call(grid.arrange, subset_2)
dev.off()

# salmon ins 25...show on its own

# plot by result category
print(module_names)



# match colors to gene name
color2gene = data.frame(unlist(colnames(input_mat_filt)),unlist(mergedColors))
color2gene$unlist.mergedColors. <- paste0("ME", color2gene$unlist.mergedColors.,sep="")
color2gene_list_data <- split(color2gene, color2gene$unlist.mergedColors.)


# cross check to see what genes in each module are differentially expressed from DESEq2 analysis

DEData_pairwise_cs<-ImportCSVs('lettuce/deseq_results/',0.05)
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
  plot_name <- paste("sunflower/wgcna/plots/change_plots/venn_DE_stages/",module_name,"_venn.png", sep="")
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
png('sunflower/wgcna/plots/change_plots/subset1_w_DE.png', width=2000, height =2200, res=300)
plot_1 <- do.call(grid.arrange, subset_1)
dev.off()


subset_2 <- plot_list_w_DE[13:24]
png('sunflower/wgcna/plots/change_plots/subset2_w_DE.png', width=2000, height =2200, res=300)
plot_2 <- do.call(grid.arrange, subset_2)
dev.off()

subset_3<- plot_list_w_DE[25:36]
png('sunflower/wgcna/plots/change_plots/subset3_w_DE.png', width=2000, height =2200, res=300)
plot_3<-do.call(grid.arrange, subset_3)
dev.off()

subset_4<- plot_list_w_DE[37:48]
png('sunflower/wgcna/plots/change_plots/subset4_w_DE.png', width=2000, height =2200, res=300)
plot_4<-do.call(grid.arrange, subset_4)
dev.off()

subset_5 <- plot_list_w_DE[49:55]
png('sunflower/wgcna/plots/change_plots/subset5_w_DE.png', width=2000, height =2200, res=300)
plot_5<-do.call(grid.arrange, subset_5)
dev.off()



# print individual plots
j<-0
for (module in seq_along(module_names)) {
  j <- j+1
  filename<-paste(module_names[module], "DE",sep='_')
  path_to_plots<-paste('sunflower/wgcna/plots/change_plots/individual_plots/',filename,".png", sep="")
  png(path_to_plots, width=1000, height =1000, res=300)
  print(plot_list_w_DE[[j]])
  dev.off()
}





# GO analysis
# match colors to gene name
color2gene = data.frame(unlist(colnames(input_mat_filt)),unlist(mergedColors))

ha412_genes<-read.table("Ha412GO_terms_interpro_inflo_aa.txt", fill=T)
colnames(ha412_genes)<-c("parent", "ontology_term")

length_test <- ha412_genes$parent
length_test2 <- unique(length_test)

ha412_genes$ontology_term<-as.character(ha412_genes$ontology_term)

test <- ha412_genes %>% dplyr::slice(match(ha412_genes$parent, color2gene$unlist.colnames.input_mat_filt..))
test2<-unique(test)
test3<-inner_join(test2, color2gene, by=c('parent'='unlist.colnames.input_mat_filt..'))
test3<-separate_rows(test3, ontology_term, sep=",")
unique_go_test3<-test3 %>% distinct()

write.csv(unique_go_test3, "wgcna/go_module_data.csv")

# find genes of interest

# CLV3
color2gene %>% filter_all(any_vars(. %in% c("g23024.t1"))) # blue

# WUS
color2gene %>% filter_all(any_vars(. %in% c("g51546.t1"))) # blue

#AGL
color2gene %>% filter_all(any_vars(. %in% c("g56393.t1"))) # blue




# UFO
color2gene %>% filter_all(any_vars(. %in% c("MSTRG.11819")))

# FT --> not present?
color2gene %>% filter_all(any_vars(. %in% c("MSTRG.8447"))) # turquoise

"MSTRG.8514" %in% DEData_pairwise_cs$result_10D_v_20D$Gene

DEData_pairwise_cs$result_30D_v_35D %>% filter_all(any_vars(. %in% c("MSTRG.8514")))

# STM
color2gene %>% filter_all(any_vars(. %in% c("MSTRG.22954")))
color2gene %>% filter_all(any_vars(. %in% c("Ha412HOChr15g00041539"))) # not present


# FD
color2gene %>% filter_all(any_vars(. %in% c("Ha412HOChr09g00022342"))) # not present
color2gene %>% filter_all(any_vars(. %in% c("MSTRG.11790"))) # red
color2gene %>% filter_all(any_vars(. %in% c("MSTRG.21496"))) # red

# LFY
color2gene %>% filter_all(any_vars(. %in% c("MSTRG.12674.2")))

# CO
color2gene %>% filter_all(any_vars(. %in% c("Ha412HOChr14g00040330")))

# AP1
color2gene %>% filter_all(any_vars(. %in% c("MSTRG.10327")))

# to see if in dataset at all
DEData_pairwise_cs$result_20D_v_30D %>% filter_all(any_vars(. %in% c("g33259.t1")))

# Reid genes of interest
color2gene %>% filter_all(any_vars(. %in% c("g62225.t1")))

# Vandana genes of interest
color2gene %>% filter_all(any_vars(. %in% c("g57369.t1")))













