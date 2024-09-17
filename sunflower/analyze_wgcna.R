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

library(tibble)
source("sunflower/Functions.R")

setwd('/home/ely67071/dev_RNAseq/')

browseVignettes("DESeq2")
# read in the cultivated metadata for wgcna
metadata<-read.csv(file='sunflower/metadata.csv')



input_mat_filt<-readRDS("sunflower/wgcna/input_matrix.RData")

#table(dimnames(input_mat_filt)[[1]]==metadata$samples)

# read in the network
netwk<-readRDS("sunflower/wgcna/netwk_bicor_blocksize_default.RData")

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





# Add rownames as a column
df <- rownames_to_column(MEs, var = "samples")

# Convert dataframe to long format
df_long <- gather(df, key = "Column_Name", value = "Value", -samples)

# Initialize an empty vector to store the new values
dev_stage <- rep(NA, nrow(df_long))

patterns_to_check <- c("10D", "20D", "30D", "35D")
for (i in seq_along(df_long$samples)) {
  for (pattern in patterns_to_check) {
    if (grepl(pattern, df_long$samples[i])) {
      # Assign the corresponding value based on the condition
      if (pattern == "10D") {
        dev_stage[i] <- "D10"
      } else if (pattern == "20D") {
        dev_stage[i] <- "D20"
      } else if (pattern == "30D") {
        dev_stage[i] <- "D30"
      } else if (pattern == "35D") {
        dev_stage[i] <- "D35"
      }
      # Break the loop once a match is found
      break
    }
  }
}

# Add the new column to the dataframe
df_long$dev_stage <- dev_stage



# do something similar to above, but with the replicates
patterns_to_check <- c("REP1", "REP2", "REP3", "_1_", "_2_", "_3_")

# Initialize an empty vector to store the new values
reps <- rep(NA, nrow(df_long))

# Loop through the values of the column
for (i in seq_along(df_long$samples)) {
  for (pattern in patterns_to_check) {
    if (grepl(pattern, df_long$samples[i])) {
      # Assign the corresponding value based on the condition
      if (pattern %in% c("REP1", "_1_")) {
        reps[i] <- "rep1"
      } else if (pattern %in% c("REP2", "_2_")) {
        reps[i] <- "rep2"
      } else if (pattern %in% c("REP3", "_3_")) {
        reps[i] <- "rep3"
      }
      # Break the loop once a match is found
      break
    }
  }
}

# Add the new column to the dataframe
df_long$reps <- reps
unique(df_long$Column_Name)

names(df_long)[2] <- "modules"
names(df_long)[3] <- "eigen_gene"


# save the DF
write.csv(df_long, file="sunflower/wgcna/All_data_kristen.csv", row.names=FALSE)



# linear regression 
# process the metadata

metadata$dev_stage<-as.factor(make.names(metadata$dev_stage))
levels(metadata$dev_stage)



rownames(MEs) <- sub("^X", "", rownames(MEs))

AllData<-merge(metadata[,c(2,3)], MEs, by.x="samples", by.y="row.names")


# think about recoding this
LR_Mod <- function(Yvar,dataset) {
  mod <- lm(Yvar ~dev_stage,
            data=dataset)
  return(mod)
}


LR_mod_results <- lapply(AllData[,c(3:58)], function(x) {LR_Mod(x,AllData)})
str(LR_mod_results)


LR_mod_ANOVA <- lapply(LR_mod_results, function(x) {as.data.frame(Anova(x,test="F",type=2))})
?Anova
LR_mod_ANOVA$MEyellow
# save f and p vals
anova_columns <- lapply(LR_mod_ANOVA, function(x) {x[c(1),c(3,4)]})



## testing out Kruskal Wallis test ##

levels(AllData$dev_stage)

LR_mod_ANOVA <- lapply(LR_mod_results, function(x) {as.data.frame(Anova(x,test="F",type=2))})



LR_Mod_np <- function(Yvar,dataset) {
  mod <- lm(rank(Yvar) ~dev_stage,
            data=dataset)
  return(mod)
}

LR_mod_results_np <- lapply(AllData[,c(3:57)], function(x) {LR_Mod_np(x,AllData)})

LR_mod_ANOVA_np <- lapply(LR_mod_results_np, function(x) {as.data.frame(Anova(x,test="F",type=2))})
anova_columns_np <- lapply(LR_mod_ANOVA_np, function(x) {x[c(1),c(3,4)]})

anova_wide_np<-lapply(anova_columns_np,function(x){as.data.frame(t(x))})
anova_widewlabels_np<-lapply(names(anova_wide_np),function(x) {anova_wide_np[[x]]$Module <- x;return(anova_wide_np[[x]])} ) 

#f vals
anova_fvals_np <- lapply(anova_widewlabels_np, function(x) {x[1,]})

# combine into 1 df
all_fvals_np <- do.call("rbind", anova_fvals_np)
colnames(all_fvals_np) <- c("dev_stage_f","Module")

anova_pvals_np <- lapply(anova_widewlabels_np, function(x) {x[2,]})
all_pvals_np <- do.call("rbind", anova_pvals_np)
colnames(all_pvals_np) <- c("dev_stage_p","Module")

anova_results_np <- merge(all_fvals_np, all_pvals_np, by="Module")


### end of testing the KruskalL Wallace test 


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
write.csv(all_results_stars_sorted, file="sunflower/wgcna/wgcna_anova_results.csv", row.names=FALSE)

# subset just the significant results from above
all_results_stars_sorted_sig <- all_results_stars_sorted[ which(all_results_stars_sorted$dev_stage_p_adj < 0.05), ]

# get a characterlist of just the significant modules names
sig_modules<-all_results_stars_sorted_sig$Module


all_data_sig<-AllData[sig_modules]
all_data_sig$dev_stage<-paste0(AllData$dev_stage)

# run t-test on all ANOVA significant results
t_test<-lapply(all_data_sig, function(x) {pairwise.t.test(x, AllData$dev_stage, p.adj="holm")})
print(t_test)


t_test_all<-lapply(AllData[3:58], function(x) {pairwise.t.test(x, AllData$dev_stage, p.adj="holm")})
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
  plotting_df<-data.frame(dev_stage=c("10D","20D","30D", "35D"), eigengene=c(module_subset$emmean.X10D, module_subset$emmean.X20D, module_subset$emmean.X30D, module_subset$emmean.X35D), se=c(module_subset$SE.X10D))
  module_size=subset(module_sizes, mergedColors ==module)
  plot_title<-paste(module_subset$Module, module_size$Freq, sep=": ")
  dev_sig<-paste("pval: ", as.character(module_subset$dev_stage_stars))
  plot_list[[i]]<-ggplot(data=plotting_df, aes(x=dev_stage, y=eigengene, group=1)) + geom_ribbon(aes(ymin = module_subset$SE.X10D/2, ymax = module_subset$SE.X10D/2), fill = "grey70") +
geom_line(size=1)+ geom_point(size=2)+ geom_errorbar(aes(x= dev_stage,ymin = eigengene-se, ymax = eigengene+se))+ ggtitle(plot_title) + annotate('text',x="35D",y=0.2,label=dev_sig, size =2)+ ylim(-0.5,0.5)+theme_bw() + theme(plot.title = element_text(size=8),axis.line=element_line(colour = "black"),
                                                                                                                                                                                                                         panel.grid.major=element_blank(),
                                                                                                                                                                                                                         panel.grid.minor = element_blank(),
                                                                                                                                                                                                                         panel.border=element_blank(),
                                                                                                                                                                                                                         panel.background = element_blank())

}


ggplot(data=plotting_df, aes(x=dev_stage, y=eigengene, group=1)) + 
  geom_line(size=1)+geom_point(size=2)+ geom_errorbar(aes(x= dev_stage,ymin = eigengene-se, ymax = eigengene+se))+ggtitle(plot_title) + annotate('text',x="35D",y=0.2,label=dev_sig, size =2)+ ylim(-0.5,0.5)+theme_bw() + theme(plot.title = element_text(size=8),axis.line=element_line(colour = "black"),
                                                                                                                                                        panel.grid.major=element_blank(),
                                                                                                                                                        panel.grid.minor = element_blank(),
                                                                                                                                                        panel.border=element_blank(),
                                                                                                                                                        panel.background = element_blank())







# print individual plots
j<-0
for (module in module_names) {
  j <- j+1
  path_to_plots<-paste('sunflower/wgcna/plots/change_plots/individual_plots/',module,"w_error.png", sep="")
  png(path_to_plots, width=1000, height =1000, res=300)
  print(plot_list[[j]])
  dev.off()
}

# plot them in subsets (55 lots seems way too large)
subset_1<-plot_list[1:12]
png('sunflower/wgcna/plots/change_plots/subset1_w_error.png', width=2000, height =2200, res=300)
plot_1 <- do.call(grid.arrange, subset_1)
dev.off()


subset_2 <- plot_list[13:24]
png('sunflower/wgcna/plots/change_plots/subset2_w_error.png', width=2000, height =2200, res=300)
plot_2 <- do.call(grid.arrange, subset_2)
dev.off()

subset_3<- plot_list[25:36]
png('sunflower/wgcna/plots/change_plots/subset3_w_error.png', width=2000, height =2200, res=300)
plot_3<-do.call(grid.arrange, subset_3)
dev.off()

subset_4<- plot_list[37:48]
png('sunflower/wgcna/plots/change_plots/subset4_w_error.png', width=2000, height =2200, res=300)
plot_4<-do.call(grid.arrange, subset_4)
dev.off()

subset_5 <- plot_list[49:56]
png('sunflower/wgcna/plots/change_plots/subset5_w_error.png', width=2000, height =2200, res=300)
plot_5<-do.call(grid.arrange, subset_5)
dev.off()


# plot by result category
print(module_names)



# match colors to gene name
color2gene = data.frame(unlist(colnames(input_mat_filt)),unlist(mergedColors))
color2gene$unlist.mergedColors. <- paste0("ME", color2gene$unlist.mergedColors.,sep="")
color2gene_list_data <- split(color2gene, color2gene$unlist.mergedColors.)


# plot individual gene expression patterns in each module (to get an idea for the variance within each module)
input_mat_filt_t<-t(input_mat_filt)
# reorder the columns
colnames(input_mat_filt_t)
col.order <- c("X10D_REP1_ATTACTCG","HA_10D_2_ACCTTGGC","HA_10D_3_ATATCTCG","X20D_REP2_TCCGGAGA","HA_20D_2_GCGCTCTA","HA_20D_3_AACAGGTT","X30D_REP2_CGCTCATT","HA_30D_2_GGTGAACC","HA_30D_3_CAACAATG","X35D_REP1_GAGATTCC","HA_35D_2_TGGTGGCA", "HA_35D_3_AGGCAGAG")
input_mat_filt_t<-input_mat_filt_t[ , col.order]

module_name<-"MEplum1"
genes_in_module <- color2gene$unlist.colnames.input_mat_filt..[color2gene$unlist.mergedColors. == module_name]

#plot CLV3 individually
# Subset the expression matrix based on the gene IDs

subset_expression <- input_mat_filt_t["g23024.t1", ]
subset_expression

# Melt the data frame for plotting with ggplot2
melted_data <- reshape2::melt(subset_expression, id.vars = "gene_id", variable.name = "dev_stage", value.name = "Expression")

# Create a new dataframe with averaged values
averaged_df <- melted_data %>%
  rownames_to_column(var = "Label") %>%  # Convert row names to a column
  mutate(Group = as.numeric(substr(gsub("[^0-9]", "", Label), 1, 2))) %>%  # Extract the numeric part from labels
  group_by(Group) %>%
  summarise(Average_Value = mean(Expression, na.rm = TRUE))

# Plotting with ggplot2
clv3<-ggplot(averaged_df, aes(x = Group, y = Average_Value)) +
  geom_line() +
  geom_point() +
  labs(title = "Expression Changes for CLV3",
       x = "dev_stage",
       y = "Average Expression Count (DESeq2 norm)") +
  scale_x_continuous(breaks = unique(averaged_df$Group))+
  theme_minimal()

png('sunflower/wgcna/plots/raw_expression/clv3_plot.png', width=2000, height =2200, res=300)
print(clv3)
dev.off()








# other random
subset_expression <- input_mat_filt_t["g2492.t1", ]
subset_expression

# Melt the data frame for plotting with ggplot2
melted_data <- reshape2::melt(subset_expression, id.vars = "gene_id", variable.name = "dev_stage", value.name = "Expression")

# Create a new dataframe with averaged values
averaged_df <- melted_data %>%
  rownames_to_column(var = "Label") %>%  # Convert row names to a column
  mutate(Group = as.numeric(substr(gsub("[^0-9]", "", Label), 1, 2))) %>%  # Extract the numeric part from labels
  group_by(Group) %>%
  summarise(Average_Value = mean(Expression, na.rm = TRUE))

# Plotting with ggplot2
clv3<-ggplot(averaged_df, aes(x = Group, y = Average_Value)) +
  geom_line() +
  geom_point() +
  labs(title = "Expression Changes for other",
       x = "dev_stage",
       y = "Average Expression Count (DESeq2 norm)") +
  scale_x_continuous(breaks = unique(averaged_df$Group))+
  theme_minimal()

print(clv3)





# WUS
# Subset the expression matrix based on the gene IDs
subset_expression <- input_mat_filt_t["g51546.t1", ]
subset_expression

# Melt the data frame for plotting with ggplot2
melted_data <- reshape2::melt(subset_expression, id.vars = "gene_id", variable.name = "dev_stage", value.name = "Expression")

# Create a new dataframe with averaged values
averaged_df <- melted_data %>%
  rownames_to_column(var = "Label") %>%  # Convert row names to a column
  mutate(Group = as.numeric(substr(gsub("[^0-9]", "", Label), 1, 2))) %>%  # Extract the numeric part from labels
  group_by(Group) %>%
  summarise(Average_Value = mean(Expression, na.rm = TRUE))

# Plotting with ggplot2
wus<-ggplot(averaged_df, aes(x = Group, y = Average_Value)) +
  geom_line() +
  geom_point() +
  labs(title = "Expression Changes for WUS",
       x = "dev_stage",
       y = "Average Expression Count (DESeq2 norm)") +
  scale_x_continuous(breaks = unique(averaged_df$Group))+
  theme_minimal()

png('sunflower/wgcna/plots/raw_expression/wus_plot.png', width=2000, height =2200, res=300)
print(wus)
dev.off()




#plot ATC individually
# Subset the expression matrix based on the gene IDs
subset_expression <- input_mat_filt_t["g24641.t1", ]
subset_expression

# Melt the data frame for plotting with ggplot2
melted_data <- reshape2::melt(subset_expression, id.vars = "gene_id", variable.name = "dev_stage", value.name = "Expression")

# Create a new dataframe with averaged values
averaged_df <- melted_data %>%
  rownames_to_column(var = "Label") %>%  # Convert row names to a column
  mutate(Group = as.numeric(substr(gsub("[^0-9]", "", Label), 1, 2))) %>%  # Extract the numeric part from labels
  group_by(Group) %>%
  summarise(Average_Value = mean(Expression, na.rm = TRUE))

# Plotting with ggplot2
atc<-ggplot(averaged_df, aes(x = Group, y = Average_Value)) +
  geom_line() +
  geom_point() +
  labs(title = "Expression Changes for ATC",
       x = "dev_stage",
       y = "Average Expression Count (DESeq2 norm)") +
  scale_x_continuous(breaks = unique(averaged_df$Group))+
  theme_minimal()

png('sunflower/wgcna/plots/raw_expression/atc_plot.png', width=2000, height =2200, res=300)
print(atc)
dev.off()


# g26711.t1

# Subset the expression matrix based on the gene IDs
subset_expression <- input_mat_filt_t["g42883.t1", ]
subset_expression

# Melt the data frame for plotting with ggplot2
melted_data <- reshape2::melt(subset_expression, id.vars = "gene_id", variable.name = "dev_stage", value.name = "Expression")

# Create a new dataframe with averaged values
averaged_df <- melted_data %>%
  rownames_to_column(var = "Label") %>%  # Convert row names to a column
  mutate(Group = as.numeric(substr(gsub("[^0-9]", "", Label), 1, 2))) %>%  # Extract the numeric part from labels
  group_by(Group) %>%
  summarise(Average_Value = mean(Expression, na.rm = TRUE))

# Plotting with ggplot2
atc<-ggplot(averaged_df, aes(x = Group, y = Average_Value)) +
  geom_line() +
  geom_point() +
  labs(title = "Expression Changes for PIN3 (g42883.t1)",
       x = "dev_stage",
       y = "Average Expression Count (DESeq2 norm)") +
  scale_x_continuous(breaks = unique(averaged_df$Group))+
  theme_minimal()

png('sunflower/wgcna/plots/raw_expression/PIN3_g42883_plot.png', width=2000, height =2200, res=300)
print(atc)
dev.off()



# cross check to see what genes in each module are differentially expressed from DESEq2 analysis

DEData_pairwise_cs<-ImportCSVs('sunflower/deseq_results/',0.05)
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
  plotting_df<-data.frame(dev_stage=c("10D","20D","30D", "35D"), eigengene=c(module_subset$emmean.X10D, module_subset$emmean.X20D, module_subset$emmean.X30D, module_subset$emmean.X35D))
  module_size=subset(module_sizes, mergedColors ==module_names[module])
  plot_title<-paste(module_subset$Module, module_size$Freq, sep=": ")
  x10v20label<-as.character(length(DE_modules_10v20[[module]]))

  x20v30label<-as.character(length(DE_modules_20v30[[module]]))
  x30v35label<-as.character(length(DE_modules_30v35[[module]]))
  print(x10v20label)
  print(x20v30label)
  print(x30v35label)
  dev_sig<-paste("pval: ", as.character(module_subset$dev_stage_stars))
  plot_list_w_DE[[i]]<-ggplot(data=plotting_df, aes(x=dev_stage, y=eigengene, group=1)) + geom_line(size=1)+ geom_point(size=2)+ ggtitle(plot_title) + annotate('text',x="35D",y=0.2,label=dev_sig, size =2)+ ylim(-0.5,0.5)+theme_bw() + theme(plot.title = element_text(size=8),axis.line=element_line(colour = "black"),
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

subset_5 <- plot_list_w_DE[49:56]
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





# GO analysis - this is matching the genes to modules, so that it can be analyzed for enrichment in GOATOOLs (python script)
# match colors to gene name
color2gene = data.frame(unlist(colnames(input_mat_filt)),unlist(mergedColors))

ha412_genes<-read.table("sunflower/Ha412GO_terms_interpro_b3.txt", fill=T)
colnames(ha412_genes)<-c("parent", "ontology_term")
ha412_genes$ontology_term<-as.character(ha412_genes$ontology_term)
# reorder the rows to match the network module data
ha412_genes_reord <- ha412_genes %>% dplyr::slice(match(ha412_genes$parent, color2gene$unlist.colnames.input_mat_filt..))
# join the two datasets to add colors
ha412_genes_reord_col<-inner_join(ha412_genes_reord, color2gene, by=c('parent'='unlist.colnames.input_mat_filt..'))
# separate each row so that there is only one GO term per row
ha412_genes_reord_sep<-separate_rows(ha412_genes_reord_col, ontology_term, sep=",")
write.csv(ha412_genes_reord_sep, "sunflower/wgcna/go_module_data.csv")

# find genes of interest

# CLV3
color2gene %>% filter_all(any_vars(. %in% c("g23024.t1"))) # brown

# WUS
color2gene %>% filter_all(any_vars(. %in% c("g51546.t1"))) # blue
color2gene %>% filter_all(any_vars(. %in% c("g62225.t1"))) # not in set

# UFO
color2gene %>% filter_all(any_vars(. %in% c("g31393.t1")))
# to see if in DE dataset at all
DEData_pairwise_cs$result_10D_v_30D %>% filter_all(any_vars(. %in% c("g31393.t1"))) # yes but not DE

# LFY
color2gene %>% filter_all(any_vars(. %in% c("g33259.t1"))) # not in set
# to see if in DE dataset at all
DEData_pairwise_cs$result_20D_v_30D%>% filter_all(any_vars(. %in% c("g33259.t1"))) # yes and DE between 20 and 30

# FT 
color2gene %>% filter_all(any_vars(. %in% c("g21886.t1"))) # blue
# g21845.t1,g21858.t1,g50267.t1 are not in network (there are a few orthologs)

# STM
color2gene %>% filter_all(any_vars(. %in% c("g1207.t1"))) # blue
color2gene %>% filter_all(any_vars(. %in% c("g12736.t1"))) # midnight blue
color2gene %>% filter_all(any_vars(. %in% c("g15656.t1"))) # greenyellow
color2gene %>% filter_all(any_vars(. %in% c("g20371.t1"))) # blue
color2gene %>% filter_all(any_vars(. %in% c("g46313.t1"))) # cyan
# g1208.t1,g2243.t1,g24510.t1 are not in network (there are a few orthologs)


# SEP1/2 orthologs
color2gene %>% filter_all(any_vars(. %in% c("g20164.t1"))) # orangered4
# g26945.t1,g32210.t1,g4990.t1,g56848.t1 are not in network (there are a few orthologs)

# AGAMOUS
color2gene %>% filter_all(any_vars(. %in% c("g60940.t1"))) # blue
# g13830.t1,g17773.t1,g32467.t1,g55563.t1,g65380.t1 are not in network (there are a few orthologs)

# FD
color2gene %>% filter_all(any_vars(. %in% c("g31331.t1"))) # blue
color2gene %>% filter_all(any_vars(. %in% c("g56610.t1"))) # light cyan

# AP1
color2gene %>% filter_all(any_vars(. %in% c("g26943.t1"))) # brown
# g24895.t1,g32872.t1,g32872.t2,g43064.t1,g52269.t1,g54628.t1 are not in network (there are a few orthologs)


# CYC2
color2gene %>% filter_all(any_vars(. %in% c("g47961.t1"))) # turquoise
color2gene %>% filter_all(any_vars(. %in% c("g49201.t1"))) # blue
color2gene %>% filter_all(any_vars(. %in% c("g957.t1"))) # turquoise
# g10697.t2 not in network (there are a few orthologs) 




# differential only


# read in the network
netwk<-readRDS("sunflower/wgcna/netwk_bicor_blocksize_default_differential_only.RData")





# columns and rows of the matrix
nGenes<-ncol(hclust_matrix)
nsamples<-nrow(hclust_matrix)
# turn the labels into colors
mergedColors=labels2colors(netwk$colors)
netwk$colors


# create a table of module sizes
module_sizes <- as.data.frame(table(mergedColors))
module_sizes$mergedColors <-paste("ME", module_sizes$mergedColors,sep="")


MEs0 <- moduleEigengenes(hclust_matrix, mergedColors)$eigengenes
MEs <- orderMEs(MEs0)
help("orderMEs")





# Add rownames as a column
df <- rownames_to_column(MEs, var = "samples")

# Convert dataframe to long format
df_long <- gather(df, key = "Column_Name", value = "Value", -samples)

# Initialize an empty vector to store the new values
dev_stage <- rep(NA, nrow(df_long))

patterns_to_check <- c("10D", "20D", "30D", "35D")
for (i in seq_along(df_long$samples)) {
  for (pattern in patterns_to_check) {
    if (grepl(pattern, df_long$samples[i])) {
      # Assign the corresponding value based on the condition
      if (pattern == "10D") {
        dev_stage[i] <- "D10"
      } else if (pattern == "20D") {
        dev_stage[i] <- "D20"
      } else if (pattern == "30D") {
        dev_stage[i] <- "D30"
      } else if (pattern == "35D") {
        dev_stage[i] <- "D35"
      }
      # Break the loop once a match is found
      break
    }
  }
}

# Add the new column to the dataframe
df_long$dev_stage <- dev_stage




# linear regression 
# process the metadata

metadata$dev_stage<-as.factor(make.names(metadata$dev_stage))
levels(metadata$dev_stage)



rownames(MEs) <- sub("^X", "", rownames(MEs))

AllData<-merge(metadata[,c(1,2)], MEs, by.x="samples", by.y="row.names")


# think about recoding this
LR_Mod <- function(Yvar,dataset) {
  mod <- lm(Yvar ~dev_stage,
            data=dataset)
  return(mod)
}


LR_mod_results <- lapply(AllData[,c(3:12)], function(x) {LR_Mod(x,AllData)})
str(LR_mod_results)


LR_mod_ANOVA <- lapply(LR_mod_results, function(x) {as.data.frame(Anova(x,test="F",type=2))})
?Anova
LR_mod_ANOVA$MEyellow
# save f and p vals
anova_columns <- lapply(LR_mod_ANOVA, function(x) {x[c(1),c(3,4)]})



## testing out Kruskal Wallis test ##

levels(AllData$dev_stage)

LR_mod_ANOVA <- lapply(LR_mod_results, function(x) {as.data.frame(Anova(x,test="F",type=2))})



#LR_Mod_np <- function(Yvar,dataset) {
#  mod <- lm(rank(Yvar) ~dev_stage,
#            data=dataset)
 # return(mod)
#}

#LR_mod_results_np <- lapply(AllData[,c(3:57)], function(x) {LR_Mod_np(x,AllData)})

#LR_mod_ANOVA_np <- lapply(LR_mod_results_np, function(x) {as.data.frame(Anova(x,test="F",type=2))})
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


### end of testing the KruskalL Wallace test 



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
write.csv(all_results_stars_sorted, file="sunflower/wgcna/wgcna_anova_results_differential_only.csv", row.names=FALSE)

# subset just the significant results from above
all_results_stars_sorted_sig <- all_results_stars_sorted[ which(all_results_stars_sorted$dev_stage_p_adj < 0.05), ]

# get a characterlist of just the significant modules names
sig_modules<-all_results_stars_sorted_sig$Module


all_data_sig<-AllData[sig_modules]
all_data_sig$dev_stage<-paste0(AllData$dev_stage)

# run t-test on all ANOVA significant results
t_test<-lapply(all_data_sig, function(x) {pairwise.t.test(x, AllData$dev_stage, p.adj="holm")})
print(t_test)


t_test_all<-lapply(AllData[3:58], function(x) {pairwise.t.test(x, AllData$dev_stage, p.adj="holm")})
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
  plotting_df<-data.frame(dev_stage=c("10D","20D","30D", "35D"), eigengene=c(module_subset$emmean.X10D, module_subset$emmean.X20D, module_subset$emmean.X30D, module_subset$emmean.X35D), se=c(module_subset$SE.X10D))
  module_size=subset(module_sizes, mergedColors ==module)
  plot_title<-paste(module_subset$Module, module_size$Freq, sep=": ")
  dev_sig<-paste("pval: ", as.character(module_subset$dev_stage_stars))
  plot_list[[i]]<-ggplot(data=plotting_df, aes(x=dev_stage, y=eigengene, group=1)) + geom_ribbon(aes(ymin = module_subset$SE.X10D/2, ymax = module_subset$SE.X10D/2), fill = "grey70") +
    geom_line(size=1)+ geom_point(size=2)+ geom_errorbar(aes(x= dev_stage,ymin = eigengene-se, ymax = eigengene+se))+ ggtitle(plot_title) + annotate('text',x="35D",y=0.2,label=dev_sig, size =2)+ ylim(-0.5,0.5)+theme_bw() + theme(plot.title = element_text(size=8),axis.line=element_line(colour = "black"),
                                                                                                                                                                                                                                     panel.grid.major=element_blank(),
                                                                                                                                                                                                                                     panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                     panel.border=element_blank(),
                                                                                                                                                                                                                                     panel.background = element_blank())
  
}



# plot them in subsets (55 lots seems way too large)
subset_1<-plot_list[1:12]
png('sunflower/wgcna/plots/change_plots/subset_differential_only.png', width=2000, height =2200, res=300)
plot_1 <- do.call(grid.arrange, subset_1)
dev.off()


hclust_matrix <- t(hclust_matrix)

# match colors to gene name
color2gene = data.frame(unlist(colnames(hclust_matrix)),unlist(mergedColors))
color2gene$unlist.mergedColors. <- paste0("ME", color2gene$unlist.mergedColors.,sep="")
color2gene_list_data <- split(color2gene, color2gene$unlist.mergedColors.)

# CLV3
color2gene %>% filter_all(any_vars(. %in% c("g23024.t1"))) # blue
