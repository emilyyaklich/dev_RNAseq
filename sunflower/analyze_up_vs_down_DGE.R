# Name: analyze DGE up v down
# Author: EY (based off of code written by E. Dittmar)
# Date: 08/03/2024 
# Version:4.1.2
# Description: Will analyze the output from the DESeq DGE 


setwd('/home/ely67071/dev_RNAseq/')
#install.packages("devtools")
# install package used to plot upset
#devtools::install_github("krassowski/complex-upset", force=TRUE)
#install.packages("lubridate")
library(dplyr)
library(ggplot2)
library(UpSetR)
library(Glimma)
library(ComplexUpset)
library(purrr)
library(naniar)
library(stringr)
source("sunflower/Functions.R")


# read in the data
DEData<-ImportCSVs('sunflower/deseq_results/pairwise',0.05)
# filter out significant results
mydataSig<-lapply(DEData,SigDEdf,PvaluesCol=7,CritP=0.05)

# different list of df for up or down expression
dataSigup<-lapply(mydataSig, MoreCritNum,column=3, critNum=0)
dataSigdown<-lapply(mydataSig, LessCritNum,column=3, critNum=0)


# (2A) Take steps to plot up and down on the same plot with "difference bar"

# add a column labeling "up" or "down" to each list of dataframes 
dataSigup<-lapply(dataSigup, transform, Direction="Up")
dataSigdown<-lapply(dataSigdown, transform, Direction="Down")

# combine the two lists of data frames 
mydataSig<-mapply(rbind, dataSigup, dataSigdown,SIMPLIFY=FALSE)

# subset the dataframes to only contain the gene ID and the direction of expression
gene_direction_only<-lapply(mydataSig, function(x) {x[c("Gene", "Direction")]})
gene_direction_only<-gene_direction_only %>% purrr::reduce(full_join,by=c('Gene'))
# combine all three columns of direction into one string with where direction 
# is outlined. For example. "Up Up Down" If the gene is not in a treatment NA is used
gene_direction_only$dif<-paste(gene_direction_only$Direction.x,gene_direction_only$Direction.y,gene_direction_only$Direction.x.x,gene_direction_only$Direction.y.y, gene_direction_only$Direction)

# loop through all of the genes present and add the genes to the "up", "down",
# or "difference" lists depending on the directionality of expression
difference=c()
up=c()
down=c()
for(i in 1:nrow(gene_direction_only)){
  if(grepl('Down',gene_direction_only[i,"dif"])&grepl("Up",gene_direction_only[i,"dif"])){
    difference<-append(difference,gene_direction_only[i,"Gene"])
  } else if (grepl('Up',gene_direction_only[i,"dif"])) {
    up<-append(up,gene_direction_only[i,"Gene"])
  } else {
    down<-append(down,gene_direction_only[i,"Gene"])
  }
}

#gene_direction_only<-names(gene_direction_only)
# for each list created in the above loop, turn it into a df and then rbind them into one df
difference<-data.frame(difference)
difference<-cbind(difference,x="Difference")
colnames(difference)[1]<-"Gene"
up<-data.frame(up)
up<-cbind(up,x="Up")
colnames(up)[1]<-"Gene"
down<-data.frame(down)
down<-cbind(down,x="Down")
colnames(down)[1]<-"Gene"
# df with Gene ID as column 1 and direction of expression as column 2 
up_or_down<-rbind(down,up,difference)

# in our initial dataframe containing all treatments, add the directionality of the gene
mydataSig$result_10D_v_20D$Direction<-up_or_down$x[match(mydataSig$result_10D_v_20D$Gene, up_or_down$Gene)]
mydataSig$result_20D_v_30D$Direction<-up_or_down$x[match(mydataSig$result_20D_v_30D$Gene, up_or_down$Gene)]
mydataSig$result_30D_v_35D$Direction<-up_or_down$x[match(mydataSig$result_30D_v_35D$Gene, up_or_down$Gene)]



# create a column in each df that contains a treatment variable

mydataSig$result_10D_v_20D<-cbind(mydataSig$result_10D_v_20D, Treatment='result_10D_v_20D')
mydataSig$result_20D_v_30D<-cbind(mydataSig$result_20D_v_30D, Treatment='result_20D_v_30D')
mydataSig$result_30D_v_35D<-cbind(mydataSig$result_30D_v_35D, Treatment='result_30D_v_35D')





# combine all the dataframes on the rows
mydataSig_full<-bind_rows(mydataSig$result_10D_v_20D, mydataSig$result_20D_v_30D, mydataSig$result_30D_v_35D)
comparisons<-c("result_10D_v_20D","result_20D_v_30D","result_30D_v_35D")
treatment<-c("Comparisons")

# one-hot encode the treatment variables
mydataSig_full<-mutate(mydataSig_full, result_10D_v_20D=ifelse(Treatment=="result_10D_v_20D", 1,0))
mydataSig_full<-mutate(mydataSig_full, result_20D_v_30D=ifelse(Treatment=="result_20D_v_30D", 1,0))
mydataSig_full<-mutate(mydataSig_full, result_30D_v_35D=ifelse(Treatment=="result_30D_v_35D", 1,0))



mydataSig_subset<-mydataSig_full[c("Gene","Direction","result_10D_v_20D","result_20D_v_30D","result_30D_v_35D")]

#sig_overlapgraph<-lapply(mydataSig_treatment[1:3], function(x) {x$Gene})

#sig_overlapgraph_colors<-lapply(mydataSig_treatment[1:3], function(x) {x$Direction})

# get the treatment names from the columns
treatments<-colnames(mydataSig_subset)[3:5]
# turn the binary to T or F
mydataSig_subset[treatments]=mydataSig_subset[treatments] ==1
# replace all F with NA
mydataSig_subset<-mydataSig_subset %>% replace_with_na_all(condition=~.x==FALSE)

# combine the columns where the gene ID is the same (keeping treatment identity intact)
mydataSig_plot<-mydataSig_subset %>% group_by(Gene)  %>% summarise_all(list(~ .[!is.na(.)][1]))
# convert the NA back to false (there has to be a better way to do this than
# converting the FALSE to NA, combining and then NA back to false...but I cannot figure it out)
mydataSig_plot[is.na(mydataSig_plot)] <- FALSE
# create a direction variable (stacked bar chart colored by column)
Direction<-colnames(mydataSig_plot)[2]
# plot the upset 
png("sunflower/plots/up_AND_down_DGE.png", width=2000, height=2000, res=300)
ComplexUpset::upset(mydataSig_plot,
                    treatments,name="Treatment",
                    base_annotations=list
                    ("Number of Intersecting Genes"=intersection_size
                      (counts=TRUE,mapping=aes(fill=Direction))
                      +scale_fill_manual(
                        values=c('Up'='#009E73','Down'='#D55E00','Difference'='#CC79A7'))+scale_y_continuous(expand=expansion(mult=c(0,0.1)))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()
                                                                                                                                                    ,axis.line=element_line(colour='black'), axis.text = element_text(color = "black", size = 12),
                                                                                                                                                    axis.title = element_text(color = "black", size = 14),
                                                                                                                                                    plot.title = element_text(color = "black", size = 16, face = "bold"))),
                    set_sizes=upset_set_size(geom=geom_bar(width = 0.4))+theme(axis.line.x = element_line(colour = 'black'),axis.ticks.x =element_line() ),
                    themes=upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=15)))))
dev.off()


# colored upset...no legend
png("sunflower/plots/up_AND_down_white_stripe_no_legend_DGE.png", width=2700, height=2100, res=300)
ComplexUpset::upset(mydataSig_plot,
                    treatments,name="Treatment",stripes="white",
                    base_annotations=list
                    ("Number of Intersecting Genes"=intersection_size
                      (counts=TRUE,mapping=aes(fill=Direction))
                      +scale_fill_manual(
                        values=c('Up'='#009E73','Down'='#D55E00','Difference'='#CC79A7'), guide="none")+scale_y_continuous(expand=expansion(mult=c(0,0.1)))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()
                                                                                                                                                                    ,axis.line=element_line(colour='black'))),
                    set_sizes=upset_set_size(geom=geom_bar(width = 0.4))+theme(axis.line.x = element_line(colour = 'black'),axis.ticks.x =element_line() ),
                    themes=upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=15)))))
dev.off()


# plot the upset, but without the stacked bar chart 
png("sunflower/plots/complexupset_DGE.png", width=2700, height=2100, res=300)
ComplexUpset::upset(mydataSig_plot,
                    treatments,name="Treatment",stripes="white",
                    base_annotations=list
                    ("Number of Intersecting Genes"=intersection_size
                      (counts=TRUE)
                      +scale_y_continuous(expand=expansion(mult=c(0,0.1)))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()
                                                                                 ,axis.line=element_line(colour='black'))),
                    set_sizes=upset_set_size(geom=geom_bar(width = 0.4))+theme(axis.line.x = element_line(colour = 'black'),axis.ticks.x =element_line() ),
                    themes=upset_modify_themes(list('intersections_matrix'=theme(text=element_text(size=15)))))


dev.off()



# look for where WUS
mydataSig_full[mydataSig_full$Gene == "g51546", ] # only 20 vs 30

# look where CLV3 is
mydataSig_full[mydataSig_full$Gene == "g23024", ] # only 20 vs 30

# LFY
mydataSig_full[mydataSig_full$Gene == "g33259", ] # only 20 vs 30


# look for where PIN1
mydataSig_full[mydataSig_full$Gene == "g50085", ] 

# look for where PIN3
mydataSig_full[mydataSig_full$Gene == "g25857", ] 






