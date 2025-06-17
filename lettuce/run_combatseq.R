# Name: run_combatseq.R
# Author: EY
# Date: 03/06/2023
# Version:4.2.1
# Description: Will run combatseq normalization

library(dplyr)
library(DESeq2)
library(Glimma)
library(sva)
library(ggplot2)
source("lettuce/Functions.R")


# read in and process data
setwd('/home/ely67071/dev_RNAseq/')

# read in the data matrix
summed_counts<-readRDS("/home/ely67071/dev_RNAseq/lettuce/gene_count_lettuce_dev_deseq.Rdata")
samples=c("CM1", "CM2" ,"CM3", "CM4", "IM1", "IM2", "IM3", "IM4","IMFM1", "IMFM2", "IMFM3", "IMFM4","TM1", "TM2", "TM3", "TM4", "VM1", "VM2", "VM3","VM4")



dev_stage <- sub("[0-9]+$", "", samples)


metadata<-data.frame(samples, dev_stage)

# create the factors of interest
metadata$dev_stage<-factor(metadata$dev_stage)

# create the model (wrt dev_stage)
summed_counts<-DESeqDataSetFromMatrix(counts(summed_counts),colData = metadata, design=~0+dev_stage)


summed_counts$samples
# pre-filter for reads where at least 3 samples have a count of 1 or higher
keep<-rowSums(counts(summed_counts)>=1)>=3
length(which(keep==1))
summed_counts_filt<-summed_counts[keep,]

# get a dataframe of the counts
count_matrix <- as.matrix(counts(summed_counts_filt))


write.csv(as.data.frame((count_matrix)), file='lettuce/raw_lettuce_counts.csv')

# sort by batch (sample group) and group (dev_stage)
batch <- c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
group <- c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5)

# adjust counds using combat seq...output is a matrix of adjusted counts 
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group)

# write to a CSV file
write.csv(as.data.frame((adjusted_counts)), file='lettuce/adjusted_counts_combatseq.csv')

# plot MDS of adjusted counts...can compare this with previous plot pre-combat seq
glimmaMDS(adjusted_counts, group=metadata)




                     

