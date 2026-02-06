# Name: plot_mds_fc.R
# Author: EY
# Date: 03/06/2023
# Version:4.2.1
# Description: read in data and plot MDS

library(dplyr)
library(DESeq2)
library(Glimma)
library(sva)
library(ggplot2)
source("sunflower/Functions.R")


# read in and process data

setwd('/home/ely67071/dev_RNAseq/')


fc_raw <- read.delim('lettuce/lettuce_dev_fc_gtf.txt', comment.char="#")


fc_counts <- fc_raw[, 7:ncol(fc_raw)]

# set gene IDs as rownames
rownames(fc_counts) <- fc_raw$Geneid

old_names <- colnames(fc_counts)[1:20]


new_col_names <- sub(".*\\.([A-Za-z0-9]+)Aligned.*", "\\1", old_names)



# apply new names
colnames(fc_counts)[1:20] <- new_col_names


samples=c("CM1", "CM2" ,"CM3", "CM4", "IM1", "IM2", "IM3", "IM4","IMFM1", "IMFM2", "IMFM3", "IMFM4","TM1", "TM2", "TM3", "TM4", "VM1", "VM2", "VM3","VM4")

dev_stage <- sub("[0-9]+$", "", samples)

metadata<-data.frame(samples, dev_stage)
write.csv(as.data.frame(metadata), file='lettuce/metadata.csv')


# create the model (wrt dev_stage)
summed_counts<-DESeqDataSetFromMatrix(fc_counts,colData = metadata, design=~0+dev_stage)


summed_counts$samples
# pre-filter for reads where at least 3 samples have a count of 1 or higher
keep<-rowSums(counts(summed_counts)>=1)>=3
length(which(keep==1))
summed_counts_filt<-summed_counts[keep,]

glimmaMDS(summed_counts_filt, group=metadata)

# get a dataframe of the counts
count_matrix <- as.matrix(counts(summed_counts_filt))

write.csv(as.data.frame((count_matrix)), file='lettuce/raw_lettuce_counts.csv')


















