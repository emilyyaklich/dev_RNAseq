# Name: plot_MDS
# Author: EY
# Date: 12/24/2022 (edited 08/21/2023)
# Version:4.2.1
# Description: will read the data matrix from load_GC_data_and_sum_reps into
# DESeq2, read in metadata,pre-filter and plotMDS

library(readxl)
library(DESeq2)
library(dplyr)
library(Glimma)
library(sva)
#install.packages("Glimma")

setwd('/home/ely67071/dev_RNAseq/')

# read in the data matrix
summed_counts<-readRDS("/home/ely67071/dev_RNAseq/lettuce/gene_count_lettuce_dev_deseq.Rdata")
dim(summed_counts)


summed_counts

samples=c("CM1", "CM2" ,"CM3", "CM4", "IM1", "IM2", "IM3", "IM4","IMFM1", "IMFM2", "IMFM3", "IMFM4","TM1", "TM2", "TM3", "TM4", "VM1", "VM2", "VM3","VM4")

dev_stage <- sub("[0-9]+$", "", samples)

metadata<-data.frame(samples, dev_stage)
write.csv(as.data.frame(metadata), file='lettuce/metadata.csv')

# create the factors of interest
metadata$dev_stage<-factor(metadata$dev_stage)



# create the model
summed_counts<-DESeqDataSetFromMatrix(counts(summed_counts),colData = metadata, design=~0+dev_stage)
summed_counts$samples
# pre-filter for reads where at least 3 samples (in ANY STAGE) have a count of 1 or higher
keep<-rowSums(counts(summed_counts)>=1)>=3
length(which(keep==1))
summed_counts_filt<-summed_counts[keep,]

counts(summed_counts)>=1

# will load the plot...need to save within the html
glimmaMDS(summed_counts_filt, groups=metadata)

# plot PCA 
#png("plots/pca_raw_data.png", res=215, width = 1200, height=1000)
#par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
#DESeq2::plotPCA((summed_counts_filt),labels=FALSE,col=dev_stage)
#legend("topright",inset=c(-0.4,0),legend=unique(dev_stage), fill=dev_stage)
#dev.off()

dds_set
