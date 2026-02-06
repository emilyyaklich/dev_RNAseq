# Name: run_combatseq_fc.R
# Author: EY
# Date: 03/06/2023 (Modified November 2025)
# Version:4.2.1
# Description: Will run combatseq normalization

library(dplyr)
library(DESeq2)
library(Glimma)
library(sva)
library(ggplot2)
source("sunflower/Functions.R")


# read in and process data

setwd('/home/ely67071/dev_RNAseq/')


fc_raw <- read.delim('sunflower/sunflower_dev_fc_gtf.txt', comment.char="#")


fc_counts <- fc_raw[, 7:ncol(fc_raw)]

# set gene IDs as rownames
rownames(fc_counts) <- fc_raw$Geneid

old_names <- colnames(fc_counts)[1:12]


new_col_names <- sub(".*STAR_output_gtf\\.([A-Za-z0-9_]+(?:_[A-Za-z0-9]+)*).*", "\\1", old_names)


# apply new names
colnames(fc_counts)[1:12] <- new_col_names


samples=c("10D_REP1_ATTACTCG", "20D_REP2_TCCGGAGA" ,"30D_REP2_CGCTCATT", "35D_REP1_GAGATTCC", 
          "HA_10D_2_ACCTTGGC", "HA_10D_3_ATATCTCG", "HA_20D_2_GCGCTCTA", 
          "HA_20D_3_AACAGGTT", "HA_30D_2_GGTGAACC", "HA_30D_3_CAACAATG", "HA_35D_2_TGGTGGCA", "HA_35D_3_AGGCAGAG")

dev_stage<-sub(".*([0-9]{2,2}D).*", "\\1",samples)
dev_stage <- as.factor(dev_stage)

metadata<-data.frame(samples, dev_stage)

# create the factors of interest
metadata$dev_stage<-factor(metadata$dev_stage)

# save metadata
write.csv(as.data.frame(metadata), file='sunflower/metadata.csv')


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

# sort by batch (sample group) and group (dev_stage)
batch <- c(1,1,1,1,2,3,2,3,2,3,2,3)
group <- c(1,2,3,4,1,1,2,2,3,3,4,4)

# adjust counds using combat seq...output is a matrix of adjusted counts 
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group)

# write to a CSV file
write.csv(as.data.frame((adjusted_counts)), file='sunflower/adjusted_counts_combatseq_fc_gtf.csv')

# plot MDS of adjusted counts...can compare this with previous plot pre-combat seq
glimmaMDS(adjusted_counts, group=metadata)

















