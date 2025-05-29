# Name: load GC data
# Author: EY
# Date: 12/9/2022 Updated: 03/06/2023
# Version:4.1.2
# Description: will load in the parsed gene count data
# This is done using DESeq

# set cwd
setwd('/scratch/ely67071/lettuce_dev_data/gene_count_data/')
library("DESeq2")

# put all of the gene_count files into a character vector in order to load them 
gene_count_files<-dir(pattern="*\\.tab$")

gene_count_files

samples=c("CM1", "CM2" ,"CM3", "CM4", "IM1", "IM2", "IM3", "IM4","IMFM1", "IMFM2", "IMFM3", "IMFM4","TM1", "TM2", "TM3", "TM4", "VM1", "VM2", "VM3","VM4")

# read the data into a dataframe
data_table<-data.frame(sampleName=samples,fileName=gene_count_files)


# load the DESeq data set
# we read in the data using the DESeqDataSetFromHTSeqCount bc our dataset is in the same format as HTSeq data
dds_set<-DESeqDataSetFromHTSeqCount(sampleTable=data_table,directory='/scratch/ely67071/lettuce_dev_data/gene_count_data/',design=~0)

# save the data set 
saveRDS(dds_set,file="/home/ely67071/dev_RNAseq/lettuce/gene_count_lettuce_dev_deseq.Rdata")


