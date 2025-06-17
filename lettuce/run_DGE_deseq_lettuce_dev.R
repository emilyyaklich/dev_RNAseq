# Name: run_DGE_sunflower_inflo_combatseq.R
# Author: EY
# Date: 08/29/2023
# Version:4.2.1
# Description: Will run the differential gene expression on the infloresence stages accounting for variation with combatseq

library(dplyr)
library(DESeq2)
library(Glimma)
library(sva)

# read in and process data
setwd('/home/ely67071/dev_RNAseq/')

source("lettuce/Functions.R")




raw_counts<-read.csv('lettuce/raw_lettuce_counts.csv', row.names=1)
metadata<-read.csv('lettuce/metadata.csv', row.names=1)
# create the factors of interest
metadata$dev_stage<-factor(metadata$dev_stage)


# create the model with adjusted counts
raw_counts_deseq<-DESeqDataSetFromMatrix((raw_counts),colData = metadata, design=~0+dev_stage)

# run the DGE analysis
DESeq_dataset_results<-DESeq(raw_counts_deseq,parallel=TRUE)



normalized_counts<-counts(DESeq_dataset_results, normalized=TRUE)
write.csv(normalized_counts, "lettuce/deseq_results/normalized_counts_size_factor.csv")




resultsNames(DESeq_dataset_results)

# save this because it will be read in to WGCNA (WGCNA takes all genes into account, but uses the normalization of DESeq)
saveRDS(DESeq_dataset_results, file='lettuce/deseq_results/deseq_dataset_results_pairwise.RData')

# set up the pairwise contrasts (10v20, 20v30, 30v35)
result_VM_v_TM<-results(DESeq_dataset_results,contrast=c("dev_stage","TM","VM"),alpha=0.05,parallel=TRUE)
result_TM_v_CM<-results(DESeq_dataset_results,contrast=c("dev_stage","CM","TM"),alpha=0.05,parallel=TRUE)
result_TM_v_IM<-results(DESeq_dataset_results,contrast=c("dev_stage","IM","TM"),alpha=0.05,parallel=TRUE)
result_CM_v_IM<-results(DESeq_dataset_results,contrast=c("dev_stage","IM","CM"),alpha=0.05,parallel=TRUE)
result_IM_v_IMFM<-results(DESeq_dataset_results,contrast=c("dev_stage","IMFM","IM"),alpha=0.05,parallel=TRUE)

# write to CSV file
write.csv(as.data.frame(result_VM_v_TM), file='lettuce/deseq_results/pairwise/result_VM_v_TM.csv')
write.csv(as.data.frame(result_TM_v_CM), file='lettuce/deseq_results/pairwise/result_TM_v_CM.csv')
write.csv(as.data.frame(result_TM_v_IM), file='lettuce/deseq_results/pairwise/result_TM_v_IM.csv')
write.csv(as.data.frame(result_CM_v_IM), file='lettuce/deseq_results/pairwise/result_CM_v_IM.csv')
write.csv(as.data.frame(result_IM_v_IMFM), file='lettuce/deseq_results/pairwise/result_IM_v_IMFM.csv')


