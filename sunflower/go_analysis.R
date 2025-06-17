# Name: GO analysis
# Author: EY 
# Date: 05/07/2023
# Version:4.1.2
# Description: Will run go analysis with result from run_DGE_deseq_cult_only.R
# need Functions.R written by ED

setwd('/home/ely67071/dev_RNAseq/sunflower/')

library(tidyr)
library(dplyr)
library(goseq)
library(GO.db)
library(matrixStats)
library(stringr)
library(data.table)
source("Functions.R")

ha412_genes<-read.table("/home/ely67071/dev_RNAseq/sunflower_old/Ha412GO_terms_interpro_b3.txt", fill=T)
colnames(ha412_genes)<-c("parent", "ontology_term")
ha412_genes$ontology_term<-as.character(ha412_genes$ontology_term)
go_mapping<-separate_rows(ha412_genes, ontology_term, sep=",")
unique_go_mapping<-go_mapping %>% distinct()
goseq_term_input<-split(unique_go_mapping$parent, unique_go_mapping$ontology_term)


#go_mapping<-go_mapping %>% drop_na()
DEData_treatment<-ImportCSVs('deseq_results/pairwise/',0.05)
all_genes<-lapply(DEData_treatment, function(x) {x$Gene})


test<-intersect(all_genes$combo, ha412_genes$parent)

# read in transcript length and make sure transcript length is the same as all genes
transcript_length<-read.table("/home/ely67071/dev_RNAseq/sunflower_old/transcript_length.csv")
tl_subset<-subset(transcript_length, ID %in% all_genes$combo)
all_genes_subset<-lapply(all_genes,function(x) {x[x %in% transcript_length$ID]})



# filter out significant results
mydataSig_treatment<-lapply(DEData_treatment,SigDEdf,PvaluesCol=7,CritP=0.05)
sig_genes<-lapply(mydataSig_treatment, function(x) {x$Gene})
sig_genes_subset<-lapply(sig_genes,function(x) {x[x %in% transcript_length$ID]})


VMvsFT<-GO_Enrichment(all_genes_subset$result_10D_v_20D, sig_genes_subset$result_10D_v_20D, transcript_length, goseq_term_input)
VMvsFT_go_separate<-separate_rows(VMvsFT,category,sep=",")
write.csv(as.data.frame(VMvsFT_go_separate), file='go_results/VMvFT_go.csv')


FTvIM<-GO_Enrichment(all_genes_subset$result_20D_v_30D, sig_genes_subset$result_20D_v_30D, transcript_length, goseq_term_input)
FTvIM_go_separate<-separate_rows(FTvIM,category,sep=",")
write.csv(as.data.frame(FTvIM_go_separate), file='go_results/FTvIM_go.csv')


IMvIMFM<-GO_Enrichment(all_genes_subset$result_30D_v_35D, sig_genes_subset$result_30D_v_35D, transcript_length, goseq_term_input)
IMvIMFM_go_separate<-separate_rows(IMvIMFM,category,sep=",")
write.csv(as.data.frame(IMvIMFM_go_separate), file='go_results/IMvIMFM_go.csv')
