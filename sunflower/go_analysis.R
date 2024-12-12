# Name: GO analysis
# Author: EY 
# Date: 05/07/2023
# Version:4.1.2
# Description: Will run go analysis with result from run_DGE_deseq_cult_only.R
# need Functions.R written by ED

setwd('/home/ely67071/sunflower_stress_analysis/')




library(tidyr)
library(dplyr)
library(goseq)
library(GO.db)
library(matrixStats)
library(stringr)
library(data.table)
source("Functions.R")

ha412_genes<-read.table("Ha412GO_terms_interpro_aa_slim.txt", fill=T)
colnames(ha412_genes)<-c("parent", "ontology_term")
ha412_genes$ontology_term<-as.character(ha412_genes$ontology_term)
go_mapping<-separate_rows(ha412_genes, ontology_term, sep=",")
unique_go_mapping<-go_mapping %>% distinct()
goseq_term_input<-split(unique_go_mapping$parent, unique_go_mapping$ontology_term)


#go_mapping<-go_mapping %>% drop_na()
DEData_treatment<-ImportCSVs('cult/treatment/',0.05)
all_genes<-lapply(DEData_treatment, function(x) {x$Gene})


test<-intersect(all_genes$combo, ha412_genes$parent)

# read in transcript length and make sure transcript length is the same as all genes
transcript_length<-read.table("transcript_length.csv")
tl_subset<-subset(transcript_length, ID %in% all_genes$combo)
all_genes_subset<-lapply(all_genes,function(x) {x[x %in% transcript_length$ID]})



# filter out significant results
mydataSig_treatment<-lapply(DEData_treatment,SigDEdf,PvaluesCol=7,CritP=0.05)
sig_genes<-lapply(mydataSig_treatment, function(x) {x$Gene})
sig_genes_subset<-lapply(sig_genes,function(x) {x[x %in% transcript_length$ID]})


high_salt_go<-GO_Enrichment(all_genes_subset$HighSalt,sig_genes_subset$HighSalt,tl_subset,goseq_term_input)
high_salt_go_separate<-separate_rows(high_salt_go,category,sep=",")
write.csv(as.data.frame(high_salt_go_separate), file='go_results/go_highsalt.csv')



low_nut_go<-GO_Enrichment(all_genes_subset$LowNut,sig_genes_subset$LowNut,tl_subset,goseq_term_input)
low_nut_go$category
low_nut_go_separate<-separate_rows(low_nut_go,category,sep=",")
write.csv(as.data.frame(low_nut_go_separate), file='go_results/go_lownut.csv')


combo_go<-GO_Enrichment(all_genes_subset$combo,sig_genes_subset$combo,tl_subset,goseq_term_input)
combo_go$category
#combo_go_separate<-separate_rows(combo_go,category,sep=",")
write.csv(as.data.frame(combo_go), file='go_results/go_combo.csv')
all.equal(combo_go,combo_go_separate)



# run MF, BP, and CC separately 
ontology_label <-AnnotationDbi::select(GO.db, keys=go_mapping$ontology_term, keytype="GOID",columns=c('ONTOLOGY'))
unique_ontology_label<-ontology_label %>% distinct()
unique_ontology_label<-na.omit(unique_ontology_label)
mf<-unique_ontology_label[unique_ontology_label$ONTOLOGY=='MF',]
bp<-unique_ontology_label[unique_ontology_label$ONTOLOGY=='BP',]
cc<-unique_ontology_label[unique_ontology_label$ONTOLOGY=='CC',]

bp_input<-c()
for(go in names(goseq_term_input)){
{if(go %in% bp$GOID) 
{bp_input<-append(bp_input, goseq_term_input[go])}
  }}

cc_input<-c()
for(go in names(goseq_term_input)){
  {if(go %in% cc$GOID) 
  {cc_input<-append(cc_input, goseq_term_input[go])}
   }}

mf_input<-c()
for(go in names(goseq_term_input)){
  {if(go %in% mf$GOID) 
  {mf_input<-append(mf_input, goseq_term_input[go])}
  }}


# MF
high_salt_go_mf<-GO_Enrichment(all_genes_subset$HighSalt,sig_genes_subset$HighSalt,tl_subset,mf_input)
high_salt_go_separate_mf<-separate_rows(high_salt_go_mf,category,sep=",")
write.csv(as.data.frame(high_salt_go_separate_mf), file='go_results/go_highsalt_mf.csv')


low_nut_go_mf<-GO_Enrichment(all_genes_subset$LowNut,sig_genes_subset$LowNut,tl_subset,mf_input)
low_nut_go_separate_mf<-separate_rows(low_nut_go_mf,category,sep=",")
write.csv(as.data.frame(low_nut_go_separate_mf), file='go_results/go_lownut_mf.csv')

combo_go_mf<-GO_Enrichment(all_genes_subset$combo,sig_genes_subset$combo,tl_subset,mf_input)
combo_go$category
combo_go_separate_mf<-separate_rows(combo_go_mf,category,sep=",")
write.csv(as.data.frame(combo_go_separate_mf), file='go_results/go_combo_mf.csv')

# BP
high_salt_go_bp<-GO_Enrichment(all_genes_subset$HighSalt,sig_genes_subset$HighSalt,tl_subset,bp_input)
high_salt_go_separate_bp<-separate_rows(high_salt_go_bp,category,sep=",")
write.csv(as.data.frame(high_salt_go_separate_bp), file='go_results/go_highsalt_bp.csv')


low_nut_go_bp<-GO_Enrichment(all_genes_subset$LowNut,sig_genes_subset$LowNut,tl_subset,bp_input)
low_nut_go_separate_bp<-separate_rows(low_nut_go_bp,category,sep=",")
write.csv(as.data.frame(low_nut_go_separate_bp), file='go_results/go_lownut_bp.csv')

combo_go_bp<-GO_Enrichment(all_genes_subset$combo,sig_genes_subset$combo,tl_subset,bp_input)
combo_go$category
combo_go_separate_bp<-separate_rows(combo_go_bp,category,sep=",")
write.csv(as.data.frame(combo_go_separate_bp), file='go_results/go_combo_bp.csv')


# CC
high_salt_go_cc<-GO_Enrichment(all_genes_subset$HighSalt,sig_genes_subset$HighSalt,tl_subset,cc_input)
high_salt_go_separate_cc<-separate_rows(high_salt_go_cc,category,sep=",")
write.csv(as.data.frame(high_salt_go_separate_cc), file='go_results/go_highsalt_cc.csv')


low_nut_go_cc<-GO_Enrichment(all_genes_subset$LowNut,sig_genes_subset$LowNut,tl_subset,cc_input)
low_nut_go_separate_cc<-separate_rows(low_nut_go_cc,category,sep=",")
write.csv(as.data.frame(low_nut_go_separate_cc), file='go_results/go_lownut_cc.csv')

combo_go_cc<-GO_Enrichment(all_genes_subset$combo,sig_genes_subset$combo,tl_subset,cc_input)
combo_go$category
combo_go_separate_cc<-separate_rows(combo_go_cc,category,sep=",")
write.csv(as.data.frame(combo_go_separate_cc), file='go_results/go_combo_cc.csv')



# find the overlap between the experimental groups 

