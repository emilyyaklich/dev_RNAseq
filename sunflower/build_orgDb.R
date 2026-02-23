# Name: build orgDB
# Author: EY 
# Date: 02/11/2026
# Version:4.2.2
# Description: build database for cluster profiler



# modified from script by E.A. Baldwin

library(AnnotationForge)
setwd('/home/ely67071/dev_RNAseq/sunflower/')


go_terms <- read.csv("GO_getter_result/Ha412HO_v2_cds.fasta.blast.besthits.tsv.genes_to_GO_terms.tsv",sep="\t")

# get rid of transcript IDs 
go_terms$gene_name <- sub("\\..*$", "", go_terms$gene_name)

s_GO <- go_terms[,c(1,3,5)]
s_GO <- unique(s_GO)
# s_GO$EVIDENCE <- "IEA" # arabidopsis evidence code is provided
colnames(s_GO) <- c("GID","GO","EVIDENCE")

makeOrgPackage(go=s_GO,
               version="0.1",
               maintainer="Emily Yaklich <ely67071@uga.edu>",
               author="Emily Yaklich <ely67071@uga.edu>",
               outputDir = "go_results/orgDB/",
               tax_id="4232",
               genus="Helianthus",
               species="annuus",
               goTable="go")

install.packages("go_results/orgDB/org.Hannuus.eg.db/", repos=NULL)
