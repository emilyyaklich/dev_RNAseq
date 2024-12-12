# Name: wgcna
# Author: EY
# Date: 06/20/2023
# Version:4.2.1
# Description: will run wgcna analysis with the sunflower inflo dev data

setwd('/home/ely67071/dev_RNAseq/')

#install.packages("WGCNA")
library(DESeq2)
library(WGCNA)
library(tidyverse)
library(BiocParallel)
source("sunflower/Functions.R")

# read in metadata
metadata<-read.csv('sunflower/metadata.csv', row.names=1)

#
metadata$samples[1:4] <- paste0("X", metadata$samples[1:4])


# read in deseq results (output from run_DGE_deseq_sunflower_inflo.R)
# note this is reading in all of the data and is NOT filtered for genes that are
# differentially expressed 
deseq <- readRDS('sunflower/deseq_results/deseq_dataset_results_pairwise_combatseq.RData')
deseq$samples



?assay

# pre-filter for reads where at least 10 samples have a count of 10 or higher
keep<-rowSums(counts(deseq)>=10)>=10
length(which(keep==1))
deseq_filt<-deseq[keep,]

# transform data into a matrix
vsd <- vst(deseq_filt,blind=TRUE)
vsd_matrix<-assay(vsd)


?vst
input_mat=t(vsd_matrix)
vsd_matrix[1:5,1:10]

# check samples
gsg<-goodSamplesGenes(input_mat,verbose=3)
gsg$allOK # this statement gives TRUE, so no further removal is necessary


# see if there are outliers to cut in teh dendrogram
sampleTree<-hclust(dist(input_mat), method="average")
par(cex=0.6)
par(mar=c(0,4,2,0))
png("sunflower/wgcna/plots/sample_clustering_outliers.png")
plot(sampleTree, main= "Sample clustering to detect outliers", sub="",xlab="",cex.lab=1.5,cex.axis=1.5,cex.main=2)
abline(h=200, col="red")
dev.off()

# remove the outlier samples
clust<-cutreeStatic(sampleTree,cutHeight = 200,minSize=10)
table(clust)
keepSamples<-(clust==1)
input_mat_filt<-input_mat[keepSamples,]
saveRDS(input_mat_filt,file='sunflower/wgcna/input_matrix.RData')


# remove the outliers from the metadata as well
#remove_outlier<-setdiff(rownames(input_mat),rownames(input_mat_filt))
#metadata<- metadata[!metadata$Plant. %in% remove_outlier,]
#write.csv(metadata, "wgcna/metadata_wgcna.csv")


# get traits of interest from metadata (aka Treatment)
samples<-rownames(input_mat_filt)
trait_rows<-match(samples, metadata$samples)
datTraits<-metadata[trait_rows,-1]
rownames(datTraits)<-metadata[trait_rows,1]
metadata[trait_rows,1]
metadata[trait_rows,-1]
# create treatment labels based after the outliers have been removed
treatment_labels<-paste(metadata$Treatment,"-",metadata$Plant.)

# see if there are global differences between samples 
sampleTree_filt<-hclust(dist(input_mat_filt), method="average")
sample_names<-rownames(input_mat_filt)

# based on stage
traitcolors<-labels2colors(metadata$dev_stage)
png("sunflower/wgcna/plots/sample_relatedness_clusters.png")
plotDendroAndColors(sampleTree_filt,traitcolors)
# no strong clusters indicating globally different groups of samples 
dev.off()

# ID power
#allowWGCNAThreads()
powers=c(c(1:10),seq(from=12, to=20, by=2))
sft=pickSoftThreshold(input_mat_filt,powerVector = powers, verbose=5)
saveRDS(sft,file='sunflower/wgcna/sft_power.RData')
sft<-readRDS("sunflower/wgcna/sft_power.RData")

# plot diagnostics
png("sunflower/wgcna/plots/soft_threshold.png")
par(mfrow=c(1,2))
cex1=0.9
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit, Signed R^2',main=paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main=paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()

# build the network
picked_power=18
temp_cor <- cor
cor <- WGCNA::cor
netwk <- blockwiseModules(input_mat_filt,power=picked_power,networkType = "signed", 
                          deepSplit=2, pamRespectsDendro = F, minModuleSize = 30, 
                          reassignThreshold = 0,mergeCutHeight = 0.25,
                          saveTOMs = T, saveTOMFileBase = "cult",numericLabels = T, verbose=3,corType = "bicor"
                          )
saveRDS(netwk,file='sunflower/wgcna/netwk_bicor_blocksize_default.RData')
netwk<-readRDS("sunflower/wgcna/netwk_bicor_blocksize_default.RData")

# look at results
mergedColors=labels2colors(netwk$colors)
png("sunflower/wgcna/plots/network_modules_bicor_blocksize_default.png")
plotDendroAndColors(netwk$dendrograms[[1]],mergedColors[netwk$blockGenes[[1]]], "Module colors", dendroLabels = FALSE,hang=0.03, addGuide = TRUE,guideHang = 0.05)
dev.off()

# print a text files of the genes/modules
module_df<-data.frame(gene_id=names(netwk$colors),colors=labels2colors(netwk$colors))
write_delim(module_df,file="sunflower/wgcna/gene_modules.txt", delim='\t')



## Try running on differentially expressed genes only 


# now analyze 
# read in the data
DEData_pairwise_cs<-ImportCSVs('sunflower/deseq_results/pairwise/',0.05)
# filter out significant results
mydataSig_pairwise_cs<-lapply(DEData_pairwise_cs,SigDEdf,PvaluesCol=7,CritP=0.05)




# Initialize an empty vector to store gene names
DE_genes <- c()

# Loop through each dataframe in the list
for (df in mydataSig_pairwise_cs) {
  # Filter rows where padj < 0.05 and extract the gene names
  filtered_genes <- df$Gene[df$padj < 0.05]
  # Append the gene names to the vector
  DE_genes <- c(DE_genes, filtered_genes)
}

# Get unique gene names
DE_genes <- unique(DE_genes)

subset_matrix <- vsd_matrix[rownames(vsd_matrix) %in% DE_genes, ]

subset_df<-as.data.frame(subset_matrix)

# Add gene names as a column
subset_df$Gene <- rownames(subset_df)

# Convert to long format
long_format <- pivot_longer(subset_df, 
                            cols = -Gene, 
                            names_to = "Condition", 
                            values_to = "Expression")


vsd_df<-as.data.frame(vsd_matrix)

# Add gene names as a column
vsd_df$Gene <- rownames(vsd_df)

# Create a matrix
hclust_matrix <- vsd_df %>% 
  dplyr::select(-Gene) %>% 
  as.matrix()



# assign rownames
rownames(hclust_matrix) <- vsd_df$Gene
hclust_matrix <- hclust_matrix[DE_genes, ]
hclust_matrix<-t(hclust_matrix)

# get traits of interest from metadata (aka Treatment)
samples<-rownames(hclust_matrix)
trait_rows<-match(samples, metadata$samples)
datTraits<-metadata[trait_rows,-1]
rownames(datTraits)<-metadata[trait_rows,1]

# create treatment labels based after the outliers have been removed
treatment_labels<-paste(metadata$Treatment,"-",metadata$Plant.)

# see if there are global differences between samples 
sampleTree_filt<-hclust(dist(hclust_matrix), method="average")
sample_names<-rownames(hclust_matrix)

# based on stage
traitcolors<-labels2colors(metadata$dev_stage)
plotDendroAndColors(sampleTree_filt,traitcolors)



# ID power
#allowWGCNAThreads()
powers=c(c(1:10),seq(from=12, to=20, by=2))
sft=pickSoftThreshold(input_mat_filt,powerVector = powers, verbose=5)


par(mfrow=c(1,2))
cex1=0.9
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit, Signed R^2',main=paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main=paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

picked_power=18
temp_cor <- cor
cor <- WGCNA::cor
netwk <- blockwiseModules(hclust_matrix,power=picked_power,networkType = "signed", 
                          deepSplit=2, pamRespectsDendro = F, minModuleSize = 30, 
                          reassignThreshold = 0,mergeCutHeight = 0.25,
                          saveTOMs = T, saveTOMFileBase = "cult",numericLabels = T, verbose=3,corType = "bicor"
)

mergedColors=labels2colors(netwk$colors)
plotDendroAndColors(netwk$dendrograms[[1]],mergedColors[netwk$blockGenes[[1]]], "Module colors", dendroLabels = FALSE,hang=0.03, addGuide = TRUE,guideHang = 0.05)

saveRDS(netwk,file='sunflower/wgcna/netwk_bicor_blocksize_default_differential_only.RData')


