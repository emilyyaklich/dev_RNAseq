# Name: clustering differential
# Author: EY 
# Date: 09/03/2024 
# Version:4.1.2
# Description: Will cluster the DE output


##### setup ####

setwd('/home/ely67071/dev_RNAseq/')

# load packages
library(tidyverse)
library(tidyr)
library(DESeq2)
source("sunflower/Functions.R")
library(DEGreport)
#install.packages("Mfuzz")
library(Mfuzz)
#install.packages("factoextra")
library(factoextra)
library(cluster)
library(dynamicTreeCut)
#install.packages("fpc")
library(fpc)

# read in metadata
metadata<-read.csv('sunflower/metadata.csv', row.names=1)

# read in deseq results (output from run_DGE_deseq_sunflower_inflo.R)
# note this is reading in all of the data and is NOT filtered for genes that are
# differentially expressed 
deseq <- readRDS('sunflower/deseq_results/deseq_dataset_results_pairwise_combatseq.RData')
deseq$samples


# transform data into a matrix
vsd <- vst(deseq,blind=TRUE)
vsd_matrix<-assay(vsd)


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



average_expression_matrix <- cbind(
  D10 = rowMeans(subset_matrix[, c(1, 5, 6)]), 
  D20 = rowMeans(subset_matrix[, c(2, 7, 8)]), 
  D30 = rowMeans(subset_matrix[, c(3, 9, 10)]), 
  D35 = rowMeans(subset_matrix[, c(4, 11, 12)]))

scaled_expression_matrix <- t(scale(t(average_expression_matrix)))


gene_dist <- dist(scaled_expression_matrix, method="euclidean")

gene_hclust <- hclust(gene_dist, method = "complete")



# The default `plot()` function can be used to produce a simple dendrogram
png("sunflower/plots/clustering_differential/hclust_tree_all_cuts.png", width=2700, height=2100, res=300)
plot(gene_hclust, labels = FALSE)
abline(h = 3.46, col = "red", lwd = 2)
abline(h = 3.1, col = "brown", lwd = 2)
abline(h = 2.9, col = "blue", lwd = 2)# add horizontal line to illustrate cutting dendrogram
abline(h = 2.7, col = "purple", lwd = 2)
abline(h = 2.3, col = "salmon", lwd = 2)
abline(h = 1.5, col = "darkgreen", lwd = 2)
abline(h = 0.4, col = "darkorange", lwd = 2)
dev.off()



png("sunflower/plots/clustering_differential/hclust_tree.png", width=2700, height=2100, res=300)
plot(gene_hclust, labels = FALSE)
abline(h = 3.1, col = "brown", lwd = 2)
dev.off()



# plot tree heights
png("sunflower/plots/clustering_differential/hclust_treeheights.png", width=2700, height=2100, res=300)
heights<-sort(gene_hclust$height, decreasing = TRUE)
plot(heights, type="p", xlab = "Number of Clusters", ylab= "Tree Cut Position")
dev.off()


gene_cluster <- cutree(gene_hclust, h=3.1)
gene_cluster_df <- enframe(gene_cluster)
# Using base R to rename columns
names(gene_cluster_df) <- c("Gene", "cluster")


average_expression_df<-as.data.frame(scaled_expression_matrix)
# Add gene names as a column
average_expression_df$Gene <- rownames(scaled_expression_matrix)


df_cluster <- average_expression_df %>% 
  inner_join(gene_cluster_df, by = "Gene")

df_long <- df_cluster %>%
  pivot_longer(cols = starts_with(c("D")), names_to = "samples", values_to = "Expression")

# Plotting
ggplot(df_long, aes(x = samples, y = Expression, group = Gene)) +
  geom_line() +
  geom_line(stat = "summary", fun = "mean", color = "brown", size = 1.5, aes(group = 1)) +
  facet_wrap(~ cluster) +
  labs(
    x = "Developmental Stage",
    y = "Averaged Expression",
    title = "Faceted Line Plot of Averaged Expression Across Developmental Stages"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


cluster_gene_count <- df_long %>%
  group_by(cluster) %>%
  summarise(number_of_genes = n_distinct(Gene))

# Create a named vector for facet labels
facet_labels <- cluster_gene_count %>%
  mutate(label = paste("Cluster", cluster, ":", number_of_genes)) %>%
  pull(label, name = cluster)

# Ensure facet_labels is a named vector
facet_labels <- setNames(facet_labels, cluster_gene_count$cluster)


# Ensure `samples` is numeric
df_long$samples <- as.numeric(gsub("D", "", df_long$samples))


# Plotting
png("sunflower/plots/clustering_differential/hclust_clusters_3_1_cut.png", width=2700, height=2100, res=300)
ggplot(df_long, aes(x = samples, y = Expression, group = Gene)) +
  geom_line() +
  geom_line(stat = "summary", fun = "mean", color = "brown", size = 1.5, aes(group = 1)) +
  facet_wrap(~ cluster, labeller = labeller(cluster = facet_labels)) +
  scale_x_continuous(
    breaks = c(10, 20, 30, 35),
    labels = c("D10", "D20", "D30", "D35")
  ) +
  labs(
    x = "Developmental Stage",
    y = "Scaled Expression"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


















# just show two clusters CLV3 is in for the lowest cut


# Specify the clusters of interest
clusters_of_interest <- c("94", "41")  # Replace with your cluster IDs

# Subset the data to include only the clusters of interest
df_subset <- df_long %>%
  filter(cluster %in% clusters_of_interest)

# Update the cluster_gene_count and facet_labels for the selected clusters
cluster_gene_count_subset <- df_subset %>%
  group_by(cluster) %>%
  summarise(number_of_genes = n_distinct(Gene))

facet_labels_subset <- cluster_gene_count_subset %>%
  mutate(label = paste("Cluster", cluster, ":", number_of_genes)) %>%
  pull(label, name = cluster)

facet_labels_subset <- setNames(facet_labels_subset, cluster_gene_count_subset$cluster)

# Ensure `samples` is numeric
df_subset$samples <- as.numeric(gsub("D", "", df_subset$samples))

# Plotting
png("sunflower/plots/clustering_differential/hclust_clusters_0_4_cut_subset.png", width=2700, height=2100, res=300)
ggplot(df_subset, aes(x = samples, y = Expression, group = Gene)) +
  geom_line() +
  geom_line(stat = "summary", fun = "mean", color = "darkorange", size = 1.5, aes(group = 1)) +
  facet_wrap(~ cluster, labeller = labeller(cluster = facet_labels_subset)) +
  scale_x_continuous(
    breaks = c(10, 20, 30, 35),
    labels = c("D10", "D20", "D30", "D35")
  ) +
  labs(
    x = "Developmental Stage",
    y = "Scaled Expression"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()





## calculate asw



# Define a range of cut heights
cut_heights <- seq(from =0, to=3.46, by=0.1)

# Function to evaluate cluster quality for a given height
evaluate_clusters_at_height <- function(cut_height, hc, dist_matrix) {
  clusters <- cutree(hc, h = cut_height)  # Cut the dendrogram at height h
  sil <- silhouette(clusters, dist_matrix)  # Compute silhouette scores
  
  # Calculate average silhouette width
  avg_sil_width <- mean(sil[, "sil_width"], na.rm = TRUE)
  
  
  return(list(avg_sil_width = avg_sil_width))
}

# Create a list to store results
cluster_quality_results <- list()

# Evaluate clusters for different heights
for (cut_height in cut_heights) {
  result <- evaluate_clusters_at_height(cut_height, gene_hclust, dist(scaled_expression_matrix))
  cluster_quality_results[[as.character(cut_height)]] <- result
  print(paste("Cut Height =", cut_height, "Average Silhouette Width =", result$avg_sil_width))
}





# test out mfuzz
object<-new("ExpressionSet", exprs=average_expression_matrix)
subset_matrix_standard <- standardise(object)
ml=mestimate(subset_matrix_standard)

# use to determine best number of clusters
Dmin(subset_matrix_standard, crange = seq(2,16,1), repeats=5, m=ml, visu=TRUE)
clusters<-mfuzz(subset_matrix_standard, c=8, m=ml)
?standardise
png("sunflower/plots/clustering_differential/mfuzz_clusters.png", width=2700, height=2100, res=300)
mfuzz.plot2(subset_matrix_standard,cl=clusters,mfrow=c(2,4), x11=FALSE)
dev.off()


as(subset_matrix_standard, "data.frame") -> fff

ddd<- acore(subset_matrix_standard, clusters,min.acore=0.10)

(ddd[[1]])->HaC1
(ddd[[2]])->HaC2
(ddd[[3]])->HaC3
(ddd[[4]])->HaC4
(ddd[[5]])->HaC5
(ddd[[6]])->HaC6
(ddd[[7]])->HaC7
(ddd[[8]])->HaC8

C1<-as.character(HaC1[,1])
C2<-as.character(HaC2[,1])
C3<-as.character(HaC3[,1])
C4<-as.character(HaC4[,1])
C5<-as.character(HaC5[,1])
C6<-as.character(HaC6[,1])
C7<-as.character(HaC7[,1])
C8<-as.character(HaC8[,1])

# see gene membership
'g51546.t1' %in% C8



# try out k-means
average_expression_df<-as.data.frame(scaled_expression_matrix)

set.seed(123)
km.res <- kmeans(average_expression_df, 8, nstart = 25)
print(km.res)

cluster<-km.res$cluster
test2<-as.data.frame(cluster)
# Add gene names as a column
test2$Gene <- rownames(test2)

# nstart - runs k-means 25 times with different random initiations (so it doesn't get stuck in local optimia). default is 1
# iter.max - max num of iterations the algorithm uses to converge. defualt is 10
fviz_nbclust(average_expression_df,kmeans,method="silhouette")

gap_stat<-clusGap(average_expression_df, FUN=kmeans, nstart=25, K.max=14, B=2, iter.max=1000)
fviz_gap_stat(gap_stat)


#pam_wrapper <- function(x,k) {
#  pam(x,k = k)
#}

#fviz_nbclust(average_expression_df,pam_wrapper,method="wss")

# Add gene names as a column
average_expression_df$Gene <- rownames(average_expression_matrix)

df_cluster_kmeans <- average_expression_df %>% 
  inner_join(test2, by = "Gene")


df_long_kmeans <- df_cluster_kmeans %>%
  pivot_longer(cols = starts_with(c("D")), names_to = "samples", values_to = "Expression")



cluster_gene_count_kmeans <- df_long_kmeans %>%
  group_by(cluster) %>%
  summarise(number_of_genes = n_distinct(Gene))

# Create a named vector for facet labels
facet_labels_kmeans <- cluster_gene_count_kmeans %>%
  mutate(label = paste("Cluster", cluster, ":", number_of_genes)) %>%
  pull(label, name = cluster)

# Ensure facet_labels is a named vector
facet_labels_kmeans <- setNames(facet_labels_kmeans, cluster_gene_count_kmeans$cluster)


# Ensure `samples` is numeric
df_long_kmeans$samples <- as.numeric(gsub("D", "", df_long_kmeans$samples))

# Plotting
png("sunflower/plots/clustering_differential/kmeans_clusters.png", width=2700, height=2100, res=300)
ggplot(df_long_kmeans, aes(x = samples, y = Expression, group = Gene)) +
  geom_line() +
  geom_line(stat = "summary", fun = "mean", color = "brown", size = 1.5, aes(group = 1)) +
  facet_wrap(~ cluster, labeller = labeller(cluster = facet_labels)) +
  scale_x_continuous(
    breaks = c(10, 20, 30, 35),
    labels = c("D10", "D20", "D30", "D35")
  ) +
  labs(
    x = "Developmental Stage",
    y = "Scaled Expression"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()





