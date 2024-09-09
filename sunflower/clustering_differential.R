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

#subset_df<-as.data.frame(subset_matrix)

# Add gene names as a column
#subset_df$Gene <- rownames(subset_df)

# Convert to long format
#long_format <- pivot_longer(subset_df, 
#                            cols = -Gene, 
#                            names_to = "Condition", 
#                            values_to = "Expression")


#vsd_df<-as.data.frame(vsd_matrix)

# Add gene names as a column
#vsd_df$Gene <- rownames(vsd_df)

# Create a matrix
#hclust_matrix <- vsd_df %>% 
#  dplyr::select(-Gene) %>% 
 # as.matrix()


# assign rownames
#rownames(hclust_matrix) <- vsd_df$Gene

#hclust_matrix <- hclust_matrix[DE_genes, ]
#dim(hclust_matrix)

#hclust_matrix <- hclust_matrix %>% 
  # transpose the matrix so genes are as columns
 # t() %>% 
  # apply scalling to each column of the matrix (genes)
  #scale() %>% 
  # transpose back so genes are as rows again
  #t()

gene_dist <- dist(scaled_expression_matrix, method="euclidean")

gene_hclust <- hclust(gene_dist, method = "complete")

# The default `plot()` function can be used to produce a simple dendrogram
png("sunflower/plots/clustering_differential/hclust_tree.png", res=215)
plot(gene_hclust, labels = FALSE)
abline(h = 8, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram
dev.off()


gene_cluster <- cutree(gene_hclust, k=8)
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


#vsd_df_cluster <- vsd_df %>% 
#  inner_join(gene_cluster_df, by = "Gene")

#head(vsd_df_cluster)



# Convert to long format
#df_long <- vsd_df_cluster %>%
 # pivot_longer(cols = starts_with(c("X", "HA")), names_to = "Sample", values_to = "Expression") %>%
  #mutate(Stage = case_when(
   # grepl("10D", Sample) ~ "10D",
    #grepl("20D", Sample) ~ "20D",
    #grepl("30D", Sample) ~ "30D",
    #grepl("35D", Sample) ~ "35D"
  #)) %>%
  #group_by(Gene, Stage, cluster) %>%
  #summarize(Average_Expression = mean(Expression), .groups = 'drop')



# Convert to long format
#df_long <- vsd_df_cluster %>%
#  pivot_longer(cols = starts_with(c("X", "HA")), names_to = "Sample", values_to = "Expression")



# Scale the expression values using z-score normalization across all data
# Scale the expression values using z-score normalization across all data
#df_long_average <- df_long %>% group_by(Gene) %>%
#  mutate(Average_Expression_scaled = scale(Average_Expression)) %>% ungroup()

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







# Plot raw expression data
ggplot(df_long, aes(x = samples, y = Expression, group = Gene)) +
  geom_line() +
  facet_wrap(~ cluster) +
  labs(
    x = "Developmental Stage",
    y = "Raw Expression",
    title = "Raw Expression Values Across Clusters"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot scaled expression data
ggplot(df_long_average, aes(x = samples, y = Average_Expression_scaled, group = Gene)) +
  geom_line() +
  facet_wrap(~ cluster) +
  labs(
    x = "Developmental Stage",
    y = "Scaled Expression",
    title = "Scaled Expression Values Across Clusters"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




df_gene_example <- df_long_average %>% filter(Gene == "g4.t1")  # Replace with a specific gene
ggplot(df_gene_example, aes(x = samples, y = Average_Expression_scaled)) +
  labs(
    x = "Developmental Stage",
    y = "Scaled Expression",
    title = "Scaled Expression for Specific Gene"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))








# Unique clusters
unique_clusters <- unique(df_long$cluster)

# Create and print plots for each cluster
plot_list <- list()  # To store plots

for (cl in unique_clusters) {
  # Filter data for the current cluster
  cluster_data <- df_long_average %>% filter(cluster == cl)
  
  # Create the plot
  p <- ggplot(cluster_data, aes(x = samples, y = mean_expression, group = Gene)) +
    geom_line() +
    geom_line(stat = "summary", fun = "mean", color = "brown", size = 1.5, aes(group = 1)) +
    labs(
      x = "Developmental Stage",
      y = "Scaled Expression",
      title = paste("Cluster", cl)
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Print the plot
  print(p)}
  


## Try mfuzz

average_expression_matrix <- cbind(
  D10 = rowMeans(subset_matrix[, c(1, 5, 6)]), 
  D20 = rowMeans(subset_matrix[, c(2, 7, 8)]), 
  D30 = rowMeans(subset_matrix[, c(3, 9, 10)]), 
  D35 = rowMeans(subset_matrix[, c(4, 11, 12)]))




object<-new("ExpressionSet", exprs=average_expression_matrix)
subset_matrix_standard <- standardise(object)
ml=mestimate(subset_matrix_standard)
clusters<-mfuzz(subset_matrix_standard, c=8, m=ml)
?standardise

mfuzz.plot(subset_matrix_standard,cl=clusters,mfrow=c(1,1))
?mfuzz.plot

subset_matrix_standard

as(subset_matrix_standard, "data.frame") -> fff

average_expression_df<-as.data.frame(scaled_expression_matrix)




# try out k-means


set.seed(123)
km.res <- kmeans(average_expression_df, 8, nstart = 25)
print(km.res)

cluster<-km.res$cluster
test2<-as.data.frame(cluster)
# Add gene names as a column
test2$Gene <- rownames(test2)


# Add gene names as a column
average_expression_df$Gene <- rownames(average_expression_matrix)

df_cluster <- average_expression_df %>% 
  inner_join(test2, by = "Gene")


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










