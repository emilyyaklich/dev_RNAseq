# Name: analyze DGE deseq sunflower inflo ruvseq
# Author: EY (based off of code written by E. Dittmar)
# Date: 08/29/2023
# Version:4.2.1
# Description: Will analyze the output from the DESeq DGE with combatseq for pairwise dev stages
# need Functions.R written by ED


setwd('/home/ely67071/dev_RNAseq/')

library(dplyr)
library(ggplot2)
library(UpSetR)
library(Glimma)
source("sunflower/Functions.R")
?DEGreport


browseVignettes("DEGreport")
#devtools::install_git("https://git@git.bioconductor.org/packages/DEGreport")
data(humanGender)
library(SummarizedExperiment)
library(ggplot2)
ma <- assays(humanGender)[[1]][1:100,]

des <- colData(humanGender)
des[["other"]] <- sample(c("a", "b"), 85, replace = TRUE)



# now analyze 
# read in the data
DEData_pairwise_cs<-ImportCSVs('sunflower/deseq_results/',0.05)
# filter out significant results
mydataSig_pairwise_cs<-lapply(DEData_pairwise_cs,SigDEdf,PvaluesCol=7,CritP=0.05)

# see which genes overlap (input into upset plot)
SigOverlap_pairwise_cs<-GeneSets(mydataSig_pairwise_cs$result_10D_v_20D[1], mydataSig_pairwise_cs$result_20D_v_30D[1],mydataSig_pairwise_cs$result_30D_v_35D[1])
names(SigOverlap_pairwise_cs)
lapply(SigOverlap_pairwise_cs,function(x) {length(x$Gene)})
SigOverlapGraph_pairwise_cs<-lapply(mydataSig_pairwise_cs, function(x) {x$Gene})

# create an upset plot of DE expression by pairwise dev_stage
png("sunflower/plots/sequential_pairwise_upset.png", res=215, width = 1800, height=1000)
upset(fromList(SigOverlapGraph_pairwise_cs),order.by="freq",nsets=13,nintersects=20, text.scale = 1.5)
dev.off()


# MA plots

deseq <- readRDS('sunflower/deseq_results/deseq_dataset_results_pairwise_combatseq.RData')
deseq$samples





# plot MA data comparing all to day 20

result_10D_v_20D_combatseq<-results(deseq,contrast=c("dev_stage","20D","10D"),alpha=0.05,parallel=TRUE)
result_20D_v_30D_combatseq<-results(deseq,contrast=c("dev_stage","30D","20D"),alpha=0.05,parallel=TRUE)
result_20D_v_35D_combatseq<-results(deseq,contrast=c("dev_stage","35D","20D"),alpha=0.05,parallel=TRUE)

#10 vs 20
ma_data <- plotMA(result_10D_v_20D_combatseq, returnData = TRUE)
# Classify Significant Points by Log2 Fold Change
ma_data$significant <- with(ma_data, ifelse(isDE & lfc > 0, "up", ifelse(isDE & lfc < 0, "down", "not_sig")))

# Create Custom MA Plot with ggplot2
ma_plot<-ggplot(ma_data, aes(x = mean, y = lfc, color = significant)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("up" = "#009E73", "down" = "#D55E00", "not_sig" = "black")) +
  labs(title = "10Dv20D",
       x = "Mean of Normalized Counts",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  scale_x_log10() +
  scale_y_continuous(limits = c(-15, 15)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
ggsave("sunflower/plots/maplot_10v20.png", plot = ma_plot, width = 10, height = 6, dpi = 300)


#20 vs 30
ma_data <- plotMA(result_20D_v_30D_combatseq, returnData = TRUE)
# Classify Significant Points by Log2 Fold Change
ma_data$significant <- with(ma_data, ifelse(isDE & lfc > 0, "up", ifelse(isDE & lfc < 0, "down", "not_sig")))

# Create Custom MA Plot with ggplot2
ma_plot<-ggplot(ma_data, aes(x = mean, y = lfc, color = significant)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("up" = "#009E73", "down" = "#D55E00", "not_sig" = "black")) +
  labs(title = "20Dv30D",
       x = "Mean of Normalized Counts",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  scale_x_log10() +
  scale_y_continuous(limits = c(-15, 15)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
ggsave("sunflower/plots/maplot_20v30.png", plot = ma_plot, width = 10, height = 6, dpi = 300)


#20 vs 35
ma_data <- plotMA(result_20D_v_35D_combatseq, returnData = TRUE)
# Classify Significant Points by Log2 Fold Change
ma_data$significant <- with(ma_data, ifelse(isDE & lfc > 0, "up", ifelse(isDE & lfc < 0, "down", "not_sig")))

# Create Custom MA Plot with ggplot2
ma_plot<-ggplot(ma_data, aes(x = mean, y = lfc, color = significant)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("up" = "#009E73", "down" = "#D55E00", "not_sig" = "black")) +
  labs(title = "20Dv35D",
       x = "Mean of Normalized Counts",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  scale_x_log10() +
  scale_y_continuous(limits = c(-15, 15)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
ggsave("sunflower/plots/maplot_20v35.png", plot = ma_plot, width = 10, height = 6, dpi = 300)


# now 30 vs 35
result_30D_v_35D_combatseq<-results(deseq,contrast=c("dev_stage","35D","30D"),alpha=0.05,parallel=TRUE)

#10 vs 20
ma_data <- plotMA(result_30D_v_35D_combatseq, returnData = TRUE)
# Classify Significant Points by Log2 Fold Change
ma_data$significant <- with(ma_data, ifelse(isDE & lfc > 0, "up", ifelse(isDE & lfc < 0, "down", "not_sig")))

# Create Custom MA Plot with ggplot2
ma_plot<-ggplot(ma_data, aes(x = mean, y = lfc, color = significant)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("up" = "#009E73", "down" = "#D55E00", "not_sig" = "black")) +
  labs(title = "30Dv35D",
       x = "Mean of Normalized Counts",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  scale_x_log10() +
  scale_y_continuous(limits = c(-15, 15)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )
ggsave("sunflower/plots/maplot_30v35.png", plot = ma_plot, width = 10, height = 6, dpi = 300)

