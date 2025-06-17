# Name: analyze DGE deseq lettuce inflo ruvseq
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
library(DESeq2)
source("lettuce/Functions.R")
?DEGreport


browseVignettes("DEGreport")
#devtools::install_git("https://git@git.bioconductor.org/packages/DEGreport")

library(SummarizedExperiment)
library(ggplot2)





# now analyze 
# read in the data
DEData_pairwise<-ImportCSVs('lettuce/deseq_results/pairwise',0.05)
# filter out significant results
mydataSig_pairwise<-lapply(DEData_pairwise,SigDEdf,PvaluesCol=7,CritP=0.05)




# see which genes overlap (input into upset plot)
SigOverlap_pairwise<-GeneSets_five(mydataSig_pairwise$result_VM_v_TM[1],mydataSig_pairwise$result_TM_v_CM[1], mydataSig_pairwise$result_TM_v_IM[1], mydataSig_pairwise$result_CM_v_IM[1],mydataSig_pairwise$result_IM_v_IMFM[1])
names(SigOverlap_pairwise)
lapply(SigOverlap_pairwise,function(x) {length(x$Gene)})
SigOverlapGraph_pairwise<-lapply(mydataSig_pairwise, function(x) {x$Gene})

# create an upset plot of DE expression by pairwise dev_stage
png("lettuce/plots/sequential_pairwise_upset_all.png", res=215, width = 1800, height=1000)
upset(fromList(SigOverlapGraph_pairwise),order.by="freq",nsets=13,nintersects=20, text.scale = 1.5)
dev.off()


# Dropping CM
mydataSig_pairwise_subset<-mydataSig_pairwise[-1]
mydataSig_pairwise_subset<-mydataSig_pairwise_subset[-2]

# see which genes overlap (input into upset plot)
SigOverlap_pairwise_subset<-GeneSets(mydataSig_pairwise$result_VM_v_TM[1], mydataSig_pairwise$result_TM_v_IM[1],mydataSig_pairwise$result_IM_v_IMFM[1])
names(SigOverlap_pairwise_subset)
lapply(SigOverlap_pairwise_subset,function(x) {length(x$Gene)})
SigOverlapGraph_pairwise_subset<-lapply(mydataSig_pairwise_subset, function(x) {x$Gene})

# create an upset plot of DE expression by pairwise dev_stage
png("lettuce/plots/sequential_pairwise_upset_subset.png", res=215, width = 1800, height=1000)
upset(fromList(SigOverlapGraph_pairwise_subset),order.by="freq",nsets=13,nintersects=20, text.scale = 1.5)
dev.off()







# MA plots

deseq <- readRDS('lettuce/deseq_results/deseq_dataset_results_pairwise.RData')
deseq$samples


# plot MA data comparing consecutive
result_VM_v_TM<-results(deseq,contrast=c("dev_stage","TM","VM"),alpha=0.05,parallel=TRUE)
result_TM_v_IM<-results(deseq,contrast=c("dev_stage","IM","TM"),alpha=0.05,parallel=TRUE)
result_IM_v_IMFM<-results(deseq,contrast=c("dev_stage","IMFM","IM"),alpha=0.05,parallel=TRUE)

#VM vs TM
ma_data <- plotMA(result_VM_v_TM, returnData = TRUE,alpha=0.05)
# Classify Significant Points by Log2 Fold Change
ma_data$significant <- with(ma_data, ifelse(isDE & lfc > 0, "up", ifelse(isDE & lfc < 0, "down", "not_sig")))
#subset(ma_data, lfc < -15 | lfc > 15)

# Create Custom MA Plot with ggplot2
ma_plot<-ggplot(ma_data, aes(x = mean, y = lfc, color = significant)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("up" = "darkgreen", "down" = "purple4", "not_sig" = "black")) +
  labs(title = "VMvTM",
       x = "Mean of Normalized Counts",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  scale_x_log10() +
  scale_y_continuous(limits = c(-15, 15)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    panel.grid = element_blank()
  )
ggsave("lettuce/plots/maplot_VMvTM.png", plot = ma_plot, dpi = 300)
# width=3, height=6 for daniel figures


#TM vs IM
ma_data <- plotMA(result_TM_v_IM, returnData = TRUE)
# Classify Significant Points by Log2 Fold Change
ma_data$significant <- with(ma_data, ifelse(isDE & lfc > 0, "up", ifelse(isDE & lfc < 0, "down", "not_sig")))

# Create Custom MA Plot with ggplot2
ma_plot<-ggplot(ma_data, aes(x = mean, y = lfc, color = significant)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("up" = "darkgreen", "down" = "purple4", "not_sig" = "black")) +
  labs(title = "TMvIM",
       x = "Mean of Normalized Counts",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  scale_x_log10() +
  scale_y_continuous(limits = c(-15, 15)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    panel.grid = element_blank()
  )
ggsave("lettuce/plots/maplot_TMvIM.png", plot = ma_plot, dpi = 300)

#TM vs IM
ma_data <- plotMA(result_IM_v_IMFM, returnData = TRUE)
# Classify Significant Points by Log2 Fold Change
ma_data$significant <- with(ma_data, ifelse(isDE & lfc > 0, "up", ifelse(isDE & lfc < 0, "down", "not_sig")))

# Create Custom MA Plot with ggplot2
ma_plot<-ggplot(ma_data, aes(x = mean, y = lfc, color = significant)) +
  geom_point(alpha = 0.4) +
  scale_color_manual(values = c("up" = "darkgreen", "down" = "purple4", "not_sig" = "black")) +
  labs(title = "IMvIMFM",
       x = "Mean of Normalized Counts",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
  scale_x_log10() +
  scale_y_continuous(limits = c(-15, 15)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    panel.grid = element_blank()
  )
ggsave("lettuce/plots/maplot_IMvIMFM.png", plot = ma_plot, dpi = 300)




