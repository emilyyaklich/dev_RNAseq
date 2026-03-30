# Name: cluster profiler multi
# Author: EY 
# Date: 02/11/2026
# Version:4.2.2
# Description: Will run go analysis using cluster profiler on multiple comparisons


#install.packages("clusterProfiler")
#install.packages("writexl")
library(clusterProfiler)
library(org.Hannuus.eg.db)
library(writexl)
library(ggplot2)


# read in the data...directory to CSV files containing DESeq2 results
DEData<-ImportCSVs('sunflower/deseq_results/pairwise',0.05)

# remove the 10v30 comparison...don't care about this for this comparison...I only want sequential 
DEData$result_10D_v_30D <- NULL

# filter: keep only significant results
mydataSig<-lapply(DEData,SigDEdf,PvaluesCol=7,CritP=0.05)


# different list of df for up or down expression
dataSigup<-lapply(mydataSig, MoreCritNum,column=3, critNum=0)
dataSigdown<-lapply(mydataSig, LessCritNum,column=3, critNum=0)



# this will be the background (all genes expressed) that the enrichment is tested with
# it is the same for all comparisons, so it doesn't matter if I chose 10v20 or 20v30, etc.
gene_universe <- DEData$result_10D_v_20D$Gene




###################################################################
#### cluster compare for up-regulated genes across dev time #######
###################################################################

### function that will run compareCluster ###

run_go_compareCluster <- function(direction = c("up", "down"),
                                  ontology = c("BP", "CC", "MF"),
                                  dataSigup,
                                  dataSigdown,
                                  gene_universe,
                                  OrgDb,
                                  output_dir = "sunflower/go_results/multi")
{
  
  direction <- match.arg(direction)
  ontology  <- match.arg(ontology)
  
  
  dataSig <- if (direction == "up") dataSigup else dataSigdown
  
  gene_list <- list(
    "10v20" = dataSig$result_10D_v_20D$Gene,
    "20v30" = dataSig$result_20D_v_30D$Gene,
    "30v35" = dataSig$result_30D_v_35D$Gene)
  
  # run compareCluster
  go_res <- compareCluster(
    geneCluster = gene_list,
    fun = "enrichGO",
    universe = gene_universe,
    keyType = "GID",
    OrgDb = OrgDb,
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff= 0.05)
  
  
  go_df <- as.data.frame(go_res@compareClusterResult)
  
  result_list <- list(go_df[go_df$Cluster == "10v20", ],go_df[go_df$Cluster == "20v30", ],go_df[go_df$Cluster == "30v35", ])
  
  # rename the lists because these will be the sheets in the spreadsheet
  names(result_list) <- c(
    paste0(direction, "_10v20_", ontology),
    paste0(direction, "_20v30_", ontology),
    paste0(direction, "_30v35_", ontology))
  
  write_xlsx(result_list,file.path(output_dir,paste0("sun_go_", direction, "_", ontology, ".xlsx")))
  
  
  return(go_res)
}

####################################################################################################
### function that will simplify compareCluster result ##############################################
#### not using as it simplifies per comparison, not between all, so can make comparisons harder ####
#simplify_go_compareCluster <- function(go_object,
#                                       direction,
#                                       ontology,
#                                       output_dir = "sunflower/go_results/multi",
#                                       cutoff = 0.7) {
  
#  go_simp <- simplify(go_object,
#                      cutoff = cutoff,
#                      by = "p.adjust",
#                      select_fun = min)
  
#  go_simp_df <- as.data.frame(go_simp@compareClusterResult)
  
#  result_list_simp <- list(go_simp_df[go_simp_df$Cluster == "10v20", ],go_simp_df[go_simp_df$Cluster == "20v30", ],go_simp_df[go_simp_df$Cluster == "30v35", ])
  
#  names(result_list_simp) <- c(
#    paste0(direction, "_10v20_", ontology),
#    paste0(direction, "_20v30_", ontology),
#    paste0(direction, "_30v35_", ontology))
  
#  write_xlsx(result_list_simp,
#             file.path(output_dir,
#                       paste0("sun_go_", direction, "_", ontology, "_simp.xlsx")))
  
  
#  return(go_simp)
#}


#################
##### BP ########
#################

# run BP enrichment on upregulated genes
go_up_bp <- run_go_compareCluster("up", "BP",dataSigup,dataSigdown,gene_universe,org.Hannuus.eg.db)

# make a non-redundant set
#go_up_bp_simp <- simplify_go_compareCluster(go_up_bp,direction = "up",ontology = "BP")


# plot the data
options(enrichplot.colours = c("hotpink", "seagreen")) 
plot <- dotplot(go_up_bp, showCategory = 15) +
  labs(title = "up-regulated enrichment analysis (BP)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("sunflower/go_results/multi/plots/sun_up_regulated_BP_enrichment_15.png", plot = plot, height = 15, width = 10, dpi = 300)





# plot the data
options(enrichplot.colours = c("hotpink", "seagreen")) 
plot <- dotplot(go_up_bp, showCategory = 10) +
  labs(title = "up-regulated enrichment analysis (BP)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("sunflower/go_results/multi/plots/sun_up_regulated_BP_enrichment_10.png", plot = plot, height = 15, width = 10, dpi = 300)



# plot the simplified/non-redundant terms
#options(enrichplot.colours = c("hotpink","seagreen")) 
#plot <- dotplot(go_up_bp_simp, showCategory = 15) +
#  labs(title = "up-regulated enrichment (BP)", x = "dev. stage") +
#  theme_minimal()

#print(plot)

#ggsave("sunflower/go_results/multi/plots/sun_up_regulated_BP_enrichment_simp_15.png", plot = plot, height = 15, width = 10, dpi = 300)



# plot the simplified/non-redundant terms
#options(enrichplot.colours = c("hotpink","seagreen")) 
#plot <- dotplot(go_up_bp_simp, showCategory = 10) +
#  labs(title = "up-regulated enrichment (BP)", x = "dev. stage") +
#  theme_minimal()

#print(plot)

#ggsave("sunflower/go_results/multi/plots/sun_up_regulated_BP_enrichment_simp_10.png", plot = plot, height = 15, width = 10, dpi = 300)



# down
go_down_bp <- run_go_compareCluster("down", "BP",dataSigup,dataSigdown,gene_universe,org.Hannuus.eg.db)


# plot the data
options(enrichplot.colours = c("firebrick2","darkgoldenrod2")) 
plot <- dotplot(go_down_bp, showCategory = 15) +
  labs(title = "down-regulated enrichment analysis (BP)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("sunflower/go_results/multi/plots/sun_down_regulated_BP_enrichment_15.png", plot = plot, height = 15, width = 10, dpi = 300)


# plot the data
options(enrichplot.colours = c("firebrick2","darkgoldenrod2")) 
plot <- dotplot(go_down_bp, showCategory = 10) +
  labs(title = "down-regulated enrichment analysis (BP)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("sunflower/go_results/multi/plots/sun_down_regulated_BP_enrichment_10.png", plot = plot, height = 15, width = 10, dpi = 300)


#################
##### MF ########
#################

go_up_mf <- run_go_compareCluster("up", "MF", dataSigup, dataSigdown, gene_universe, org.Hannuus.eg.db)

# plot the data
options(enrichplot.colours = c("hotpink", "seagreen")) 
plot <- dotplot(go_up_mf, showCategory = 15) +
  labs(title = "up-regulated enrichment analysis (MF)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("sunflower/go_results/multi/plots/sun_up_regulated_MF_enrichment_15.png", plot = plot, height = 15, width = 10, dpi = 300)




# plot the data
options(enrichplot.colours = c("hotpink", "seagreen")) 
plot <- dotplot(go_up_mf, showCategory = 10) +
  labs(title = "up-regulated enrichment analysis (MF)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("sunflower/go_results/multi/plots/sun_up_regulated_MF_enrichment_10.png", plot = plot, height = 15, width = 10, dpi = 300)


# down

go_down_mf <- run_go_compareCluster("down", "MF", dataSigup, dataSigdown, gene_universe, org.Hannuus.eg.db)

# plot the data
options(enrichplot.colours = c("firebrick2", "darkgoldenrod2")) 
plot <- dotplot(go_down_mf, showCategory = 15) +
  labs(title = "down-regulated enrichment analysis (MF)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("sunflower/go_results/multi/plots/sun_down_regulated_MF_enrichment_15.png", plot = plot, height = 15, width = 10, dpi = 300)



# plot the data
options(enrichplot.colours = c("firebrick2", "darkgoldenrod2")) 
plot <- dotplot(go_down_mf, showCategory = 10) +
  labs(title = "down-regulated enrichment analysis (MF)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("sunflower/go_results/multi/plots/sun_down_regulated_MF_enrichment_10.png", plot = plot, height = 15, width = 10, dpi = 300)


#################
##### CC ########
#################

go_up_cc <- run_go_compareCluster("up", "CC", dataSigup, dataSigdown, gene_universe, org.Hannuus.eg.db)


# plot the data
options(enrichplot.colours = c("hotpink", "seagreen")) 
plot <- dotplot(go_up_cc, showCategory = 15) +
  labs(title = "up-regulated enrichment analysis (CC)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("sunflower/go_results/multi/plots/sun_up_regulated_CC_enrichment_15.png", plot = plot, height = 15, width = 10, dpi = 300)


# plot the data
options(enrichplot.colours = c("hotpink", "seagreen")) 
plot <- dotplot(go_up_cc, showCategory = 10) +
  labs(title = "up-regulated enrichment analysis (CC)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("sunflower/go_results/multi/plots/sun_up_regulated_CC_enrichment_10.png", plot = plot, height = 15, width = 10, dpi = 300)


# down

go_down_cc <- run_go_compareCluster("down", "CC", dataSigup, dataSigdown, gene_universe, org.Hannuus.eg.db)

# plot the data
options(enrichplot.colours = c("firebrick2", "darkgoldenrod2")) 
plot <- dotplot(go_down_cc, showCategory = 15) +
  labs(title = "down-regulated enrichment analysis (CC)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("sunflower/go_results/multi/plots/sun_down_regulated_CC_enrichment_15.png", plot = plot, height = 15, width = 10, dpi = 300)


# plot the data
options(enrichplot.colours = c("firebrick2", "darkgoldenrod2")) 
plot <- dotplot(go_down_cc, showCategory = 10) +
  labs(title = "down-regulated enrichment analysis (CC)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("sunflower/go_results/multi/plots/sun_down_regulated_CC_enrichment_10.png", plot = plot, height = 15, width = 10, dpi = 300)








