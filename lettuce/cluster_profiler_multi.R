# Name: cluster profiler multi
# Author: EY 
# Date: 02/11/2026
# Version:4.2.2
# Description: Will run go analysis using cluster profiler on multiple comparisons


#install.packages("clusterProfiler")
#install.packages("writexl")
library(clusterProfiler)
library(org.Lsativa.eg.db)
library(writexl)
library(ggplot2)


setwd('/home/ely67071/dev_RNAseq/')

# read in the data...directory to CSV files containing DESeq2 results
DEData<-ImportCSVs('lettuce/deseq_results/pairwise',0.05)

# remove the comparisons with CM...I don't care about those right now
DEData$result_CM_v_IM <- NULL
DEData$result_TM_v_CM <- NULL
# also do not care about VM vs IM right now...only want sequential pairwise
DEData$result_VM_v_IM <- NULL

# filter: keep only significant results
mydataSig<-lapply(DEData,SigDEdf,PvaluesCol=7,CritP=0.05)


# different list of df for up or down expression
dataSigup<-lapply(mydataSig, MoreCritNum,column=3, critNum=0)
dataSigdown<-lapply(mydataSig, LessCritNum,column=3, critNum=0)



# this will be the background (all genes expressed) that the enrichment is tested with
# it is the same for all comparisons, so it doesn't matter if I chose 10v20 or 20v30, etc.
gene_universe <- DEData$result_VM_v_TM$Gene




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
                                  output_dir = "lettuce/go_results/multi")
                                 {
  
  direction <- match.arg(direction)
  ontology  <- match.arg(ontology)
  

  dataSig <- if (direction == "up") dataSigup else dataSigdown

  gene_list <- list(
    "VMvTM" = dataSig$result_VM_v_TM$Gene,
    "TMvIM" = dataSig$result_TM_v_IM$Gene,
    "IMvIMFM" = dataSig$result_IM_v_IMFM$Gene)
  
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
    qvalueCutoff = 0.05)
  
 
  go_df <- as.data.frame(go_res@compareClusterResult)
  
  result_list <- list(go_df[go_df$Cluster == "VMvTM", ],go_df[go_df$Cluster == "TMvIM", ],go_df[go_df$Cluster == "IMvIMFM", ])

  # rename the lists because these will be the sheets in the spreadsheet
  names(result_list) <- c(
      paste0(direction, "_VMvTM_", ontology),
      paste0(direction, "_TMvIM_", ontology),
      paste0(direction, "_IMvIMFM_", ontology))
  
  write_xlsx(result_list,file.path(output_dir,paste0("lsat_go_", direction, "_", ontology, ".xlsx")))
  
  
  return(go_res)
}


### function that will simplify compareCluster result ###

simplify_go_compareCluster <- function(go_object,
                                       direction,
                                       ontology,
                                       output_dir = "sunflower/go_results/multi",
                                       cutoff = 0.8) {
  
  go_simp <- simplify(go_object,
                      cutoff = cutoff,
                      by = "p.adjust",
                      select_fun = min)
  
  go_simp_df <- as.data.frame(go_simp@compareClusterResult)
  
  result_list_simp <- list(go_simp_df[go_simp_df$Cluster == "VMvTM", ],go_simp_df[go_simp_df$Cluster == "TMvIM", ],go_simp_df[go_simp_df$Cluster == "IMvIMFM", ])
  
  names(result_list_simp) <- c(
    paste0(direction, "_VMvTM_", ontology),
    paste0(direction, "_TMvIM_", ontology),
    paste0(direction, "_IMvIMFM_", ontology))
  
  write_xlsx(result_list_simp,
             file.path(output_dir,
                       paste0("go_", direction, "_", ontology, "_simp.xlsx")))
  
  
  return(go_simp)
}


#################
##### BP ########
#################

# run BP enrichment on upregulated genes
go_up_bp <- run_go_compareCluster("up", "BP",dataSigup,dataSigdown,gene_universe,org.Lsativa.eg.db)



# plot the data
options(enrichplot.colours = c("hotpink","seagreen")) 
plot <- dotplot(go_up_bp, showCategory = 15) +
  labs(title = "up-regulated enrichment analysis (BP)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("lettuce/go_results/multi/plots/lsat_up_regulated_BP_enrichment_15.png", plot = plot, height = 15, width = 10, dpi = 300)



# plot the data
options(enrichplot.colours = c("hotpink","seagreen")) 
plot <- dotplot(go_up_bp, showCategory = 10) +
  labs(title = "up-regulated enrichment analysis (BP)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("lettuce/go_results/multi/plots/lsat_up_regulated_BP_enrichment_10.png", plot = plot, height = 15, width = 10, dpi = 300)





# down
go_down_bp <- run_go_compareCluster("down", "BP",dataSigup,dataSigdown,gene_universe,org.Lsativa.eg.db)


# plot the data
options(enrichplot.colours = c("firebrick2","darkgoldenrod2")) 
plot <- dotplot(go_down_bp, showCategory = 15) +
  labs(title = "down-regulated enrichment analysis (BP)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("lettuce/go_results/multi/plots/lsat_down_regulated_BP_enrichment_15.png", plot = plot, height = 15, width = 10, dpi = 300)


options(enrichplot.colours = c("firebrick2","darkgoldenrod2")) 
plot <- dotplot(go_down_bp, showCategory = 10) +
  labs(title = "down-regulated enrichment analysis (BP)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("lettuce/go_results/multi/plots/lsat_down_regulated_BP_enrichment_10.png", plot = plot, height = 15, width = 10, dpi = 300)


#################
##### MF ########
#################

go_up_mf <- run_go_compareCluster("up", "MF", dataSigup, dataSigdown, gene_universe, org.Lsativa.eg.db)

# plot the data
options(enrichplot.colours = c("hotpink", "seagreen")) 
plot <- dotplot(go_up_mf, showCategory = 15) +
  labs(title = "up-regulated enrichment analysis (MF)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("lettuce/go_results/multi/plots/lsat_up_regulated_MF_enrichment_15.png", plot = plot, height = 15, width = 10, dpi = 300)


options(enrichplot.colours = c("hotpink", "seagreen")) 
plot <- dotplot(go_up_mf, showCategory = 15) +
  labs(title = "up-regulated enrichment analysis (MF)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("lettuce/go_results/multi/plots/lsat_up_regulated_MF_enrichment_10.png", plot = plot, height = 15, width = 10, dpi = 300)


# down

go_down_mf <- run_go_compareCluster("down", "MF", dataSigup, dataSigdown, gene_universe, org.Lsativa.eg.db)

# plot the data
options(enrichplot.colours = c("firebrick2", "darkgoldenrod2")) 
plot <- dotplot(go_down_mf, showCategory = 15) +
  labs(title = "down-regulated enrichment analysis (MF)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("lettuce/go_results/multi/plots/lsat_down_regulated_MF_enrichment_15.png", plot = plot, height = 15, width = 10, dpi = 300)

options(enrichplot.colours = c("firebrick2", "darkgoldenrod2")) 
plot <- dotplot(go_down_mf, showCategory = 10) +
  labs(title = "down-regulated enrichment analysis (MF)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("lettuce/go_results/multi/plots/lsat_down_regulated_MF_enrichment_10.png", plot = plot, height = 15, width = 10, dpi = 300)


#################
##### CC ########
#################

go_up_cc <- run_go_compareCluster("up", "CC", dataSigup, dataSigdown, gene_universe, org.Lsativa.eg.db)


# plot the data
options(enrichplot.colours = c("hotpink", "seagreen")) 
plot <- dotplot(go_up_cc, showCategory = 15) +
  labs(title = "up-regulated enrichment analysis (CC)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("lettuce/go_results/multi/plots/lsat_up_regulated_CC_enrichment_15.png", plot = plot, height = 15, width = 10, dpi = 300)


options(enrichplot.colours = c("hotpink", "seagreen")) 
plot <- dotplot(go_up_cc, showCategory = 10) +
  labs(title = "up-regulated enrichment analysis (CC)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("lettuce/go_results/multi/plots/lsat_up_regulated_CC_enrichment_10.png", plot = plot, height = 15, width = 10, dpi = 300)

# down

go_down_cc <- run_go_compareCluster("down", "CC", dataSigup, dataSigdown, gene_universe, org.Lsativa.eg.db)

# plot the data
options(enrichplot.colours = c("firebrick2", "darkgoldenrod2")) 
plot <- dotplot(go_down_cc, showCategory = 15) +
  labs(title = "down-regulated enrichment analysis (CC)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("lettuce/go_results/multi/plots/lsat_down_regulated_CC_enrichment_15.png", plot = plot, height = 15, width = 10, dpi = 300)

options(enrichplot.colours = c("firebrick2", "darkgoldenrod2")) 
plot <- dotplot(go_down_cc, showCategory = 10) +
  labs(title = "down-regulated enrichment analysis (CC)", x = "dev. stage") +
  theme_minimal()

print(plot)

ggsave("lettuce/go_results/multi/plots/lsat_down_regulated_CC_enrichment_10.png", plot = plot, height = 15, width = 10, dpi = 300)


