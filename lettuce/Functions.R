#########################
##### DE GENE SETS ######
#########################

NumSigDE <- function(vector) {
  sig <- which(vector < 0.05)
  return (length(sig))
}

SigDE <- function(vector) {
  sig <- which(vector < 0.05)
  return (sig)
}

LessCritNum <- function(dataset, column, critNum) {
  less <- which(dataset[,column] < critNum)
  data_less <- dataset[less,]
  return (data_less)
}

MoreCritNum <- function(dataset, column, critNum) {
  more <- which(dataset[,column] > critNum)
  data_more <- dataset[more,]
  return (data_more)
}


VennNumbers <- function (dataset1, dataset2, dataset3) {
  name1 <- deparse(substitute(dataset1))
  name2 <- deparse(substitute(dataset2))
  name3 <- deparse(substitute(dataset3))
  sig1_2 <- intersect(dataset1, dataset2)
  sig2_3 <- intersect(dataset2, dataset3)
  sig1_3 <- intersect(dataset1, dataset3)
  sig1_2only <- setdiff(sig1_2, dataset3)
  sig2_3only <- setdiff(sig2_3, dataset1)
  sig1_3only <- setdiff(sig1_3, dataset2)
  all <- intersect(sig1_2, dataset3)
  sig1only <- setdiff(dataset1, union(dataset2, dataset3))
  sig2only <- setdiff(dataset2, union(dataset1, dataset3))
  sig3only <- setdiff(dataset3, union(dataset1, dataset2))
  #print(paste0("Common in all three = ", length(all)))
  #print(paste0("In ", name1, " and ", name2, " only = ", length(sig1_2only)))
  #print(paste0("In ", name2, " and ", name2, " only = ", length(sig2_3only)))
  #print(paste0("In ", name1, " and ", name3, " only = ", length(sig1_3only)))
  #print(paste0(name1, " only = ", length(sig1only)))
  #print(paste0(name2, " only = ", length(sig2only)))
  #print(paste0(name3, " only = ", length(sig3only)))
  my_list <- data.frame("Category" = c("InCommonAll", paste0(name1, name2, "Only"), paste0(name2, name3, "Only"), paste0(name1, name3, "Only"), paste0(name1, "Only"), paste0(name2, "Only"), paste0(name3, "Only")), 
                        "Number" = c(length(all), length(sig1_2only), length(sig2_3only), length(sig1_3only), length(sig1only), length(sig2only), length(sig3only)))
  #my_list <- data.frame("Category" = c("InCommonAll", paste0(names[1], names[2], "Only"), paste0(names[2], names[3], "Only"), paste0(names[1], names[3], "Only"), paste0(names[1], "Only"), paste0(names[2], "Only"), paste0(names[3], "Only")), 
  #"Number" = c(length(all), length(sig1_2only), length(sig2_3only), length(sig1_3only), length(sig1only), length(sig2only), length(sig3only)))
  return (my_list)
}

SigDEdf <- function(dataset, PvaluesCol, CritP) {
  sig <- which(dataset[,PvaluesCol] < CritP)
  Sigdf <- dataset[sig,]
  return (Sigdf)
}

GeneSets <- function (dataset1, dataset2, dataset3) {
  name1 <- deparse(substitute(dataset1))
  name2 <- deparse(substitute(dataset2))
  name3 <- deparse(substitute(dataset3))
  sig1_2 <- intersect(dataset1, dataset2)
  sig2_3 <- intersect(dataset2, dataset3)
  sig1_3 <- intersect(dataset1, dataset3)
  sig1_2only <- setdiff(sig1_2, dataset3)
  #name(sig1_2only) <- paste0(name1, name2, "Only")
  sig2_3only <- setdiff(sig2_3, dataset1)
  #name(sig2_3only) <- paste0(name2, name3, "Only")
  sig1_3only <- setdiff(sig1_3, dataset2)
  #name(sig1_3only) <- paste0(name1, name3, "Only") 
  InCommonAll <- intersect(sig1_2, dataset3)
  sig1only <- setdiff(dataset1, union(dataset2, dataset3))
  #name(sig1only) <- paste0(name1, "Only")
  sig2only <- setdiff(dataset2, union(dataset1, dataset3))
  #name(sig2only) <- paste0(name2, "Only")
  sig3only <- setdiff(dataset3, union(dataset1, dataset2))
  #name(sig3only) <- paste0(name3, "Only")
  Indices <- list(InCommonAll, sig1_2only, sig2_3only, sig1_3only, sig1only, sig2only, sig3only)
  names(Indices) <- c("InCommonAll", paste0(name1, name2, "Only"), paste0(name2, name3, "Only"), paste0(name1, name3, "Only"), paste0(name1, "Only"), paste0(name2, "Only"), paste0(name3, "Only"))
  return (Indices)
}

GeneSets_only2 <- function (dataset1, dataset2) {
  name1 <- deparse(substitute(dataset1))
  name2 <- deparse(substitute(dataset2))
  sig1_and2 <-intersect(dataset1, dataset2)
  sig1only <- setdiff(dataset1,(dataset2))
  #name(sig1only) <- paste0(name1, "Only")
  sig2only <- setdiff(dataset2,(dataset1))
  #name(sig3only) <- paste0(name3, "Only")
  Indices <- list(sig1_and2,sig1only, sig2only)
  names(Indices) <- c(paste0(name1, name2, "Only"), paste0(name1, "Only"), paste0(name2, "Only"))
  return (Indices)
}

GeneSets_four <- function (dataset1, dataset2, dataset3, dataset4) {
  name1 <- deparse(substitute(dataset1))
  name2 <- deparse(substitute(dataset2))
  name3 <- deparse(substitute(dataset3))
  name4 <- deparse(substitute(dataset4))
  
  sig1_2 <- intersect(dataset1, dataset2)
  sig2_3 <- intersect(dataset2, dataset3)
  sig1_3 <- intersect(dataset1, dataset3)
  sig3_4 <- intersect(dataset3, dataset4)
  sig1_4 <- intersect(dataset1, dataset4)
  sig2_4 <- intersect(dataset2, dataset4)
  
  sig1_2only <- setdiff(sig1_2, dataset3)
  #name(sig1_2only) <- paste0(name1, name2, "Only")
  sig2_3only <- setdiff(sig2_3, dataset1)
  #name(sig2_3only) <- paste0(name2, name3, "Only")
  sig1_3only <- setdiff(sig1_3, dataset2)
  #name(sig1_3only) <- paste0(name1, name3, "Only") 
  sig3_4only <- setdiff(sig3_4, union(dataset1, dataset2))
  sig1_4only <- setdiff(sig1_4, union(dataset2, dataset3))
  sig2_4only <- setdiff(sig2_4, union(dataset1, dataset3))
  
  
  
  
  
  InCommonAll <- intersect(sig1_2, sig3_4)
  sig1only <- setdiff(dataset1, union(dataset2,union(dataset3, dataset4)))
  #name(sig1only) <- paste0(name1, "Only")
  sig2only <- setdiff(dataset2, union(dataset1,union(dataset3, dataset4)))
  #name(sig2only) <- paste0(name2, "Only")
  sig3only <- setdiff(dataset3, union(dataset1,union(dataset2, dataset4)))
  #name(sig3only) <- paste0(name3, "Only")
  sig4only <- setdiff(dataset4, union(dataset1,union(dataset2, dataset3)))
  
  
  
  Indices <- list(InCommonAll, sig1_2only, sig2_3only, sig1_3only, sig3_4only, sig1_4only, sig2_4only, sig1only, sig2only, sig3only, sig4only)
  names(Indices) <- c("InCommonAll", paste0(name1, name2, "Only"), paste0(name2, name3, "Only"), paste0(name1, name3, "Only"), paste0(name3, name4, "Only"),paste0(name1, name4, "Only"), paste0(name2, name4, "Only"),
                      paste0(name1, "Only"), paste0(name2, "Only"), paste0(name3, "Only"), paste0(name4, "Only"))
  return (Indices)
}





GeneSets_five <- function(dataset1, dataset2, dataset3, dataset4, dataset5) {
  name1 <- deparse(substitute(dataset1))
  name2 <- deparse(substitute(dataset2))
  name3 <- deparse(substitute(dataset3))
  name4 <- deparse(substitute(dataset4))
  name5 <- deparse(substitute(dataset5))
  
  # Pairwise intersections
  sig1_2 <- intersect(dataset1, dataset2)
  sig2_3 <- intersect(dataset2, dataset3)
  sig3_4 <- intersect(dataset3, dataset4)
  sig4_5 <- intersect(dataset4, dataset5)
  sig1_3 <- intersect(dataset1, dataset3)
  sig1_4 <- intersect(dataset1, dataset4)
  sig1_5 <- intersect(dataset1, dataset5)
  sig2_4 <- intersect(dataset2, dataset4)
  sig2_5 <- intersect(dataset2, dataset5)
  sig3_5 <- intersect(dataset3, dataset5)
  
  # Pairwise-only
  sig1_2only <- setdiff(sig1_2, union(dataset3, union(dataset4, dataset5)))
  sig2_3only <- setdiff(sig2_3, union(dataset1, union(dataset4, dataset5)))
  sig3_4only <- setdiff(sig3_4, union(dataset1, union(dataset2, dataset5)))
  sig4_5only <- setdiff(sig4_5, union(dataset1, union(dataset2, dataset3)))
  sig1_3only <- setdiff(sig1_3, union(dataset2, union(dataset4, dataset5)))
  sig1_4only <- setdiff(sig1_4, union(dataset2, union(dataset3, dataset5)))
  sig1_5only <- setdiff(sig1_5, union(dataset2, union(dataset3, dataset4)))
  sig2_4only <- setdiff(sig2_4, union(dataset1, union(dataset3, dataset5)))
  sig2_5only <- setdiff(sig2_5, union(dataset1, union(dataset3, dataset4)))
  sig3_5only <- setdiff(sig3_5, union(dataset1, union(dataset2, dataset4)))
  
  # In common across all five datasets
  InCommonAll <- Reduce(intersect, list(dataset1, dataset2, dataset3, dataset4, dataset5))
  
  # Unique to each dataset
  sig1only <- setdiff(dataset1, union(dataset2, union(dataset3, union(dataset4, dataset5))))
  sig2only <- setdiff(dataset2, union(dataset1, union(dataset3, union(dataset4, dataset5))))
  sig3only <- setdiff(dataset3, union(dataset1, union(dataset2, union(dataset4, dataset5))))
  sig4only <- setdiff(dataset4, union(dataset1, union(dataset2, union(dataset3, dataset5))))
  sig5only <- setdiff(dataset5, union(dataset1, union(dataset2, union(dataset3, dataset4))))
  
  # Collect results
  Indices <- list(
    InCommonAll,
    sig1_2only, sig2_3only, sig3_4only, sig4_5only,
    sig1_3only, sig1_4only, sig1_5only,
    sig2_4only, sig2_5only, sig3_5only,
    sig1only, sig2only, sig3only, sig4only, sig5only
  )
  
  names(Indices) <- c(
    "InCommonAll",
    paste0(name1, name2, "Only"),
    paste0(name2, name3, "Only"),
    paste0(name3, name4, "Only"),
    paste0(name4, name5, "Only"),
    paste0(name1, name3, "Only"),
    paste0(name1, name4, "Only"),
    paste0(name1, name5, "Only"),
    paste0(name2, name4, "Only"),
    paste0(name2, name5, "Only"),
    paste0(name3, name5, "Only"),
    paste0(name1, "Only"),
    paste0(name2, "Only"),
    paste0(name3, "Only"),
    paste0(name4, "Only"),
    paste0(name5, "Only")
  )
  
  return(Indices)
}


## Import files

ImportCSVs <- function (DirPath, critP) {
  my_files <- list.files(path = DirPath, pattern = "*.csv", full.names = TRUE)
  my_data <- lapply(my_files, read.csv)
  names(my_data) <- gsub("\\.csv$", "", my_files)
  names(my_data) <- gsub(DirPath, "", names(my_data))
  names(my_data) <- gsub("/", "", names(my_data))
  for (i in seq_along(my_data)) {
    name <- names(my_data[i])
    colnames(my_data[[i]])[1] <- "Gene"
    my_data[[i]] <- my_data[[i]][order(my_data[[i]]$Gene),]
  }
  return(my_data)
}

PairwiseCats <- function (Overlap, Diff, BothUp, BothDown) {
  Unconditional <- setdiff(Overlap, Diff)
  Up_Unconditional <- intersect(Unconditional, BothUp)
  Down_Unconditional <- intersect(Unconditional, BothDown)
  Conditional <- intersect(Overlap, Diff)
  Up_Conditional <- intersect(Conditional, BothUp)
  Down_Conditional <- intersect(Conditional, BothDown)
  DiffDir_Conditional <- setdiff(Conditional, union(Up_Conditional, Down_Conditional))
  my_list <- list(Unconditional_Up = Up_Unconditional,
                  Unconditional_Down = Down_Unconditional,
                  Conditional_Up = Up_Conditional,
                  Conditional_Down = Down_Conditional,
                  Conditional_Diff_Dir = DiffDir_Conditional)
  return(my_list)
}

Prioritized <- function (Diff, OverlapUp, BothUp_Col, OverlapDown, BothDown_Col) {
  Prioritized_Up <- lapply(OverlapUp, function(x) {intersect(x$Gene, Diff)})
  Prioritized_Down <- lapply(OverlapDown, function(x) {intersect(x$Gene, Diff)})
  my_list <- list(All_Up = Prioritized_Up[[1]],
                  Both_Up = Prioritized_Up[[BothUp_Col]],
                  All_Down = Prioritized_Down[[1]],
                  Both_Down = Prioritized_Down[[BothDown_Col]])
  return(my_list)
}

DirectionDf <- function(dataset, Yvar, Xvar) {
  Dir1 <- dataset[which(dataset[,Xvar] > 0 &
                          dataset[,Yvar] < 0),]
  Dir2 <- dataset[which(dataset[,Xvar] < 0 &
                          dataset[,Yvar] > 0),]
  DirGenes <- union(Dir1$Gene, Dir2$Gene)
  ReducedGenes <- setdiff(dataset[which(abs(dataset[,Xvar]) >
                                          abs(dataset[,Yvar])),"Gene"], 
                          DirGenes)
  ReducedData <- dataset[which(dataset$Gene %in% ReducedGenes),]
  IncreasedGenes <- setdiff(dataset[which(abs(dataset[,Xvar]) <
                                            abs(dataset[,Yvar])),"Gene"], 
                            DirGenes)
  IncreasedData <- dataset[which(dataset$Gene %in% IncreasedGenes),]
  Yname <- deparse(substitute(Yvar))
  Xname <- deparse(substitute(Xvar))
  MyList <- list("Down" = Dir1,
                 "Up" = Dir2,
                 "Reduced" = ReducedData,
                 "Increased" = IncreasedData)
  return (MyList)
}

# a function to return lists of shared genes that are significantly different (and aren't)
Compare_shared <- function(SharedOverlap, PairwiseDiffs) {
  Shared_same <- setdiff(SharedOverlap, PairwiseDiffs)
  Shared_diff <- intersect(SharedOverlap, PairwiseDiffs)
  Shared_lists <- list(Shared_same, Shared_diff)
  names(Shared_lists) <- c("Shared_same", "Shared_different")
  return (Shared_lists)
}

# function to compare direction in 2 columns:

Treat_compare <- function(dataframe, col1, col2, col1_name, col2_name) {
  col1_Up <- MoreCritNum(dataframe, col1, 0)
  col1_Down <- LessCritNum(dataframe, col1, 0)
  col1_Up_col2Down <- LessCritNum(col1_Up, col2, 0)
  col1_Down_col2Up <- MoreCritNum(col1_Down, col2, 0)
  SameDir <- dataframe[which(!dataframe$Gene %in% union(col1_Up_col2Down$Gene,
                                                        col1_Down_col2Up$Gene)),]
  col2_Red <- SameDir[which(abs(SameDir[,col1]) >
                              abs(SameDir[,col2])),]
  col2_Inc <- SameDir[which(abs(SameDir[,col1]) <
                              abs(SameDir[,col2])),]
  col1_Up_col2Red <- MoreCritNum(col2_Red, col1, 0)
  col1_Down_col2Red <- LessCritNum(col2_Red, col1, 0)
  col1_Up_col2Inc <- MoreCritNum(col2_Inc, col1, 0)
  col1_Down_col2Inc <- LessCritNum(col2_Inc, col1, 0)
  Magnitudes <- list(col1_Up_col2Red, col1_Down_col2Red, col1_Up_col2Inc, col1_Down_col2Inc, col1_Up_col2Down, col1_Down_col2Up)
  names(Magnitudes) <- c(paste0("Up_", col1_name, "_reduced_", col2_name), paste0("Down_", col1_name, "_reduced_", col2_name),
                         paste0("Up_", col1_name, "_increased_", col2_name), paste0("Down_", col1_name, "_increased_", col2_name),
                         paste0("Up_", col1_name, "_OppositeDir_", col2_name), paste0("Down_", col1_name, "_OppositeDir_", col2_name))
  return (Magnitudes)
}
Treatment_compare <- function(Sig_df1, Sig_df2, Diff_df, Name_treat1, Name_treat2) {
  SigSet1_only <- setdiff(Sig_df1$Gene, Sig_df2$Gene)
  SharedOverlap <- intersect(Sig_df1$Gene, Sig_df2$Gene)
  Shared_same <- setdiff(SharedOverlap, Diff_df$Gene)
  Shared_diff <- intersect(SharedOverlap, Diff_df$Gene)
  Shared_df <- merge(Sig_df1, Sig_df2, by="Gene")
  Shared_diff_df <- Shared_df[which(Shared_df$Gene %in% Shared_diff),]
  treat1_Up <- MoreCritNum(Shared_diff_df, "log2FoldChange.x", 0)
  treat1_Down <- LessCritNum(Shared_diff_df, "log2FoldChange.x", 0)
  treat1Up_treat2Down <- LessCritNum(treat1_Up, "log2FoldChange.y", 0)
  treat1Down_treat2Up <- MoreCritNum(treat1_Down, "log2FoldChange.y", 0)
  DiffDir <- union(treat1Up_treat2Down$Gene, treat1Down_treat2Up$Gene)
  SameDir <- Shared_diff_df[which(!Shared_diff_df$Gene %in% DiffDir),]
  treat2_Red <- SameDir[which(abs(SameDir[,"log2FoldChange.x"]) >
                                abs(SameDir[,"log2FoldChange.y"])),]
  treat2_Inc <- SameDir[which(abs(SameDir[,"log2FoldChange.x"]) <
                                abs(SameDir[,"log2FoldChange.y"])),]
  Categories <- list(SigSet1_only, Shared_same, treat2_Red$Gene, treat2_Inc$Gene, DiffDir)
  names(Categories) <- c(paste0(Name_treat1, "_only"), 
                         paste0(Name_treat2, "_", Name_treat1, "_Same"), 
                         paste0(Name_treat2, "_reduced"),
                         paste0(Name_treat2, "_increased"),
                         paste0(Name_treat2,"_Diff_direction"))
  return (Categories)
}

df_from_List <- function(list, labels, list_name) {
  df <- data.frame(unlist(list, recursive = TRUE, use.names = TRUE))
  colnames(df) <- c("Number")
  df$Stress <- c(list_name)
  df_labels <- cbind(df, labels)
  df_total <- sum(df_labels$Number)
  df_labels$Prop <- df_labels$Number / df_total
  return (df_labels)
}

#########################
##### REGRESSION #####
#########################

Predictdf <- function(dataset, geneset, Yvar, Xvar) {
  dfSubset <- dataset[which(dataset$Gene %in% geneset),]
  Mod <- lm(dfSubset[,Yvar] ~ dfSubset[,Xvar])
  sum<-summary(Mod)
  residSE <- sum$sigma
  dfPredict <- as.data.frame(predict(Mod, level=0.95, interval= 'prediction'))
  dfConfid <- as.data.frame(predict(Mod, level=0.95, interval= 'confidence'))
  colnames(dfConfid) <- c("CI_fit", "CI_lower", "CI_upper")
  resid <- resid(Mod)
  dfSubsetPredict <- cbind(dfSubset, resid, dfPredict, dfConfid)
  dfSubsetPredict$Interval <- ifelse(dfSubsetPredict[,Yvar] < 
                                       dfSubsetPredict$lwr |
                                       dfSubsetPredict[,Yvar] > 
                                       dfSubsetPredict$upr, "Outside", "Within")
  dfSubsetPredict$ResidSE_Interval <- ifelse(abs(dfSubsetPredict$resid) > 
                                               residSE, "Outside", "Within")
  dfSubsetPredict$CI <- ifelse(dfSubsetPredict[,Yvar] < 
                                       dfSubsetPredict$CI_lower |
                                       dfSubsetPredict[,Yvar] > 
                                       dfSubsetPredict$CI_upper, "Outside", "Within")
  return (dfSubsetPredict)
}

#########################
##### GO ENRICHMENT #####
#########################

### Function for GO Term Enrichment

#LengthTable must have Gene ID as first column + Length as 2nd column
#GO_Terms must have GeneID ("Parent") as first column, and GO Term ("Ontology_term") in 2nd column

GO_Enrichment <- function (AllGenes, DE_Genes, LengthTable, GO_Terms) {
  common <- intersect(AllGenes, LengthTable[,1])
  LengthTable <- as.data.frame(LengthTable[LengthTable[,1] %in% common, ]) #"all" genes
  LengthTable <- LengthTable[order(LengthTable[,1]),] #order by gene ID
  gene.vector=as.integer(AllGenes %in% DE_Genes)
  names(gene.vector)=AllGenes
  pwf <- nullp(gene.vector, bias.data = LengthTable[,2])
  GO.wall <- goseq(pwf, gene2cat = GO_Terms)
  enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05]
  enrichedGO_table <- GO.wall[is.element(GO.wall$category, enriched.GO),]
  return(enrichedGO_table)
}

#########################
##### COEXPRESSION ######
#########################

# function to vary "deepSplit" and "pamStage" paramter:

Module_ID <- function(dataframe, GeneTree, TOMdiss, split, pam, minModuleSize) {
dynamicMods = cutreeDynamic(dendro = GeneTree, distM = TOMdiss,
              method = "hybrid",
              #cutHeight = , default for "hybrid" is 99% of range between 5th percentile and joining heights of dendrogram
                            deepSplit = split,  #from 0-4
                            pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize,
                            pamStage = pam)
Mods_Genes <- as.data.frame(table(dynamicMods))
UnAssigned <- ifelse(Mods_Genes[1,1] == 0, Mods_Genes[1,2], "NA")
Mods_Genes_noZero <- subset(Mods_Genes, Mods_Genes[,1]!=0)
dynamicColors = labels2colors(dynamicMods)
MEList = moduleEigengenes(dataframe, colors = dynamicColors)
datVA <- MEList$varExplained
df <- data.frame(matrix(unlist(datVA ), nrow=length(datVA), byrow=T))
df_noMin <- subset(df, df[,1] != min(df[,1]))  # for values that don't include "Grey" modules
ValueList <- data.frame("NumMods" = length(Mods_Genes_noZero$Freq), "NumUnAssigned" = UnAssigned, 
  "MeanGeneNum" = mean(Mods_Genes_noZero[,2]), "MinGeneNum" = min(Mods_Genes_noZero[,2]),  
  "MaxGeneNum" = max(Mods_Genes_noZero[,2]), "MeanVarExp" = mean(df[,1]), "MeanVarExp_noMin" = mean(df_noMin[,1]),
  "MinVarExp" = min(df[,1]), "MinVarExp_noMin" = min(df_noMin[,1]),  "MaxVarExp" = max(df[,1]))
return(ValueList)
}

Module_Cluster <- function(MEDissThres, datExpr, dynamicColors) {
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  mergedMEs = merge$newMEs
  datVA <- moduleEigengenes(datExpr, merge$colors)$varExplained
  df <- data.frame(matrix(unlist(datVA ), nrow=length(datVA), byrow=T))
  ModGeneNums <- as.data.frame(table(merge$colors))
  Metric <- c("NumMods", "MeanGeneNum", "MinGeneNum", "MaxGeneNum",
              "MeanVarExp", "MinVarExp", "MaxVarExp")
  Value <- c(length(mergedMEs), mean(ModGeneNums$Freq), min(ModGeneNums$Freq), max(ModGeneNums$Freq),
             mean(df[,1]), min(df[,1]), max(df[,1]))
  Result <- data.frame(Metric, Value)
  return(Result)
}

# function to vary "deepSplit" and "pamStage" paramters for a consensus tree:
Module_ID_multiExpr <- function(dataframe, GeneTree, TOMdiss, split, pam, minModuleSize) {
dynamicMods = cutreeDynamic(dendro = GeneTree, distM = TOMdiss,
              method = "hybrid",
              #cutHeight = , default for "hybrid" is 99% of range between 5th percentile and joining heights of dendrogram
                            deepSplit = split,  #from 0-4
                            pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize,
                            pamStage = pam)
Mods_Genes <- as.data.frame(table(dynamicMods))
UnAssigned <- ifelse(Mods_Genes[1,1] == 0, Mods_Genes[1,2], "NA")
Mods_Genes_noZero <- subset(Mods_Genes, Mods_Genes[,1]!=0)
dynamicColors = labels2colors(dynamicMods)
MEList = multiSetMEs(dataframe, colors = NULL, universalColors = dynamicColors)
datVA1 <- MEList[[1]]$varExplained
df1 <- data.frame(matrix(unlist(datVA1), nrow=length(datVA1), byrow=T))
datVA2 <- MEList[[2]]$varExplained
df2 <- data.frame(matrix(unlist(datVA2), nrow=length(datVA2), byrow=T))
df_all <- cbind(df1, df2)
df_all$AveVarExp <- (df_all[,1] + df_all[,2]) / 2
df_noMin <- subset(df_all, df_all$AveVarExp != min(df_all$AveVarExp)) # for values that don't include "Grey" modules
ValueList <- data.frame("NumMods" = length(Mods_Genes_noZero$Freq), "NumUnAssigned" = UnAssigned, 
  "MeanGeneNum" = mean(Mods_Genes_noZero[,2]), "MinGeneNum" = min(Mods_Genes_noZero[,2]),  
  "MaxGeneNum" = max(Mods_Genes_noZero[,2]), "MeanVarExp" = mean(df_all$AveVarExp), "MeanVarExp_noMin" = mean(df_noMin$AveVarExp),
  "MinVarExp" = min(df_all$AveVarExp), "MinVarExp_noMin" = min(df_noMin$AveVarExp),  "MaxVarExp" = max(df_all$AveVarExp))
return(ValueList)
}

Module_Cluster_multiExpr <- function(MEDissThres, multiExpr, unmergedLabels) {
  merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = MEDissThres, verbose = 3)
  datVA1 <- merge$newMEs[[1]]$varExplained
  df1 <- data.frame(matrix(unlist(datVA1), nrow=length(datVA1), byrow=T))
  datVA2 <- merge$newMEs[[2]]$varExplained
  df2 <- data.frame(matrix(unlist(datVA2), nrow=length(datVA2), byrow=T))
  df_all <- cbind(df1, df2)
  df_all$AveVarExp <- (df_all[,1] + df_all[,2]) / 2
  Metric <- c("NumMods", "MeanVarExp", "MinVarExp", "MaxVarExp")
  Value <- c(length(df_all$AveVarExp), mean(df_all$AveVarExp), min(df_all$AveVarExp), max(df_all$AveVarExp))
  Result <- data.frame(Metric, Value)
  return(Result)
}

# Function to look at Module-Trait Correlations
Module_Trait_Corr <- function(TraitData, model, modEigenvalues, NumSamples) {
  datTraits <- as.data.frame(model.matrix(model, data=TraitData))
  moduleTraitCor = cor(modEigenvalues, datTraits, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, NumSamples)
  result = list(Traits = datTraits, ModTraitCorr = moduleTraitCor, ModTraitP = moduleTraitPvalue)
  return(result)
}

# Function to plot raw data as a barplot
base_barplot <- function(var1, module_name){
  barplot(var1, main = module_name)
}

# Function to correct p-values using the Holm or Hochberg procedure
Holm_Hochberg <- function(pvals_column, df){
  pOrder <- df[order(df[,pvals_column], decreasing = FALSE),
               c("Module", pvals_column)]
  colnames(pOrder) <- c("Module", "Pvals")
  pOrder$Order <- seq(1,length(pOrder$Module), by=1)
  N <- length(pOrder$Module)
  pOrder$alpha <- 0.05 / (N - pOrder$Order + 1)
  pOrder$Sig <- ifelse(pOrder$Pvals < pOrder$alpha,
                         "Significant", "N.S.")
  return(pOrder)
}

# summarize connectivity for all genes in each module
ConnectivitySummary <- function(Yvar, Group, YvarLabel) {
  NumGenes <- aggregate(Yvar, by=list(Group), length)
  mean <- aggregate(Yvar, by=list(Group), mean)
  sd <- aggregate(Yvar, by=list(Group), sd)
  SummaryStats <- merge(NumGenes, 
                        merge(mean, sd, by="Group.1"), 
                        by="Group.1")
  colnames(SummaryStats) <- c("ModuleColor", "NumGenes", paste0(YvarLabel,"_Mean"), paste0(YvarLabel,"_SD"))
  SummaryStats[,paste0(YvarLabel, "_SE")] <- SummaryStats[,paste0(YvarLabel,"_SD")] / (sqrt(SummaryStats$NumGenes))
  SummaryStats$Module <- paste0("ME", SummaryStats$ModuleColor)
  return(SummaryStats)
}

#########################
##### GO ENRICHMENT #####
#########################

### Function for GO Term Enrichment

#LengthTable must have Gene ID as first column + Length as 2nd column
#GO_Terms must have GeneID ("Parent") as first column, and GO Term ("Ontology_term") in 2nd column

GO_Enrichment <- function (AllGenes, DE_Genes, LengthTable, GO_Terms) {
  common <- intersect(AllGenes, LengthTable[,1])
  LengthTable <- as.data.frame(LengthTable[LengthTable[,1] %in% common, ]) #"all" genes
  LengthTable <- LengthTable[order(LengthTable[,1]),] #order by gene ID
  gene.vector=as.integer(AllGenes %in% DE_Genes)
  names(gene.vector)=AllGenes
  pwf <- nullp(gene.vector, bias.data = LengthTable[,2])
  GO.wall <- goseq(pwf, gene2cat = GO_Terms)
  enriched.GO=GO.wall$category[p.adjust(GO.wall$over_represented_pvalue, method="BH")<.05]
  enrichedGO_table <- GO.wall[is.element(GO.wall$category, enriched.GO),]
  return(enrichedGO_table)
}

#########################
### PLOTTING MODULES ####
#########################

# this function plots the overall lsmeans in each treatment for each module
Plot_ModOverallMeans <- function(ModuleName, dataframe, color, title){
  means <- as.data.frame(t(dataframe[which(dataframe$Module==ModuleName),c("Module", 
                                                                           "emmean.Control", "emmean.LowNut",
                                                                           "emmean.Salt", "emmean.Combo")]))
  se <- as.data.frame(t(dataframe[which(dataframe$Module==ModuleName),c("Module", "SE.Control",
                                                                        "SE.LowNut", "SE.Salt", "SE.Combo")]))
  pvals <- as.data.frame(t(dataframe[which(dataframe$Module==ModuleName),c("Module", "Difference_p.DE_Nut",
                                                                           "Difference_p.DE_Salt","Difference_p.DE_Combo")]))
  Mod_df <- data.frame(Treatment=c("Control", "Nutrient", "Salt", "Combo"), 
                       DE=as.numeric(means[2:5,1]),
                       SE=as.numeric(se[2:5,1]),
                       pval=as.numeric(c(NA, pvals[2:4,1])))
  Mod_df$Treatment <- factor(Mod_df$Treatment, levels=c("Control", "Nutrient", "Salt", "Combo"))
  Mod_df$plabs <- ifelse(Mod_df$pval < 0.0001, "***",
                         ifelse(Mod_df$pval < 0.001, "**",
                                ifelse(Mod_df$pval < 0.01, "*",
                                       ifelse(Mod_df$pval == "NA", "NA",
                                              paste0("p=",round(Mod_df$pval, digits=2))))))
  scale <- abs(max(Mod_df$DE) - min(Mod_df$DE))
  Mod_df$plabPositions <- ifelse(Mod_df$pval < 0.01,
                                 Mod_df$DE + sign(Mod_df$DE) * Mod_df$SE +
                                   sign(Mod_df$DE) * (scale/8),              
                                 Mod_df$DE + sign(Mod_df$DE) * Mod_df$SE +
                                   sign(Mod_df$DE) * (scale/5))
  ymin <- min(Mod_df$DE) - (2*max(Mod_df$SE) + (scale/5))
  #ymin <- min(Mod_df$DE) - scale/5
  ymin_plot <- ifelse(ymin > 0, 0, ymin)
  ymax <- max(Mod_df$DE) + (2*max(Mod_df$SE) + (scale/5))
  #ymax <- max(Mod_df$DE) + scale/5
  ymax_plot <- ifelse(ymax < 0, 0, ymax)
  p <- ggplot(data = Mod_df, aes(x=Treatment, y=DE)) +
    geom_bar(stat="identity", fill=color, alpha=0.7, color="black") +
    geom_errorbar(aes(ymin=DE-SE, ymax=DE+SE), width=0.05) +
    #geom_linerange(aes(ymin=DE-SE, ymax=DE+SE)) +
    theme_minimal() +
    ggtitle(title) +
    geom_text(aes(label=plabs, y=plabPositions)) +
    #vjust="outward") +
    ylim(ymin_plot, ymax_plot)
  return(p)
}
