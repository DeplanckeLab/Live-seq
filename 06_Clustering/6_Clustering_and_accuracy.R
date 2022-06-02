################################################################
#                                                              #
#                     CLUSTERING ACCURACY                      #
#                                                              #
################################################################


### Author: Pernille Rainer, pernille.rainer@epfl.ch
### Date: 01-02.2022
### Datasets: scRNA-seq and Live-seq 
### Goal: Compute clustering accuracy (ARI, Barplot comparison and clustree)

library(rprojroot)
root_dir <- find_root(has_file("Live-seq.RProj"))
source(paste0(root_dir, "/utils.R"))
require(Seurat); require(clustree); require(ggplot2); require(gridExtra); require(dplyr)

##---------------------------------------------##
##------------------Load data------------------##
##---------------------------------------------##

f_subsample <- function(all, CT , B, Sampling, correct_batch = T ){
  myseu <- subset(all, subset = Cell_type %in% CT & (Batch %in% B)  & (Species != "Human") & (sampling_type == Sampling))
  if(length(B) > 1 & correct_batch){
    var_reg <- c("Batch", "nCount_RNA", "nFeature_RNA")}else{
      var_reg <- c("nCount_RNA", "nFeature_RNA")}
  myseu <- myseu %>% FindVariableFeatures() %>% NormalizeData() %>% ScaleData(vars.to.regress = var_reg) %>% RunPCA()
  myseu <- myseu %>% JackStraw(num.replicate = 100) %>% ScoreJackStraw(dims = 1:20)
  return(myseu)
}

## Seurat objects 
scRNAseq_only <- readRDS(paste0(root_dir, "/data/scRNAseq_only.rds"))
Liveseq_only <- readRDS(paste0(root_dir, "/data/Liveseq_only.rds"))

##---------------------------------------------##
##--------------Compute clustering-------------##
##---------------------------------------------##
# Cell types were clustered on the entire data set, ASPCs were clustered independently 

### General clustering:
Liveseq_only <- Liveseq_only %>% FindClusters(res = 0.2)
scRNAseq_only <- scRNAseq_only %>% FindClusters(res = 0.2)

### Clustering on ASPCs:
# scRNAseq:
## !! Do not correct for batch when clustering on ASPCs from all batches:
# one batch contains only not treated cells 
sc_aspcs_all <- f_subsample(scRNAseq_only, CT = "ASPC", B = c("5Reseq", "6",  "7", "scRNA","8_8", "9_1"), Sampling = "scRNA", correct_batch = F)
sc_aspcs_all <- sc_aspcs_all %>%  RunTSNE(dims = 1:10, perplexity = 10) %>% FindNeighbors(dims =  1:10) %>% FindClusters(res = 0.2)
levels(sc_aspcs_all$integrated_snn_res.0.2) <- c("ASPC_not_treated", "ASPC_DMIR_treated")

scRNAseq_only$MODIFIED_integrated_snn_res.0.2 <- scRNAseq_only$integrated_snn_res.0.2
levels(scRNAseq_only$MODIFIED_integrated_snn_res.0.2) <- c("Raw264.7_not_treated", "IBA_not_treated", "Raw264.7_LPS_treated", "ASPC")
scRNAseq_only$MODIFIED_integrated_snn_res.0.2 <- as.character(scRNAseq_only$MODIFIED_integrated_snn_res.0.2)
scRNAseq_only$MODIFIED_integrated_snn_res.0.2[colnames(sc_aspcs_all)] <- as.character(sc_aspcs_all$integrated_snn_res.0.2)

# Live-seq:
## Same as scRNA-seq
## !! Do not correct for batch when clustering on ASPCs from all batches:
# one batch contains only not treated cells 
live_aspcs_all <- f_subsample(Liveseq_only, CT = "ASPC", B = c("5Reseq", "6", "7", "8_8", "9_4"), Sampling = "Live_seq", correct_batch = F)
live_aspcs_all <- live_aspcs_all %>%  RunTSNE(dims = 1:10, perplexity = 10) %>% 
  FindNeighbors(dims =  1:10) %>% FindClusters(res = 0.6)
levels(live_aspcs_all$RNA_snn_res.0.6) <- c("ASPC_not_treated", "ASPC_DMIR_treated")

Liveseq_only$MODIFIED_RNA_snn_res.0.2 <- Liveseq_only$RNA_snn_res.0.2
levels(Liveseq_only$MODIFIED_RNA_snn_res.0.2) <- c("Raw264.7_not_treated", "ASPCs", "Raw264.7_LPS_treated", "IBA_not_treated")
Liveseq_only$MODIFIED_RNA_snn_res.0.2 <- as.character(Liveseq_only$MODIFIED_RNA_snn_res.0.2)
Liveseq_only$MODIFIED_RNA_snn_res.0.2[colnames(live_aspcs_all)] <- as.character(live_aspcs_all$RNA_snn_res.0.6)

saveRDS(scRNAseq_only, paste0(root_dir,"/data/scRNAseq_only.rds"))
saveRDS(Liveseq_only, paste0(root_dir,"/data/Liveseq_only.rds"))

##---------------------------------------------##
##--------Calculate Adjusted Rand Index--------##
##---------------------------------------------##

## Live-seq:
  # Extended data Figure 2 f (first two panels)
grid.arrange(DimPlot(Liveseq_only, group.by = "MODIFIED_RNA_snn_res.0.2", reduction = "tsne",
                     cols = unname(celltype_treatment_colors[levels(as.factor(Liveseq_only$MODIFIED_RNA_snn_res.0.2))])) + theme(legend.position = "none"),
             DimPlot(Liveseq_only, group.by = "celltype_treatment", reduction = "tsne",
                     cols = unname(c(celltype_treatment_colors[levels(as.factor(Liveseq_only$celltype_treatment))]))) + theme(legend.position = "none"),
             nrow = 1,
             top = paste0("ARI = ", round(mclust::adjustedRandIndex(Liveseq_only$celltype_treatment, Liveseq_only$MODIFIED_RNA_snn_res.0.2), 2)))
mclust::adjustedRandIndex(Liveseq_only$celltype_treatment, Liveseq_only$MODIFIED_RNA_snn_res.0.2) #0.8132462

## scRNA-seq:
  # Extended data Figure 2 o (first two panels)
grid.arrange(DimPlot(scRNAseq_only, group.by = "MODIFIED_integrated_snn_res.0.2", reduction = "tsne", 
                     cols = unname(celltype_treatment_colors[levels(as.factor(scRNAseq_only$MODIFIED_integrated_snn_res.0.2))])) + theme(legend.position = "none"),
             DimPlot(scRNAseq_only, group.by = "celltype_treatment", reduction = "tsne",
                     cols = unname(c(celltype_treatment_colors[levels(as.factor(scRNAseq_only$celltype_treatment))]))) + theme(legend.position = "none"),
             nrow = 1,
             top = paste0("ARI = ", round(mclust::adjustedRandIndex(scRNAseq_only$celltype_treatment, scRNAseq_only$MODIFIED_integrated_snn_res.0.2), 2)))
mclust::adjustedRandIndex(scRNAseq_only$celltype_treatment, scRNAseq_only$MODIFIED_integrated_snn_res.0.2) #0.9522529

##---------------------------------------------##
##-------------Barplot comparison--------------##
##---------------------------------------------##

## Barplot showing the number of cells per cell state/type in each cluster

## Live-seq:
  # Extended data figure 2 g
df <- reshape2::melt(table(Liveseq_only$celltype_treatment, paste0(Liveseq_only$MODIFIED_RNA_snn_res.0.2, "-clust")))
df$Var1 <- factor(as.character(df$Var1), levels = names(celltype_treatment_colors))
df$Var2 <- factor(as.character(df$Var2), levels = paste0(names(celltype_treatment_colors), "-clust"))
ggplot(df, aes(Var2, value, fill = Var1)) + geom_bar(stat = "identity") + scale_fill_manual(values = unname(celltype_treatment_colors)) + 
  mashaGgplot2Theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Number of cells")

## scRNA-seq:
  # Extended data figure 2 q
df <- reshape2::melt(table(scRNAseq_only$celltype_treatment, paste0(scRNAseq_only$MODIFIED_integrated_snn_res.0.2, "-clust")))
df$Var1 <- factor(as.character(df$Var1), levels = names(celltype_treatment_colors))
df$Var2 <- factor(as.character(df$Var2), levels = paste0(names(celltype_treatment_colors), "-clust"))
ggplot(df, aes(Var2, value, fill = Var1)) + geom_bar(stat = "identity") + scale_fill_manual(values = unname(celltype_treatment_colors)) + 
  mashaGgplot2Theme + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylab("Number of cells")

##---------------------------------------------##
##------------------Accuracy-------------------##
##---------------------------------------------##
## 
## Preparation
Liveseq_only$MODIFIED_RNA_snn_res.0.2 <- factor(Liveseq_only$MODIFIED_RNA_snn_res.0.2, levels = levels(Liveseq_only$celltype_treatment))
scRNAseq_only$MODIFIED_integrated_snn_res.0.2 <- factor(scRNAseq_only$MODIFIED_integrated_snn_res.0.2, levels = levels(scRNAseq_only$celltype_treatment))

Liveseq_only$Celltype_simple <- as.factor(Liveseq_only$Cell_type); levels(Liveseq_only$Celltype_simple) <- c("ASPC", "IBA", "RAW", "RAW")
Liveseq_only$MODIFIED_RNA_snn_res.0.2_CellType <- Liveseq_only$MODIFIED_RNA_snn_res.0.2; levels(Liveseq_only$MODIFIED_RNA_snn_res.0.2_CellType) <- c(rep("ASPC", 2), "IBA", rep("RAW", 2))

scRNAseq_only$Celltype_simple <- as.factor(scRNAseq_only$Cell_type); levels(scRNAseq_only$Celltype_simple) <- c("ASPC", "IBA", "RAW")
scRNAseq_only$MODIFIED_integrated_snn_res.0.2_CellType <- scRNAseq_only$MODIFIED_integrated_snn_res.0.2; levels(scRNAseq_only$MODIFIED_integrated_snn_res.0.2_CellType) <- c(rep("ASPC", 2), "IBA", rep("RAW", 2))

# ## 1. Accuracy for Cell type classification
#   #Live-seq
# t_live_CT <- table(Liveseq_only$Celltype_simple, Liveseq_only$MODIFIED_RNA_snn_res.0.2_CellType)
# Correct_Live_CT <- sum(diag(t_live_CT))/sum(t_live_CT) #0.9897959
# Wrong_Live_CT <- sum(t_live_CT[row(t_live_CT) != col(t_live_CT)])/sum(t_live_CT) #0.01020408 
#   #scRNA-seq
# t_sc_CT <- table(scRNAseq_only$Celltype_simple, scRNAseq_only$MODIFIED_integrated_snn_res.0.2_CellType)
# Correct_sc_CT <- sum(diag(t_sc_CT))/sum(t_sc_CT) #0.9981949
# Wrong_sc_CT <- sum(t_sc_CT[row(t_sc_CT) != col(t_sc_CT)])/sum(t_sc_CT) #0.001805054

## 2. Accuracy for cell type treatment classification 
## (§ Live-seq enables the stratification of cell types and states, ¶ 4)
  #Live-seq
t_live_noIBA <- table(Liveseq_only$celltype_treatment[Liveseq_only$celltype_treatment != "IBA_not_treated"], Liveseq_only$MODIFIED_RNA_snn_res.0.2[Liveseq_only$celltype_treatment != "IBA_not_treated"])
Correct_Live_noIBA <- sum(diag(t_live_noIBA))/sum(t_live_noIBA) #0.924
Wrong_Live_noIBA <- sum(t_live_noIBA[row(t_live_noIBA) != col(t_live_noIBA)])/sum(t_live_noIBA) #0.076 
  #scRNA-seq
t_sc_noIBA <- table(scRNAseq_only$celltype_treatment[scRNAseq_only$celltype_treatment != "IBA_not_treated"], scRNAseq_only$MODIFIED_integrated_snn_res.0.2[scRNAseq_only$celltype_treatment != "IBA_not_treated"])
Correct_sc_noIBA <- sum(diag(t_sc_noIBA))/sum(t_sc_noIBA) #0.9800499
Wrong_sc_noIBA <- sum(t_sc_noIBA[row(t_sc_noIBA) != col(t_sc_noIBA)])/sum(t_sc_noIBA) #0.01995012 

##---------------------------------------------##
##------------------Clustree-------------------##
##---------------------------------------------##
scRNAseq_only@meta.data <- scRNAseq_only@meta.data[, -grep("snn_res", colnames(scRNAseq_only@meta.data))]
Liveseq_only@meta.data <- Liveseq_only@meta.data[, -grep("snn_res", colnames(Liveseq_only@meta.data))]

for(i in seq(0.2, 1.6, 0.2)){
  Liveseq_only <- Liveseq_only %>% FindClusters(res = i)
  scRNAseq_only <- scRNAseq_only %>% FindClusters(res = i)
}

library(clustree)
clustree::clustree(Liveseq_only, suffix = "RNA_snn_res.") #Extended data Figure 2 h
clustree::clustree(scRNAseq_only, suffix = "integrated_snn_res.") #Extended data Figure 2 q
#5.7x6.08

# Plot clustering
# p1<- lapply( colnames(Liveseq_only@meta.data)[grep("RNA_snn_res", colnames(Liveseq_only@meta.data))], 
#              function (i) DimPlot(Liveseq_only, group.by = i, reduction = "tsne" ))
# marrangeGrob(p1, nrow = 2, ncol = 1)
# 
# p1<- lapply( colnames(scRNAseq_only@meta.data)[grep("integrated_snn_res", colnames(scRNAseq_only@meta.data))], 
#              function (i) DimPlot(scRNAseq_only, group.by = i, reduction = "tsne" ))
# marrangeGrob(p1, nrow = 2, ncol = 1)
