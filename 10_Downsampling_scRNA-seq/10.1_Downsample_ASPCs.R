################################################################
#                                                              #
#                    Comparison DE results                     #
#                 between Live and scRNA-seq                   #
#                         DOWNSAMPLE                           #
#                                                              #
################################################################


### Author: Pernille Rainer pernille.rainer@epfl.ch
### Date: 01-02.2022
### Datasets: scRNA-seq and Live seq 
### Goal: Compare DE results obtained with scRNA and Live-seq 
###       For Live-seq manuscript
###       Downsample scRNAseq data so that match complexity of Live-seq -> perform DE

library(Seurat)
library(edgeR)
library(ggplot2)
library(dplyr)
library(plyr)
library(gridExtra)

root_dir <- find_root(has_file("Live-seq.RProj"))
source(paste0(root_dir,"/utils/myFunctions_DEacrossCT.R"))
source(paste0(root_dir,"/utils/myFunctions_DEwithinCT.R"))

##---------------------------------------------##
##------------------FUNCTION-------------------##
##---------------------------------------------##
# Compare density of nCount and nFeat sc and live-seq
density_plots <- function(live, sc, sc_matrix = T, tittle = ""){

  if(sc_matrix){
    sc <- data.frame(nCount_RNA = colSums(sc),
                     nFeature_RNA = colSums(sc > 0))
  }
  
  df_lib <- data.frame(lib = c(live$nCount_RNA, sc$nCount_RNA),
                       which = c(rep("live", ncol(live)), rep("sc", length(sc$nCount_RNA))))
  mu <- ddply(df_lib, "which", summarise, grp.mean=mean(lib))
  
  p_lib <- ggplot(df_lib, aes(x=lib, col = which))+
    geom_density(alpha=0.4) + 
    geom_vline(data=mu, aes(xintercept=grp.mean, color=which),
               linetype="dashed") + 
    mashaGgplot2Theme + 
    scale_color_manual(values = c("#2888C9","#CD6226"))
  
  df_feat <- data.frame(feat = c(live$nFeature_RNA, sc$nFeature_RNA),
                        which = c(rep("live", ncol(live)), rep("sc", length(sc$nCount_RNA))))
  mu <- ddply(df_feat, "which", summarise, grp.mean=mean(feat))
  
  p_feat <- ggplot(df_feat, aes(x=feat, col = which))+
    geom_density(alpha=0.4) + 
    geom_vline(data=mu, aes(xintercept=grp.mean, color=which),
               linetype="dashed") + 
    mashaGgplot2Theme + 
    scale_color_manual(values = c("#2888C9","#CD6226"))
  grid.arrange(p_lib, p_feat, ncol = 2, top = tittle)
}

# Sample lib.size values from vector x (count vector of a cell)
upOrdown_sample <- function(x, lib.size){
  
  # Prepare data.to.sample so that it repeats each ens.id the number of corresponding counts
  if(is.null(dim(x))){
    d <- data.frame(ens_id = names(x), freq = x)
  }else{d <- data.frame(ens_id = rownames(x), freq = x[,1]) }
  data.to.sample <- d %>% filter(freq != 0)
  data.to.sample <- unlist(sapply(1:nrow(data.to.sample), function(i) rep(data.to.sample[i, "ens_id"], data.to.sample[i, "freq"])))
  
  #Sampling
  data.sampled <- sample(data.to.sample, lib.size, replace = T)
  
  # Recreate count vector
  data.table.count <- table(data.sampled)
  data.output <- data.frame(row.names = d$ens_id, value = rep(0, nrow(d)))
  data.output[names(data.table.count), "value"] <- data.table.count
  
  return(data.output)
}

#  Downsample until reaching same n of feat as paired cell
#' f_downsample
#' @description Downsample until reaching same n of feat as paired cell
#' @param i which cell number
#' @param dataset_to_sample single-cell data to downsample
#' @param nCount_origSC_all_cells number of counts in sc data sets to set initial lib.size to sample
#' @param cell_ord_sc order of cell in sc to match paired cells
#' @param order of cells in live-seq to match paired cells
f_downsample <- function(i, dataset_to_sample = subset_sc_aspc, 
                         nCount_origSC_all_cells = nCount_subset_aspc,
                         cell_ord_sc = sc_cell_ord , cell_order_l = l_cell_ord){
  cat(paste("CELL", i, ": goal nFeat", cell_order_l[i])); cat("\n")
  
  # initialize first lib.size_test (the library size we will sample)
  # proportion of original lib.size depending on diff of nFeat between sc cells and paired live-seq cell
  lib.size_test <- round((1-(cell_ord_sc[i] - cell_order_l[i])/cell_ord_sc[i])*nCount_origSC_all_cells[names(cell_ord_sc[i])])
  
  counter  <- 0
  ok <- F
  
  while(counter < 50 & ok == F){ 
    # Sample from sc data
    out <- upOrdown_sample(dataset_to_sample[, names(cell_ord_sc[i])], lib.size_test)
    # Calculate nFeat in down-sampled data
    out_nFeat <- sum(out > 0)  
    cat("\r", "iteration:", counter," nFeat =", out_nFeat, sep = " ")
    
    # If absolute difference between goal nFeat and reached nFeat < 5 Stop otherwise reiterate
    if(abs(out_nFeat - cell_order_l[i]) < 5){ ok <- T}
    lib.size_test <- round((1-(out_nFeat - cell_order_l[i])/out_nFeat)*lib.size_test)
    counter <- counter + 1
  } # Max 50 iterations
  cat("\n")
  
  return(out)
}

##---------------------------------------------##
##------------------LOAD DATA------------------##
##---------------------------------------------##

live_aspcs <- readRDS(paste0(root_dir,"/data/Liveseq_ASPCs_seu.rds"))
sc_aspcs <- readRDS(paste0(root_dir, "/data/scRNAseq_ASPCs_seu.rds"))

##---------------------------------------------##
##-------------------EXPLORE-------------------##
##---------------------------------------------##

p_orig <- density_plots(live_aspcs, sc_aspcs, sc_matrix = F, tittle = "Original data")

##---------------------------------------------##
##-------------DOWNSAMPLE scRNA-seq------------##
##---------------------------------------------##

## Down-sample single-cell RNA-seq data to have similar distribution of nFeat as Live-seq data
# 1. Order the cells based on nFeat in Live-seq and scRNA-seq data and pair cells based on this metric
# 2. Downsample each scRNA-seq data to have the same nFeat as the paired live-seq cells
# 3. Up sample the down-sample scRNA-seq to reach same lib size as live-seq (based on paired cell once again)

## APPROX SAME NUMBER CELLs
table(live_aspcs$celltype_treatment)
# ASPC_DMIR_treated  ASPC_not_treated 
# 37                43 
table(sc_aspcs$celltype_treatment)
# ASPC_DMIR_treated  ASPC_not_treated 
# 35                44 

## PREPARE FOR DOWNSAMPLE
# Subset sc data to remove EGFP and mCherry
subset_sc_aspc <- GetAssayData(sc_aspcs, slot = "count")
subset_sc_aspc <- subset_sc_aspc[!rownames(subset_sc_aspc) %in% c("EGFP", "mCherry"),]

# Calculate nCount and nFeat of scRNA-seq data without EGFP and mCherry "genes"
nCount_subset_aspc <- colSums(subset_sc_aspc)
nFeat_subset_aspc <-  colSums(subset_sc_aspc > 0)

# Order cells by nFeat
sc_cell_ord <- sort(nFeat_subset_aspc)
l_cell_ord <- sort(colSums(GetAssayData(live_aspcs, slot = "counts")[ !rownames(live_aspcs) %in% c("EGFP", "mCherry"), ] > 0)) 

## DOWNSAMPLE
# Downsample and recreate count matrix
t <- lapply(1:length(sc_cell_ord), function(i) f_downsample(i,dataset_to_sample = subset_sc_aspc, nCount_origSC_all_cells = nCount_subset_aspc))
down_sc_aspc <- do.call(cbind,t)
colnames(down_sc_aspc) <- names(sc_cell_ord)

# Compute meta of downsample data
down_meta_aspcs <- data.frame(row.names = colnames(down_sc_aspc), nCount = colSums(down_sc_aspc), nFeat = colSums(down_sc_aspc > 0))
down_meta_aspcs$prop <- down_meta_aspcs$nCount/nCount_subset_aspc[rownames(down_meta_aspcs)]

p_down <- density_plots(live_aspcs, down_sc_aspc, tittle = "Down-sampled sc")

## UPSAMPLE
# Calculate nCount of live-seq cells without considering EGFP and mCherry
nCount_l <- colSums(GetAssayData(live_aspcs, slot = "counts")[ !rownames(live_aspcs) %in% c("EGFP", "mCherry"), ]) 

# Up sample
up_sc_aspc <- do.call(cbind, lapply(1:ncol(down_sc_aspc), function(x) upOrdown_sample(down_sc_aspc[,x, drop = F], nCount_l[names(l_cell_ord)[x]])))
colnames(up_sc_aspc) <- colnames(down_sc_aspc)

p_up <- density_plots(live_aspcs, up_sc_aspc, tittle = "Up-scaled sc")

## CREATE SEU OBJ OF DS scRNAseq 
sc_aspcs_ds <- CreateSeuratObject(up_sc_aspc)
sc_aspcs_ds$celltype_treatment <- sc_aspcs$celltype_treatment[colnames(sc_aspcs_ds)]
sc_aspcs_ds$Batch <- sc_aspcs_ds$orig.ident

grid.arrange(p_orig,  p_down, p_up, ncol = 1)
# 8.58x6.2 density_comparisons_datasets

##---------------------------------------------##
##-------------------DE EdgeR------------------##
##---------------------------------------------##

## Perform DE on down-sampled sc data
# We want to see how the results stands compare to full scRNA-seq data so 
# we perform the DE on the genes that are tested for DE in the scRNA-seq data
g <- select_genes(sc_aspcs, "celltype_treatment")

sc_DownSampled_edgeR_res_aspc <- f_edgeR_CT(sc_aspcs_ds, genes_to_keep = g)

##---------------------------------------------##
##----------------CALCULATE PCT----------------##
##---------------------------------------------##

# Genes are defined as DE only if expressed in at least 15% cells of a group with at least 2 counts
df_l_NT <- as.matrix(GetAssayData(sc_aspcs_ds, slot = "counts")[,sc_aspcs_ds$celltype_treatment == "ASPC_not_treated"])
df_l_DMI <- as.matrix(GetAssayData(sc_aspcs_ds, slot = "counts")[,sc_aspcs_ds$celltype_treatment == "ASPC_DMIR_treated"])
prop_NT <- rowSums(df_l_NT > 2)/ncol(df_l_NT)
prop_DMI <- rowSums(df_l_DMI > 2)/ncol(df_l_DMI)

sc_DownSampled_edgeR_res_aspc$table$pct_NT <- prop_NT[rownames(sc_DownSampled_edgeR_res_aspc$table)]
sc_DownSampled_edgeR_res_aspc$table$pct_DMI <- prop_DMI[rownames(sc_DownSampled_edgeR_res_aspc$table)]
sc_DownSampled_edgeR_res_aspc$table$high_pct <- F
sc_DownSampled_edgeR_res_aspc$table$high_pct[sc_DownSampled_edgeR_res_aspc$table$pct_NT > 0.15 | sc_DownSampled_edgeR_res_aspc$table$pct_DMI > 0.15] <- T

saveRDS(sc_DownSampled_edgeR_res_aspc, paste0(root_dir,"/data/DOWNSAMPLEDscRNAeq_edgeR_res_ASPCs_NTvsDMIRTreated.rds"))

##---------------------------------------------##
##----------------PLOT logFC COMP---------------##
##---------------------------------------------##
liveSeq_edgeR_aspc <- readRDS(paste0(root_dir,"/data/liveseq_edgeR_res_ASPCs_NTvsDMIRTreated.rds"))

## Extended data figure 2 x (right)
P_t <- f_plot_CT(liveSeq_edgeR_aspc$table, sc_DownSampled_edgeR_res_aspc$table) #10845

P_t$lm_AllGenes 
# Adjusted R-squared:  0.3268 
P_t$lm_DEsc 
# Adjusted R-squared:  0.7281 
P_t$lm_DEboth 
# Adjusted R-squared:  0.9337
P_t$pearson_AllGenes 
# [1] 0.5717463
P_t$pearson_DEsc 
# [1] 0.8535672
P_t$pearson_DEboth 
# [1] 0.9664163
