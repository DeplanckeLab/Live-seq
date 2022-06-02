################################################################
#                                                              #
#                    Comparison DE results                     #
#                 between Live and scRNA-seq                   #
#                    DOWNSAMPLE & BOOTSTRAP                    #
#                                                              #
################################################################


### Author: Pernille
### Date: 14.12.2021
### Datasets: scRNA-seq and Live seq 
### Goal: Compare DE results obtained with scRNA and Live-seq 
###       For Live-seq manuscript
###       Downsample scRNAseq data -> perform DE

library(Seurat)
library(edgeR)
library(ggplot2)
library(dplyr)
library(plyr)

root_dir <- find_root(has_file("Live-seq.RProj"))
source(paste0(root_dir,"/utils/myFunctions_DEacrossCT.R"))
source(paste0(root_dir,"/8.DEwithinCT/myFunctions_DEwithinCT.R"))

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
f_downsample <- function(i, dataset_to_sample = subset_sc_raw, 
                         nCount_origSC_all_cells = nCount_subset_raw,
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

live_raw <- readRDS(paste0(root_dir, "/data/Liveseq_RAW_seu.rds"))
sc_raw <- readRDS(paste0(root_dir, "/data/scRNAseq_RAW_seu.rds"))

##---------------------------------------------##
##-------------------EXPLORE-------------------##
##---------------------------------------------##

o_orig <- density_plots(live_raw, sc_raw, sc_matrix = F, tittle = "Original data")

##---------------------------------------------##
##-------------DOWNSAMPLE scRNA-seq------------##
##---------------------------------------------##
## Down-sample single-cell RNA-seq data to have similar distribution of nFeat as Live-seq data
# 1. Order the cells based on nFeat in Live-seq and scRNA-seq data and pair cells based on this metric
# 2. Downsample each scRNA-seq data to have the same nFeat as the paired live-seq cells
# 3. Up sample the down-sample scRNA-seq to reach same lib size as live-seq (based on paired cell once again)

table(live_raw$celltype_treatment, live_raw$Batch)
#                       6  7 8_8 9_4
# Raw264.7_not_treated  8 20  61  12
# Raw264.7_LPS_treated  4 16   9  10
table(live_raw$celltype_treatment)
# Raw264.7_not_treated Raw264.7_LPS_treated 
# 102                   50 
table(sc_raw$celltype_treatment, sc_raw$Batch)
#                       8_8 9_1 scRNA
# Raw264.7_not_treated  16  46    95
# Raw264.7_LPS_treated  15  46    88

## PREPARE DATA
# Sample cells from sc to have matching number of cells with live-seq
sample_NT_cells <- sample(colnames(sc_raw)[sc_raw$treatment == "not_treated"], 102)
sample_LPS_cells <- sample(colnames(sc_raw)[sc_raw$treatment == "LPS_treated"], 50)

# table(sc_raw$celltype_treatment[c(sample_NT_cells, sample_LPS_cells)], sc_raw$Batch[c(sample_NT_cells, sample_LPS_cells)])
# # 8_8 9_1 scRNA
# # Raw264.7_not_treated  12  33    56
# # Raw264.7_LPS_treated   4  11    24

# Subset sc data to remove EGFP and mCherry
subset_sc_raw <- GetAssayData(sc_raw, slot = "counts")[, c(sample_NT_cells, sample_LPS_cells)]
subset_sc_raw <- subset_sc_raw[!rownames(subset_sc_raw) %in% c("EGFP", "mCherry"),]

# Compute nCount and nFeat sc without EGFP and mCherry
nCount_subset_raw <- colSums(subset_sc_raw)
nFeat_subset_raw <-  colSums(subset_sc_raw > 0)

# Order cells by nFeat
sc_cell_ord <- sort(nFeat_subset_raw)
l_cell_ord <- sort(colSums(GetAssayData(live_raw, slot = "counts")[ !rownames(live_raw) %in% c("EGFP", "mCherry"), ] > 0)) 

## DOWNSAMPLE
# Downsample and recreate count matrix
t <- lapply(1:length(sc_cell_ord), function(i) f_downsample(i, dataset_to_sample = subset_sc_raw, nCount_origSC_all_cells = nCount_subset_raw))
down_sc_raw <- do.call(cbind,t)
colnames(down_sc_raw) <- names(sc_cell_ord)

# Compute meta of downsample data
down_meta_raw <- data.frame(row.names = colnames(down_sc_raw), nCount = colSums(down_sc_raw), nFeat = colSums(down_sc_raw > 0))
down_meta_raw$prop <- down_meta_raw$nCount/nCount_subset_raw[rownames(down_meta_raw)]

p_down <- density_plots(live_raw, down_sc_raw, tittle = "Down-sampled sc")

## UPSAMPLE
# Calculate nCount of live-seq cells without considering EGFP and mCherry
nCount_l <- colSums(GetAssayData(live_raw, slot = "counts")[ !rownames(live_raw) %in% c("EGFP", "mCherry"), ]) 

# Up sample
up_sc_raw <- do.call(cbind, lapply(1:ncol(down_sc_raw), function(x) upOrdown_sample(down_sc_raw[,x, drop = F], nCount_l[names(l_cell_ord)[x]])))
colnames(up_sc_raw) <- colnames(down_sc_raw)

p_up <- density_plots(live_raw, up_sc_raw, tittle = "Up scaled sc")

## CREATE SEU DOWNSAMPLE
sc_raw_ds <- CreateSeuratObject(up_sc_raw)
sc_raw_ds$celltype_treatment <- sc_raw$celltype_treatment[colnames(sc_raw_ds)]
sc_raw_ds$Batch <- sc_raw_ds$orig.ident

##---------------------------------------------##
##-------------------DE EdgeR------------------##
##---------------------------------------------##
## Perform DE on down-sampled sc data
# We want to see how the results stands compare to full scRNA-seq data so 
# we perform the DE on the genes that are tested for DE in the scRNA-seq data
g <- select_genes(sc_raw, "celltype_treatment")

sc_DownSampled_edgeR_res_raw <- f_edgeR_CT(sc_raw_ds, genes_to_keep = g)

##---------------------------------------------##
##----------------CALCULATE PCT----------------##
##---------------------------------------------##

df_l_NT <- as.matrix(GetAssayData(sc_raw_ds, slot = "counts")[,sc_raw_ds$celltype_treatment == "Raw264.7_not_treated"])
df_l_LPS <- as.matrix(GetAssayData(sc_raw_ds, slot = "counts")[,sc_raw_ds$celltype_treatment == "Raw264.7_LPS_treated"])
prop_NT <- rowSums(df_l_NT > 2)/ncol(df_l_NT)
prop_LPS <- rowSums(df_l_LPS > 2)/ncol(df_l_LPS)

sc_DownSampled_edgeR_res_raw$table$pct_NT <- prop_NT[rownames(sc_DownSampled_edgeR_res_raw$table)]
sc_DownSampled_edgeR_res_raw$table$pct_LPS <- prop_LPS[rownames(sc_DownSampled_edgeR_res_raw$table)]
sc_DownSampled_edgeR_res_raw$table$high_pct <- F
sc_DownSampled_edgeR_res_raw$table$high_pct[sc_DownSampled_edgeR_res_raw$table$pct_NT > 0.15 | sc_DownSampled_edgeR_res_raw$table$pct_LPS > 0.15] <- T

saveRDS(sc_DownSampled_edgeR_res_raw, paste0(root_dir,"/data/DOWNSAMPLEDscRNAeq_edgeR_res_RAW_NTvsLPSTreated.rds"))

##---------------------------------------------##
##----------------PLOT logFC COMP---------------##
##---------------------------------------------##
liveSeq_edgeR_raw <- readRDS(paste0(root_dir,"data/liveseq_edgeR_res_RAW-NTvsLPSTreated.rds"))

## Extended data figure 2 x (left)
P_t <- f_plot_CT(liveSeq_edgeR_raw$table, sc_DownSampled_edgeR_res_raw$table) #9964

P_t$lm_AllGenes 
# Adjusted R-squared:  0.1272 
P_t$lm_DEsc 
# Adjusted R-squared:  0.6102 
P_t$lm_DEboth 
# Adjusted R-squared:  0.7487
P_t$pearson_AllGenes 
# [1] 0.3568188
P_t$pearson_DEsc 
# [1] 0.7816351
P_t$pearson_DEboth 
# [1] 0.8661536

