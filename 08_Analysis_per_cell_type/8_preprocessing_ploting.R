################################################################
#                                                              #
#                    FIGURES PER CELL TYPES                    #
#                                                              #
################################################################


### Author: Pernille Rainer, pernille.rainer@epfl.ch
### Date: 01-02.2022
### Datasets: scRNA-seq and Live seq per cell types
### Goal: Compute figures Live-seq manuscript per cell types:
# 1. tSNE per cell types colored by metadata and clustering
# 2. tSNE colored per extracted volumes 
# 3. Avg expression per cell type - correlation between sc and live

library(Seurat)
library(edgeR)
library(ggplot2)
library(dplyr)

## CT Cell type 
## SM Sampling Method

##---------------------------------------------##
##------------------LOAD DATA------------------##
##---------------------------------------------##

root_dir <- find_root(has_file("Live-seq.RProj"))
source(paste0(root_dir, "/utils/utils.R"))

Liveseq_all <- readRDS(paste0(root_dir, "/data/Liveseq_all_fluoresecene_added.3_7.5h.RDS"))
Liveseq_all$uniquely.mapped.rate <- Liveseq_all$uniquely.mapped / Liveseq_all$input.reads

scRNAseq_all <- readRDS("~/SVRAW1/wchen/data_analysis/Live_seq/final_analysis_V3/scRNAseq_only/scRNAseq.rds")
scRNAseq_all$uniquely.mapped.rate <- scRNAseq_all$uniquely.mapped / scRNAseq_all$input.reads

##---------------------------------------------##
##-----------PREPARE SEU OBJ PER CT------------##
##---------------------------------------------##
f_prepare <- function(all, CT , B, Sampling ){
  myseu <- subset(all, subset = Cell_type %in% CT & (Batch %in% B)  & (Species != "Human") & (sampling_type == Sampling))
  myseu <- subset(x = myseu, subset = nFeature_RNA > 1000 & percent.mt < 30 & uniquely.mapped.rate > 0.3)
  if(length(B) > 1){
    var_reg <- c("Batch", "nCount_RNA", "nFeature_RNA")}else{
      var_reg <- c("nCount_RNA", "nFeature_RNA")}
  myseu <- myseu %>% FindVariableFeatures() %>% NormalizeData() %>% ScaleData(vars.to.regress = var_reg) %>% RunPCA()
  myseu <- myseu %>% JackStraw(num.replicate = 100) %>% ScoreJackStraw(dims = 1:20)
  return(myseu)
}

# Live-seq - Raw cells
live_raw <- f_prepare(Liveseq_all, CT = c("Raw264.7_G9", "Raw264.7"), B = c("5Reseq", "6", "7", "8_8", "9_4"), Sampling = "Live_seq")
live_raw <- live_raw %>% RunTSNE(dims = 1:10, perplexity = 10) %>% FindNeighbors(dims =  1:10) %>% FindClusters(res = 0.3)
live_raw$celltype_treatment <- factor(live_raw$celltype_treatment, levels = c("Raw264.7_not_treated", "Raw264.7_LPS_treated"))

# Live-seq - ASPCs (use only 9_4 as only batch with both DMIR and non-treated cells)
live_aspcs <- f_prepare(Liveseq_all, CT = "ASPC", B = "9_4", Sampling = "Live_seq")
live_aspcs <- live_aspcs %>%  RunTSNE(dims = 1:10, perplexity = 10) %>% FindNeighbors(dims =  1:10) %>% FindClusters(res = 0.6)
live_aspcs$celltype_treatment <- factor(live_aspcs$celltype_treatment, levels = c("ASPC_not_treated", "ASPC_DMIR_treated"))

# scRNA-seq - RAW cells
sc_raw <- f_prepare(scRNAseq_all, CT = c("Raw264.7_G9", "Raw264.7"), B = c("5Reseq", "6",  "7", "scRNA","8_8", "9_1"), Sampling = "scRNA")
sc_raw <- sc_raw %>% RunTSNE(dims = 1:10, perplexity = 10) %>% FindNeighbors(dims =  1:10) %>% FindClusters(res = 0.1)
sc_raw$celltype_treatment <- factor(sc_raw$celltype_treatment, levels = c("Raw264.7_not_treated", "Raw264.7_LPS_treated"))

# scRNA-seq - ASPCs (use only 9_1 as only batch with both DMIR and non-treated cells)
sc_aspcs <- f_prepare(scRNAseq_all, CT = "ASPC", B = c("9_1"), Sampling = "scRNA")
sc_aspcs <- sc_aspcs %>% RunTSNE(dims = 1:10, perplexity = 10) %>% FindNeighbors(dims =  1:10) %>% FindClusters(res = 0.2)
sc_aspcs$celltype_treatment <- factor(sc_aspcs$celltype_treatment, levels = c("ASPC_not_treated", "ASPC_DMIR_treated"))

### Live-seq - IBA cells (for Volume extracted plot only)
live_iba <- f_prepare(Liveseq_all, CT = "IBA", B = c("5Reseq", "6","7", "8_8","9_4"), Sampling = "Live_seq")
live_iba <- live_iba %>% RunTSNE(dims = 1:10, perplexity = 10)


saveRDS(live_aspcs, 
        file = paste0(root_dir, "/data/Liveseq_ASPCs_seu.Rds"))
saveRDS(sc_aspcs, 
        file =  paste0(root_dir, "/data/scRNAseq_ASPCs_seu.Rds"))

saveRDS(live_raw, 
        file =  paste0(root_dir, "/data/Liveseq_RAW_seu.Rds"))
saveRDS(sc_raw, 
        file =  paste0(root_dir, "/data/scRNAseq_RAW_seu.Rds"))

saveRDS(live_iba, 
        file =  paste0(root_dir, "/data/Liveseq_IBA_seu.Rds"))

##---------------------------------------------##
##------------1. PLOT tSNEs METADATA-----------##
##---------------------------------------------##
f_tsne_plots <- function(myseu, clust_res = 0.3){
  p <- DimPlot(myseu, reduction = "tsne", group.by = "celltype_treatment", cols = unname(celltype_treatment_colors[levels(myseu$celltype_treatment)])) + theme(legend.position = "none")
  p1 <- DimPlot(myseu, reduction = "tsne", group.by = paste0("RNA_snn_res.", clust_res), cols = unname(celltype_treatment_colors[levels(myseu$celltype_treatment)])) + theme(legend.position = "none")
  p2 <- FeaturePlot(myseu, reduction = "tsne", feature = "nFeature_RNA") + theme(legend.position = "right")
  p3 <- FeaturePlot(myseu, reduction = "tsne", feature = "nCount_RNA") + theme(legend.position = "right")
  p4 <- DimPlot(myseu, reduction = "tsne", group.by = "Batch") + theme(legend.position = "none")
  
  return(grid.arrange(p, p2+ theme(legend.position = "none"), p1, p3+ theme(legend.position = "none"), p4, p2, p3, ncol = 2))
}

#Extended Data figure 2 i, j
P <- f_tsne_plots(live_aspcs, clust_res = 0.6)
P <- f_tsne_plots(live_raw, clust_res = 0.3)

#P <- f_tsne_plots(sc_raw, clust_res = 0.1)
#P <- f_tsne_plots(sc_aspcs, clust_res = 0.2)

##---------------------------------------------##
##--------2. PLOT tSNEs VOLUME EXTRACTED-------##
##---------------------------------------------##

## Extended data Figure 3c
extracted_vol <- read.csv(paste0(root_dir, "/data/InfoContent_Updated_VX-ASPC-9_4.csv"), header = T, row.names = 1)

live_aspcs$Volume.extract <- extracted_vol[colnames(live_aspcs), "Volume.extract"]
live_raw$Volume.extract <- extracted_vol[colnames(live_raw), "Volume.extract"]
live_iba$Volume.extract <- extracted_vol[colnames(live_iba), "Volume.extract"]

Liveseq_only <- readRDS(paste0(rootdir, "/data/Liveseq_only.rds"))
Liveseq_only$Volume.extract <- extracted_vol[colnames(Liveseq_only), "Volume.extract"]

p1 <- FeaturePlot(live_aspcs, reduction = "tsne", features = "Volume.extract")
p2 <- FeaturePlot(live_raw, reduction = "tsne", features = "Volume.extract")
p3 <- FeaturePlot(live_iba, reduction = "tsne", features = "Volume.extract")
p4 <- FeaturePlot(Liveseq_only, reduction = "tsne", features = "Volume.extract")

grid.arrange( p4 + ggtitle("ALL"), p1 + ggtitle("ASPCs"), p3 + ggtitle("IBA"), p2 + ggtitle("RAW"), ncol = 4)
#11.7x2.613 

##---------------------------------------------##
##--------3. Corr AVG EXP PER CT PER SM--------##
##---------------------------------------------##

## (Not shown)

#Correlation of average expression between sampling methods and calculated per cell type and treatments
cal_avg <- function(myseu){
  d <- as.data.frame(t(GetAssayData(myseu)))
  d$Celltype_treatment <- myseu$celltype_treatment
  d <- d %>% group_by(Celltype_treatment) %>% dplyr::summarise(across(.cols = everything(), mean, na.rm = TRUE))
  return(d)
}
avg_exp_plot <- function(myseu_sc, myseu_live){
  live <- reshape2::melt(cal_avg(myseu_live), id = "Celltype_treatment")
  colnames(live) <- c("Celltype_treatment", "variable", "live_avg_exp")
  live$variable <- as.character(live$variable)
  sc <- reshape2::melt(cal_avg(myseu_sc), id = "Celltype_treatment")
  colnames(sc) <- c("Celltype_treatment", "variable", "sc_avg_exp")
  sc$variable <- as.character(sc$variable)
  
  to_keep <- intersect(sc$variable, live$variable)
  to_keep <- to_keep[!to_keep %in% c("EGFP", "mCherry", gene.blacklist$ensembl_gene_id)]
  
  live <- live %>% subset(variable %in% to_keep)
  sc <- sc %>% subset(variable %in% to_keep)
  d <- left_join(live, sc, by = c("variable", "Celltype_treatment"))

  p <- ggplot(d, aes(x = sc_avg_exp, y = live_avg_exp, col = Celltype_treatment)) + 
    geom_point(size = 0.8) + scale_color_manual(values = unname(celltype_treatment_colors[levels(d$Celltype_treatment)])) + 
    mashaGgplot2Theme + geom_abline(slope = 1, intercept = 0, linetype = "dashed", col = "gray") + 
    xlab("Avg. log norm. exp. scRNA-seq") + ylab("Avg. log norm. exp. Live-seq ")

  cor_res <- c(cor(d[d$Celltype_treatment == levels(d$Celltype_treatment)[1], "sc_avg_exp"], d[d$Celltype_treatment == levels(d$Celltype_treatment)[1], "live_avg_exp"]),
               cor(d[d$Celltype_treatment == levels(d$Celltype_treatment)[2], "sc_avg_exp"], d[d$Celltype_treatment == levels(d$Celltype_treatment)[2], "live_avg_exp"]))
  names(cor_res) <- levels(d$Celltype_treatment)
  
  return(list(p = p, cor_res = cor_res))
}

p_aspcs <- avg_exp_plot(sc_aspcs, live_aspcs)
p_raw <- avg_exp_plot(sc_raw, live_raw)
grid.arrange(p_aspcs$p + theme(legend.position = "none"), p_raw$p + theme(legend.position = "none"), ncol = 2)
#4.59
p_aspcs$cor_res
# ASPC_not_treated ASPC_DMIR_treated 
# 0.9696420         0.9536276 
p_raw$cor_res
# Raw264.7_not_treated Raw264.7_LPS_treated 
# 0.9599684            0.9525767 


