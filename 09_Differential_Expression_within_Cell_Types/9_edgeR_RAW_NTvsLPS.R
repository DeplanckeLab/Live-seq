
################################################################
#                                                              #
#                    Comparison DE results                     #
#                    of Live and scRNA-seq                     #
#                  RAW cells - Non vs LPS-treated              #
#                                                              #
################################################################


### Author: Pernille
### Date: 6.09.2021
### Datasets: scRNA-seq and Live seq 
### Goal: Compare DE results obtained with scRNA and Live-seq 
###       For Live-seq manuscript

library(Seurat)
library(edgeR)
library(ggplot2)
library(dplyr)

root_dir <- find_root(has_file("Live-seq.RProj"))
source(paste0(root_dir, "/utils/utils.R"))
source(paste0(root_dir, "/utils/myFunctions_DEacrossCT.R"))
source(paste0(root_dir, "/utils/myFunctions_DEwithinCT.R"))

##---------------------------------------------##
##------------------LOAD DATA------------------##
##---------------------------------------------##
live_raw <- readRDS(paste0(root_dir, "/data/Liveseq_RAW_seu.rds"))
live_raw <- subset(live_raw, Batch != "5Reseq")
sc_raw <- readRDS(paste0(root_dir, "/data/scRNAseq_RAW_seu.rds"))

##---------------------------------------------##
##---------------RUN DE ANALYSIS---------------##
##---------------------------------------------##

# We want to see how live-seq stands compare to scRNA-seq data so 
# we perform the DE on the genes that are tested for DE in the scRNA-seq data
genes_to_test <- select_genes(sc_raw, "celltype_treatment")

liveSeq_edgeR_raw <- f_edgeR_CT(live_raw, genes_to_test)
scRNA_edgeR_res_raw <- f_edgeR_CT(sc_raw, genes_to_test)

##---------------------------------------------##
##----------------CALCULATE PCT----------------##
##---------------------------------------------##

# Genes are defined as DE only if expressed with 2 counts in at least 15% cells in one group
df_l_NT <- as.matrix(GetAssayData(live_raw, slot = "counts")[,live_raw$celltype_treatment == "Raw264.7_not_treated"])
df_l_LPS <- as.matrix(GetAssayData(live_raw, slot = "counts")[,live_raw$celltype_treatment == "Raw264.7_LPS_treated"])
prop_NT <- rowSums(df_l_NT > 2)/ncol(df_l_NT)
prop_LPS <- rowSums(df_l_LPS > 2)/ncol(df_l_LPS)

liveSeq_edgeR_raw$table$pct_NT <- prop_NT[rownames(liveSeq_edgeR_raw$table)]
liveSeq_edgeR_raw$table$pct_LPS <- prop_LPS[rownames(liveSeq_edgeR_raw$table)]
liveSeq_edgeR_raw$table$high_pct <- F
liveSeq_edgeR_raw$table$high_pct[liveSeq_edgeR_raw$table$pct_NT > 0.15 | liveSeq_edgeR_raw$table$pct_LPS > 0.15] <- T

saveRDS(liveSeq_edgeR_raw, "~/SVRAW1/prainer/FluidFM_exp5_Wanze/Participation_Manuscript_2021-22/Figures_Liveseq/DE-EdgeR_and_Downsampling/FULLData/RAW/liveseq_edgeR_res_RAW-Not-vs-LPSTreated.rds")
saveRDS(scRNA_edgeR_res_raw, "~/SVRAW1/prainer/FluidFM_exp5_Wanze/Participation_Manuscript_2021-22/Figures_Liveseq/DE-EdgeR_and_Downsampling/FULLData/RAW/scRNAeq_edgeR_res_RAW-Not-vs-LPSTreated.rds")

saveRDS(liveSeq_edgeR_raw, paste0(root_dir, "/data/liveseq_edgeR_res_RAW_NTvsLPSTreated.rds"))
saveRDS(scRNA_edgeR_res_raw, paste0(root_dir, "/data/scRNAeq_edgeR_res_RAW_NTvsLPSTreated.rds"))

##---------------------------------------------##
##-------------PLOT COMPARAISON FC-------------##
##---------------------------------------------##

## Extended Data Figure 2 w (left)
P_t <- f_plot_CT(DE_res_live = liveSeq_edgeR_raw$table, DE_res_sc = scRNA_edgeR_res_raw$table) #9964

#5.03x4.1
# P_t$lm_AllGenes
# Adjusted R-squared:  0.1998
# P_t$lm_DEsc
# Adjusted R-squared:  0.482 
# P_t$lm_DEboth
# Adjusted R-squared:  0.7413 
# P_t$pearson_AllGenes
# 0.4470811
# P_t$pearson_DEsc
# 0.694485
# P_t$pearson_DEboth
# 0.8648645
