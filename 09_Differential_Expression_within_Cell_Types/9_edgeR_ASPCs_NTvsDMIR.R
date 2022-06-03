################################################################
#                                                              #
#                    Comparison DE results                     #
#                 between Live and scRNA-seq                   #
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

root_dir <- rprojroot::find_root(rprojroot::is_rstudio_project)
source(paste0(root_dir, "/utils/utils.R"))
source(paste0(root_dir, "/utils/myFunctions_DEacrossCT.R"))
source(paste0(root_dir, "/utils/myFunctions_DEwithinCT.R"))

##---------------------------------------------##
##------------------LOAD DATA------------------##
##---------------------------------------------##
live_aspcs <-
  readRDS(paste0(root_dir, "/data/Liveseq_ASPCs_seu.rds"))
sc_aspcs <-
  readRDS(paste0(root_dir, "/data/scRNAseq_ASPCs_seu.rds"))

##---------------------------------------------##
##---------------RUN DE ANALYSIS---------------##
##---------------------------------------------##

# We want to see how live-seq stands compare to scRNA-seq data so
# we perform the DE on the genes that are tested for DE in the scRNA-seq data
genes_to_test <- select_genes(sc_aspcs, "celltype_treatment")

liveSeq_edgeR_aspc <- f_edgeR_CT(live_aspcs, genes_to_test)
scRNA_edgeR_res_aspc <- f_edgeR_CT(sc_aspcs, genes_to_test)

##---------------------------------------------##
##----------------CALCULATE PCT----------------##
##---------------------------------------------##

# Genes are defined as DE only if expressed with 2 counts in at least 15% cells in one group
df_l_NT <-
  as.matrix(GetAssayData(live_aspcs, slot = "counts")[, live_aspcs$celltype_treatment == "ASPC_not_treated"])
df_l_DMI <-
  as.matrix(GetAssayData(live_aspcs, slot = "counts")[, live_aspcs$celltype_treatment == "ASPC_DMIR_treated"])
prop_NT <- rowSums(df_l_NT > 2) / ncol(df_l_NT)
prop_DMI <- rowSums(df_l_DMI > 2) / ncol(df_l_DMI)

liveSeq_edgeR_aspc$table$pct_NT <-
  prop_NT[rownames(liveSeq_edgeR_aspc$table)]
liveSeq_edgeR_aspc$table$pct_DMI <-
  prop_DMI[rownames(liveSeq_edgeR_aspc$table)]
liveSeq_edgeR_aspc$table$high_pct <- F
liveSeq_edgeR_aspc$table$high_pct[liveSeq_edgeR_aspc$table$pct_NT > 0.15 |
                                    liveSeq_edgeR_aspc$table$pct_DMI > 0.15] <- T

saveRDS(
  liveSeq_edgeR_aspc,
  paste0(
    root_dir,
    "/data/liveseq_edgeR_res_ASPCs_NTvsDMIRTreated.rds"
  )
)
saveRDS(
  scRNA_edgeR_res_aspc,
  paste0(
    root_dir,
    "/data/scRNAeq_edgeR_res_ASPCs_NTvsDMIRTreated.rds"
  )
)

##---------------------------------------------##
##-------------PLOT COMPARAISON FC-------------##
##---------------------------------------------##

## Extended data Figure 2 w (right)
P_t <-
  f_plot_CT(liveSeq_edgeR_aspc$table, scRNA_edgeR_res_aspc$table) #10845
# P_t$lm_AllGenes
# Adjusted R-squared:  0.3805
# P_t$lm_DEsc
# Adjusted R-squared:  0.6896
# P_t$lm_DEboth
# Adjusted R-squared:  0.9168
# P_t$pearson_AllGenes
# 0.6168892
# P_t$pearson_DEsc
# 0.8305017
# P_t$pearson_DEboth
# 0.9575881
