################################################################
#                                                              #
#                    DIFFERENTIAL EXPRESSION                   #
#                   for each cluster vs Rest                   #
#                                                              #
################################################################

### Author: Pernille Rainer pernille.rainer@epfl.ch  
### Date: 6.09.2021
### Datasets: scRNA-seq and Live seq 
### Goal: Compute DE for each cluster (corresponding to a cell type/state) versus the rest 
###       Done for both Live-seq and scRNA-seq data
###       For Live-seq manuscript

## DE are performed per sampling method of one group (i.e., cluster) versus All the rest. The logFC obtained with the two sampling methods are then compared 
## CLUSTERING SELECTED: 
## Live-seq: MODIFIED_RNA_snn_res.0.2.Rds
## scRNA-seq: MODIFIED_integrated_snn_res.0.2

library(Seurat)
library(edgeR)
library(ggplot2)
library(dplyr)
library(dynutils)
library(data.table)

##---------------------------------------------##
##------------------LOAD DATA------------------##
##---------------------------------------------##
root_dir <- find_root(has_file("Live-seq.RProj"))
source(paste0(root_dir, "/utils/utils.R"))
source(paste0(root_dir, "/utils/myFunctions_DEacrossCT.R"))

# Seurat object
scRNAseq_only <- readRDS(paste0(root_dir, "/data/scRNAseq_only.Rds"))
Liveseq_only <- readRDS(paste0(root_dir, "/data/Liveseq_only.Rds"))

##---------------------------------------------##
##---------------RUN DE ANALYSIS---------------##
##---------------------------------------------##
## Select genes for DE: Express in at least 15% of one group with at least 2 counts 
  # Select genes from the scRNAseq data
genes_to_test <- select_genes(scRNAseq_only, "MODIFIED_integrated_snn_res.0.2")

## Compute DE
liveseq_edgeR_res_t <- f_edgeR(Liveseq_only, "MODIFIED_RNA_snn_res.0.2", genes_to_test )
scRNA_edgeR_res_t <- f_edgeR(scRNAseq_only, "MODIFIED_integrated_snn_res.0.2", genes_to_test)

## Adjust p-value
liveseq_edgeR_res_t <- lapply(liveseq_edgeR_res_t, function(x) f_adjPval(x))
scRNA_edgeR_res_t <- lapply(scRNA_edgeR_res_t, function(x) f_adjPval(x))

##---------------------------------------------##
##----------------CALCULATE PCT----------------##
##---------------------------------------------##

## Calculate the percentage of cells expressing each genes 

# Genes are defined as DE only if expressed with 2 counts in at least 15% cells in one group
liveseq_edgeR_res <- lapply(names(liveseq_edgeR_res_t), function(x) calculate_pcts(x, Liveseq_only, liveseq_edgeR_res_t[[x]]))
names(liveseq_edgeR_res) <- names(liveseq_edgeR_res_t)

scRNA_edgeR_res <- lapply(names(scRNA_edgeR_res_t), function(x) calculate_pcts(x, scRNAseq_only, scRNA_edgeR_res_t[[x]]))
names(scRNA_edgeR_res) <- names(scRNA_edgeR_res_t)

rm(scRNA_edgeR_res_t, liveseq_edgeR_res_t)

##---------------------------------------------##
##--------------------IS DE--------------------##
##---------------------------------------------##

# Define which genes are DE: padj < 0.05 & abs(logFC) > 1 & high_pct == T
liveseq_edgeR_res <- lapply(liveseq_edgeR_res, f_is_de)
scRNA_edgeR_res <- lapply(scRNA_edgeR_res, f_is_de)

liveseq_edgeR_res_filt <- lapply(liveseq_edgeR_res, function(x) return(x$table %>% filter(is_de == T)))
scRNA_edgeR_res_filt <- lapply(scRNA_edgeR_res, function(x) return(x$table %>% filter(is_de == T)))

names(liveseq_edgeR_res_filt) <- c("ASPC_NT", "ASPC_DMIR", "IBA", "RAW_Mock", "RAW_LPS")
liveseq_edgeR_res_filt <- lapply(names(liveseq_edgeR_res_filt), function(x) f_add(x, liveseq_edgeR_res_filt))
liveseq_edgeR_res_filt <- rbindlist(liveseq_edgeR_res_filt)
liveseq_edgeR_res_filt <- as.data.frame(liveseq_edgeR_res_filt)

names(scRNA_edgeR_res_filt) <- c("ASPC_NT", "ASPC_DMIR", "IBA", "RAW_Mock", "RAW_LPS")
scRNA_edgeR_res_filt <- lapply(names(scRNA_edgeR_res_filt), function(x) f_add(x, scRNA_edgeR_res_filt))
scRNA_edgeR_res_filt <- rbindlist(scRNA_edgeR_res_filt)
scRNA_edgeR_res_filt <- as.data.frame(scRNA_edgeR_res_filt)

saveRDS(liveseq_edgeR_res, file = paste0(root_dir, "/data/Liveseq_edgeR_Clust-vs-Rest_MODIFIEDres0.2.Rds"))
saveRDS(scRNA_edgeR_res, file = paste0(root_dir, "/data/scRNAseq_edgeR_Clust-vs-Rest_MODIFIEDres0.2.Rds"))

##---------------------------------------------##
##-------------PLOT COMPARAISON FC-------------##
##---------------------------------------------##
## Figure 2 e 
P_t <- lapply(names(liveseq_edgeR_res), 
              function(x) f_plot(liveseq_edgeR_res[[x]]$table,
                                 scRNA_edgeR_res[[x]]$table))

names(P_t) <- names(liveseq_edgeR_res)

P_t$de_ASPC_NT$p + ggtitle("ASPCs not treated")
P_t$de_ASPC_DMIR$p + ggtitle("ASPCs + DMIR")

P_t$de_Raw_NT$p + ggtitle("Raw264.7 not treated")
P_t$de_Raw_LPS$p + ggtitle("Raw264.7 + LPS")

P_t$de_IBA$p + ggtitle("IBA")

grid.arrange(P_t$de_Raw_NT$p + theme(legend.position = "none"), 
             P_t$de_Raw_LPS$p + theme(legend.position = "none"), nrow = 1)
#logFC_scRNAseq-vs-liveseq_RAW-NT-vsRest_RAW-LPS-vsRest.pdf
#5.17x2.43


##---------------------------------------------##
##--------BARPLOT COMPARISON N DE GENES--------##
##---------------------------------------------##

liveseq_edgeR_res_filt <- as.data.frame(liveseq_edgeR_res_filt)

comp_de <- lapply(names(table(liveseq_edgeR_res_filt$cat)), f_comp)
comp_de <- rbindlist(comp_de)
comp_de$de_cat <- "None"; comp_de$de_cat[comp_de$inlive] <- "live"; comp_de$de_cat[comp_de$insc] <- "sc"
comp_de$de_cat[comp_de$inlive & comp_de$insc] <- "Both";
comp_de$de_cat <- factor(comp_de$de_cat, levels = c("live", "sc", "Both"))
comp_de$cat <- factor(comp_de$cat, levels = c("ASPC_NT", "ASPC_DMIR", "IBA", "RAW_Mock", "RAW_LPS")); levels(comp_de$cat) <- c("ASPC_Pre", "ASPC_Post", "IBA", "RAW_Mock", "RAW_LPS")
## Extended data Figure 2 v 
ggplot(comp_de, aes(x = cat, y = value, fill = de_cat)) + geom_bar(stat = "identity") + ylab("Number of DE genes") + xlab("") + 
  scale_fill_manual(values = c("#B9E3BB", "#7ACCC3", "#2C8BBD")) + mashaGgplot2Theme + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
#4.05x3.98


##---------------------------------------------##
##----HEATMAP top 20 DE per Cell type/states---##
##---------------------------------------------##

## Live-seq - Extended data Figure 2 k
g <- c(rownames(liveseq_edgeR_DEgenes$de_ASPC_NT[1:20,]),
           rownames(liveseq_edgeR_DEgenes$de_ASPC_DMIR[1:20,]),
           rownames(liveseq_edgeR_DEgenes$de_IBA[1:20,]),
           rownames(liveseq_edgeR_DEgenes$de_Raw_NT[1:20,]),
           rownames(liveseq_edgeR_DEgenes$de_Raw_LPS[1:20,]))
g <- data.annot[g,]
g <- g[!rownames(g) %in% c("ENSMUSG00000053279", "ENSMUSG00000030116"),]
g <- g[!duplicated(g$Ensembl),]
d <- GetAssayData(Liveseq_only, slot = "scale.data")[g$Ensembl, ]
rownames(d) <- data.annot[rownames(d), "Name"]
d <- t(dynutils::scale_quantile(t(d)))
d <- d[,order(Liveseq_only$MODIFIED_RNA_snn_res.0.2)]
my_annot_c <- data.frame(row.names = colnames(d),
                         cellType = Liveseq_only@meta.data[colnames(d), "celltype_treatment"],
                         clust = Liveseq_only@meta.data[colnames(d), "MODIFIED_RNA_snn_res.0.2"])

## scRNA-seq  - Extended data Figure 2r
g <- c(rownames(scRNA_edgeR_DEgenes$de_ASPC_NT[1:20,]),
       rownames(scRNA_edgeR_DEgenes$de_ASPC_DMIR[1:20,]),
       rownames(scRNA_edgeR_DEgenes$de_IBA[1:20,]),
       rownames(scRNA_edgeR_DEgenes$de_Raw_NT[1:20,]),
       rownames(scRNA_edgeR_DEgenes$de_Raw_LPS[1:20,]))
g <- data.annot[g,]
g <- g[!rownames(g) %in% c("ENSMUSG00000053279", "ENSMUSG00000030116"),]
g <- g[!duplicated(g$Ensembl),]

d <- GetAssayData(scRNAseq_only, slot = "data")[g$Ensembl, ]
rownames(d) <- data.annot[rownames(d), "Name"]
d <- t(dynutils::scale_quantile(t(d)))
d <- d[,order(scRNAseq_only$MODIFIED_integrated_snn_res.0.2)]
my_annot_c <- data.frame(row.names = colnames(d),
                         cellType = scRNAseq_only@meta.data[colnames(d), "celltype_treatment"],
                         clust = scRNAseq_only@meta.data[colnames(d), "MODIFIED_integrated_snn_res.0.2"])


## PLOT
my_colour = list(
  cellType = celltype_treatment_colors,
  clust = celltype_treatment_colors
)

library(viridis)
pheatmap::pheatmap(d, #scale = "row", 
                   #clustering_method = "ward.D2",
                   cluster_rows = F, cluster_cols = F,
                   annotation_colors = my_colour,
                   #color = colorRampPalette(brewer.pal(n = 7, name = "RdYlBu"))(100),
                   color = viridis(100),
                   annotation = my_annot_c,
                   #annotation_row = my_annot_r,
                   border_color = "NA",
                   cutree_col = 2,
                   fontsize_row = 7)


