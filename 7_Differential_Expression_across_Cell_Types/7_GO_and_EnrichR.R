################################################################
#                                                              #
#                      GO TERMS ANALYSIS                       #
#       OF DE GENES CELL TYPE/STATES & SAMPLING METHOD         #
#                                                              #
################################################################

### Author: Wanze Chen wz.chen@siat.ac.cn, adapted by Pernille Rainer pernille.rainer@epfl.ch
### Date: 01-02.2022
### Datasets: 
### Goal: enrichR on DE genes for BP GO terms and Mouse Cell Atlas, of each cluster versus Rest per sampling method

library(enrichR)
library(dplyr)
library(reshape2)
library(pheatmap)
library(gplots)
library(viridis)

root_dir <- find_root(has_file("Live-seq.RProj"))
source(paste0(root_dir, "/utils/utils.R"))

##---------------------------------------------##
##-------------------FUNCTION------------------##
##---------------------------------------------##

heatmap_and_save <- function(enr_res, dbp, n_top = 10, Term_to_select = NULL){
  enrich.DBP <- lapply(enr_res, function(x) return(x[[dbp]]))
  if(!is.null(Term_to_select)){
    enrich.DBP.top <- Term_to_select
  }else{
    enrich.DBP.top <- unique(unlist(lapply(enrich.DBP, function(x) return(x$Term[1:n_top]))))
  }
  
  enrich.DBP <- rbindlist(enrich.DBP)
  
  enrich.DBP.dcast <-  dcast(enrich.DBP, Term ~ cat, value.var = "Adjusted.P.value")
  enrich.DBP.dcast[ is.na(enrich.DBP.dcast) ]  <- 1  ## the NA value is filled with 1
  
  enrich.DBP.dcast <- as.data.frame(enrich.DBP.dcast)
  rownames(enrich.DBP.dcast) <- enrich.DBP.dcast$Term
  enrich.DBP.dcast <- enrich.DBP.dcast[,-1]
  
  enrich_top <- enrich.DBP.dcast[enrich.DBP.top,]
  
  enrich_top <- -log10(enrich_top)
  
  quantile_breaks <- function(xs, n = 10) {
    breaks <- quantile(as.matrix(xs), probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
  }
  
  mat_breaks <- quantile_breaks(enrich_top, n = 400)
  p <- pheatmap(
    mat               = enrich_top,
    color             = viridis((length(mat_breaks))),
    breaks            = mat_breaks,
    border_color      = NA,
    show_colnames     = TRUE,
    show_rownames     = TRUE,
    # cluster_cols      = mat_cluster_cols,
    # annotation_row = rownames(enrich_top5),
    cluster_rows      = F,
    cluster_cols = F,
    # annotation_col    = mat_col,
    # annotation_colors = mat_colors,
    drop_levels       = TRUE,
    fontsize          = 7
    # main              = "Quantile Color Scale"
    
  )
  
  return(list(p = p, res_top = enrich_top))
}

##---------------------------------------------##
##------------------LOAD DATA------------------##
##---------------------------------------------##

liveseq_edgeR_res <- readRDS( file = paste0(root_dir, "/data/Liveseq_edgeR_Clust-vs-Rest.Rds"))
scRNA_edgeR_res <- readRDS(file = paste0(root_dir, "/data/scRNAseq_edgeR_Clust-vs-Rest.Rds"))

##---------------------------------------------##
##----------------COMPUTE GO RES---------------##
##---------------------------------------------##
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015", "Mouse_Gene_Atlas")

### LIVE-SEQ
# Select top 100 de genes
top100_degenes <- lapply(liveseq_edgeR_res, function(x) {x$table[x$table$is_de == T,"Name"][1:100]})
# enrichR
enrichR_res_live <- lapply(top100_degenes, function(x) enrichr(x, dbs))
# Add category (which group of cells)
enrichR_res_live <- lapply(names(enrichR_res_live), function(x) {out <- enrichR_res_live[[x]];
                                                                 out <- lapply(out, function(xx) {xx$cat = x; return(xx)}); 
                                                                 return(out)})

### scRNA-SEQ
top100_degenes <- lapply(scRNA_edgeR_res, function(x) {x$table[x$table$is_de == T,"Name"][1:100]})

enrichR_res_sc <- lapply(top100_degenes, function(x) enrichr(x, dbs))
enrichR_res_sc <- lapply(names(enrichR_res_sc), function(x) {out <- enrichR_res_sc[[x]];
                                                             out <- lapply(out, function(xx) {xx$cat = x; return(xx)}); 
                                                             return(out)})

##---------------------------------------------##
##-----------------PLOT & SAVE-----------------##
##---------------------------------------------##

## COMPUTE the heatmap

## Extended data figure 2 m
live_MGA <- heatmap_and_save(enrichR_res_live, dbp = "Mouse_Gene_Atlas", n_top = 6)

## Extended data figure 2 t
sc_MGA <- heatmap_and_save(enrichR_res_sc, dbp = "Mouse_Gene_Atlas", n_top = 6)

  # Select 5 among top 10 GO terms of each group (cluster)
## Extended data figure 2 l
live_BP_f <- heatmap_and_save(enrichR_res_live, dbp = "GO_Biological_Process_2015", 
                              Term_to_select = c("regulation of endothelial cell proliferation (GO:0001936)", "regulation of response to wounding (GO:1903034)",
                                                 "positive regulation of ion transport (GO:0043270)", "negative regulation of fibroblast growth factor receptor signaling pathway (GO:0040037)",
                                                 "regulation of epithelial cell proliferation (GO:0050678)",
                                                 
                                                 "extracellular matrix organization (GO:0030198)", "extracellular structure organization (GO:0043062)",
                                                 "regulation of vasculature development (GO:1901342)", "regulation of developmental growth (GO:0048638)",
                                                 "regulation of epithelial cell proliferation (GO:0050678)",
                                                 
                                                 "stem cell proliferation (GO:0072089)", "Wnt signaling pathway (GO:0016055)", "angiogenesis (GO:0001525)",
                                                 "hematopoietic stem cell proliferation (GO:0071425)", "angiogenesis involved in wound healing (GO:0060055)",
                                                 
                                                 "regulation of B cell receptor signaling pathway (GO:0050855)", "leukocyte activation (GO:0045321)",
                                                 "antigen receptor-mediated signaling pathway (GO:0050851)", "immune response-activating cell surface receptor signaling pathway (GO:0002429)",
                                                 "lymphocyte activation (GO:0046649)",
                                                 
                                                 "cellular response to lipopolysaccharide (GO:0071222)", "positive regulation of cytokine production (GO:0001819)",
                                                 "response to molecule of bacterial origin (GO:0002237)", "inflammatory response (GO:0006954)",
                                                 "chemotaxis (GO:0006935)"))

## Extended data figure 2 s
sc_BP_f <- heatmap_and_save(enrichR_res_sc, dbp = "GO_Biological_Process_2015", 
                            Term_to_select = c("regulation of endothelial cell proliferation (GO:0001936)", "regulation of response to wounding (GO:1903034)",
                                               "regulation of ion homeostasis (GO:2000021)", "positive regulation of endothelial cell proliferation (GO:0001938)",
                                               "regulation of vasculature development (GO:1901342)",
                                               
                                               "regulation of homeostatic process (GO:0032844)", "extracellular matrix organization (GO:0030198)", 
                                               "extracellular structure organization (GO:0043062)", "cellular response to lipid (GO:0071396)", 
                                               "regulation of developmental growth (GO:0048638)",
                                               
                                               "stem cell proliferation (GO:0072089)", "hematopoietic stem cell proliferation (GO:0071425)",
                                               "angiogenesis (GO:0001525)","angiogenesis involved in wound healing (GO:0060055)",
                                               "positive regulation of transforming growth factor beta production (GO:0071636)",
                                               
                                               "regulation of B cell receptor signaling pathway (GO:0050855)", "antigen receptor-mediated signaling pathway (GO:0050851)", 
                                               "leukocyte activation (GO:0045321)", "lymphocyte activation (GO:0046649)", "cell activation involved in immune response (GO:0002263)",
                                               
                                               "myeloid leukocyte activation (GO:0002274)",
                                               "inflammatory response (GO:0006954)", "cellular response to lipopolysaccharide (GO:0071222)", 
                                               "response to molecule of bacterial origin (GO:0002237)", "positive regulation of cytokine production (GO:0001819)"
                                               ))



