################################################################
#                                                              #
#                      GO TERMS ANALYSIS                       #
#       OF GENES DETECTED ONLY BY ONE SAMPLING METHODS         #
#                                                              #
################################################################

### Author: Pernille Rainer pernille.rainer@epfl.ch
### Date: 01-02.2022
### Datasets: scRNA-seq and Live seq DE genes per cell types 
### Goal: Identify GO Terms of genes detected only by live-seq or scRNA-seq to find any potential bias

library(topGO)
library(Mapping)
library(ggplot2)
library(dplyr)

root_dir <- find_root(has_file("Live-seq.RProj"))

##---------------------------------------------##
##-------------------FUNCTION------------------##
##---------------------------------------------##
# GO enrichment
GOenrichment <- function(selectedGenes, allGenesList, ont = "BP", 
                         topNodes = 10, algo_Ranks = "eLimFisher",
                         Mapping = "org.Mm.eg.db"){
  require(topGO)
  require(Mapping)
  # algo_Ranks <- "eLimFisher" or "classicFisher"
  allGenesList_bg <- rep(1, length(allGenesList))
  names(allGenesList_bg) <- allGenesList
  allGenesList_bg[selectedGenes] <- 0.01
  
  tg.1 <- new("topGOdata", description = "GO analysis",
              ontology =  ont, allGenes = allGenesList_bg,
              geneSel = function (x) {return (x < 0.05)},
              annot = annFUN.org ,
              nodeSize = topNodes , # minimum number of genes in a GO categorie
              ID = "ENSEMBL", mapping = Mapping)
  #GO.resKS <- runTest(tg.1, algorithm = "classic", statistic = "ks")
  GO.resF <- runTest(tg.1, algorithm = "classic", statistic = "fisher")
  #GO.resKS.elim <- runTest(tg.1, algorithm = "elim", statistic = "ks")
  GO.resF.elim <- runTest(tg.1, algorithm = "elim", statistic = "fisher")
  
  result <- GenTable(tg.1, classicFisher = GO.resF, eLimFisher = GO.resF.elim, orderBy = algo_Ranks,
                     ranksOf = algo_Ranks, topNodes = topNodes)
  return(result)
}
# Dot plot
dotPlot_GO <- function(d, ref){
  y <- d
  y$eLimFisher <- as.numeric(y$eLimFisher); y$classicFisher <- as.numeric(y$classicFisher)
  y <- y[order(y[, colnames(y) == ref], decreasing = F),]
  y$Term <- factor(y$Term, levels = y$Term)
  c <- colnames(y); c[grep(ref, colnames(y))] <- "Pvalue"
  colnames(y) <- c
  
  p <- ggplot(y, aes(x = Pvalue, y = Term)) + 
    geom_point(aes(size = Pvalue, color = "red")) +
    theme_bw(base_size = 14) + scale_size(trans = 'reverse')  + theme(legend.position = "none")
  
  return(p)
}

##---------------------------------------------##
##------------------LOAD DATA------------------##
##---------------------------------------------##

liveSeq_edgeR_aspc <- readRDS(paste0(root_dir, "/data/liveseq_edgeR_res_ASPCs_NTvsDMIRTreated.Rds"))
liveSeq_edgeR_raw <- readRDS(paste0(root_dir, "/data/liveseq_edgeR_res_RAW-NTvsLPSTreated.Rds"))

##---------------------------------------------##
##----------------COMPUTE GO RES---------------##
##---------------------------------------------##

GO_BP_aspcs_onlyL <- GOenrichment(selectedGenes = rownames(liveSeq_edgeR_aspc$table)[liveSeq_edgeR_aspc$table$Detected_as_DE_sc == "Only_Liveseq"], allGenesList = rownames(liveSeq_edgeR_aspc$table), ont = "BP" )
GO_CC_aspcs_onlyL <- GOenrichment(selectedGenes = rownames(liveSeq_edgeR_aspc$table)[liveSeq_edgeR_aspc$table$Detected_as_DE_sc == "Only_Liveseq"], allGenesList = rownames(liveSeq_edgeR_aspc$table), ont = "CC" )

GO_BP_aspcs_onlysc <- GOenrichment(selectedGenes = rownames(liveSeq_edgeR_aspc$table)[liveSeq_edgeR_aspc$table$Detected_as_DE_sc == "Only_sc"], allGenesList = rownames(liveSeq_edgeR_aspc$table), ont = "BP" )
GO_CC_aspcs_onlysc <- GOenrichment(selectedGenes = rownames(liveSeq_edgeR_aspc$table)[liveSeq_edgeR_aspc$table$Detected_as_DE_sc == "Only_sc"], allGenesList = rownames(liveSeq_edgeR_aspc$table), ont = "CC" )

GO_BP_raw_onlyL <- GOenrichment(selectedGenes = rownames(liveSeq_edgeR_raw$table)[liveSeq_edgeR_raw$table$Detected_as_DE_sc == "Only_Liveseq"], allGenesList = rownames(liveSeq_edgeR_raw$table), ont = "BP" )
GO_CC_raw_onlyL <- GOenrichment(selectedGenes = rownames(liveSeq_edgeR_raw$table)[liveSeq_edgeR_raw$table$Detected_as_DE_sc == "Only_Liveseq"], allGenesList = rownames(liveSeq_edgeR_raw$table), ont = "CC" )

GO_BP_raw_onlysc <- GOenrichment(selectedGenes = rownames(liveSeq_edgeR_raw$table)[liveSeq_edgeR_raw$table$Detected_as_DE_sc == "Only_sc"], allGenesList = rownames(liveSeq_edgeR_raw$table), ont = "BP" )
GO_CC_raw_onlysc <- GOenrichment(selectedGenes = rownames(liveSeq_edgeR_raw$table)[liveSeq_edgeR_raw$table$Detected_as_DE_sc == "Only_sc"], allGenesList = rownames(liveSeq_edgeR_raw$table), ont = "CC" )

##---------------------------------------------##
##---------------------PLOTS-------------------##
##---------------------------------------------##
p_CC_a_l <- dotPlot_GO(GO_CC_aspcs_onlyL, ref = "eLimFisher")
p_CC_r_l <- dotPlot_GO(GO_CC_raw_onlyL, ref = "eLimFisher")

p_BP_a_l <- dotPlot_GO(GO_BP_aspcs_onlyL, ref = "eLimFisher")
p_BP_r_l <- dotPlot_GO(GO_BP_raw_onlyL, ref = "eLimFisher")

p_CC_a_sc <- dotPlot_GO(GO_CC_aspcs_onlysc, ref = "eLimFisher")
p_CC_r_sc <- dotPlot_GO(GO_CC_raw_onlysc, ref = "eLimFisher")

p_BP_a_sc <- dotPlot_GO(GO_BP_aspcs_onlysc, ref = "eLimFisher")
p_BP_r_sc <- dotPlot_GO(GO_BP_raw_onlysc, ref = "eLimFisher")

grid.arrange(p_CC_a_l, p_CC_r_l,
             p_BP_a_l, p_BP_r_l,
             p_CC_a_sc, p_CC_r_sc,
             p_BP_a_sc, p_BP_r_sc,
             ncol = 2)
## Extended Data Figure 2 z
