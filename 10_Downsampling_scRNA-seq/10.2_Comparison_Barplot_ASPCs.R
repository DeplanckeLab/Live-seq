################################################################
#                                                              #
#      Compare DE results of Down-sampled and full data        #
#                                                              #
################################################################


### Author: Pernille Rainer pernille.rainer@epfl.ch
### Date: 01-02.2022
## Goal: Compare number of common genes between Live-seq vs scRNA-seq OR Live-seq vs Downsampled scRNA-seq

library(Seurat)
library(edgeR)
library(ggplot2)
library(dplyr)
library(ggbreak)
root_dir <- rprojroot::find_root(rprojroot::is_rstudio_project)

##---------------------------------------------##
##------------------LOAD DATA------------------##
##---------------------------------------------##
liveSeq_edgeR_aspc <-
  readRDS(paste0(
    root_dir,
    "/data/liveseq_edgeR_res_ASPCs_NTvsDMIRTreated.rds"
  ))
scRNA_edgeR_res_aspc <-
  readRDS(paste0(
    root_dir,
    "/data/scRNAeq_edgeR_res_ASPCs_NTvsDMIRTreated.rds"
  ))
sc_DownSampled_edgeR_res_aspc <-
  readRDS(
    paste0(
      root_dir,
      "/data/DOWNSAMPLEDscRNAeq_edgeR_res_ASPCs_NTvsDMIRTreated.rds"
    )
  )

##---------------------------------------------##
##------------------COMPARE--------------------##
##---------------------------------------------##

## DEFINE DE GENES FROM EACH CATEGORIES
de_ls <-
  rownames(filter(
    as.data.frame(liveSeq_edgeR_aspc),
    padj < 0.05 & abs(logFC) > 1 & high_pct == T
  ))
de_both <-
  intersect(rownames(filter(
    as.data.frame(liveSeq_edgeR_aspc),
    padj < 0.05 & abs(logFC) > 1 & high_pct == T
  )),
  rownames(filter(
    as.data.frame(sc_DownSampled_edgeR_res_aspc),
    padj < 0.05 & abs(logFC) > 1 & high_pct == T
  )))
de_sc_ds <-
  rownames(filter(
    as.data.frame(sc_DownSampled_edgeR_res_aspc),
    padj < 0.05 & abs(logFC) > 1 & high_pct == T
  ))
liveSeq_edgeR_aspc$table$Detected_as_DE_DOWNSAMPLEsc <- "None"
liveSeq_edgeR_aspc$table[de_ls[!de_ls %in% de_both], "Detected_as_DE_DOWNSAMPLEsc"] <-
  "Only_Liveseq"
liveSeq_edgeR_aspc$table[de_sc_ds[!de_sc_ds %in% de_both], "Detected_as_DE_DOWNSAMPLEsc"] <-
  "Only_scDownsampled"
liveSeq_edgeR_aspc$table[de_both, "Detected_as_DE_DOWNSAMPLEsc"] <-
  "Both"

de_ls <-
  rownames(filter(
    as.data.frame(liveSeq_edgeR_aspc),
    padj < 0.05 & abs(logFC) > 1 & high_pct == T
  ))
de_both <-
  intersect(rownames(filter(
    as.data.frame(liveSeq_edgeR_aspc),
    padj < 0.05 & abs(logFC) > 1 & high_pct == T
  )),
  rownames(filter(
    as.data.frame(scRNA_edgeR_res_aspc), padj < 0.05 & abs(logFC) > 1
  )))
de_sc_ds <-
  rownames(filter(as.data.frame(scRNA_edgeR_res_aspc), padj < 0.05 &
                    abs(logFC) > 1))
liveSeq_edgeR_aspc$table$Detected_as_DE_sc <- "None"
liveSeq_edgeR_aspc$table[de_ls[!de_ls %in% de_both], "Detected_as_DE_sc"] <-
  "Only_Liveseq"
liveSeq_edgeR_aspc$table[de_sc_ds[!de_sc_ds %in% de_both], "Detected_as_DE_sc"] <-
  "Only_sc"
liveSeq_edgeR_aspc$table[de_both, "Detected_as_DE_sc"] <- "Both"

## CREATE DATA.FRAME FOR PLOT
df <-
  reshape2::melt(
    table(
      liveSeq_edgeR_aspc$table$Detected_as_DE_DOWNSAMPLEsc,
      liveSeq_edgeR_aspc$table$Detected_as_DE_sc
    )
  )
df$Var1 <-
  factor(
    as.character(df$Var1),
    levels = c("None", "Both", "Only_Liveseq", "Only_scDownsampled")
  )
df$Var2 <-
  factor(as.character(df$Var2),
         levels = c("None", "Both", "Only_Liveseq", "Only_sc"))

## PLOT
## Extended data Figures 2 y (right)
p <-
  ggplot(df, aes(x = Var2, y = value, fill = Var1)) + geom_bar(stat = "identity") + mashaGgplot2Theme +
  scale_y_break(c(1500, 8500)) + theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  scale_fill_manual(values = c("gray80",  "#2C8BBD",  "#B9E3BB", "#7ACCC3"))

df %>% subset(Var2 == "Only_sc") %>% mutate(percent = value / sum(value))
# Var1    Var2 value   percent
# 13               Both Only_sc     0 0.0000000
# 14               None Only_sc  1165 0.8158263
# 15       Only_Liveseq Only_sc     0 0.0000000
# 16 Only_scDownsampled Only_sc   263 0.1841737


##---------------------------------------------##
##-----------------IN NUMBERS-------...--------##
##---------------------------------------------##

### DOWNSAMPLE
## N DE genes per categories:
tab <-
  table(liveSeq_edgeR_aspc$table$Detected_as_DE_DOWNSAMPLEsc)[c("Both", "Only_Liveseq", "Only_scDownsampled")]
# Both       Only_Liveseq Only_scDownsampled
# 256                341                302
## In percent:
tab / sum(tab) * 100
#     Both       Only_Liveseq Only_scDownsampled
# 28.47608           37.93103           33.59288

### FULL
## N DE genes per categories:
tab <-
  table(liveSeq_edgeR_aspc$table$Detected_as_DE_sc)[c("Both", "Only_Liveseq", "Only_sc")]
# Both Only_Liveseq      Only_sc
# 410          187         1428
## In percent:
tab / sum(tab) * 100
#      Both Only_Liveseq      Only_sc
# 20.246914     9.234568    70.518519
