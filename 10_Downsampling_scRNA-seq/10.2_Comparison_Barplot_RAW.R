################################################################
#                                                              #
#      Compare DE results of Down-sampled and full data        #
#                                                              #
################################################################


### Author: Pernille
### Date: 01-02.2022

library(Seurat)
library(edgeR)
library(ggplot2)
library(dplyr)
library(ggbreak)

root_dir <- rprojroot::find_root(rprojroot::is_rstudio_project)

##---------------------------------------------##
##------------------LOAD DATA------------------##
##---------------------------------------------##
liveSeq_edgeR_raw <-
  readRDS(paste0(root_dir, "/data/liveseq_edgeR_res_RAW-NTvsLPSTreated.rds"))
scRNA_edgeR_res_raw <-
  readRDS(paste0(root_dir, "/data/scRNAeq_edgeR_res_RAW-NTvsLPSTreated.rds"))
sc_DownSampled_edgeR_res_raw <-
  readRDS(paste0(
    root_dir,
    "/data/DOWNSAMPLEDscRNAeq_edgeR_res_RAW_NTvsLPSTreated.rds"
  ))

##---------------------------------------------##
##------------------COMPARE--------------------##
##---------------------------------------------##

## DEFINE DE GENES EACH CATEGORIES
de_ls <-
  rownames(filter(
    as.data.frame(liveSeq_edgeR_raw),
    padj < 0.05 & abs(logFC) > 1 & high_pct == T
  ))
de_sc_ds <-
  rownames(filter(
    as.data.frame(sc_DownSampled_edgeR_res_raw),
    padj < 0.05 & abs(logFC) > 1 & high_pct == T
  ))
de_both <-
  intersect(rownames(filter(
    as.data.frame(liveSeq_edgeR_raw),
    padj < 0.05 & abs(logFC) > 1 & high_pct == T
  )),
  rownames(filter(
    as.data.frame(sc_DownSampled_edgeR_res_raw),
    padj < 0.05 & abs(logFC) > 1
  )))

liveSeq_edgeR_raw$table$Detected_as_DE_DOWNSAMPLEsc <- "None"
liveSeq_edgeR_raw$table[de_ls, "Detected_as_DE_DOWNSAMPLEsc"] <-
  "Only_Liveseq"
liveSeq_edgeR_raw$table[de_sc_ds, "Detected_as_DE_DOWNSAMPLEsc"] <-
  "Only_scDownsampled"
liveSeq_edgeR_raw$table[de_both, "Detected_as_DE_DOWNSAMPLEsc"] <-
  "Both"

de_ls <-
  rownames(filter(
    as.data.frame(liveSeq_edgeR_raw),
    padj < 0.05 & abs(logFC) > 1 & high_pct == T
  ))
de_both <-
  intersect(rownames(filter(
    as.data.frame(liveSeq_edgeR_raw),
    padj < 0.05 & abs(logFC) > 1 & high_pct == T
  )),
  rownames(filter(
    as.data.frame(scRNA_edgeR_res_raw), padj < 0.05 & abs(logFC) > 1
  )))
de_sc <-
  rownames(filter(as.data.frame(scRNA_edgeR_res_raw), padj < 0.05 &
                    abs(logFC) > 1))

liveSeq_edgeR_raw$table$Detected_as_DE_sc <- "None"
liveSeq_edgeR_raw$table[de_ls, "Detected_as_DE_sc"] <-
  "Only_Liveseq"
liveSeq_edgeR_raw$table[de_sc, "Detected_as_DE_sc"] <- "Only_sc"
liveSeq_edgeR_raw$table[de_both, "Detected_as_DE_sc"] <- "Both"

## PREPARE DATA.FRAME FOR PLOT
df <-
  reshape2::melt(
    table(
      liveSeq_edgeR_raw$table$Detected_as_DE_DOWNSAMPLEsc,
      liveSeq_edgeR_raw$table$Detected_as_DE_sc
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
# Extended data Figure 2 y (left)
p <-
  ggplot(df, aes(x = Var2, y = value, fill = Var1)) + geom_bar(stat = "identity") + mashaGgplot2Theme +
  scale_y_break(c(1300, 8000)) + theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  scale_fill_manual(values = c("gray80",  "#2C8BBD",  "#B9E3BB", "#7ACCC3"))

df %>% subset(Var2 == "Only_sc") %>% mutate(percent = value / sum(value))
# Var1    Var2 value   percent
# 13               Both Only_sc     0 0.0000000
# 14               None Only_sc   988 0.7835052
# 15       Only_Liveseq Only_sc     0 0.0000000
# 16 Only_scDownsampled Only_sc   273 0.2164948


##---------------------------------------------##
##-----------------IN NUMBERS------------------##
##---------------------------------------------##

# (Not shown)
### DOWNSAMPLE
## N DE genes per categories:
tab <-
  table(liveSeq_edgeR_raw$table$Detected_as_DE_DOWNSAMPLEsc)[c("Both", "Only_Liveseq", "Only_scDownsampled")]
# Both       Only_Liveseq Only_scDownsampled
# 183                402                380
## In percent:
tab / sum(tab) * 100
#     Both       Only_Liveseq Only_scDownsampled
# 18.96373           41.65803           39.37824

### FULL
## N DE genes per categories:
tab <-
  table(liveSeq_edgeR_raw$table$Detected_as_DE_sc)[c("Both", "Only_Liveseq", "Only_sc")]
# Both Only_Liveseq      Only_sc
# 270          315         1261
## In percent:
tab / sum(tab) * 100
#     Both Only_Liveseq      Only_sc
# 14.62622     17.06392     68.30986
