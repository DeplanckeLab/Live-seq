packages <-
  c(
    "Seurat",
    "ggplot2",
    "tidyverse",
    "cowplot",
    "ggpubr",
    "dendextend",
    "doSNOW",
    "viridis",
    "data.table",
    "dynutils"
  )
install.packages(packages)

BiocManager::install("edgeR")