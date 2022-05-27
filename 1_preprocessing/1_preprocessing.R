library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(scater)
library(reshape2)
library(Matrix)
# source('~/NAS2/wchen/data_analysis/Resource/R_Funcitons/My_Rfunctions.R', local=TRUE)
# Set working directory
setwd("~/SVFASRAW/wchen/data_analysis/Live_seq/final_analysis_V3/Code_github/")
# save packages versions 
writeLines(capture.output(sessionInfo()), "1_preprocessing/sessionInfo.txt")

## input data
# read count_matrix
count.all <- read.csv("count.final.csv", row.names = 1)
dim(count.all)
# read meta data
meta.all <- read.csv("meta.final.csv", row.names = 1)
rownames(meta.all) <- meta.all$sample_ID
# check the cell sample_ID
all(colnames(count.all) == rownames(meta.all))


## Create Seurat object
Seu.all <- CreateSeuratObject(counts =  count.all, meta.data = meta.all , min.cells  = 2, min.features = 1, 
                           project = "Seu.all")
plot.text( c("number of cells =", ncol(Seu.all) ))  # 2148



## Add feature meta
# Load Ensembl 87 annotation, EGFP and mCherry gene was included
geneName <- read.table(file = "mouseGeneTable87_mCherry_EGFP.txt", sep = "\t", header = T, row.names = 1 )
# reduce the duplicate records
geneName.uni <- geneName[!duplicated(geneName$ensembl_gene_id),]
# check all the feature of Liveseq_all object is included in the geneName.uni
all (rownames(Seu.all) %in% geneName.uni$ensembl_gene_id)
# order the geneName.uni as the feature of Liveseq object
geneName.uni <- geneName.uni[ match(rownames(Seu.all), geneName.uni$ensembl_gene_id)  , ]
all(geneName.uni$ensembl_gene_id == rownames(Seu.all))   ## must be true here
# add feature meta to the Liveseq object
Seu.all@assays$RNA@meta.features <- geneName.uni
rownames(Seu.all@assays$RNA@meta.features ) <- geneName.uni$ensembl_gene_id
# the feature meta rowname shall the same as the count row name
all(rownames(Seu.all@assays$RNA@meta.features ) == rownames(Seu.all))  ## must be true here



## QC 
VlnPlot(Seu.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rRNA", "percent.protein"), group.by = "sampling_type" , ncol = 3)
FeatureScatter(Seu.all, feature1 = "nFeature_RNA", feature2 = "nCount_RNA", group.by = "sampling_type")


FeatureScatter(Seu.all, feature1 = "percent.rRNA", feature2 = "nFeature_RNA", group.by = "sampling_type")
FeatureScatter(Seu.all, feature1 = "percent.mt", feature2 = "nFeature_RNA", group.by = "sampling_type")
FeatureScatter(Seu.all, feature1 = "percent.protein", feature2 = "nFeature_RNA", group.by = "sampling_type")

FeatureScatter(Seu.all, feature1 = "percent.rRNA", feature2 = "percent.protein", group.by = "sampling_type")
FeatureScatter(Seu.all, feature1 = "percent.mt", feature2 = "percent.protein", group.by = "sampling_type")
FeatureScatter(Seu.all, feature1 = "percent.mt", feature2 = "percent.rRNA", group.by = "sampling_type")


## save Liveseq_all seurat object
saveRDS(Seu.all, "1_preprocessing/Seu.all.RDS")





