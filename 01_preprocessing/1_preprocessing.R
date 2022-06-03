################################################################
#                                                              #
#                         Preprocessing                        #
#                                                              #
################################################################

### Questions: https://github.com/DeplanckeLab/Live-seq/issues
### Date: 2022-03-06
### Datasets: scRNA-seq and Live-seq 
### Goal: Preprocessing of scRNA-seq and Live-seq data

library(rprojroot)
root_dir <- find_root(has_file("Live-seq.RProj"))

source(paste0(root_dir, "/utils/utils.R"))

library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(scater)
library(reshape2)
library(Matrix)

# download count matrix
file <- paste0(root_dir, "/data/GSE141064_count.final.csv.gz")
if (!file.exists(file)) {
  download.file(
    "https://0-www-ncbi-nlm-nih-gov.brum.beds.ac.uk/geo/download/?acc=GSE141064&format=file&file=GSE141064%5Fcount%2Efinal%2Ecsv%2Egz",
    file
  )
}

# save packages versions 
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

## input data
# read count_matrix
count.all <- read.csv(paste0(root_dir, "/data/GSE141064_count.final.csv.gz"), row.names = 1)
dim(count.all)
# read meta data
meta.all <- read.csv(paste0(root_dir, "/data/meta.final.csv"), row.names = 1)
rownames(meta.all) <- meta.all$sample_ID
# check the cell sample_ID
all(colnames(count.all) == rownames(meta.all))


## Create Seurat object
Seu.all <- CreateSeuratObject(counts =  count.all, meta.data = meta.all , min.cells  = 2, min.features = 1, 
                           project = "Seu.all")
# plot.text( c("number of cells =", ncol(Seu.all) ))  # 2148



## Add feature meta
# Load Ensembl 87 annotation, EGFP and mCherry gene was included
geneName <- read.table(file = paste0(root_dir, "/data/mouseGeneTable87_mCherry_EGFP.txt"), sep = "\t", header = T, row.names = 1 )
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
saveRDS(Seu.all, paste0(root_dir,"/01_preprocessing/Seu.all.rds"))

