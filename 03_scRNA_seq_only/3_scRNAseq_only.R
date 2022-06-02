library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(scater)
library(reshape2)
library(Matrix)
library(ggplot2)
library(cowplot)
# Set working directory
setwd("~/SVFASRAW/wchen/data_analysis/Live_seq/final_analysis_V3/Code_github/")

###### version control
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
#######################

Seu.all <- readRDS("1_preprocessing/Seu.all.RDS")



### subset of scRNA
scRNA <- subset(Seu.all, subset = (sampling_type == "scRNA" & celltype_treatment %in% c("ASPC_not_treated", "ASPC_DMIR_treated", "IBA_not_treated", "Raw264.7_not_treated", "Raw264.7_LPS_treated") )) 
dim(scRNA)
### remove 20 genes in black list, which are derived from the 0 pg input RNA negative control. 
gene.blacklist <- read.csv("gene.blacklist.csv")
data.count <- as.matrix(scRNA@assays$RNA@counts) 
data.count <- data.count[ !rownames(data.count) %in% gene.blacklist$ensembl_gene_id, ]

##### remove ribosomal protein genes
gene.ribo <- scRNA@assays$RNA@meta.features$ensembl_gene_id[grep("^Rp[s,l]", scRNA@assays$RNA@meta.features$external_gene_name) ]
gene.ribo <- as.character(gene.ribo)
data.count <- data.count[ !rownames(data.count) %in% gene.ribo, ]


scRNA <- CreateSeuratObject(counts =data.count, meta.data = scRNA@meta.data )

## add the meta feature
scRNA@assays$RNA@meta.features <- Seu.all@assays$RNA@meta.features[ rownames(scRNA), ]

## add uniquely.mapped.rate and intron.mapped.rate in meta data
scRNA@meta.data$uniquely.mapped.rate <- scRNA@meta.data$uniquely.mapped / scRNA@meta.data$input.reads
scRNA@meta.data$intron.mapped.rate <- scRNA@meta.data$nCount_RNA / scRNA@meta.data$uniquely.mapped


## QC before filtering
dim(scRNA)
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rRNA", "percent.protein", "uniquely.mapped.rate", "intron.mapped.rate"), group.by = "celltype_treatment" , ncol = 4)

## filter and QC

scRNA <- subset(x = scRNA, subset = nFeature_RNA > 1000 & percent.mt < 30 & uniquely.mapped.rate > 0.3)
dim(scRNA)

## plot the cells passing filtering
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rRNA", "percent.protein", "uniquely.mapped.rate", "intron.mapped.rate"), group.by = "celltype_treatment" , ncol = 4)




p1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "celltype_treatment")
p2 <- FeatureScatter(scRNA, feature1 = "cDNA_concentration", feature2 = "nFeature_RNA", group.by = "celltype_treatment")
p3 <- FeatureScatter(scRNA, feature1 = "percent.rRNA", feature2 = "nFeature_RNA", group.by = "celltype_treatment")
p4 <- FeatureScatter(scRNA, feature1 = "percent.mt", feature2 = "nFeature_RNA", group.by = "celltype_treatment")
p5 <- FeatureScatter(scRNA, feature1 = "percent.protein", feature2 = "nFeature_RNA", group.by = "celltype_treatment")
p6 <- FeatureScatter(scRNA, feature1 = "percent.rRNA", feature2 = "percent.protein", group.by = "celltype_treatment")
p7 <- FeatureScatter(scRNA, feature1 = "percent.mt", feature2 = "percent.protein", group.by = "celltype_treatment")
p8 <- FeatureScatter(scRNA, feature1 = "percent.mt", feature2 = "percent.rRNA", group.by = "celltype_treatment")

plot_grid(p1, p2, p3,p4, p5, p6,p7,p8)



## data normalization 
scRNA <- NormalizeData(scRNA, verbose = FALSE)


### variable genes selection method 
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 500, 
                                    verbose = FALSE)

# Identify the most highly variable genes
top30 <- head(VariableFeatures(scRNA), 30)
plot1 <- VariableFeaturePlot(scRNA)
plot2 <- LabelPoints(plot = plot1, points = top30,  labels = scRNA@assays$RNA@meta.features[top30, "external_gene_name"] , repel = TRUE)
plot2

####    evaluate  cell cycle effects  ####
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- readRDS("s.genes.mouse.rds")    
g2m.genes <- readRDS("g2m.genes.mouse.rds")    

scRNA <- CellCycleScoring(scRNA, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


## scale all the genes
scRNA <- ScaleData(scRNA, verbose = FALSE)
# scRNA <- ScaleData(scRNA, verbose = FALSE,vars.to.regress = c("nFeature_RNA") )

scRNA <- RunPCA(scRNA, npcs = 20, verbose = FALSE)
scRNA <- JackStraw(scRNA, num.replicate = 100)
scRNA <- ScoreJackStraw(scRNA, dims = 1:20)
JackStrawPlot(scRNA, dims = 1:20)
ElbowPlot(scRNA)


scRNA <- RunTSNE(scRNA, reduction = "pca", dims = 1:10, verbose = FALSE, seed.use = 1000)
scRNA <- RunUMAP(scRNA, reduction = "pca", dims = 1:10)

scRNA <- FindNeighbors(scRNA, dims = 1:10 )  # , k.param =15
scRNA <- FindClusters(scRNA, resolution = 0.2)     ## resolution 0.3 
scRNA <- FindClusters(scRNA, resolution = 0.4)     ## resolution 0.4
scRNA <- FindClusters(scRNA, resolution = 0.5)     ## resolution 0.4
scRNA <- FindClusters(scRNA, resolution = 0.3)     ## resolution 0.3 
# levels(Idents(scRNA)) <- c(0,2,1)

p1 <- DimPlot(scRNA, reduction = "tsne")
p2 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype_treatment")
p3 <- DimPlot(scRNA, reduction = "tsne", group.by = "Batch")
p4<- FeaturePlot(scRNA, reduction = "tsne", features  = "nFeature_RNA" )
plot_grid(p1,p2,p3, p4)

saveRDS(scRNA  ,"3_scRNA_seq_only/scRNAseq.rds")






################  batchs correction using seurat intergration  ################
# 9_1 has different primers and enzyme batches, and indeed show most distant.
scRNA$Batch
scRNA$Batch_9_1 <- "others"
scRNA$Batch_9_1[ which(scRNA$Batch == "9_1")] <- "9_1"

scRNA.list <- SplitObject(scRNA, split.by = "Batch_9_1") 

scRNA.list[[1]]@meta.data <-  scRNA@meta.data[ colnames(scRNA.list[[1]]) , ]
scRNA.list[[2]]@meta.data <-  scRNA@meta.data[ colnames(scRNA.list[[2]]) , ]
# scRNA.list[[3]]@meta.data <-  scRNA@meta.data[ colnames(scRNA.list[[3]]) , ]


for (i in 1:length(scRNA.list )) {
  
  ## remove unwanted cell based on the different criteria
  scRNA.list [[i]] <- subset(  scRNA.list [[i]], subset = nFeature_RNA > 1000 & percent.mt < 30 & uniquely.mapped.rate > 0.3)
  scRNA.list [[i]] <- NormalizeData(scRNA.list [[i]], verbose = FALSE)
  scRNA.list [[i]] <- FindVariableFeatures(scRNA.list [[i]], selection.method = "vst", nfeatures = 500, 
                                                 verbose = FALSE)
  ## remove EGFP and mCherry from the variable gene, as it biases.
  # VariableFeatures(scRNA.list [[i]]) <- VariableFeatures(scRNA.list [[i]]) [!VariableFeatures(scRNA.list [[i]]) %in% c("EGFP", "mCherry") ]
  
}

## integrate Live-seq and scRNA-seq data
reference.list <- scRNA.list[c("others", "9_1")]
scRNA.anchors <- FindIntegrationAnchors(object.list = reference.list, reference = c(2),dims = 1:10, k.anchor = 5, k.filter = 50, k.score = 60 ) # k.anchor = 5, k.filter = 5 ,k.score = 30, 
## 
scRNA.integrated <- IntegrateData(anchorset = scRNA.anchors, features.to.integrate = rownames(scRNA.list[[1]]) , dims = 1:10)

all(rownames(scRNA.integrated@assays$RNA@meta.features) == rownames( scRNA@assays$RNA@meta.features))
scRNA.integrated@assays$RNA@meta.features <- scRNA@assays$RNA@meta.features

############################################


#######################################
DefaultAssay(scRNA.integrated) <- "RNA"


scRNA.integrated <- NormalizeData(scRNA.integrated, verbose = FALSE)


# scRNA.RNA <- NormalizeData(scRNA.RNA, verbose = FALSE, normalization.method = "CLR")
### variable genes selection method "mean.var.plot"  

scRNA.integrated <- FindVariableFeatures(scRNA.integrated, selection.method = "vst", nfeatures = 500, 
                                               verbose = FALSE)  ## 

## remove EGFP and mCherry from the variable gene, as it biases.
# VariableFeatures(scRNA) <- VariableFeatures(scRNA) [!VariableFeatures(scRNA) %in% c("EGFP", "mCherry") ]


# Identify the 10 most highly variable genes
top30 <- head(VariableFeatures(scRNA.integrated), 30)
plot1 <- VariableFeaturePlot(scRNA.integrated)
plot2 <- LabelPoints(plot = plot1, points = top30,  labels = scRNA.integrated@assays$RNA@meta.features[top30, "external_gene_name"] , repel = TRUE)
plot2

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- readRDS("s.genes.mouse.rds")    
g2m.genes <- readRDS("g2m.genes.mouse.rds")    

scRNA.integrated <- CellCycleScoring(scRNA.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)





# Run the standard workflow for visualization and clustering
scRNA.integrated <- ScaleData(scRNA.integrated, verbose = FALSE)
scRNA.integrated <- RunPCA(scRNA.integrated, npcs  = 20, verbose = FALSE)
scRNA.integrated <- JackStraw(scRNA.integrated, num.replicate = 100)
scRNA.integrated <- ScoreJackStraw(scRNA.integrated, dims = 1:20)
JackStrawPlot(scRNA.integrated, dims = 1:20)
ElbowPlot(scRNA.integrated)

scRNA.integrated <- RunUMAP(scRNA.integrated, reduction = "pca", dims = 1:10)  # , n.neighbors = 20
scRNA.integrated <- RunTSNE(scRNA.integrated, reduction = "pca", dims = 1:10)
scRNA.integrated <- FindNeighbors(scRNA.integrated, dims = 1:10)
scRNA.integrated <- FindClusters(scRNA.integrated, resolution = 0.2)  ## ", algorithm =4"  Leiden algorithm outperform a bit here, 3 cells misassign instead of 5
scRNA.integrated <- FindClusters(scRNA.integrated, resolution = 0.4)  ## ", algorithm =4"  Leiden algorithm outperform a bit here, 3 cells misassign instead of 5
scRNA.integrated <- FindClusters(scRNA.integrated, resolution = 0.5)  ## ", algorithm =4"  Leiden algorithm outperform a bit here, 3 cells misassign instead of 5
scRNA.integrated <- FindClusters(scRNA.integrated, resolution = 0.7)  ## ", algorithm =4"  Leiden algorithm outperform a bit here, 3 cells misassign instead of 5
scRNA.integrated <- FindClusters(scRNA.integrated, resolution = 0.6)  ## ", algorithm =4"  Leiden algorithm outperform a bit here, 3 cells misassign instead of 5
scRNA.integrated <- FindClusters(scRNA.integrated, resolution = 0.3)  ## ", algorithm =4"  Leiden algorithm outperform a bit here, 3 cells misassign instead of 5

### evaluate the clustering accuracy 
scRNA.integrated@meta.data$RNA.ident <- Idents(scRNA.integrated)
ggplot(scRNA.integrated@meta.data,aes(x=RNA.ident,fill=celltype_treatment))+
  geom_bar()+
  geom_text(aes(label=stat(count)),stat="count",position=position_stack(0.5)) +
  xlab("cluster")


p1 <- DimPlot(scRNA.integrated, reduction = "umap")
p2 <- DimPlot(scRNA.integrated, reduction = "umap", group.by = "celltype_treatment")
p3 <- DimPlot(scRNA.integrated, reduction = "umap", group.by = "Batch")
p4 <- FeaturePlot(scRNA.integrated, reduction  = "umap", features = "nFeature_RNA")
plot_grid(p1, p2, p3,p4, align = "hv", axis = "lrtb")

p1 <- DimPlot(scRNA.integrated, reduction = "pca")
p2 <- DimPlot(scRNA.integrated, reduction = "pca", group.by = "celltype_treatment")
p3 <- DimPlot(scRNA.integrated, reduction = "pca", group.by = "Batch")
p4 <- FeaturePlot(scRNA.integrated, reduction  = "pca", features = "nFeature_RNA")
plot_grid(p1, p2, p3,p4, align = "hv", axis = "lrtb")

p1 <- DimPlot(scRNA.integrated, reduction = "tsne")
p2 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "celltype_treatment")
p3 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "Batch")
p4 <- FeaturePlot(scRNA.integrated, reduction  = "tsne", features = "nFeature_RNA")
p5 <- DimPlot(scRNA.integrated, reduction  = "tsne", group.by = "Phase")
plot_grid(p1, p2, p3,p4,p5, align = "hv", axis = "lrtb")


p1 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "RNA_snn_res.0.3")
p2 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "RNA_snn_res.0.4")
p3 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "RNA_snn_res.0.5")
p4 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "RNA_snn_res.0.6")
p5 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "RNA_snn_res.0.7")
plot_grid(p1, p2, p3,p4,p5, align = "hv", axis = "lrtb")


#######################################
DefaultAssay(scRNA.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA.integrated <- ScaleData(scRNA.integrated, verbose = FALSE)
scRNA.integrated <- RunPCA(scRNA.integrated, npcs  = 20, verbose = FALSE)
scRNA.integrated <- JackStraw(scRNA.integrated, num.replicate = 100)
scRNA.integrated <- ScoreJackStraw(scRNA.integrated, dims = 1:20)
JackStrawPlot(scRNA.integrated, dims = 1:20)
ElbowPlot(scRNA.integrated)

scRNA.integrated <- RunUMAP(scRNA.integrated, reduction = "pca", dims = 1:10) # , n.neighbors = 20
scRNA.integrated <- RunTSNE(scRNA.integrated, reduction = "pca", dims = 1:10)
scRNA.integrated <- FindNeighbors(scRNA.integrated, dims = 1:10)

scRNA.integrated <- FindClusters(scRNA.integrated, resolution = 0.2)  ## ", algorithm =4"  Leiden algorithm outperform a bit here, 3 cells misassign instead of 5
scRNA.integrated <- FindClusters(scRNA.integrated, resolution = 0.3)  ## ", algorithm =4"  Leiden algorithm outperform a bit here, 3 cells misassign instead of 5
scRNA.integrated <- FindClusters(scRNA.integrated, resolution = 0.4)  ## ", algorithm =4"  Leiden algorithm outperform a bit here, 3 cells misassign instead of 5
scRNA.integrated <- FindClusters(scRNA.integrated, resolution = 0.5)  ## ", algorithm =4"  Leiden algorithm outperform a bit here, 3 cells misassign instead of 5
scRNA.integrated <- FindClusters(scRNA.integrated, resolution = 0.7)  ## ", algorithm =4"  Leiden algorithm outperform a bit here, 3 cells misassign instead of 5
scRNA.integrated <- FindClusters(scRNA.integrated, resolution = 0.6)  ## ", algorithm =4"  Leiden algorithm outperform a bit here, 3 cells misassign instead of 5
scRNA.integrated <- FindClusters(scRNA.integrated, resolution = 0.1)  ## ", algorithm =4"  Leiden algorithm outperform a bit here, 3 cells misassign instead of 5

DimPlot(scRNA.integrated, reduction = "tsne")

### evaluate the clustering accuracy 
scRNA.integrated@meta.data$interated.ident <- Idents(scRNA.integrated)
ggplot(scRNA.integrated@meta.data,aes(x=interated.ident,fill=celltype_treatment))+
  geom_bar()+
  geom_text(aes(label=stat(count)),stat="count",position=position_stack(0.5)) +
  xlab("cluster")



### similar to Live-seq. ASPC of scRNA also need to be further clustered
scRNA.integrated_ASPC <- subset(scRNA.integrated, Cell_type == "ASPC")
scRNA.integrated_ASPC@active.assay
## normalization
scRNA.integrated_ASPC <- NormalizeData(scRNA.integrated_ASPC, verbose = FALSE)
## find MVGs
scRNA.integrated_ASPC <- FindVariableFeatures(scRNA.integrated_ASPC, selection.method = "vst",  nfeatures = 500, verbose = FALSE)
scRNA.integrated_ASPC <- CellCycleScoring(scRNA.integrated_ASPC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
## scale all the genes
scRNA.integrated_ASPC <- ScaleData(scRNA.integrated_ASPC, verbose = FALSE, features = rownames(scRNA.integrated_ASPC))
## run PCA
scRNA.integrated_ASPC <- RunPCA(scRNA.integrated_ASPC, npcs = 20, verbose = FALSE)
## evaluate PCA components
p1<- DimPlot(scRNA.integrated_ASPC, reduction = "pca", dims = c(1,2), group.by = "celltype_treatment")
p2<- DimPlot(scRNA.integrated_ASPC, reduction = "pca", dims = c(3,4), group.by = "celltype_treatment")
p3<- DimPlot(scRNA.integrated_ASPC, reduction = "pca", dims = c(5,6), group.by = "celltype_treatment")
p4<- DimPlot(scRNA.integrated_ASPC, reduction = "pca", dims = c(7,8), group.by = "celltype_treatment")
p5<- DimPlot(scRNA.integrated_ASPC, reduction = "pca", dims = c(9,10), group.by = "celltype_treatment")
p6<- DimPlot(scRNA.integrated_ASPC, reduction = "pca", dims = c(11,12), group.by = "celltype_treatment")
plot_grid(p1,p2,p3,p4,p5,p6)
## select pca accordingly for tsne and clustering
scRNA.integrated_ASPC <- RunTSNE(scRNA.integrated_ASPC, reduction = "pca", dims = 1:6, verbose = FALSE, seed.use = 1000)
scRNA.integrated_ASPC <- FindNeighbors(scRNA.integrated_ASPC, dims = 1:6 )  # , k.param =15
scRNA.integrated_ASPC <- FindClusters(scRNA.integrated_ASPC, resolution = 0.3)     ## resolution 0.3 
## also try cluster with k.mean
cluster.kmeans <- kmeans(t(scRNA.integrated_ASPC@assays$RNA@data[VariableFeatures(scRNA.integrated_ASPC) ,]),  2)
all(colnames(scRNA.integrated_ASPC) == names(cluster.kmeans$cluster ))
scRNA.integrated_ASPC$kmean.2 <- cluster.kmeans$cluster 

## plots
p1 <- DimPlot(scRNA.integrated_ASPC, reduction = "tsne")
p2 <- DimPlot(scRNA.integrated_ASPC, reduction = "tsne", group.by = "celltype_treatment")
p3 <- DimPlot(scRNA.integrated_ASPC, reduction = "tsne", group.by = "Batch")
p4<- FeaturePlot(scRNA.integrated_ASPC, reduction = "tsne", features  = "nFeature_RNA" )
p5 <- DimPlot(scRNA.integrated_ASPC, reduction = "tsne", group.by = "kmean.2")
plot_grid(p1,p2,p3, p4, p5)

## get the subcluster cell name
subcluster0 <- names(Idents(scRNA.integrated_ASPC) [Idents(scRNA.integrated_ASPC) == "0"])
subcluster1 <- names(Idents(scRNA.integrated_ASPC) [Idents(scRNA.integrated_ASPC) == "1"])

## rename the clusters, with subclustering result
scRNA.integrated$subcluster <- as.character(Idents(scRNA.integrated))
scRNA.integrated$subcluster[subcluster0] <- "4"
scRNA.integrated$subcluster[subcluster1] <- "5"

scRNA.integrated$subcluster <- factor(scRNA.integrated$subcluster )
DimPlot(scRNA.integrated, reduction = "tsne", group.by = "subcluster")

## reorder the level for plotting
levels(scRNA.integrated$subcluster) <- c(4,3,5,1,2)
scRNA.integrated$subcluster <- factor(scRNA.integrated$subcluster, levels = c(1,2,3,4,5))
DimPlot(scRNA.integrated, reduction = "tsne", group.by = "subcluster")
Idents(scRNA.integrated) <- scRNA.integrated$subcluster 


p1 <- DimPlot(scRNA.integrated, reduction = "umap")
p2 <- DimPlot(scRNA.integrated, reduction = "umap", group.by = "celltype_treatment")
p3 <- DimPlot(scRNA.integrated, reduction = "umap", group.by = "Batch")
p4 <- FeaturePlot(scRNA.integrated, reduction  = "umap", features = "nFeature_RNA")
plot_grid(p1, p2, p3,p4, align = "hv", axis = "lrtb")

p1 <- DimPlot(scRNA.integrated, reduction = "pca")
p2 <- DimPlot(scRNA.integrated, reduction = "pca", group.by = "celltype_treatment")
p3 <- DimPlot(scRNA.integrated, reduction = "pca", group.by = "Batch")
p4 <- FeaturePlot(scRNA.integrated, reduction  = "pca", features = "nFeature_RNA")
plot_grid(p1, p2, p3,p4, align = "hv", axis = "lrtb")

p1 <- DimPlot(scRNA.integrated, reduction = "tsne")
p2 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "celltype_treatment")
p3 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "Batch")
p4 <- FeaturePlot(scRNA.integrated, reduction  = "tsne", features = "nFeature_RNA")
p5 <- DimPlot(scRNA.integrated, reduction  = "tsne", group.by = "Phase")
plot_grid(p1, p2, p3,p4,p5, align = "hv", axis = "lrtb")


p1 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.3")
p2 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.4")
p3 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.5")
p4 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.6")
p5 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.7")
plot_grid(p1, p2, p3,p4,p5, align = "hv", axis = "lrtb")


saveRDS(scRNA.integrated, "3_scRNA_seq_only/scRNAseq_batch_corrected.rds")

# scRNA.integrated <-  readRDS("scRNAseq_batch_corrected.rds")


##################### DE genes  ##########
##################### Note that the edgR package are used to incorporate the batch effect, refer to the edgR analysis   ##########
## remove the pseudogene from DE analysis
geneName <- read.table(file = "mouseGeneTable87_mCherry_EGFP.txt", sep = "\t", header = T, row.names = 1 )
gene.pseudogene <- geneName[endsWith(as.character(geneName$gene_biotype),  "pseudogene"),]
gene.pseudoRemoved <- subset(scRNA@assays$RNA@meta.features, !(ensembl_gene_id %in% gene.pseudogene$ensembl_gene_id))


# For performing differential expression after integration, we switch back to the original
# data

DefaultAssay(scRNA.integrated) <- "integrated"

# # DE base on celltype_treatment
# Idents(scRNA.integrated) <- scRNA.integrated$celltype_treatment




scRNA.markers <- FindAllMarkers(scRNA.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(gene.pseudoRemoved))

# scRNA.markers.allGene <- FindAllMarkers(scRNA.integrated, only.pos = FALSE, min.pct = 0, logfc.threshold = 0)

# cluster2.markers <- FindMarkers(scRNA, ident.1 = 0,  ident.2 =2 , min.pct = 0.25, only.pos = T)
# cluster2.markers$gene <- rownames(cluster2.markers)

scRNA.markers$genesymbol <- scRNA@assays$RNA@meta.features[ scRNA.markers$gene , "external_gene_name"]
# scRNA.markers.allGene$genesymbol <- scRNA.integrated@assays$RNA@meta.features[ scRNA.markers.allGene$gene , "external_gene_name"]

# cluster2.markers$genesymbol <- scRNA@assays$RNA@meta.features[ cluster2.markers$gene , "external_gene_name"]


top10 <- scRNA.markers %>% subset(genesymbol != "mCherry" & genesymbol != "EGFP") %>% group_by(cluster) %>% top_n(n = -20, wt = p_val_adj)
top10 <- subset(top10, p_val_adj < 0.05)

DoHeatmap(scRNA.integrated, features = top10$gene)  + 
  scale_y_discrete(labels= rev(top10$genesymbol))

write.csv(scRNA.markers, "3_scRNA_seq_only/DEs.scRNA.csv")
# write.csv(scRNA.markers.allGene, "DEs.scRNA.allGenes.csv")




DefaultAssay(scRNA.integrated) <- "integrated"
p1 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "celltype_treatment" )
p2 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "Batch")

p3 <- FeaturePlot(scRNA.integrated, reduction = "tsne", features = symbol.to.ensembl("Zbtb16"))
p4 <- FeaturePlot(scRNA.integrated, reduction = "tsne", features = symbol.to.ensembl("Rora"))
p5 <- FeaturePlot(scRNA.integrated, reduction = "tsne", features = symbol.to.ensembl("Slc10a6"))
p6 <- FeaturePlot(scRNA.integrated, reduction = "tsne", features = symbol.to.ensembl("Dpep1"))
plot_grid(p1, p2,p3,p4,p5,p6,  align = "hv", axis = "lrtb")



