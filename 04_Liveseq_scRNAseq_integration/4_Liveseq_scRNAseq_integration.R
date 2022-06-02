library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(scater)
library(reshape2)
library(Matrix)
library(cowplot)
# Set working directory
setwd("~/SVFASRAW/wchen/data_analysis/Live_seq/final_analysis_V3/Code_github/")


Seu.all <- readRDS("1_preprocessing/Seu.all.RDS")

### subset of exp5,6,7,scRNA
Twodata <- subset(Seu.all, subset = ( (sampling_type %in% c("Live_seq", "scRNA")) & celltype_treatment %in% c("ASPC_not_treated", "ASPC_DMIR_treated", "IBA_not_treated", "Raw264.7_not_treated", "Raw264.7_LPS_treated") )) 
### remove 20 genes in black list, which are derived from the 0 pg input RNA negative control. 
gene.blacklist <- read.csv("gene.blacklist.csv")
data.count <- as.matrix(Twodata@assays$RNA@counts) 
data.count <- data.count[ !rownames(data.count) %in% gene.blacklist$ensembl_gene_id, ]

##### remove ribosomal protein genes
gene.ribo <- Twodata@assays$RNA@meta.features$ensembl_gene_id[grep("^Rp[s,l]", Twodata@assays$RNA@meta.features$external_gene_name) ]
gene.ribo <- as.character(gene.ribo)
data.count <- data.count[ !rownames(data.count) %in% gene.ribo, ]


Twodata <- CreateSeuratObject(counts =data.count, meta.data = Twodata@meta.data )

## add the meta feature
Twodata@assays$RNA@meta.features <- Twodata@assays$RNA@meta.features[ rownames(Twodata), ]

## add uniquely.mapped.rate and intron.mapped.rate in meta data
Twodata@meta.data$uniquely.mapped.rate <- Twodata@meta.data$uniquely.mapped / Twodata@meta.data$input.reads
Twodata@meta.data$intron.mapped.rate <- Twodata@meta.data$nCount_RNA / Twodata@meta.data$uniquely.mapped

## split object by sampling_type
Twodata.list <- SplitObject(Twodata, split.by = "sampling_type") 

Twodata.list[[1]]@meta.data <-  Twodata@meta.data[ colnames(Twodata.list[[1]]) , ]
Twodata.list[[2]]@meta.data <-  Twodata@meta.data[ colnames(Twodata.list[[2]]) , ]


for (i in 1:length(Twodata.list )) {
  
    ## remove unwanted cell based on the different criteria
  Twodata.list [[i]] <- subset(  Twodata.list [[i]], subset = nFeature_RNA > 1000 & percent.mt < 30 & uniquely.mapped.rate > 0.3)
  Twodata.list [[i]] <- NormalizeData(Twodata.list [[i]], verbose = FALSE)
  Twodata.list [[i]] <- FindVariableFeatures(Twodata.list [[i]], selection.method = "vst", nfeatures = 500, 
                                                 verbose = FALSE)
  ## remove EGFP and mCherry from the variable gene, as it biases.
  # VariableFeatures(Twodata.list [[i]]) <- VariableFeatures(Twodata.list [[i]]) [!VariableFeatures(Twodata.list [[i]]) %in% c("EGFP", "mCherry") ]
  
}




## integrate Live-seq and scRNA-seq data
reference.list <- Twodata.list[c("Live_seq", "scRNA")]
Twodata.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:10)
Twodata.integrated <- IntegrateData(anchorset = Twodata.anchors, dims = 1:10)

all(rownames(Twodata.integrated@assays$RNA@meta.features) == rownames( Twodata@assays$RNA@meta.features))
Twodata.integrated@assays$RNA@meta.features <- Twodata@assays$RNA@meta.features


## analysis on un-intergrated data
DefaultAssay(Twodata.integrated) <- "RNA"
Twodata.RNA <- NormalizeData(Twodata.integrated, verbose = FALSE)


# Twodata.RNA <- NormalizeData(Twodata.RNA, verbose = FALSE, normalization.method = "CLR")
### variable genes selection method "mean.var.plot"  

Twodata.RNA <- FindVariableFeatures(Twodata.RNA, selection.method = "vst", nfeatures = 500, verbose = FALSE)   

# Identify the  most highly variable genes
top30 <- head(VariableFeatures(Twodata.RNA), 30)
plot1 <- VariableFeaturePlot(Twodata.RNA)
plot2 <- LabelPoints(plot = plot1, points = top30,  labels = Twodata.RNA@assays$RNA@meta.features[top30, "external_gene_name"] , repel = TRUE)
plot2


Twodata.RNA <- ScaleData(Twodata.RNA, verbose = FALSE, features = rownames(Twodata.RNA))
Twodata.RNA <- RunPCA(Twodata.RNA, npcs = 50, verbose = FALSE)
Twodata.RNA <- JackStraw(Twodata.RNA, num.replicate = 100)
Twodata.RNA <- ScoreJackStraw(Twodata.RNA, dims = 1:20)
JackStrawPlot(Twodata.RNA, dims = 1:20)
ElbowPlot(Twodata.RNA)
Twodata.RNA <- RunTSNE(Twodata.RNA,reduction = "pca", dims = 1:10 , verbose = FALSE)
Twodata.RNA <- RunUMAP(Twodata.RNA, reduction = "pca", dims = 1:10,  n.neighbors = 5)
Twodata.RNA <- FindNeighbors(Twodata.RNA, dims = 1:10 )
Twodata.RNA <- FindClusters(Twodata.RNA, resolution = 0.1, verbose = T)     ## how stable the clusters are? 

p1 <- DimPlot(Twodata.RNA, reduction = "tsne", group.by = "sampling_type")
p2 <- DimPlot(Twodata.RNA, reduction = "tsne", group.by = "celltype_treatment", label = TRUE, repel = TRUE) + 
  NoLegend()
p3 <- FeaturePlot(Twodata.RNA, reduction  = "tsne", features = "nFeature_RNA")
p4 <- DimPlot(Twodata.RNA, reduction = "tsne", label = TRUE, repel = TRUE) + NoLegend()
plot_grid(p1, p2, p3,p4)



### order the cluster
levels(Idents(Twodata.RNA)) <- c(4, 3, 5, 1, 2)
Idents(Twodata.RNA) <- factor(Idents(Twodata.RNA), levels = c(1,2,3,4,5))

# Twodata.RNA@meta.data$celltype_treatment <- droplevels(Twodata.RNA@meta.data$celltype_treatment)
Twodata.RNA@meta.data$celltype_treatment <- factor(Twodata.RNA@meta.data$celltype_treatment, levels = c("ASPC_not_treated", "ASPC_DMIR_treated","IBA_not_treated","Raw264.7_not_treated", "Raw264.7_LPS_treated"))


### evaluate the cluster accuracy
Twodata.RNA@meta.data$RNA.ident <- Idents(Twodata.RNA)
ggplot(Twodata.RNA@meta.data,aes(x=RNA.ident,fill=celltype_treatment))+
  geom_bar(stat="count")+
  geom_text(aes(label=stat(count)),stat="count",position=position_stack(0.5)) +
  xlab("cluster")


p1 <- DimPlot(Twodata.RNA, reduction = "umap", group.by = "sampling_type")
p2 <- DimPlot(Twodata.RNA, reduction = "umap", group.by = "celltype_treatment", label = TRUE, repel = TRUE) + 
  NoLegend()
p3 <- FeaturePlot(Twodata.RNA, reduction  = "umap", features = "nFeature_RNA")
p4 <- DimPlot(Twodata.RNA, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
plot_grid(p1, p2, p3,p4)

p1 <- DimPlot(Twodata.RNA, reduction = "tsne", group.by = "sampling_type")
p2 <- DimPlot(Twodata.RNA, reduction = "tsne", group.by = "celltype_treatment", label = TRUE, repel = TRUE) + 
  NoLegend()
p3 <- FeaturePlot(Twodata.RNA, reduction  = "tsne", features = "nFeature_RNA")
p4 <- DimPlot(Twodata.RNA, reduction = "tsne", label = TRUE, repel = TRUE) + NoLegend()
plot_grid(p1, p2, p3,p4)


p1 <- DimPlot(Twodata.RNA, reduction = "pca", group.by = "sampling_type")
p2 <- DimPlot(Twodata.RNA, reduction = "pca", group.by = "celltype_treatment", label = TRUE, repel = TRUE) + 
  NoLegend()
p3 <- FeaturePlot(Twodata.RNA, reduction  = "pca", features = "nFeature_RNA", max.cutoff = 2000)
p4 <- DimPlot(Twodata.RNA, reduction = "pca", label = TRUE, repel = TRUE) + NoLegend()
plot_grid(p1, p2, p3,p4)



## mCherry and Tnf scatter plot


## functin for genesymbol and ensemble name conversation
gene.info <- read.table(file = "mouseGeneTable87_mCherry_EGFP.txt", sep = "\t", header = T, row.names = 1 )
symbol.to.ensembl <- function(x) {
  
  df <- subset(gene.info, external_gene_name == x) 
  gene <- as.character(df$ensembl_gene_id)
  gene <- gene[!duplicated(gene)]
  if(length(gene)>1) warning("more two ensemble names")
  return(gene[1])
}
ensembl.to.symbol <- function(x) {
  
  df <- subset(gene.info, ensembl_gene_id == x) 
  gene <- as.character(df$external_gene_name)
  gene <- gene[!duplicated(gene)]
  if(length(gene)>1) warning("more two gene symbol")
  return(gene[1])
} 



FeatureScatter(subset(Twodata.RNA, Cell_type == "Raw264.7_G9"), feature1 = "mCherry", feature2 = symbol.to.ensembl("Tnf") , group.by = "treatment")  + ylab("Tnf")
## mCherry and Tnf scatter plot of Live-seq sample only
FeatureScatter(subset(Twodata.RNA, Cell_type == "Raw264.7_G9" & sampling_type == "Live_seq"), feature1 = "mCherry", feature2 = symbol.to.ensembl("Tnf") , group.by = "treatment")  + ylab("Tnf")

VlnPlot(Twodata.RNA, features =symbol.to.ensembl("mCherry"), group.by = "mCherry_TNFa" ) + ggtitle("mCherry in all samples")
VlnPlot(subset(Twodata.RNA, treatment=="LPS_treated" ), features =symbol.to.ensembl("mCherry"), group.by = "mCherry_TNFa" ) + 
    ggtitle("mCherry LPS_treated")
VlnPlot(subset(Twodata.RNA, treatment=="not_treated" ), features =symbol.to.ensembl("mCherry"), group.by = "mCherry_TNFa" ) + 
    ggtitle("mCherry not_treated")



VlnPlot(Twodata.RNA, features =symbol.to.ensembl("Tnf"), group.by = "mCherry_TNFa" ) + ggtitle("Tnf in all samples")
VlnPlot(subset(Twodata.RNA, treatment=="LPS_treated" ), features =symbol.to.ensembl("Tnf"), group.by = "mCherry_TNFa" ) + 
    ggtitle("Tnf LPS_treated")
VlnPlot(subset(Twodata.RNA, treatment=="not_treated" ), features =symbol.to.ensembl("Tnf"), group.by = "mCherry_TNFa" ) + 
  ggtitle("Tnf not_treated")

saveRDS(Twodata.RNA, "4_Liveseq_scRNAseq_integration/Intergrated_data.RNA.rds")



#########  switch to integrated assay   #######


# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(Twodata.integrated) <- "integrated"


# Run the standard workflow for visualization and clustering
Twodata.integrated <- ScaleData(Twodata.integrated, verbose = FALSE)
Twodata.integrated <- RunPCA(Twodata.integrated, npcs  = 10, verbose = FALSE)
Twodata.integrated <- JackStraw(Twodata.integrated, num.replicate = 100)
Twodata.integrated <- ScoreJackStraw(Twodata.integrated, dims = 1:10)
JackStrawPlot(Twodata.integrated, dims = 1:10)
ElbowPlot(Twodata.integrated)

Twodata.integrated <- RunUMAP(Twodata.integrated, reduction = "pca", dims = 1:10)  ## , n.neighbors = 5
Twodata.integrated <- RunTSNE(Twodata.integrated, reduction = "pca", dims = 1:10, seed.use = 1)
Twodata.integrated <- FindNeighbors(Twodata.integrated, dims = 1:10)
Twodata.integrated <- FindClusters(Twodata.integrated, resolution = 0.2, verbose = T)  



p1 <- DimPlot(Twodata.integrated, reduction = "tsne", group.by = "sampling_type")
p2 <- DimPlot(Twodata.integrated, reduction = "tsne", group.by = "celltype_treatment")
p3 <- FeaturePlot(Twodata.integrated, reduction  = "tsne", features = "nFeature_RNA", max.cutoff = 5000)
p4 <- DimPlot(Twodata.integrated, reduction = "tsne")
plot_grid(p1, p2, p3,p4)


### order the cluster

levels(Idents(Twodata.integrated)) <- c(4, 5, 3, 1,2)
Idents(Twodata.integrated) <- factor(Idents(Twodata.integrated), levels = c(1,2,3,4,5))

# Twodata.integrated@meta.data$celltype_treatment <- droplevels(Twodata.integrated@meta.data$celltype_treatment)
Twodata.integrated@meta.data$celltype_treatment <- factor(Twodata.integrated@meta.data$celltype_treatment, levels = c("ASPC_not_treated","ASPC_DMIR_treated","IBA_not_treated","Raw264.7_not_treated", "Raw264.7_LPS_treated"))

### evaluate the clustering accuracy 
Twodata.integrated@meta.data$interated.ident <- Idents(Twodata.integrated)
ggplot(Twodata.integrated@meta.data,aes(x=interated.ident,fill=celltype_treatment))+
  geom_bar()+
  geom_text(aes(label=stat(count)),stat="count",position=position_stack(0.5)) +
  xlab("cluster")



p1 <- DimPlot(Twodata.integrated, reduction = "umap", group.by = "sampling_type")
p2 <- DimPlot(Twodata.integrated, reduction = "umap", group.by = "celltype_treatment")
p3 <- FeaturePlot(Twodata.integrated, reduction  = "umap", features = "nFeature_RNA")
p4 <- DimPlot(Twodata.integrated, reduction = "umap")

plot_grid(p1, p2, p3,p4)



p1 <- DimPlot(Twodata.integrated, reduction = "tsne", group.by = "sampling_type")
p2 <- DimPlot(Twodata.integrated, reduction = "tsne", group.by = "celltype_treatment")
p3 <- FeaturePlot(Twodata.integrated, reduction  = "tsne", features = "nFeature_RNA", max.cutoff = 5000)
p4 <- DimPlot(Twodata.integrated, reduction = "tsne")
plot_grid(p1, p2, p3,p4)


p1 <- DimPlot(Twodata.integrated, reduction = "pca", group.by = "sampling_type")
p2 <- DimPlot(Twodata.integrated, reduction = "pca", group.by = "celltype_treatment")
p3 <- FeaturePlot(Twodata.integrated, reduction  = "pca", features = "nFeature_RNA", max.cutoff = 5000)
p4 <- DimPlot(Twodata.integrated, reduction = "pca")
plot_grid(p1, p2, p3,p4)



#### plot genes of interest
### 
symbol.to.ensembl("Tnf")
ensembl.to.symbol("ENSMUSG00000024401")

genesToPlot <-  c("Grb14", "Thy1", "Gpx3","Sfrp2" ,"Mif","Tnf" )

plist <- lapply(genesToPlot, function(x) FeaturePlot(object = Twodata.integrated, reduction = "tsne", features = symbol.to.ensembl(x) ))
plot_grid(plotlist=plist)



## save integrated data
saveRDS(Twodata.integrated, "4_Liveseq_scRNAseq_integration/Intergrated_data.integrated.rds")


#### label the double extraction samples

Twodata.integrated$double_extraction <- factor(Twodata.integrated$double_extraction)

plist <- list()
for (i in levels(Twodata.integrated$double_extraction)) {
  cells.selected <- names( Twodata.integrated$double_extraction[ Twodata.integrated$double_extraction == i] ) ## select double extraction cells
  treatment <- Twodata.integrated$treatment[ Twodata.integrated$double_extraction == i]           ## get the LPS treatment info of double extraction cells
  double_extraction_order <- Twodata.integrated$double_extraction_order[ Twodata.integrated$double_extraction == i]  ## get the sampling order (1st or 2nd sampling)
  ilabel <- paste(cells.selected, treatment, double_extraction_order,    collapse = " ")
  ## plot double extraction sample
  # icolor <- rep("grey", ncol(Twodata))
  # icolor[Twodata$double_extraction == i & Twodata$double_extraction_order == 1 ] <- "#619CFF"
  # icolor[Twodata$double_extraction == i & Twodata$double_extraction_order == 2 ] <- "#00BA38"
  
  
  p<- DimPlot(Twodata.integrated, reduction = "tsne", group.by = "double_extraction", cells.highlight = cells.selected,  sizes.highlight = 1.5, pt.size = 0.3) + 
    ggtitle(ilabel ) 
  plist[[i]] <- p
}

p1 <- DimPlot(Twodata.integrated, reduction = "tsne" ,label = TRUE,label.size = 2 ,repel = T, pt.size = 0.3) + NoLegend()
plist[[1]] <- p1
plot_grid(plotlist = plist[1:16], align = "hv", axis = "lrtb")
plot_grid(plotlist = plist[17:32], align = "hv", axis = "lrtb")



dev.off()


## read data if needed
# Twodata.integrated <- readRDS("4_Liveseq_scRNAseq_integration/Intergrated_data.integrated.rds")

## Down_Sample_Matrix function
# if the sample size is smaill then expected depth, return the original sample 
Down_Sample_Matrix <- function (expr_mat, expected_depth)
{
  down_sample <- function(x, expected_depth) {
    prob <- expected_depth/sum(x)
    # if the sample size is smaill then expected depth, return the original sample 
    if (prob >1) { return(x)} 
    return(unlist(lapply(x, function(y) {
      rbinom(1, y, prob)
    })))
  }
  down_sampled_mat <- apply(expr_mat, 2, function(x)  down_sample(x,expected_depth ) )
  return(down_sampled_mat)
}

## pbmc.data.filter is the count matrix; Meta_data is the metadata of the count matrix; expected_depth in the number of reads need;
# cells with reads smaller than 0.8*expected_depth are removed
Down_Sample_Seurat_integrated <-    function( pbmc.data.filter, Meta_data , split_Feature , expected_depth)
{
  print(paste0("start downsampling to ",expected_depth, " reads", sep=""))
  pbmc.data.down <-  Down_Sample_Matrix( pbmc.data.filter,  expected_depth )
  print("downsampling matrix completed")
  # # remove the cells with read <= 0.8*expected_depth
  # pbmc.data.down <- pbmc.data.down[,  colSums(pbmc.data.down) > 0.8*expected_depth]
  
  # creat seurate object of the downsampled data
  pbmcDown <- CreateSeuratObject(counts = pbmc.data.down, min.cells  = 1, min.features = 1, project = "downsampling")
  
  ## check cells are in the Meta_data or not
  if(!all(rownames(pbmcDown@meta.data) == rownames(Meta_data) )) {
    warning("meta data does not match")
    if (!all(rownames(pbmcDown@meta.data) %in% rownames(Meta_data) )) { 
      stop("some cells are missing in Meta_data") }
    Meta_data <- Meta_data[rownames(pbmcDown@meta.data), ]
  }
  
  Meta_data$sample_name <- factor(Meta_data$sample_name, levels = Meta_data$sample_name)
  # change the order of the level of meta.data$sampling_type, for plotting
  # Meta_data$sampling_type <- factor(Meta_data$sampling_type, levels = c( "Negative_control", "Live_seq"   , "Bulk_raw_2ng", "Bulk_raw_500ng", "PC_Hela_5pg", "PC_IBA_1pg" , "PC_IBA_500pg"  , "PC_IBA_5pg" ,   "PC_Raw_5pg"  ,  "AFM"  ))
  Meta_data$sampling_type <- factor(Meta_data$sampling_type)
  Meta_data$sampling_type <- droplevels(Meta_data$sampling_type )
  # check the redundancy of the meta data, rename the original name with prefix "ori_'
  prefix <-  ifelse( colnames(Meta_data) %in% colnames(pbmcDown@meta.data), "ori_", "")
  print(prefix)
  colnames(Meta_data) <- paste(prefix, colnames(Meta_data), sep = "" )
  print(colnames(Meta_data))
  # combine the meta data
  print("start combine the meta data")
  pbmcDown@meta.data <- cbind ( pbmcDown@meta.data, Meta_data)
  pbmcDown@meta.data$cellType_condition <- paste(pbmcDown@meta.data$Cell_type, pbmcDown@meta.data$LPS_treatment, sep = "+")
  print("meta data combined sucessfully")
  
  
  
  pbmcDown.list <- SplitObject(pbmcDown, split.by = split_Feature) 
  
  pbmcDown.list[[1]]@meta.data <-  pbmcDown@meta.data[ colnames(pbmcDown.list[[1]]) , ]
  pbmcDown.list[[2]]@meta.data <-  pbmcDown@meta.data[ colnames(pbmcDown.list[[2]]) , ]
  
  
  for (i in 1:length(pbmcDown.list )) {
    
    ## remove unwanted cell based on the different criteria
    pbmcDown.list [[i]] <- subset(  pbmcDown.list [[i]], subset = nFeature_RNA > 1000 & percent.mt < 30 & uniquely.mapped.rate > 0.3)
    pbmcDown.list [[i]] <- NormalizeData(pbmcDown.list [[i]], verbose = FALSE)
    pbmcDown.list [[i]] <- FindVariableFeatures(pbmcDown.list [[i]], selection.method = "vst", nfeatures = 500, 
                                                   verbose = FALSE)
    ## remove EGFP and mCherry from the variable gene, as it biases.
    # VariableFeatures(Twodata.list [[i]]) <- VariableFeatures(Twodata.list [[i]]) [!VariableFeatures(Twodata.list [[i]]) %in% c("EGFP", "mCherry") ]
    
  }
  
  
  
  
  ## integrate Live-seq and scRNA-seq data
  reference.list <- pbmcDown.list[c("Live_seq", "scRNA")]
  pbmcDown.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:10)
  pbmcDown.integrated <- IntegrateData(anchorset = pbmcDown.anchors, dims = 1:10)
  
  all(rownames(pbmcDown.integrated@assays$RNA@meta.features) == rownames( pbmcDown@assays$RNA@meta.features))
  pbmcDown.integrated@assays$RNA@meta.features <- pbmcDown@assays$RNA@meta.features
  
  
  DefaultAssay(pbmcDown.integrated) <- "integrated"
  
  
  # Run the standard workflow for visualization and clustering
  pbmcDown.integrated <- ScaleData(pbmcDown.integrated, verbose = FALSE)
  pbmcDown.integrated <- RunPCA(pbmcDown.integrated, npcs  = 10, verbose = FALSE)
  # pbmcDown.integrated <- JackStraw(pbmcDown.integrated, num.replicate = 100)
  # pbmcDown.integrated <- ScoreJackStraw(pbmcDown.integrated, dims = 1:10)
  # JackStrawPlot(pbmcDown.integrated, dims = 1:10)
  # ElbowPlot(pbmcDown.integrated)
  
  pbmcDown.integrated <- RunUMAP(pbmcDown.integrated, reduction = "pca", dims = 1:10, n.neighbors = 5)
  pbmcDown.integrated <- RunTSNE(pbmcDown.integrated, reduction = "pca", dims = 1:10, seed.use = 2)
  pbmcDown.integrated <- FindNeighbors(pbmcDown.integrated, dims = 1:10)
  pbmcDown.integrated <- FindClusters(pbmcDown.integrated, resolution = 0.2, verbose = T)  ## Leiden algorithm outperform a bit here, 3 cells misassign instead of 5
  return(pbmcDown.integrated)
  
  # 
  # 
  # ### order the cluster
  # levels(pbmcDown.integrated) <- c(2,3,1,0)
  # Idents(Twodata.integrated) <- factor(Idents(Twodata.integrated), levels = c(0,1,2,3))
  
}


library(doSNOW)

nb_cores <- 4
cl <- makeCluster(nb_cores)
registerDoSNOW(cl)
sampling.size <- c(500000, 250000, 100000, 50000, 25000, 10000, 5000, 2500)

strt<-Sys.time()
## as the print() in the slave thread could not be shown. Creat a file to record the info from the slave threads. 
## Use the sink() function inside the foreach loop (slave thread) 
writeLines(c(""), "log.foreach.txt")
downsample.list <- foreach(i = 1:length(sampling.size), .packages = "Seurat", .export = c("Down_Sample_Matrix", "Down_Sample_Seurat_integrated"), .verbose=T) %dopar% {
  sink("log.foreach.txt", append=TRUE)
  Down_Sample_Seurat_integrated(Twodata.integrated@assays$RNA@counts, Twodata.integrated@meta.data, "sampling_type", sampling.size[i] )
}

print(Sys.time()-strt)
stopCluster(cl)

## plot the downsampled data

plist <- lapply( downsample.list, function(x) {
  x$idents <- Idents(x)
  p <- ggplot(x@meta.data,aes(x=idents,fill=celltype_treatment))+
    geom_bar()+
    geom_text(aes(label=stat(count)),stat="count",position=position_stack(0.5), size=3) +
    xlab("Cluster") +
    # scale_fill_manual(values =  c("#dfc27d", "lightgoldenrod4", "#80cdc1","#018571") )+
    ggtitle(mean(x$nCount_RNA))
  return(p)
})

plot_grid(plotlist = plist)



## plot the downsampled data

downsample.new <- downsample.list
for (i in 1:length(downsample.new)) {
  p <- DimPlot(downsample.new[[i]], reduction = "tsne", group.by = "celltype_treatment",label = TRUE,label.size = 2 ,repel = T, pt.size = 0.3) + 
    NoLegend() +
    ggtitle(mean(downsample.new[[i]]@meta.data$nCount_RNA))
}

plot_grid(plotlist = plist)

# save
saveRDS(downsample.list, "4_Liveseq_scRNAseq_integration/downsample.list.rds")




