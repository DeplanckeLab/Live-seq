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
library(viridis)
# source('~/NAS2/wchen/data_analysis/Resource/R_Funcitons/My_Rfunctions.R', local=TRUE)
# Set working directory
setwd("~/SVFASRAW/wchen/data_analysis/Live_seq/final_analysis_V3/Code_github/")
# save packages versions 
writeLines(capture.output(sessionInfo()), "2_Live_seq_only/sessionInfo.txt")


## make Live-seq object
# read all Seu.all
Seu.all <- readRDS("1_preprocessing/Seu.all.RDS")
# subset Liveseq
Liveseq <- subset(Seu.all, subset = (sampling_type == "Live_seq" & nFeature_RNA > 1000) )
dim(Liveseq)
# remove 20 genes in black list, which are derived from the 0 pg input RNA negative control. 
gene.blacklist <- read.csv("gene.blacklist.csv")
data.count <- as.matrix(Liveseq@assays$RNA@counts)
data.count <- data.count[ !rownames(data.count) %in% gene.blacklist$ensembl_gene_id, ]
# remove ribosomal protein genes
gene.ribo <- Liveseq@assays$RNA@meta.features$ensembl_gene_id[grep("^Rp[s,l]", Liveseq@assays$RNA@meta.features$external_gene_name) ]
gene.ribo <- as.character(gene.ribo)
data.count <- data.count[ !rownames(data.count) %in% gene.ribo, ]
dim(data.count)    
data.count <- data.count[ !rownames(data.count) %in% "EGFP", ]
dim(data.count)    ## 30081   512
# creat seurat object
Liveseq <- CreateSeuratObject(counts =data.count, meta.data = Liveseq@meta.data )
# add the meta feature
Liveseq@assays$RNA@meta.features <- Seu.all@assays$RNA@meta.features[ rownames(Liveseq), ]
# add uniquely.mapped.rate and intron.mapped.rate in meta data
Liveseq@meta.data$uniquely.mapped.rate <- Liveseq@meta.data$uniquely.mapped / Liveseq@meta.data$input.reads
Liveseq@meta.data$intron.mapped.rate <- Liveseq@meta.data$nCount_RNA / Liveseq@meta.data$uniquely.mapped
# Liveseq_sub@meta.data$celltype_treatment  <-  droplevels(Liveseq_sub@meta.data$celltype_treatment  )
Liveseq@meta.data$celltype_treatment  <- factor(Liveseq@meta.data$celltype_treatment , levels = c("ASPC_not_treated", "ASPC_DMIR_treated", "IBA_not_treated","Raw264.7_not_treated", "Raw264.7_LPS_treated"))




## QC 
VlnPlot(Liveseq, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rRNA", "percent.protein", "uniquely.mapped.rate", "intron.mapped.rate"), group.by = "celltype_treatment" , ncol = 3)
VlnPlot(Liveseq, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rRNA", "percent.protein", "uniquely.mapped.rate", "intron.mapped.rate"), group.by = "Batch" , ncol = 3)

# linear regression on input reads and nGene
fit.lm <- lm(nFeature_RNA ~ input.reads, Liveseq@meta.data)
summary(fit.lm)

getP.lm <- function(x) {
  summary.x <- summary(x)
  pf(summary.x$fstatistic[1],summary.x$fstatistic[2],summary.x$fstatistic[3],lower.tail=FALSE) 
}

ggplot(Liveseq@meta.data,aes(input.reads, nFeature_RNA)) +
        geom_point(size=0.1) + geom_smooth(method = lm) +
        labs(title = paste( round(fit.lm$coefficients[2]*1000000, 4) , "(1M scale up) x+", round(fit.lm$coefficients[1], 4), "\n R2=", summary(fit.lm)$adj.r.squared , "\n P(F test)", getP.lm(fit.lm)   )) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        scale_y_continuous(expand = c(0, 0)) +
        scale_x_continuous(expand = c(0, 0)) +
        xlab("Input reads") +
        ylab("nGene")





## data normalization 
Liveseq <- NormalizeData(Liveseq, verbose = FALSE)

## scale data
Liveseq <- FindVariableFeatures(Liveseq, selection.method = "vst",  nfeatures = 500, verbose = FALSE)

## Identify the 30 most highly variable genes
top30 <- head(VariableFeatures(Liveseq), 30)
plot1 <- VariableFeaturePlot(Liveseq)
plot2 <- LabelPoints(plot = plot1, points = top30,  labels = Liveseq@assays$RNA@meta.features[top30, "external_gene_name"] , repel = TRUE)
plot2

## cell cycle score
s.genes <- readRDS("s.genes.mouse.rds")
g2m.genes <- readRDS("g2m.genes.mouse.rds")

Liveseq <- CellCycleScoring(Liveseq, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


## scale all the genes
Liveseq <- ScaleData(Liveseq, verbose = FALSE, features = rownames(Liveseq))

## run PCA
Liveseq <- RunPCA(Liveseq, npcs = 20, verbose = FALSE)
Liveseq <- JackStraw(Liveseq, num.replicate = 100)
Liveseq <- ScoreJackStraw(Liveseq, dims = 1:20)
JackStrawPlot(Liveseq, dims = 1:20)
ElbowPlot(Liveseq)

## run tsne, umap 
Liveseq <- RunTSNE(Liveseq, reduction = "pca", dims = 1:10, verbose = FALSE, seed.use = 1000)
Liveseq <- RunUMAP(Liveseq, reduction = "pca", dims = 1:10)

## clustering
Liveseq <- FindNeighbors(Liveseq, dims = 1:10 )  
Liveseq <- FindClusters(Liveseq, resolution = 0.2)     ## resolution 0.3 


p1 <- DimPlot(Liveseq, reduction = "tsne")
p2 <- DimPlot(Liveseq, reduction = "tsne", group.by = "celltype_treatment")
p3 <- DimPlot(Liveseq, reduction = "tsne", group.by = "Batch")
p4<- FeaturePlot(Liveseq, reduction = "tsne", features  = "nFeature_RNA" )
plot_grid(p1,p2,p3, p4)

p1 <- DimPlot(Liveseq, reduction = "umap")
p2 <- DimPlot(Liveseq, reduction = "umap", group.by = "celltype_treatment")
p3 <- DimPlot(Liveseq, reduction = "umap", group.by = "Batch")
p4<- FeaturePlot(Liveseq, reduction = "umap", features  = "nFeature_RNA" )
plot_grid(p1,p2,p3,p4)


p1 <- DimPlot(Liveseq, reduction = "pca")
p2 <- DimPlot(Liveseq, reduction = "pca", group.by = "celltype_treatment")
p3 <- DimPlot(Liveseq, reduction = "pca", group.by = "Batch")
p4<- FeaturePlot(Liveseq, reduction = "pca", features  = "nFeature_RNA" )
plot_grid(p1,p2,p3,p4)


p1<- DimPlot(Liveseq, reduction = "pca", dims = c(3,4), group.by = "celltype_treatment")
p2<- DimPlot(Liveseq, reduction = "pca", dims = c(5,6), group.by = "celltype_treatment")
p3<- DimPlot(Liveseq, reduction = "pca", dims = c(7,8), group.by = "celltype_treatment")
p3<- DimPlot(Liveseq, reduction = "pca", dims = c(9,10), group.by = "celltype_treatment")
plot_grid(p1,p2,p3,p4)


### further clustering of ASPC
Liveseq_ASPC <- subset(Liveseq, Cell_type == "ASPC")
dim(Liveseq)
## normalization
Liveseq_ASPC <- NormalizeData(Liveseq_ASPC, verbose = FALSE)
## find MVGs
Liveseq_ASPC <- FindVariableFeatures(Liveseq_ASPC, selection.method = "vst",  nfeatures = 500, verbose = FALSE)
Liveseq_ASPC <- CellCycleScoring(Liveseq_ASPC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
## scale all the genes
Liveseq_ASPC <- ScaleData(Liveseq_ASPC, verbose = FALSE, features = rownames(Liveseq_ASPC))
## run PCA
Liveseq_ASPC <- RunPCA(Liveseq_ASPC, npcs = 20, verbose = FALSE)
## evaluate PCA components
p1<- DimPlot(Liveseq_ASPC, reduction = "pca", dims = c(1,2), group.by = "celltype_treatment")
p2<- DimPlot(Liveseq_ASPC, reduction = "pca", dims = c(3,4), group.by = "celltype_treatment")
p3<- DimPlot(Liveseq_ASPC, reduction = "pca", dims = c(5,6), group.by = "celltype_treatment")
p4<- DimPlot(Liveseq_ASPC, reduction = "pca", dims = c(7,8), group.by = "celltype_treatment")
p5<- DimPlot(Liveseq_ASPC, reduction = "pca", dims = c(9,10), group.by = "celltype_treatment")
p6<- DimPlot(Liveseq_ASPC, reduction = "pca", dims = c(11,12), group.by = "celltype_treatment")
plot_grid(p1,p2,p3,p4,p5,p6)
## select pca accordingly for tsne and clustering
Liveseq_ASPC <- RunTSNE(Liveseq_ASPC, reduction = "pca", dims = 1:8, verbose = FALSE, seed.use = 1000)
Liveseq_ASPC <- FindNeighbors(Liveseq_ASPC, dims = 1:8 )  # , k.param =15
Liveseq_ASPC <- FindClusters(Liveseq_ASPC, resolution = 0.3)     ## resolution 0.3 
## also try cluster with k.mean
cluster.kmeans <- kmeans(t(Liveseq_ASPC@assays$RNA@data[VariableFeatures(Liveseq_ASPC) ,]),  2)
all(colnames(Liveseq_ASPC) == names(cluster.kmeans$cluster ))
Liveseq_ASPC$kmean.2 <- cluster.kmeans$cluster 

## plots
p1 <- DimPlot(Liveseq_ASPC, reduction = "tsne")
p2 <- DimPlot(Liveseq_ASPC, reduction = "tsne", group.by = "celltype_treatment")
p3 <- DimPlot(Liveseq_ASPC, reduction = "tsne", group.by = "Batch")
p4<- FeaturePlot(Liveseq_ASPC, reduction = "tsne", features  = "nFeature_RNA" )
p5 <- DimPlot(Liveseq_ASPC, reduction = "tsne", group.by = "kmean.2")
plot_grid(p1,p2,p3, p4, p5)

## get the subcluster cell name
subcluster0 <- names(Idents(Liveseq_ASPC) [Idents(Liveseq_ASPC) == "0"])
subcluster1 <- names(Idents(Liveseq_ASPC) [Idents(Liveseq_ASPC) == "1"])

## rename the clusters, with subclustering result
Liveseq$subcluster <- as.character(Idents(Liveseq))
Liveseq$subcluster[subcluster0] <- "4"
Liveseq$subcluster[subcluster1] <- "5"

Liveseq$subcluster <- factor(Liveseq$subcluster )

levels(Liveseq$subcluster) <- c(4,5,3,2,1)
# levels(Liveseq$subcluster) <- c(2,1,3,4,5)

Liveseq$subcluster <- factor(Liveseq$subcluster, levels = c(1,2,3,4,5))
Idents(Liveseq) <- Liveseq$subcluster 


## plot
p1 <- DimPlot(Liveseq, reduction = "tsne", group.by = "subcluster")
p2 <- DimPlot(Liveseq, reduction = "tsne", group.by = "celltype_treatment")
p3 <- DimPlot(Liveseq, reduction = "tsne", group.by = "Batch")
plot_grid(p1,p2,p3)

## save object
saveRDS(Liveseq, "2_Live_seq_only/Liveseq.rds")



#########  DE genes per cluster  ##########
## remove the pseudogene from DE analysis
geneName <- read.table(file = "mouseGeneTable87_mCherry_EGFP.txt", sep = "\t", header = T, row.names = 1 )
gene.pseudogene <- geneName[endsWith(as.character(geneName$gene_biotype),  "pseudogene"),]
gene.pseudoRemoved <- subset(Liveseq@assays$RNA@meta.features, !(ensembl_gene_id %in% gene.pseudogene$ensembl_gene_id))

# find markers for every cluster compared to all remaining cells, report only the positive ones
Liveseq.markers <- FindAllMarkers(Liveseq, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(gene.pseudoRemoved))

# add gene symbol
Liveseq.markers$genesymbol <- Liveseq@assays$RNA@meta.features[ Liveseq.markers$gene , "external_gene_name"]

# heatmap with top 10 DE genes
top10.ASPC <- Liveseq.markers %>% subset(cluster==1 & avg_log2FC > 0 )%>% top_n(n = -10, wt = p_val_adj)
top10.ASPC_DMIR <- Liveseq.markers %>% subset(cluster==2 & avg_log2FC > 0 )%>% top_n(n = -10, wt = p_val_adj)
top10.IBA <- Liveseq.markers %>% subset(cluster==3 & avg_log2FC > 0 )%>% top_n(n = -10, wt = p_val_adj)
top10.RAW <- Liveseq.markers %>% subset(cluster==4 & avg_log2FC > 0 ) %>% top_n(n = -10, wt = p_val_adj)
top10.RAW_LPS <- Liveseq.markers %>% subset(cluster==5 & avg_log2FC > 0 )%>% top_n(n = -10, wt = p_val_adj)
top10 <- rbind(top10.ASPC,top10.ASPC_DMIR, top10.IBA,top10.RAW,top10.RAW_LPS  )
all(!duplicated(top10$genesymbol))
p <-DoHeatmap(Liveseq, features = top10$gene) + 
  scale_y_discrete(labels= rev(top10$genesymbol)) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) 
p

# save DE genes
write.csv(Liveseq.markers, "2_Live_seq_only/DEs.Liveseq_only.csv")
          
          


#########  DE genes per condition (celltype_treatment)  ##########          
# find markers for every cluster compared to all remaining cells, report only the positive ones

Idents(Liveseq) <- Liveseq$celltype_treatment
Liveseq_condition.markers <- FindAllMarkers(Liveseq, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, features = rownames(gene.pseudoRemoved))

# add genesymbol
Liveseq_condition.markers$genesymbol <- Liveseq@assays$RNA@meta.features[ Liveseq$gene , "external_gene_name"]

# plot heatmap of top 10 DE genes
top10.ASPC <- Liveseq_condition.markers %>% subset(cluster=="ASPC_not_treated" & avg_log2FC > 0 )%>% top_n(n = -10, wt = p_val_adj)
top10.ASPC_DMIR <- Liveseq_condition.markers %>% subset(cluster=="ASPC_DMIR_treated" & avg_log2FC > 0 )%>% top_n(n = -10, wt = p_val_adj)
top10.IBA <- Liveseq_condition.markers %>% subset(cluster=="IBA_not_treated" & avg_log2FC > 0 )%>% top_n(n = -10, wt = p_val_adj)
top10.RAW <- Liveseq_condition.markers %>% subset(cluster=="Raw264.7_not_treated" & avg_log2FC > 0 ) %>% top_n(n = -10, wt = p_val_adj)
top10.RAW_LPS <- Liveseq_condition.markers %>% subset(cluster=="Raw264.7_LPS_treated" & avg_log2FC > 0 )%>% top_n(n = -10, wt = p_val_adj)
top10 <- rbind(top10.ASPC,top10.ASPC_DMIR, top10.IBA,top10.RAW,top10.RAW_LPS  )
all(!duplicated(top10$genesymbol))

p <-DoHeatmap(Liveseq, features = top10$gene) + 
  scale_y_discrete(labels= rev(top10$genesymbol)) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) 
p
# save DE 
write.csv(Liveseq_condition.markers, "2_Live_seq_only/DEs.Liveseq_per_celltype_treatment.csv")




############### down sample the Live-seq data ############### 

## Down_Sample_Matrix function
# if the sample size is smaill then expected depth, return the original sample 
Down_Sample_Matrix <- function (expr_mat, expected_depth)
{
  # down sample function of each colume
  down_sample <- function(x, expected_depth) {
    # establish the probability
    prob <- expected_depth/sum(x)
    # if the sample size is smaill then expected depth, return the original sample 
    if (prob >1) { return(x)} 
    return(unlist(lapply(x, function(y) {
      rbinom(1, y, prob)
    })))
  }
  # down sample the whole matirx
  down_sampled_mat <- apply(expr_mat, 2, function(x)  down_sample(x,expected_depth ) )
  return(down_sampled_mat)
}

## pbmc.data.filter is the count matrix; Meta_data is the metadata of the count matrix; expected_depth in the number of reads need;
# cells with reads smaller than 0.8*expected_depth are removed
Down_Sample_Seurat <-    function( pbmc.data.filter, Meta_data , expected_depth)
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
  Meta_data$sampling_type <- factor(Meta_data$sampling_type )
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
  # Normalizing the data
  pbmcDown <- NormalizeData(pbmcDown, verbose = FALSE)

  # Detection and plot the variable genes across the single cells
  pbmcDown <- FindVariableFeatures(pbmcDown, selection.method = "vst", nfeatures = 500, 
                                      verbose = FALSE)

  # pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))
  pbmcDown <- ScaleData(pbmcDown, verbose = FALSE, features = rownames(pbmcDown))
  pbmcDown <- RunPCA(pbmcDown, npcs = 20, verbose = FALSE)
  # pbmcDown <- JackStraw(pbmcDown, num.replicate = 100)
  # pbmcDown <- ScoreJackStraw(pbmcDown, dims = 1:20)
  # JackStrawPlot(pbmcDown, dims = 1:20)
  # ElbowPlot(pbmcDown)
  pbmcDown <- RunTSNE(pbmcDown, reduction = "pca", dims = 1:10, verbose = FALSE, perplexity =15)
  pbmcDown <- RunUMAP(pbmcDown, reduction = "pca", dims = 1:10,  n.neighbors = 5)
  
  pbmcDown <- FindNeighbors(pbmcDown, dims = 1:10, k.param = 15 )
  pbmcDown <- FindClusters(pbmcDown, resolution = 0.25)     ## resolution 0.3 
  return(pbmcDown)
}

## use tue doSNOW for parallel processing
library(doSNOW)
# assign core number and register
nb_cores <- 5
cl <- makeCluster(nb_cores)
registerDoSNOW(cl)
sampling.size <- c(500000, 250000, 100000, 50000, 25000, 10000, 5000, 2500, 1000, 500)

# recode starting time
strt<-Sys.time()
## as the print() in the slave thread could not be shown. Creat a file to record the info from the slave threads. 
## Use the sink() function inside the foreach loop (slave thread) 
writeLines(c(""), "log.foreach.txt")
downsample.list <- foreach(i = 1:length(sampling.size), .packages = "Seurat", .export = c("Down_Sample_Matrix", "Down_Sample_Seurat"), .verbose=T) %dopar% {
  sink("log.foreach.txt", append=TRUE)
  Down_Sample_Seurat(Liveseq@assays$RNA@counts[, colnames(Liveseq)], Liveseq@meta.data, sampling.size[i] )
}

print(Sys.time()-strt)
stopCluster(cl)

# check some downsample results

p1 <- DimPlot(downsample.list[[2]], reduction = "tsne" )
p2 <- DimPlot(downsample.list[[2]] ,reduction = "tsne" , group.by = "celltype_treatment")
plot_grid(p1, p2)

p1 <- DimPlot(downsample.list[[5]], reduction = "tsne" )
p2 <- DimPlot(downsample.list[[5]] ,reduction = "tsne" , group.by = "celltype_treatment")
plot_grid(p1, p2)

## refine clustering if needed
# downsample.list[[1]] <-   FindClusters(downsample.list[[1]], resolution = 0.3)     ## resolution 0.3 
# downsample.list[[2]] <-   FindClusters(downsample.list[[2]], resolution = 0.3)     ## resolution 0.3 
# downsample.list[[3]] <-   FindClusters(downsample.list[[3]], resolution = 0.3)     ## resolution 0.3 
# downsample.list[[4]] <-   FindClusters(downsample.list[[4]], resolution = 0.3)     ## resolution 0.3 
# downsample.list[[5]] <-   FindClusters(downsample.list[[5]], resolution = 0.3)     ## resolution 0.3 

downsample.new <- downsample.list

pdf("2_Live_seq_only/downsampling.intermediate.pdf", width = 10)
for (i in 1:length(downsample.new)) {
  objectToPlot <- downsample.new[[i]] 
  # p1 <- DimPlot(objectToPlot, reduction = "umap")
  # p2 <- DimPlot(objectToPlot, reduction = "umap", group.by = "celltype_treatment")
  # print(plot_grid(p1,p2))

  p1 <- DimPlot(objectToPlot, reduction = "tsne")
  p2 <- DimPlot(objectToPlot, reduction = "tsne", group.by = "celltype_treatment")
  print(plot_grid(p1,p2))
  # 
  # p1 <- DimPlot(objectToPlot, reduction = "pca")
  # p2 <- DimPlot(objectToPlot, reduction = "pca", group.by = "celltype_treatment")
  # print(plot_grid(p1,p2))
  
  objectToPlot@meta.data$ident <- Idents(objectToPlot)
  p <- ggplot(objectToPlot@meta.data,aes(x=ident,fill=celltype_treatment))+
    geom_bar()+
    geom_text(aes(label=stat(count)),stat="count",position=position_stack(0.5), size=3) +
    xlab("cluster")
  print(p)
  # color.celltype <- c("#dfc27d", "lightgoldenrod4", "#80cdc1","#018571") # 
  
  plist <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rRNA", "percent.protein", "uniquely.mapped.rate", "intron.mapped.rate"), function(x) VlnPlot(object = objectToPlot, features = x,group.by ="celltype_treatment"))
  print(plot_grid(plotlist = plist))

}

dev.off()

saveRDS(downsample.new, "2_Live_seq_only/downsample.new.rds")
# downsample.list <- readRDS("downsample.list.rds")
