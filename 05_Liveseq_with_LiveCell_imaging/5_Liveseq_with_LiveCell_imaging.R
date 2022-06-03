################################################################
#                                                              #
#               Live-seq with LiveCell imaging                 #
#                                                              #
################################################################

### Questions: https://github.com/DeplanckeLab/Live-seq/issues
### Date: 2022-03-06
### Datasets: Live-seq and scRNA-seq only
### Goal: Use Live-seq to predict how a cell will respond to TNF stimulation

root_dir <- rprojroot::find_root(rprojroot::is_rstudio_project)

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

## load Seu.all
Seu.all <- readRDS(file.path(root_dir, "01_preprocessing/Seu.all.rds"))

## function for genesymbol and ensemble name conversation
gene.info <-
  read.table(
    file = paste0(root_dir, "/data/mouseGeneTable87_mCherry_EGFP.txt"),
    sep = "\t",
    header = T,
    row.names = 1
  )
symbol.to.ensembl <- function(x) {
  df <- subset(gene.info, external_gene_name == x)
  gene <- as.character(df$ensembl_gene_id)
  gene <- gene[!duplicated(gene)]
  if (length(gene) > 1)
    warning("more two ensemble names")
  return(gene[1])
}
ensembl.to.symbol <- function(x) {
  df <- subset(gene.info, ensembl_gene_id == x)
  gene <- as.character(df$external_gene_name)
  gene <- gene[!duplicated(gene)]
  if (length(gene) > 1)
    warning("more two gene symbol")
  return(gene[1])
}

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######

###### find the most variable genes from Raw cell without LPS treatment, both from Live-seq and scRNA-seq

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######

## Liveseq data of RAW cells not treated
## Batch 9_4 are experimental a bit distinct from others (new batches of primers and enzymes), but can also be included with overall consistent result
Liveseq_G9Mock <-
  subset(Seu.all, subset = ((sampling_type == "Live_seq") &
                              (treatment == "not_treated") &
                              (Cell_type == "Raw264.7_G9") &
                              (Batch %in% c("5Reseq", "6",  "7", "scRNA", "8_8"))
  ))
## add uniquely.mapped.rate and intron.mapped.rate in meta data
Liveseq_G9Mock@meta.data$uniquely.mapped.rate <-
  Liveseq_G9Mock@meta.data$uniquely.mapped / Liveseq_G9Mock@meta.data$input.reads
Liveseq_G9Mock@meta.data$intron.mapped.rate <-
  Liveseq_G9Mock@meta.data$nCount_RNA / Liveseq_G9Mock@meta.data$uniquely.mapped
Liveseq_G9Mock <-
  subset(
    x = Liveseq_G9Mock,
    subset = nFeature_RNA > 1000 &
      percent.mt < 30 &
      uniquely.mapped.rate > 0.3 &
      percent.rRNA < 30 & intron.mapped.rate > 0.5
  )


### remove 211 genes in black list, which are derived from the 0 pg input RNA negative control.
gene.blacklist <- read.csv(paste0(root_dir, "/data/gene.blacklist.csv"))
data.count <- as.matrix(Liveseq_G9Mock@assays$RNA@counts)
data.count <-
  data.count[!rownames(data.count) %in% gene.blacklist$ensembl_gene_id,]

##### remove ribosomal protein genes
gene.ribo <-
  Liveseq_G9Mock@assays$RNA@meta.features$ensembl_gene_id[grep("^Rp[s,l]",
                                                               Liveseq_G9Mock@assays$RNA@meta.features$external_gene_name)]
gene.ribo <- as.character(gene.ribo)
data.count <- data.count[!rownames(data.count) %in% gene.ribo,]

### creat seurat object
Liveseq_G9Mock <-
  CreateSeuratObject(counts = data.count, meta.data = Liveseq_G9Mock@meta.data)

Liveseq_G9Mock@assays$RNA@meta.features <-
  Seu.all@assays$RNA@meta.features[rownames(Liveseq_G9Mock),]



## QC
VlnPlot(
  Liveseq_G9Mock,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt",
    "percent.rRNA",
    "percent.protein",
    "uniquely.mapped.rate",
    "intron.mapped.rate"
  ),
  group.by = "celltype_treatment" ,
  ncol = 3
)


plist <-
  lapply(c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt",
    "percent.rRNA",
    "percent.protein",
    "uniquely.mapped.rate",
    "intron.mapped.rate",
    "input.reads"
  ), function(i) {
    # p <- myVlnPlot(Liveseq_sub, x, "sampling_type", cols = c("grey","#00BA38","#619CFF"))
    p <-
      ggplot(Liveseq_G9Mock@meta.data,
             aes_string(x = "sampling_type", y = i)) +
      geom_violin(trim = FALSE,  fill = "gray") +
      stat_summary(
        fun.data = mean_sdl,
        fun.args = list(mult = 1),
        geom = "pointrange",
        color = viridis(10)[1]
      ) +
      theme_classic()
    p <- p + geom_jitter(shape = 16, position = position_jitter(0.2))
    # p <- p+ geom_boxplot(width = 0.2)
    # p<- p + geom_violin(draw_quantiles = 0.5, alpha = 0)  ## alpha for transparent
    # p<- p + stat_summary(fun.y=mean, geom="point",  size=2, color="red")
    return(p)
  })

plot_grid(
  plotlist = plist,
  ncol = 8,
  align = "hv",
  axis = "tblr"
)
# ggsave("05_Liveseq_with_LiveCell_imaging/QC_noGroup_before_filtering_Liveseq_G9Mock.pdf",width = 10, height = 1.6, useDingbats=FALSE )


## data normalization
Liveseq_G9Mock <- NormalizeData(Liveseq_G9Mock, verbose = FALSE)
## scale data
Liveseq_G9Mock <-
  FindVariableFeatures(
    Liveseq_G9Mock,
    selection.method = "vst",
    nfeatures = 500,
    verbose = FALSE
  )

MVG_Liveseq_G9Mock <-
  data.frame(ensembl_ID = (VariableFeatures(Liveseq_G9Mock)))
MVG_Liveseq_G9Mock$gene_symbol <-
  sapply(as.character(MVG_Liveseq_G9Mock$ensembl_ID),
         ensembl.to.symbol)
write.csv(MVG_Liveseq_G9Mock,
          paste0(root_dir, "/05_Liveseq_with_LiveCell_imaging/MVG_Liveseq_G9Mock.csv"))

# Identify the  most highly variable genes
top30 <- head(VariableFeatures(Liveseq_G9Mock), 30)
plot1 <- VariableFeaturePlot(Liveseq_G9Mock)
plot2 <-
  LabelPoints(
    plot = plot1,
    points = top30,
    labels = Liveseq_G9Mock@assays$RNA@meta.features[top30, "external_gene_name"] ,
    repel = TRUE
  )
plot2

## scale all the genes
Liveseq_G9Mock <-
  ScaleData(Liveseq_G9Mock,
            verbose = FALSE,
            features = rownames(Liveseq_G9Mock))

# Liveseq_sub <- ScaleData(Liveseq_sub, verbose = FALSE,vars.to.regress = c("nFeature_RNA") )

Liveseq_G9Mock <- RunPCA(Liveseq_G9Mock, npcs = 20, verbose = FALSE)
Liveseq_G9Mock <- JackStraw(Liveseq_G9Mock, num.replicate = 100)
Liveseq_G9Mock <- ScoreJackStraw(Liveseq_G9Mock, dims = 1:20)
JackStrawPlot(Liveseq_G9Mock, dims = 1:20)
ElbowPlot(Liveseq_G9Mock)


Liveseq_G9Mock <-
  RunTSNE(
    Liveseq_G9Mock,
    reduction = "pca",
    dims = 1:10,
    verbose = FALSE,
    perplexity = 15,
    seed.use = 3
  )
Liveseq_G9Mock <-
  RunUMAP(
    Liveseq_G9Mock,
    reduction = "pca",
    dims = 1:10,
    n.neighbors = 5
  )

Liveseq_G9Mock <-
  FindNeighbors(Liveseq_G9Mock, dims = 1:10, k.param = 15)
Liveseq_G9Mock <-
  FindClusters(Liveseq_G9Mock, resolution = 0.3)     ## resolution 0.3

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- readRDS("s.genes.mouse.rds")
g2m.genes <- readRDS("g2m.genes.mouse.rds")
# add cell cycle score
Liveseq_G9Mock <-
  CellCycleScoring(
    Liveseq_G9Mock,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE
  )

# check cell cycle and mCherry.log.slope
FeatureScatter(Liveseq_G9Mock, feature1 = "S.Score", feature2 = "mCherry.log.slope")
# ggsave("05_Liveseq_with_LiveCell_imaging/S.score_mCherry.slope.pdf")

# check cell cycle and mCherry.log.slope
FeatureScatter(Liveseq_G9Mock, feature1 = "G2M.Score", feature2 = "mCherry.log.slope")
# ggsave("05_Liveseq_with_LiveCell_imaging/G2M.score_mCherry.slope.pdf")


# dimplot with cell cycle
p1 <-
  DimPlot(Liveseq_G9Mock, reduction = "tsne", group.by = "Phase")
p2 <-
  FeaturePlot(Liveseq_G9Mock, reduction = "tsne", features = "S.Score")
p3 <-
  FeaturePlot(Liveseq_G9Mock, reduction = "tsne", features = "G2M.Score")
p4 <- DimPlot(Liveseq_G9Mock, reduction = "pca", group.by = "Phase")
p5 <-
  FeaturePlot(
    Liveseq_G9Mock,
    reduction = "pca",
    dims = c(1, 2),
    features = "S.Score"
  )
p6 <-
  FeaturePlot(
    Liveseq_G9Mock,
    reduction = "pca",
    dims = c(1, 2),
    features = "G2M.Score"
  )

plot_grid(p1, p2, p3, p4, p5, p6)
# ggsave("05_Liveseq_with_LiveCell_imaging/Liveseq_G9Mock_cellCycle.pdf", width = 10, height = 6)


p1 <- DimPlot(Liveseq_G9Mock, reduction = "tsne")
p2 <-
  DimPlot(Liveseq_G9Mock, reduction = "tsne", group.by = "celltype_treatment")
plot_grid(p1, p2)

p1 <- DimPlot(Liveseq_G9Mock, reduction = "pca")
p2 <-
  DimPlot(Liveseq_G9Mock, reduction = "pca", group.by = "celltype_treatment")
plot_grid(p1, p2)

# p <- FeaturePlot(Liveseq_G9Mock, features = symbol.to.ensembl("Tnf"))
# HoverLocator(p)

# save object
saveRDS(Liveseq_G9Mock,
        "05_Liveseq_with_LiveCell_imaging/Liveseq_G9Mock.rds")
# Liveseq_G9Mock <- readRDS("Liveseq_G9Mock.rds")




#####################
## scRNA data of RAW cells not treated
#####################

## Batch 9_4 are experimental a bit distinct from others (new batches of primers and enzymes), but can also be included with overall consistent result
scRNA_G9Mock <-
  subset(
    Seu.all,
    subset =  (sampling_type == "scRNA") &
      (treatment == "not_treated") &
      (Cell_type == "Raw264.7_G9") &
      (Batch %in% c("5Reseq", "6",  "7", "scRNA", "8_8"))
  )
## add uniquely.mapped.rate and intron.mapped.rate in meta data
scRNA_G9Mock@meta.data$uniquely.mapped.rate <-
  scRNA_G9Mock@meta.data$uniquely.mapped / scRNA_G9Mock@meta.data$input.reads
scRNA_G9Mock@meta.data$intron.mapped.rate <-
  scRNA_G9Mock@meta.data$nCount_RNA / scRNA_G9Mock@meta.data$uniquely.mapped
scRNA_G9Mock <-
  subset(
    x = scRNA_G9Mock,
    subset = nFeature_RNA > 1000 &
      percent.mt < 30 &
      uniquely.mapped.rate > 0.3 &
      percent.rRNA < 30 & intron.mapped.rate > 0.5
  )


### remove 211 genes in black list, which are derived from the 0 pg input RNA negative control.
gene.blacklist <- read.csv(file.path(root_dir, "data/gene.blacklist.csv"))
data.count <- as.matrix(scRNA_G9Mock@assays$RNA@counts)
data.count <-
  data.count[!rownames(data.count) %in% gene.blacklist$ensembl_gene_id,]

##### remove ribosomal protein genes
gene.ribo <-
  scRNA_G9Mock@assays$RNA@meta.features$ensembl_gene_id[grep("^Rp[s,l]",
                                                             scRNA_G9Mock@assays$RNA@meta.features$external_gene_name)]
gene.ribo <- as.character(gene.ribo)
data.count <- data.count[!rownames(data.count) %in% gene.ribo,]

### creat seurat object
scRNA_G9Mock <-
  CreateSeuratObject(counts = data.count, meta.data = scRNA_G9Mock@meta.data)


scRNA_G9Mock@assays$RNA@meta.features <-
  Seu.all@assays$RNA@meta.features[rownames(scRNA_G9Mock),]

## add uniquely.mapped.rate and intron.mapped.rate in meta data
scRNA_G9Mock@meta.data$uniquely.mapped.rate <-
  scRNA_G9Mock@meta.data$uniquely.mapped / scRNA_G9Mock@meta.data$input.reads
scRNA_G9Mock@meta.data$intron.mapped.rate <-
  scRNA_G9Mock@meta.data$nCount_RNA / scRNA_G9Mock@meta.data$uniquely.mapped

## QC
VlnPlot(
  scRNA_G9Mock,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt",
    "percent.rRNA",
    "percent.protein",
    "uniquely.mapped.rate",
    "intron.mapped.rate"
  ),
  group.by = "celltype_treatment" ,
  ncol = 3
)


plist <-
  lapply(c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt",
    "percent.rRNA",
    "percent.protein",
    "uniquely.mapped.rate",
    "intron.mapped.rate",
    "input.reads"
  ), function(i) {
    # p <- myVlnPlot(Liveseq_sub, x, "sampling_type", cols = c("grey","#00BA38","#619CFF"))
    p <-
      ggplot(scRNA_G9Mock@meta.data, aes_string(x = "sampling_type", y = i)) +
      geom_violin(trim = FALSE,  fill = "gray") +
      stat_summary(
        fun.data = mean_sdl,
        fun.args = list(mult = 1),
        geom = "pointrange",
        color = viridis(10)[1]
      ) +
      theme_classic()
    p <- p + geom_jitter(shape = 16, position = position_jitter(0.2))
    # p <- p+ geom_boxplot(width = 0.2)
    # p<- p + geom_violin(draw_quantiles = 0.5, alpha = 0)  ## alpha for transparent
    # p<- p + stat_summary(fun.y=mean, geom="point",  size=2, color="red")
    return(p)
  })

plot_grid(
  plotlist = plist,
  ncol = 8,
  align = "hv",
  axis = "tblr"
)

# ggsave("05_Liveseq_with_LiveCell_imaging/QC_noGroup_scRNA_G9Mock.pdf",width = 10, height = 1.6, useDingbats=FALSE )



## data normalization
scRNA_G9Mock <- NormalizeData(scRNA_G9Mock, verbose = FALSE)
## scale data
scRNA_G9Mock <-
  FindVariableFeatures(
    scRNA_G9Mock,
    selection.method = "vst",
    nfeatures = 500,
    verbose = FALSE
  )

MVG_scRNA_G9Mock <-
  data.frame(ensembl_ID = (VariableFeatures(scRNA_G9Mock)))
MVG_scRNA_G9Mock$gene_symbol <-
  sapply(as.character(MVG_scRNA_G9Mock$ensembl_ID),
         ensembl.to.symbol)
write.csv(MVG_scRNA_G9Mock,
          "05_Liveseq_with_LiveCell_imaging/MVG_scRNA_G9Mock.csv")

## remove EGFP and mCherry from the variable gene, as it biases.
# VariableFeatures(Liveseq_sub) <- VariableFeatures(Liveseq_sub) [!VariableFeatures(Liveseq_sub) %in% c("EGFP", "mCherry") ]
# Identify the 10 most highly variable genes
top30 <- head(VariableFeatures(scRNA_G9Mock), 30)
top30 <- c(top30, "EGFP", "ENSMUSG00000024927")

plot1 <- VariableFeaturePlot(scRNA_G9Mock)
plot2 <-
  LabelPoints(
    plot = plot1,
    points = top30,
    labels = as.character(scRNA_G9Mock@assays$RNA@meta.features[top30, "external_gene_name"]) ,
    repel = TRUE
  )
plot2

## scale all the genes
scRNA_G9Mock <-
  ScaleData(scRNA_G9Mock,
            verbose = FALSE,
            features = rownames(scRNA_G9Mock))

# Liveseq_sub <- ScaleData(Liveseq_sub, verbose = FALSE,vars.to.regress = c("nFeature_RNA") )

scRNA_G9Mock <- RunPCA(scRNA_G9Mock, npcs = 20, verbose = FALSE)
scRNA_G9Mock <- JackStraw(scRNA_G9Mock, num.replicate = 100)
scRNA_G9Mock <- ScoreJackStraw(scRNA_G9Mock, dims = 1:20)
JackStrawPlot(scRNA_G9Mock, dims = 1:20)
ElbowPlot(scRNA_G9Mock)


scRNA_G9Mock <-
  RunTSNE(
    scRNA_G9Mock,
    reduction = "pca",
    dims = 1:10,
    verbose = FALSE,
    perplexity = 15,
    seed.use = 3
  )
scRNA_G9Mock <-
  RunUMAP(
    scRNA_G9Mock,
    reduction = "pca",
    dims = 1:10,
    n.neighbors = 5
  )

scRNA_G9Mock <-
  FindNeighbors(scRNA_G9Mock, dims = 1:10, k.param = 15)
scRNA_G9Mock <-
  FindClusters(scRNA_G9Mock, resolution = 0.3)     ## resolution 0.3



p1 <- DimPlot(scRNA_G9Mock, reduction = "tsne")
p2 <-
  DimPlot(scRNA_G9Mock, reduction = "tsne", group.by = "celltype_treatment")
plot_grid(p1, p2)

p1 <- DimPlot(scRNA_G9Mock, reduction = "pca")
p2 <-
  DimPlot(scRNA_G9Mock, reduction = "pca", group.by = "celltype_treatment")
plot_grid(p1, p2)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase


scRNA_G9Mock <-
  CellCycleScoring(
    scRNA_G9Mock,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE
  )

p1 <- DimPlot(scRNA_G9Mock, reduction = "tsne")
p2 <- DimPlot(scRNA_G9Mock, reduction = "tsne", group.by = "Phase")
plot_grid(p1, p2)

p1 <- DimPlot(scRNA_G9Mock, reduction = "pca")
p2 <- DimPlot(scRNA_G9Mock, reduction = "pca", group.by = "Phase")
plot_grid(p1, p2)

# HoverLocator(p1)

saveRDS(scRNA_G9Mock,
        "05_Liveseq_with_LiveCell_imaging/scRNA_G9Mock.rds")


## save the MVG of both scRNA-seq and Live-seq
MVG_scRNA_G9Mock <-
  read.csv("05_Liveseq_with_LiveCell_imaging/MVG_scRNA_G9Mock.csv",
           row.names = 1)
dim(MVG_scRNA_G9Mock)
MVG_Liveseq_G9Mock <-
  read.csv("05_Liveseq_with_LiveCell_imaging/MVG_Liveseq_G9Mock.csv",
           row.names = 1)
dim(MVG_Liveseq_G9Mock)
MVG_all_G9Mock <- rbind(MVG_scRNA_G9Mock, MVG_Liveseq_G9Mock)
dim(MVG_all_G9Mock)
MVG_all_G9Mock <-
  MVG_all_G9Mock[!duplicated(MVG_all_G9Mock$ensembl_ID) ,]
dim(MVG_all_G9Mock)
write.csv(MVG_all_G9Mock,
          "05_Liveseq_with_LiveCell_imaging/MVG_all_G9Mock.csv")


###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######

###### test the correclation between mCherry fluorescene and most variable genes

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######



#### subset the recorded cells only

Liveseq_recordCell <-
  subset(Seu.all,
         mCherry.log.intercept > 0 &
           treatment == "not_treated" & Batch == "8_8")

### remove 211 genes in black list, which are derived from the 0 pg input RNA negative control.
gene.blacklist <- read.csv(file.path(root_dir, "data/gene.blacklist.csv"))
data.count <- as.matrix(Liveseq_recordCell@assays$RNA@counts)
data.count <-
  data.count[!rownames(data.count) %in% gene.blacklist$ensembl_gene_id,]

##### remove ribosomal protein genes
gene.ribo <-
  Liveseq_recordCell@assays$RNA@meta.features$ensembl_gene_id[grep("^Rp[s,l]",
                                                                   Liveseq_recordCell@assays$RNA@meta.features$external_gene_name)]
gene.ribo <- as.character(gene.ribo)
data.count <- data.count[!rownames(data.count) %in% gene.ribo,]

### Liveseq_recordCell seurat object
Liveseq_recordCell <-
  CreateSeuratObject(counts = data.count, meta.data = Liveseq_recordCell@meta.data)


## add the meta feature
Liveseq_recordCell@assays$RNA@meta.features <-
  Seu.all@assays$RNA@meta.features[rownames(Liveseq_recordCell),]

## add uniquely.mapped.rate and intron.mapped.rate in meta data
Liveseq_recordCell@meta.data$uniquely.mapped.rate <-
  Liveseq_recordCell@meta.data$uniquely.mapped / Liveseq_recordCell@meta.data$input.reads
Liveseq_recordCell@meta.data$intron.mapped.rate <-
  Liveseq_recordCell@meta.data$nCount_RNA / Liveseq_recordCell@meta.data$uniquely.mapped

## QC
VlnPlot(
  Liveseq_recordCell,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt",
    "percent.rRNA",
    "percent.protein",
    "uniquely.mapped.rate",
    "intron.mapped.rate"
  ),
  group.by = "celltype_treatment" ,
  ncol = 3
)

qplot(Liveseq_recordCell$nFeature_RNA)

## normalization
Liveseq_recordCell <-
  NormalizeData(Liveseq_recordCell, verbose = FALSE)
## scale data
Liveseq_recordCell <-
  FindVariableFeatures(
    Liveseq_recordCell,
    selection.method = "vst",
    nfeatures = 500,
    verbose = FALSE
  )


# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase


Liveseq_recordCell <-
  CellCycleScoring(
    Liveseq_recordCell,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = TRUE
  )

saveRDS(Liveseq_recordCell,
        file.path(root_dir, "05_Liveseq_with_LiveCell_imaging/Liveseq_recordCell.rds"))


## cell cycle and slope
p1 <-
  FeatureScatter(Liveseq_recordCell, feature1 = "S.Score", feature2 = "mCherry.log.slope")
p2 <-
  FeatureScatter(Liveseq_recordCell, feature1 = "G2M.Score", feature2 = "mCherry.log.slope")
p3 <-
  FeatureScatter(Liveseq_recordCell, feature1 = "S.Score", feature2 = "mCherry.log.intercept")
p4 <-
  FeatureScatter(Liveseq_recordCell, feature1 = "G2M.Score", feature2 = "mCherry.log.intercept")

plot_grid(p1, p2, p3, p4)

ggsave(
  "05_Liveseq_with_LiveCell_imaging/lm.cellcycle.slope_intercept.3_7.5h.pdf",
  width = 5.5,
  height = 4
)

## the top MVG
top30 <- head(VariableFeatures(Liveseq_recordCell), 30)
plot1 <- VariableFeaturePlot(Liveseq_recordCell)
plot2 <-
  LabelPoints(
    plot = plot1,
    points = top30,
    labels = Liveseq_recordCell@assays$RNA@meta.features[top30, "external_gene_name"] ,
    repel = TRUE
  )
plot2


## transform the data
df <-
  as.data.frame(t(as.matrix(Liveseq_recordCell@assays$RNA@data)))

rownames(df) == colnames(Liveseq_recordCell)

df <-
  cbind(df, Liveseq_recordCell@meta.data[, c("mCherry.log.intercept", "mCherry.log.slope", "mCherry.AUC")])


## linear regression model function
lmMod <- function(exp.trans, targetGene) {
  name.ori <- names(exp.trans)
  exp.trans.renamed <- exp.trans
  names(exp.trans.renamed) <-
    paste("Gene", 1:length(name.ori), sep = "_")
  index.target <- which(name.ori == targetGene)
  
  cof.list <- list()
  strt <- Sys.time()
  ### linear regressoin, take ~ 13 min for this loop
  for (i in 1:length(name.ori)) {
    # sink("log.foreach.txt", append=TRUE)
    X <-
      paste("Gene_",
            index.target,
            " ~ ",
            colnames(exp.trans.renamed)[i],
            sep = "")
    cof <- lm(X , data = exp.trans.renamed)
    cof.list[[i]] <- cof
  }
  print(Sys.time() - strt)
  
  # function to calculate the P value
  lmp <- function (modelobject) {
    if (class(modelobject) != "lm")
      stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    if (is.null(f) == TRUE) {
      warning("this column does not contain fstatistic value")
      return(0)
    }
    p <-
      pf(f[1], f[2], f[3], lower.tail = F)   ## get the p from the fstatistic
    #   p = summary(modelobject)$coefficients[,4]["x"]  ## better use this one if you have varience = 0, in that case this funciont return NA while the function above return 0.
    attributes(p) <- NULL
    return(p)
  }
  
  ## get the p value
  lm.p <- sapply(cof.list, lmp)
  # barplot(-log10(sort(lm.p+(1e-300))), main = "p of lm()")  ## add a pseudo p on the top to avoid infinite value for 0
  
  ## get the R2
  lm.r2 <- sapply(cof.list, function(x)
    summary(x)$r.squared)
  # barplot(sort(lm.r2), ylim = c(0,1), main = "R2 of lm()")
  
  ## get the slope
  lm.slope <- sapply(cof.list, function(x) {
    if (nrow(summary(x)$coefficients) < 2)
      return(NA)
    else
      return(summary(x)$coefficients[2, 1])
  })
  
  # creat the data.frame with all the info
  lm.data <- data.frame(
    ensemble_ID = colnames(exp.trans),
    gene_symbol = sapply(colnames(exp.trans), ensembl.to.symbol),
    p_value = lm.p,
    ### make sure the colume have the same length
    r2 = lm.r2,
    ### make sure the colume have the same length
    slope = lm.slope
  )
  
  lm.data$p_value_adjusted <-
    p.adjust(lm.data$p_value, method = "fdr")
  lm.data$expression_sum <- colSums(exp.trans)
  ## save the lm info
  write.csv(
    lm.data,
    paste(
      "05_Liveseq_with_LiveCell_imaging/lm.",
      deparse(substitute(exp.trans)),
      ".",
      targetGene,
      "3_7.5h.csv",
      sep = ""
    )
  )
  
  lm.data <-
    lm.data %>% arrange(desc(r2)) %>% subset(expression_sum > 3)
  p1 <- ggplot(data = lm.data, aes(p_value, r2)) + geom_point()
  p2 <-
    qplot(
      exp.trans[, as.character(lm.data$ensemble_ID[1])],
      exp.trans[, targetGene],
      geom = c("point", "smooth"),
      method = "lm",
      se = FALSE
    ) + xlab(lm.data$gene_symbol[1])
  p3 <-
    qplot(
      exp.trans[, as.character(lm.data$ensemble_ID[2])],
      exp.trans[, targetGene],
      geom = c("point", "smooth"),
      method = "lm",
      se = FALSE
    ) + xlab(lm.data$gene_symbol[2])
  p4 <-
    qplot(
      exp.trans[, as.character(lm.data$ensemble_ID[3])],
      exp.trans[, targetGene],
      geom = c("point", "smooth"),
      method = "lm",
      se = FALSE
    ) + xlab(lm.data$gene_symbol[3])
  p5 <-
    qplot(
      exp.trans[, as.character(lm.data$ensemble_ID[4])],
      exp.trans[, targetGene],
      geom = c("point", "smooth"),
      method = "lm",
      se = FALSE
    ) + xlab(lm.data$gene_symbol[4])
  p6 <-
    qplot(
      exp.trans[, as.character(lm.data$ensemble_ID[5])],
      exp.trans[, targetGene],
      geom = c("point", "smooth"),
      method = "lm",
      se = FALSE
    ) + xlab(lm.data$gene_symbol[5])
  p7 <-
    qplot(
      exp.trans[, as.character(lm.data$ensemble_ID[6])],
      exp.trans[, targetGene],
      geom = c("point", "smooth"),
      method = "lm",
      se = FALSE
    ) + xlab(lm.data$gene_symbol[6])
  p8 <-
    qplot(
      exp.trans[, as.character(lm.data$ensemble_ID[7])],
      exp.trans[, targetGene],
      geom = c("point", "smooth"),
      method = "lm",
      se = FALSE
    ) + xlab(lm.data$gene_symbol[7])
  p9 <-
    qplot(
      exp.trans[, as.character(lm.data$ensemble_ID[8])],
      exp.trans[, targetGene],
      geom = c("point", "smooth"),
      method = "lm",
      se = FALSE
    ) + xlab(lm.data$gene_symbol[8])
  p10 <-
    qplot(
      exp.trans[, as.character(lm.data$ensemble_ID[9])],
      exp.trans[, targetGene],
      geom = c("point", "smooth"),
      method = "lm",
      se = FALSE
    ) + xlab(lm.data$gene_symbol[9])
  p11 <-
    qplot(
      exp.trans[, as.character(lm.data$ensemble_ID[10])],
      exp.trans[, targetGene],
      geom = c("point", "smooth"),
      method = "lm",
      se = FALSE
    ) + xlab(lm.data$gene_symbol[10])
  p12 <-
    qplot(
      exp.trans[, as.character(lm.data$ensemble_ID[11])],
      exp.trans[, targetGene],
      geom = c("point", "smooth"),
      method = "lm",
      se = FALSE
    ) + xlab(lm.data$gene_symbol[11])
  
  plot_grid(p1,
            p2,
            p3,
            p4,
            p5,
            p6,
            p7,
            p8,
            p9,
            p10,
            p11,
            p12,
            align = "hv",
            axis = "tblr")
  ggsave(
    paste(
      "05_Liveseq_with_LiveCell_imaging/plot_",
      deparse(substitute(exp.trans)),
      "." ,
      targetGene,
      "_lm.3_7.5h.pdf",
      sep = ""
    ),
    width = 12,
    height = 10
  )
  
}

genes_of_interest <-
  c("mCherry.log.intercept", "mCherry.log.slope", "mCherry.AUC")

##############      linear regression on all the genes    ################

library(tibble)

symbol_to_gene <- Liveseq_recordCell@assays$RNA@meta.features %>% rownames_to_column("gene") %>% select(mgi_symbol, gene) %>% deframe()
gene_to_symbol <- Liveseq_recordCell@assays$RNA@meta.features %>% rownames_to_column("gene") %>% select(gene, mgi_symbol) %>% deframe()

X <- as.matrix(Matrix::t(as.matrix(Liveseq_recordCell@assays$RNA@data)[MVG_all_G9Mock$ensembl_ID, ]))
X <- X[, apply(X, 2, sd) > 0.5]
ncol(X)

X <- cbind(X, as.matrix(Liveseq_recordCell@meta.data[, c("S.Score", "G2M.Score")]))

# determine the output (the output_id variable will also be used downstream for output files, etc)
output_id <- "slope"
# output_id <- "intercept"

if (output_id == "slope") {
  y <- Liveseq_recordCell@meta.data$mCherry.log.slope
} else if (output_id == "intercept") {
  y <- Liveseq_recordCell@meta.data$mCherry.log.intercept
} else {
  stop()
}

model <- function(x, y) {
  lm <- lm(y ~ x)
  tibble(
    intercept = lm$coefficients["x"],
    pval = summary(lm)$coefficients[,4]["x"],
    coef = lm$coefficients["x"],
    r2 = summary(lm)$r.squared
  )
}
scores <- purrr::map_dfr(colnames(X), function(gene) {
  x <- X[, gene]
  model(x, y) %>% mutate(gene = gene)
})
scores <- scores %>% mutate(symbol = gene_to_symbol[gene])

scores$fdr <- p.adjust(scores$pval, method = "fdr")

scores %>% arrange(fdr)


# Bootstrapping p-value ---------------------------------------------------
library(furrr)

# calculate scores for samples
n_samples <- 500

# this will take a couple of minutes depending on the number of samples
# we use furrr::future_map here to run things in parallel (using the 5 workers we 'planned' earlier)
resamples_file <- file.path(root_dir, paste0("05_Liveseq_with_LiveCell_imaging/resamples_" , output_id, ".rds"))

if (!file.exists(resamples_file)){
  resamples <- furrr::future_map_dfr(1:n_samples, function(sample_ix) {
    library(Matrix)
    set.seed(sample_ix)
    X_resampled <- X[sample(nrow(X), replace = FALSE),]
    map_dfr(colnames(X_resampled), function(gene) {
      x <- X_resampled[, gene]
      model(x, y) %>% mutate(gene = gene)
    }) %>% mutate(sample_ix = sample_ix)
  }, .options = furrr_options(seed = TRUE))
  saveRDS(resamples, resamples_file)
}
resamples <- readRDS(resamples_file)

resamples <- resamples %>% filter(gene %in% colnames(X))
dim(resamples)

# visualize the distribution for all genes with one gene on top
resamples %>% ggplot(aes(coef, y=..count../sum(..count..))) + # plot normalized counts
  geom_histogram() + # plot for all genes
  geom_histogram(fill = "red", alpha = 0.5, data = function(data){filter(data, gene == symbol_to_gene["Nfkbia"])}) # plot for one gene

# create an empirical distribution for each gene
# essentially just a vector containing all r2 :-)
empirical <- resamples %>% group_by(gene) %>% summarize(coef = list(coef))

# calculate the p-value by checking how many samples have a higher r2 than the actual observed r2
scores2 <- scores %>% 
  left_join(empirical, "gene", suffix = c("", "_empirical")) %>% 
  mutate(
    pval = pmap_dbl(., function(coef, coef_empirical, ...) {
      mean(abs(coef_empirical) >= abs(coef)) *0.99 + 0.01
    })
  ) %>% 
  # filter(abs(coef) > 0) %>% 
  # filter(gene %in% rownames(seu@assays$RNA@meta.features)[log2(seu@assays$RNA@meta.features$vst.mean) > 5])
  identity()
scores2


# plot the p-value distribution
# if all null hypotheses are true, this should look like a uniform distribution between 0 and 1
# if some null hypotheses are false, this distribution will be skewed towards zero
ggplot(scores2, aes(pval)) + geom_histogram()

# correct for multiple testing
scores$fdr <- p.adjust(scores$pval, method = "fdr")
scores2$fdr <- p.adjust(scores2$pval, method = "fdr")

scores2 %>% arrange(fdr) %>% filter(fdr > 0)

qplot(scores2$pval, scores2$fdr)
qplot(scores$pval, scores$fdr)

scores$symbol <- gene_to_symbol[scores$gene]

# volcano plot
plot <- ggplot(scores, aes(coef, fdr)) + 
  geom_point() + 
  ggrepel::geom_label_repel(aes(label = symbol), data = function(data) {data %>% arrange(fdr, -abs(coef)) %>% head(5)}, nudge_y = 0.05) +
  scale_y_continuous(limits = c(0, 1)) +
  ggtitle(output_id) +
  theme_minimal()
plot

# save scores
final_scores <- left_join(
  scores %>% select(pval, fdr, r2, coef, gene, symbol),
  scores2 %>% select(pval, fdr, gene, symbol),
  by = c("gene", "symbol"),
  suffix = c("_lm", "_bootstrap")
) %>% arrange(fdr_lm)

file <- file.path(root_dir, paste0("05_Liveseq_with_LiveCell_imaging/", output_id, "_scores.csv"))
readr::write_excel_csv2(final_scores, file)
openxlsx::write.xlsx(final_scores, file)

wb <- openxlsx::createWorkbook()
openxlsx::saveWorkbook(wb, file.path(root_dir, paste0("05_Liveseq_with_LiveCell_imaging/scores.xlsx")))
wb <- openxlsx::loadWorkbook(file.path(root_dir, paste0("05_Liveseq_with_LiveCell_imaging/scores.xlsx")))
openxlsx::removeWorksheet(wb, output_id)
openxlsx::addWorksheet(wb, output_id)
openxlsx::writeData(wb, output_id, final_scores)

openxlsx::saveWorkbook(wb, file.path(root_dir, paste0("05_Liveseq_with_LiveCell_imaging/scores.xlsx")), overwrite = TRUE)
