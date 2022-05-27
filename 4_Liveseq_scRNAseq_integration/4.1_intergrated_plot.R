library(Seurat)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(scater)
library(reshape2)
library(Matrix)
library(cowplot)
source('~/NAS2/wchen/data_analysis/Resource/R_Funcitons/My_Rfunctions.R', local=TRUE)
# Set working directory
setwd("~/SVFASRAW/wchen/data_analysis/Live_seq/final_analysis_V3/Code_github/")

## read Twodata.RNA object
Twodata.RNA <- readRDS("4_Liveseq_scRNAseq_integration/Intergrated_data.RNA.rds")


#### define color code
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = gg_color_hue(5)
cols   ## "#F8766D" "#00BA38" "#619CFF"
color.batch <- cols

color.celltype <- c("lightgoldenrod4","#dfc27d" , "steelblue1"  , "#80cdc1","#018571") # ,darkslateblue,  tan3
color.samplingtype  <- c("#1c88cc", "#cc601c")

###################### correlation test of scRNA and Liveseq data  #####

raw.data.sub <-as.matrix(Twodata.RNA@assays$RNA@data)

cells.scRNA <- base::subset(Twodata.RNA@meta.data, sampling_type == "scRNA" ) [,"sample_ID"]
cells.scRNA.ASPC <- base::subset(Twodata.RNA@meta.data, sampling_type == "scRNA" &  Cell_type == c("ASPC") &  treatment == "not_treated") [,"sample_ID"]
cells.scRNA.ASPCtreated <- base::subset(Twodata.RNA@meta.data, sampling_type == "scRNA" &  Cell_type == c("ASPC") &  treatment == "DMIR_treated") [,"sample_ID"]

cells.scRNA.IBA <- base::subset(Twodata.RNA@meta.data, sampling_type == "scRNA" &  Cell_type == c("IBA") ) [,"sample_ID"]
cells.scRNA.RawLpsN <- base::subset(Twodata.RNA@meta.data, sampling_type == "scRNA" &  Cell_type %in% c("Raw264.7", "Raw264.7_G9") &  treatment == "not_treated") [,"sample_ID"]
cells.scRNA.RawLpsP <- base::subset(Twodata.RNA@meta.data, sampling_type == "scRNA" &  Cell_type %in% c("Raw264.7", "Raw264.7_G9") &  treatment == "LPS_treated") [,"sample_ID"]

cells.Liveseq <- base::subset(Twodata.RNA@meta.data, sampling_type == "Live_seq" ) [,"sample_ID"]
cells.Liveseq.ASPC <- base::subset(Twodata.RNA@meta.data, sampling_type == "Live_seq" &  Cell_type == c("ASPC") &  treatment == "not_treated") [,"sample_ID"]
cells.Liveseq.ASPCtreated <- base::subset(Twodata.RNA@meta.data, sampling_type == "Live_seq" &  Cell_type == c("ASPC") &  treatment == "DMIR_treated") [,"sample_ID"]

cells.Liveseq.IBA <- base::subset(Twodata.RNA@meta.data, sampling_type == "Live_seq" &  Cell_type == c("IBA") ) [,"sample_ID"]
cells.Liveseq.RawLpsN <- base::subset(Twodata.RNA@meta.data, sampling_type == "Live_seq" &  Cell_type %in% c("Raw264.7", "Raw264.7_G9") &  treatment == "not_treated") [,"sample_ID"]
cells.Liveseq.RawLpsP <- base::subset(Twodata.RNA@meta.data, sampling_type == "Live_seq" &  Cell_type %in% c("Raw264.7", "Raw264.7_G9") &  treatment == "LPS_treated") [,"sample_ID"]


countsum <-  data.frame( cells.scRNA=  rowMeans(raw.data.sub[, cells.scRNA ] ),
                         cells.scRNA.ASPC = rowMeans(raw.data.sub[, cells.scRNA.ASPC ] ) ,
                         cells.scRNA.ASPCtreated = rowMeans(raw.data.sub[, cells.scRNA.ASPCtreated ] ) ,
                         cells.scRNA.IBA=  rowMeans(raw.data.sub[, cells.scRNA.IBA ] ) ,
                         cells.scRNA.RawLpsN =  rowMeans(raw.data.sub[, cells.scRNA.RawLpsN ] ) ,
                         cells.scRNA.RawLpsP =  rowMeans(raw.data.sub[, cells.scRNA.RawLpsP ] ) ,
                         cells.Liveseq  = rowMeans( raw.data.sub[,cells.Liveseq] ) ,
                         cells.Liveseq.ASPC = rowMeans( raw.data.sub[,cells.Liveseq.ASPC ] ) ,
                         cells.Liveseq.ASPCtreated = rowMeans( raw.data.sub[,cells.Liveseq.ASPCtreated ] ) ,
                         cells.Liveseq.IBA  = rowMeans( raw.data.sub[,cells.Liveseq.IBA ] ) ,
                         cells.Liveseq.RawLpsN =  rowMeans(raw.data.sub[, cells.Liveseq.RawLpsN ] ),
                         cells.Liveseq.RawLpsP =  rowMeans(raw.data.sub[, cells.Liveseq.RawLpsP ] )

)

write.csv(countsum, "4_Liveseq_scRNAseq_integration/mimicBulk_gene_expression_normalized.csv")



## ggpairs plot the correlation  
## http://www.sthda.com/english/wiki/ggally-r-package-extension-to-ggplot2-for-correlation-matrix-and-survival-plots-r-software-and-data-visualization

library(GGally)
library(viridis)
GGscatterPlot <- function(data, mapping, ..., 
                          method = "pearson") {
  
  #Get correlation coefficient
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  
  cor <- cor(x, y, method = method)
  #Assemble data frame
  df <- data.frame(x = x, y = y)
  # PCA
  nonNull <- x!=0 & y!=0
  dfpc <- prcomp(~x+y, df[nonNull,])
  df$cols <- predict(dfpc, df)[,1]
  # Define the direction of color range based on PC1 orientation:
  dfsum <- x+y
  colDirection <- ifelse(dfsum[which.max(df$cols)] < 
                           dfsum[which.min(df$cols)],
                         1,
                         -1)
  #Get 2D density for alpha
  dens2D <- MASS::kde2d(df$x, df$y)
  df$density <- fields::interp.surface(dens2D , 
                                       df[,c("x", "y")])
  
  if (any(df$density==0)) {
    mini2D = min(df$density[df$density!=0]) #smallest non zero value
    df$density[df$density==0] <- mini2D
  }
  #Prepare plot
  pp <- ggplot(df, aes(x=x, y=y, color = cols, alpha = 1/density)) +
    ggplot2::geom_point(shape=16, show.legend = FALSE) +
    ggplot2::scale_color_viridis_c(direction = colDirection) +
    #                scale_color_gradient(low = "#0091ff", high = "#f0650e") +
    ggplot2::scale_alpha(range = c(.05, .6)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, col="darkred") +
    ggplot2::geom_label(
      data = data.frame(
        xlabel = min(x, na.rm = TRUE),
        ylabel = max(y, na.rm = TRUE),
        lab = round(cor, digits = 3)),
      mapping = ggplot2::aes(x = xlabel, 
                             y = ylabel, 
                             label = lab),
      hjust = 0, vjust = 1,
      size = 3, fontface = "bold",
      inherit.aes = FALSE # do not inherit anything from the ...
    ) +
    theme_minimal()
  
  return(pp)
}


countsum.plot <- countsum[,c(2,3,4,5,6,8,9,10,11,12)]
### filter lowly expressed genes
countsum.plot <- countsum.plot[ rowMeans(countsum.plot) >0.01 , ]


GGally::ggpairs(countsum.plot,
                lower = list(continuous = GGscatterPlot),
                upper = list(continuous = wrap("cor", method= "spearman")))
# ggsave("4_Liveseq_scRNAseq_integration/Livesesq_scRNAseq_cor.pdf", width = 10, height = 10, useDingbats=F)


## evaluate_cluster_RNA
p <- ggplot(Twodata.RNA@meta.data,aes(x=RNA.ident,fill=celltype_treatment))+
  geom_bar(stat="count")+
  geom_text(aes(label=stat(count)),stat="count",position=position_stack(0.5)) +
  xlab("cluster") + scale_fill_manual(values = color.celltype)
p
# ggsave("4_Liveseq_scRNAseq_integration/evaluate_cluster_RNA.pdf", width = 3, height = 2)

## dimplot
p1 <- DimPlot(Twodata.RNA, reduction = "umap", pt.size = 0.3, label = TRUE, repel = TRUE, cols = color.celltype ) + NoLegend()

p2 <- DimPlot(Twodata.RNA, reduction = "umap", group.by = "celltype_treatment", pt.size = 0.3, label = TRUE, repel = TRUE, cols = color.celltype ) + 
  NoLegend()

p3 <- DimPlot(Twodata.RNA, reduction = "umap", group.by = "sampling_type", pt.size = 0.3, cols = c("#1c88cc", "#cc601c")) 

p4 <- FeaturePlot(Twodata.RNA, reduction  = "umap", features = "nFeature_RNA", pt.size = 0.3)

p5 <- FeaturePlot(Twodata.RNA, reduction  = "umap", features = "nCount_RNA", pt.size = 0.3)


plot_grid(p1, p3, p2,p4,p5, ncol = 2, align = "hv", axis = "lrtb")    ### axis indicates the alignment of l: left, r: right, t: top, b: bottom
ggsave("4_Liveseq_scRNAseq_integration/cluster_umap_RNA.pdf", width = 5.4, height = 6, useDingbats=F )



p1 <- DimPlot(Twodata.RNA, reduction = "tsne", pt.size = 0.05, label = TRUE, repel = TRUE, cols = color.celltype ) + NoLegend()

p2 <- DimPlot(Twodata.RNA, reduction = "tsne", group.by = "celltype_treatment", pt.size = 0.05, label = TRUE, repel = TRUE, cols = color.celltype ) + 
  NoLegend()

p3 <- DimPlot(Twodata.RNA, reduction = "tsne", group.by = "sampling_type", pt.size = 0.05, cols = c("#1c88cc", "#cc601c")) 

p4 <- FeaturePlot(Twodata.RNA, reduction  = "tsne", features = "nFeature_RNA", pt.size = 0.05)

p5 <- FeaturePlot(Twodata.RNA, reduction  = "tsne", features = "nCount_RNA", pt.size = 0.05)

plot_grid(p1, p3, p2,p4,p5, ncol = 2, align = "hv", axis = "lrtb")    ### axis indicates the alignment of l: left, r: right, t: top, b: bottom
# ggsave("4_Liveseq_scRNAseq_integration/cluster_tsne_RNA.pdf", width = 5.4, height = 6, useDingbats=F )

## plot cell cycle
s.genes <- readRDS("s.genes.mouse.rds")    
g2m.genes <- readRDS("g2m.genes.mouse.rds")    

Twodata.RNA <- CellCycleScoring(Twodata.RNA, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

p6 <- DimPlot(Twodata.RNA, reduction = "tsne", group.by = "Phase", pt.size = 0.05) 

p7 <- FeaturePlot(Twodata.RNA, reduction  = "tsne", features = "S.Score", pt.size = 0.05)

p8 <- FeaturePlot(Twodata.RNA, reduction  = "tsne", features = "G2M.Score", pt.size = 0.05)

p9 <- DimPlot(Twodata.RNA, reduction = "tsne", group.by = "Batch", pt.size = 0.05) 

plot_grid(p6,p7,p8,p9, ncol = 2, align = "hv", axis = "lrtb")    ### axis indicates the alignment of l: left, r: right, t: top, b: bottom

# ggsave("4_Liveseq_scRNAseq_integration/cluster_tsne_RNA2.pdf", width = 5.4, height = 4, useDingbats=F )



p1 <- DimPlot(Twodata.RNA, reduction = "pca", pt.size = 0.3, label = TRUE, repel = TRUE, cols = color.celltype ) + NoLegend()

p2 <- DimPlot(Twodata.RNA, reduction = "pca", group.by = "celltype_treatment", pt.size = 0.3, label = TRUE, repel = TRUE, cols = color.celltype ) + 
  NoLegend()

p3 <- DimPlot(Twodata.RNA, reduction = "pca", group.by = "sampling_type", pt.size = 0.3, cols = c("#1c88cc", "#cc601c")) 

p4 <- FeaturePlot(Twodata.RNA, reduction  = "pca", features = "nFeature_RNA", pt.size = 0.3)

p5 <- FeaturePlot(Twodata.RNA, reduction  = "pca", features = "nCount_RNA", pt.size = 0.3)

plot_grid(p1, p3, p2,p4,p5, ncol = 2, align = "hv", axis = "lrtb")    ### axis indicates the alignment of l: left, r: right, t: top, b: bottom
# ggsave("4_Liveseq_scRNAseq_integration/cluster_pca_RNA.pdf", width = 5.4, height = 6, useDingbats=F )


## mCherry and Tnf scatter plot
p1 <- FeatureScatter(subset(Twodata.RNA, Cell_type == "Raw264.7_G9"), feature1 = "mCherry", feature2 = symbol.to.ensembl("Tnf") , group.by = "treatment", pt.size = 0.3,cols = color.celltype[3:2])  + ylab("Tnf")
p2 <- FeatureScatter(subset(Twodata.RNA, Cell_type == "Raw264.7_G9" & sampling_type == "Live_seq"), feature1 = "mCherry", feature2 = symbol.to.ensembl("Tnf") , group.by = "treatment", pt.size = 0.3, cols = color.celltype[3:2])  + ylab("Tnf")

plot_grid(p1,p2, align = "vh", axis = "lrtb")
# ggsave("4_Liveseq_scRNAseq_integration/Tnf_mcherry_cor_RNA.pdf", width=5.4, height = 2, useDingbats=F)





#########  switch to integrated assay   #######



Twodata.integrated <- readRDS("4_Liveseq_scRNAseq_integration/Intergrated_data.integrated.rds")
DefaultAssay(Twodata.integrated) <- "RNA"

s.genes <- readRDS("s.genes.mouse.rds")    
g2m.genes <- readRDS("g2m.genes.mouse.rds")    

Twodata.integrated <- CellCycleScoring(Twodata.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# Idents(Twodata.integrated) <- Twodata.integrated$integrated_snn_res.0.2 


# switch to integrated assay. The variable features of this assay are automatically set during
DefaultAssay(Twodata.integrated) <- "integrated"


## remove the pseudogene from DE analysis
geneName <- read.table(file = "mouseGeneTable87_mCherry_EGFP.txt", sep = "\t", header = T, row.names = 1 )
gene.pseudogene <- geneName[endsWith(as.character(geneName$gene_biotype),  "pseudogene"),]
gene.pseudoRemoved <- subset(Twodata.integrated@assays$RNA@meta.features, !(ensembl_gene_id %in% gene.pseudogene$ensembl_gene_id))




### evaluate the clustering accuracy 
Twodata.integrated@meta.data$interated.ident <- Twodata.integrated$integrated_snn_res.0.2
Idents(Twodata.integrated) <- Twodata.integrated$integrated_snn_res.0.2
p1 <- ggplot(Twodata.integrated@meta.data,aes(x=interated.ident,fill=celltype_treatment))+
  geom_bar()+
  geom_text(aes(label=stat(count)),stat="count",position=position_stack(0.5)) +
  xlab("cluster") + scale_fill_manual(values = color.celltype)
p1
# ggsave("4_Liveseq_scRNAseq_integration/evaluate_cluster_intergrated.pdf", width = 3, height = 2)




p1 <- DimPlot(Twodata.integrated, reduction = "tsne", pt.size = 0.15, label = TRUE, repel = TRUE, cols = color.celltype) + NoLegend()

p2 <- DimPlot(Twodata.integrated, reduction = "tsne", group.by = "celltype_treatment", pt.size = 0.15, label = TRUE, repel = TRUE, cols = color.celltype) + 
  NoLegend()

p3 <- DimPlot(Twodata.integrated, reduction = "tsne", group.by = "sampling_type", pt.size = 0.15, cols = c("#1c88cc", "#cc601c")) 

p4 <- FeaturePlot(Twodata.integrated, reduction  = "tsne", features = "nFeature_RNA", pt.size = 0.15)

p5 <- FeaturePlot(Twodata.integrated, reduction  = "tsne", features = "nCount_RNA", pt.size = 0.15)

p6 <- DimPlot(Twodata.integrated, reduction = "tsne", group.by = "Phase", pt.size = 0.15) 

p7 <- FeaturePlot(Twodata.integrated, reduction  = "tsne", features = "S.Score", pt.size = 0.15)

p8 <- FeaturePlot(Twodata.integrated, reduction  = "tsne", features = "G2M.Score", pt.size = 0.15)

p9 <- DimPlot(Twodata.integrated, reduction = "tsne", group.by = "Batch", pt.size = 0.15) 


plot_grid(p1, p3, p2,p4,p5,p6 , ncol = 2, align = "hv", axis = "lrtb")    ### axis indicates the alignment of l: left, r: right, t: top, b: bottom
# ggsave("4_Liveseq_scRNAseq_integration/cluster_tsne_intergrated.pdf", width = 6.4, height = 6.6, useDingbats=F )

plot_grid(p7, p8, p9, ncol = 2, align = "hv", axis = "lrtb")    ### axis indicates the alignment of l: left, r: right, t: top, b: bottom
# ggsave("4_Liveseq_scRNAseq_integration/cluster_tsne_intergrated2.pdf", width = 6.4, height = 4.4, useDingbats=F )




## culstomized tsne plot
tSNE.embedding <- Twodata.integrated@reductions$tsne@cell.embeddings

all(rownames(tSNE.embedding) == rownames(Twodata.integrated@meta.data))
tSNE.embedding  <- cbind(tSNE.embedding, Twodata.integrated@meta.data) 

p1 <- ggplot(tSNE.embedding %>% filter(sampling_type == "Live_seq"), aes(tSNE_1, tSNE_2, color=celltype_treatment)) +
  geom_point(size = 0.5) +
  scale_color_manual(values=color.celltype) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

p2 <- ggplot(tSNE.embedding %>% filter(sampling_type == "scRNA"), aes(tSNE_1, tSNE_2, color=celltype_treatment)) +
  geom_point(size = 0.3) +
  scale_color_manual(values=color.celltype) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot_grid(p1,p2)

ggsave("4_Liveseq_scRNAseq_integration/tSNE_intergrated_different_visilizaiton.pdf", width = 7.5, height = 2, useDingbats=F  )




p1 <- DimPlot(Twodata.integrated, reduction = "pca", pt.size = 0.15, label = TRUE, repel = TRUE, cols = color.celltype) + NoLegend()

p2 <- DimPlot(Twodata.integrated, reduction = "pca", group.by = "celltype_treatment", pt.size = 0.15, label = TRUE, repel = TRUE, cols = color.celltype) + 
  NoLegend()

p3 <- DimPlot(Twodata.integrated, reduction = "pca", group.by = "sampling_type", pt.size = 0.15, cols = c("#1c88cc", "#cc601c")) 

p4 <- FeaturePlot(Twodata.integrated, reduction  = "pca", features = "nFeature_RNA", pt.size = 0.15)

p5 <- FeaturePlot(Twodata.integrated, reduction  = "pca", features = "nCount_RNA", pt.size = 0.15)

plot_grid(p1, p3, p2,p4,p5, ncol = 2, align = "hv", axis = "lrtb")    ### axis indicates the alignment of l: left, r: right, t: top, b: bottom
ggsave("4_Liveseq_scRNAseq_integration/cluster_pca_intergrated.pdf", width = 5.4, height = 6, useDingbats=F )





### plot the DE genes
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

genesToPlot <- c("Dpep1", "Thy1", "Gpx3","Sfrp2" ,"Mycn","Tnf" )

plist <- lapply(genesToPlot, function(x) {
    p <- FeaturePlot(Twodata.integrated, reduction  = "tsne", features = symbol.to.ensembl(x), pt.size = 0.15) + ggtitle(x)
    return(p)
  }
)

plot_grid(plotlist=plist, align = "hv", axis = "lrtb", ncol = 2)
# ggsave("4_Liveseq_scRNAseq_integration/gene_expression_tsne.pdf", width = 4.7, height = 6, useDingbats=FALSE)



#### label the double extraction samples

Twodata.integrated$double_extraction <- factor(Twodata.integrated$double_extraction)

# plot list of the double extraction
plist <- list()
for (i in levels(Twodata.integrated$double_extraction)) {
  cells.selected <- names( Twodata.integrated$double_extraction[ Twodata.integrated$double_extraction == i] ) ## select double extraction cells
  treatment <- Twodata.integrated$treatment[ Twodata.integrated$double_extraction == i]           ## get the LPS treatment info of double extraction cells
  double_extraction_order <- Twodata.integrated$double_extraction_order[ Twodata.integrated$double_extraction == i]  ## get the sampling order (1st or 2nd sampling)
  # ilabel <- paste(cells.selected, treatment, double_extraction_order,    collapse = " ")
  ilabel <- paste(Twodata.integrated@meta.data[cells.selected, "original_sample_name"], treatment, double_extraction_order,    collapse = " ")
  ## plot double extraction sample
  # icolor <- rep("grey", ncol(Liveseq_sub))
  # icolor[Liveseq_sub$double_extraction == i & Liveseq_sub$double_extraction_order == 1 ] <- "#619CFF"
  # icolor[Liveseq_sub$double_extraction == i & Liveseq_sub$double_extraction_order == 2 ] <- "#00BA38"
  
  
  p<- DimPlot(Twodata.integrated, reduction = "tsne", group.by = "double_extraction", cells.highlight = cells.selected,  sizes.highlight = 1.5, pt.size = 0.3) + 
    ggtitle(ilabel ) 
  plist[[i]] <- p
}

# add the celltype_treatment plot in the first list
p1 <- DimPlot(Twodata.integrated, reduction = "tsne", cols = color.celltype,label = TRUE,label.size = 2 ,repel = T, pt.size = 0.3) + NoLegend()
plist[[1]] <- p1

plot_grid(plotlist = plist[1:16], align = "hv", axis = "lrtb")
# ggsave("4_Liveseq_scRNAseq_integration/double_extraction1.tsne.pdf", width = 12, height = 8,  useDingbats=FALSE)

plot_grid(plotlist = plist[17:32], align = "hv", axis = "lrtb")
# ggsave("4_Liveseq_scRNAseq_integration/double_extraction2.tsne.pdf", width = 12, height = 8,  useDingbats=FALSE)

plot_grid(plotlist = plist[33:48], align = "hv", axis = "lrtb")
# ggsave("4_Liveseq_scRNAseq_integration/double_extraction3.tsne.pdf", width = 12, height = 8,  useDingbats=FALSE)

plot_grid(plotlist = plist[49:length(plist)], align = "hv", axis = "lrtb", cols = 4)
# ggsave("4_Liveseq_scRNAseq_integration/double_extraction4.tsne.pdf", width = 12, height = 4,  useDingbats=FALSE)



### downsamping data for integration ##$$

downsample.list <- readRDS("4_Liveseq_scRNAseq_integration/downsample.list.rds")
## evulate the clustering in downsampled data
plist <- lapply( downsample.list, function(x) {
  x$idents <- Idents(x)
  p <- ggplot(x@meta.data,aes(x=idents,fill=celltype_treatment))+
    geom_bar()+
    geom_text(aes(label=stat(count)),stat="count",position=position_stack(0.5), size=3) +
    xlab("Cluster") +
    scale_fill_manual(values =  color.celltype )+
    ggtitle(mean(x$nCount_RNA))
  return(p)
})

plot_grid(plotlist = plist)

## adjust the level accordingly
downsample.list <- lapply( downsample.list, function(x) {
  levels(Idents(x)) <- c(4, 5, 3, 1,2)
  Idents(x) <- factor(Idents(x), levels = c(1,2,3,4,5) )
  x$celltype_treatment <- factor(x$celltype_treatment, levels = c("ASPC_not_treated", "ASPC_DMIR_treated","IBA_not_treated","Raw264.7_not_treated", "Raw264.7_LPS_treated") )
  return(x)
})

plist <- lapply( downsample.list[1:7], function(x) {
  x$idents <- Idents(x)
  p <- ggplot(x@meta.data,aes(x=idents,fill=celltype_treatment))+
    geom_bar()+
    geom_text(aes(label=stat(count)),stat="count",position=position_stack(0.5), size=3) +
    xlab("Cluster") +
    scale_fill_manual(values =  color.celltype )+
    ggtitle(mean(x$nCount_RNA))
  return(p)
})

plot_grid(plotlist = plist)
# ggsave("4_Liveseq_scRNAseq_integration/evaluate_clustering_downsampling.pdf", width = 9, height = 6)


plist <- lapply( downsample.list, function(x) {
 p<- DimPlot(x, reduction = "tsne", group.by = "celltype_treatment",  pt.size = 0.2,cols = color.celltype) 
})
plot_grid(plotlist = plist)
# ggsave("4_Liveseq_scRNAseq_integration/clustering_tsne_downsampling.pdf", width = 15, height = 7.5, useDingbats=FALSE)


plist <- lapply( downsample.list, function(x) {
  p<- DimPlot(x, reduction = "tsne", group.by = "sampling_type", pt.size = 0.2,cols = color.samplingtype) 
})
plot_grid(plotlist = plist)
# ggsave("4_Liveseq_scRNAseq_integration/clustering_tsne_downsampling_samplingtype.pdf", width = 12, height = 7.5,useDingbats=FALSE )









