################################################################
#                                                              #
#                         Plot scRNA-seq                       #
#                                                              #
################################################################

### Questions: https://github.com/DeplanckeLab/Live-seq/issues
### Date: 2022-03-06
### Datasets: scRNA-seq
### Goal: Making some summary plots of the scRNA-seq data

library(rprojroot)
root_dir <- find_root(has_file("Live-seq.RProj"))

library(Seurat)
library(dplyr)
library(grid)
library(gridExtra)
library(scater)
library(reshape2)
library(Matrix)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(viridis)

#######################

scRNA <- readRDS("03_scRNA_seq_only/scRNAseq_only.rds")

mean(scRNA$nFeature_RNA)    ### nFeature mean = 8328.863

### save sample_name for RNA velocity analysis
cells.scRNA.exp7 <- as.character(subset(scRNA, Batch == "scRNA")$sample_name)
cells.scRNA.exp8_8 <- as.character(subset(scRNA, Batch == "8_8")$sample_name)
cells.scRNA.exp9_1 <- as.character(subset(scRNA, Batch == "9_1")$sample_name)

write.csv(cells.scRNA.exp7, "03_scRNA_seq_only/cells.scRNA.exp7.csv")
write.csv(cells.scRNA.exp8_8, "03_scRNA_seq_only/cells.scRNA.exp8_8.csv")
write.csv(cells.scRNA.exp9_1, "03_scRNA_seq_only/cells.scRNA.exp9_1.csv")

#### get the default color of ggplot2. n is the number of colors needed. 
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = gg_color_hue(3)
cols   ## "#F8766D" "#00BA38" "#619CFF"

color.batch <- cols
color.celltype <- c("lightgoldenrod4","#dfc27d" , "steelblue1"  , "#80cdc1","#018571") # ,darkslateblue,  tan3
scRNA$celltype_treatment <- factor(scRNA$celltype_treatment , levels =c("ASPC_not_treated","ASPC_DMIR_treated", "IBA_not_treated","Raw264.7_not_treated" ,"Raw264.7_LPS_treated" ))

## plot the number of samples in each step
df <- data.frame(catalog=c("Total",  "Passing QC"), number=c(573,  ncol(scRNA) ))
df$catalog <- factor(df$catalog , levels = df$catalog )
p <- ggplot(df, aes(catalog, number)) + geom_bar(stat = "identity") + xlab(label = "") + ylim(0,400)
p <- p + theme(
  axis.ticks = element_line(colour = "black", size = 0.2),
  axis.text.x = element_text(angle = 45, colour = "black", hjust = 1),
  axis.text.y = element_text(colour = "black")) +
  stat_summary(geom = 'text', label = df$number, fun.y = max, vjust = -1)
p <- p +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,700))
p
# ggsave("03_scRNA_seq_only/sample_number.pdf", width = 1, height = 2, useDingbats=FALSE)


## plot the basic QC, group by batch or sample type

plist <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rRNA", "percent.protein", "uniquely.mapped.rate", "intron.mapped.rate", "input.reads"), function(i) { 
  # p <- myVlnPlot(Liveseq_sub, x, "sampling_type", cols = c("grey","#00BA38","#619CFF"))
  p <- ggplot(scRNA@meta.data, aes_string(x="sampling_type", y=i)) + 
    geom_violin(trim=FALSE,  fill="gray")+
    stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color=viridis(10)[1]) +
    theme_classic()
  p <- p + geom_jitter(shape=16, position=position_jitter(0.2))

  return(p)})                 

plot_grid(plotlist = plist,  ncol = 8, align = "hv", axis = "tblr")

# ggsave("03_scRNA_seq_only/QC_noGroup.pdf", width = 12, height = 1.6, useDingbats=FALSE )





ylab.metrics <- c("Number of reads", "Fraction of input reads", "Fraction of uniquely mapped reads", "Total count of all genes", "Number of detected genes", "Fraction of nCount (%)", "Fraction of nCount (%)","Fraction of nCount (%)")
names(ylab.metrics) <-  c("input.reads","uniquely.mapped.rate","intron.mapped.rate",  "nCount_RNA","nFeature_RNA", "percent.mt", "percent.rRNA", "percent.protein")
title.metrics <- c("Input reads", "Uniquely mapped reads", "Reads mapped in exon", "nCount", "nGene", "Percent MT","Percent rRNA", "Percent Protein")
names(title.metrics) <-  c("input.reads","uniquely.mapped.rate","intron.mapped.rate",  "nCount_RNA","nFeature_RNA", "percent.mt", "percent.rRNA", "percent.protein")

plist <- lapply(c("input.reads","uniquely.mapped.rate","intron.mapped.rate",  "nCount_RNA","nFeature_RNA", "percent.mt", "percent.rRNA", "percent.protein" ), function(x) {
  p <- VlnPlot(object = scRNA, features = x,group.by = "celltype_treatment", cols = color.celltype )
  p <- p + scale_x_discrete(labels= c("ASPC_Pre", "ASPC_Post", "IBA", "RAW_Mock", "RAW_LPS"))
  p <- p + xlab("")
  p <- p + ylab(ylab.metrics[x]) + ggtitle(title.metrics[x])
  return(p)
})
plot_grid(plotlist = plist, align = "hv", axis = "lrtb")
# ggsave("03_scRNA_seq_only/QC_celltype.pdf", width = 5, height = 6, useDingbats=FALSE )


## plot scatter plot 
p1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "celltype_treatment", pt.size = 0.3, cols = color.celltype) 
p2 <- FeatureScatter(scRNA, feature1 = "cDNA_concentration", feature2 = "nFeature_RNA", group.by = "celltype_treatment", pt.size = 0.3, cols = color.celltype) 

plot_grid(p1,p2, align = "hv", axis = "lrtb")
# ggsave("03_scRNA_seq_only/scatter_cDNA,nCount_nFeature_byCelltype.pdf", width = 6, height = 1.5, useDingbats=FALSE)


## save the name of each cell type and treatment
ASPC.cell <- colnames(subset(scRNA, celltype_treatment == "ASPC_not_treated"))
ASPCtreated.cell <- colnames(subset(scRNA, celltype_treatment == "ASPC_DMIR_treated"))
IBA.cell <- colnames(subset(scRNA, celltype_treatment == "IBA_not_treated"))
RAWLpsN.cell <- colnames(subset(scRNA, celltype_treatment %in% c("Raw264.7_not_treated")))
RAWLpsP.cell <- colnames(subset(scRNA, celltype_treatment %in% c("Raw264.7_LPS_treated")))

## save the name of each cell type and treatment
save(ASPC.cell,ASPCtreated.cell, IBA.cell,RAWLpsN.cell,RAWLpsP.cell, file = "03_scRNA_seq_only/cells.types.RData")
write(c(ASPC.cell, " ", ASPCtreated.cell, " ",  IBA.cell," " ,RAWLpsN.cell," ", RAWLpsP.cell), "03_scRNA_seq_only/cells.types.txt", sep=" ")


### compare the cluster with ground trust.
### two categorical variable stacked barplot with label  

Idents(scRNA) <- scRNA$RNA_snn_res.0.2  ### adjust the res accordingly 
scRNA@meta.data$ident <- scRNA$RNA_snn_res.0.3
p <- ggplot(scRNA@meta.data,aes(x=ident,fill=celltype_treatment))+
  geom_bar()+
  geom_text(aes(label=stat(count)),stat="count",position=position_stack(0.5), size=3) +
  xlab("cluster") + scale_fill_manual(values = color.celltype)
p
ggsave("03_scRNA_seq_only/evaluate_cluster_RNA.res.0.3.pdf", width = 2.6, height = 1.8)


## umap plots
DimPlot(scRNA, reduction = "pca", label = TRUE,label.size = 2 ,repel = T, pt.size = 0.3)
p1 <- DimPlot(scRNA, reduction = "umap", label = TRUE,label.size = 2 ,repel = T, pt.size = 0.3) + NoLegend()

p2 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype_treatment",label = TRUE,label.size = 2 ,repel = T, pt.size = 0.3, cols = color.celltype) + NoLegend()

p3 <- FeaturePlot(scRNA, reduction  = "umap", features = "nFeature_RNA", pt.size = 0.3)

plot_grid(p1, p2, p3, align = "hv", axis = "lrtb", ncol = 3)
# ggsave("03_scRNA_seq_only/cluster_umap_RNA.res.0.3.pdf", width = 7.1, height = 1.8, useDingbats=FALSE)




## tsne plots
p1 <- DimPlot(scRNA, reduction = "tsne", label = TRUE,label.size = 2 ,repel = T, pt.size = 0.3) + NoLegend()

p2 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype_treatment",label = TRUE,label.size = 2 ,repel = T, pt.size = 0.3, cols = color.celltype) + NoLegend()

p3 <- FeaturePlot(scRNA, reduction  = "tsne", features = "nFeature_RNA", pt.size = 0.3)

plot_grid(p1, p2, p3, align = "hv", axis = "lrtb", ncol = 3)
ggsave("03_scRNA_seq_only/cluster_tsne_RNA.res.0.3.pdf", width = 7.1, height = 1.8, useDingbats=FALSE)



## tsne plots with differet RNA_snn_res.
p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "RNA_snn_res.0.2", label = TRUE,label.size = 2 ,repel = T, pt.size = 0.3) + NoLegend()

p2 <- DimPlot(scRNA, reduction = "tsne", group.by = "RNA_snn_res.0.3", label = TRUE,label.size = 2 ,repel = T, pt.size = 0.3) + NoLegend()

p3 <- DimPlot(scRNA, reduction = "tsne", group.by = "RNA_snn_res.0.4", label = TRUE,label.size = 2 ,repel = T, pt.size = 0.3) + NoLegend()

p4 <- DimPlot(scRNA, reduction = "tsne", group.by = "RNA_snn_res.0.5", label = TRUE,label.size = 2 ,repel = T, pt.size = 0.3) + NoLegend()

p5 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype_treatment",label = TRUE,label.size = 2 ,repel = T, pt.size = 0.3, cols = color.celltype) + NoLegend()


plot_grid(p1, p2, p3,p4,p5, align = "hv", axis = "lrtb", ncol = 3)
# ggsave("03_scRNA_seq_only/cluster_tsne_different_res.pdf", width = 7.1, height = 3.8, useDingbats=FALSE)




## pca plots
p1 <- DimPlot(scRNA, reduction = "pca", cols = color.celltype,label = TRUE,label.size = 2 ,repel = T, pt.size = 0.3) + NoLegend()

p2 <- DimPlot(scRNA, reduction = "pca", group.by = "celltype_treatment",label = TRUE,label.size = 2 ,repel = T, pt.size = 0.3, cols = color.celltype) + NoLegend()

p3 <- FeaturePlot(scRNA, reduction  = "pca", features = "nFeature_RNA", pt.size = 0.3)

plot_grid(p1, p2, p3, align = "hv", axis = "lrtb", ncol = 3)
# ggsave("03_scRNA_seq_only/cluster_pca_RNA.res.0.3.pdf", width = 7.1, height = 1.8, useDingbats=FALSE)



####################### plot the batch corrected scRNA data

scRNA.integrated <- readRDS("03_scRNA_seq_only/scRNAseq_batch_corrected.rds")

## reorder the celltype_treatment factor 
scRNA.integrated$celltype_treatment <- factor(scRNA.integrated$celltype_treatment , levels = c("ASPC_not_treated", "ASPC_DMIR_treated", "IBA_not_treated","Raw264.7_not_treated", "Raw264.7_LPS_treated"))

## tsne plots
p1 <- DimPlot(scRNA.integrated, reduction = "tsne", label = TRUE,label.size = 2 ,repel = T, pt.size = 0.15) + NoLegend()

p2 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.6", label = TRUE,label.size = 2 ,repel = T, pt.size = 0.15) + NoLegend()

p3 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.1", label = TRUE,label.size = 2 ,repel = T, pt.size = 0.15) + NoLegend()

p4 <- DimPlot(scRNA.integrated, reduction = "tsne", label = TRUE,label.size = 2 ,repel = T, pt.size = 0.15, cols = color.celltype) + NoLegend()

p5 <- DimPlot(scRNA.integrated, reduction = "tsne", group.by = "celltype_treatment",label = TRUE,label.size = 2 ,repel = T, pt.size = 0.15, cols = color.celltype) + NoLegend()

p6 <- FeaturePlot(scRNA.integrated, reduction  = "tsne", features = "nFeature_RNA", pt.size = 0.15)

p7 <- DimPlot(scRNA.integrated, reduction  = "tsne",  group.by = "Phase",pt.size = 0.15)

p8 <- FeaturePlot(scRNA.integrated, reduction  = "tsne",  features = "G2M.Score",pt.size = 0.15)

p9 <- FeaturePlot(scRNA.integrated, reduction  = "tsne",  features = "S.Score",pt.size = 0.15)

plot_grid(p1, p2, p3,p4,p5,p6,p7,p8,p9, align = "hv", axis = "lrtb", ncol = 3)
# ggsave("03_scRNA_seq_only/cluster_tsne_batchCorrected.pdf", width = 9.6, height = 7.2, useDingbats=FALSE)


### compare the cluster with ground trust.
### two categorical variable stacked barplot with label  

scRNA.integrated$ident <- Idents(scRNA.integrated)

p <- ggplot(scRNA.integrated@meta.data,aes(x=ident,fill=celltype_treatment))+
  geom_bar()+
  geom_text(aes(label=stat(count)),stat="count",position=position_stack(0.5), size=3) +
  xlab("cluster") + scale_fill_manual(values=color.celltype) +
  ylim(c(0,200))
## remove the grid
p <- p +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,200))+
  ylab("Number of cells") +
  xlab("Cluster")
p
# ggsave("03_scRNA_seq_only/evaluate_cluster_batchCorrected.pdf", width = 3, height = 2)

##### plot the percentage
p <- ggplot(scRNA.integrated@meta.data,aes(x=ident,fill=celltype_treatment))+
  geom_bar(position="fill", stat="count")+
  geom_text(aes(label=stat(count)),stat="count",position=position_stack(0.5), size=3) +
  xlab("cluster") + scale_fill_manual(values=color.celltype) +
  ylim(c(0,1))
## remove the grid
p <- p +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))+
  ylab("Number of cells") +
  xlab("Cluster") 
p

# ggsave("03_scRNA_seq_only/evaluate_cluster_batchCorrected_percent.pdf", width = 3, height = 2)






## plot heatmap of top DE of scRNA-seq data

scRNA.markers <- read.csv("03_scRNA_seq_only/DEs.scRNA.csv")
## the batch corretion incorporated DE analyis, refer to scRNA.markers.celltype.edgeR.csv below
# scRNA.markers <- read.csv("../LiveSeq_vsscRNAseq_9.09.21/scRNA.markers.celltype.edgeR.csv")



## reorder the cluster
# Liveseq_sub.markers$cluster <- factor(Liveseq_sub.markers$cluster, levels = c(0,2,1))
# top30 <- Liveseq_sub.markers %>% group_by(cluster) %>% top_n(n = -30, wt = p_val_adj)

# top10 <- Liveseq_sub.markers %>% group_by(cluster) %>% subset(p_val_adj < 0.05) %>% top_n(n = 10, wt = avg_logFC)
# top10 <- Liveseq_sub.markers %>% subset(genesymbol != "mCherry" & genesymbol != "EGFP") %>% group_by(cluster)  %>% dplyr::arrange(p_val_adj, .by_group = TRUE) %>% top_n(n = -10, wt = p_val_adj) 

## get the top DE genes

top10.ASPC <- scRNA.markers %>% subset(genesymbol != "mCherry" & genesymbol != "EGFP") %>% subset(cluster=="1" & avg_log2FC > 0 & pct.1/pct.2 > 2)%>% top_n(n = -20, wt = p_val_adj)
top10.ASPC_DMIR <- scRNA.markers %>% subset(genesymbol != "mCherry" & genesymbol != "EGFP") %>% subset(cluster=="2" & avg_log2FC > 0 & pct.1/pct.2 > 2)%>% top_n(n = -20, wt = p_val_adj)
top10.IBA <- scRNA.markers %>% subset(genesymbol != "mCherry" & genesymbol != "EGFP") %>% subset(cluster=="3" & avg_log2FC > 0 & pct.1/pct.2 > 2)%>% top_n(n = -20, wt = p_val_adj)
top10.RAW <- scRNA.markers %>% subset(genesymbol != "mCherry" & genesymbol != "EGFP") %>% subset(cluster=="4" & avg_log2FC > 0 & pct.1/pct.2 > 2) %>% top_n(n = -20, wt = p_val_adj)
top10.RAW_LPS <- scRNA.markers %>% subset(genesymbol != "mCherry" & genesymbol != "EGFP") %>% subset(cluster=="5" & avg_log2FC > 0 & pct.1/pct.2 > 2)%>% top_n(n = -20, wt = p_val_adj)
top10 <- rbind(top10.ASPC,top10.ASPC_DMIR, top10.IBA,top10.RAW,top10.RAW_LPS  )
top10 <- top10 %>% filter(!duplicated(genesymbol))     

# top10 <- top10  %>% group_by(cluster) %>%  top_n(n = -8, wt = p_val_adj)

head(top10)
all(!duplicated(top10$genesymbol))

DefaultAssay(scRNA.integrated) <- "RNA"
## scale all the data to avoid some genes missing for heatmap
scRNA.integrated <- ScaleData(scRNA.integrated, features = rownames(scRNA.integrated) ,verbose = FALSE)


# top10.IBA <- Liveseq_sub.markers %>% subset(cluster==2 & avg_logFC > 0 )%>% top_n(n = -10, wt = p_val_adj)
# top10.RAW <- Liveseq_sub.markers %>% subset(cluster==2 & avg_logFC < 0 ) %>% top_n(n = -10, wt = p_val_adj)
# top10.RAW_LPS <- Liveseq_sub.markers %>% subset(cluster==1 & avg_logFC > 0 )%>% top_n(n = -10, wt = p_val_adj)
# top10 <- rbind(top10.IBA,top10.RAW,top10.RAW_LPS  )   
p <- DoHeatmap(scRNA.integrated, features = as.character(top10$gene) , group.colors = color.celltype ) + 
  scale_y_discrete(labels= rev(top10$genesymbol)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_x_discrete(labels=NULL)  +
  scale_fill_gradientn(colors = viridis(3)) 
p
# ggsave("03_scRNA_seq_only/heatmap_scRNA_batch_corrected_DE.pdf", width = 5, height = 5, useDingbats=FALSE)

## plot the color bar of heatmap
## Retrieve values for axis labels in ggplot2

p <- DoHeatmap(scRNA.integrated, features = as.character(top10$gene) ) + NoLegend() + 
  scale_y_discrete(labels= rev(top10$genesymbol)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

sampleorder.x <- ggplot_build(p)$layout$panel_params[[1]]$x$limits

color.x <- rep("white", length(sampleorder.x))

load("03_scRNA_seq_only/cells.types.RData")
color.x[match(ASPC.cell, sampleorder.x)]  <- color.celltype[1]
color.x[match(ASPCtreated.cell, sampleorder.x)]  <- color.celltype[2]
color.x[match(IBA.cell, sampleorder.x)]  <- color.celltype[3]
color.x[match( RAWLpsN.cell, sampleorder.x)]  <- color.celltype[4]
color.x[match(RAWLpsP.cell, sampleorder.x)]  <-  color.celltype[5]
color.x

df <- data.frame(sample.name =sampleorder.x, color.x = color.x )
df$sample.name <- factor(df$sample.name, levels = df$sample.name)


ggplot(df, aes(x=sample.name, y=1, fill=sample.name))+
  geom_bar(stat="identity", color=color.x)+
  scale_fill_manual(values=color.x) + 
  theme_nothing()
# ggsave("03_scRNA_seq_only/heatmap_scRNA_batch_corrected_DE_colorBar.pdf", width = 5, height = 0.5, useDingbats=FALSE )




#### Dimplot the DE genes

## functin for genesymbol and ensemble name conversation
gene.info <- read.table(file = file.path(root_dir, "data/mouseGeneTable87_mCherry_EGFP.txt"), sep = "\t", header = T, row.names = 1 )
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

FeaturePlot(object = scRNA.integrated, reduction = "tsne", features = symbol.to.ensembl("Tnf") )

# Dimplot the DE gene
genesToPlot <- c("Dpep1", "Thy1", "Gpx3","Sfrp2" ,"Mycn","Tnf" )
plist <- lapply(genesToPlot, function(x) FeaturePlot(object = scRNA.integrated, reduction = "umap", features = symbol.to.ensembl(x) ))
plot_grid(plotlist=plist)
# ggsave("03_scRNA_seq_only/gene_expression_umap.pdf", width = 7.2, height = 4, useDingbats=FALSE)

plist <- lapply(genesToPlot, function(x) FeaturePlot(object = scRNA.integrated, reduction =  "tsne", features = symbol.to.ensembl(x)  ))
plot_grid(plotlist=plist)
# ggsave("03_scRNA_seq_only/gene_expression_tsne.pdf", width = 7.2, height = 4, useDingbats=FALSE)





