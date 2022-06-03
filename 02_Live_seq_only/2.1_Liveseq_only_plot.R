################################################################
#                                                              #
#                         Plot Live-seq                        #
#                                                              #
################################################################

### Questions: https://github.com/DeplanckeLab/Live-seq/issues
### Date: 2022-03-06
### Datasets: Live-seq
### Goal: Making some summary plots of the Live-seq data

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

## load Liveseq object
Liveseq <-
  readRDS(file.path(root_dir, "02_Live_seq_only/Liveseq_only.rds"))

## mean of nFeature = 4112.35
mean(Liveseq$nFeature_RNA)

#### get the default color of ggplot2. n is the number of colors needed.
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = gg_color_hue(5)
cols   ## "#F8766D" "#00BA38" "#619CFF"

color.batch <- cols
## color for celle_treatment
color.celltype <-
  c("lightgoldenrod4",
    "#dfc27d" ,
    "steelblue1"  ,
    "#80cdc1",
    "#018571") # ,darkslateblue,  tan3

# Liveseq@meta.data$celltype_treatment  <-  droplevels(Liveseq@meta.data$celltype_treatment  )
# Liveseq@meta.data$celltype_treatment  <- factor(Liveseq@meta.data$celltype_treatment , levels = c("ASPC_not_treated", "IBA_not_treated","Raw264.7_not_treated", "Raw264.7_LPS_treated"))


## plot the number of samples in each step
df <-
  data.frame(
    catalog = c("Total", "Library preparation", "Passing QC"),
    number = c(641, 588, ncol(Liveseq))
  )
df$catalog <- factor(df$catalog , levels = df$catalog)
p <-
  ggplot(df, aes(catalog, number)) + geom_bar(stat = "identity") + xlab(label = "") + ylim(0, 500)
p <- p + theme(
  axis.ticks = element_line(colour = "black", size = 0.2),
  axis.text.x = element_text(
    angle = 45,
    colour = "black",
    hjust = 1
  ),
  axis.text.y = element_text(colour = "black")
) +
  stat_summary(
    geom = 'text',
    label = df$number,
    fun.y = max,
    vjust = -1
  )

p <-
  p +  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 800))
p <- p + ylab("Number of cells")
p
# ggsave("sample_number.pdf", width = 1.3, height = 2, useDingbats=FALSE)


## plot the basic QC, group by batch or sample type
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
    # p <- myVlnPlot(Liveseq, x, "sampling_type", cols = c("grey","#00BA38","#619CFF"))
    p <-
      ggplot(Liveseq@meta.data, aes_string(x = "sampling_type", y = i)) +
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
# ggsave("QC_noGroup.pdf",width = 10, height = 1.6, useDingbats=FALSE )



## plot the metrics
# per celltype
ylab.metrics <-
  c(
    "Number of reads",
    "Fraction of input reads",
    "Fraction of uniquely mapped reads",
    "Total count of all genes",
    "Number of detected genes",
    "Fraction of nCount (%)",
    "Fraction of nCount (%)",
    "Fraction of nCount (%)"
  )
names(ylab.metrics) <-
  c(
    "input.reads",
    "uniquely.mapped.rate",
    "intron.mapped.rate",
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt",
    "percent.rRNA",
    "percent.protein"
  )
title.metrics <-
  c(
    "Input reads",
    "Uniquely mapped reads",
    "Reads mapped in exon",
    "nCount",
    "nGene",
    "Percent MT",
    "Percent rRNA",
    "Percent Protein"
  )
names(title.metrics) <-
  c(
    "input.reads",
    "uniquely.mapped.rate",
    "intron.mapped.rate",
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt",
    "percent.rRNA",
    "percent.protein"
  )

plist <-
  lapply(c(
    "input.reads",
    "uniquely.mapped.rate",
    "intron.mapped.rate",
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt",
    "percent.rRNA",
    "percent.protein"
  ), function(x) {
    p <-
      VlnPlot(
        object = Liveseq,
        features = x,
        group.by = "celltype_treatment",
        cols = color.celltype
      )
    p <-
      p + scale_x_discrete(labels = c("ASPC_Pre", "ASPC_Post", "IBA", "RAW_Mock", "RAW_LPS"))
    p <- p + xlab("")
    p <- p + ylab(ylab.metrics[x]) + ggtitle(title.metrics[x])
    return(p)
  })
plot_grid(plotlist = plist,
          align = "hv",
          axis = "lrtb")
# ggsave("QC_celltype.pdf", width = 5, height = 6, useDingbats=FALSE )

# per batch
plist <-
  lapply(c(
    "input.reads",
    "uniquely.mapped.rate",
    "intron.mapped.rate",
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt",
    "percent.rRNA",
    "percent.protein"
  ), function(x) {
    p <-
      VlnPlot(
        object = Liveseq,
        features = x,
        group.by = "Batch",
        cols = color.batch
      )
    p <-
      p + scale_x_discrete(labels = c("ASPC_Pre", "ASPC_Post", "IBA", "RAW_Mock", "RAW_LPS"))
    p <- p + xlab("")
    p <- p + ylab(ylab.metrics[x]) + ggtitle(title.metrics[x])
    return(p)
  })
# plist <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rRNA", "percent.protein", "uniquely.mapped.rate", "intron.mapped.rate"), function(x) myVlnPlot(Liveseq, x,"Batch", cols = color.batch ))
plot_grid(plotlist = plist)
# ggsave("QC_batch.pdf", width = 5, height = 4.5, useDingbats=FALSE )

## plot scatter plot
p1 <-
  FeatureScatter(
    Liveseq,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    group.by = "Batch",
    pt.size = 0.3,
    cols = color.batch
  )
# p1 <- mytheme(p1)

p2 <-
  FeatureScatter(
    Liveseq,
    feature1 = "cDNA_concentration",
    feature2 = "nFeature_RNA",
    group.by = "Batch",
    pt.size = 0.3,
    cols = color.batch
  )
# p2 <- mytheme(p2)

p3 <-
  FeatureScatter(
    Liveseq,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    group.by = "celltype_treatment",
    pt.size = 0.3,
    cols = color.celltype
  )
# p3 <- mytheme(p3)

p4 <-
  FeatureScatter(
    Liveseq,
    feature1 = "cDNA_concentration",
    feature2 = "nFeature_RNA",
    group.by = "celltype_treatment",
    pt.size = 0.3,
    cols = color.celltype
  )
# p4 <- mytheme(p4)

plot_grid(p1, p2, p3, p4, align = "hv", axis = "lrtb")

# ggsave("scatter_cDNA,nCount_nFeature.pdf", width = 7, height = 4, useDingbats=FALSE)


# predict the amount of input RNA based on the cDNA yield with experimental stardard curve
# y = 5.651x + 7.911   (y is the cDNA yield (pg) and x is the input RNA pg) , this is calculated by smart2_optimization.R, based on stardar curve
# the volume to elute cDNA is 15 ul and the concentration has already normalized to control (0 pg RNA) in each experiment. So only the slope is taken.

Liveseq$input.RNA <- (Liveseq$cDNA_concentration * 15) / 5.651

ggplot(Liveseq@meta.data, aes(input.RNA, nFeature_RNA)) + geom_point() + geom_smooth(method = lm)
## remove the outliers with input.RNA > 5
fit.inputRNA <-
  lm(nFeature_RNA ~ input.RNA   ,
     data = Liveseq@meta.data %>% filter(input.RNA < 5))

ggplot(Liveseq@meta.data %>% filter(input.RNA < 5),
       aes(input.RNA, nFeature_RNA)) + geom_point() + geom_smooth(method = lm) +
  annotate(
    "text",
    x = 2,
    y = 9000,
    label = paste(fit.inputRNA$coefficients[2], "x+",    fit.inputRNA$coefficients[1])
  )
# ggsave("inputRNA.nGene.curve.pdf")


## explore the cells misclustered
Liveseq$idents <- Idents(Liveseq)
misClusterASPC_Mock <-
  Liveseq@meta.data %>% filter(idents != "1" &
                                 (celltype_treatment == "ASPC_not_treated"))
misClusterASPC_DMIR <-
  Liveseq@meta.data %>% filter(idents != "2" &
                                 (celltype_treatment == "ASPC_DMIR_treated"))
misClusterIBA <-
  Liveseq@meta.data %>% filter(idents != "3" &
                                 celltype_treatment == "IBA_not_treated")
misClusterRAW_Mock <-
  Liveseq@meta.data %>% filter(idents != "4" &
                                 celltype_treatment == "Raw264.7_not_treated")
misClusterRAW_LPS <-
  Liveseq@meta.data %>% filter(idents != "5" &
                                 celltype_treatment == "Raw264.7_LPS_treated")

misCLuster.all <-
  rbind(
    misClusterASPC_Mock,
    misClusterASPC_DMIR,
    misClusterIBA,
    misClusterRAW_Mock,
    misClusterRAW_LPS
  )


Liveseq$is.clusterOK <- "yes"
Liveseq$is.clusterOK [Liveseq$sample_name %in% misCLuster.all$sample_name] <-
  "no"

# cells_misClustered per nFeature
p.clusterOK.nFeature <-
  t.test(Liveseq@meta.data[Liveseq@meta.data$is.clusterOK == "no" , "nFeature_RNA"] , Liveseq@meta.data[Liveseq@meta.data$is.clusterOK == "yes" , "nFeature_RNA"])

ggplot(Liveseq@meta.data, aes(is.clusterOK, nFeature_RNA)) + geom_violin() + geom_jitter()


p <- ggplot(Liveseq@meta.data, aes(is.clusterOK, nFeature_RNA))  +
  geom_violin(trim = TRUE,  fill = "gray") +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "pointrange",
    color = viridis(10)[1]
  ) +
  theme_classic() +
  annotate(
    "text",
    x = 1,
    y = 9000,
    label = paste("t test, p ",  p.clusterOK.nFeature$p.value)
  )

p <- p + geom_jitter(shape = 16, position = position_jitter(0.2))
p
# ggsave("cells_misClustered_nFeature.pdf", width = 2, height = 3, useDingbats=FALSE)


# cells_misClustered per input.RNA
p.clusterOK.input.RNA <-
  t.test(Liveseq@meta.data[Liveseq@meta.data$is.clusterOK == "no" , "input.RNA"] , Liveseq@meta.data[Liveseq@meta.data$is.clusterOK == "yes" , "input.RNA"])

p <- ggplot(Liveseq@meta.data, aes(is.clusterOK, input.RNA))  +
  geom_violin(trim = TRUE,  fill = "gray") +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "pointrange",
    color = viridis(10)[1]
  ) +
  theme_classic() +
  ylim(0, 8) +
  annotate(
    "text",
    x = 1,
    y = 7,
    label = paste("t test, p ",  p.clusterOK.input.RNA$p.value)
  )

p <- p + geom_jitter(shape = 16, position = position_jitter(0.2))
p
# ggsave("cells_misClustered_input.RNA.pdf", width = 2, height = 3, useDingbats=FALSE)



### Hierarchical clustering
library(dendextend)
Liveseq.rename <-
  RenameCells(Liveseq, new.names = as.character(Liveseq@meta.data$original_sample_name))
Liveseq.rename <-
  RenameCells(Liveseq, new.names = as.character(Liveseq@meta.data$sample_name))

dend <-
  Liveseq.rename@assays$RNA@data[VariableFeatures(Liveseq.rename),] %>% # get the most variable gene expression data
  as.matrix %>% # change to matrix
  t %>%        # transform
  scale %>% # Scale the data
  dist %>% # calculate a distance matrix,
  hclust(method = "ward.D2") %>% # Hierarchical clustering
  as.dendrogram # Turn the object into a dendrogram.

# arrange(df, factor(colume.name, levels = "target order" ) can be used for order the df according to the column in a customized order
df.temp <-  Liveseq.rename@meta.data %>%
  dplyr::arrange(factor(sample_name, levels = labels(dend)))  # labels(dend) to get the label in order

mycolor <-
  droplevels.factor(df.temp$celltype_treatment)   ## color according to the level
levels(mycolor) <- color.celltype

# Change color and size for labels
par()              # view current settings
opar <- par()      # make a copy of current settings
par(cex = 0.5) # red x and y labels

pdf(
  "02_Live_seq_only/Liveseq_hierachical_clustering.pdf",
  width = 7,
  height = 3
)
dend %>%
  set("labels_col", as.character(mycolor)) %>% # change color
  set("labels_cex", 0.5) %>% # Change size
  plot(main = "Hierarchical clustering") # plot
dev.off()


## plot the color bar
df <-
  data.frame(sample.name = as.character(seq(from = 1, to = length(mycolor))), mycolor = as.character(mycolor))
df$sample.name <- factor(df$sample.name, levels = df$sample.name)
ggplot(df, aes(x = sample.name, y = 1, fill = sample.name)) +
  geom_bar(stat = "identity", color = as.character(mycolor)) +
  scale_fill_manual(values = as.character(mycolor)) +
  theme_nothing()
ggsave("02_Live_seq_only/Liveseq_hierachical_clustering_colorBar.pdf")


par(opar)          # restore original settings



##### save the cell name of each batch and celltype_treatmment
ASPC.cell <-
  colnames(subset(Liveseq, celltype_treatment == "ASPC_not_treated"))
ASPCtreated.cell <-
  colnames(subset(Liveseq, celltype_treatment == "ASPC_DMIR_treated"))
IBA.cell <-
  colnames(subset(Liveseq, celltype_treatment == "IBA_not_treated"))
RAWLpsN.cell <-
  colnames(subset(Liveseq, celltype_treatment %in% c("Raw264.7_not_treated")))
RAWLpsP.cell <-
  colnames(subset(Liveseq, celltype_treatment %in% c("Raw264.7_LPS_treated")))

batch5.cell <- colnames(subset(Liveseq, Batch == "5Reseq"))
batch6.cell <- colnames(subset(Liveseq, Batch == "6"))
batch7.cell <- colnames(subset(Liveseq, Batch == "7"))
batch8_8.cell <- colnames(subset(Liveseq, Batch == "8_8"))
batch9_4.cell <- colnames(subset(Liveseq, Batch == "9_4"))

## save the name of each cell type and treatment
save(
  ASPC.cell,
  ASPCtreated.cell,
  IBA.cell,
  RAWLpsN.cell,
  RAWLpsP.cell,
  batch5.cell,
  batch6.cell,
  batch7.cell,
  batch8_8.cell,
  batch9_4.cell ,
  file = "02_Live_seq_only/cells.types.batch5,6,7,8_8,9_4.RData"
)
write(
  c(
    batch5.cell,
    " ",
    batch6.cell,
    " ",
    batch7.cell,
    " ",
    batch8_8.cell,
    " ",
    batch9_4.cell
  ),
  "02_Live_seq_only/cells.batch5,6,7,8_8,9_4.txt",
  sep = " "
)
write(
  c(
    ASPC.cell,
    " ",
    ASPCtreated.cell,
    " ",
    IBA.cell,
    " " ,
    RAWLpsN.cell,
    " ",
    RAWLpsP.cell
  ),
  "02_Live_seq_only/cells.types5,6,7,8_8,9_4.txt",
  sep = " "
)


### compare the cluster with ground trust.
### two categorical variable stacked barplot with label
Liveseq$ident <- Idents(Liveseq)

# with absolute cell number
p <- ggplot(Liveseq@meta.data, aes(x = ident, fill = celltype_treatment)) +
  geom_bar() +
  geom_text(
    aes(label = stat(count)),
    stat = "count",
    position = position_stack(0.5),
    size = 3
  ) +
  xlab("cluster") + scale_fill_manual(values = color.celltype) +
  ylim(c(0, 150))
## remove the grid
p <-
  p +  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) +
  ylab("Number of cells") +
  xlab("Cluster")
p
# ggsave("02_Live_seq_only/Liveseq_evaluate_cluster.pdf", width = 3, height = 2)

# with percentage
p <- ggplot(Liveseq@meta.data, aes(x = ident, fill = celltype_treatment)) +
  geom_bar(position = "fill", stat = "count") +
  geom_text(
    aes(label = stat(count)),
    stat = "count",
    position = position_stack(0.5),
    size = 3
  ) +
  xlab("cluster") + scale_fill_manual(values = color.celltype) +
  ylim(c(0, 1))
## remove the grid
p <-
  p +  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  ylab("Number of cells") +
  xlab("Cluster")

p
# ggsave("02_Live_seq_only/Liveseq_evaluate_cluster_percent.pdf", width = 3, height = 2)





## umap plots
p1 <-
  DimPlot(
    Liveseq,
    reduction = "umap",
    cols = color.celltype,
    label = TRUE,
    label.size = 2 ,
    repel = T,
    pt.size = 0.3
  ) + NoLegend()

p2 <-
  DimPlot(
    Liveseq,
    reduction = "umap",
    group.by = "celltype_treatment",
    label = TRUE,
    label.size = 2 ,
    repel = T,
    pt.size = 0.3,
    cols = color.celltype
  ) + NoLegend()

p3 <-
  FeaturePlot(Liveseq,
              reduction  = "umap",
              features = "nFeature_RNA",
              pt.size = 0.3)

p4 <-
  DimPlot(
    Liveseq,
    reduction  = "umap",
    group.by = "Batch",
    cols = color.batch,
    pt.size = 0.3
  )

p5 <-
  FeaturePlot(Liveseq,
              reduction  = "umap",
              features = "nCount_RNA",
              pt.size = 0.3)

plot_grid(p1, p2, p3, p4, p5, align = "hv", axis = "lrtb")
# ggsave("02_Live_seq_only/cluster_umap.pdf", width = 8.2, height = 4, useDingbats=F)

## tsne plots
p1 <-
  DimPlot(
    Liveseq,
    reduction = "tsne",
    cols = color.celltype,
    label = TRUE,
    label.size = 2 ,
    repel = T,
    pt.size = 0.3
  ) + NoLegend()

p2 <-
  DimPlot(
    Liveseq,
    reduction = "tsne",
    group.by = "celltype_treatment",
    label = TRUE,
    label.size = 2 ,
    repel = T,
    pt.size = 0.3,
    cols = color.celltype
  ) + NoLegend()

p3 <-
  FeaturePlot(Liveseq,
              reduction  = "tsne",
              features = "nFeature_RNA",
              pt.size = 0.3)

p4 <-
  DimPlot(
    Liveseq,
    reduction  = "tsne",
    group.by = "Batch",
    cols = color.batch,
    pt.size = 0.3
  )

p5 <-
  FeaturePlot(Liveseq,
              reduction  = "tsne",
              features = "nCount_RNA",
              pt.size = 0.3)

p6 <-
  DimPlot(
    Liveseq,
    reduction  = "tsne",
    group.by = "Cell_type",
    cols = c("wheat3", "steelblue1", "darkcyan", "cyan3"),
    pt.size = 0.3
  )
plot_grid(p1, p2, p3, p4, p5, p6, align = "hv", axis = "lrtb")
# ggsave("02_Live_seq_only/cluster_tsne.pdf", width = 9.2, height = 4, useDingbats=F)




## pca plots

p1 <-
  DimPlot(
    Liveseq,
    reduction = "pca",
    cols = color.celltype,
    label = TRUE,
    label.size = 2 ,
    repel = T,
    pt.size = 0.3
  ) + NoLegend()

p2 <-
  DimPlot(
    Liveseq,
    reduction = "pca",
    group.by = "celltype_treatment",
    label = TRUE,
    label.size = 2 ,
    repel = T,
    pt.size = 0.3,
    cols = color.celltype
  ) + NoLegend()

p3 <-
  FeaturePlot(Liveseq,
              reduction  = "pca",
              features = "nFeature_RNA",
              pt.size = 0.3)

p4 <-
  DimPlot(
    Liveseq,
    reduction  = "pca",
    group.by = "Batch",
    cols = color.batch,
    pt.size = 0.3
  )

p5 <-
  FeaturePlot(Liveseq,
              reduction  = "tsne",
              features = "nCount_RNA",
              pt.size = 0.3)

plot_grid(p1, p2, p3, p4, p5, align = "hv", axis = "lrtb")
# ggsave("02_Live_seq_only/cluster_pca.pdf", width = 8.2, height = 4, useDingbats=F)


### cell cycle plot

p2 <-
  DimPlot(
    Liveseq,
    reduction = "tsne",
    group.by = "Phase",
    label = TRUE,
    label.size = 2 ,
    repel = T,
    pt.size = 0.3
  )
p3 <-
  FeaturePlot(Liveseq,
              reduction  = "tsne",
              features = "S.Score",
              pt.size = 0.3)
p4 <-
  FeaturePlot(Liveseq,
              reduction  = "tsne",
              features = "G2M.Score",
              pt.size = 0.3)

plot_grid(p2, p3, p4, align = "hv", axis = "lrtb")
# ggsave("02_Live_seq_only/cluster_cellCycle.pdf", width = 5.5, height = 4, useDingbats=F)



## plot heatmap of top DE of Live-seq data
Liveseq.markers <-
  read.csv(file.path(root_dir, "02_Live_seq_only/DEs.Liveseq_only.csv"))
# Liveseq.markers <- read.csv("../LiveSeq_vsscRNAseq_9.09.21/Liveseq.markers.celltype.edgeR.csv")  ## DE from edgeR


######## get the top genes
top10.ASPC <-
  Liveseq.markers %>% subset(cluster == "1" &
                               avg_log2FC > 0) %>% top_n(n = 20, wt = avg_log2FC)
top10.ASPC_DMIR <-
  Liveseq.markers %>% subset(cluster == "2" &
                               avg_log2FC > 0) %>% top_n(n = 20, wt = avg_log2FC)
top10.IBA <-
  Liveseq.markers %>% subset(cluster == "3" &
                               avg_log2FC > 0) %>% top_n(n = 20, wt = avg_log2FC)
top10.RAW <-
  Liveseq.markers %>% subset(cluster == "4" &
                               avg_log2FC > 0) %>% top_n(n = 20, wt = avg_log2FC)
top10.RAW_LPS <-
  Liveseq.markers %>% subset(cluster == "5" &
                               avg_log2FC > 0) %>% top_n(n = 20, wt = avg_log2FC)
top10 <-
  rbind(top10.ASPC,
        top10.ASPC_DMIR,
        top10.IBA,
        top10.RAW,
        top10.RAW_LPS)

top10 <- top10 %>% filter(genesymbol != c("mCherry"))

top10 <- top10 %>% filter(!duplicated(genesymbol))

top10 <-
  top10  %>% group_by(cluster) %>%  top_n(n = 8, wt = avg_log2FC)
head(top10)
all(!duplicated(top10$genesymbol))




p <-
  DoHeatmap(Liveseq,
            features = as.character(top10$gene),
            group.colors = color.celltype) +
  scale_y_discrete(labels = rev(top10$genesymbol)) +
  scale_fill_gradientn(colors = viridis(3))
p
# ggsave("heatmap_Liveseq_DE.pdf", width = 5, height = 5, useDingbats=FALSE )

## plot the color bar of heatmap
## Retrieve values for axis labels in ggplot2
p <-
  DoHeatmap(Liveseq, features = as.character(top10$gene)) + NoLegend() +
  scale_y_discrete(labels = rev(top10$genesymbol)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
sampleorder.x <- ggplot_build(p)$layout$panel_params[[1]]$x$limits

color.x <- rep("white", length(sampleorder.x))

# load cells of each celltype_treament
load("02_Live_seq_only/cells.types.batch5,6,7,8_8,9_4.RData")
color.x[match(ASPC.cell, sampleorder.x)]  <- color.celltype[1]
color.x[match(ASPCtreated.cell, sampleorder.x)]  <-
  color.celltype[2]
color.x[match(IBA.cell, sampleorder.x)]  <- color.celltype[3]
color.x[match(RAWLpsN.cell, sampleorder.x)]  <- color.celltype[4]
color.x[match(RAWLpsP.cell, sampleorder.x)]  <- color.celltype[5]
color.x
match(sampleorder.x, RAWLpsN.cell)

## careful about the factor, the order mater
df <- data.frame(sample.name = sampleorder.x, cell.color = color.x)
df$sample.name <- factor(df$sample.name, levels = df$sample.name)
df$cell.color <- as.character(df$cell.color)

ggplot(df, aes(x = sample.name, y = 1, fill = sample.name)) +
  geom_bar(stat = "identity", color = color.x) +
  scale_fill_manual(values = color.x) +
  theme_nothing()

# ggsave("heatmap_Liveseq_DE_colorBar.pdf", width = 5, height = 0.5, useDingbats=FALSE )




#### Dimplot the DE genes

## functin for genesymbol and ensemble name conversation
gene.info <-
  read.table(
    file = file.path(root_dir, "data/mouseGeneTable87_mCherry_EGFP.txt"),
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

FeaturePlot(object = Liveseq,
            reduction = "umap",
            features = symbol.to.ensembl("Tnf"))

# Dimplot the DE gene
genesToPlot <- c("Dpep1", "Thy1", "Gpx3", "Sfrp2" , "Mycn", "Tnf")
plist <-
  lapply(genesToPlot, function(x)
    FeaturePlot(
      object = Liveseq,
      reduction = "umap",
      features = symbol.to.ensembl(x)
    ))
plot_grid(plotlist = plist)
# ggsave("gene_expression_umap.pdf", width = 7.2, height = 4, useDingbats=FALSE)

plist <-
  lapply(genesToPlot, function(x)
    FeaturePlot(
      object = Liveseq,
      reduction =  "tsne",
      features = symbol.to.ensembl(x)
    ))
plot_grid(plotlist = plist)
# ggsave("gene_expression_tsne.pdf", width = 7.2, height = 4, useDingbats=FALSE)



#### label the double extraction samples

Liveseq$double_extraction <- factor(Liveseq$double_extraction)

plist <- list()
for (i in levels(Liveseq$double_extraction)) {
  cells.selected <-
    names(Liveseq$double_extraction[Liveseq$double_extraction == i]) ## select double extraction cells
  treatment <-
    Liveseq$treatment[Liveseq$double_extraction == i]           ## get the LPS treatment info of double extraction cells
  double_extraction_order <-
    Liveseq$double_extraction_order[Liveseq$double_extraction == i]  ## get the sampling order (1st or 2nd sampling)
  ilabel <-
    paste(Liveseq@meta.data[cells.selected, "original_sample_name"],
          treatment,
          double_extraction_order,
          collapse = " ")
  
  p <-
    DimPlot(
      Liveseq,
      reduction = "tsne",
      group.by = "double_extraction",
      cells.highlight = cells.selected,
      sizes.highlight = 1.5,
      pt.size = 0.3
    ) +
    ggtitle(ilabel)
  plist[[i]] <- p
}

p1 <-
  DimPlot(
    Liveseq,
    reduction = "tsne",
    cols = color.celltype,
    label = TRUE,
    label.size = 2 ,
    repel = T,
    pt.size = 0.3
  ) + NoLegend()
plist[[1]] <- p1

## plot every 16 pannels
plot_grid(plotlist = plist[1:16],
          align = "hv",
          axis = "lrtb")
# ggsave("02_Live_seq_only/double_extraction1.tsne.pdf", width = 12, height = 8,  useDingbats=FALSE)

plot_grid(plotlist = plist[17:32],
          align = "hv",
          axis = "lrtb")
# ggsave("02_Live_seq_only/double_extraction2.tsne.pdf", width = 12, height = 8,  useDingbats=FALSE)

plot_grid(plotlist = plist[33:48],
          align = "hv",
          axis = "lrtb")
# ggsave("02_Live_seq_only/double_extraction3.tsne.pdf", width = 12, height = 8,  useDingbats=FALSE)

plot_grid(
  plotlist = plist[49:length(plist)],
  align = "hv",
  axis = "lrtb",
  cols = 4
)
# ggsave("02_Live_seq_only/double_extraction4.tsne.pdf", width = 12, height = 4,  useDingbats=FALSE)




### evaluate the sequencing depth effect with downsampling data

downsample.new <- readRDS("02_Live_seq_only/downsample.new.rds")


## plot the QC of downsampled data
df <- downsample.new[[3]]@meta.data
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
    # p <- myVlnPlot(Liveseq, x, "sampling_type", cols = c("grey","#00BA38","#619CFF"))
    p <- ggplot(df, aes_string(x = "sampling_type", y = i)) +
      geom_violin(trim = T,  fill = "gray") +
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

# ggsave("02_Live_seq_only/downsampled.QC.pdf",width = 10, height = 1.6, useDingbats=FALSE )


plist <- lapply(1:length(downsample.new), function(i) {
  p <-
    DimPlot(
      downsample.new[[i]],
      reduction = "tsne",
      group.by = "celltype_treatment",
      label = TRUE,
      label.size = 2 ,
      repel = T,
      pt.size = 0.3,
      cols = color.celltype
    ) +
    NoLegend() +
    ggtitle(mean(downsample.new[[i]]@meta.data$nCount_RNA))
})

plot_grid(plotlist = plist)
# ggsave("02_Live_seq_only/Liveseq_cluster_downsampling.pdf", width = 8, height = 6)
