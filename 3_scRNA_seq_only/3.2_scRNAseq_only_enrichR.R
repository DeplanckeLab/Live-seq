# install.packages("enrichR", lib = "~/NAS2/wchen/data_analysis/Resource/R_libs/")
library(enrichR, lib.loc = "~/NAS2/wchen/data_analysis/Resource/R_libs/")
library(dplyr)
library(reshape2)
library(pheatmap)
library(gplots)
library(viridis)
setwd("~/SVFASRAW/wchen/data_analysis/Live_seq/final_analysis_V3/Code_github/")

## list the enrichr database
dbs <- listEnrichrDbs()
head(dbs)

## select the database needed
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015", "Mouse_Gene_Atlas")


############### for scRNA-seq data  ##########
###########################################
scRNA.markers <- read.csv("3_scRNA_seq_only/DEs.scRNA.csv")

## the batch corretion incorporated DE analyis, refer to scRNA.markers.celltype.edgeR.csv below
# scRNA.markers <- read.csv("../LiveSeq_vsscRNAseq_9.09.21/scRNA.markers.celltype.edgeR.csv")

# DE.liveseq.0to2 <- read.csv("Liveseq_only/DEs0to2.csv")


# top100.liveseq.c0 <-  DE.liveseq.all %>% subset(cluster == 0) %>% arrange(desc(avg_logFC)) %>% top_n(100, wt = avg_logFC)
# top100.liveseq.c0 <-  DE.liveseq.0to2 %>% subset(avg_logFC >0) %>% arrange(desc(avg_logFC)) %>% top_n(100, wt = avg_logFC)
top100.liveseq.c1 <-  scRNA.markers %>% subset(cluster == 1 & pct.1/pct.2 > 3) %>% arrange(desc(avg_log2FC)) %>% top_n(100, wt = avg_log2FC)
top100.liveseq.c2 <-  scRNA.markers %>% subset(cluster == 2 & pct.1/pct.2 > 3) %>% arrange(desc(avg_log2FC)) %>% top_n(100, wt = avg_log2FC)
top100.liveseq.c3 <-  scRNA.markers %>% subset(cluster == 3 & pct.1/pct.2 > 3) %>% arrange(desc(avg_log2FC)) %>% top_n(100, wt = avg_log2FC)
top100.liveseq.c4 <-  scRNA.markers %>% subset(cluster == 4 & pct.1/pct.2 > 3) %>% arrange(desc(avg_log2FC)) %>% top_n(100, wt = avg_log2FC)
top100.liveseq.c5 <-  scRNA.markers %>% subset(cluster == 5 & pct.1/pct.2 > 3) %>% arrange(desc(avg_log2FC)) %>% top_n(100, wt = avg_log2FC)



geneset.c1 <-as.character( top100.liveseq.c1$genesymbol)
geneset.c2 <-as.character( top100.liveseq.c2$genesymbol)
geneset.c3 <-as.character( top100.liveseq.c3$genesymbol)
geneset.c4 <-as.character( top100.liveseq.c4$genesymbol)
geneset.c5 <-as.character( top100.liveseq.c5$genesymbol)


enriched.c1 <- enrichr(geneset.c1, dbs)
enriched.c1 <- lapply(enriched.c1, function(x) {x$cluster = "1"; return(x)} )

enriched.c2 <- enrichr(geneset.c2, dbs)
enriched.c2 <- lapply(enriched.c2, function(x) {x$cluster = "2"; return(x)} )

enriched.c3 <- enrichr(geneset.c3, dbs)
enriched.c3 <- lapply(enriched.c3, function(x) {x$cluster = "3"; return(x)} )

enriched.c4 <- enrichr(geneset.c4, dbs) 
enriched.c4 <- lapply(enriched.c4, function(x) {x$cluster = "4"; return(x)} )

enriched.c5 <- enrichr(geneset.c5, dbs) 
enriched.c5 <- lapply(enriched.c5, function(x) {x$cluster = "5"; return(x)} )


colnames(enriched.c2[["GO_Biological_Process_2015"]]) == colnames(enriched.c1[["GO_Biological_Process_2015"]])

######
enrich.BP <-  rbind(enriched.c1[["GO_Biological_Process_2015"]], enriched.c2[["GO_Biological_Process_2015"]], enriched.c3[["GO_Biological_Process_2015"]], enriched.c4[["GO_Biological_Process_2015"]],enriched.c5[["GO_Biological_Process_2015"]] )



# enrich.BP.top <-  enrich.BP %>% group_by(cluster) %>% top_n(n=-10, wt= Adjusted.P.value)

# enrich.BP.top.dcast <-  dcast(enrich.BP.top, Term ~ cluster, value.var = "Adjusted.P.value")

enrich.BP.dcast <-  dcast(enrich.BP, Term ~ cluster, value.var = "Adjusted.P.value")
enrich.BP.dcast[ is.na(enrich.BP.dcast) ]  <- 1  ## the NA value is filled with 1

### function to get the row index of the top 5 terms of colume 2-5  (clusters)
index_top5 <- lapply(2:6, function(col_index) { 
  as.numeric(rownames(enrich.BP.dcast[enrich.BP.dcast[[col_index]] %in% sort(enrich.BP.dcast[[col_index]], decreasing = F)[1:5], ]))
})

enrich_top5 <- enrich.BP.dcast[ c(index_top5[[1]], index_top5[[2]], index_top5[[3]], index_top5[[4]], index_top5[[5]]), ]
enrich_top5 <- subset(enrich_top5, !duplicated(Term))

rownames(enrich_top5 ) <- enrich_top5$Term


enrich_top5 <- as.matrix(enrich_top5[,2:6] )
enrich_top5 <- -log10(enrich_top5)
# subset(enrich.BP.dcast, top_n()
# 
# head( enrich.BP.dcast  %>% top_n(n=-20, wt=0), n = 10)

# rownames(enrich.BP.top.dcast ) <- enrich.BP.top.dcast$Term
# enrich.BP.top.dcast <- as.matrix(enrich.BP.top.dcast[,2:4] )
# enrich.BP.top.dcast[ is.na(enrich.BP.top.dcast) ]  <- 1

# pdf("Liveseq_only/enrichr.GoBP.pdf", width = 6, height = 4)

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(enrich_top5, n = 400)


library(dendsort)

# ## sort the row and column by cluster

# sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
# 
# mat_cluster_cols <- hclust(dist(t(enrich_top5)))
# plot(mat_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")
# 
# mat_cluster_cols <- sort_hclust(mat_cluster_cols)
# plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")
# 
# mat_cluster_rows <- sort_hclust(hclust(dist(enrich_top5)))
# plot(mat_cluster_rows, main = "Sorted Dendrogram", xlab = "", sub = "")
# dev.off()
pdf("3_scRNA_seq_only/enrichr.GoBP.pdf", width = 6, height = 3.6)
pheatmap(
  mat               = enrich_top5,
  color             = viridis(length(mat_breaks) - 1),
  breaks            = mat_breaks,
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  # cluster_cols      = mat_cluster_cols,
  # annotation_row = rownames(enrich_top5),
  cluster_rows      = F,
  cluster_cols = F,
  # annotation_col    = mat_col,
  # annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 7
  # main              = "Quantile Color Scale"
  
)

# pheatmap(  enrich_top5, scale = "none", color = viridis(100) , fontsize = 7,   breaks            = mat_breaks,)
dev.off()

write.csv(enrich_top5 , "3_scRNA_seq_only/enrichr.GoBP.csv")








#########


enrich.MGA <-  rbind(enriched.c1[["Mouse_Gene_Atlas"]], enriched.c2[["Mouse_Gene_Atlas"]], enriched.c3[["Mouse_Gene_Atlas"]], enriched.c4[["Mouse_Gene_Atlas"]],enriched.c5[["Mouse_Gene_Atlas"]] )




# enrich.BP.top <-  enrich.BP %>% group_by(cluster) %>% top_n(n=-10, wt= Adjusted.P.value)

# enrich.BP.top.dcast <-  dcast(enrich.BP.top, Term ~ cluster, value.var = "Adjusted.P.value")

enrich.MGA.dcast <-  dcast(enrich.MGA, Term ~ cluster, value.var = "Adjusted.P.value")
enrich.MGA.dcast[ is.na(enrich.MGA.dcast) ]  <- 1  ## the NA value is filled with 1

### function to get the row index of the top 5 terms of colume 2-4  (clusters)
index_top5 <- lapply(2:6, function(col_index) { 
  as.numeric(rownames(enrich.MGA.dcast[enrich.MGA.dcast[[col_index]] %in% sort(enrich.MGA.dcast[[col_index]], decreasing = F)[1:5], ]))
})

enrich_top5 <- enrich.MGA.dcast[ unique(c(index_top5[[1]], index_top5[[2]], index_top5[[3]], index_top5[[4]], index_top5[[5]])), ]
enrich_top5 <- subset(enrich_top5, !duplicated(Term))


rownames(enrich_top5 ) <- enrich_top5$Term
enrich_top5 <- as.matrix(enrich_top5[,2:6] )
enrich_top5 <- -log10(enrich_top5)
# subset(enrich.MGA.dcast, top_n()
# 
# head( enrich.MGA.dcast  %>% top_n(n=-20, wt=0), n = 10)

# rownames(enrich.MGA.top.dcast ) <- enrich.MGA.top.dcast$Term
# enrich.MGA.top.dcast <- as.matrix(enrich.MGA.top.dcast[,2:4] )
# enrich.MGA.top.dcast[ is.na(enrich.MGA.top.dcast) ]  <- 1



quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(enrich_top5, n = 400)


library(dendsort)

sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))

mat_cluster_cols <- hclust(dist(t(enrich_top5)))
plot(mat_cluster_cols, main = "Unsorted Dendrogram", xlab = "", sub = "")

mat_cluster_cols <- sort_hclust(mat_cluster_cols)
plot(mat_cluster_cols, main = "Sorted Dendrogram", xlab = "", sub = "")

mat_cluster_rows <- sort_hclust(hclust(dist(enrich_top5)))
plot(mat_cluster_rows, main = "Sorted Dendrogram", xlab = "", sub = "")
dev.off()
pdf("3_scRNA_seq_only/enrichr.MGA.pdf", width = 2.8, height = 3)
pheatmap(
  mat               = enrich_top5,
  color             = viridis(length(mat_breaks) - 1),
  breaks            = mat_breaks,
  border_color      = NA,
  show_colnames     = TRUE,
  show_rownames     = TRUE,
  # cluster_cols      = mat_cluster_cols,
  cluster_rows      = F,
  cluster_cols      = F,
  # annotation_col    = mat_col,
  # annotation_colors = mat_colors,
  drop_levels       = TRUE,
  fontsize          = 7
  # main              = "Quantile Color Scale"
  
)

# pheatmap(  enrich_top5, scale = "none", color = viridis(100) , fontsize = 7,   breaks            = mat_breaks,)
dev.off()

write.csv(enrich_top5 , "3_scRNA_seq_only/enrichr.MGA.csv")








