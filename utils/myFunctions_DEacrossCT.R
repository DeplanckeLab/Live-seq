
# select genes for DE: Express in at least 15% of one group with at least 2 counts 
select_genes <- function(seu_scRNA, clust_col_name){
  counts_diffexp <- seu_scRNA@assays$RNA@counts
  counts_diffexp <- counts_diffexp[!rownames(counts_diffexp) %in% c("EGFP", "mCherry", gene.blacklist$ensembl_gene_id),]
  
  metadata_diffexp <- seu_scRNA@meta.data[, c(clust_col_name, "Batch")]
  colnames(metadata_diffexp) <- c("Clustering", "Batch")
  
  #Filter out genes that are expressed in less than 5% of the data with counts <2
  tokeep <- filterByExpr(counts_diffexp, min.count = 2, min.prop = 0.05)
  
  counts_diffexp <- counts_diffexp[tokeep,]
  
  #Filter out if not express in at least 15% of cells in one of the 5 state with at least 1 count
  for_filt <- as.data.frame(t(as.matrix(counts_diffexp)))
  for_filt$info <- as.character(metadata_diffexp[rownames(for_filt), "Clustering"])
  colSums_perGroup <- as.data.frame(for_filt  %>% group_by(info) %>%
                                      dplyr::summarise(across(everything(), ~sum(.>2))))
  rownames(colSums_perGroup) <- colSums_perGroup$info
  colSums_perGroup <- dplyr::select(colSums_perGroup, -info)
  ncells_pergroup <- table(metadata_diffexp$Clustering)
  colSums_perGroup <- colSums_perGroup/ncells_pergroup[rownames(colSums_perGroup)]
  to_filt <- colnames(colSums_perGroup[,colSums(colSums_perGroup < 0.15) == 2])
  
  counts_diffexp <- counts_diffexp[!rownames(counts_diffexp) %in% to_filt, ]
  
  return(rownames(counts_diffexp))
}

f_edgeR <- function(seu, clust_col_name, genes_to_keep){
  counts_diffexp <- seu@assays$RNA@counts[rownames(seu) %in% genes_to_keep, ]
  
  metadata_diffexp <- seu@meta.data[, c(clust_col_name, "Batch")]
  colnames(metadata_diffexp) <- c("Clustering", "Batch")
  
  # 
  # do the actual modelling
  # https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLF.R
  dge <- edgeR::DGEList(counts_diffexp)
  dge <- edgeR::calcNormFactors(dge)
  #metadata_diffexp$cdr <- scale(Matrix::colMeans(counts_diffexp > 0))
  designmatrixm <- model.matrix(~0 + Clustering + Batch, metadata_diffexp)
  
  colnames(designmatrixm) <- c(substring(colnames(designmatrixm)[1:5], 11),
                               colnames(designmatrixm)[6:ncol(designmatrixm)])
  dge <- edgeR::estimateDisp(dge, design = designmatrixm)
  fit <- edgeR::glmQLFit(dge, design = designmatrixm)
  ave.con <- makeContrasts("ASPC_NT" = ASPC_not_treated - (ASPC_DMIR_treated + IBA_not_treated + Raw264.7_not_treated + Raw264.7_LPS_treated)/4,
                           "ASPC_DMIR" = ASPC_DMIR_treated - (ASPC_not_treated + IBA_not_treated + Raw264.7_not_treated + Raw264.7_LPS_treated)/4,
                           "IBA" =  IBA_not_treated - (ASPC_not_treated + ASPC_DMIR_treated + Raw264.7_not_treated + Raw264.7_LPS_treated)/4,
                           "Raw_NT" = Raw264.7_not_treated - (ASPC_not_treated + ASPC_DMIR_treated + IBA_not_treated + Raw264.7_LPS_treated)/4, 
                           "Raw_LPS" = Raw264.7_LPS_treated - (ASPC_not_treated + ASPC_DMIR_treated + IBA_not_treated + Raw264.7_not_treated)/4, 
                           levels=designmatrixm)
  
  de_ASPC_NT <- edgeR::glmQLFTest(fit, contrast = ave.con[, "ASPC_NT"])
  de_ASPC_DMIR <- edgeR::glmQLFTest(fit, contrast = ave.con[, "ASPC_DMIR"])
  de_IBA <- edgeR::glmQLFTest(fit, contrast = ave.con[, "IBA"])
  de_Raw_NT <- edgeR::glmQLFTest(fit, contrast = ave.con[, "Raw_NT"])
  de_Raw_LPS <- edgeR::glmQLFTest(fit, contrast = ave.con[, "Raw_LPS"])
  
  return(list(de_ASPC_NT = de_ASPC_NT,
              de_ASPC_DMIR = de_ASPC_DMIR,
              de_IBA = de_IBA,
              de_Raw_NT = de_Raw_NT,
              de_Raw_LPS = de_Raw_LPS))
}
f_adjPval <- function(x){
  x$table$padj <- p.adjust(x$table$PValue, method = "BH")
  x$table$Name <- data.annot[rownames(x$table), "Name"]
  x$table <- x$table[order(x$table$logFC, decreasing = T),]
  return(x)
}

# Calculate the percentage of cells expressing each genes 
calculate_pcts <- function(group, seu, de_res){
  convert <- c(de_ASPC_NT = "ASPC_not_treated",  de_ASPC_DMIR = "ASPC_DMIR_treated", 
               de_IBA = "IBA_not_treated",
               de_Raw_NT = "Raw264.7_not_treated", de_Raw_LPS = "Raw264.7_LPS_treated")
  ident.1 <- convert[group]
  ident.2 <- convert[names(convert) != group]
  
  pct.1 <- rowSums(as.matrix(seu@assays$RNA@counts[,seu$celltype_treatment == ident.1] > 1))/sum(seu$celltype_treatment == ident.1)
  pct.2 <- rowSums(as.matrix(seu@assays$RNA@counts[,seu$celltype_treatment %in% ident.2] > 1))/sum(seu$celltype_treatment %in% ident.2)
  
  de_res$table$pct.1 <- pct.1[rownames(de_res$table)]
  de_res$table$pct.2 <- pct.2[rownames(de_res$table)]
  de_res$table$high_pct <- F
  de_res$table$high_pct[de_res$table$pct.1 > 0.15 | de_res$table$pct.2 > 0.15] <- T
  
  return(de_res)
}

f_is_de <- function(de_res){
  de_res$table$is_de <- F
  de_res$table[rownames(de_res$table %>% filter(padj < 0.05 & abs(logFC) > 1 & high_pct == T)), "is_de"] <- T
  return(de_res)
}

f_add <- function(which_d, d){
  out <- d[[which_d]]
  out$cat <- which_d
  out$Ensembl <- rownames(out)
  return(out)
}

f_plot <- function(DE_res_live, DE_res_sc){
  
  my_genes <- intersect(rownames(DE_res_live), rownames(DE_res_sc))
  cat(length(my_genes))
  cat("\n")
  my_genes <- my_genes[!my_genes %in% c("EGFP", "mCherry", gene.blacklist$ensembl_gene_id)]
  df <- data.frame(FC_Liveseq = DE_res_live[my_genes, "logFC"],
                   FC_scRNA = DE_res_sc[my_genes, "logFC"],
                   ens = my_genes,
                   gene = data.annot[my_genes, "Name"])
  rownames(df) <- df$ens
  
  df$DE <- "Rest"
  df$DE[df$ens %in% rownames(filter(as.data.frame(DE_res_live), padj < 0.05 & abs(logFC) >1 & high_pct == T))] <- "DE Liveseq"
  df$DE[df$ens %in% rownames(filter(as.data.frame(DE_res_sc), padj < 0.05 & abs(logFC) > 1 & high_pct == T))] <- "DE scRNAseq"
  df$DE[df$ens %in% intersect(rownames(filter(as.data.frame(DE_res_sc), padj < 0.05 & abs(logFC) >1 & high_pct == T)),
                              rownames(filter(as.data.frame(DE_res_live), padj < 0.05 & abs(logFC) >1 & high_pct == T)))] <- "DE both"
  df$DE <- factor(df$DE,levels = c("Rest", "DE Liveseq", "DE scRNAseq", "DE both"))
  
  df <- df[order(df$DE),]
  p <- ggplot(df, aes(x = FC_scRNA, y = FC_Liveseq, label = gene, col = DE)) + 
    geom_point(size = 0.3) + 
    scale_color_manual(values = c( "gray80", "#B9E3BB", "#7ACCC3", "#2C8BBD")) + 
    geom_hline(yintercept = -1, linetype = "dashed", col = "gray47") + 
    geom_hline(yintercept = 1, linetype = "dashed", col = "gray47") + 
    geom_vline(xintercept = -1, linetype = "dashed", col = "gray47") +
    geom_vline(xintercept = 1, linetype = "dashed", col = "gray47") +
    mashaGgplot2Theme + 
    stat_smooth(data = df, 
                mapping = aes(x = FC_scRNA, y = FC_Liveseq),
                method='lm', formula= y~x, col = "gray60", se = T) +
    stat_smooth(data = df %>% filter(DE %in% c("DE scRNAseq", "DE both")), 
                mapping = aes(x = FC_scRNA, y = FC_Liveseq),
                method='lm', col = "#8dd3c7", se = T, fill = "gray80") +
    stat_smooth(data = df %>% filter(DE %in% c("DE both")), 
                mapping = aes(x = FC_scRNA, y = FC_Liveseq),
                method='lm', col = "#253494", se = T, fill = "gray80") + 
    #geom_abline(slope = 1, intercept = 0) +
    xlab("logFC based on scRNA-seq") + 
    ylab("logFC based on Live-seq")
  
  lm_ALL <- summary(lm(FC_Liveseq~FC_scRNA, data = df))
  lm_DEsc <- summary(lm(FC_Liveseq~FC_scRNA, data = filter(df, DE %in% c("DE both", "DE scRNAseq"))))
  lm_DEboth <- summary(lm(FC_Liveseq~FC_scRNA, data = filter(df, DE %in% c("DE both"))))
  
  return(list(p = p,
              lm_AllGenes=lm_ALL,
              lm_DEsc=lm_DEsc,
              lm_DEboth=lm_DEboth))
  
}

xtab_set <- function(live,sc){
  both    <-  union(live,sc)
  inlive     <-  both %in% live
  insc     <-  both %in% sc
  return(table(inlive,insc))
}
f_comp <- function(x){
  out <- reshape2::melt(xtab_set((liveseq_edgeR_res_filt %>% filter(cat == x))[,"Ensembl"], 
                                 (scRNA_edgeR_res_filt %>% filter(cat == x))[,"Ensembl"]))
  out$cat <- x
  return(out)
}


