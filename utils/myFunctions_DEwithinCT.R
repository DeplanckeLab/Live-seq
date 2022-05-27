f_edgeR_CT <- function(seu, genes_to_keep){
  
  counts_diffexp <- seu@assays$RNA@counts[rownames(seu) %in% genes_to_keep, ]
  metadata_diffexp <- seu@meta.data[, c("celltype_treatment", "Batch")]
  
  # do the actual modelling
  # https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLF.R
  dge <- edgeR::DGEList(counts_diffexp)
  dge <- edgeR::calcNormFactors(dge)
  #metadata_diffexp$cdr <- scale(Matrix::colMeans(counts_diffexp > 0))
  if(length(table(metadata_diffexp$Batch)) == 1){
    designmatrixm <- model.matrix(~0 + celltype_treatment, metadata_diffexp)}else{
      designmatrixm <- model.matrix(~0 + celltype_treatment + Batch, metadata_diffexp)}
  
  dge <- edgeR::estimateDisp(dge, design = designmatrixm)
  fit <- edgeR::glmQLFit(dge, design = designmatrixm)
  
  colnames(designmatrixm)[grep("not_treated",  colnames(designmatrixm))] <- "NOTTreated"
  colnames(designmatrixm)[grep("_treated",  colnames(designmatrixm))] <- "Treated"
  
  ave.con <- makeContrasts("TreatedvsNOT" = Treated - NOTTreated, 
                           levels=designmatrixm)
  
  out <- edgeR::glmQLFTest(fit, contrast = ave.con[, "TreatedvsNOT"])
  out$table$padj <- p.adjust(out$table$PValue, method = "BH")
  return(out)
}

f_plot_CT <- function(DE_res_live, DE_res_sc){
  
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
  df$DE[df$ens %in% rownames(filter(as.data.frame(DE_res_live), padj < 0.05 & abs(logFC) > 1 & high_pct == T))] <- "DE Liveseq"
  df$DE[df$ens %in% rownames(filter(as.data.frame(DE_res_sc), padj < 0.05 & abs(logFC) > 1))] <- "DE scRNAseq"
  df$DE[df$ens %in% intersect(rownames(filter(as.data.frame(DE_res_sc), padj < 0.05 & abs(logFC) > 1 )),
                              rownames(filter(as.data.frame(DE_res_live), padj < 0.05 & abs(logFC) > 1 & high_pct == T)))] <- "DE both"
  df$DE <- factor(df$DE,levels = c("Rest", "DE Liveseq", "DE scRNAseq", "DE both"))
  
  df <- df[order(df$DE),]
  p <- ggplot(df, aes(x = FC_scRNA, y = FC_Liveseq, label = gene, col = DE)) + 
    geom_point(size = 0.8) + 
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
  
  pearson_ALL <- cor(df$FC_Liveseq, df$FC_scRNA)
  d <- filter(df, DE %in% c("DE both", "DE scRNAseq"))
  pearson_DEsc <- cor(d$FC_Liveseq, d$FC_scRNA)
  d <- filter(df, DE %in% c("DE both"))
  pearson_DEboth <- cor(d$FC_Liveseq, d$FC_scRNA)
  
  return(list(p = p,
              lm_AllGenes=lm_ALL,
              lm_DEsc=lm_DEsc,
              lm_DEboth=lm_DEboth,
              pearson_AllGenes = pearson_ALL,
              pearson_DEsc = pearson_DEsc,
              pearson_DEboth = pearson_DEboth))
  
}
