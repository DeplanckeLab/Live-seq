root_dir <- find_root(has_file("Live-seq.RProj"))

require(ggplot2)

#aesthetics 
celltype_treatment_colors = c(
  "ASPC_not_treated"="#8a804aff",
  "ASPC_DMIR_treated"="#dfc17aff",
  "IBA_not_treated"="#72b3e3ff",
  "Raw264.7_not_treated"="#75bfb3ff",
  "Raw264.7_LPS_treated"="#0c8772ff"
)
sampling_type_colors <- c("Live_seq"="#2289caff", "scRNA"="#cd6228ff")

mashaGgplot2Theme <- list(
  theme_classic(base_size = 14) + 
    theme(text = element_text(size = 14)) +
    theme(axis.line.x = element_line(colour = 'black', size = 0.5,
                                     linetype = 'solid'),
          axis.line.y = element_line(colour = 'black', size=0.5,
                                     linetype ='solid'),
          panel.grid.minor = element_line(colour = "white", size = 0.5,
                                          linetype = 2))
)

data.annot <- read.table(paste0(root_dir, "/data/Mus_musculus.GRCm38.100_data.annot.txt"), header = T)
rownames(data.annot) <- data.annot$Ensembl
gene.blacklist <- read.csv(paste0(root_dir, "/data/gene.blacklist.csv"))

