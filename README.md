# Live-seq

## Abstract
Live-seq is a single-cell transcriptome profiling approach that preserves cell viability during RNA extraction, using fluidic force microscopy. Therefore, Live-seq can address a broad range of biological questions by transforming scRNA-seq from an end-point to a temporal analysis approach.

## Citation
W. Chen, O. Guillaume-Gentil, et al., [Genome-wide molecular recording using Live-seq](https://www.biorxiv.org/content/10.1101/2021.03.24.436752v1), *in press*, 2022

## Pipeline
This repository contains the pipeline to be able to reproduce all figures of the paper.

**00. Dependencies:** To install the dependencies for the project. You can find the package versions used in the paper in [sessionInfo.txt](sessionInfo.txt)

**01. Preprocessing:** Create a Seurat object from the count matrix generated after Live-seq and single-cell sequencing. Some basic single-cell QCs.

**02. Live-seq:** Subset the Seurat object to Live-seq cells only. Generate Live-seq-specific QCs. Downsampling of Live-seq.

**03. scNRA-seq:** Subset the Seurat object to single-cell RNA-seq (scRNA-seq) cells only. Generate scRNA-seq-specific QCs. Downsampling of scRNA-seq.

**04. Live-seq scNRA-seq integration:** Integration of both layers of data using [Seurat's integration pipeline](https://satijalab.org/seurat/articles/integration_introduction.html).

**05. Live-seq with live cell imaging:** Connecting the transcriptomic profile with a downstream phenotypic response, i.e. TNF upregulation upon LPS treatment

**06. Clustering:** Compute clustering accuracy (ARI, Barplot comparison and clustree)

**07. Differential expression across cell types:** Compute DE for each cluster (corresponding to a cell type/state) versus the rest for both Live-seq and scRNA-seq data. EnrichR on DE genes for BP GO terms and Mouse Cell Atlas, of each cluster versus Rest per sampling method

**08. Analysis per cell type:** Compute figures Live-seq manuscript per cell types: 1. tSNE per cell types colored by metadata and clustering, 2. tSNE colored per extracted volumes, 3. Avg expression per cell type - correlation between sc and live

**09. Differential expression within cell type:** Compare DE results obtained with scRNA and Live-seq. Identify GO Terms of genes detected only by live-seq or scRNA-seq to find any potential bias.

**10. Downsampling scRNA-seq:** Downsample scRNAseq data so that match complexity of Live-seq -> perform DE. Compare DE results obtained with scRNA and Live-seq. Downsampling scripts for scRNA-seq. (see **2.** for downsampling of Live-seq data). Compare number of common genes between Live-seq vs scRNA-seq OR Live-seq vs Downsampled scRNA-seq.

**./utils** Some utility functions used across the scripts

## Data
This repository also contains the data used in the pipeline.

The main count matrix is downloadable on GEO: [GSE141064](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141064).
In the pipeline, it's automatically downloaded and parsed in the first script: **01. Preprocessing**

**[./data](data):** 
  - *mouseGeneTable87_mCherry_EGFP.txt* and *Mus_musculus.GRCm38.100_data.annot.txt*: The gene annotation information (gene name, biotype, exon length, ...)
  - *gene.blacklist.csv*: List of 20 genes that are blacklisted and removed from many analyses, which are derived from the 0 pg input RNA negative control. 
  - *InfoContent_Updated_VX-ASPC-9_4.csv* and *meta.final.csv*: Sample (cell) metadata
  - *log.foreach.txt*: Log of clustering script
