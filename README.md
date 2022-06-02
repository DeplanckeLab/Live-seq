# Live-seq

## Abstract
Live-seq is a single-cell transcriptome profiling approach that preserves cell viability during RNA extraction, using fluidic force microscopy. Therefore, Live-seq can address a broad range of biological questions by transforming scRNA-seq from an end-point to a temporal analysis approach.

## Citation
W. Chen, O. Guillaume-Gentil, et al., [Genome-wide molecular recording using Live-seq](https://www.biorxiv.org/content/10.1101/2021.03.24.436752v1), in press, 2022

## Pipeline
This repository contains the pipeline to be able to reproduce all figures of the paper.

**1. Preprocessing:** Create a Seurat object from the count matrix generated after Live-seq and single-cell sequencing. Some basic single-cell QCs.
**2. Live-seq:** Subset the Seurat object to Live-seq cells only. Generate Live-seq-specific QCs. Downsampling of Live-seq.
**3. scNRA-seq:** Subset the Seurat object to single-cell RNA-seq (scRNA-seq) cells only. Generate scRNA-seq-specific QCs. Downsampling of scRNA-seq.
**4. Live-seq scNRA-seq integration:** Integration of both layers of data using [Seurat's integration pipeline](https://satijalab.org/seurat/articles/integration_introduction.html).
**5. Live-seq with LiveCell imaging:** ??
**6. Clustering:**
**7. Differential expression across cell types:**
**8. Analysis per cell type:**
**9. Differential expression within cell type:**
**10. Downsampling scRNA-seq:** Downsampling scripts for scRNA-seq. (see **2.** for downsampling of Live-seq data)

## Data
This repository also contains the data used in the pipeline.

**- data:** 
  - mouseGeneTable87_mCherry_EGFP.txt and Mus_musculus.GRCm38.100_data.annot.txt: The gene annotation information (gene name, biotype, exon length, ...)
  - gene.blacklist.csv: List of 20 genes that are blacklisted and removed from many analyses, which are derived from the 0 pg input RNA negative control. 
  - mouseGeneTable87_mCherry_EGFP.txt: The count matrix obtained after preprocessing the sequencing data, i.e. alignment and gene annotation of the reads
  - 
