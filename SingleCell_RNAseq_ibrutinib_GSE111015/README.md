## Single-Cell RNA-Seq Data Analysis Pipeline

This repository contains an analysis pipeline for single-cell RNA-Seq data, focusing on the study of MYC activation in cancer samples.

### Key Steps:
1. **Load Necessary Libraries:**
   - Load various libraries including `Seurat`, `ggpubr`, `tidyverse`, `data.table`, `rlist`, `plyr`, `pacman`, `UCell`, and `liana`.

2. **Set Working Directory and Load Data:**
   - Set the working directory and load data from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111015).
   - Load the data and create a Seurat object.

3. **Process Metadata:**
   - Process and join metadata to the Seurat object.

4. **Save and Load Seurat Object:**
   - Save the Seurat object and reload it, setting the default assay to RNA.
   - Update the Seurat object.

5. **Analyze Metadata:**
   - Analyze and summarize metadata, generating a violin plot of library size per tissue.

6. **Quality Control:**
   - Add mitochondrial percentage to the Seurat object and generate histograms for quality control.
   - Subset low-quality cells and filter out low-quality genes.

7. **Preprocess Data:**
   - Normalize data, find variable features, scale data, run PCA, and UMAP.
   - Find neighbors and clusters with optimal resolution.

8. **Generate Feature Plots and Annotate Clusters:**
   - Generate feature plots and annotate clusters manually.

9. **Add MYC Activation Signature:**
   - Load and add MYC activation signature to the Seurat object.

### Time Point 1 Analysis (Before Ibrutinib):
1. **Subset Data for Time Point 1:**
   - Subset the data for the first time point (T1) and preprocess.
   - Annotate MYC activation categories.

2. **Generate Plots:**
   - Generate UMAP plots colored by MYC activation category and patient ID.

### LIANA Analysis:
1. **Perform LIANA Analysis:**
   - Analyze cell-cell interactions using LIANA and generate a heatmap.

### MYC Under Ibrutinib:
1. **Subset Data for B Cells:**
   - Subset the data for B cells and add metadata.
   - Generate dot plots to visualize MYC activation under ibrutinib treatment.

### Dependencies:
- Seurat, ggpubr, tidyverse, data.table, rlist, plyr, pacman, UCell, liana, dplyr

### Data Sources:
- Data can be downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111015).
- Metadata file should be provided in the working directory.

### Description
The script and the required files to analyze MYC in cancer samples using single-cell RNAseq data.

