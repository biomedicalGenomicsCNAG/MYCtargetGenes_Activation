## Single-Cell RNA-Seq Data Analysis Pipeline

This repository contains an analysis pipeline for single-cell RNA-Seq data, focusing on the study of MYC activation in cancer samples.

### Key Steps:
1. **Load Necessary Libraries:**
   - Load various libraries including `data.table`, `Seurat`, `ggpubr`, `tidyverse`, `pacman`, `UCell`, and `liana`.

2. **Load Main Object:**
   - Load the main Seurat object from [Zenodo](https://zenodo.org/records/6631966).

3. **Update Seurat Object and Preprocess:**
   - Normalize data, find variable features, scale data, run PCA, and UMAP.
   - Find neighbors and clusters with optimal resolution.

4. **Plot UMAP and Annotate Clusters:**
   - Generate UMAP plots and manually annotate clusters using known marker genes.

5. **Cell Cycle Scoring:**
   - Score cell cycle stages using predefined gene sets.

6. **Add MYC Activation Signature:**
   - Load and add MYC activation signature to the Seurat object.

7. **Generate Supplementary Figures:**
   - Create violin plots and boxplots to visualize MYC activation.

8. **CLL Phase Analysis:**
   - Subset data for CLL phase, preprocess, and annotate clusters.
   - Categorize MYC activation levels and generate relevant plots.

9. **RS Phase Analysis:**
   - Subset data for RS phase, preprocess, and annotate clusters.
   - Categorize MYC activation levels and generate relevant plots.

10. **LIANA - Cell-Cell Interactions:**
    - Analyze cell-cell interactions using LIANA and generate heatmaps and dot plots.

### Dependencies:
- data.table, Seurat, ggpubr, tidyverse, pacman, UCell, liana

### Data Sources:
- The main Seurat object can be downloaded from [Zenodo](https://zenodo.org/records/6631966).
- Cell cycle scoring gene sets can be downloaded from [Dropbox](https://www.dropbox.com/s/3dby3bjsaf5arrw/cell_cycle_vignette_files.zip?dl=1).

### Description
The script and the required files to analyze MYC in cancer samples using single-cell RNAseq data.
