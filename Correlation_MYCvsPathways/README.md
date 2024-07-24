## Correlation of MYC target gene activation score with other pathways

This repository contains a script for correlation analysis using bulk RNA-seq data from ICGC and RS datasets.

### Key Steps:
1. **Load Necessary Libraries:**
   - Load libraries including `Seurat`, `ggpubr`, `tidyverse`, `GSVA`, and `ggplot2`.

2. **ICGC Analysis:**
   - Generate the `df` table using the script `Bulk_RNAseq_ICGC.R`.
   - Load hallmark gene sets using `msigdbr`.
   - Calculate signatures and correlations for each pathway using GSVA.
   - Read and calculate correlation with KEGG_BCR pathway.
   - Update the results table with the correlations and p-values.

3. **RS Analysis:**
   - Generate the `data_vvt` object using the script `Bulk_RNAseq_RT_MYC.R`.
   - Load hallmark gene sets using `msigdbr`.
   - Calculate signatures and correlations for each pathway using GSVA.
   - Read and calculate correlation with KEGG_BCR pathway.
   - Update the results table with the correlations and p-values.

4. **Combine and Visualize Results:**
   - Combine results from ICGC and RS analyses.
   - Create a dot plot with size representing the correlation and color representing the p-value.
   - Adjust p-values for multiple testing and filter significant pathways.

5. **Generate Final Plots:**
   - Create a filtered dot plot highlighting significant correlations.

### Prerequisites:
- The `df` table must be generated using the script `Bulk_RNAseq_ICGC.R`. The script and related files are available [here](https://github.com/biomedicalGenomicsCNAG/MYCtargetGenes_Activation/tree/main/Bulk_RNAseq_ICGC).
- The `data_vvt` object must be produced using the script `Bulk_RNAseq_RT_MYC.R`. The script and related files are available [here](https://github.com/biomedicalGenomicsCNAG/MYCtargetGenes_Activation/tree/main/Bulk_RNAseq_RT).

### Required Files:
- All required files are located within this folder.

### Data Sources:
- Hallmark gene sets are loaded using the `msigdbr` package.

### Description
This pipeline provides tools for the analysis of bulk RNA-seq data using GSVA, from preprocessing to advanced statistical analyses and visualization.

To run this script, ensure the following prerequisites are met:
- The `df` table must be generated using the script `Bulk_RNAseq_ICGC.R`. The script and related files are available [here](https://github.com/biomedicalGenomicsCNAG/MYCtargetGenes_Activation/tree/main/Bulk_RNAseq_ICGC).
- The `data_vvt` object must be produced using the script `Bulk_RNAseq_RT_MYC.R`. The script and related files are available [here](https://github.com/biomedicalGenomicsCNAG/MYCtargetGenes_Activation/tree/main/Bulk_RNAseq_RT).
