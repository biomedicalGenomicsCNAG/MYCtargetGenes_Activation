
RNA-Seq Data Analysis Pipeline
This repository contains an analysis pipeline for RNA-Seq data, focusing on the study of MYC activation in cancer samples.

Key Steps:
Data Loading and Preprocessing:

Download RNA-Seq count data and metadata.
Preprocess data by removing duplicates and normalizing with DESeq2.
Batch Effect Correction:

Correct for batch effects using ComBat.
Scoring and Gene Set Variation Analysis:

Score samples based on MYC target gene activation using GSVA.
Survival and Differential Expression Analysis:

Perform survival analysis to correlate MYC activation with patient outcomes.
Conduct differential expression analysis to identify significant genes.
Visualization:

Generate various plots to visualize MYC activation levels, survival curves, and differential expression results.
GSEA and Mutation Analysis:

Perform Gene Set Enrichment Analysis (GSEA).
Analyze the impact of specific gene mutations on MYC activation.
Visualizations:
Boxplots and scatter plots to visualize gene expression and activation scores.
Heatmaps for clustering and visualization of differentially expressed genes.
Survival curves illustrating the impact of MYC activation on patient outcomes.
Dependencies:
DESeq2, dplyr, tidyverse, data.table, ggplot2, GSVA, sva, survminer, survival, ComplexHeatmap, limma, fgsea, msigdbr, clusterProfiler, enrichplot
Data Sources:
RNA-Seq counts and metadata can be downloaded from CLLmap.
