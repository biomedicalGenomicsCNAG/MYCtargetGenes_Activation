# MYC Target Gene Activation in Chronic Lymphocytic Leukemia and Richter Transformation: Links to Aggressiveness and Tumor Microenvironment Interactions

This repository contains code and links to the used data accompanying the study.

## Abstract

This study investigates the activation of MYC target genes in Chronic Lymphocytic Leukemia (CLL) and Richter Transformation (RT) to understand their links to disease aggressiveness and interactions within the tumor microenvironment.

## MYC Target Gene Activation Signature

The genes can be found [here](https://github.com/biomedicalGenomicsCNAG/MYCtargetGenes_Activation/tree/main/MYCtargetgenes_signature).

## Analyses

- **Bulk RNA-seq Analysis**
  - Code to reproduce the ICGC analysis can be found [here](https://github.com/biomedicalGenomicsCNAG/MYCtargetGenes_Activation/tree/main/Bulk_RNAseq_ICGC).
  - Code to reproduce the RT analysis can be found [here](https://github.com/biomedicalGenomicsCNAG/MYCtargetGenes_Activation/tree/main/Bulk_RNAseq_RT).

- **Single-Cell RNA-seq Analysis**
  - Analysis of single-cell RNA-seq data from Nadeu et al. can be found [here](https://github.com/biomedicalGenomicsCNAG/MYCtargetGenes_Activation/tree/main/SingleCell_RNAseq_Nadeu_RT).
  - Analysis of single-cell RNA-seq data under ibrutinib treatment can be found [here](https://github.com/biomedicalGenomicsCNAG/MYCtargetGenes_Activation/tree/main/SingleCell_RNAseq_ibrutinib_GSE111015).

## Data Availability

- The CLL phase files were downloaded from the CLL-map project and are available at [CLLmap](https://cllmap.org/downloads.html) under the gene counts for the ICGC project.
- The bulk RNA-seq data from RT samples in the RS phase are available through a kallisto table, which can be downloaded from [GitHub](https://github.com/ferrannadeu/RichterTransformation/tree/main/bulkRNA-seq/kallisto).
- The scRNA-seq expression object from the Nadeu et al. study, along with their metadata, is available on [Zenodo](https://zenodo.org/records/6631966).
- The single-cell RNA-seq data under ibrutinib treatment were downloaded from GEO [GSE111015](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111015) as reported by Rendeiro et al.


## About

This repository provides comprehensive tools for the analysis of MYC target gene activation in CLL and RT using both bulk and single-cell RNA-seq data.

### License

This project is licensed under the MIT License.
