## RNA-Seq Data Analysis Pipeline

This repository contains an analysis pipeline for RNA-Seq data, focusing on the study of MYC activation in cancer samples.

### Key Steps:
1. **Load Necessary Libraries:**
   - Load various libraries including `tximport`, `DESeq2`, `ggplot2`, `ComplexHeatmap`, and others.

2. **Load Gene and Transcript Data:**
   - Load gene and transcript data from Ensembl.

3. **Load Metadata:**
   - Load metadata and prepare file names for transcript data.

4. **Import Transcript Data Using tximport:**
   - Import transcript data using `tximport`.

5. **Prepare DESeq2 Dataset:**
   - Prepare and normalize the DESeq2 dataset.

6. **Annotate Data with Gene Symbols:**
   - Annotate the dataset with gene symbols and remove duplicates.

7. **Load Gene Signature:**
   - Load the MYC target gene signature.

8. **GSVA Analysis:**
   - Perform Gene Set Variation Analysis (GSVA) to score samples based on MYC target gene activation.

9. **Heatmap:**
   - Prepare data for heatmap visualization and generate the heatmap.

10. **Boxplot:**
    - Generate a boxplot to visualize MYC activation scores.

11. **Correlation with MYC Expression:**
    - Perform correlation analysis between MYC expression and activation scores.

### Dependencies:
- tximport, DESeq2, ggplot2, ComplexHeatmap, circlize, ggpubr, data.table, dplyr, GSVA

### Data Sources:
- Gene and transcript data from Ensembl.
- Metadata and transcript data from RNAseq annotation files.
- The original dataset can be downloaded from [GitHub](https://github.com/ferrannadeu/RichterTransformation/tree/main/bulkRNA-seq/kallisto).

### Description
The script and the required files to analyze MYC in RT using bulk RNAseq data.
