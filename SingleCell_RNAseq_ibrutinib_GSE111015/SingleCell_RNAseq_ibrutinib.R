# Load necessary libraries
library(Seurat)
library(ggpubr)
library(tidyverse)
library(data.table)
library(rlist)
library(plyr)
library(pacman)
library(UCell)
library(liana)
library(dplyr)

# Set working directory and load data
data_dir <- "CLL_GSE111015_ibr" #add your folder, download the file from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111015
setwd(data_dir)
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
CLL_ibr <- Read10X(data.dir = data_dir)
CLL_ibr <- CreateSeuratObject(counts = CLL_ibr)

# Process metadata
meta1 <- names(Idents(CLL_ibr))
meta1 <- str_split_fixed(meta1, "_", 2)
meta1 <- as.data.frame(meta1)
names(meta1) <- c("V1", "sample_name")
annotation <- read.csv("metadata.csv")
metadata <- join(meta1, annotation, by = "sample_name")

# Add metadata to Seurat object
CLL_ibr <- AddMetaData(CLL_ibr, metadata = metadata$sample_name, col.name = "Sample_id")
CLL_ibr <- AddMetaData(CLL_ibr, metadata = metadata$project_id, col.name = "Project_id")
CLL_ibr <- AddMetaData(CLL_ibr, metadata = metadata$patient_id, col.name = "Patient_id")
CLL_ibr <- AddMetaData(CLL_ibr, metadata = metadata$Diagnosis, col.name = "Diagnosis")
CLL_ibr <- AddMetaData(CLL_ibr, metadata = metadata$timepoint, col.name = "Timepoint")
CLL_ibr <- AddMetaData(CLL_ibr, metadata = metadata$Status, col.name = "Status")
CLL_ibr <- AddMetaData(CLL_ibr, metadata = metadata$sex, col.name = "Sex")
CLL_ibr <- AddMetaData(CLL_ibr, metadata = metadata$ighv_mutation, col.name = "Category")
CLL_ibr <- AddMetaData(CLL_ibr, metadata = metadata$number_of_prior_treatments, col.name = "number_of_prior_treatments")
CLL_ibr <- AddMetaData(CLL_ibr, metadata = metadata$response, col.name = "response")
CLL_ibr <- AddMetaData(CLL_ibr, metadata = metadata$Tissue, col.name = "Tissue")
CLL_ibr <- AddMetaData(CLL_ibr, metadata = metadata$Genetics_ibrutinib, col.name = "Genetics")

# Save Seurat object
saveRDS(CLL_ibr, file = paste0("CLL_ibr_2023.rds"))

# Load Seurat object and set default assay
seurat <- readRDS('CLL_ibr_2023.rds')
DefaultAssay(seurat) <- "RNA"

# Update Seurat object
p_unload(SeuratDisk)
p_unload(Seurat)
p_load(Seurat) # Using version 3.2.0
seurat <- UpdateSeuratObject(seurat)

# Analyze metadata
meta <- seurat@meta.data
unique(meta$Patient_id)

lib_complexity_df <- seurat@meta.data %>%
  group_by(Project_id, as.character(Patient_id)) %>%
  summarise(mean_n_features = round(mean(nFeature_RNA), 2), sd_n_features = round(sd(nFeature_RNA), 2)) %>%
  arrange(Project_id)
print(lib_complexity_df)

# Generate violin plot of library size per tissue
lib_size_per_tissue <- seurat@meta.data %>%
  ggplot(aes(Project_id, nFeature_RNA, fill = Project_id)) +
  geom_violin() +
  scale_y_log10() +
  labs(title = "", x = "", y = "Number of Detected Genes") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.title = element_blank()
  )
lib_size_per_tissue

# Add mitochondrial percentage to Seurat object
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

# Choose cut-offs for filtering
min_lib_size <- 1000
max_lib_size <- 12000
min_n_genes <- 100
max_n_genes <- 3000
max_pct_mt <- 15
min_cells <- 5

# Generate histograms for quality control
lib_size_hist <- seurat@meta.data %>%
  ggplot(aes(x = nCount_RNA)) +
  geom_histogram() +
  geom_vline(xintercept = min_lib_size, linetype = "dashed", color = "red") +
  geom_vline(xintercept = max_lib_size, linetype = "dashed", color = "red") +
  scale_x_continuous(limits = c(0, 30000))
lib_size_hist

n_genes_hist <- seurat@meta.data %>%
  ggplot(aes(x = nFeature_RNA)) +
  geom_histogram() +
  geom_vline(xintercept = min_n_genes, linetype = "dashed", color = "red") +
  geom_vline(xintercept = max_n_genes, linetype = "dashed", color = "red")
n_genes_hist

pct_mt_hist <- seurat@meta.data %>%
  ggplot(aes(x = percent.mt)) +
  geom_histogram() +
  geom_vline(xintercept = max_pct_mt, linetype = "dashed", color = "red") +
  scale_x_continuous(limits = c(0, 100))
pct_mt_hist

# Choose cut-offs for filtering
min_lib_size <- 1000
max_lib_size <- 12000
min_n_genes <- 300
max_n_genes <- 3000
max_pct_mt <- 15
min_cells <- 5

# Subset low quality cells
is_low_quality <- seurat$nCount_RNA < min_lib_size | seurat$nCount_RNA > max_lib_size | seurat$nFeature_RNA < min_n_genes | seurat$nFeature_RNA > max_n_genes
table(is_low_quality)
seurat$keep_cells <- !is_low_quality
Idents(seurat) <- "keep_cells"
seurat <- subset(seurat, idents = TRUE)

# Filter out low quality genes
n_cells <- Matrix::rowSums(seurat[["RNA"]]@counts > 0)
kept_genes <- rownames(seurat)[n_cells > min_cells]
table(n_cells > min_cells)

seurat <- subset(seurat, features = kept_genes)
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:30)
optimal_resolution <- 0.06
seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:30) %>% FindClusters(resolution = optimal_resolution)

# Generate feature plots
FeaturePlot(seurat, features = c("CD79A", "CD3D", "IL7R", "NKG7", "LYZ"), pt.size = 0.45)

# Annotate clusters
seurat$annotation_level1 <- seurat$seurat_clusters
annotation_level1 <- c("B cells", "T cells", "B cells", "Myeloid cells", "B cells", "NK cells", "Plasma cells")
levels(seurat$annotation_level1) <- annotation_level1
Idents(seurat) <- "annotation_level1"

# Add MYC activation signature
signature <- fread('MYCtargetgenes.txt')
colnames(signature) = "MYC_activation"
signature = as.data.frame(signature)
markers <- list(MYC_activation = signature)
seurat <- AddModuleScore_UCell(obj = seurat, features = markers)

################################################################################
# Time point 1 analysis, before ibrutinib

# Subset data for time point 1
Idents(seurat) = seurat$Timepoint
testdata <- subset(seurat, idents = c("T1"))
testdata <- NormalizeData(testdata)
testdata <- FindVariableFeatures(testdata, selection.method = "vst", nfeatures = 2000)
testdata <- ScaleData(testdata)
testdata <- RunPCA(testdata)
testdata <- RunUMAP(testdata, reduction = "pca", dims = 1:30)
optimal_resolution <- 0.06
testdata <- FindNeighbors(testdata, reduction = "pca", dims = 1:30) %>% FindClusters(resolution = optimal_resolution)

# Annotate MYC activation category
Idents(testdata) = testdata$annotation_level1
testdata$MYC_activation_category <- case_when(
  testdata$annotation_level1 == "B cells" ~ case_when(
    testdata$MYC_activation_UCell <= quantile(testdata$MYC_activation_UCell, 0.25) ~ "low MYC_TGA",
    testdata$MYC_activation_UCell <= quantile(testdata$MYC_activation_UCell, 0.75) ~ "intermediate MYC_TGA",
    TRUE ~ "high MYC_TGA"
  ),
  TRUE ~ as.character(testdata$annotation_level1)
)

# Generate Supplementary Figure 4a
cols <- c("#F8766D", "#00BA38", "#5254A3", "#FD8D3C", "#619CFF", "pink", "purple")
names(cols) = levels(Idents(testdata))
cols2 <- c("gold3", "brown", "gray22", "pink4")
DimPlot(testdata, group.by = c("MYC_activation_category"), cols = cols)
DimPlot(testdata, group.by = c("Patient_id"), cols = cols2)

################################################################################
# LIANA Analysis

unique(testdata$MYC_activation_category)
Idents(testdata) = testdata$MYC_activation_category
liana_test <- liana_wrap(testdata)
liana_test <- liana_test %>% liana_aggregate()

# Filter LIANA results for heatmap
filtered_liana_heatmap <- liana_test %>%
  filter(
    !(source == "Myeloid cells" & target == "Myeloid cells") &
      !(source == "Myeloid cells" & target == "T cells") &
      !(target == "Myeloid cells" & source == "T cells") &
      !(source == "T cells" & target == "Myeloid cells") &
      !(source == "Myeloid cells" & target == "NK cells") &
      !(target == "Myeloid cells" & source == "NK cells") &
      !(source == "T cells" & target == "NK cells") &
      !(target == "T cells" & source == "NK cells") &
      !(source == "NK cells" & target == "NK cells") &
      !(source == "T cells" & target == "T cells")
  ) %>%
  ungroup()

# Generate Supplementary Figure 4b
heat_freq(filtered_liana_heatmap, pallette = c("blue", "white", "violetred2"))


################################################################################
# MYC under Ibrutinib
Idents(seurat) = seurat$annotation_level1
seurat_B <- subset(seurat, idents = c("B cells"))
seurat_B <- AddMetaData(object = seurat_B, metadata = paste(seurat_B@meta.data$Patient_id, seurat_B@meta.data$Timepoint,  sep = "_"), col.name = "donor")

# Supplementary figure 4c
DotPlot(seurat_B, features= "MYC_activation_UCell", group.by ="donor")
