# Load necessary libraries
library(data.table)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(pacman)
library(UCell)
library(liana)

################################################################################
# Load main object, https://zenodo.org/records/6631966 
main_obj <- readRDS('3.seurat_filtered.rds')
unique(main_obj$donor_id)

# Update Seurat object and preprocess
p_unload(SeuratDisk)
p_unload(Seurat)
p_load(Seurat)
main_obj <- UpdateSeuratObject(main_obj)
main_obj <- NormalizeData(main_obj)
main_obj <- FindVariableFeatures(main_obj, selection.method = "vst", nfeatures = 2000)
main_obj <- ScaleData(main_obj)
main_obj <- RunPCA(main_obj)
main_obj <- RunUMAP(main_obj, reduction = "pca", dims = 1:20)
optimal_resolution <- 0.1
main_obj <- FindNeighbors(main_obj, dims = 1:20) %>% FindClusters(resolution = optimal_resolution)

# Plot UMAP
DimPlot(main_obj, label = T)

# Manual annotation
features <- c("LYZ", "CD14", "FCGR3A", "CD3D", "NKG7", "CD79A", "MZB1", "CD19", "CD5")
DotPlot(main_obj, features = features) + RotatedAxis()

# Name the clusters
main_obj$annotation_level1 <- main_obj$seurat_clusters
annotation_level1 <- c("B cells", "B cells", "B cells", "T cells", "B cells", "B cells", "B cells", "B cells", "Myeloid cells")
levels(main_obj$annotation_level1) <- annotation_level1
Idents(main_obj) <- "annotation_level1"

# Cell cycle scoring
#downloaded from https://www.dropbox.com/s/3dby3bjsaf5arrw/cell_cycle_vignette_files.zip?dl=1
exp.mat <- read.table(file = "cell_cycle_vignette_files\\nestorawa_forcellcycle_expressionMatrix.txt", header = TRUE, as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
main_obj <- CellCycleScoring(main_obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Add MYC activation signature
signature <- fread('MYCtargetgenes.txt')
colnames(signature) = "MYC_activation"
signature = as.data.frame(signature)
markers <- list(MYC_activation = signature)
main_obj <- AddModuleScore_UCell(obj = main_obj, features = markers)

# Generate Supplementary Figure 3e
Idents(main_obj) = main_obj$annotation_level1
sub <- subset(main_obj, idents = c("T cells", "Myeloid cells"), invert = TRUE)
group_means <- tapply(sub$MYC_activation_UCell, sub$sample_description_FN, mean)
print(group_means)
reference_group <- names(sort(group_means, decreasing = TRUE)[1])
VlnPlot(sub, features = "MYC_activation_UCell", group.by = c("sample_description_FN"), cols = c("gray92", "gray62", "gray52", "gray42", "gray32", "gray12"), sort = TRUE, pt.size = 0, y.max = max(main_obj$MYC_activation_UCell) + 0.02) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_compare_means(method = "anova", label.y = 40) + 
  stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test", ref.group = reference_group)

################################################################################
# CLL Phase Analysis

# Subset for CLL phase and preprocess
Idents(main_obj) = main_obj$sample_description_FN
main_obj_CLL <- subset(main_obj, idents = c("RS"), invert = TRUE)
main_obj_CLL <- NormalizeData(main_obj_CLL)
main_obj_CLL <- FindVariableFeatures(main_obj_CLL, selection.method = "vst", nfeatures = 2000)
main_obj_CLL <- ScaleData(main_obj_CLL)
main_obj_CLL <- RunPCA(main_obj_CLL)
main_obj_CLL <- RunUMAP(main_obj_CLL, reduction = "pca", dims = 1:30)
optimal_resolution <- 0.06
main_obj_CLL <- FindNeighbors(main_obj_CLL, reduction = "pca", dims = 1:30) %>% FindClusters(resolution = optimal_resolution)

# Annotate and categorize MYC activation
Idents(main_obj_CLL) = main_obj_CLL$annotation_level1
main_obj_CLL$MYC_activation_category <- case_when(
  main_obj_CLL$annotation_level1 == "B cells" ~ case_when(
    main_obj_CLL$MYC_activation_UCell <= quantile(main_obj_CLL$MYC_activation_UCell, 0.25) ~ "low MYC_TGA",
    main_obj_CLL$MYC_activation_UCell <= quantile(main_obj_CLL$MYC_activation_UCell, 0.75) ~ "intermediate MYC_TGA",
    TRUE ~ "high MYC_TGA"
  ),
  TRUE ~ as.character(main_obj_CLL$annotation_level1)
)

# Generate Figure 2a
cols <- c("#F8766D", "#00BA38", "#619CFF", "pink", "purple")
DimPlot(main_obj_CLL, group.by = "MYC_activation_category", cols = cols, raster = FALSE)

# Generate DotPlot for Supplementary Figure 3a
DotPlot(main_obj_CLL, feature = c("CD79A", "CD19", "CD3D", "CD68", "CD14", "FCGR3A")) & scale_color_viridis_c(option = "C") & RotatedAxis()

# Further subset and annotate for cell cycle analysis
Idents(main_obj_CLL) = main_obj_CLL$annotation_level1
main_obj_CLL_B <- subset(main_obj_CLL, idents = c("B cells"))
main_obj_CLL_B <- AddMetaData(object = main_obj_CLL_B, metadata = paste(main_obj_CLL_B@meta.data$MYC_activation_category, main_obj_CLL_B@meta.data$donor_id, sep = "_"), col.name = "MYC_donor")
main_obj_CLL_B2 <- CellCycleScoring(main_obj_CLL_B, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Generate Supplementary Figure 3c
FeatureScatter(main_obj_CLL_B2, feature1 = "S.Score", feature2 = "MYC_activation_UCell", group.by = "MYC_activation_category")

################################################################################
# RS Phase Analysis

# Subset for RS phase and preprocess
Idents(main_obj) = main_obj$sample_description_FN
main_obj_RS <- subset(main_obj, idents = c("RS"))
main_obj_RS <- NormalizeData(main_obj_RS)
main_obj_RS <- FindVariableFeatures(main_obj_RS, selection.method = "vst", nfeatures = 2000)
main_obj_RS <- ScaleData(main_obj_RS)
main_obj_RS <- RunPCA(main_obj_RS)
main_obj_RS <- RunUMAP(main_obj_RS, reduction = "pca", dims = 1:30)
optimal_resolution <- 0.06
main_obj_RS <- FindNeighbors(main_obj_RS, reduction = "pca", dims = 1:30) %>% FindClusters(resolution = optimal_resolution)

# Annotate and categorize MYC activation
Idents(main_obj_RS) = main_obj_RS$annotation_level1
main_obj_RS$MYC_activation_category <- case_when(
  main_obj_RS$annotation_level1 == "B cells" ~ case_when(
    main_obj_RS$MYC_activation_UCell <= quantile(main_obj_RS$MYC_activation_UCell, 0.25) ~ "low MYC_TGA",
    main_obj_RS$MYC_activation_UCell <= quantile(main_obj_RS$MYC_activation_UCell, 0.75) ~ "intermediate MYC_TGA",
    TRUE ~ "high MYC_TGA"
  ),
  TRUE ~ as.character(main_obj_RS$annotation_level1)
)

# Generate Figure 2b
cols <- c("#F8766D", "#00BA38", "#619CFF", "pink", "purple")
DimPlot(main_obj_RS, group.by = "MYC_activation_category", cols = cols, raster = FALSE)

# Generate DotPlot for Supplementary Figure 3b
DotPlot(main_obj_RS, feature = c("CD79A", "CD19", "CD3D", "CD68", "CD14", "FCGR3A")) & scale_color_viridis_c(option = "C") & RotatedAxis()

# Further subset and annotate for cell cycle analysis
Idents(main_obj_RS) = main_obj_RS$annotation_level1
main_obj_RS_B <- subset(main_obj_RS, idents = c("B cells"))
main_obj_RS_B <- AddMetaData(object = main_obj_RS_B, metadata = paste(main_obj_RS_B@meta.data$MYC_activation_category, main_obj_RS_B@meta.data$donor_id, sep = "_"), col.name = "MYC_donor")
main_obj_RS_B2 <- CellCycleScoring(main_obj_RS_B, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Generate Supplementary Figure 3d
FeatureScatter(main_obj_RS_B2, feature1 = "S.Score", feature2 = "MYC_activation_UCell", group.by = "MYC_activation_category")

################################################################################
# LIANA - Cell-Cell Interactions

# CLL Phase
unique(main_obj_CLL$MYC_activation_category)
Idents(main_obj_CLL) = main_obj_CLL$MYC_activation_category
liana_test <- liana_wrap(main_obj_CLL)
liana_test <- liana_test %>% liana_aggregate()

# Filter LIANA results
filtered_liana_test <- liana_test %>%
  group_by(source, ligand.complex, receptor.complex) %>%
  mutate(std_dev_natmi = sd(natmi.edge_specificity, na.rm = TRUE), count = n()) %>%
  filter(count == 1 | std_dev_natmi > 0.1) %>%
  ungroup()

filtered_liana_test_target <- liana_test %>%
  group_by(target, ligand.complex, receptor.complex) %>%
  mutate(std_dev_natmi = sd(natmi.edge_specificity, na.rm = TRUE), count = n()) %>%
  filter(count == 1 | std_dev_natmi > 0.1) %>%
  ungroup()

# Generate heatmap for Figure 2c
filtered_liana_heatmap <- liana_test %>%
  filter(!(source == "Myeloid cells" & target == "Myeloid cells") &
           !(target == "Myeloid cells" & source == "Myeloid cells") &
           !(source == "Myeloid cells" & target == "T cells") &
           !(target == "Myeloid cells" & source == "T cells") &
           !(source == "T cells" & target == "Myeloid cells") &
           !(target == "T cells" & source == "Myeloid cells")) %>%
  ungroup()
heat_freq(filtered_liana_heatmap, pallette = c("blue", "white", "violetred2"))

# Generate dot plots for Figures 2e and Supplementary Figure 3f
filtered_liana_test %>%
  liana_dotplot(target_groups = c("high MYC_TGA", "intermediate MYC_TGA", "low MYC_TGA"), source_groups = c("Myeloid cells"), ntop = 10) +
  scale_x_discrete(guide = guide_axis(angle = 90)) + theme(strip.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) + theme(axis.text.x = element_text(size = 13))

filtered_liana_test %>%
  liana_dotplot(target_groups = c("high MYC_TGA", "intermediate MYC_TGA", "low MYC_TGA"), source_groups = c("T cells"), ntop = 10) +
  scale_x_discrete(guide = guide_axis(angle = 90)) + theme(strip.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) + theme(axis.text.x = element_text(size = 13))

filtered_liana_test_target %>%
  liana_dotplot(target_groups = c("Myeloid cells"), source_groups = c("high MYC_TGA", "intermediate MYC_TGA", "low MYC_TGA"), ntop = 10) +
  scale_x_discrete(guide = guide_axis(angle = 90)) + theme(strip.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) + theme(axis.text.x = element_text(size = 13))

filtered_liana_test_target %>%
  liana_dotplot(target_groups = c("T cells"), source_groups = c("high MYC_TGA", "intermediate MYC_TGA", "low MYC_TGA"), ntop = 10) +
  scale_x_discrete(guide = guide_axis(angle = 90)) + theme(strip.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) + theme(axis.text.x = element_text(size = 13))

filtered_liana_test_target %>%
  liana_dotplot(target_groups = c("high MYC_TGA", "intermediate MYC_TGA", "low MYC_TGA"), source_groups = c("high MYC_TGA", "intermediate MYC_TGA", "low MYC_TGA"), ntop = 10) +
  scale_x_discrete(guide = guide_axis(angle = 90)) + theme(strip.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) + theme(axis.text.x = element_text(size = 13))

filtered_liana_test %>%
  liana_dotplot(target_groups = c("high MYC_TGA", "intermediate MYC_TGA", "low MYC_TGA"), source_groups = c("high MYC_TGA", "intermediate MYC_TGA", "low MYC_TGA"), ntop = 10) +
  scale_x_discrete(guide = guide_axis(angle = 90)) + theme(strip.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) + theme(axis.text.x = element_text(size = 13))

################################################################################
# RS Phase

unique(main_obj_RS$MYC_activation_category)
Idents(main_obj_RS) = main_obj_RS$MYC_activation_category
liana_test <- liana_wrap(main_obj_RS)
liana_test <- liana_test %>% liana_aggregate()

# Filter LIANA results for RS phase
filtered_liana_test <- liana_test %>%
  group_by(source, ligand.complex, receptor.complex) %>%
  mutate(std_dev_natmi = sd(natmi.edge_specificity, na.rm = TRUE), count = n()) %>%
  filter(count == 1 | std_dev_natmi > 0.1) %>%
  ungroup()

filtered_liana_test_target <- liana_test %>%
  group_by(target, ligand.complex, receptor.complex) %>%
  mutate(std_dev_natmi = sd(natmi.edge_specificity, na.rm = TRUE), count = n()) %>%
  filter(count == 1 | std_dev_natmi > 0.1) %>%
  ungroup()

# Generate heatmap for Figure 2d
filtered_liana_heatmap <- liana_test %>%
  filter(!(source == "Myeloid cells" & target == "Myeloid cells") &
           !(target == "Myeloid cells" & source == "Myeloid cells") &
           !(source == "Myeloid cells" & target == "T cells") &
           !(target == "Myeloid cells" & source == "T cells") &
           !(source == "T cells" & target == "Myeloid cells") &
           !(target == "T cells" & source == "Myeloid cells")) %>%
  ungroup()
heat_freq(filtered_liana_heatmap, pallette = c("blue", "white", "violetred2"))

# Generate dot plots for Figures 2f-g and Supplementary Figure 3g
filtered_liana_test %>%
  liana_dotplot(target_groups = c("high MYC_TGA", "intermediate MYC_TGA", "low MYC_TGA"), source_groups = c("Myeloid cells"), ntop = 10) +
  scale_x_discrete(guide = guide_axis(angle = 90)) + theme(strip.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) + theme(axis.text.x = element_text(size = 13))

filtered_liana_test %>%
  liana_dotplot(target_groups = c("high MYC_TGA", "intermediate MYC_TGA", "low MYC_TGA"), source_groups = c("T cells"), ntop = 10) +
  scale_x_discrete(guide = guide_axis(angle = 90)) + theme(strip.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) + theme(axis.text.x = element_text(size = 13))

filtered_liana_test_target %>%
  liana_dotplot(target_groups = c("Myeloid cells"), source_groups = c("high MYC_TGA", "intermediate MYC_TGA", "low MYC_TGA"), ntop = 10) +
  scale_x_discrete(guide = guide_axis(angle = 90)) + theme(strip.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) + theme(axis.text.x = element_text(size = 13))

filtered_liana_test_target %>%
  liana_dotplot(target_groups = c("T cells"), source_groups = c("high MYC_TGA", "intermediate MYC_TGA", "low MYC_TGA"), ntop = 10) +
  scale_x_discrete(guide = guide_axis(angle = 90)) + theme(strip.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) + theme(axis.text.x = element_text(size = 13))

filtered_liana_test_target %>%
  liana_dotplot(target_groups = c("high MYC_TGA", "intermediate MYC_TGA", "low MYC_TGA"), source_groups = c("high MYC_TGA", "intermediate MYC_TGA", "low MYC_TGA"), ntop = 10) +
  scale_x_discrete(guide = guide_axis(angle = 90)) + theme(strip.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) + theme(axis.text.x = element_text(size = 13))

filtered_liana_test %>%
  liana_dotplot(target_groups = c("high MYC_TGA", "intermediate MYC_TGA", "low MYC_TGA"), source_groups = c("high MYC_TGA", "intermediate MYC_TGA", "low MYC_TGA"), ntop = 10) +
  scale_x_discrete(guide = guide_axis(angle = 90)) + theme(strip.text.x = element_text(size = 10)) + theme(axis.text.y = element_text(size = 10)) + theme(axis.text.x = element_text(size = 13))
