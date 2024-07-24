# Load necessary libraries
library(DESeq2)
library(dplyr)
library(tidyverse)
library(data.table)
library(ggplot2)
library(GSVA)
library(sva)
library(survminer)
library(survival)
library(ComplexHeatmap)
library(limma)
library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)

# Load data
# Data can be downloaded from https://cllmap.org/downloads.html
data <- fread("cllmap_rnaseq_counts_full.tsv")
data <- as.data.frame(data)

# Calculate the number of counts
colnames(data)
row.names(data) <- data[, 1]
data <- data[, -c(1:2)]
x <- t(data %>% summarize_if(is.numeric, sum, na.rm = TRUE))

# DESeq2 normalization
metadata <- fread("metadata.txt")
metadata <- as.data.frame(metadata)
metadata <- metadata %>% filter(cohort == "ICGC")
data <- data[!duplicated(data$Description), ]
row.names(data) <- data$Description
data <- data[, -c(1, 2)]
colData <- metadata
countData <- data

# Convert participant IDs to character vector
participants <- as.character(colData$participant_id)
matching_cols <- intersect(participants, colnames(countData))
countData_sub <- countData[, matching_cols]
colData <- colData %>% filter(participant_id %in% matching_cols)

# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = countData_sub, colData = colData, design = ~ tumor_molecular_subtype)
dds <- DESeq(dds)
res <- results(dds)
resultsNames(dds)
vsd <- vst(dds)
vvt <- assay(vsd)
vvt <- as.data.frame(vvt)
data <- as.matrix(vvt)

# Batch effect correction
x <- as.data.frame(x)
x$v2 <- row.names(x)
x <- x[, c(2, 1)]
df <- as.data.frame(x, stringsAsFactors = FALSE)
colnames(df) <- c("Sample", "Reads")
num_bins <- 10
breaks <- quantile(df$Reads, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)
df$Bin <- cut(df$Reads, breaks = breaks, labels = FALSE, include.lowest = TRUE)
annotation <- merge(colData, df, by.x = "participant_id", by.y = "Sample")
ctype <- annotation
row.names(ctype) <- ctype$participant_id
modcombat <- model.matrix(~tumor_molecular_subtype, data = ctype)
batch <- annotation$Bin
edata <- as.matrix(vvt)
combat_edata <- ComBat(dat = edata, batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)

# Scoring
markers <- list()
markers$MYC_285genes <- as.data.frame(fread('MYCtargetgenes.txt'))
data_vvt <- as.data.frame(combat_edata)
data_vvt$gene <- row.names(data_vvt)
signature <- markers$MYC_285genes
signature <- signature %>% filter(final_common_genes %in% data_vvt$gene)
data_vvt_plot <- data_vvt %>% filter(gene %in% signature$final_common_genes)
scores <- gsva(as.matrix(combat_edata), markers$MYC_285genes, method = c("ssgsea"), mx.diff = FALSE)
row.names(scores) <- "MYC_activation"
df <- data.frame(t(scores))
df$Sample_ID <- row.names(df)

write.csv(df, "df.csv") # needed for the correlation analysis
# Load additional metadata
cyto <- fread("cytogenetics_metadata.csv")
annotation_all <- merge(annotation, cyto, by = "participant_id")
genetic <- fread("genetics_metadata.csv")
annotation_all <- merge(annotation_all, genetic, by.x = "participant_id", by.y = "Tumor_Sample_Barcode")

# Add NOTCH1 mutation data
notch1_mutated_ids <- annotation_all %>% filter(Hugo_Symbol == "NOTCH1")
annotation_all$NOTCH1 <- ifelse(annotation_all$participant_id %in% notch1_mutated_ids$participant_id, "yes", "no")
df <- merge(df, annotation_all, by.x = "Sample_ID", by.y = "participant_id")

# Categorize MYC activation levels
df <- df %>%
  mutate(MYC_activation_category = case_when(
    MYC_activation <= quantile(MYC_activation, 0.25) ~ "low",
    MYC_activation <= quantile(MYC_activation, 0.75) ~ "intermediate",
    TRUE ~ "high"
  ))

# Figure 1e
df %>%
  ggplot(aes(x = tumor_type, y = MYC_activation)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1), aes(col = MYC_activation_category)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme_minimal() +
  theme(legend.position = "none")

# Figure 1a
df %>%
  filter(!is.na(IGHV_mut)) %>%
  ggplot(aes(x = IGHV_mut, y = MYC_activation)) +
  geom_boxplot() +
  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_minimal()

# Figure 1b
df %>%
  ggplot(aes(x = as.character(tri_12), y = MYC_activation)) +
  geom_boxplot() +
  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_minimal() +
  scale_x_discrete(breaks = c("0", "1"), labels = c("non-trisomy 12", "trisomy 12"))

# Supplementary Figure 1c
df %>%
  filter(!is.na(IGHV_mut)) %>%
  ggplot(aes(x = as.character(tri_12), y = MYC_activation)) +
  geom_boxplot() +
  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_minimal() +
  scale_x_discrete(breaks = c("0", "1"), labels = c("non-trisomy 12", "trisomy 12")) +
  facet_grid(~tumor_molecular_subtype)

# Survival analysis
df_survival <- df %>%
  filter(!is.na(IGHV_mut)) %>%
  filter(treatment_status_at_samp == "Untreated")

fit <- survfit(Surv(ffs_days, ffs) ~ MYC_activation, data = df_survival)
ggsurvplot(fit, data = df_survival, size = 1, conf.int = TRUE, pval = TRUE, risk.table = TRUE,
           risk.table.col = "strata", risk.table.height = 0.25, ggtheme = theme_bw())

# Supplementary Figure 2f
df_survival <- df %>%
  filter(treatment_status_at_samp == "Untreated") %>%
  filter(MYC_activation_category %in% c("low", "high"))

fit <- survfit(Surv(os_days, os) ~ MYC_activation, data = df_survival)
ggsurvplot(fit, data = df_survival, size = 1, conf.int = TRUE, pval = TRUE, risk.table = TRUE,
           risk.table.col = "strata", risk.table.height = 0.25, ggtheme = theme_bw())

# Correlation with MYC expression
MYC_data <- as.data.frame(combat_edata)
MYC_data$gene <- row.names(MYC_data)
MYC_data <- MYC_data %>% filter(gene == "MYC")
MYC_data <- MYC_data[, -ncol(MYC_data)]
MYC_data <- data.frame(t(MYC_data))
MYC_data$Sample_ID <- row.names(MYC_data)

MYC_data <- merge(MYC_data, metadata, by.x = "Sample_ID", by.y = "participant_id")
MYC_table <- merge(MYC_data, df[, c(1, 2, 123)], by = "Sample_ID")

cor.test(MYC_table$MYC, MYC_table$MYC_activation)

# Supplementary Figure 2d
ggscatter(MYC_table, x = "MYC", y = "MYC_activation") +
  stat_cor(method = "pearson", label.x = 7.3, label.y = 26.5) +
  theme_minimal()

# Differential expression analysis
dfexp <- as.matrix(combat_edata)
df_compare <- df %>%
  filter(!is.na(IGHV_mut)) %>%
  filter(treatment_status_at_samp == "Untreated") %>%
  filter(cll_epitype != "unclassified") %>%
  filter(MYC_activation_category %in% c("high", "low"))

list <- as.character(df_compare$Sample_ID)
dfexp <- combat_edata[, list]
pheno <- factor(df_compare$MYC_activation_category)

phenoMat <- model.matrix(~pheno)
colnames(phenoMat) <- sub("^pheno", "", colnames(phenoMat))

fit <- lmFit(object = dfexp, design = phenoMat)
fit <- eBayes(fit)

degCLLs <- topTable(fit, number = nrow(dfexp), adjust.method = "fdr", sort.by = "p")
sign.table <- as.data.frame(degCLLs) %>% filter(adj.P.Val <= 0.05) %>% filter(abs(logFC) >= 1)
sign.table$ID <- row.names(sign.table)

dfexp <- as.data.frame(combat_edata[, list])
dfexp$ID <- row.names(dfexp)

sign.table.full <- merge(sign.table, dfexp, by.x = "ID", by.y = "ID")
z.mat <- t(scale(t(sign.table.full[, -c(1:7)]), center = TRUE, scale = TRUE))

design_data_ALL <- df_compare
column_ha3 <- HeatmapAnnotation(groups = design_data_ALL$MYC_activation_category,
                                IGHV = design_data_ALL$IGHV_mut, tris12 = as.character(design_data_ALL$tri_12),
                                col = list(groups = c("low" = "#619CFF", "high" = "#F8766D"),
                                           IGHV = c("mutated" = "yellow", "unmutated" = "orange"),
                                           tris12 = c("1" = 'gray3', "0" = "gray60")),
                                na_col = "gray90")

cl <- kmeans(z.mat, centers = 2)$cluster

# Figure 2g
Heatmap(as.matrix(z.mat), clustering_distance_columns = "euclidean", row_split = cl,
        clustering_method_columns = "complete", show_column_names = FALSE,
        top_annotation = column_ha3)

cut.info <- as.data.frame(cl)
cut.info$rowID <- row.names(cut.info)
sign.table.full$rowID <- row.names(sign.table.full)
sign.table.full.2 <- merge(cut.info, sign.table.full, by = "rowID")

# GSEA analysis
original_gene_list <- sign.table.full.2$logFC
names(original_gene_list) <- sign.table.full.2$ID

gene_list <- na.omit(original_gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)

msigdbr_df <- msigdbr(species = "human", category = "H")
pathwaysH <- split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)

fgseaRes <- fgsea(pathways = pathwaysH, gene_list)
gos <- fgseaRes %>% filter(pval <= 0.05)

table <- as.data.frame(gos)
table$log <- -log10(table$pval)

# Figure 2h
ggplot(table, aes(x = log, y = pathway, color = -NES, size = 6)) +
  geom_point(alpha = 0.8) +
  xlab("-log10(pvalue)") +
  ylab("Hallmark database") +
  theme_classic()

# HMCN1 mutation analysis
HMCN1_mutated_ids <- annotation_all %>% filter(Hugo_Symbol == "HMCN1")
annotation_all$HMCN1 <- ifelse(annotation_all$participant_id %in% HMCN1_mutated_ids$participant_id, "yes", "no")
df <- merge(df, annotation_all, by.x = "Sample_ID", by.y = "participant_id")

# Supplementary Figure 2b
df %>%
  filter(!is.na(HMCN1)) %>%
  ggplot(aes(x = HMCN1, y = MYC_activation)) +
  geom_boxplot() +
  stat_compare_means() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_minimal()
