# Load necessary libraries
library(tximport)
library(DESeq2)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(data.table)
library(dplyr)
library(GSVA)

# Load gene and transcript data
ensembl <- read.table("ensembl_genes_transcripts_v100.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tx2gene <- ensembl[, c(4, 2)]

# Load metadata
sampleTable <- read.table("RNAseq_annotation.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sampleTable$fileName <- paste0("kalisto/", sampleTable$Sample, "_abundance.tsv")

# Import transcript data using tximport
txi <- tximport(files = sampleTable$fileName, type = "kallisto", txIn = TRUE, txOut = FALSE, tx2gene = tx2gene)

# Prepare DESeq2 dataset
sampleTable$Case <- factor(as.character(sampleTable$Case), levels = unique(sampleTable$Case))
sampleTable$Diagnosis <- factor(sampleTable$Diagnosis, levels = c("CLL", "RT"))
dds_CLLRT <- DESeqDataSetFromTximport(txi = txi, colData = sampleTable, design = ~ Case + Diagnosis)
keep <- rowSums(counts(dds_CLLRT)) >= 10
dds_CLLRT <- dds_CLLRT[keep,]
dds_CLLRT <- DESeq(dds_CLLRT)
vsd <- vst(dds_CLLRT, blind = FALSE)
vsd_CLLRT <- assay(vsd)
colnames(vsd_CLLRT) <- colData(dds_CLLRT)[, "Sample"]
vsd_CLLRT <- as.data.frame(vsd_CLLRT)

# Annotate data with gene symbols
vsd_CLLRT$EnsemblGeneStableID <- rownames(vsd_CLLRT)
vsd_CLLRT$HGNC.symbol <- ensembl$HGNC.symbol[match(vsd_CLLRT$EnsemblGeneStableID, ensembl$Gene.stable.ID.version)]
vsd_CLLRT <- vsd_CLLRT[!duplicated(vsd_CLLRT$HGNC.symbol),]
rownames(vsd_CLLRT) <- vsd_CLLRT$HGNC.symbol
vsd_CLLRT <- vsd_CLLRT[, -c(13:14)]

# Load gene signature
signature <- fread('MYCtargetgenes')
signature <- as.data.frame(signature)

# GSVA analysis
scores <- gsva(as.matrix(vsd_CLLRT), signature, method = c("ssgsea"), mx.diff = FALSE)
scores <- scores[, c(1, 3, 5, 7, 9, 11, 2, 4, 6, 8, 10, 12)]
df <- data.frame(Group = c(rep("Pre treatment", 6), rep("Resistance/Richter", 6)), Score = scores)
df$paired <- c(1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6)

################################################################################
# Heatmap

# Prepare data for heatmap
data_vvt <- vsd_CLLRT
data_vvt$gene <- rownames(data_vvt)
signature <- signature %>% filter(final_common_genes %in% data_vvt$gene)
data_vvt_plot <- data_vvt %>% filter(gene %in% signature$final_common_genes)

# Scale data for heatmap
z.mat <- t(scale(t(data_vvt_plot[, -c(13)]), center = TRUE, scale = TRUE))
z.mat <- na.omit(z.mat)

# Heatmap annotation
column_ha3 <- HeatmapAnnotation(Group = df$Group, score = df$Score)

# Generate heatmap
Heatmap(as.matrix(z.mat), clustering_distance_columns = "maximum", show_row_names = FALSE,
        clustering_method_columns = "complete", show_column_names = TRUE, column_names_gp = gpar(fontsize = 7),
        top_annotation = column_ha3)

################################################################################
# Boxplot

# Update group names
df$Group <- gsub("Pre treatment", "pre-RS", df$Group)
df$Group <- gsub("Resistance/Richter", "RS", df$Group)

# Figure 1c
ggplot(df, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot() +
  geom_point(aes(fill = Group, group = paired), size = 5, shape = 21) +
  geom_line(aes(group = paired), size = 2, color = 'gray', alpha = 0.6) +
  stat_compare_means(label.x = 1.2, label.y = 6.8) +
  scale_fill_manual(values = c("white", "white"), drop = FALSE) +
  theme_minimal() +
  labs(y = "MYC_activation") +
  theme(legend.position = "none")

################################################################################
# Correlation with MYC expression

df_MYC <- cbind(df, df1$MYC)
cor(df_MYC$Score, df_MYC$`df1$MYC`)

# Supplementary figure 2e
ggscatter(df_MYC, x = "df1$MYC", y = "Score") +
  stat_cor(method = "pearson", label.x = 8.3, label.y = 5.5) +
  theme_minimal() +
  labs(x = "MYC", y = "MYC_activation")

ggplot(data = df_MYC, mapping = aes(x = Score, y = `df1$MYC`)) +
  geom_point(color = "black", size = 3) +
  stat_cor(method = "pearson")
