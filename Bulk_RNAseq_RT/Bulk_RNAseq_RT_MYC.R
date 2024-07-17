library("tximport")
library("DESeq2")
library("org.Hs.eg.db")
library("ggplot2")
library("ggrepel")
library("reshape2")
library("gridExtra")
library("ComplexHeatmap")
library("clusterProfiler")
library("circlize")
library("ggpubr")

getwd()
analysis.path=getwd()  ;analysis.path
analysis.path <- "E:\\cnag\\MCL_CLL\\analysis\\bulk_Richter" 
setwd(analysis.path)
getwd()

gc()
list.files()

ensembl <- read.table(paste0("ensembl_genes_transcripts_v100.txt"), header = T, sep = "\t", stringsAsFactors = F)
tx2gene <- ensembl[,c(4,2)]


# PCA

# sampleTable
sampleTable <- read.table(paste0( "RNAseq_annotation.txt"), sep = "\t", header = T, stringsAsFactors = F)
sampleTable$fileName <- paste0("kalisto/", sampleTable$Sample, "_abundance.tsv")
# Import using tximport package
txi <- tximport(files = sampleTable$fileName, type = "kallisto", txIn = TRUE, txOut = FALSE, tx2gene = tx2gene)
# DESeq
sampleTable$Case <- factor(as.character(sampleTable$Case), levels = unique(sampleTable$Case))
sampleTable$Diagnosis <- factor(sampleTable$Diagnosis, levels = c("CLL", "RT"))
dds_CLLRT <- DESeqDataSetFromTximport(txi = txi, colData = sampleTable, design = ~ Case + Diagnosis)
keep <- rowSums(counts(dds_CLLRT)) >= 10
dds_CLLRT <- dds_CLLRT[keep,]
dds_CLLRT <- DESeq(dds_CLLRT)
vsd <- vst(dds_CLLRT, blind=FALSE)
vsd_CLLRT <- assay(vsd)
colnames(vsd_CLLRT) <- colData(dds_CLLRT)[,c("Sample")]
vsd_CLLRT = as.data.frame(vsd_CLLRT)


vsd_CLLRT$EnsemblGeneStableID <- rownames(vsd_CLLRT)
vsd_CLLRT$HGNC.symbol <- ensembl$HGNC.symbol[match(vsd_CLLRT$EnsemblGeneStableID, ensembl$Gene.stable.ID.version)]

vsd_CLLRT <- vsd_CLLRT[!duplicated(vsd_CLLRT$HGNC.symbol),]
row.names(vsd_CLLRT) = vsd_CLLRT$HGNC.symbol
vsd_CLLRT =vsd_CLLRT[,-c(13:14)]
# Load the gene signature as a list of gene symbols

signature <- fread('E:\\cnag\\BCLLatlas\\data\\MYCtargets_285MYC3ou4NoNotch.txt')
signature= as.data.frame(signature)

scores <- gsva(as.matrix(vsd_CLLRT), signature,method=c("ssgsea"), mx.diff=FALSE)



scores=scores[,c(1,3,5,7,9,11,2,4,6,8,10,12)]
df <- data.frame(Group = c(rep("Pre treatment", 6), rep("Resistance/Richter", 6)), Score = scores)
df$paired= c(1,2,3,4,5,6,1,2,3,4,5,6)



#########################################
####heatmap
data_vvt=vsd_CLLRT
data_vvt$gene= row.names(data_vvt)
signature = signature %>% filter(final_common_genes %in% data_vvt$gene)
data_vvt_plot= data_vvt %>% filter (gene %in% signature$final_common_genes)
colnames(data_vvt_plot)


library(ComplexHeatmap)
library(gplots)

df$Sample_ID= row.names(df)
design_data_ALL<-df 


colnames(data_vvt_plot)

z.mat <- t(scale(t(data_vvt_plot[,-c(13)]), center=TRUE, scale=TRUE))
head(z.mat)

column_ha3 = HeatmapAnnotation(Group = design_data_ALL$Group, #No_of_reads = anno_points(design_data_ALL$reads),
                               #IGHV= design_data_ALL$IGHV, tris12= design_data_ALL$tris12,
                               score=design_data_ALL$Score
)

z.mat <- na.omit(z.mat)
#"euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
Heatmap(as.matrix(z.mat), clustering_distance_columns = "maximum", show_row_names =F,
        clustering_method_columns = "complete", show_column_names = T, column_names_gp = gpar(fontsize = 7),
        top_annotation = column_ha3)


####################################################################################


ggplot(df, aes(x = Group, y = Score)) + geom_boxplot()


# Replace values in the "Group" column
df$Group <- gsub("Pre treatment", "pre-RS", df$Group)
df$Group <- gsub("Resistance/Richter", "RS", df$Group)

# Print the updated dataframe
print(df)

# Figure 1c
ggplot(df, aes(x = Group, y = Score, fill=Group)) +
  geom_boxplot()+
  geom_point(aes(fill=Group,group=paired),size=5,shape=21)+
  geom_line(aes(group=paired), size=2, color='gray', alpha=0.6) + 
  stat_compare_means(label.x = 1.2, 
                     label.y = 6.8)+ 
  ggplot2::scale_fill_manual(values = c("white", "white"), drop = FALSE)+ theme_minimal()+
  labs(y= "MYC_activation") +
  theme(legend.position = "none")


################################
### correlation with MYC expression

# Supplementary figure 2e
df_MYC= cbind(df, df1$MYC)
cor(df_MYC$Score, df_MYC$`df1$MYC`)
library(ggpubr)
ggplot(data = df_MYC, mapping = aes(x = Score, y = `df1$MYC`)) + 
  geom_point(color= "black", size=3) +stat_cor(method="pearson")

ggscatter(df_MYC, x = "df1$MYC", y = "Score", #color = "Group",
          #add = "reg.line",  # Add regressin line
          #add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          #conf.int = TRUE # Add confidence interval
)  + stat_cor(method = "pearson", label.x = 8.3, label.y = 5.5)+ theme_minimal() +
  labs(x="MYC", y="MYC_activation")