library(data.table)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(GSVA)
library(ggplot2)




############################################## ICGC
############################
head(df) # the df table has been generated in Bulk_RNAseq_ICGC

library(msigdbr)
h_gene_sets = msigdbr(species = "human", category = "H")
unique(h_gene_sets$gs_name)

# Step 1: Create an empty results table
results_table2 <- data.frame(Pathway = character(), Correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)
# Step 2: Calculate signatures and correlations for each pathway
for (pathway in unique(h_gene_sets$gs_name)) {
  # Subset h_gene_sets for the current pathway
  pathway_data <- h_gene_sets[h_gene_sets$gs_name == pathway, ]
  
  # Extract gene symbols for the current pathway
  gene_symbols <- pathway_data$human_gene_symbol
  your_count_matrix= count_matrix
  gene_symbols= as.data.frame(gene_symbols)
  
  # Calculate signatures using gsva
  signature <- gsva(as.matrix(your_count_matrix),gene_symbols, method = c("ssgsea"), mx.diff = FALSE)
  signature <- data.frame(t(signature))
 
  
  # Check if row names are in correct order and reorder if needed
  if (!identical(rownames(signature), rownames(df))) {
    signature <- signature[rownames(df), ]
  }
  
  
  # Calculate correlation
  correlation_result <- cor.test( df$Score, signature$gene_symbols)
  
  # Save results to the table
  result_row <- data.frame(
    Pathway = pathway,
    Correlation = correlation_result$estimate,
    p_value = correlation_result$p.value
  )
  
  results_table2 <- rbind(results_table2, result_row)
}

# Print the results



# Read and calculate correlation with KEGG_BCR
signature_bcr <- fread('HALLMARK_BCR_KEGG.txt')
signature_bcr <- as.data.frame(signature_bcr)
scores_bcr <- gsva(as.matrix(count_matrix), signature_bcr, method = c("ssgsea"), mx.diff = FALSE)
scores_bcr = as.data.frame(t(scores_bcr))

if (!identical(rownames(scores_bcr), rownames(df))) {
  scores_bcr <- scores_bcr[rownames(df), ]
}

correlation_result_bcr <- cor.test(scores_bcr$gene, df$Score)


# Extract correlation and p-value for KEGG_BCR
cor_bcr <- correlation_result_bcr$estimate
p_value_bcr <- correlation_result_bcr$p.value

# Add the correlations and p-values to the results_table
results_table2 <- rbind(
  results_table2,
  data.frame(
    Pathway = "KEGG_BCR",  # Replace with the actual pathway name
    Correlation = cor_bcr,
    p_value = p_value_bcr
  )
)

# Print the updated results_table
print(results_table2)



#################################### RS
#############################################
df <- fread("df_RS.txt")
row.names(df)= df$V1

count_matrix=data_vvt #this table has been produced in Bulk_RNAseq_RT_MYC.R


library(GSVA)
library(dplyr)


library(msigdbr)
h_gene_sets = msigdbr(species = "human", category = "H")
unique(h_gene_sets$gs_name)


# Step 1: Create an empty results table
results_table3 <- data.frame(Pathway = character(), Correlation = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Step 2: Calculate signatures and correlations for each pathway
for (pathway in unique(h_gene_sets$gs_name)) {
  # Subset h_gene_sets for the current pathway
  pathway_data <- h_gene_sets[h_gene_sets$gs_name == pathway, ]
  
  # Extract gene symbols for the current pathway
  gene_symbols <- pathway_data$human_gene_symbol
  your_count_matrix= count_matrix
  gene_symbols= as.data.frame(gene_symbols)
  
  # Calculate signatures using gsva
  signature <- gsva(as.matrix(your_count_matrix),gene_symbols, method = c("ssgsea"), mx.diff = FALSE)
  signature=signature[,c(1,3,5,7,9,11,2,4,6,8,10,12)]
  signature <- data.frame(Group = c(rep("Pre treatment", 6), rep("Resistance/Richter", 6)), Score = signature)
  signature$paired= c(1,2,3,4,5,6,1,2,3,4,5,6)
  
  
  # Check if row names are in correct order and reorder if needed
  if (!identical(rownames(signature), rownames(df))) {
    signature <- signature[rownames(df), ]
  }
  
  
  # Calculate correlation
  correlation_result <- cor.test( df$Score, signature$Score,)
  
  # Save results to the table
  result_row <- data.frame(
    Pathway = pathway,
    Correlation = correlation_result$estimate,
    p_value = correlation_result$p.value
  )
  
  results_table3 <- rbind(results_table3, result_row)
}



# Read and calculate correlation with KEGG_BCR
signature_bcr <- fread('HALLMARK_BCR_KEGG.txt')
signature_bcr <- as.data.frame(signature_bcr)
scores_bcr <- gsva(as.matrix(count_matrix), signature_bcr, method = c("ssgsea"), mx.diff = FALSE)

scores_bcr=scores_bcr[,c(1,3,5,7,9,11,2,4,6,8,10,12)]
scores_bcr <- data.frame(Group = c(rep("Pre treatment", 6), rep("Resistance/Richter", 6)), Score = scores_bcr)
scores_bcr$paired= c(1,2,3,4,5,6,1,2,3,4,5,6)

if (!identical(rownames(scores_bcr), rownames(df))) {
  scores_bcr <- scores_bcr[rownames(df), ]
}

correlation_result_bcr <- cor.test(scores_bcr$Score, df$Score)


# Extract correlation and p-value for KEGG_BCR
cor_bcr <- correlation_result_bcr$estimate
p_value_bcr <- correlation_result_bcr$p.value

# Add the correlations and p-values to the results_table
results_table3 <- rbind(
  results_table3,
  data.frame(
    Pathway = "KEGG_BCR",  # Replace with the actual pathway name
    Correlation = cor_bcr,
    p_value = p_value_bcr
  )
)

# Print the updated results_table




# 1. Add a "Group" column to both tables
results_table2$Group <- "1.CLL"
results_table3$Group <- "2.RS"
# 2. Check if the order of "Pathway" in results_table2 and results_table is the same
if (identical(results_table2$Pathway, results_table3$Pathway)) {
  # 3. Combine the two tables one after the other, adding "_CLL" and "_Venetoclax" to each row in the "Pathway" column
  combined_table <- rbind( results_table2[-51,], results_table3[-51,]
  )
  
  # 4. Create a dotplot with size representing the correlation and color representing the p-value
  library(ggplot2)
  
  ggplot(combined_table, aes(x = Pathway, y = Group,  size = abs(Correlation),color = Correlation)) + #size = -log10(p_value),
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal() + coord_flip() 
   
} else {
  cat("Pathway order is not the same in results_table2 and results_table.\n")
}



combined_table$FDR <- p.adjust(combined_table$p_value, method = "fdr")




# Assuming combined_table is your data frame
filtered_table <- combined_table %>%
  group_by(Pathway) %>%
  filter(any(abs(Correlation) >= 0.6))

filtered_table <- combined_table %>%
  group_by(Pathway) %>%
  filter(any(FDR < 0.01)) %>%
  filter(any(abs(Correlation) >= 0.4))


nrow(filtered_table)/2
nrow(combined_table)/2

library(ggplot2)

ggplot(filtered_table, aes(x = Pathway, y = Group, color = Correlation, shape = abs(Correlation) > 0.7)) +
  geom_point(aes(size = 2)) +
  scale_colour_gradient2(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.2)) +
  guides(colour = guide_coloursteps(show.limits = TRUE)) +
  geom_text(aes(label = ifelse(FDR < 0.01, "*", "")),
            position = position_nudge(x = 0.1, y = 0.1),
            size = 3, hjust = 0, vjust = 0.5, color = "black") +
  theme_minimal() +
  coord_flip()


