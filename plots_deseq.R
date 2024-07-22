library(ggbiplot)
pca_Tissue <- function(dds, meta, Tissue){
  # Perform variance stabilizing transformation
  object <- vst(dds)
  
  # Calculate variance for each gene
  rv <- rowVars(assay(object))
  
  # Select top n genes by variance
  ntop <- 500
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # Perform PCA on the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # Prepare PCA data for plotting
  pca_data <- as.data.frame(pca$x)
  pca_data$Tissue <- Tissue
  pca_data$Sample <- meta$Sample
  
  # Return PCA data and loadings
  return(list(pca_data = pca_data, loadings = pca$rotation))
}

# Define a list of tissues and corresponding DESeq2 objects
tissues <- c("Heart", "Liver", "Lungs","Pect","Gut1","Gut2","Gut3")  # Add all tissues
dds_list <- list(dds_Heart, dds_Liver, dds_Lungs,dds_Pect,dds_Gut1,dds_Gut2,dds_Gut3)  # List of DESeq2 objects for each tissue

# Initialize empty lists to store PCA results
pca_results <- list()

# Run PCA for each tissue and store results
for (i in seq_along(tissues)) {
  tissue <- tissues[i]
  dds <- dds_list[[i]]
  meta <- colData(dds)  # Assuming colData contains sample metadata
  
  pca_results[[tissue]] <- pca_Tissue(dds, meta, tissue)
}


# Combine PCA results into one data frame
combined_pca_data <- bind_rows(lapply(pca_results, `[[`, "pca_data"))

# Plot combined PCA
ggplot(combined_pca_data, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 3) +
  labs(title = "PCA of All Tissues",
       x = "PC1",
       y = "PC2") +
  theme_minimal() +
  scale_color_viridis_d() +my_theme2


rld <- rlog(dds)
plotPCA(rld)


# Convert results to a data frame
res_df <- as.data.frame(res_Heart_ND)

# Create a column for -log10(p-value)
res_df$logPval <- -log10(res_df$pvalue)

# Create a column for significance (optional)
res_df$significant <- ifelse(res_df$pvalue < 0.05, "Significant", "Not Significant")

# For volcano plot: filter or adjust as needed
# Ensure gene names are included if required (e.g., rownames(res_heart))
res_df$Gene <- rownames(res_df)
res_df$color <- with(res_df, ifelse(
  significant == "Significant" & abs(log2FoldChange) > 0.58, 
  "#FDE725FF",  # Yellow for significant and log2FoldChange > 0.58 or < -0.58
  "#21908CFF"   # Blue for all other cases
))
ggplot(res_df, aes(x = log2FoldChange, y = logPval, color = color)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_identity() +  # Use colors defined in the data frame
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "right") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = '#440154FF') +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = '#440154FF') + my_theme2
# Plot volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = logPval, color = significant)) +
  geom_point(alpha = 0.4, size = 2) +
  scale_color_manual(values = c("Significant" = "#FDE725FF", "Not Significant" ="#21908CFF")) +  # Choose an option like "viridis", "magma", etc.
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "right") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = '#440154FF') +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = '#440154FF') +
  my_theme2

('D' = "#FDE725FF", 'T' = '#21908CFF',
  "N"="#440154FF")
viridis_colors <- viridis(256, option = "viridis")
head(viridis_colors)




res_df$color <- with(res_df, ifelse(
  log2FoldChange > 0.58, "#FDE725FF",  # Blue for log2FoldChange > 0.58
  ifelse(log2FoldChange < -0.58, "#6A0D91FF",  # Purple for log2FoldChange < -0.58
         NA)  # No color for other points
))
ggplot(res_df, aes(x = log2FoldChange, y = logPval, color = color)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c("#FDE725FF" = "#FDE725FF", "#6A0D91FF" = "#6A0D91FF"),
    na.value = "grey",  # Color for NA (non-significant points)
    breaks = c("#FDE725FF", "#6A0D91FF", "grey"),
    labels = c("Upregulated", "Downregulated", "None")
  ) +  # Use grey for colorless points (optional)
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "right") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = '#440154FF') +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = '#440154FF')+
  my_theme2





