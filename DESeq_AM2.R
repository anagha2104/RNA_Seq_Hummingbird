# Intro-------------------
####  Understanding DESeq for non-brain tissues
# Hellooooo
# Packages-------------

# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# BiocManager::install('DESeq2')
# BiocManager::install('EnhancedVolcano') 
# BiocManager::install('apeglm')                                                                        
# BiocManager::install("clusterProfiler")
# BiocManager::install("ComplexHeatmap")
# if (!require("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
# BiocManager::install("biomaRt")
# BiocManager::install("pathview")
# ComplexHeatmap, clusterProfiler (used BiocManager for these also cause R version problems)
# Enhanced Volcano plot help: 
# https:\\\\bioconductor.org\\packages\\devel\\bioc\\vignettes\\EnhancedVolcano\\inst\\doc\\EnhancedVolcano.html
# if (!("msigdbr" %in% installed.packages())) {
#   # Install this package if it isn't installed yet
#   BiocManager::install("msigdbr", update = FALSE)
# }
# if (!("ggupset" %in% installed.packages())) {
#   # Install this package if it isn't installed yet
#   BiocManager::install("ggupset", update = FALSE)
# }

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(EnhancedVolcano)
library(dplyr)
library(here)
library(viridis)
library(gridExtra)
library(apeglm)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(dendextend) # For making gene trees for heatmaps
library(ComplexHeatmap)
library(cowplot) # for plot_grid function
library(enrichR)
library(scales) #show_col(viridis_pal()(20))
library(DOSE)
library(pathview)
library(matrixStats)#cv
library(AnnotationDbi)
library(biomaRt)
library(org.Hs.eg.db)
library(org.Gg.eg.db)
library(PCAtools)
library(Rfast)
library(msigdbr)
library(ggupset)
library(AnnotationForge)
library(GenomeInfoDb)
library(GO.db)

# Read data in------------------
data <- read.csv(here("DESeq_Data_Mar2022_all tissues", "all tissues", "final", "RNASeq_rawcounts_data.csv"))
meta <- read.csv(here("DESeq_Data_Mar2022_all tissues", "all tissues", "final", "Meta_use_this.csv"))
# star_alignment <- read.csv(here("DESeq_Data_Mar2022_all tissues", "all tissues", "final", "star_alignment_plot.csv"))
# star <- read.csv(here("DESeq_Data_Mar2022_all tissues", "all tissues", "final", "star_gene_counts.csv"))
# foldchange <- read.csv(here("DESeq_Data_Mar2022_all tissues", "all tissues", "final", "TopFoldChanges.csv"))
# norm_counts_df <- read.csv(here("DESeq_Data_Mar2022_all tissues", "all tissues", "final", "RNASeq_NormCounts.csv"))


# Setting my themes--------------------
my_theme <- theme_classic(base_size = 30) + 
  theme(panel.border = element_rect(colour = "black", fill=NA)) + theme(legend.key.height = unit(2, "line"))

my_theme2 <- theme_classic(base_size = 20) + 
  theme(panel.border = element_rect(colour = "black", fill=NA)) + theme(legend.key.height = unit(2, "line"))

# Virids colors
my_gradient <- c("#823de9", "#7855ce", "#6e6eb2", "#648697", "#599e7c", "#4fb760", "#45cf45")
my_col_rainbows <- c("#f94144", "#f3722c", "#f8961e", "#f9844a",
                     "#f9c74f", "#90be6d", "#43aa8b", "#4d908e", "#577590", "#277da1")
mycols <- c("orange", "navy", "springgreen4", "mediumorchid", "#760431", "plum1", "springgreen")

# # removing outliers----------
# outlier_samples <- meta %>% filter(meta$BirdID=="CAAN10" |  meta$BirdID=="CAAN04") %>% pull(X)
# data <- data %>% dplyr::select(-outlier_samples)
# meta <- meta %>% filter(meta$BirdID!="CAAN10" & meta$BirdID!="CAAN04")



# Making meta2 and getting lists of tissue specific samples----------
meta2 <- meta %>%
  rename(Sample = X)

meta2$Metabolic_State <- as.factor(meta2$Metabolic_State)
meta2$Tissue_State <- paste0(meta2$Tissue, "_", meta2$Metabolic_State)
meta2$Tissue_State <- as.factor(meta2$Tissue_State)

# samples_summary <- summary(meta2$Tissue_State) %>% data.frame() %>% rownames_to_column("Sample_type")
# write.csv(samples_summary, here("../samples_summary.csv"), row.names = FALSE)


Gut1_samples <- meta2$Sample[meta2$Tissue=="Gut1"]
Gut2_samples <- meta2$Sample[meta2$Tissue=="Gut2"]
Gut3_samples <- meta2$Sample[meta2$Tissue=="Gut3"]
Liver_samples <- meta2$Sample[meta2$Tissue=="Liver"]
Heart_samples <- meta2$Sample[meta2$Tissue=="Heart"]
Lungs_samples <- meta2$Sample[meta2$Tissue=="Lungs"]
Pect_samples <- meta2$Sample[meta2$Tissue=="Pect"]



# Make data proper for deseq----------------
# doing it before subsetting saves steps
rownames(meta) <- unique(meta$X)
# rownames(meta) <- meta$X
rownames(data) <- unique(data$X)
data <- subset(data, select = -c(X))
meta <- subset(meta, select = -c(X))
meta <- meta[match(colnames(data), rownames(meta)),]

# ## Pre-filtering
# ## Taking out rows that have counts less than 10 reads total - recommended, not required, step
# keep <- rowSums(data) >=10
# data <- data[keep,]

# Make tissue based subsets  ----
# (6\\10\\2024, Aravind H helped with this. It is better to get dispersion estimates for a 
# single tissue and not across because some genes may not be expressed in all tissues)
meta_Gut1 <- filter(meta[meta$Tissue=="Gut1",])
meta_Gut2 <- filter(meta[meta$Tissue=="Gut2",])
meta_Gut3 <- filter(meta[meta$Tissue=="Gut3",])
meta_Heart <- filter(meta[meta$Tissue=="Heart",])
meta_Liver <- filter(meta[meta$Tissue=="Liver",])
meta_Lungs <- filter(meta[meta$Tissue=="Lungs",])
meta_Pect <- filter(meta[meta$Tissue=="Pect",])

data_Gut1 <- filter(data[,Gut1_samples])
data_Gut2 <- filter(data[,Gut2_samples])
data_Gut3 <- filter(data[,Gut3_samples])
data_Heart <- filter(data[,Heart_samples])
data_Liver <- filter(data[,Liver_samples])
data_Lungs <- filter(data[,Lungs_samples])
data_Pect <- filter(data[,Pect_samples])



# Checking if all ok-------
## Check that all column names from the data set are present as rownames in the metadata file
all(colnames(data) %in% rownames(meta))
## Check that the cols are in the same order in the data as the rownames in the metadata
all(colnames(data) == rownames(meta))



# one fn to rule them all :)---------
Tissue_Specific_DESeq <- function(data,meta,Tissue){
  ## Pre-filtering
  ## Taking out rows that have counts less than 10 reads total - recommended, not required, step
  keep <- rowSums(data) >=10
  data <- data[keep,]
  

  # Create DESeq2Dataset object and running DESeq()
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~Metabolic_State)
  
 
  ## Set the factor level to compare against (in our case, normothermy)
  dds$Metabolic_State <- factor(dds$Metabolic_State, levels=c("N","T","D"))
  dds$Metabolic_State <- relevel(dds$Metabolic_State, ref="N")
  
  ## Run DESeq
  dds <- DESeq(dds)
  
  # res_ND <- results(dds, contrast = c("Metabolic_State","D","N"))
  res_ND <- lfcShrink(dds, coef = "Metabolic_State_D_vs_N",type = "apeglm") %>% data.frame()
  # res_NT <- results(dds, contrast = c("Metabolic_State","T","N"))
  res_NT <- lfcShrink(dds, coef = "Metabolic_State_T_vs_N",type = "apeglm") %>% data.frame()
  
  # normalized counts
  normalized_counts <- counts(dds, normalized=TRUE) %>% 
    data.frame()

  # Pulling significant genes
  upreg_ND <-  res_ND %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble() %>%
    filter(log2FoldChange>0.58 & padj < 0.05 & baseMean > 10) %>%
    arrange(padj) %>%
    pull(gene)

  downreg_ND <-  res_ND %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble() %>%
    filter(log2FoldChange < -0.58 & padj < 0.05 & baseMean > 10) %>%
    arrange(padj) %>%
    pull(gene)

  upreg_NT <-  res_NT %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble() %>%
    filter(log2FoldChange>0.58 & padj < 0.05 & baseMean > 10) %>%
    arrange(padj) %>%
    pull(gene)

  downreg_NT <-  res_NT %>%
    data.frame() %>%
    rownames_to_column(var="gene") %>%
    as_tibble() %>%
    filter(log2FoldChange < -0.58 & padj < 0.05 & baseMean > 10) %>%
    arrange(padj) %>%
    pull(gene)

  
  #creating dynamic variables
  dds_name <- paste0("dds_", Tissue)
  res_ND_name <- paste0("res_", Tissue,"_ND")
  res_NT_name <- paste0("res_", Tissue,"_NT")
  nor_name <- paste0("norm_", Tissue)
  upreg_ND_name <- paste0("upreg", Tissue,"_ND")
  upreg_NT_name <- paste0("upreg", Tissue,"_NT")
  downreg_ND_name <- paste0("downreg", Tissue,"_ND")
  downreg_NT_name <- paste0("downreg", Tissue,"_NT")
  
  # Assign the objects globally
  assign(dds_name, dds, envir = .GlobalEnv)
  assign(res_ND_name, res_ND, envir = .GlobalEnv)
  assign(res_NT_name, res_NT, envir = .GlobalEnv)
  assign(nor_name, normalized_counts, envir = .GlobalEnv)
  assign(upreg_NT_name, upreg_NT, envir = .GlobalEnv)
  assign(upreg_ND_name, upreg_ND, envir = .GlobalEnv)
  assign(downreg_NT_name, downreg_NT, envir = .GlobalEnv)
  assign(downreg_ND_name, downreg_ND, envir = .GlobalEnv)
}


Tissue_Specific_DESeq(data_Heart,meta_Heart,Tissue = "Heart")
Sys.sleep(15)
Tissue_Specific_DESeq(data_Gut1,meta_Gut1,Tissue = "Gut1")
Sys.sleep(15)
Tissue_Specific_DESeq(data_Gut2,meta_Gut2,Tissue = "Gut2")
Sys.sleep(15)
Tissue_Specific_DESeq(data_Liver,meta_Liver,Tissue = "Liver")
Sys.sleep(15)
Tissue_Specific_DESeq(data_Lungs,meta_Lungs,Tissue = "Lungs")
Sys.sleep(15)
Tissue_Specific_DESeq(data_Pect,meta_Pect,Tissue = "Pect")
Sys.sleep(15)
Tissue_Specific_DESeq(data_Gut3,meta_Gut3,Tissue = "Gut3")

# #res shrinking
# res_Heart_ND2 <- lfcShrink(dds_Heart,coef="Metabolic_State_D_vs_N",res=res_Heart_ND)
# res_Heart_ND3 <- lfcShrink(dds_Heart,coef="Metabolic_State_D_vs_N",type = "apeglm")

# used to write a bunch of tables so I dont have run everything again
res_Heart_ND <- res_Heart_ND %>% rownames_to_column("gene")
norm_Heart <- norm_Heart %>% rownames_to_column("gene")
write.csv(res_Heart_ND, "res_Heart_ND.csv", row.names = FALSE)
write.csv(norm_Heart, "norm_Heart.csv", row.names = FALSE)
save(upregHeart_ND, downregHeart_ND, file = here("lists.RData"))


# Following lines can be used to copy a list to clipboard and then excel
df <- data.frame(my_column = unlist(upregHeart_ND))
write.table(df, file = "clipboard-128", sep = "/t", row.names = FALSE, col.names = FALSE)

# Convert the list to a dataframe
df <- data.frame(my_column = unlist(upregHeart_ND))

# Add double quotes around each element
df$my_column <- paste0('"', df$my_column, '"')

# Convert the column to a comma-separated string
comma_separated_string <- paste(df$my_column, collapse = ",")

# Copy the string to the clipboard
writeClipboard(comma_separated_string)

# disp graph and data quality related graphs---------
plotDispEsts(
  dds_Heart,
  CV = T,
  genecol = "black",
  fitcol = "#FDE725FF",
  finalcol = "#440154FF") 

res <- lfcShrink(dds_Heart, coef = "Metabolic_State_D_vs_N",type = "apeglm")
plotMA(res,ylim=c(-2,2),
       main = "MA plot heart",
       colNonSig = "#440154FF",
       colSig = "#FDE725FF",
       colLine = "black")


# rld <- rlog(dds_Heart, blind=T)
# rld_mat <- assay(rld) 
# rld_cor <- cor(rld_mat)    ## cor() is a base R function
# 
# head(rld_cor)
# pheatmap(rld_cor, annotation = meta, color = heat.colors, border_color=NA, fontsize = 10, 
#          fontsize_row = 10, height=20)
# heat.colors <- viridis_pal()(20)
# pheatmap(rld_cor, annotation = meta, color = heat.colors, border_color=NA, fontsize = 10, 
#          fontsize_row = 10, height=20)

# volcano plots------------
res_df <- as.data.frame(res_Heart_ND)
res_df$logPval <- -log10(res_df$pvalue)
res_df$significant <- ifelse(res_df$pvalue < 0.05, "Significant", "Not Significant")
res_df$Gene <- rownames(res_df)
res_df$color <- with(res_df, ifelse(
  log2FoldChange > 0.58, "#FDE725FF",  
  ifelse(log2FoldChange < -0.58, "#6A0D91FF",  
         NA)  
))
ggplot(res_df, aes(x = log2FoldChange, y = logPval, color = color)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c("#FDE725FF" = "#FDE725FF", "#6A0D91FF" = "#6A0D91FF"),
    na.value = "grey", 
    breaks = c("#FDE725FF", "#6A0D91FF", "grey"),
    labels = c("Upregulated", "Downregulated", "None")
  ) +  guides(color = guide_legend(title = "Gene Expression Status"))+
  theme(axis.title.x = element_text(size = 6),  
        axis.title.y = element_text(size = 6),
        axis.text.x = element_text(size = 4),    
        axis.text.y = element_text(size = 6),
        legend.position = "right",
        legend.title = element_text(size = 8),  
        legend.text = element_text(size = 6)) +
  labs(title = "Volcano Plot Heart",
       x = "Log2 Fold Change",
       y = "-Log10(p-value)")  +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = 'black') +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = 'black')+
  my_theme2 

# bar plots----------
tissues <- c("Liver", "Lungs", "Heart", "Pect","Gut1","Gut2","Gut3")
upregulated <- c(47,61,33,68,100,59,15)
downregulated <- c(-83,-63,-29,-49,-284,-98,-34)

upreg_downreg_ND <- data.frame(
  Tissue = rep(tissues, 2),
  Count = c(upregulated, downregulated),
  Regulation = rep(c("Upregulated", "Downregulated"), each = length(tissues))
)

ggplot(upreg_downreg_ND, aes(x = Tissue, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "identity", width = 0.5) +  # Set width to make bars thinner
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  # Add horizontal line at y=0
  labs(title = "Number of upreg and downreg genes ND",
       x = "Tissue",
       y = "Number of Genes") +
  scale_fill_viridis_d(option = "viridis") +
  theme_minimal() +
  scale_y_continuous(labels = abs, 
                     breaks = seq(-max(data$Count), max(data$Count), by = 20)) +
  my_theme2+  theme(plot.title = element_text(size = 17),
  axis.title.x = element_text(size = 12),
  axis.title.y = element_text(size = 12),
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 8))+
  scale_y_continuous(limits = c(min(upreg_downreg_ND$Count), 
                                max(upreg_downreg_ND$Count)))

# Number_of_upreg_downreg_heart_ND
# for NT now
tissues <- c("Liver", "Lungs", "Heart", "Pect","Gut1","Gut2","Gut3")
upregulated <- c(14,21,4,13,24,11,0)
downregulated <- c(-17,-41,-3,-10,-67,-18,-3)

upreg_downreg_NT <- data.frame(
  Tissue = rep(tissues, 2),
  Count = c(upregulated, downregulated),
  Regulation = rep(c("Upregulated", "Downregulated"), each = length(tissues))
)


ggplot(upreg_downreg_NT, aes(x = Tissue, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "identity", width = 0.5) +  # Set width to make bars thinner
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +  # Add horizontal line at y=0
  labs(title = "Number of upreg and downreg genes NT",
       x = "Tissue",
       y = "Number of Genes") +
  scale_fill_viridis_d(option = "viridis") +
  theme_minimal() +
  scale_y_continuous(labels = abs, 
                     breaks = seq(-max(data$Count), max(data$Count), by = 20)) +
  my_theme2+  theme(plot.title = element_text(size = 17),
                    axis.title.x = element_text(size = 12),
                    axis.title.y = element_text(size = 12),
                    axis.text.x = element_text(size = 8),
                    axis.text.y = element_text(size = 8),
                    legend.title = element_text(size = 12),
                    legend.text = element_text(size = 8))+
  scale_y_continuous(limits = c(min(upreg_downreg_NT$Count), 
                                max(upreg_downreg_NT$Count)))










all_upreg_ND <-c(upregGut1_ND,upregGut2_ND,upregGut3_ND,upregHeart_ND,upregLiver_ND,upregLiver_ND,upregLungs_ND,upregPect_ND)
all_downreg_ND<- c(downregGut1_ND,downregGut2_ND,downregGut3_ND,downregHeart_ND,downregLiver_ND,downregLiver_ND,downregLungs_ND,downregPect_ND)
all_upreg_NT <-c(upregGut1_NT,upregGut2_NT,upregGut3_NT,upregHeart_NT,upregLiver_NT,upregLiver_NT,upregLungs_NT,upregPect_NT)
all_downreg_NT <- c(downregGut1_NT,downregGut2_NT,downregGut3_NT,downregHeart_NT,downregLiver_NT,downregLiver_NT,downregLungs_NT,downregPect_NT)


# Coeffcient of varience------------
# norm_Heart %>% data.matrix() %>% rowVars() %>% cbind( )
# cbind(norm_Heart, var = apply(df, 1, function(x) var(na.omit(x))))
# var(norm_Heart)
heart_N_samples <- meta2 %>%
  filter(Metabolic_State=="N"&Tissue=="Heart") %>% pull(Sample)
heart_D_samples <- meta2 %>%
  filter(Metabolic_State=="D"&Tissue=="Heart") %>% pull(Sample)
heart_T_samples <- meta2 %>% 
  filter(Metabolic_State=="T"&Tissue=="Heart") %>% pull(Sample)
norm_Heart[,heart_N_samples] %>% data.matrix() %>% rowVars()


norm_Heart$Dmean <- rowMeans(norm_Heart[,heart_D_samples]) %>% data.frame() 
norm_Heart$Dsd <- apply(norm_Heart[,heart_D_samples], 1, sd, na.rm=TRUE) %>% data.frame() 
norm_Heart$Dcv <- norm_Heart$Dsd/norm_Heart$Dmean

norm_Heart$Nmean <- rowMeans(norm_Heart[,heart_N_samples]) %>% data.frame() 
norm_Heart$Nsd <- apply(norm_Heart[,heart_N_samples], 1, sd, na.rm=TRUE) %>% data.frame() 
norm_Heart$Ncv <- norm_Heart$Nsd/norm_Heart$Nmean

norm_Heart$Tmean <- rowMeans(norm_Heart[,heart_T_samples]) %>% data.frame() 
norm_Heart$Tsd <- apply(norm_Heart[,heart_T_samples], 1, sd, na.rm=TRUE) %>% data.frame() 
norm_Heart$Tcv <- norm_Heart$Tsd/norm_Heart$Tmean

list_heart <-  filter(norm_Heart,Tcv > 0.40 | Dcv>0.40 | Ncv>0.40) %>%
  rownames_to_column(var="gene") %>% pull(gene)

Liver_N_samples <- meta2 %>%
  filter(Metabolic_State=="N"&Tissue=="Liver") %>% pull(Sample)
Liver_D_samples <- meta2 %>%
  filter(Metabolic_State=="D"&Tissue=="Liver") %>% pull(Sample)
Liver_T_samples <- meta2 %>% 
  filter(Metabolic_State=="T"&Tissue=="Liver") %>% pull(Sample)


norm_Liver$Dmean <- rowMeans(norm_Liver[,Liver_D_samples]) %>% data.frame() 
norm_Liver$Dsd <- apply(norm_Liver[,Liver_D_samples], 1, sd, na.rm=TRUE) %>% data.frame() 
norm_Liver$Dcv <- norm_Liver$Dsd/norm_Liver$Dmean

norm_Liver$Nmean <- rowMeans(norm_Liver[,Liver_N_samples]) %>% data.frame() 
norm_Liver$Nsd <- apply(norm_Liver[,Liver_N_samples], 1, sd, na.rm=TRUE) %>% data.frame() 
norm_Liver$Ncv <- norm_Liver$Nsd/norm_Liver$Nmean

norm_Liver$Tmean <- rowMeans(norm_Liver[,Liver_T_samples]) %>% data.frame() 
norm_Liver$Tsd <- apply(norm_Liver[,Liver_T_samples], 1, sd, na.rm=TRUE) %>% data.frame() 
norm_Liver$Tcv <- norm_Liver$Tsd/norm_Liver$Tmean

list_Liver <-  filter(norm_Liver,Tcv > 0.40 | Dcv>0.40 | Ncv>0.40) %>%
  rownames_to_column(var="gene") %>% pull(gene)

list_Liver %>% data.frame() %>% view()
length(unique(c(list_Liver,list_heart)))


Pect_N_samples <- meta2 %>%
  filter(Metabolic_State=="N"&Tissue=="Pect") %>% pull(Sample)
Pect_D_samples <- meta2 %>%
  filter(Metabolic_State=="D"&Tissue=="Pect") %>% pull(Sample)
Pect_T_samples <- meta2 %>% 
  filter(Metabolic_State=="T"&Tissue=="Pect") %>% pull(Sample)


norm_Pect$Dmean <- rowMeans(norm_Pect[,Pect_D_samples]) %>% data.frame() 
norm_Pect$Dsd <- apply(norm_Pect[,Pect_D_samples], 1, sd, na.rm=TRUE) %>% data.frame() 
norm_Pect$Dcv <- norm_Pect$Dsd/norm_Pect$Dmean

norm_Pect$Nmean <- rowMeans(norm_Pect[,Pect_N_samples]) %>% data.frame() 
norm_Pect$Nsd <- apply(norm_Pect[,Pect_N_samples], 1, sd, na.rm=TRUE) %>% data.frame() 
norm_Pect$Ncv <- norm_Pect$Nsd/norm_Pect$Nmean

norm_Pect$Tmean <- rowMeans(norm_Pect[,Pect_T_samples]) %>% data.frame() 
norm_Pect$Tsd <- apply(norm_Pect[,Pect_T_samples], 1, sd, na.rm=TRUE) %>% data.frame() 
norm_Pect$Tcv <- norm_Pect$Tsd/norm_Pect$Tmean

list_Pect <-  filter(norm_Pect,Tcv > 0.40 | Dcv>0.40 | Ncv>0.40) %>%
  rownames_to_column(var="gene") %>% pull(gene)

list_Pect %>% data.frame() %>% view()
length(unique(c(list_Pect,list_heart,list_Liver)))

norm_Heart$Dcv <- as.numeric(norm_Heart$Dcv)
norm_Heart$Dmean <- as.numeric(norm_Heart$Dmean)
hist(norm_Heart$Dcv)
ggplot(norm_Heart, aes(x=Dcv)) + 
  geom_histogram()

Dcv_vector <- unlist(norm_Heart$Dcv)
Dcv_numeric <- as.numeric(Dcv_vector)
Dmean_vector <- unlist(norm_Heart$Dmean)
Dmean_numeric <- as.numeric(Dmean_vector)
Ncv_vector <- unlist(norm_Heart$Ncv)
Ncv_numeric <- as.numeric(Ncv_vector)
Nmean_vector <- unlist(norm_Heart$Nmean)
Nmean_numeric <- as.numeric(Nmean_vector)

plot(log(Dmean_numeric),Dcv_numeric)

ggplot(norm_Heart, aes(x = Dcv, y = Dmean)) +
  geom_point()



# Heatmaps------------

annotation_liver <- meta_Liver %>% rownames_to_column(var="Sample") %>% 
  dplyr::select(Sample, Metabolic_State) %>% 
  mutate(Metabolic_State = fct_relevel(Metabolic_State, c("N", "T", "D"))) %>%
  arrange(Metabolic_State)
 
anno_colors <- list( Metabolic_State = c(N = "#f9c74f", T = "violet", D = "#599e7c"))



pheatmap(
  mat               = norm_Liver[upregLiver_ND,],
  # color             = inferno(10),
  border_color      = NA,
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  annotation_col    = annotation_liver,
  annotation_colors = anno_colors,
  drop_levels       = TRUE,
  fontsize          = 14,
  main              = "Liver upreg Heatmap"
)  


sampleDistMatrix_heart <- as.matrix(sampleDists_heart)
#rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
#colnames(sampleDistMatrix) <- NULL
anno_colors <- list(Metabolic_State = c(N = "#f9c74f", D = "#43aa8b", T = "#277da1"))
colors <- rev(brewer.pal(9, "BrBG"))
meta_heart2 <- meta_Heart
meta_heart2$Metabolic_State <- factor(meta_heart2$Metabolic_State, levels = c("N", "T", "D"))
meta_heart2 <- meta_heart2 %>%
  arrange(Metabolic_State)
annotation <- data.frame(Metabolic_State = factor (meta_heart2$Metabolic_State, levels = c("N", "T", "D")))
pheatmap(sampleDistMatrix_heart,
         clustering_distance_rows=sampleDists_heart,
         clustering_distance_cols=sampleDists_heart,
         cluster_rows = F,
         cluster_cols=F,
         col=colors,
         annotation_row = annotation,
         annotation_col = annotation,
         annotation_colors = anno_colors,
         show_colnames = T,
         show_rownames = T)


  
# PCA------------
 
pca_Tissue <- function(dds,meta,Tissue){
  # change x and y in biplot for diff pcs and include in loadings
  
  # Run variance stabilizing transformation on the counts.
  object <- vst(dds)
  
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # Top n genes by variance to keep.
  ntop <- 500
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # Loadings for the first two PCs.
  loadings <- pca$rotation[, c(3,6)]
  loadings_inside <- loadings %>% data.frame() %>% 
    rownames_to_column(var="gene")
  
  loadings_name <- paste0("loadings_", Tissue)
  assign(loadings_name, loadings_inside, envir = .GlobalEnv)
  
  p <- pca(assay(object), metadata = meta, removeVar = 0.1)
  
  biplot(p,lab = rownames(p$metadata), 
         x = 'PC1', y = 'PC2',
         colby = 'Metabolic_State',
         # showLoadings = TRUE,
         max.overlaps = 30,legendPosition = "right",
         colkey = c('D' = "#FDE725FF", 'T' = '#21908CFF',
                    "N"="#440154FF"),
         pointSize=4, title=paste("PCA",Tissue),
         axisLabSize=13,gridlines.major = F,gridlines.minor = F) 
  
  
}

pca_Tissue(dds_Liver,meta_Liver,"Liver")
pca_Tissue(dds_Heart,meta_Heart,"Heart")
pca_Tissue(dds_Pect,meta_Pect,"Pect")


# box plots and lists-----------

all_upreg_ND <- c(upregGut1_ND, upregGut2_ND, upregGut3_ND, upregHeart_ND,  upregLiver_ND, upregLungs_ND, upregPect_ND)

all_downreg_ND <- c(downregGut1_ND, downregGut2_ND, downregGut3_ND, downregHeart_ND, downregLiver_ND, downregLungs_ND, downregPect_ND)

# Create the source vectors
upreg_sources <- c(rep('upregGut1_ND', length(upregGut1_ND)), 
                   rep('upregGut2_ND', length(upregGut2_ND)), 
                   rep('upregGut3_ND', length(upregGut3_ND)), 
                   rep('upregHeart_ND', length(upregHeart_ND)), 
                   rep('upregLiver_ND', length(upregLiver_ND)), 
                   rep('upregLungs_ND', length(upregLungs_ND)), 
                   rep('upregPect_ND', length(upregPect_ND)))

downreg_sources <- c(rep('downregGut1_ND', length(downregGut1_ND)), 
                     rep('downregGut2_ND', length(downregGut2_ND)), 
                     rep('downregGut3_ND', length(downregGut3_ND)), 
                     rep('downregHeart_ND', length(downregHeart_ND)), 
                     rep('downregLiver_ND', length(downregLiver_ND)), 
                     rep('downregLungs_ND', length(downregLungs_ND)), 
                     rep('downregPect_ND', length(downregPect_ND)))


# Create the DataFrames
upreg_df <- data.frame(Value = all_upreg_ND, Source = upreg_sources)
downreg_df <- data.frame(Value = all_downreg_ND, Source = downreg_sources)

# Display the DataFrames
print(upreg_df)
print(downreg_df)

# Example function to pad data frames with NA rows to match the maximum number of rows
pad_with_na <- function(df, max_rows) {
  if (nrow(df) < max_rows) {
    pad_rows <- max_rows - nrow(df)
    df <- bind_rows(df, data.frame(matrix(NA, nrow = pad_rows, ncol = ncol(df))))
    names(df) <- names(df)
  }
  return(df)
}

# Check the number of rows in each data frame
nrows_list <- c(nrow(norm_Heart), nrow(norm_Liver), nrow(norm_Lungs), nrow(norm_Pect),
                nrow(norm_Gut1), nrow(norm_Gut2), nrow(norm_Gut3))

# Find the maximum number of rows
max_rows <- max(nrows_list)

# Pad data frames with NA rows to match the maximum number of rows
norm_Heart_padded <- pad_with_na(norm_Heart, max_rows)
norm_Liver_padded <- pad_with_na(norm_Liver, max_rows)
norm_Lungs_padded <- pad_with_na(norm_Lungs, max_rows)
norm_Pect_padded <- pad_with_na(norm_Pect, max_rows)
norm_Gut1_padded <- pad_with_na(norm_Gut1, max_rows)
norm_Gut2_padded <- pad_with_na(norm_Gut2, max_rows)
norm_Gut3_padded <- pad_with_na(norm_Gut3, max_rows)

# Example function to ensure all data frames have a 'Gene' column
prepare_df <- function(df) {
  df %>% 
    data.frame() %>% 
    rownames_to_column(var = "Gene")
}

# Prepare data frames by adding 'Gene' column
norm_Heart_prepared <- prepare_df(norm_Heart)
norm_Liver_prepared <- prepare_df(norm_Liver)
norm_Lungs_prepared <- prepare_df(norm_Lungs)
norm_Pect_prepared <- prepare_df(norm_Pect)
norm_Gut1_prepared <- prepare_df(norm_Gut1)
norm_Gut2_prepared <- prepare_df(norm_Gut2)
norm_Gut3_prepared <- prepare_df(norm_Gut3)

# Perform a full join on 'Gene' to combine all data frames
combined_data <- norm_Heart_prepared %>%
  full_join(norm_Liver_prepared, by = "Gene") %>%
  full_join(norm_Lungs_prepared, by = "Gene") %>%
  full_join(norm_Pect_prepared, by = "Gene") %>%
  full_join(norm_Gut1_prepared, by = "Gene") %>%
  full_join(norm_Gut2_prepared, by = "Gene") %>%
  full_join(norm_Gut3_prepared, by = "Gene")

# Restore the 'Gene' column as row names and remove the 'Gene' column
rownames(combined_data) <- combined_data$Gene
combined_data <- combined_data %>% select(-Gene)
combined_data <- combined_data %>%
  column_to_rownames(var = "Gene") 
if ("Gene" %in% colnames(combined_data)) {
  combined_data <- combined_data %>% dplyr::select(-Gene)
}

combined_data <- bind_cols(
  norm_Heart %>% data.frame(),
  norm_Liver %>% data.frame() ,
  norm_Lungs %>% data.frame(), 
  norm_Pect %>% data.frame(),
  norm_Gut1%>% data.frame(),
  norm_Gut2 %>% data.frame(),
  norm_Gut3 %>% data.frame() )
 
datlong <- combined_data %>% data.frame() %>% rownames_to_column(var = "Gene") %>% 
  gather(key = 'Sample', value= 'counts',-Gene) %>% 
  left_join(., meta2, by="Sample") %>%
  mutate(Metabolic_State = 
           fct_relevel(Metabolic_State, 
                       "N", "T", "D")) 

# Boxplot\\Violin plot 
# just change gene name and plot
datlong %>%
  filter(Gene == "PER2") %>%
  ggplot(., aes(y=counts, x=Metabolic_State)) +
  geom_boxplot() + 
  my_theme2 + facet_wrap(.~Tissue, scales = "free") +
  ggtitle("PER2 gene expression across metabolic states") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene counts") + xlab("Metabolic state") +my_theme2+
  scale_fill_viridis_d()

# Plot the boxplot with custom colors
datlong %>%
  filter(Gene == "RSRP1") %>%
  ggplot(aes(y = counts, x = Metabolic_State, fill = Metabolic_State)) +
  geom_boxplot() +
  my_theme2 +
  facet_wrap(. ~ Tissue, scales = "free") +
  ggtitle("RSRP1 Gene Expression Across Metabolic States") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Gene Counts") + 
  xlab("Metabolic State") +
  scale_fill_viridis_d()+
  theme(axis.title.x = element_text(size = 13),  
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 8),    
        axis.text.y = element_text(size = 8),
        legend.position = "right",
        legend.title = element_text(size = 12),  
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 17))


res_Lungs_ND_tb <- res_Lungs_ND %>% data.frame() %>% view()


res_Lungs_ND <- rownames_to_column(res_Lungs_ND, "gene")
res_Lungs_ND <- filter(res_Lungs_ND, padj != "NA" )
write.csv(res_Lungs_ND,"C:\\Users\\anagh\\Documents\\RNASeq\\res_Lungs_ND_tb.csv")

"C:\Users\anagh\Documents\RNASeq\res_Heart_ND_tb.csv"



# new stuff with cluster profiler----------

df <- res_Heart_ND_tb %>% mutate(diffexpressed = case_when(
  log2FoldChange > 0.58 & padj < 0.05 ~ 'UP',
  log2FoldChange < -0.58 & padj < 0.05 ~ 'DOWN',
  padj > 0.05 ~ 'NO',
  log2FoldChange < 0.58 & log2FoldChange > -0.58 ~'NO'
  
))#Add base mean ka threshold also.


df_filtered <- subset(df, diffexpressed != "NO")

deg_Heart_res_list <- split(df_filtered,df_filtered$diffexpressed)

# ORA old-------------
  
all_genes <- as.character(res_Lungs_ND_tb$gene)

## Extract significant results
sig <- dplyr::filter(res_Lungs_ND_tb, padj < 0.05)

sig_genes <- as.character(sig$gene)

ego <- enrichGO(gene = all_genes, 
                keyType = "SYMBOL",
                OrgDb = org.Gg.eg, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)

go_annotations <- enrichGO(gene         = sig_genes,
                           OrgDb        = org.Gg.eg.db,
                           universe = all_genes,
                           keyType      = "SYMBOL",
                           ont          = "ALL",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05)

cluster_summary <- data.frame(go_annotations)

ego <- enrichGO(gene = sig_genes, 
                universe = all_genes,
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
cluster_summary2 <- data.frame(ego)
dotplot(go_annotations, showCategory=10)

go_annotations <- enrichplot::pairwise_termsim(go_annotations)
emapplot(go_annotations, showCategory = 50)
OE_foldchanges <- sig$log2FoldChange

names(OE_foldchanges) <- sig$gene
cnetplot(go_annotations, 
         showCategory = 5, 
         color.params=list(foldChange=OE_foldchanges),
         vertex.label.font=6)

goplot(go_annotations)
buildGOmap(all_genes)

# annotationdbi and biomart----


# listMarts(mart = NULL, host="www.ensembl.org", path="/biomart/martservice",
#           includeHosts = FALSE, archive = FALSE, verbose = FALSE)
ensembl <- useMart("ensembl")


datasets <- listDatasets(ensembl)
view(datasets)
Arctic_ground_squ <- useDataset("uparryii_gene_ensembl", mart = ensembl)
Arctic_ground_squ_attr<- listAttributes(Arctic_ground_squ)
view(Arctic_ground_squ)
go_annotations_squ <- getBM(attributes = c("external_gene_name",
                                           'ensembl_gene_id', 
                                           "description",
                                           "entrezgene_id", 
                                           'go_id', 
                                           'name_1006',
                                           'namespace_1003',
                                           "definition_1006"),
                          filters = 'external_gene_name',
                          values = all_heart_genes,
                          mart = Arctic_ground_squ)

Zebra_finch <- useDataset("tguttata_gene_ensembl", mart = ensembl)
Zebra_finch_attr <- listAttributes(Zebra_finch)
view(Zebra_finch_attr)
Golden_eagle <- useDataset("acchrysaetos_gene_ensembl", mart = ensembl)
Golden_eagle_attr <- listAttributes(Golden_eagle)
view(Golden_eagle_attr)
Kakapo <- useDataset("shabroptila_gene_ensembl", mart = ensembl)
Kakapo_attr <- listAttributes(Kakapo)
view(Kakapo_attr)
Great_tit <- useDataset("pmajor_gene_ensembl", mart = ensembl)
Great_tit_attr <- listAttributes(Great_tit)
view(Great_tit_attr)
Collared_flycatcher <- useDataset("falbicollis_gene_ensembl", mart = ensembl)
Collared_flycatcher_attr <- listAttributes(Collared_flycatcher)
view(Collared_flycatcher_attr)
Canary <- useDataset("scanaria_gene_ensembl", mart = ensembl)
Canary_attr <- listAttributes(Canary)
view(Canary_attr)
Medium_ground_finch <- useDataset("gfortis_gene_ensembl", mart = ensembl)
Medium_ground_finch_attr <- listAttributes(Medium_ground_finch)
view(Medium_ground_finch_attr)


go_annotations <- getBM(attributes = c("external_gene_name",'ensembl_gene_id', 'go_id', 'name_1006', 'namespace_1003',"definition_1006" ,"ensembl_transcript_id"),
                        filters = 'external_gene_name',
                        values = res_Heart_ND$gene,
                        mart = Zebra_finch)

go_annotations_2 <- getBM(attributes = c("external_gene_name",'ensembl_gene_id', "description","entrezgene_id"),
                        filters = 'external_gene_name',
                        values = all_heart_genes,
                        mart = Zebra_finch)


go_annotations_3 <- getBM(attributes = c("external_gene_name",
                                         "description",
                                         'go_id', 
                                         'name_1006',
                                         'namespace_1003',
                                         "definition_1006"),
                          filters = 'external_gene_name',
                          values = all_heart_genes,
                          mart = Zebra_finch)
go_annotations_3 <- go_annotations_3 %>% mutate_all(na_if,"")
na_indices <- which(is.na(go_annotations_3$go_id))
values_col2 <- go_annotations_3$external_gene_name[na_indices]

# view(go_annotations)
# downregHeart_ND %>% data.frame() %>% view()

all_heart_genes <- norm_Heart %>% rownames_to_column(var="gene") %>% pull(gene)

# List of all gene symbols I have just in case there are repitiions
# all_gene_symbols <- unique(res_Heart_ND$gene)

# List of gene symbols retrieved in go_annotations
retrieved_gene_symbols <- unique(go_annotations_squ$external_gene_name)




# Find missing gene symbols
missing_gene_symbols <- setdiff(all_heart_genes, retrieved_gene_symbols)
items_starting_with_LOC <- Filter(function(x) startsWith(x, "LOC"), all_heart_genes)
length(missing_gene_symbols)-length(items_starting_with_LOC) 
print(missing_gene_symbols)
actually_missing_gene_symbols <- setdiff(missing_gene_symbols,items_starting_with_LOC)
sum(downregHeart_ND %in% actually_missing_gene_symbols)


df <- data.frame(my_column = unlist(items_starting_with_L))
write.table(df, file = "clipboard-128", sep = "/t", row.names = FALSE, col.names = FALSE)
# items_starting_with_L

# Check if missing gene symbols exist in the dataset
# if (length(missing_gene_symbols) > 0) {
  missing_gene_info <- getBM(attributes = c('external_gene_name', 'ensembl_gene_id'),
                             filters = 'external_gene_name',
                             values = "GCNA",
                             mart = Zebra_finch)
  
  print(missing_gene_info)
# } else {
#   print("All gene symbols have been found in the dataset.")
# }
 
  
  
  # Convert the list to a dataframe
  df <- data.frame(my_column = unlist(all_heart_genes))
  
  # Add double quotes around each element
  df$my_column <- paste0('"', df$my_column, '"')
  
  # Convert the column to a comma-separated string
  comma_separated_string <- paste(df$my_column, collapse = ",")
  
  # Copy the string to the clipboard
  writeClipboard(comma_separated_string)
  
  items_starting_with_LOC <- Filter(function(x) startsWith(x, "LOC"), all_heart_genes)
  print(items_starting_with_LOC)  
  
# ORA--------- 
  res_Heart_ND_nNA <- filter(res_Heart_ND, padj != "NA" ) # %>% rownames_to_column(var="gene")
  res_ids <- left_join(res_Heart_ND_nNA, go_annotations_2, by=c("gene"="external_gene_name")) 
  
  
  allOE_genes <- as.character(res_ids$gene)
  
  ## Extract significant results
  sigOE <- dplyr::filter(res_Heart_ND, padj < 0.05)
  
  sigOE_genes <- as.character(sigOE$entrezgene_id)
  
  ego <- enrichGO(gene = sigOE_genes, 
                  universe = allOE_genes$ensembl_gene_id,
                  OrgDb = Zebra_finch, 
                  ont = "BP", 
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)
  
  ego <- enricher(gene = sigOE_genes,
                  universe = allOE_genes$ensembl_gene_id,
                  TERM2GENE = gene2go,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)
  
  ego <- enrichGO(gene = sigOE_genes, 
                  universe = all_heart_genes,
                  keyType = "SYMBOL",
                  OrgDb = org.Gg.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)
  
  
  
  term2gene <- data.frame(res_ids$ensembl_gene_id,res_ids$gene)
  enricher(
    gene=sigOE_genes,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe = allOE_genes,
    minGSSize = 0,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    gson = NULL,
    TERM2GENE=term2gene,
    TERM2NAME = NA
  )
  
  heart_DEG_ND <- c(upregLungs_ND,downregLungs_ND)
  heart_DEG_ND <- c(upregHeart_ND,downregHeart_ND)
  items_starting_with_LOC <- Filter(function(x) startsWith(x, "LOC"), heart_DEG_ND)
  heart_DEG_ND <-  heart_DEG_ND[! heart_DEG_ND %in% items_starting_with_LOC ]
  go_annotations_2 <- na.omit(go_annotations_2)
  in_both <- intersect(go_annotations_2$external_gene_name,heart_DEG_ND)
  
  sig <- go_annotations_2 %>% filter(external_gene_name %in% in_both) %>% pull(entrezgene_id)
  ekeg <- enrichKEGG(
    gene=sig,
    organism = "tgu",
    keyType = "kegg",
    pvalueCutoff = 2)
  
  ekeg
  
  sum(is.na(res_ids$ensembl_gene_id))
  
# vIRUS VIRUS VIRUS !!!!!! GENE GON----

# GSEA--------------
  msigdbr_species()
  mm_hallmark_sets <- msigdbr(
    species = "Homo sapiens", # Replace with species name relevant to your data
    category = "H"
  )

  mm_hallmark_sets %>% data.frame() %>% view()
  keytypes(org.Hs.eg.db)
 
  res_Heart_ND <- res_Heart_ND %>% rownames_to_column(var="gene")

  # Let's create a named vector ranked based on the log2 fold change values
  lfc_vector <- res_Heart_ND$log2FoldChange
  names(lfc_vector) <- res_Heart_ND$gene
  lfc_vector <- sort(lfc_vector, decreasing = TRUE)
  head(lfc_vector)
  
  # Set the seed so our results are reproducible:
  set.seed(2020)

  gsea_results <- GSEA(
    geneList = lfc_vector, # Ordered ranked gene list
    minGSSize = 25, # Minimum gene set size
    maxGSSize = 500, # Maximum gene set set
    pvalueCutoff = 1, # p-value cutoff
    eps = 0, # Boundary for calculating the p value
    seed = TRUE, # Set seed to make results reproducible
    pAdjustMethod = "BH", # Benjamini-Hochberg correction
    TERM2GENE = 
    )

  missing_gene_symbols <- setdiff(res_Heart_ND$gene,mm_hallmark_sets$gene_symbol )

# Annotation forge----------
  
  makeOrgPackageFromNCBI(version = "0.1",
                         author = "Some One <so@someplace.org>",
                         maintainer = "Some One <so@someplace.org>",
                         outputDir = ".",
                         tax_id = "59729",
                         genus = "Taeniopygia",
                         species = "guttata")



  makeOrgPackageFromNCBI(version = "0.1", 
                         author = "W. Zac Stephens <zac.stephens@path.utah.edu>",
                         maintainer = "W. Zac Stephens <zac.stephens@path.utah.edu>",
                         outputDir = "../",
                         tax_id = "237561",
                         genus = "Candida", 
                         species =  "albicans")







































































