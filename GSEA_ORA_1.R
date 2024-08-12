
# -------
# remove libraries not needed
library(clusterProfiler)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(here)
library(viridis)
library(enrichR)
library(scales) #show_col(viridis_pal()(20))
library(AnnotationDbi)
library(biomaRt)
library(org.Hs.eg.db)
library(org.Gg.eg.db)
library(PCAtools)
library(msigdbr)
library(ggupset)
library(AnnotationForge)
library(AnnotationHub)
library(UpSetR)

# Read data in------------------
res_Heart_ND <- read.csv(here("res_Heart_ND.csv"))
norm_Heart <- read.csv(here("norm_Heart.csv"))
load(here("lists.RData"))
# remeber to add upreg and downreg lists as databases too.

# all genes and sig genes-----
all_heart_genes <- norm_Heart %>% pull(gene)
# -----
annot = read.csv(file = "ZF_gonad_annGeneNames.csv")
annot_2 = read.csv(file = "ZF_PTR_annGeneNames.csv")


common_genes <- intersect(all_heart_genes,annot$Gene.Name)
not_found <- setdiff(all_heart_genes,annot$Gene.Name)

# lets see how much comes out of biomart
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
human <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
human_attr<- listAttributes(human)
ensembl_anno <- getBM(attributes = c("external_gene_name",
                                           'ensembl_gene_id',
                                     "ensembl_transcript_id",
                                     "transcript_is_canonical",
                                     "ensembl_gene_id_version"),
                            filters = 'external_gene_name',
                            values = all_heart_genes,
                            mart = human)

ensembl_anno_ref <- ensembl_anno %>% 
  filter(!is.na(ensembl_gene_id)) %>% 
  select(-c(ensembl_gene_id_version,ensembl_transcript_id,transcript_is_canonical) ) %>% unique()
  




#  right now focusing Heart ND only--------



# biomart----------

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

# go_annotations_tgu_sig <- getBM(attributes = c("external_gene_name",
#                                            'go_id', 
#                                            'name_1006',
#                                            'namespace_1003',
#                                            "definition_1006"),
#                             filters = 'external_gene_name',
#                             values = c(upregHeart_ND,downregHeart_ND),
#                             mart = Zebra_finch)
# 
# go_annotations_tgu_all <- getBM(attributes = c("external_gene_name",
#                                                'go_id', 
#                                                'name_1006',
#                                                'namespace_1003',
#                                                "definition_1006"),
#                                 filters = 'external_gene_name',
#                                 values = all_heart_genes_10,
#                                 mart = Zebra_finch)
# 
# for (i in 1:10) {
#   
# }
# 
# 
# retrieved_gene_symbols <- unique(go_annotations_tgu_all$external_gene_name)
# missing_gene_symbols <- setdiff(all_heart_genes, retrieved_gene_symbols)


# Determine the number of elements in each part
n <- length(all_heart_genes)
group_size <- ceiling(n / 10)
# Create a factor that repeats 1 to 4, enough times to cover the list length
grouping <- rep(1:10, each = group_size)[1:n]
# Split the list into 4 parts
split_list <- split(all_heart_genes, grouping)
for (i in 1:10) {
  assign(paste0("all_heart_genes_", i), split_list[[i]])
}

# List of gene lists
gene_lists <- list(all_heart_genes_1,
                   all_heart_genes_2,
                   all_heart_genes_3,
                   all_heart_genes_4,
                   all_heart_genes_5,
                   all_heart_genes_6,
                   all_heart_genes_7,
                   all_heart_genes_8,
                   all_heart_genes_9,
                   all_heart_genes_10)

                   
# making a list that can store results as I run the loop
go_annotations_results <- list()
# I had to loop 1:2,3:4 so on because it would not run otherwise.
for (i in 1:10) {
  go_annotations_results[[i]] <- getBM(attributes = c("external_gene_name",
                                                      'go_id', 
                                                      'name_1006',
                                                      'namespace_1003',
                                                      "definition_1006"),
                                       filters = 'external_gene_name',
                                       values = gene_lists[[i]],
                                       mart = Zebra_finch)
  assign(paste0("go_annotations_tgu_all_", i), go_annotations_results[[i]])
}


combined_df <- do.call(rbind, list(go_annotations_tgu_all_1, go_annotations_tgu_all_2, 
                                   go_annotations_tgu_all_3, go_annotations_tgu_all_4, 
                                   go_annotations_tgu_all_5, go_annotations_tgu_all_6, 
                                   go_annotations_tgu_all_7, go_annotations_tgu_all_8, 
                                   go_annotations_tgu_all_9, go_annotations_tgu_all_10))

# 
# gene_names_external_format <- getBM(attributes = "external_gene_name",
#                                      filters = 'external_gene_name',
#                                      values = all_heart_genes,
#                                      mart = Zebra_finch)

# Check for the number of missing genes
# retrieved_gene_symbols <- unique(combined_df$external_gene_name)
# missing_gene_symbols <- setdiff(all_heart_genes, retrieved_gene_symbols)
# length(missing_gene_symbols)
# # 4741 is a lot a lot! ok. list 5 nothing annotated! X0 :=/
# # List 6 also has suspiciously low number of rows
# 
# retrieved_gene_symbols_6 <- unique(go_annotations_tgu_all_6$external_gene_name)
# missing_gene_symbols_6 <- setdiff(all_heart_genes_6, retrieved_gene_symbols)
# length(missing_gene_symbols_6) 
# # damnnn...991 are missing of 1479. why?
# all of list 5 i.e 1479 genes are missing. 2271 from other tables


# go_annotations_tgu_mis <- getBM(attributes = c("external_gene_name",
#                                                'go_id', 
#                                                'name_1006',
#                                                'namespace_1003',
#                                                "definition_1006"),
#                                 filters = 'external_gene_name',
#                                 values = missing_gene_symbols,
#                                 mart = Zebra_finch)

# intresting it annotated FOXP2 and prnp
# kuch to gadbad hein....these 2 genes are 
# there in combined so why is it here in missing_gene_symbols????

# Okayy missing list have intresting problem only 
# 10,042 get external_gene_names :(  so 4739 genes dont get external names. but there are fewer missing genes
# this list of external_gene_names is also the same as retrived.
# 4739 is the diff b/w all_heart_genes and retrived. so why is missing list longer?!
# hhh I see cause missing list also has PRNP prnp wala diff. ok

combined_df <- combined_df %>% mutate_all(~na_if(., ""))
combined_df_ref <- combined_df %>% filter(!is.na(go_id))

length(unique(combined_df_ref$external_gene_name))
retrieved_gene_symbols <- unique(combined_df_ref$external_gene_name)
missing_gene_symbols <- setdiff(all_heart_genes, retrieved_gene_symbols)
length(missing_gene_symbols)


go_annotations_upa_mis <- getBM(attributes = c("external_gene_name",
                                               'go_id', 
                                               'name_1006',
                                               'namespace_1003',
                                               "definition_1006"),
                                filters = 'external_gene_name',
                                values = missing_gene_symbols,
                                mart = Arctic_ground_squ)




# go_annotations_tgu_mis  is actuaclly arctic  ground squ  
#  refine it!
go_annotations_upa_mis<- go_annotations_upa_mis %>% mutate_all(~na_if(., ""))
go_annotations_upa_mis_ref <- go_annotations_upa_mis %>% filter(!is.na(go_id))

retrieved_fm_gene_symbols <- unique(go_annotations_upa_mis_ref$external_gene_name)
missing_fm_gene_symbols <- setdiff(missing_gene_symbols, retrieved_fm_gene_symbols)
length(missing_fm_gene_symbols)

items_starting_with_LOC <- Filter(function(x) startsWith(x, "LOC"), missing_fm_gene_symbols)
length(items_starting_with_LOC)

I_can_fix_missing <- setdiff(missing_fm_gene_symbols,items_starting_with_LOC)

# ok 1022 missing genes. 16069 to begin with.
# A bunch removed becuse they were outliers. 
# 14781 then. of that. 9660 annotated by tgu.
# Of the remaining 5121, 1357 were annotated by upa.
# Of the 3764 I realized LOC is 2742! 
# the remaining 1022 need to somehow be annotated ig.
# marmot annoatated a 100 genes!


marmot <- useDataset("mmmarmota_gene_ensembl", mart = ensembl)
marmot_attr<- listAttributes(marmot)
# view(Arctic_ground_squ)

go_annotations_mmm_mis <- getBM(attributes = c("external_gene_name",
                                               'go_id', 
                                               'name_1006',
                                               'namespace_1003',
                                               "definition_1006"),
                                filters = 'external_gene_name',
                                values = I_can_fix_missing,
                                mart = marmot)

go_annotations_mmm_mis<- go_annotations_mmm_mis %>% mutate_all(~na_if(., ""))
go_annotations_mmm_mis_ref <- go_annotations_mmm_mis %>% filter(!is.na(go_id))

retrieved_fm_gene_symbols_2 <- unique(go_annotations_mmm_mis_ref$external_gene_name)
missing_fm_gene_symbols_2 <- setdiff(I_can_fix_missing, retrieved_fm_gene_symbols_2)
length(missing_fm_gene_symbols_2)

human <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
human_attr<- listAttributes(human)

go_annotations_hsa_mis <- getBM(attributes = c("external_gene_name",
                                               'go_id', 
                                               'name_1006',
                                               'namespace_1003',
                                               "definition_1006"),
                                filters = 'external_gene_name',
                                values = missing_fm_gene_symbols_2,
                                mart = human)

go_annotations_hsa_mis<- go_annotations_hsa_mis %>% mutate_all(~na_if(., ""))
go_annotations_hsa_mis_ref <- go_annotations_hsa_mis %>% filter(!is.na(go_id))

retrieved_fm_gene_symbols_3 <- unique(go_annotations_hsa_mis_ref$external_gene_name)
missing_fm_gene_symbols_3 <- setdiff(missing_fm_gene_symbols_2, retrieved_fm_gene_symbols_3)
length(missing_fm_gene_symbols_3)

combined_df <- do.call(rbind, list(go_annotations_mmm_mis_ref,
                                   go_annotations_upa_mis_ref,
                                   combined_df_ref,
                                   go_annotations_hsa_mis_ref))
# #### turned out to be useless 0 annotations. 
# # lets use mouse 
# mmusculus_gene_ensembl
# 
# 
# mouse <- useDataset("mmusculus_gene_ensembl", mart = ensembl)
# mouse_attr<- listAttributes(mouse)
# # view(Arctic_ground_squ)
# 
# go_annotations_mmu_mis <- getBM(attributes = c("external_gene_name",
#                                                'go_id', 
#                                                'name_1006',
#                                                'namespace_1003',
#                                                "definition_1006"),
#                                 filters = 'external_gene_name',
#                                 values = I_can_fix_missing,
#                                 mart = mouse)
# 
# go_annotations_mmu_mis<- go_annotations_mmu_mis %>% mutate_all(~na_if(., ""))
# go_annotations_mmu_mis_ref <- go_annotations_mmu_mis %>% filter(!is.na(go_id))
# 
# retrieved_fm_gene_symbols_2 <- unique(go_annotations_mmu_mis_ref$external_gene_name)
# missing_fm_gene_symbols_2 <- setdiff(I_can_fix_missing, retrieved_fm_gene_symbols_2)
# length(missing_fm_gene_symbols_2)


# combined_df <- do.call(rbind, list(go_annotations_mmm_mis_ref,
#                                    go_annotations_upa_mis_ref,
#                                    combined_df_ref))


retrieved_gene_symbols <- unique(combined_df$external_gene_name)
missing_gene_symbols <- setdiff(all_heart_genes, retrieved_gene_symbols)
length(missing_gene_symbols)

items_starting_with_LOC <- Filter(function(x) startsWith(x, "LOC"), all_heart_genes)
length(items_starting_with_LOC)
length(missing_gene_symbols)-length(items_starting_with_LOC)

# ok 3664 are missing of which 2742 LOCs I can't do anything.
# But, only 922 remaining so yeah!!

view(combined_df)

length(setdiff(c(upregHeart_ND,downregHeart_ND),retrieved_gene_symbols))

go_annotations_tgu_sig <- getBM(attributes = c("external_gene_name",
                                               'go_id', 
                                               'name_1006',
                                               'namespace_1003',
                                               "definition_1006"),
                                filters = 'external_gene_name',
                                values = c(upregHeart_ND,downregHeart_ND),
                                mart = Zebra_finch)
go_annotations_tgu_sig<- go_annotations_tgu_sig %>% mutate_all(~na_if(., ""))
go_annotations_tgu_sig_ref <- go_annotations_tgu_sig %>% filter(!is.na(go_id))

retrieved_fm_gene_symbols_4 <- unique(go_annotations_tgu_sig_ref$external_gene_name)
missing_fm_gene_symbols_4 <- setdiff(c(upregHeart_ND,downregHeart_ND), retrieved_fm_gene_symbols_4)
length(missing_fm_gene_symbols_4)

go_annotations_mmm_sig <- getBM(attributes = c("external_gene_name",
                                               'go_id', 
                                               'name_1006',
                                               'namespace_1003',
                                               "definition_1006"),
                                filters = 'external_gene_name',
                                values = missing_fm_gene_symbols_4,
                                mart = marmot)

go_annotations_mmm_sig<- go_annotations_mmm_sig %>% mutate_all(~na_if(., ""))
go_annotations_mmm_sig_ref <- go_annotations_mmm_sig %>% filter(!is.na(go_id))

retrieved_fm_gene_symbols_5 <- unique(go_annotations_mmm_sig_ref$external_gene_name)
missing_fm_gene_symbols_5 <- setdiff(missing_fm_gene_symbols_4, retrieved_fm_gene_symbols_5)
length(missing_fm_gene_symbols_5)

go_annotations_hsa_sig <- getBM(attributes = c("external_gene_name",
                                               'go_id', 
                                               'name_1006',
                                               'namespace_1003',
                                               "definition_1006"),
                                filters = 'external_gene_name',
                                values = missing_fm_gene_symbols_5,
                                mart = human)

go_annotations_hsa_sig<- go_annotations_hsa_sig %>% mutate_all(~na_if(., ""))
go_annotations_hsa_sig_ref <- go_annotations_hsa_sig %>% filter(!is.na(go_id))

retrieved_fm_gene_symbols_6 <- unique(go_annotations_hsa_sig_ref$external_gene_name)
missing_fm_gene_symbols_6 <- setdiff(missing_fm_gene_symbols_5, retrieved_fm_gene_symbols_6)
length(missing_fm_gene_symbols_6)



combined_df_sig <- do.call(rbind, list(go_annotations_mmm_sig_ref,
                                       go_annotations_tgu_sig_ref,
                                       go_annotations_hsa_sig_ref))

view(combined_df_sig)

# go_annotations_upa_sig <- getBM(attributes = c("external_gene_name",
#                                                'go_id', 
#                                                'name_1006',
#                                                'namespace_1003',
#                                                "definition_1006"),
#                                 filters = 'external_gene_name',
#                                 values = missing_fm_gene_symbols_3,
#                                 mart = Arctic_ground_squ)
# 
# # go_annotations_upa_sig<- go_annotations_upa_sig %>% mutate_all(~na_if(., ""))
# go_annotations_upa_sig %>%
#   mutate(across(where(is.character), ~ na_if(., "")))
# go_annotations_upa_sig_ref <- go_annotations_upa_sig %>% filter(!is.na(go_id))
# 
# retrieved_fm_gene_symbols_2 <- unique(go_annotations_upa_sig_ref$external_gene_name)
# missing_fm_gene_symbols_4 <- setdiff(missing_fm_gene_symbols_3, retrieved_fm_gene_symbols_2)
# length(missing_fm_gene_symbols_4)

# 
# combined_df_sig <- do.call(rbind, list(go_annotations_mmm_sig_ref,
#                                        go_annotations_tgu_sig_ref))
# 
# 
# view(combined_df_sig)

# 
# as.factor(combined_df$go_id)
# as.factor(combined_df_sig$go_id)
# 
# combined_df$go_id <- factor(combined_df$go_id)
# value_counts <- table(combined_df$go_id)
# 
# combined_df_sig$go_id <- factor(combined_df_sig$go_id)
# value_counts_sig <- table(combined_df_sig$go_id)
# value_counts <- value_counts %>% data.frame()
# value_counts_sig  <- value_counts_sig %>% data.frame()
# 
# hsapiens_gene_ensembl
# 
# 
# human <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
# human_attr<- listAttributes(human)
# 
# go_annotations_hsa_mis <- getBM(attributes = c("external_gene_name",
#                                                'go_id', 
#                                                'name_1006',
#                                                'namespace_1003',
#                                                "definition_1006"),
#                                 filters = 'external_gene_name',
#                                 values = missing_fm_gene_symbols_2,
#                                 mart = human)
# 
# go_annotations_hsa_mis<- go_annotations_hsa_mis %>% mutate_all(~na_if(., ""))
# go_annotations_hsa_mis_ref <- go_annotations_hsa_mis %>% filter(!is.na(go_id))
# 
# retrieved_fm_gene_symbols_4 <- unique(go_annotations_hsa_mis_ref$external_gene_name)
# missing_fm_gene_symbols_4 <- setdiff(missing_fm_gene_symbols_2, retrieved_fm_gene_symbols_4)
# length(missing_fm_gene_symbols_4)
# 
# combined_df <- do.call(rbind, list(go_annotations_mmm_mis_ref,
#                                    go_annotations_upa_mis_ref,
#                                    combined_df_ref,
#                                    go_annotations_hsa_mis_ref))
# 
# 
# 
# go_annotations_hsa_sig <- getBM(attributes = c("external_gene_name",
#                                                'go_id', 
#                                                'name_1006',
#                                                'namespace_1003',
#                                                "definition_1006"),
#                                 filters = 'external_gene_name',
#                                 values = missing_fm_gene_symbols_3,
#                                 mart = human)
# 
# go_annotations_hsa_sig<- go_annotations_hsa_sig %>% mutate_all(~na_if(., ""))
# go_annotations_hsa_sig_ref <- go_annotations_hsa_sig %>% filter(!is.na(go_id))
# 
# retrieved_fm_gene_symbols_5 <- unique(go_annotations_hsa_sig_ref$external_gene_name)
# missing_fm_gene_symbols_5 <- setdiff(missing_fm_gene_symbols_3, retrieved_fm_gene_symbols_5)
# length(missing_fm_gene_symbols_5)
# 
# 
# 
# combined_df_sig <- do.call(rbind, list(go_annotations_mmm_sig_ref,
#                                        go_annotations_tgu_sig_ref,
#                                        go_annotations_hsa_sig_ref))
# 


combined_df$go_id <- factor(combined_df$go_id)
value_counts <- table(combined_df$go_id)

combined_df_sig$go_id <- factor(combined_df_sig$go_id)
value_counts_sig <- table(combined_df_sig$go_id)
value_counts <- value_counts %>% data.frame()
value_counts_sig  <- value_counts_sig %>% data.frame()

merged_df <- merge(value_counts, value_counts_sig, by = "Var1", all.x = TRUE, suffixes = c("_bg", "_sig"))

merged_df$rest <- nrow(combined_df)-merged_df$Freq_bg

num_go_sig <- nrow(combined_df_sig)

merged_df_ref <-  merged_df %>% filter(!is.na(Freq_sig))
   
# dhyper(1,207,146178,598)

merged_df$p_value <- mapply(function(x, m, k, n) {
  phyper(x - 1, m, k - m, n, lower.tail = FALSE)
}, merged_df$count_chosen, merged_df$count_bag, sum(merged_df$count_bag), sum(merged_df$count_chosen))


merged_df$p_value <- apply(merged_df, 1, function(row) {

  dhyper(merged_df$Freq_sig, merged_df$Freq_bg, merged_df$rest, num_go_sig)
})
merged_df_ref$rest
merged_df_ref$p_value <- apply(merged_df_ref, 1, function(row) {
  white_chosen <- as.numeric(row["Freq_sig"])
  white_bag <- as.numeric(row["Freq_bg"])
  black_bag <- as.numeric(row["rest"])
  balls_chosen <- num_go_sig
  
  dhyper(white_chosen, white_bag, black_bag, balls_chosen)
})

# p.bonferroni <- p.adjust(p.values, method = "bonferroni")
merged_df_ref$p_adj <- p.adjust(merged_df_ref$p_value, method = "BH")

# write.csv(merged_df_ref, file = here("pval_go_heart.csv"), row.names = FALSE)
# write.csv(combined_df, file = here("annotated_heart_nd.csv"), row.names = FALSE)
# write.csv(combined_df_ref, file = here("annotated_heart_nd_sig.csv"), row.names = FALSE)

go_id_sig <- merged_df_ref %>% filter(p_adj<0.05) %>% pull(Var1)

sig_go_fn <- combined_df_sig[combined_df_sig$go_id %in% go_id_sig, ] %>% 
  subset(select = -external_gene_name) %>%  unique()
view(sig_go_fn)
write.csv(sig_go_fn, file = here("go_functions_sig_heart.csv"), row.names = FALSE)





# enrichKEGG------



ekeg <- enrichKEGG(
  gene=c(upregHeart_ND,downregHeart_ND),
  organism = "cpea",
  universe=all_heart_genes,
  keyType = "kegg",
  pvalueCutoff = 0.1)

# ----
hub <- AnnotationHub()
query(hub, c("clostridium","orgdb"))

# enricher-----
combined_df <- read.csv(here("annotated_heart_nd.csv"))
combined_df_sig <- read.csv(here("annotated_heart_nd_sig.csv"))

term2gene <- data.frame(combined_df$go_id ,combined_df$external_gene_name)
goi <- combined_df_sig$external_gene_name
x <- enricher(goi, TERM2GENE = term2gene)
go_id_sig <- x$ID
combined_df_sig[go_id == go_id_sig,]

go_sig_info <- combined_df_sig %>% filter(go_id %in% go_id_sig ) %>%
  select(-external_gene_name)%>%  unique()



