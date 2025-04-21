
# GSE57065
geo_data <- getGEO("GSE57065", GSEMatrix = TRUE)

# Retrieve the expression data
expression_data <- exprs(geo_data[[1]]) # already normalized

# Retrieve the phenotype data
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c(2, 35, 37)]
sub_pheno <- sub_pheno %>% 
  rename_with(~"Label", matches("sapsii"))
sub_pheno <- sub_pheno %>% 
  mutate(Label = case_when(
    Label == "NA" ~ "Healthy",
    Label == "SAPSII-High" ~ "Septic",
    Label == "SAPSII-Low" ~ "Septic",
    TRUE ~ Label
  ))
write.csv(sub_pheno, "sub_phenoGSE57065.csv", row.names = FALSE)

# create PCAs
pca_data <- prcomp(t(expression_data), scale. = TRUE)

# extract PCAs
pca_df <- as.data.frame(pca_data$x)
pca_df$Sample <- rownames(pca_df)
pca_df$batch <- pheno_data$`collection time:ch1`

ggplot(pca_df, aes(x = PC1, y = PC2, colour = batch))+
  geom_point(size = 3)+
  labs(title = "PCA:",
       x = "PC1",
       y = "PC1")+
  theme_minimal() #  there is clustering between healthy and other patients, but it can be biological differences

# Retrieve the feature data to map the genes symbols
feature_data <- fData(geo_data[[1]])
sub_feature <- feature_data[, c(1, 11)]

# Merge expression data to sub feature data to map gene symbols for each probe
expression_with_symbols <- merge(expression_data, sub_feature, by.x = "row.names", by.y = "ID", all.x = TRUE)
write.csv(expression_with_symbols, "expression_with_symbolsGSE75065.csv", row.names = FALSE)

# genes of interest names
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

# Search for our genes of interest
present_genes <- genes_of_interest[genes_of_interest %in% expression_with_symbols$`Gene Symbol`]
missing_genes <- genes_of_interest[!genes_of_interest %in% expression_with_symbols$`Gene Symbol`]

write.csv(present_genes, "present_genesGSE57065.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE57065.csv", row.names = FALSE)

# Filter our genes of interest
filtered_genes <- expression_with_symbols %>%
  filter(`Gene Symbol` %in% present_genes)

# Make unique value for different probe of each gene
filtered_genes$`Gene Symbol` <- make.unique(as.character(filtered_genes$`Gene Symbol`))

# Remove the probe id column
filtered_genes <- filtered_genes[, -1]

# Transpose teh data to make the gene symbols as header
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = `Gene Symbol`,
    values_from = Value
  )

# Merge the transposed data to sub pheno data to identify the septic shock and healthy ones
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

################# filter patients of sampling in 0hour
sepsis_hr0 <- merged_transposed %>% 
  filter(Label == "Healthy" | (Label == "Septic" & `collection time:ch1` == "0 hr" ) )
sepsis_hr0 <- sepsis_hr0[, -101]

write.csv(sepsis_hr0, "sepsis_hr0_GSE57065.csv", row.names = FALSE)

################# filter patients of sampling in 24hour
sepsis_hr24 <- merged_transposed %>% 
  filter(Label == "Healthy" | (Label == "Septic" & `collection time:ch1` == "24 hr" ) )
sepsis_hr24 <- sepsis_hr24[, -101]

write.csv(sepsis_hr24, "sepsis_hr24_GSE57065.csv", row.names = FALSE)

################# filter patients of sampling in 48hour
sepsis_hr48 <- merged_transposed %>% 
  filter(Label == "Healthy" | (Label == "Septic" & `collection time:ch1` == "48 hr" ) )
sepsis_hr48 <- sepsis_hr48[, -101]

write.csv(sepsis_hr48, "sepsis_hr48_GSE57065.csv", row.names = FALSE)

################################################################################
############################# Selected Genes ###################################
################################################################################
selected_genes <- c( "TREM1", "PLAUR", "IL6", "HMGB1", "ELANE", "CXCL8", "IL10", "IL1B",
                     "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
                     "C5AR1", "LCN2", "CXCL10", "IFNG", "CCL25", "CX3CR1", "P2RX7",
                     "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88",
                     "NLRP3", "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1",
                     "IL1R2", "CD177", "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3",
                     "GATA3", "CCR2", "HIF1A", "LY96", "CXCR4", "CXCL12" )


# Search for our genes of interest
present_genes <- selected_genes[selected_genes %in% expression_with_symbols$`Gene Symbol`]
missing_genes <- selected_genes[!selected_genes %in% expression_with_symbols$`Gene Symbol`]

write.csv(present_genes, "present_genesGSE57065.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE57065.csv", row.names = FALSE)

# Filter our genes of interest
filtered_genes <- expression_with_symbols %>%
  filter(`Gene Symbol` %in% present_genes)

# Make unique value for different probe of each gene
filtered_genes$`Gene Symbol` <- make.unique(as.character(filtered_genes$`Gene Symbol`))

# Remove the probe id column
filtered_genes <- filtered_genes[, -1]

# Transpose teh data to make the gene symbols as header
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = `Gene Symbol`,
    values_from = Value
  )

# Merge the transposed data to sub pheno data to identify the septic shock and healthy ones
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

################# filter patients of sampling in 0hour
sepsis_hr0 <- merged_transposed %>% 
  filter(Label == "Healthy" | (Label == "Septic" & `collection time:ch1` == "0 hr" ) )
sepsis_hr0 <- sepsis_hr0[, -94]

write.csv(sepsis_hr0, "sepsis_hr0selected_GSE57065.csv", row.names = FALSE)

################# filter patients of sampling in 24hour
sepsis_hr24 <- merged_transposed %>% 
  filter(Label == "Healthy" | (Label == "Septic" & `collection time:ch1` == "24 hr" ) )
sepsis_hr24 <- sepsis_hr24[, -94]

write.csv(sepsis_hr24, "sepsis_hr24selected_GSE57065.csv", row.names = FALSE)

################# filter patients of sampling in 48hour
sepsis_hr48 <- merged_transposed %>% 
  filter(Label == "Healthy" | (Label == "Septic" & `collection time:ch1` == "48 hr" ) )
sepsis_hr48 <- sepsis_hr48[, -94]

write.csv(sepsis_hr48, "sepsis_hr48selected_GSE57065.csv", row.names = FALSE)

###############################################################################
############################################################### all top genes

top_genes <- c("CD177", "ITGAM", "ARG1", "C3AR1", "FCGR1A",
               "S100A8", "S100A9", "S100A12", "HMGB1",
               "GATA3", "SOCS3", "MAPK14", "IL10", "BCL2",
               "IL1R2", "MYD88", "HIF1A")

# Search for our genes of interest
present_genes <- top_genes[top_genes %in% expression_with_symbols$`Gene Symbol`]
missing_genes <- top_genes[!top_genes %in% expression_with_symbols$`Gene Symbol`]

# Filter our genes of interest
filtered_genes <- expression_with_symbols %>%
  filter(`Gene Symbol` %in% present_genes)

# Make unique value for different probe of each gene
filtered_genes$`Gene Symbol` <- make.unique(as.character(filtered_genes$`Gene Symbol`))

# Remove the probe id column
filtered_genes <- filtered_genes[, -1]

# Transpose teh data to make the gene symbols as header
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = `Gene Symbol`,
    values_from = Value
  )

# Merge the transposed data to sub pheno data to identify the septic shock and healthy ones
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

sepsis_hr48 <- merged_transposed %>% 
  filter(Label == "Healthy" | (Label == "Septic" & `collection time:ch1` == "48 hr" ) )
sepsis_hr48 <- sepsis_hr48[, -39]

write.csv(sepsis_hr48, "sepsis_hr48alltop_GSE57065.csv", row.names = FALSE)

################################################################################
############################################################# top cluster 1

cluster1 <- c("CD177", "ITGAM", "C3AR1", "FCGR1A")

# Search for our genes of interest
present_genes <- cluster1[cluster1 %in% expression_with_symbols$`Gene Symbol`]
missing_genes <- cluster1[!cluster1 %in% expression_with_symbols$`Gene Symbol`]

# Filter our genes of interest
filtered_genes <- expression_with_symbols %>%
  filter(`Gene Symbol` %in% present_genes)

# Make unique value for different probe of each gene
filtered_genes$`Gene Symbol` <- make.unique(as.character(filtered_genes$`Gene Symbol`))

# Remove the probe id column
filtered_genes <- filtered_genes[, -1]

# Transpose teh data to make the gene symbols as header
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = `Gene Symbol`,
    values_from = Value
  )

# Merge the transposed data to sub pheno data to identify the septic shock and healthy ones
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

sepsis_hr48 <- merged_transposed %>% 
  filter(Label == "Healthy" | (Label == "Septic" & `collection time:ch1` == "48 hr" ) )
sepsis_hr48 <- sepsis_hr48[, -6]

write.csv(sepsis_hr48, "sepsis_hr48cluster1_GSE57065.csv", row.names = FALSE)

################################################################################
############################################################# top cluster 2

cluster2 <- c("S100A8", "S100A9", "S100A12", "HMGB1")

# Search for our genes of interest
present_genes <- cluster2[cluster2 %in% expression_with_symbols$`Gene Symbol`]
missing_genes <- cluster2[!cluster2 %in% expression_with_symbols$`Gene Symbol`]

# Filter our genes of interest
filtered_genes <- expression_with_symbols %>%
  filter(`Gene Symbol` %in% present_genes)

# Make unique value for different probe of each gene
filtered_genes$`Gene Symbol` <- make.unique(as.character(filtered_genes$`Gene Symbol`))

# Remove the probe id column
filtered_genes <- filtered_genes[, -1]

# Transpose teh data to make the gene symbols as header
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = `Gene Symbol`,
    values_from = Value
  )

# Merge the transposed data to sub pheno data to identify the septic shock and healthy ones
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

sepsis_hr48 <- merged_transposed %>% 
  filter(Label == "Healthy" | (Label == "Septic" & `collection time:ch1` == "48 hr" ) )
sepsis_hr48 <- sepsis_hr48[, -11]

write.csv(sepsis_hr48, "sepsis_hr48cluster2_GSE57065.csv", row.names = FALSE)

################################################################################
###################################################################. cluster3
cluster3 <- c("GATA3", "SOCS3", "MAPK14", "IL10", "BCL2",
              "IL1R2", "MYD88", "HIF1A", "ARG1")

# Search for our genes of interest
present_genes <- cluster3[cluster3 %in% expression_with_symbols$`Gene Symbol`]
missing_genes <- cluster3[!cluster3 %in% expression_with_symbols$`Gene Symbol`]

# Filter our genes of interest
filtered_genes <- expression_with_symbols %>%
  filter(`Gene Symbol` %in% present_genes)

# Make unique value for different probe of each gene
filtered_genes$`Gene Symbol` <- make.unique(as.character(filtered_genes$`Gene Symbol`))

# Remove the probe id column
filtered_genes <- filtered_genes[, -1]

# Transpose teh data to make the gene symbols as header
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = `Gene Symbol`,
    values_from = Value
  )

# Merge the transposed data to sub pheno data to identify the septic shock and healthy ones
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

sepsis_hr48 <- merged_transposed %>% 
  filter(Label == "Healthy" | (Label == "Septic" & `collection time:ch1` == "48 hr" ) )
sepsis_hr48 <- sepsis_hr48[, -26]

write.csv(sepsis_hr48, "sepsis_hr48cluster3_GSE57065.csv", row.names = FALSE)

################################################################################
################################################################# last top genes


last_cluster <- c("CD177", "S100A8", "S100A9", "GATA3", "S100A12", "ARG1","ITGAM",
                  "SOCS3", "C3AR1", "IL1R2", "IL10", "FCGR1A", "MMP8", "MAPK14")

# Search for our genes of interest
present_genes <- last_cluster[last_cluster %in% expression_with_symbols$`Gene Symbol`]
missing_genes <- last_cluster[!last_cluster %in% expression_with_symbols$`Gene Symbol`]

# Filter our genes of interest
filtered_genes <- expression_with_symbols %>%
  filter(`Gene Symbol` %in% present_genes)

# Make unique value for different probe of each gene
filtered_genes$`Gene Symbol` <- make.unique(as.character(filtered_genes$`Gene Symbol`))

# Remove the probe id column
filtered_genes <- filtered_genes[, -1]

# Transpose teh data to make the gene symbols as header
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = `Gene Symbol`,
    values_from = Value
  )

# Merge the transposed data to sub pheno data to identify the septic shock and healthy ones
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

sepsis_hr0 <- merged_transposed %>% 
  filter(Label == "Healthy" | (Label == "Septic" & `collection time:ch1` == "0 hr" ) )
sepsis_hr0 <- sepsis_hr0[, -30]

write.csv(sepsis_hr0, "sepsis_hr0lasttop_GSE57065.csv", row.names = FALSE)

