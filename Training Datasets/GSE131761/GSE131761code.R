
# Retrieve the geo data as alist
geo_data <- getGEO("GSE131761", GSEMatrix = TRUE)

# REtrieve the expression data
expression_data <- exprs(geo_data[[1]])

# Retrieve the feature data
feature_data <- fData(geo_data[[1]])
sub_feature <- feature_data[, c(1, 7)]
 
# Retrieve the phenotype data
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c(2, 37)]
sub_pheno <- sub_pheno %>% 
  rename_with(~"Label", matches("diagnosis"))

# investigate the batch effect
library(ggplot2)
pca_data <- prcomp(t(expression_data))  # Transpose: genes as columns
pca_df <- data.frame(
  PC1 = pca_data$x[,1],
  PC2 = pca_data$x[,2],
  Batch = pheno_data$`age:ch1`  # Adjust to match the batch column in our pheno_data
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA of Expression Data Before Batch Correction") +
  labs(color = "Batch")

summary(aov(PC1 ~ pheno_data$`age:ch1`, data = pca_df)) # more than 0.005
summary(aov(PC2 ~ pheno_data$`age:ch1`, data = pca_df)) # more than 0.005 no need for batch correction


# Merge expression data an fub feature
expression_with_symbols <- merge(expression_data, sub_feature, by.x = "row.names", by.y = "ID", all.x = TRUE)
write.csv(expression_with_symbols, "expression_with_dataGSE131761.csv", row.names = FALSE)

# our genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "IL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

# Search for our genes of interest
present_genes <- genes_of_interest[genes_of_interest %in% expression_with_symbols$GENE_SYMBOL]
missing_genes <- genes_of_interest[!genes_of_interest %in% expression_with_symbols$GENE_SYMBOL] # "FCGR1A" "IFNA1" 

write.csv(present_genes, "present_genesGSE131761.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE131761.csv", row.names = FALSE)

# Filter our genes of interest
filtered_genes <- expression_with_symbols %>% 
  filter(GENE_SYMBOL %in% genes_of_interest)

# Make unique value of each probe id
filtered_genes$GENE_SYMBOL <- make.unique(as.character(filtered_genes$GENE_SYMBOL))

# Remove the probe ids columns
filtered_genes <- filtered_genes[, -1]

# Transpose the data
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -GENE_SYMBOL,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = GENE_SYMBOL,
    values_from = Value
  )

# Merge to phenotype data to filter setic patients and control ones
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# filter septic shock and control patients
septic_data <- merged_data %>% 
  filter(Label %in% c ("control patient", "septic shock"))

write.csv(septic_data, "septic_dataGSE131761.csv", row.names = FALSE)


#################################################### Selected genes

# our genes of interest
genes_of_interest <- c(
  "CALCA", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "IL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNB1", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A",
  "LY96", "CXCR4", "TLR5"
)
# Search for our genes of interest
present_genes <- genes_of_interest[genes_of_interest %in% expression_with_symbols$GENE_SYMBOL]
missing_genes <- genes_of_interest[!genes_of_interest %in% expression_with_symbols$GENE_SYMBOL] # "FCGR1A" "IFNA1" 

write.csv(present_genes, "present_genesGSE131761.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE131761.csv", row.names = FALSE)

# Filter our genes of interest
filtered_genes <- expression_with_symbols %>% 
  filter(GENE_SYMBOL %in% genes_of_interest)

# Make unique value of each probe id
filtered_genes$GENE_SYMBOL <- make.unique(as.character(filtered_genes$GENE_SYMBOL))

# Remove the probe ids columns
filtered_genes <- filtered_genes[, -1]

# Transpose the data
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -GENE_SYMBOL,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = GENE_SYMBOL,
    values_from = Value
  )

# Merge to phenotype data to filter setic patients and control ones
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# filter septic shock and control patients
septic_data <- merged_data %>% 
  filter(Label %in% c ("control patient", "septic shock"))

write.csv(septic_data, "septic_selectedGSE131761.csv", row.names = FALSE)

################################################################## all top genes
top_genes <- c("CD177", "ITGAM", "ARG1", "C3AR1", "FCGR1A",
               "S100A8", "S100A9", "S100A12", "HMGB1",
               "GATA3", "SOCS3", "MAPK14", "IL10", "BCL2",
               "IL1R2", "MYD88", "HIF1A")

# Search for our genes of interest
present_genes <- top_genes[top_genes %in% expression_with_symbols$GENE_SYMBOL]
missing_genes <- top_genes[!top_genes %in% expression_with_symbols$GENE_SYMBOL] # "FCGR1A"  


# Filter our genes of interest
filtered_genes <- expression_with_symbols %>% 
  filter(GENE_SYMBOL %in% top_genes)

# Make unique value of each probe id
filtered_genes$GENE_SYMBOL <- make.unique(as.character(filtered_genes$GENE_SYMBOL))

# Remove the probe ids columns
filtered_genes <- filtered_genes[, -1]

# Transpose the data
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -GENE_SYMBOL,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = GENE_SYMBOL,
    values_from = Value
  )

# Merge to phenotype data to filter setic patients and control ones
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# filter septic shock and control patients
septic_data <- merged_data %>% 
  filter(Label %in% c ("control patient", "septic shock"))

write.csv(septic_data, "septic_alltopGSE131761.csv", row.names = FALSE)

############################################################# top cluster 1

cluster1 <- c("CD177", "ITGAM", "C3AR1", "FCGR1A")

# Search for our genes of interest
present_genes <- cluster1[cluster1 %in% expression_with_symbols$GENE_SYMBOL]
missing_genes <- cluster1[!cluster1 %in% expression_with_symbols$GENE_SYMBOL] # "FCGR1A"  


# Filter our genes of interest
filtered_genes <- expression_with_symbols %>% 
  filter(GENE_SYMBOL %in% cluster1)

# Make unique value of each probe id
filtered_genes$GENE_SYMBOL <- make.unique(as.character(filtered_genes$GENE_SYMBOL))

# Remove the probe ids columns
filtered_genes <- filtered_genes[, -1]

# Transpose the data
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -GENE_SYMBOL,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = GENE_SYMBOL,
    values_from = Value
  )

# Merge to phenotype data to filter setic patients and control ones
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# filter septic shock and control patients
septic_data <- merged_data %>% 
  filter(Label %in% c ("control patient", "septic shock"))

write.csv(septic_data, "septic_cluster1_GSE131761.csv", row.names = FALSE)

################################################################# cluster2
cluster2 <- c("S100A8", "S100A9", "S100A12", "HMGB1")

# Search for our genes of interest
present_genes <- cluster2[cluster2 %in% expression_with_symbols$GENE_SYMBOL]
missing_genes <- cluster2[!cluster2 %in% expression_with_symbols$GENE_SYMBOL]  


# Filter our genes of interest
filtered_genes <- expression_with_symbols %>% 
  filter(GENE_SYMBOL %in% cluster2)

# Make unique value of each probe id
filtered_genes$GENE_SYMBOL <- make.unique(as.character(filtered_genes$GENE_SYMBOL))

# Remove the probe ids columns
filtered_genes <- filtered_genes[, -1]

# Transpose the data
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -GENE_SYMBOL,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = GENE_SYMBOL,
    values_from = Value
  )

# Merge to phenotype data to filter setic patients and control ones
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# filter septic shock and control patients
septic_data <- merged_data %>% 
  filter(Label %in% c ("control patient", "septic shock"))

write.csv(septic_data, "septic_cluster2_GSE131761.csv", row.names = FALSE)

################################################################ cluster 3

cluster3 <- c("GATA3", "SOCS3", "MAPK14", "IL10", "BCL2",
              "IL1R2", "MYD88", "HIF1A", "ARG1")
# Search for our genes of interest
present_genes <- cluster3[cluster3 %in% expression_with_symbols$GENE_SYMBOL]
missing_genes <- cluster3[!cluster3 %in% expression_with_symbols$GENE_SYMBOL]  


# Filter our genes of interest
filtered_genes <- expression_with_symbols %>% 
  filter(GENE_SYMBOL %in% cluster3)

# Make unique value of each probe id
filtered_genes$GENE_SYMBOL <- make.unique(as.character(filtered_genes$GENE_SYMBOL))

# Remove the probe ids columns
filtered_genes <- filtered_genes[, -1]

# Transpose the data
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -GENE_SYMBOL,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = GENE_SYMBOL,
    values_from = Value
  )

# Merge to phenotype data to filter setic patients and control ones
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# filter septic shock and control patients
septic_data <- merged_data %>% 
  filter(Label %in% c ("control patient", "septic shock"))

write.csv(septic_data, "septic_cluster3_GSE131761.csv", row.names = FALSE)

################################################################################
############################################################### last top genes
last_cluster <- c("CD177", "S100A8", "S100A9", "GATA3", "S100A12", "ARG1","ITGAM",
                  "SOCS3", "C3AR1", "IL1R2", "IL10", "FCGR1A", "MMP8", "MAPK14")

# Search for our genes of interest
present_genes <- last_cluster[last_cluster %in% expression_with_symbols$GENE_SYMBOL]
missing_genes <- last_cluster[!last_cluster %in% expression_with_symbols$GENE_SYMBOL]  


# Filter our genes of interest
filtered_genes <- expression_with_symbols %>% 
  filter(GENE_SYMBOL %in% present_genes)

# Make unique value of each probe id
filtered_genes$GENE_SYMBOL <- make.unique(as.character(filtered_genes$GENE_SYMBOL))

# Remove the probe ids columns
filtered_genes <- filtered_genes[, -1]

# Transpose the data
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -GENE_SYMBOL,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = GENE_SYMBOL,
    values_from = Value
  )

# Merge to phenotype data to filter setic patients and control ones
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# filter septic shock and control patients
septic_data <- merged_data %>% 
  filter(Label %in% c ("control patient", "septic shock"))

write.csv(septic_data, "septic_lasttop_GSE131761.csv", row.names = FALSE)
