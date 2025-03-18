
# read the series matrix file
seriesmatrix <- read.delim("GSE28750_series_matrix.txt", header = TRUE, stringsAsFactors = FALSE, skip = 76)
write.csv(seriesmatrix, "seriesmatrixGSE28750.csv", row.names = TRUE)

# Retrieve the geo dataset as list
geo_data <- getGEO("GSE28750", GSEMatrix = TRUE)

# Retrieve the expression data
expression_data <- exprs(geo_data[[1]])

# Phenotypic Data
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c(2,35)]


# feature data
feature_data <- fData(geo_data[[1]])
sub_feature <- feature_data[, c(1, 11)]

# merge the expression data with gene symbols
expression_with_symbols <- merge(expression_data, sub_feature, by.x = "row.names", by.y = "ID", all.x = TRUE )
write.csv(expression_with_symbols, "expression_with_symbolGSE28750.csv", row.names = TRUE)

# genes of interest names
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

# find the present and missing genes
present_genes <- genes_of_interest[genes_of_interest %in% expression_with_symbols$`Gene Symbol`]
write.csv(present_genes, "present_genesGSE28750.csv", row.names = TRUE)

missing_genes <- genes_of_interest[!genes_of_interest %in% expression_with_symbols$`Gene Symbol`]
write.csv(missing_genes, "missing_genesGSE28750.csv", row.names = TRUE)


###########################################################################
####################### prepare data for Random forest#####################
###########################################################################


# filter out our genes of interest from expression data
filtered_genes <- expression_with_symbols %>%
   filter(`Gene Symbol` %in% genes_of_interest)

filtered_genes <- filtered_genes[, -1]

# Make unique of all probe ids for each gene
filtered_genes$`Gene Symbol` <- make.unique(as.character(filtered_genes$`Gene Symbol`))

# Transpose the data
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

# Merge transposed data to sub pheno data to identify the sepsis patients
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Exclude the post surgical samples form data and keep only SEPSIS and HEALTHY
sepsis_data <- merged_transposed %>%
  filter(`health status:ch1` != "POST_SURGICAL")

# Rename the target colomn's name
sepsis_data <- sepsis_data %>%
  rename_with(~"Label", matches("health"))
# Save the ready dataset for RF
write.csv(sepsis_data, "sepsis_data_GSE28750.csv", row.names = FALSE)

################################################################################
################################## 60 Genes ####################################

genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A",
  "LY96", "CXCR4", "TRIL", "TLR5", "CXCL13"
)

# filter out our genes of interest from expression data
filtered_genes_60 <- expression_with_symbols %>%
  filter(`Gene Symbol` %in% genes_of_interest)

filtered_genes_60 <- filtered_genes_60[, -1]

# Make unique of all probe ids for each gene
filtered_genes_60$`Gene Symbol` <- make.unique(as.character(filtered_genes_60$`Gene Symbol`))

# Transpose the data
transposed_data_60 <- filtered_genes_60 %>%
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = `Gene Symbol`,
    values_from = Value
  )

# Merge transposed data to sub pheno data to identify the sepsis patients
merged_transposed_60 <- merge(transposed_data_60, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Exclude the post surgical samples form data and keep only SEPSIS and HEALTHY
sepsis_data_60 <- merged_transposed_60 %>%
  filter(`health status:ch1` != "POST_SURGICAL")

# Rename the target colomn's name
sepsis_data_60 <- sepsis_data_60 %>%
  rename_with(~"Label", matches("health"))

# Save the ready dataset for RF
write.csv(sepsis_data_60, "sepsis_data_60genes_GSE28750.csv", row.names = FALSE)

##################################################################### last cluster


last_cluster <- c("CD177", "S100A8", "S100A9", "GATA3", "S100A12", "ARG1","ITGAM",
                  "SOCS3", "C3AR1", "IL1R2", "IL10", "FCGR1A", "MMP8", "MAPK14")

# filter out our genes of interest from expression data
filtered_genes <- expression_with_symbols %>%
  filter(`Gene Symbol` %in% last_cluster)

filtered_genes <- filtered_genes[, -1]

# Make unique of all probe ids for each gene
filtered_genes$`Gene Symbol` <- make.unique(as.character(filtered_genes$`Gene Symbol`))

# Transpose the data
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

# Merge transposed data to sub pheno data to identify the sepsis patients
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Exclude the post surgical samples form data and keep only SEPSIS and HEALTHY
sepsis_data <- merged_transposed %>%
  filter(`health status:ch1` != "POST_SURGICAL")

# Rename the target colomn's name
sepsis_data <- sepsis_data %>%
  rename_with(~"Label", matches("health"))

# Save the ready dataset for RF
write.csv(sepsis_data, "sepsis_lasttop_GSE28750.csv", row.names = FALSE)
