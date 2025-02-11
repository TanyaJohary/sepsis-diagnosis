
# GSE9960

# Retrive geo datas as list
geo_data <- getGEO("GSE9960", GSEMatrix = TRUE)

# Retrieve the expression data
expression_data <- exprs(geo_data[[1]]) # data is already normalized by RMA

# Retrieve the phenotype data
pheno_data <- pData(geo_data[[1]])
colnames(pheno_data)
sub_pheno <- pheno_data[, c("geo_accession", "characteristics_ch1")]
write.csv(sub_pheno, "sub_phenoGSE9960.csv", row.names = TRUE)

# Retrieve the feature data for gene symbols
feature_data <- fData(geo_data[[1]])
colnames(feature_data)
sub_feature <- feature_data[, c("ID", "Gene Symbol")]
write.csv(sub_feature, "sub_featureGSE9960.csv", row.names = TRUE)


# Merge the expression data and sub_feature data
expression_with_symbols <- merge(expression_data, sub_feature, by.x = "row.names", by.y = "ID", all.x = TRUE)
write.csv(expression_with_symbols, "expression_with_symbolsGSE9960.csv", row.names = TRUE)

# search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)


present_genes <- genes_of_interest[genes_of_interest %in% expression_with_symbols$`Gene Symbol`]
missing_genes <- genes_of_interest[!genes_of_interest %in% expression_with_symbols$`Gene Symbol`]

# save the present and missing genes for one sample
write.csv(present_genes, "present_genesGSE9960.csv", row.names = TRUE, quote = FALSE)
write.csv(missing_genes, "missing_genesGSE9960.csv", row.names = TRUE, quote = FALSE)


# Filter out the expression of genes of interest
library(dplyr)
filtered_gene_expression <- expression_with_symbols %>%
  filter(`Gene Symbol` %in% genes_of_interest)

# Remove the column contains the probbe ids
filtered_gene_expression <- filtered_gene_expression[,-1]

# make unique value for each probe 
filtered_gene_expression$`Gene Symbol` <- make.unique(as.character(filtered_gene_expression$`Gene Symbol`))

# Transpose the data to make header by genes names
transposed_data <- filtered_gene_expression %>%
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = `Gene Symbol`,
    values_from = Value
  )

# Merge with pheno data to iddentify the sepsis patients and controls
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Rename the target column as Label
sepsis_data <- merged_transposed %>%
  rename_with(~"Label", matches("characteristics"))

write.csv(sepsis_data, "sepsis_dataGSE9960.csv", row.names = FALSE)

################################################################################
################################# 60 Genes #####################################
################################################################################

genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A",
  "LY96", "CXCR4", "TRIL", "TLR5", "CXCL13"
)

# Filter out the expression of genes of interest
library(dplyr)
filtered_gene_expression <- expression_with_symbols %>%
  filter(`Gene Symbol` %in% genes_of_interest)

# Remove the column contains the probbe ids
filtered_gene_expression <- filtered_gene_expression[,-1]

# make unique value for each probe 
filtered_gene_expression$`Gene Symbol` <- make.unique(as.character(filtered_gene_expression$`Gene Symbol`))

# Transpose the data to make header by genes names
transposed_data <- filtered_gene_expression %>%
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = `Gene Symbol`,
    values_from = Value
  )

# Merge with pheno data to iddentify the sepsis patients and controls
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Rename the target column as Label
sepsis_data <- merged_transposed %>%
  rename_with(~"Label", matches("characteristics"))

write.csv(sepsis_data, "sepsis_data_60_GSE9960.csv", row.names = FALSE)

################################################################################
################################# Selected Genes #####################################
################################################################################

genes_of_interest <- c(
  "CALCA", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNB1", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A",
  "LY96", "CXCR4", "TLR5"
)

# Filter out the expression of genes of interest
library(dplyr)
filtered_gene_expression <- expression_with_symbols %>%
  filter(`Gene Symbol` %in% genes_of_interest)

# Remove the column contains the probbe ids
filtered_gene_expression <- filtered_gene_expression[,-1]

# make unique value for each probe 
filtered_gene_expression$`Gene Symbol` <- make.unique(as.character(filtered_gene_expression$`Gene Symbol`))

# Transpose the data to make header by genes names
transposed_data <- filtered_gene_expression %>%
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = `Gene Symbol`,
    values_from = Value
  )

# Merge with pheno data to iddentify the sepsis patients and controls
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Rename the target column as Label
sepsis_data <- merged_transposed %>%
  rename_with(~"Label", matches("characteristics"))

write.csv(sepsis_data, "sepsis_data_selected_GSE9960.csv", row.names = FALSE)

