
# Retrieve the dataset as list
geo_data <- getGEO("GSE13015", GSEMatrix = TRUE) # dataset has 2 platforms GSE13015-GPL6106 and GSE13015-GPL6947

# Retrieve all data related to first GPL 6106
expression_6106 <-exprs(geo_data[[1]])
feature_6106 <- fData(geo_data[[1]])
pheno_6106 <- pData(geo_data[[1]])


# expression data is not well normalized, we want to make sure
library(limma)
normalized_6106 <- normalizeBetweenArrays(expression_6106, method = "quantile")
boxplot(normalized_6106, main = "Boxplot of Quantile Normalized Data")

# Merge the expression data with fetaure data to map the genes symbols
sub_feature_6106 <- feature_6106[, c(1, 7)]
merged_6106 <- merge(normalized_6106, sub_feature_6106, by.x = "row.names", by.y = "ID", all.x = TRUE)
write.csv(merged_6106, "merged_6106.csv", row.names = FALSE)

# Search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELA2", "IL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "CIAS1",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

present_genes <- genes_of_interest[genes_of_interest %in% merged_6106$Symbol] # gene NLRP3 exist as "CIAS1", ELANE as "ELA2"
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_6106$Symbol]

# Save the present and missing genes
write.csv(present_genes, "present_genesGSE13015-6106.csv", row.names = TRUE)
write.csv(missing_genes, "missing_genesGSE13015-6106.csv", row.names = TRUE)

# Filter the genes of interest
filtered_gene_6106 <- merged_6106 %>%
  filter(Symbol %in% genes_of_interest)

# Remove the column probe id and make unique the symbol column to transpose
filtered_gene_6106 <- filtered_gene_6106[,-1]
filtered_gene_6106$Symbol <- make.unique(as.character(filtered_gene_6106$Symbol))

# Transpose the data and make genes symbol as colnames
transposed_6106 <- filtered_gene_6106 %>%
  pivot_longer(
    cols = -Symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Symbol,
    values_from = Value
  )

# Merge the transposed data with pheno data to identift sepsis and healthy ones
sub_pheno_6106 <- pheno_6106[, c(2, 40)]
transposed_merged_6106 <- merge(transposed_6106, sub_pheno_6106, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Define the target column
transposed_merged_6106 <- transposed_merged_6106 %>%
  rename_with(~"Label", matches("Illness"))

write.csv(transposed_merged_6106, "transposed_merged_6106.csv", row.names = FALSE)

# Simplify the label values and filter out the sepsis patients and healthy control one, discards the diabetes patients
transposed_merged_6106 <- transposed_merged_6106 %>%
  mutate(Label = case_when(
    grepl("Sepsis", Label) ~ "Sepsis",
    grepl("Control/healthy", Label) ~ "healthy",
    TRUE ~ Label
  )) 

sepsis_limma_6106 <- transposed_merged_6106 %>%
  filter(Label %in% c("Sepsis", "healthy"))

write.csv(sepsis_limma_6106, "sepsis_limma_6106.csv", row.names = FALSE)


######################## GPL6947


# Retrieve all data related to first GPL 6106
expression_6947 <-exprs(geo_data[[2]])
feature_6947 <- fData(geo_data[[2]])
pheno_6947 <- pData(geo_data[[2]])


# expression data is not well normalized, we want to make sure
library(limma)
normalized_6947 <- normalizeBetweenArrays(expression_6947, method = "quantile")
boxplot(normalized_6947, main = "Boxplot of Quantile Normalized Data")

# Merge the expression data with fetaure data to map the genes symbols
sub_feature_6947 <- feature_6947[, c(1, 14)]
merged_6947 <- merge(normalized_6947, sub_feature_6947, by.x = "row.names", by.y = "ID", all.x = TRUE)
write.csv(merged_6947, "merged_6947.csv", row.names = FALSE)

# Search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELA2", "IL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

present_genes <- genes_of_interest[genes_of_interest %in% merged_6947$Symbol] #  ELANE as "ELA2", CXCL8 as IL8
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_6947$Symbol]

# Save the present and missing genes
write.csv(present_genes, "present_genesGSE13015-6947.csv", row.names = TRUE)
write.csv(missing_genes, "missing_genesGSE13015-6947.csv", row.names = TRUE)

# Filter the genes of interest
filtered_gene_6947 <- merged_6947 %>%
  filter(Symbol %in% genes_of_interest)

# Remove the column probe id and make unique the symbol column to transpose
filtered_gene_6947 <- filtered_gene_6947[,-1]
filtered_gene_6947$Symbol <- make.unique(as.character(filtered_gene_6947$Symbol))

# Transpose the data and make genes symbol as colnames
transposed_6947 <- filtered_gene_6947 %>%
  pivot_longer(
    cols = -Symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Symbol,
    values_from = Value
  )

# Merge the transposed data with pheno data to identift sepsis and healthy ones
sub_pheno_6947 <- pheno_6947[, c(2, 40)]
transposed_merged_6947 <- merge(transposed_6947, sub_pheno_6947, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Define the target column
transposed_merged_6947 <- transposed_merged_6947 %>%
  rename_with(~"Label", matches("Illness"))

write.csv(transposed_merged_6106, "transposed_merged_6106.csv", row.names = FALSE)

# Simplify the label values and filter out the sepsis patients and healthy control one, discards the diabetes patients
transposed_merged_6947 <- transposed_merged_6947 %>%
  mutate(Label = case_when(
    grepl("Sepsis", Label) ~ "Sepsis",
    grepl("Control/healthy", Label) ~ "healthy",
    TRUE ~ Label
  )) 

sepsis_limma_6947 <- transposed_merged_6947 %>%
  filter(Label %in% c("Sepsis", "healthy"))

write.csv(sepsis_limma_6947, "sepsis_limma_6947.csv", row.names = FALSE)








