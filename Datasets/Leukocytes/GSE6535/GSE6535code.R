

# GSE6535

# Retrieve the geo data
geo_data <- getGEO("GSE6535", GSEMatrix = TRUE)

# Retrieve the expression data
expression_data <- exprs(geo_data[[1]])

# Retrieve the phenotype data
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c("geo_accession", "characteristics_ch1")]
write.csv(sub_pheno, "sub_phenoGSE6535.csv", row.names = FALSE)

# Retrieve feature data of gene annotation
feature_data <- fData(geo_data[[1]])
colnames(feature_data) # does not hhave gene symbols

# Search for gene symbols based on GenBank IDs
library(org.Hs.eg.db)
library(AnnotationDbi)

# Map GenBank IDs to gene symbols
gene_mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = feature_data$GB_ACC,
  columns = c("SYMBOL"),
  keytype = "ACCNUM"
)

# Merge expression data and sub feature data
expression_with_symbols <- merge(expression_data, gene_mapping, by.x = "row.names", by.y = "ACCNUM", all.x = TRUE)
write.csv(expression_with_symbols, "expression_with_symbolsGSE6535.csv", row.names = TRUE)

# Search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

present_genes <- genes_of_interest[genes_of_interest %in% expression_with_symbols$SYMBOL]
missing_genes <- genes_of_interest[!genes_of_interest %in% expression_with_symbols$SYMBOL] # "IFNA1" "HLA.DRA" "CCR2"  

# Save the present and missing genes
write.csv(present_genes, "present_genesGSE6535.csv", row.names = TRUE)
write.csv(missing_genes, "missing_genesGSE6535.csv", row.names = TRUE)


# search for aliases genes
alias_genes <- c(
  "IFNA1", "IFNA2", "IFNB1", "IFNA13", "IFN", "IFL", "LeIF A", "IFB", "IFNB", "IFG", "IFI", "IFNA@", "IFN-AlphaD",
  "HLA-DRA", "HLA-DRA1", "HLA_DR", "HLADR", "HLADRA",
  "CCR2", "CC-CKR-2", "MCP-1-R", "CMKBR2", "CD192", "CKR2", "FLJ78302", "CCR-2", "MCP-1 Receptor", "CCR2A",
  "CCR2B", "CKR2A", "CKR2B", "PCLUD"
)

# reload the expression with symbol and search for aliases
present_alias_genes <- alias_genes[alias_genes %in% expression_with_symbols$SYMBOL]


# Filter the expression data for genes of interest
filtered_genes <- expression_with_symbols %>%
  filter(SYMBOL %in% genes_of_interest) 

# Remove the accession number
filtered_genes <- filtered_genes[, -1]

# Make unique value for each peobe
filtered_genes$SYMBOL <- make.unique(as.character(filtered_genes$SYMBOL))

# Transpose the data
transposed_data <- filtered_genes %>%
  pivot_longer(
    cols = -SYMBOL,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = SYMBOL,
    values_from = Value
  )

# Merge to phenotype data
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Rename the label
sepsis_data <- merged_transposed %>%
  rename_with(~"Label", matches("characteristics"))

# simplify the target columns value
sepsis_data <- sepsis_data %>%
  mutate(Label = case_when(
    Label == "Control Patient" ~ "Control",
    Label == "Sepsis Patient" ~ "Sepsis",
    TRUE ~ Label
  ))

# calculate the missing values
sum(is.na(sepsis_data))
rowSums(is.na(sepsis_data))
mean(is.na(sepsis_data)) * 100

# Remove the Na rows
sepsis_data_clean <- sepsis_data[rowSums(is.na(sepsis_data)) == 0, ]

write.csv(sepsis_data_clean, "sepsis_dataGSE6535.csv", row.names = FALSE)


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


present_genes <- genes_of_interest[genes_of_interest %in% expression_with_symbols$SYMBOL]
missing_genes <- genes_of_interest[!genes_of_interest %in% expression_with_symbols$SYMBOL]


# filter the genes of interest expression data

# Load dplyr package
library(dplyr)

# Filter the expression data for genes of interest
filtered_genes <- expression_with_symbols %>%
  filter(SYMBOL %in% genes_of_interest) 

# Remove the accession number
filtered_genes <- filtered_genes[, -1]

# Make unique value for each peobe
filtered_genes$SYMBOL <- make.unique(as.character(filtered_genes$SYMBOL))

# Transpose the data
transposed_data <- filtered_genes %>%
  pivot_longer(
    cols = -SYMBOL,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = SYMBOL,
    values_from = Value
  )

# Merge to phenotype data
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Rename the label
sepsis_data <- merged_transposed %>%
  rename_with(~"Label", matches("characteristics"))

# simplify the target columns value
sepsis_data <- sepsis_data %>%
  mutate(Label = case_when(
    Label == "Control Patient" ~ "Control",
    Label == "Sepsis Patient" ~ "Sepsis",
    TRUE ~ Label
  ))

# calculate the missing values
sum(is.na(sepsis_data))
rowSums(is.na(sepsis_data))
mean(is.na(sepsis_data)) * 100

# Remove the Na rows
sepsis_data_clean <- sepsis_data[rowSums(is.na(sepsis_data)) == 0, ]

write.csv(sepsis_data_clean, "sepsis_60_dataGSE6535.csv", row.names = FALSE)

################################################################################
################################# Selected genes ###############################
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


present_genes <- genes_of_interest[genes_of_interest %in% expression_with_symbols$SYMBOL]
missing_genes <- genes_of_interest[!genes_of_interest %in% expression_with_symbols$SYMBOL]


# filter the genes of interest expression data

# Load dplyr package
library(dplyr)

# Filter the expression data for genes of interest
filtered_genes <- expression_with_symbols %>%
  filter(SYMBOL %in% genes_of_interest) 

# Remove the accession number
filtered_genes <- filtered_genes[, -1]

# Make unique value for each peobe
filtered_genes$SYMBOL <- make.unique(as.character(filtered_genes$SYMBOL))

# Transpose the data
transposed_data <- filtered_genes %>%
  pivot_longer(
    cols = -SYMBOL,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = SYMBOL,
    values_from = Value
  )

# Merge to phenotype data
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Rename the label
sepsis_data <- merged_transposed %>%
  rename_with(~"Label", matches("characteristics"))

# simplify the target columns value
sepsis_data <- sepsis_data %>%
  mutate(Label = case_when(
    Label == "Control Patient" ~ "Control",
    Label == "Sepsis Patient" ~ "Sepsis",
    TRUE ~ Label
  ))

# calculate the missing values
sum(is.na(sepsis_data))
rowSums(is.na(sepsis_data))
mean(is.na(sepsis_data)) * 100

# Remove the Na rows
sepsis_data_clean <- sepsis_data[rowSums(is.na(sepsis_data)) == 0, ]

write.csv(sepsis_data_clean, "sepsis_selected_dataGSE6535.csv", row.names = FALSE)

