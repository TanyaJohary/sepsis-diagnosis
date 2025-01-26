

################################################################################
############################# Prepare data for RF###############################
################################################################################

# Retrieve the geo data list
geo_data <- getGEO("GSE100159", GSEMatrix = TRUE)

# Retrieve the expression dataset
expression_data <- exprs(geo_data[[1]])

# Retrieve the phenodata 
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c(2,32)]

# Retrieve the feature data 
feature_data <- fData(geo_data[[1]])
sub_feature <- feature_data[, c(1,13)]

# Mereg the expression data and sub feature to map genes symbols
merged_data <- merge(expression_data, sub_feature, by.x = "row.names", by.y = "ID", all.x = TRUE)

# Search for our genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELA2", "IL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

present_genes <- genes_of_interest[genes_of_interest %in% merged_data$Symbol] # ELANE as "ELA2", CXCL8 as "IL8"
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_data$Symbol]

# Save the present and missing genes
write.csv(present_genes, "present_genesGSE100159.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE100159.csv", row.names = FALSE)

# Filter our data for genes of interest
filtered_genes <- merged_data %>%
  filter(Symbol %in% genes_of_interest)

# make unique value for each probe of each gene
filtered_genes$Symbol <- make.unique(as.character(filtered_genes$Symbol))

# Remove the probe id column
filtered_genes <- filtered_genes[,-1]

# Transpose the data to have the genes symbols as header in table
transposed_data <- filtered_genes %>%
  pivot_longer(
    cols = -Symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Symbol,
    values_from = Value
  )

# Merge transposed data to pheno data to identift sepsis patients
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Filter rows where all contains sepsis and control without recovery patients
sepsis_data <- merged_transposed %>%
  filter(`sample class:ch1` != "Recovery")

# simplify the values of target column
sepsis_data <- sepsis_data %>%
  rename_with(~"Label", matches("sample class"))

sepsis_data <- sepsis_data %>%
  mutate(Label = case_when(
    grepl("melioidosis", Label) ~ "Sepsis",
    grepl("other infection", Label) ~ "Sepsis",
    TRUE ~ Label
  ))

table(sepsis_data$Label) # 12=Control, 33=Sepsis
# Save the prepared data for downstream analysis
write.csv(sepsis_data, "sepsis_dataGSE100159.csv", row.names = FALSE)


