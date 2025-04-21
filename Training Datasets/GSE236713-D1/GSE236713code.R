


# Retrieve the geo data as list
geo_data <- getGEO("GSE236713", GSEMatrix = TRUE)

# Retrieve the expression data
expression_data <- exprs(geo_data[[1]])

# Retrieve the phenotype data
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c(2, 8, 42)]

# retrieve the feature data to map genes symbols
feature_data <- fData(geo_data[[1]]) 
sub_feature <- feature_data[, c("ID", "GENE_SYMBOL")]

# Merge expression data to feature data to mapp the genes symbols
merged_data <- merge(expression_data, sub_feature, by.x = "row.names", by.y = "ID", all.x = TRUE)


# genes of interest names
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "IL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

present_genes <- genes_of_interest[genes_of_interest %in% merged_data$GENE_SYMBOL]
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_data$GENE_SYMBOL]

write.csv(present_genes, "present_genesGSE236713.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE236713.csv", row.names = FALSE)

# Filter out our genes of intereset
filtered_genes <- merged_data %>%
  filter(GENE_SYMBOL %in% genes_of_interest)

# Make unique values for each probe
filtered_genes$GENE_SYMBOL <- make.unique(as.character(filtered_genes$GENE_SYMBOL))
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

# Merge to sub pheno data
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# filter sepsis patients
sepsis_filtered <- merged_transposed %>%
  filter(`disease:ch1` %in% c("Sepsis", "Control"))
write.csv(sepsis_filtered, "sepsis_filtered.csv")


# Filter dataset for Control samples and Sepsis samples with Day 1 label
sepsis_day1 <- sepsis_filtered %>%
  filter(`disease:ch1` == "Control" | (`disease:ch1` == "Sepsis" & `source_name_ch1` == "Peripheral Blood, Day: 1"))

sepsis_D1 <- sepsis_day1 %>%
  rename_with(~"Label", matches("disease:ch1"))
sepsis_D1 <- sepsis_D1[, -78]
write.csv(sepsis_D1, "sepsis_D1_dataGSE236713.csv", row.names = FALSE)

################################# Day2
# Filter dataset for Control samples and Sepsis samples with Day 2 label
sepsis_day2 <- sepsis_filtered %>%
  filter(`disease:ch1` == "Control" | (`disease:ch1` == "Sepsis" & `source_name_ch1` == "Peripheral Blood, Day: 2"))

sepsis_D2 <- sepsis_day2 %>%
  rename_with(~"Label", matches("disease:ch1"))
sepsis_D2 <- sepsis_D2[, -78]  

write.csv(sepsis_D2, "sepsis_D2_dataGSE236713.csv", row.names = FALSE)

################################ Day5
# Filter dataset for Control samples and Sepsis samples with Day 5 label
sepsis_day5 <- sepsis_filtered %>%
  filter(`disease:ch1` == "Control" | (`disease:ch1` == "Sepsis" & `source_name_ch1` == "Peripheral Blood, Day: 5"))

sepsis_D5 <- sepsis_day5 %>%
  rename_with(~"Label", matches("disease:ch1"))
sepsis_D5 <- sepsis_D5[, -78]  

write.csv(sepsis_D5, "sepsis_D5_dataGSE236713.csv", row.names = FALSE)

