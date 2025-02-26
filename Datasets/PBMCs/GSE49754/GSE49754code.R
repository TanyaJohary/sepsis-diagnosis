
# GSE49754

# Retrive geo datas as list
geo_data <- getGEO("GSE49754", GSEMatrix = TRUE)

# Retrieve the expression data
expression_data <- exprs(geo_data[[1]]) 

# Retrive phenotypr data of patients
pheno_data <- pData(geo_data[[1]])
colnames(pheno_data)
sub_pheno <- pheno_data[, c("geo_accession", "activators:ch1")]
write.csv(sub_pheno, "sub_phenoGSE49754.csv", row.names = TRUE)

# Retrieve the feature data for gene symbols
feature_data <- fData(geo_data[[1]]) 
colnames(feature_data)
sub_feature <- feature_data[, c("ID", "Symbol")]
write.csv(sub_feature, "sub_featureGSE49754.csv", row.names = TRUE)


# Merge expression data with gene symbols and probe_id
expression_with_symbols <- merge(expression_data, sub_feature, by.x = "row.names", by.y = "ID", all.x = TRUE)
write.csv(expression_with_symbols, "expression_with_symbolsGSE49754.csv", row.names = TRUE)


# search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "IL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)


present_genes <- genes_of_interest[genes_of_interest %in% expression_with_symbols$Symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% expression_with_symbols$Symbol] # "CRP" "LBP" "IFNA1" "IFNA2" "IFNB1" "CD177"

# save the present and missing genes for one sample
write.csv(present_genes, "present_genesGSE49754.csv", row.names = TRUE, quote = FALSE)
write.csv(missing_genes, "missing_genesGSE49754.csv", row.names = TRUE, quote = FALSE)

# search for aliases genes
alias_genes <- c(
  "IFNA1", "IFNA2", "IFNB1", "IFNA13", "IFN", "IFL", "LeIF A", "IFB", "IFNB", "IFG", "IFI",
  "PTX1", "C-Reactive Protein", "Pentraxin 1",
  "BPIFD2", "Lipopolysaccharide Binding Protein",
  "CD177", "PRV1", "NB1", "HNA2A", "HNA-2a", "NB1 GP", "PRV-1", "NB1 Glycoprotein"
)

# reload the expression with symbol and search for aliases
present_alias_genes <- alias_genes[alias_genes %in% c(expression_with_symbols$Symbol, expression_with_symbols$Synonyms)]


# Filter our genes of interest from dataset
filtered_genes <- expression_with_symbols %>%
  filter(Symbol %in% genes_of_interest)

# remove the probe_id column
filtered_genes <- filtered_genes[,-1]

# make unique value for each probe
filtered_genes$Symbol <- make.unique(as.character(filtered_genes$Symbol))

# Transpose the data
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

# Merge transposed data to pheno data 
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# rename the target column
merged_transposed <- merged_transposed %>%
  rename_with(~"Label", matches("activators"))

# Filter the sepsis and uninfected samples
sepsis_data <- merged_transposed %>%
  filter(Label %in% c("Septic plasma", "Uninfected plasma"))

# Rename the levels of label column to simplify
sepsis_data <- sepsis_data %>%
  mutate(Label = case_when(
    Label == "Septic plasma" ~ "Sepsis",
    Label == "Uninfected plasma" ~ "Control",
    TRUE ~ Label
  ))
write.csv(sepsis_data, "sepsis_dataGSE49754.csv", row.names = FALSE)


#################################################################################
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

present_genes <- genes_of_interest[genes_of_interest %in% expression_with_symbols$Symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% expression_with_symbols$Symbol] # "CRP"   "LBP"   "CXCL8" "IFNA1" "IFNA2" "IFNB1" "CD177" "TRIL" 

# Filter our genes of interest from dataset
filtered_genes <- expression_with_symbols %>%
  filter(Symbol %in% genes_of_interest)

# remove the probe_id column
filtered_genes <- filtered_genes[,-1]

# make unique value for each probe
filtered_genes$Symbol <- make.unique(as.character(filtered_genes$Symbol))

# Transpose the data
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

# Merge transposed data to pheno data 
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# rename the target column
merged_transposed <- merged_transposed %>%
  rename_with(~"Label", matches("activators"))

# Filter the sepsis and uninfected samples
sepsis_data <- merged_transposed %>%
  filter(Label %in% c("Septic plasma", "Uninfected plasma"))

# Rename the levels of label column to simplify
sepsis_data <- sepsis_data %>%
  mutate(Label = case_when(
    Label == "Septic plasma" ~ "Sepsis",
    Label == "Uninfected plasma" ~ "Control",
    TRUE ~ Label
  ))
write.csv(sepsis_data, "sepsis_60_dataGSE49754.csv", row.names = FALSE)


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

present_genes <- genes_of_interest[genes_of_interest %in% expression_with_symbols$Symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% expression_with_symbols$Symbol] # "LBP"   "CXCL8" "IFNA1" "IFNB1" "CD177"

# Filter our genes of interest from dataset
filtered_genes <- expression_with_symbols %>%
  filter(Symbol %in% genes_of_interest)

# remove the probe_id column
filtered_genes <- filtered_genes[,-1]

# make unique value for each probe
filtered_genes$Symbol <- make.unique(as.character(filtered_genes$Symbol))

# Transpose the data
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

# Merge transposed data to pheno data 
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# rename the target column
merged_transposed <- merged_transposed %>%
  rename_with(~"Label", matches("activators"))

# Filter the sepsis and uninfected samples
sepsis_data <- merged_transposed %>%
  filter(Label %in% c("Septic plasma", "Uninfected plasma"))

# Rename the levels of label column to simplify
sepsis_data <- sepsis_data %>%
  mutate(Label = case_when(
    Label == "Septic plasma" ~ "Sepsis",
    Label == "Uninfected plasma" ~ "Control",
    TRUE ~ Label
  ))
write.csv(sepsis_data, "sepsis_selected_dataGSE49754.csv", row.names = FALSE)
