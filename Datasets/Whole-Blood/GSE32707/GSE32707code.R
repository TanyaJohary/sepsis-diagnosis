
data <- read.delim("GSE32707_non_normalized.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Retrieve the dataset from GEO
geo_data <- getGEO("GSE32707", GSEMatrix = TRUE)

# Access Expression Data
expression_data <- exprs(geo_data[[1]])
meta_data <- Meta(geo_data[[1]]) # not provided
geo_data[[1]]@annotation

# Preview the metadata
pheno_data <- pData(geo_data[[1]])
head(metadata)
write.csv(pheno_data, "pheno_dataGSE32707.csv", row.names = TRUE)

# Preview the feature data
feature_data <- fData(geo_data[[1]])
head(feature_data)  

# Extract relevant columns from feature_data
gene_mapping <- feature_data[, c("ID", "Symbol")]

# Rename for clarity
colnames(gene_mapping) <- c("Probe_ID", "Gene_Symbol") 


# Merge Gene Mapping with Expression Data
merged_data <- merge(expression_data, gene_mapping, by.x = "row.names", by.y = "Probe_ID", all.x = TRUE)
 

#save the expression data with symbol
write.csv(merged_data, "Merged_Expression_Data32707.csv", row.names = FALSE)

# search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "IL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

present_genes <- genes_of_interest[genes_of_interest %in% merged_data$Gene_Symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_data$Gene_Symbol]

# search for aliases genes
alias_genes <- c(
  "IFNA1", "IFNA2", "IFNB1",
  "ELANE", "HLE", "PMN-E", "ELA2", "HNE", "NE",
  "MDNCF", "NAP-1", "GCP-1", "SCYB8", "LYNAP", "IL-8", "IL8"
)

# save the present and missing genes
write.csv(present_genes, "present_genesGSE32707.csv", row.names = FALSE, quote = FALSE)
write.csv(missing_genes, "missing_genesGSE32707.csv", row.names = FALSE, quote = FALSE)

# filter out the genes of interest
filtered_genes_data <- merged_data %>%
  filter(Gene_Symbol %in% genes_of_interest)

write.csv(filtered_genes_data, "filtered_genes_data_bysymbol.csv", row.names = FALSE)

# Transpose data
transposed_data <- t(filtered_genes_data)
transposed_data <- transposed_data[-1,]

# Make gene symbol as colnames
colnames(transposed_data) <- transposed_data[nrow(transposed_data), ]

# Get geo accession number and patients status from pheno data
sub_pheno <- pheno_data[, c(2, 8)]

# Merge transposed data with sub pheno
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "row.names", by.y = "geo_accession", all.x = TRUE)
merged_transposed <- merged_transposed[-1,]

write.csv(merged_transposed, "merged_transposedGSE32707.csv", row.names = FALSE)

# Make unique the colnames to prevent duplication
colnames(merged_transposed) <- make.unique(colnames(merged_transposed))

# Define the target column as label
merged_transposed <- merged_transposed %>%
  rename_with(~"Label", matches("source"))

#################### D0-sepsis
# filter out the sepsis patients for Day 0 and control ones
filtered_sepsis_D0 <- merged_transposed %>%
  filter(Label %in% c("Sepsis Day 0", "untreated"))

# Change the value of target column for D0
filtered_sepsis_day0 <- filtered_sepsis_D0 %>%
  mutate(Label = case_when(
    Label == "Sepsis Day 0" ~ "Sepsis",
    Label == "untreated" ~ "Control",
    TRUE ~ Label
  )) %>%
  rename(Sample = Row.names)

write.csv(filtered_sepsis_day0, "filtered_sepsD0_GSE32707.csv", row.names = FALSE)  


###################### D7-sepsis
# filter out the sepsis patients for Day 7 and control ones
filtered_sepsis_D7 <- merged_transposed %>%
  filter(Label %in% c("Sepsis Day 7", "untreated"))

# # Change the value of target column for D7
filtered_sepsis_day7 <- filtered_sepsis_D7 %>%
  mutate(Label = case_when(
    Label == "Sepsis Day 7" ~ "Sepsis",
    Label == "untreated" ~ "Contrl",
    TRUE ~ Label
  )) %>%
  rename(Sample = Row.names)

write.csv(filtered_sepsis_day7, "filtered_sepsD7_GSE32707.csv", row.names = FALSE)


####################### D0 and D7 sepsis
# filter out the sepsis patients for Day 0, 7 and control ones
filtered_sepsis_all <- merged_transposed %>%
  filter(Label %in% c("Sepsis Day 0", "Sepsis Day 7", "untreated"))

# Change the value of target column for D7
filtered_sepsis_D07 <- filtered_sepsis_all %>%
  mutate(Label = case_when(
    Label == "Sepsis Day 0" ~ "Sepsis",
    Label == "Sepsis Day 7" ~ "Sepsis",
    Label == "untreated" ~ "Control",
    TRUE ~ Label
  )) %>%
  rename(Sample = Row.names)
  
write.csv(filtered_sepsis_D07, "filtered_sepsis_all_GSE32707.csv", row.names = FALSE)  




























