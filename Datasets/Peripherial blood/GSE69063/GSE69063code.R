
geo_data <- getGEO("GSE69063", GSEMatrix = TRUE)

# get expression data
expression_data <- exprs(geo_data[[1]])

# get the phenotype data of patients
pheno_data <- pData(geo_data[[1]])
colnames(pheno_data)
sub_pheno <- pheno_data[, c("geo_accession", "title", "characteristics_ch1.2", "characteristics_ch1.3")]
write.csv(sub_pheno, "sub_phenoGSE69063.csv", row.names = TRUE)


# get annotation data of genes
feature_data <- fData(geo_data[[1]])
colnames(feature_data)

# Connect to Ensembl biomart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Map Entrez IDs to Gene Symbols
mapping <- getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
                 filters = "entrezgene_id",
                 values = feature_data$ENTREZ_GENE_ID,
                 mart = ensembl)

# rename the coulmn in mapping data to merge easier
colnames(mapping)[colnames(mapping) == "entrezgene_id"] <- "ENTREZ_GENE_ID"

# Merge the mappging file to feature data
feature_data_merged <- merge(feature_data, mapping, by = "ENTREZ_GENE_ID", all.x = TRUE)
write.csv(feature_data_merged, "feature_data_with_gene_symbols.csv", row.names = FALSE)

# Merge the expression data with gene symbols
expression_with_symbols <- merge(expression_data, feature_data_merged, by.x = "row.names", by.y = "ID", all.x = TRUE)
write.csv(expression_with_symbols, "expression_with_symbolsGSE69063.csv", row.names = TRUE)


# search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

present_genes <- genes_of_interest[genes_of_interest %in% expression_with_symbols$hgnc_symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% expression_with_symbols$hgnc_symbol]

write.csv(present_genes, "present_genesGSE69063.csv", row.names = TRUE, quote = FALSE)
write.csv(missing_genes, "missing_genesGSE69063.csv", row.names = TRUE, quote = FALSE)


# search for aliases genes
alias_genes <- c(
  "ELANE", "HLE", "PMN-E", "ELA2", "HNE", "NE",
  "MDNCF", "NAP-1", "GCP-1", "SCYB8", "LYNAP", "IL-8", "IL8",
  "FCGR1", "FCG1", "CD64", "CD64A", "IGFR1", "FCRI", "FcRI"
)

# reload the expression with symbol and search for aliases
present_alias_genes <- alias_genes[alias_genes %in% expression_with_symbols$hgnc_symbol]

# Read the expression with symbol data
expression <- read.csv("expression_with_symbolsGSE69063.csv")

present_genes <- genes_of_interest[genes_of_interest %in% expression$hgnc_symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% expression$hgnc_symbol]

# filter out the genes of interest from expression data
filtered_genes <- expression %>%
  filter(hgnc_symbol %in% genes_of_interest)

# delete the extra columns X and probe id and entrez gene id
library(dplyr)

filtered_genes <- dplyr::select(filtered_genes, -c(X, Row.names, ENTREZ_GENE_ID))

# transpose the data to have genes symbol as column names
transposed_data <- filtered_genes %>%
  pivot_longer(
    cols = -hgnc_symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = hgnc_symbol,
    values_from = Value
  )

# Read the sub_phenotype data to merge with expression data and identify the sepsis patients
sub_pheno <- read.csv("sub_phenoGSE69063.csv")
sub_pheno <- dplyr::select(sub_pheno, -c(X, characteristics_ch1.2, characteristics_ch1.3))

merged_trasposed_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Filter rows where title contains "sepsis patients" or "healthy control"
filtered_sepsis_data <- merged_trasposed_data %>%
  filter(grepl("sepsis patient", title, ignore.case = TRUE) | 
           grepl("healthy control", title, ignore.case = TRUE))

# We have 3 timeline of sampling from each sepsis patient
######### separating each timeline T0 
T0_sepsis_data <- filtered_sepsis_data %>%
  filter(grepl("T0", title, ignore.case = TRUE) |
           grepl("healthy control", title, ignore.case = TRUE))

# Simplify the 'title' column values
T0_sepsis_data <- T0_sepsis_data %>%
  mutate(title = case_when(
    grepl("healthy control", title, ignore.case = TRUE) ~ "healthy control",
    grepl("sepsis patient", title, ignore.case = TRUE) ~ "sepsis",
    TRUE ~ title  # Retain the original value if no match is found
  ))
 
# Change the column title to label
T0_sepsis_labeled <- T0_sepsis_data %>%
  rename_with(~ "Label", matches("title"))

# Save the dataset for T0
write.csv(T0_sepsis_labeled, "T0_sepsis_labeledGSE69063.csv", row.names = FALSE)

########## separate the T1 dataset
T1_filtered_data <- filtered_sepsis_data %>%
  filter(grepl("T1", title, ignore.case = TRUE) |
           grepl("healthy control", title, ignore.case = TRUE))

# Simplify the 'title' column values
T1_sepsis_data <- T1_filtered_data %>%
  mutate(title = case_when(
    grepl("healthy control", title, ignore.case = TRUE) ~ "healthy control",
    grepl("sepsis patient", title, ignore.case = TRUE) ~ "sepsis",
    TRUE ~ title
  ))

# Change the column title to label
T1_sepsis_labeled <- T1_sepsis_data %>%
  rename_with(~ "Label", matches("title"))

# Save the T1 dataset
write.csv(T1_sepsis_labeled, "T1_sepsis_labeledGSE69063.csv", row.names = FALSE)


########### separate T2 dataset
T2_filtered_data <- filtered_sepsis_data %>%
  filter(grepl("T2", title, ignore.case = TRUE) |
           grepl("healthy control", title, ignore.case = TRUE))

# Simplify the 'title' column values
T2_sepsis_data <- T2_filtered_data %>%
  mutate(title = case_when(
    grepl("healthy control", title, ignore.case = TRUE) ~ "healthy control",
    grepl("sepsis patient", title, ignore.case = TRUE) ~ "sepsis",
    TRUE ~ title
  ))

# Change the column title to label in T2 dataset
T2_sepsis_labeled <- T2_sepsis_data %>%
  rename_with(~ "Label", matches("title"))

# Save the T2 dataset
write.csv(T2_sepsis_labeled, "T2_sepsis_labeledGSE69063.csv", row.names = FALSE)


x <- read.csv("T1_sepsis_labeledGSE69063.csv")
