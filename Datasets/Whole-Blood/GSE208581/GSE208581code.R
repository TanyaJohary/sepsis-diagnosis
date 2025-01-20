# Read the raw data because the mRNA count data has too many out ranged values that random forest does not work properly
raw_data <- read.delim("GSE208581_raw_counts.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

if (!requireNamespace("vsn", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("vsn")
}
library(vsn)

# Save Gene IDs
gene_ids <- raw_data$Ensembl_Gene_ID

# Remove Gene IDs from the dataset (keep only numeric data)
raw_data <- raw_data[, -1]

# Convert to matrix and normalize using VSN
raw_matrix <- as.matrix(raw_data)

# Apply variance-stabilizing normalization
vsn_model <- vsn2(raw_matrix)  # Fit the VSN model
vsn_normalized <- predict(vsn_model, raw_matrix)  # Apply the model to the data

# To confirm the normalization was successful
meanSdPlot(raw_matrix, main = "Before VSN Normalization")
meanSdPlot(vsn_normalized, main = "After VSN Normalization")

# Convert normalized matrix to a data frame
vsn_normalized_df <- as.data.frame(vsn_normalized)

# Add Gene IDs back to the normalized data
vsn_normalized_df$Ensembl_Gene_ID <- gene_ids

# Reorder columns to place Gene IDs first
vsn_normalized_df <- vsn_normalized_df[, c(ncol(vsn_normalized_df), 1:(ncol(vsn_normalized_df) - 1))]

# Save to a CSV file
write.csv(vsn_normalized_df, "GSE208581_vsn_normalized_with_ids.csv", row.names = FALSE)

all(rownames(vsn_normalized) == gene_ids)  # Should return TRUE



# Use the Ensembl database to retrieve gene information:
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

gene_id <- vsn_normalized_df$Ensembl_Gene_ID

# Retrieve Gene Symbols
gene_symbols <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = gene_id,
  mart = ensembl
)

# Merge the retrieved gene symbols with your original dataset:
vsn_norm_with_symbols <- merge(vsn_normalized_df, gene_symbols, 
                                      by.x = "Ensembl_Gene_ID", 
                                      by.y = "ensembl_gene_id", 
                                      all.x = TRUE)

# Save the updated dataset with gene symbols:
write.csv(vsn_norm_with_symbols, "vsn_norm_with_symbols.csv", row.names = FALSE)




# search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

present_genes <- genes_of_interest[genes_of_interest %in% vsn_norm_with_symbols$hgnc_symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% vsn_norm_with_symbols$hgnc_symbol]

# search for aliases genes
alias_genes <- c(
  "IFNA1", "IFNA2", "IFNB1",
  "ELANE", "HLE", "PMN-E", "ELA2", "HNE", "NE",
  "MDNCF", "NAP-1", "GCP-1", "SCYB8", "LYNAP", "IL-8", "IL8",
  "MRP8", "P8", "S100-A8", "MRP-8", "CGLA", "CP-10", "CAGA", "CFAG",
  "MRP-14", "MRP14", "P14", "S100-A9", "MAC387", "60B8AG", "LIAG", "CGLB", "CAGB", "CFAG", "MIF"
)

write.csv(present_genes, "present_genesGSE208581.csv", row.names = FALSE, quote = FALSE)
write.csv(missing_genes, "missing_genesGSE208581.csv", row.names = FALSE, quote = FALSE)


# search for gene by ensemble gene id
missing_gene_id <- c("ENSG00000143546", "ENSG00000163220") # genes "S100A8" "S100A9" are absent by symbols and ensembl ids
present_gene_id <- missing_gene_id[missing_gene_id %in% vsn_norm_with_symbols$Ensembl_Gene_ID]

# filter our genes of interest from dataset
filtered_vsn_genes <- vsn_norm_with_symbols %>%
  filter(hgnc_symbol %in% present_genes)

# Remove the ensembl gene id column
filtered_vsn_genes <- dplyr::select(filtered_vsn_genes, -Ensembl_Gene_ID)

# Transpose data to make genes symbol as colnames
transposed_filtered_vsn <- filtered_vsn_genes %>%
  pivot_longer(
    cols = -hgnc_symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = hgnc_symbol,
    values_from = Value
  )

# extract phenotype data
gse <- getGEO("GSE208581", GSEMatrix = TRUE)
pheno_data <- pData(gse[[1]])
sub_phenodata <- pheno_data[, c(18, 46)]
write.csv(sub_phenodata, "pheno_dataGSE280581.csv", row.names = TRUE)

# Merge the sub_pheno data to transposed data
merged_transposed_vsn <- merge(transposed_filtered_vsn, sub_phenodata, by.x = "Sample", by.y = "description", all.x = TRUE)

# Filter out just Sepsis+ patients and SIRS- as control group
sepsis_vsn_data <- merged_transposed_vsn %>%
  filter(`condition:ch1` %in% c("Sepsis+", "SIRS-"))

# define a target column as label
sepsis_vsn_data <- sepsis_vsn_data %>%
  rename_with(~"Label", matches("condition"))

# Save the ready sepsis data
write.csv(sepsis_vsn_data, "sepsis_vsn_dataGSE208581.csv", row.names = FALSE)
