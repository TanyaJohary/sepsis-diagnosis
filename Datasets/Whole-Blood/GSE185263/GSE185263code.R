

##################################################################
############prepare dataset for random forest#####################
##################################################################

expression_data <- read.csv("GSE185263_raw_counts.csv.gz")

install.packages("BiocManager")
BiocManager::install("vsn")
library(vsn)

# Extract Numeric Data for Normalization
count_matrix <- as.matrix(expression_data[, -1])  # Extract all columns except the first
rownames(count_matrix) <- expression_data$ensembl_gene_id  # Set gene IDs as rownames

# Apply VSN to Normalize the Data
vsn_normalized <- vsn::justvsn(count_matrix)

# Convert Normalized Data Back to a Data Frame
vsn_data <- as.data.frame(vsn_normalized)
vsn_data$ensembl_gene_id <- rownames(vsn_normalized)

# Save or Inspect the Normalized Data
head(vsn_data)  # Inspect the first few rows
write.csv(vsn_data, "vsn_normalized_data.csv", row.names = FALSE)

# Verify that the row names of normalized_data match the original ensembl_gene_id
all(rownames(vsn_data) == expression_data$ensembl_gene_id)

# quality control of normalization
library(vsn)
meanSdPlot(vsn_normalized)

# Use Biomart to get gene symbols from gene ensembl id
library(biomaRt)
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

gene_ids <- rownames(vsn_data)

gene_symbols <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                      filters = "ensembl_gene_id", 
                      values = gene_ids, 
                      mart = ensembl)

# Merge gene symbols with your data
vsn_with_symbols <- merge(vsn_data, gene_symbols, 
                                      by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE)

write.csv(vsn_with_symbols, "vsn_with_symbols185263.csv", row.names = FALSE)

# Search for genes of interest by gene symbols
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

present_genes <- genes_of_interest[genes_of_interest %in% vsn_with_symbols$hgnc_symbol] # all the genes are present
missing_genes <- genes_of_interest[!genes_of_interest %in% vsn_with_symbols$hgnc_symbol]

# List of Ensembl IDs of genes of interest
genes_ensembl_id <- c(
  "ENSG00000110680", "ENSG00000132693", "ENSG00000124731", "ENSG00000011422", "ENSG00000129988",
  "ENSG00000136244", "ENSG00000189403", "ENSG00000197561", "ENSG00000169429", "ENSG00000136634",
  "ENSG00000125538", "ENSG00000143546", "ENSG00000163220", "ENSG00000163221", "ENSG00000170458",
  "ENSG00000150337", "ENSG00000169896", "ENSG00000171860", "ENSG00000197405", "ENSG00000148346",
  "ENSG00000169245", "ENSG00000111537", "ENSG00000172724", "ENSG00000131142", "ENSG00000168329",
  "ENSG00000089041", "ENSG00000163661", "ENSG00000121858", "ENSG00000118113", "ENSG00000100985",
  "ENSG00000204287", "ENSG00000232810", "ENSG00000172936", "ENSG00000162711", "ENSG00000137462",
  "ENSG00000136869", "ENSG00000148400", "ENSG00000171791", "ENSG00000188389", "ENSG00000108691",
  "ENSG00000118520", "ENSG00000115590", "ENSG00000204936", "ENSG00000102837", "ENSG00000112062",
  "ENSG00000162692", "ENSG00000090339", "ENSG00000184557", "ENSG00000107485", "ENSG00000126353",
  "ENSG00000121807", "ENSG00000100644", "ENSG00000188379", "ENSG00000197919", "ENSG00000171855"
)

present_id <- genes_ensembl_id[genes_ensembl_id %in% vsn_with_symbols$ensembl_gene_id]. # all the genes are present
missing_id <- genes_ensembl_id[!genes_ensembl_id %in% vsn_with_symbols$ensembl_gene_id]

# save the present and missing genes
write.csv(present_genes, "present_genesGSE185263.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE185263.csv", row.names = FALSE)

# filter out our genes of interest
filtered_vsn_data <- vsn_with_symbols %>%
  filter(ensembl_gene_id %in% genes_ensembl_id)

# Remove the ensembl id from data
filtered_vsn_data <- dplyr::select(filtered_vsn_data, -c(1, 394))

# Transpose the data to genes symbols become the column names
transposed_data <- filtered_vsn_data %>%
  pivot_longer(
    cols = -hgnc_symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = hgnc_symbol,
    values_from = Value
  )

# Retrieve the pheno data to merge with transposed data to identify sepsis patients
gse_data <- getGEO("GSE185263", GSEMatrix = TRUE)

pheno <- pData(gse_data[[1]])
sub_pheno <- pheno[, c(1, 49, 50)]

# Merge sub pheno data to transposed adat set based on sample title
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "title", all.x = TRUE)
write.csv(merged_transposed, "merged_transposedGSE185263.csv", row.names = FALSE)

# Remove the extra column for have just label as sepsis and healthy
sepsis_data <- dplyr::select(merged_transposed, -c(58))

# Define the label column
sepsis_labeled_data <- sepsis_data %>%
  rename_with(~"Label", matches("disease state"))
write.csv(sepsis_labeled_data, "sepsis_labeled_dataGSE185263.csv", row.names = FALSE)

# Define the column for outcome of disease (prognostics)
sepsis_prognostic <- dplyr::select(merged_transposed, -c(57))
sepsis_prognos_labeled <- sepsis_prognostic %>%
  rename_with(~"Label", matches("hospital"))
write.csv(sepsis_prognos_labeled, "sepsis_prognos_labeledGSE185263.csv", row.names = FALSE)
