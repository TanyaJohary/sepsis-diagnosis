

################################################################################
############### Deseq2 normalized dataset for random forest ####################
################################################################################
library(DESeq2)
geo_data <- getGEO("GSE185263", GSEMatrix = TRUE)

expression_data <- exprs(geo_data[[1]]) # empty

pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c(2, 49)]

feature_data <- fData(geo_data[[1]]) # empty

raw_count <- read.delim("GSE185263_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
row.names(raw_count) <- raw_count$GeneID
raw_count <- raw_count[, -1]

coldata <- sub_pheno %>%
  rename_with(~"Label", matches("disease"))

# Ensure sample names match exactly
all(row.names(coldata) %in% colnames(raw_count))
all(colnames(raw_count) %in% rownames(coldata)) 

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
   countData = raw_count,
   colData = coldata,
   design = ~ Label
 )

# Filter lowly expressed genes
dds <- dds[rowSums(counts(dds)) > 10, ] 

# Apply Variance Stabilizing Transformation
vsd <- vst(dds, blind = TRUE)

# Extract normalized data
normalized_counts <- assay(vsd)

# Save filtered & normalized data
write.csv(normalized_counts, "normalized_counts_VST_GSE185263.csv", row.names = TRUE)

# Read the annotataion data
annotation_data <- read.delim("Human.GRCh38.p13.annot.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sub_annot <- annotation_data[, c(1, 2)] 

# Merge annotation data to normalized data to map the genes symbols
merged_normalized <- merge(normalized_counts, sub_annot, by.x = "row.names", by.y = "GeneID", all.x = TRUE)
write.csv(merged_normalized, "normalized_with_symbolsGSE185263.csv", row.names = FALSE)

# Search for genes of interest by gene symbols
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

present_genes <- genes_of_interest[genes_of_interest %in% merged_normalized$Symbol] 
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_normalized$Symbol] # IFNA1, IFNA2 missing

write.csv(present_genes, "present_genesGSE185263.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE185263.csv", row.names = FALSE)

# Filter our genes of interest from normalized data
filtered_genes <- merged_normalized %>%
  filter(Symbol %in% present_genes)

# Remove the probe id column 
filtered_genes <- filtered_genes[, -1]

# Transpose the data to make the gene symbols as header
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

# Merge to pheno data to identify the sepsis patients
sepsis_data <- merge(transposed_data, coldata, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)
write.csv(sepsis_data, "sepsis_dataGSE185263.csv", row.names = FALSE)

################################################################################
################################ 60 Genes ######################################
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

present_genes <- genes_of_interest[genes_of_interest %in% merged_normalized$Symbol] 
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_normalized$Symbol]

# Filter our genes of interest from normalized data
filtered_genes <- merged_normalized %>%
  filter(Symbol %in% present_genes)

# Remove the probe id column 
filtered_genes <- filtered_genes[, -1]

# Transpose the data to make the gene symbols as header
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

# Merge to pheno data to identify the sepsis patients
sepsis_data <- merge(transposed_data, coldata, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)
write.csv(sepsis_data, "sepsis_60_dataGSE185263.csv", row.names = FALSE)

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

present_genes <- genes_of_interest[genes_of_interest %in% merged_normalized$Symbol] 
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_normalized$Symbol]

# Filter our genes of interest from normalized data
filtered_genes <- merged_normalized %>%
  filter(Symbol %in% present_genes)

# Remove the probe id column 
filtered_genes <- filtered_genes[, -1]

# Transpose the data to make the gene symbols as header
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

# Merge to pheno data to identify the sepsis patients
sepsis_data <- merge(transposed_data, coldata, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)
write.csv(sepsis_data, "sepsis_selected_dataGSE185263.csv", row.names = FALSE)
