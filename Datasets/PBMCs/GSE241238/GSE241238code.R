
################################################################################
#######################Normalize the raw data for RF############################
################################################################################

Raw_data <- read.delim("GSE241238_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

row.names(Raw_data) <- Raw_data$GeneID 
Raw_data <- Raw_data[, -1]

# Retrieve the geo data as list
geo_data <- getGEO("GSE241238", GSEMatrix = TRUE)

# Retrieve the phenotype data
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c("geo_accession", "condition:ch1")]
sub_pheno <- sub_pheno %>%
  rename_with(~"Label", matches("condition"))

sub_pheno$Label <- as.factor(sub_pheno$Label) 

all(rownames(sub_pheno) %in% colnames(Raw_data))
all(colnames(Raw_data) %in% rownames(sub_pheno))

# Prepare the Deseq2 list
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = Raw_data,
                              colData = sub_pheno,
                              design = ~ Label)

# Run DESeq to Estimate Size Factors
dds <- DESeq(dds)

#Extract Normalized Counts
vsd <- vst(dds, blind = FALSE)

#Convert to DataFrame
normalized_data <- assay(vsd)

# Read the annotation data
annotation <- read.delim("Human.GRCh38.p13.annot (2).tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sub_annot <- annotation[, c(1,2)]

# Merge the annotataion data to normalized data to map the genes
merged_norm <- merge(normalized_data, sub_annot, by.x = "row.names", by.y = "GeneID", all.x = TRUE)
write.csv(merged_norm, "normalized_with_symbolsGSE241238.csv", row.names = FALSE)

################################################################################
################################################################################
################################################################################

# Search for our original genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)


present_genes <- genes_of_interest[genes_of_interest %in% merged_norm$Symbol]
write.csv(present_genes, "present_genesGSE241238.csv", row.names = TRUE)

missing_genes <- genes_of_interest[!genes_of_interest %in% merged_norm$Symbol]
write.csv(missing_genes, "missing_genesGSE241238.csv", row.names = TRUE)

# Filter our genes

filtered_genes <- merged_norm %>%
  filter(Symbol %in% genes_of_interest)

filtered_genes <- filtered_genes[, -1]

# Transpose the data
transposed_norm <- filtered_genes %>%
  pivot_longer(
    cols = -Symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Symbol,
    values_from = Value
  )

# Merge the phenotype data to identify the sepsis patients
sepsis_data <- merge(transposed_norm, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)
write.csv(sepsis_data, "sepsis_dataGSE241238.csv", row.names = FALSE)

################################################################################
################################  60 Genes. ####################################
################################################################################

# Search for our original genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A",
  "LY96", "CXCR4", "TRIL", "TLR5", "CXCL13"
)


present_genes <- genes_of_interest[genes_of_interest %in% merged_norm$Symbol]
write.csv(present_genes, "present_genesGSE241238.csv", row.names = TRUE)

missing_genes <- genes_of_interest[!genes_of_interest %in% merged_norm$Symbol]
write.csv(missing_genes, "missing_genesGSE241238.csv", row.names = TRUE)

# Filter our genes

filtered_genes <- merged_norm %>%
  filter(Symbol %in% genes_of_interest)

filtered_genes <- filtered_genes[, -1]

# Transpose the data
transposed_norm <- filtered_genes %>%
  pivot_longer(
    cols = -Symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Symbol,
    values_from = Value
  )

# Merge the phenotype data to identify the sepsis patients
sepsis_data_60 <- merge(transposed_norm, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)
write.csv(sepsis_data_60, "sepsis_data_60_GSE241238.csv", row.names = FALSE)

################################################################################
##############################Selected Genes####################################
################################################################################

# Search for our original genes of interest
genes_of_interest <- c(
  "CALCA", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNB1", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A",
  "LY96", "CXCR4", "TLR5"
)


present_genes <- genes_of_interest[genes_of_interest %in% merged_norm$Symbol]
write.csv(present_genes, "present_genesGSE241238.csv", row.names = TRUE)

missing_genes <- genes_of_interest[!genes_of_interest %in% merged_norm$Symbol]
write.csv(missing_genes, "missing_genesGSE241238.csv", row.names = TRUE)

# Filter our genes

filtered_genes <- merged_norm %>%
  filter(Symbol %in% genes_of_interest)

filtered_genes <- filtered_genes[, -1]

# Transpose the data
transposed_norm <- filtered_genes %>%
  pivot_longer(
    cols = -Symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Symbol,
    values_from = Value
  )

# Merge the phenotype data to identify the sepsis patients
sepsis_data_selected <- merge(transposed_norm, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)
write.csv(sepsis_data_selected, "sepsis_data_selected_GSE241238.csv", row.names = FALSE)

