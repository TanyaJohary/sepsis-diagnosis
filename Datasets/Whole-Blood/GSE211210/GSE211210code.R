# Retrieve the geo_data list
geo_data <- getGEO("GSE211210", GSEMatrix = TRUE) 

expression_data <- exprs(geo_data[[1]]). # empty
feature_data <- fData(geo_data[[1]]).  # empty
pheno_data <- pData(geo_data[[1]])

# Load the supplementary data
supplementary_data <- read.delim("GSE211210_gene_fpkm_all_samples.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Connect to Ensembl's BioMart database
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Use the getBM function to fetch the gene symbols for the given Ensembl IDs
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = supplementary_data$Gene.ID,
  mart = ensembl
)


# Merge the retrieved gene symbols with the supplementary data
merged_data <- merge(supplementary_data, gene_mapping, by.x = "Gene.ID", by.y = "ensembl_gene_id", all.x = TRUE)

write.csv(annotated_data, "Annotated_Supplementary_DataGSE211210.csv", row.names = FALSE)

# search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)


present_genes <- genes_of_interest[genes_of_interest %in% merged_data$hgnc_symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_data$hgnc_symbol]

# save the present and missing genes
write.csv(present_genes, "present_genesGSE211210.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE211210.csv", row.names = FALSE)


# Filter our genes of interest
filtered_genes <- merged_data %>%
  filter(hgnc_symbol %in% genes_of_interest)

# Remove the gene id column
filtered_genes <- filtered_genes[, -1]

# transpose the data
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

# Merge sub_pheno data to expression data
sub_pheno <- pheno_data[, c(1, 43)]
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "title", all.x = TRUE)

# Define the target column as label
merged_transposed <- merged_transposed %>%
  rename_with(~"Label", matches("subject"))

# Simplify the levels
sepsis_data <- merged_transposed %>%
  mutate(Label = case_when(
    grepl("patients", Label) ~ "Sepsis",
    grepl("health", Label) ~ "healthy",
    TRUE ~ Label
  ))

write.csv(sepsis_data, "sepsis_dataGSE211210.csv", row.names = FALSE)

################################################################################
##################### Use Raw data to normalize by Deseq2#######################
################################################################################

# Download supplementary files (often contains raw counts)
raw_data <- read.delim("GSE211210_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, sep = "\t")

# Retrive the meta data
geo_data <- getGEO("GSE211210", GSEMatrix = TRUE)
feature_data <- fData(geo_data[[1]])
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c(1,2)]

# prepare the phenodata for Deseq2
sub_pheno <- sub_pheno %>%
  mutate(title= case_when(
    grepl("CON", title) ~ "Control",
    grepl("SEP", title) ~ "Sepsis"
  ))

# Prepare raw data for deseq2
rownames(raw_data) <- raw_data$GeneID
raw_data <- raw_data[, -1]

# DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(
  countData = raw_data,
  colData = sub_pheno,
  design = ~title
)

# Perform Deseq normalization
dds <- DESeq(dds)

# Obtain normalized count
count_normalized <- counts(dds , normalized = TRUE)
write.csv(count_normalized, "normalized_counts_Deseq2GSE211210.csv", row.names = FALSE)

# Load the annotation file
annotation_data <- read.delim("Human.GRCh38.p13.annot.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sub_annot <- annotation_data[, c(1,2)]

# Merge annotation data to normalized data to map teh genes
merged_normalized <- merge(count_normalized, sub_annot, by.x = "row.names", by.y = "GeneID", all.x = TRUE)
write.csv(merged_normalized, "merged_normalized_Deseq2GSE211210.csv", row.names = FALSE)

# Search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)


present_genes <- genes_of_interest[genes_of_interest %in% merged_normalized$Symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_normalized$Symbol]

# save the present and missing genes
write.csv(present_genes, "present_genesGSE211210.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE211210.csv", row.names = FALSE)

# Filter our genes of interest
filtered_norm <- merged_normalized %>%
  filter(Symbol %in% genes_of_interest)

# Transpose the data
filtered_norm <- filtered_norm[,-1]
transposed_norm <- filtered_norm %>%
  pivot_longer(
    cols = -Symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Symbol,
    values_from = Value
  )

# Merge to the sub pheno data to identify sepsis and healthy controls
sepsis_normalized <- merge(transposed_norm, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Define the target column as Label
sepsis_normalized <- sepsis_normalized %>%
  rename_with(~"Label", matches("title"))
write.csv(sepsis_normalized, "sepsis_Deseq_normGSE211210.csv", row.names = FALSE)  

