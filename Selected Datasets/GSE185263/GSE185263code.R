

################################################################################
############### Deseq2 normalized dataset for random forest ####################
################################################################################
# Retrieve the geo data as a list
geo_data <- getGEO("GSE185263", GSEMatrix = TRUE)

# Retrieve the expression data
expression_data <- exprs(geo_data[[1]]) # empty
raw_count <- read.delim("GSE185263_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
row.names(raw_count) <- raw_count$GeneID
raw_count <- raw_count[, -1]

# Retrieve the phenotype data for dataset
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c(2, 47, 49)]
coldata <- sub_pheno %>%
  rename_with(~"Label", matches("disease"))
# Ensure Label is a factor
coldata$Label <- as.factor(coldata$Label)

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
vsd <- vst(dds, blind = FALSE)

# Extract normalized data
normalized_counts <- assay(vsd)

# Save filtered & normalized data
write.csv(normalized_counts, "normalized_counts_VST_GSE185263.csv", row.names = TRUE)

# Check for Batch Effects
# Perform PCA on normalized counts
pca_data <- prcomp(t(normalized_counts))  # Transpose so samples are rows

# Create a PCA dataframe with metadata
pca_df <- data.frame(
  PC1 = pca_data$x[,1],
  PC2 = pca_data$x[,2],
  Batch = coldata$`collection location:ch1`  
)

# Plot PCA to see batch effect
ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA of Expression Data Before Batch Correction") +
  labs(color = "Collection Location") #Since batch effects exist, apply ComBat from the sva package to correct 

library(sva)

# Ensure collection_location is a factor
coldata$`collection location:ch1` <- as.factor(coldata$`collection location:ch1`)

# Model matrix (preserve biological condition)
mod <- model.matrix(~ Label, data = coldata)

# Apply ComBat batch correction
corrected_counts <- ComBat(
  dat = normalized_counts, 
  batch = coldata$`collection location:ch1`,  
  mod = mod
)

# Save corrected data
write.csv(corrected_counts, "batch_corrected_VST_GSE185263.csv", row.names = TRUE)

# Perform PCA on batch-corrected data
pca_corrected <- prcomp(t(corrected_counts))

# Create a PCA dataframe with metadata
pca_df_corrected <- data.frame(
  PC1 = pca_corrected$x[,1],
  PC2 = pca_corrected$x[,2],
  Batch = coldata$`collection location:ch1`
)

# Plot PCA after batch correction
library(ggplot2)
ggplot(pca_df_corrected, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA After Batch Effect Correction") +
  labs(color = "Collection Location")

corrected_data <- as.data.frame(corrected_counts)

# Read the annotataion data
feature_data <- fData(geo_data[[1]]) # empty
annotation_data <- read.delim("Human.GRCh38.p13.annot.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sub_annot <- annotation_data[, c(1, 2)] 

# Merge annotation data to normalized data to map the genes symbols
corrected_with_symbols <- merge(corrected_data, sub_annot, by.x = "row.names", by.y = "GeneID", all.x = TRUE)
write.csv(corrected_with_symbols, "corrected_with_symbolsGSE185263.csv", row.names = FALSE)

# Search for genes of interest by gene symbols
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
  
)

present_genes <- genes_of_interest[genes_of_interest %in% corrected_with_symbols$Symbol] 
missing_genes <- genes_of_interest[!genes_of_interest %in% corrected_with_symbols$Symbol] # IFNA1, IFNA2 missing

write.csv(present_genes, "present_genesGSE185263.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE185263.csv", row.names = FALSE)

# Filter our genes of interest from normalized data
filtered_genes <- corrected_with_symbols %>%
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
sepsis_data <- sepsis_data[, -63]
write.csv(sepsis_data, "sepsis_dataGSE185263.csv", row.names = FALSE)


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


############################################################# All top genes

top_genes <- c("CD177", "ITGAM", "ARG1", "C3AR1", "FCGR1A",
               "S100A8", "S100A9", "S100A12", "HMGB1",
               "GATA3", "SOCS3", "MAPK14", "IL10", "BCL2",
               "IL1R2", "MYD88", "HIF1A")

present_genes <- top_genes[top_genes %in% corrected_with_symbols$Symbol] 
missing_genes <- top_genes[!top_genes %in% corrected_with_symbols$Symbol]

# Filter our genes of interest from normalized data
filtered_genes <- corrected_with_symbols %>%
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
sepsis_data <- sepsis_data[, -19]
write.csv(sepsis_data, "sepsis_alltops_GSE185263.csv", row.names = FALSE)

############################################################# top cluster 1

cluster1 <- c("CD177", "ITGAM", "C3AR1", "FCGR1A")

# Filter our genes of interest from normalized data
filtered_genes <- corrected_with_symbols %>%
  filter(Symbol %in% cluster1)

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
sepsis_data <- sepsis_data[, -7]
write.csv(sepsis_data, "sepsis_cluster1_GSE185263.csv", row.names = FALSE)

################################################################# cluster2
cluster2 <- c("S100A8", "S100A9", "S100A12", "HMGB1")

# Filter our genes of interest from normalized data
filtered_genes <- corrected_with_symbols %>%
  filter(Symbol %in% cluster2)

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
sepsis_data <- sepsis_data[, -6]
write.csv(sepsis_data, "sepsis_cluster2_GSE185263.csv", row.names = FALSE)

################################################################ cluster 3

cluster3 <- c("GATA3", "SOCS3", "MAPK14", "IL10", "BCL2",
              "IL1R2", "MYD88", "HIF1A", "ARG1")

# Filter our genes of interest from normalized data
filtered_genes <- corrected_with_symbols %>%
  filter(Symbol %in% cluster3)

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
sepsis_data <- sepsis_data[, -10]
write.csv(sepsis_data, "sepsis_cluster3_GSE185263.csv", row.names = FALSE)

##################################################################### last cluster
corrected_with_symbols <- read.csv("corrected_with_symbolsGSE185263.csv")

last_cluster <- c("CD177", "S100A8", "S100A9", "GATA3", "S100A12", "ARG1","ITGAM",
                   "SOCS3", "C3AR1", "IL1R2", "IL10", "FCGR1A", "MMP8", "MAPK14")


present_genes <- last_cluster[last_cluster %in% corrected_with_symbols$Symbol] 
missing_genes <- last_cluster[!last_cluster %in% corrected_with_symbols$Symbol]

# Filter our genes of interest from normalized data
filtered_genes <- corrected_with_symbols %>%
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
sepsis_data <- sepsis_data[, -16]
write.csv(sepsis_data, "sepsis_last_GSE185263.csv", row.names = FALSE)



