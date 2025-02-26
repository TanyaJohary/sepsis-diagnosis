
# GSE133822
geo_data <- getGEO("GSE133822", GSEMatrix = TRUE)


# Get genes expression data
expression_data <- exprs(geo_data[[1]]) #empty
raw_data <- read.delim("GSE133822_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(raw_data) <- raw_data$GeneID
raw_data <- raw_data[, -1]

# Get phenotype data
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c(2, 45)]
sub_pheno <- sub_pheno %>%
  rename_with(~"Label", matches("patient"))

# Prepare for Deseq2 noemalization
all(rownames(sub_pheno) %in% colnames(raw_data)) # FALSE
setdiff(rownames(sub_pheno), colnames(raw_data)) # lot of samples
setdiff(colnames(raw_data), rownames(sub_pheno)) # zero

# Find the common samples between sub_pheno and raw_data
common_samples <- intersect(rownames(sub_pheno), colnames(raw_data))

# Subset sub_pheno and raw_data to only include common samples
sub_pheno <- sub_pheno[common_samples, ]
raw_data <- raw_data[, common_samples]

all(rownames(sub_pheno) %in% colnames(raw_data)) # TRUE
all(colnames(raw_data) %in% rownames(sub_pheno)) # TRUE

# create Deseq object
dds <- DESeqDataSetFromMatrix(countData = raw_data,
                              colData = sub_pheno,
                              design = ~ 1)

# Normalization
dds <- DESeq(dds)

# Apply a Variance Stabilizing Transformation (VST)
# Using blind = TRUE prevents the transformation from using the experimental design.
vsd <- vst(dds, blind = TRUE)

norm_data <- assay(vsd)

# load the annotation data to map the genes symbol
annotation_data <- read.delim("Human.GRCh38.p13.annot.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sub_annot <- annotation_data[, c(1, 2)]

# Merge the anootation data and normalized expression data
norm_with_symbol <- merge(norm_data, sub_annot, by.x = "row.names", by.y = "GeneID", all.x = TRUE)
write.csv(norm_with_symbol, "normalized_with_symbolGSE133822.csv", row.names = FALSE)

# Search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

present_genes <- genes_of_interest[genes_of_interest %in% norm_with_symbol$Symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% norm_with_symbol$Symbol]

# save the presents and missging genes of interest
write.csv(present_genes, "present_genesGSE133822.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE133822.csv", row.names = FALSE)


# filter out the goal genes expression data
filtered_genes <- norm_with_symbol %>%
  filter(Symbol %in% genes_of_interest)

# Remove the probe id column
filtered_genes <- filtered_genes[, -1]

# Treanspose data
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


# merge the transposed data with sub pheno data
merged_transoped <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Filter out the sepsis and healthy ones
filtered_sepsis <- merged_transoped %>%
  filter(Label %in% c( "Sepsis", "Healthy"))
write.csv(filtered_sepsis, "sepsis_dataGSE133822.csv", row.names = FALSE)

################################################################################
################################### 60 Genes  ##################################
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

present_genes <- genes_of_interest[genes_of_interest %in% norm_with_symbol$Symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% norm_with_symbol$Symbol]

# save the presents and missging genes of interest
write.csv(present_genes, "present_genesGSE133822.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE133822.csv", row.names = FALSE)


# filter out the goal genes expression data
filtered_genes <- norm_with_symbol %>%
  filter(Symbol %in% genes_of_interest)

# Remove the probe id column
filtered_genes <- filtered_genes[, -1]

# Treanspose data
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


# merge the transposed data with sub pheno data
merged_transoped <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Filter out the sepsis and healthy ones
filtered_sepsis <- merged_transoped %>%
  filter(Label %in% c( "Sepsis", "Healthy"))
write.csv(filtered_sepsis, "sepsis_60_dataGSE133822.csv", row.names = FALSE)



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


present_genes <- genes_of_interest[genes_of_interest %in% norm_with_symbol$Symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% norm_with_symbol$Symbol]

# save the presents and missging genes of interest
write.csv(present_genes, "present_genesGSE133822.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE133822.csv", row.names = FALSE)


# filter out the goal genes expression data
filtered_genes <- norm_with_symbol %>%
  filter(Symbol %in% genes_of_interest)

# Remove the probe id column
filtered_genes <- filtered_genes[, -1]

# Treanspose data
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


# merge the transposed data with sub pheno data
merged_transoped <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Filter out the sepsis and healthy ones
filtered_sepsis <- merged_transoped %>%
  filter(Label %in% c( "Sepsis", "Healthy"))
write.csv(filtered_sepsis, "sepsis_selected_dataGSE133822.csv", row.names = FALSE)
