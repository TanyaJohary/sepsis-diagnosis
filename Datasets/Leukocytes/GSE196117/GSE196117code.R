
# GSE196117
# get geo data list
geo_data <- getGEO("GSE196117", GSEMatrix = TRUE)

# Retrieve the expression data 
expression_data <- exprs(geo_data[[1]]) #empty
raw_data <- read.delim("GSE196117_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(raw_data) <- raw_data$GeneID
raw_data <- raw_data[, -1]

# Retrieve the phenotype data of patients
pheno_data <- pData(geo_data[[1]])
sub_pheno_data <- pheno_data[, c(2, 46 )]
sub_pheno_data <- sub_pheno_data %>%
  rename_with(~"Label", matches("condition"))

all(rownames(sub_pheno_data) %in% colnames(raw_data))
all(colnames(raw_data) %in% rownames(sub_pheno_data))

# Create the Deseq object
dds <- DESeqDataSetFromMatrix(
  countData = raw_data,
  colData = sub_pheno_data,
  design =~ Label
)

# amke normalization
# Keep genes with at least 10 counts in total
dds <- dds[ rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds) 

# variance stabilizing
vsd <- vst(dds, blind = FALSE)

normalized_data <- as.data.frame(assay(vsd))

# Retrieve the feature data
feature_data <- fData(geo_data[[1]]) #empty
annotation_data <- read.delim("Human.GRCh38.p13.annot.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sub_annoatation <- annotation_data[, c(1, 2)]


# Merge  expression data with annotation data
merged_normalized <- merge(normalized_data, sub_annoatation, by.x = "row.names", by.y = "GeneID", all.x = TRUE)
write.csv(merged_normalized, "normalizedDeseq_with_symbolsGSE196117.csv", row.names = FALSE)

# search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

present_genes <- genes_of_interest[genes_of_interest %in% merged_normalized$Symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_normalized$Symbol] # "IFNA1" "IFNA2"

write.csv(present_genes, "present_genesGSE196117.csv", row.names = TRUE)
write.csv(missing_genes, "missing_genesGGSE196117.csv", row.names = TRUE)

# Filter our genes of interst
filtered_genes <- merged_normalized %>%
  filter(Symbol %in% genes_of_interest)

# Remove the probe id column
filtered_genes <- filtered_genes[, -1]

# Transpose the data frame
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

# Merge to phenotype data to find the sepsis patients and classify based on sampling day and no-treatment cases
sub_pheno2 <- pheno_data[, c(2, 14, 15, 46)] 
sub_pheno2 <- sub_pheno2 %>%
  rename_with(~"Label", matches("condition"))

merged_transposed <- merge(transposed_data, sub_pheno2, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

################################ Filter the Day1 patients sampling
sepsis_dataD1 <- merged_transposed %>%
  filter(Label == "control" | 
           (Label == "sepsis" & characteristics_ch1.4 == "treatment: Placebo" & characteristics_ch1.5 == "time: treatment day 1"))

sepsis_D1_data <- sepsis_dataD1[, -c(55, 56)]
write.csv(sepsis_D1_data, "sepsis_D1_dataGSE196117.csv", row.names = FALSE)  

################################ Filter the Day5 patients sampling
sepsis_dataD5 <- merged_transposed %>%
  filter(Label == "control" | 
           (Label == "sepsis" & characteristics_ch1.4 == "treatment: Placebo" & characteristics_ch1.5 == "time: treatment day 5"))

sepsis_D5_data <- sepsis_dataD5[, -c(55, 56)]
write.csv(sepsis_D5_data, "sepsis_D5_dataGSE196117.csv", row.names = FALSE) 



