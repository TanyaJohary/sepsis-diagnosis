################################################################################
############################Prepare data for RF#################################
#########################with already normalized TMM############################
################################################################################

# Retrieve the geo data from GEO
geo_data <- getGEO("GSE60424", GSEMatrix = TRUE)

# Retrieve the expression data
expression_data <- exprs(geo_data[[1]]) # empty
feature_data <- fData(geo_data[[1]]) # empty

#Check the phenoData for information about the samples
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c("title", "characteristics_ch1.4")]


#read normalized count data
expression_data <- read.delim("GSE60424_GEOSubmit_FC1to11_normalized_counts.txt", 
                              header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Connect to Ensembl biomart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Query for gene symbols
gene_mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = expression_data$genenames,
                     mart = ensembl)

# Merge annotations with expression data
expression_with_symbolGSE60424 <- merge(expression_data, gene_mapping, 
                                   by.x = "genenames", by.y = "ensembl_gene_id", 
                                   all.x = TRUE)


write.csv(expression_with_symbolGSE60424, "expression_with_symbolGSE60424.csv", row.names = TRUE)

# search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)


present_genes <- genes_of_interest[genes_of_interest %in% expression_with_symbolGSE60424$hgnc_symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% expression_with_symbolGSE60424$hgnc_symbol]

# save the present and missing genes
write.csv(present_genes, "present_genesSE60424.csv", row.names = TRUE, quote = FALSE)
write.csv(missing_genes, "missing_genesSE60424.csv", row.names = TRUE, quote = FALSE)


# Filter genes of interest from expression data
filtered_genes <- expression_with_symbolGSE60424 %>%
  filter(hgnc_symbol %in% genes_of_interest)

# remove the ensemble gene id
filtered_genes <- filtered_genes[, -1]

# Transpose the data to make genes symbols  colnames
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

# Merge the sub_pheno to identify sepsis and cotrols
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "title", all.x = TRUE)

# Filter sepsis patients and healthy controls among other patients
filtered_sepsis <- merged_transposed %>%
  filter(characteristics_ch1.4 %in% c("diseasestatus: Sepsis", "diseasestatus: Healthy Control"))

# Update the characteristics column by simple label
filtered_sepsis <- filtered_sepsis %>%
  mutate(characteristics_ch1.4 = case_when(
    characteristics_ch1.4 == "diseasestatus: Sepsis" ~ "Sepsis",
    characteristics_ch1.4 == "diseasestatus: Healthy Control" ~ "healthy",
    TRUE ~ characteristics_ch1.4  # Keep other values unchanged
  ))

# Rename the target column
sepsis_data <- filtered_sepsis %>%
  rename_with(~"Label", matches("characteristics_ch1.4"))
write.csv(sepsis_data, "sepsis_dataGSE60424.csv", row.names = FALSE)


################################################################################
############################Prepare data for RF#################################
#########################with normalization by Deseq############################
################################################################################

# Read the raw dataset
raw_data <- read.delim("GSE60424_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(raw_data) <- raw_data$GeneID
raw_data <- raw_data[, -1]

# Retrieve the geo data from GEO
geo_data <- getGEO("GSE60424", GSEMatrix = TRUE)

#Check the phenoData for information about the samples and prepare it for coldata in Deseq object
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c("geo_accession", "characteristics_ch1.4")]
sub_pheno <- sub_pheno %>%
 rename_with(~"Label", matches("characteristics_ch1.4"))
sub_pheno <- sub_pheno %>%
  mutate(Label = case_when(
    Label == "diseasestatus: Sepsis" ~ "Sepsis",
    Label == "diseasestatus: Healthy Control" ~ "healthy",
    Label == "diseasestatus: Type 1 Diabetes" ~ "Diabetes",
    Label == "diseasestatus: MS posttreatment" ~ "MS",
    Label == "diseasestatus: MS pretreatment" ~ "MS",
    Label == "diseasestatus: ALS" ~ "ALS",
    TRUE ~ Label
  ))
sub_pheno$Label <- as.factor(sub_pheno$Label)
sub_pheno$Label <- relevel(sub_pheno$Label, ref = "healthy")
all(rownames(sub_pheno) %in% colnames(raw_data))
all(colnames(raw_data) %in% rownames(sub_pheno))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = raw_data,
  colData = sub_pheno,
  design =~ Label
)

# make normalization
dds <- DESeq(dds)

# Apply Variance Stabilizing Transformation (VST)
vsd <- vst(dds, blind = FALSE)

# Get normalized counts
normalized_data <- assay(vsd)

# Remove the gene with very low expression
mean_expression <- rowMeans(normalized_data)
filtered_norm <- normalized_data[mean_expression > 1, ] 

# Read the annotation data for dataset
annotation_data <- read.delim("Human.GRCh38.p13.annot.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sub_annot <- annotation_data[, c(1,2)]


# Merge the normalized data and sub annotation data to map the genes
merged_norm <- merge(filtered_norm, sub_annot, by.x = "row.names", by.y = "GeneID", all.x = TRUE)
write.csv(merged_norm, "normalized_DEseq_with_SymbolGSE60424.csv", row.names = FALSE)

################################################################################
############################## Genes of interest ###############################
################################################################################

# search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

present_genes <- genes_of_interest[genes_of_interest %in% merged_norm$Symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_norm$Symbol]

# Filter our genes of interest
filtered_genes <- merged_norm %>%
  filter(Symbol %in% present_genes)
filtered_genes <- filtered_genes[, -1]

# Transpose the data
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

# Merge to phenotype data to filter the sepsis patients
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Filter sepsis and helathy samples
sepsis_data <- merged_transposed %>%
  filter(Label %in% c("Sepsis", "healthy"))
write.csv(sepsis_data, "sepsis_Deseq_GSE60424.csv", row.names = FALSE)

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

present_genes <- genes_of_interest[genes_of_interest %in% merged_norm$Symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_norm$Symbol]

# Filter our genes of interest
filtered_genes <- merged_norm %>%
  filter(Symbol %in% present_genes)
filtered_genes <- filtered_genes[, -1]

# Transpose the data
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

# Merge to phenotype data to filter the sepsis patients
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Filter sepsis and helathy samples
sepsis_data <- merged_transposed %>%
  filter(Label %in% c("Sepsis", "healthy"))
write.csv(sepsis_data, "sepsis_60Deseq_GSE60424.csv", row.names = FALSE)

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


present_genes <- genes_of_interest[genes_of_interest %in% merged_norm$Symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_norm$Symbol]

# Filter our genes of interest
filtered_genes <- merged_norm %>%
  filter(Symbol %in% present_genes)
filtered_genes <- filtered_genes[, -1]

# Transpose the data
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

# Merge to phenotype data to filter the sepsis patients
merged_transposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Filter sepsis and helathy samples
sepsis_data <- merged_transposed %>%
  filter(Label %in% c("Sepsis", "healthy"))
write.csv(sepsis_data, "sepsis_selected_Deseq_GSE60424.csv", row.names = FALSE)

