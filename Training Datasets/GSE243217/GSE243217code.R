
#GSE243217

# Retrieve geo data as list
geo_data <- getGEO("GSE243217", GSEMatrix = TRUE)

# Retrieve the phenotype data
pheno_data <- pData(geo_data[[1]])
colnames(pheno_data)
sub_pheno <- pheno_data[, c("geo_accession", "disease state:ch1")]
sub_pheno <- sub_pheno%>%
  rename_with(~"Label", matches("disease"))

sub_pheno$Label <- gsub("[^a-zA-Z0-9_.]", "_", sub_pheno$Label)  # Replace spaces/special chars with '_'
sub_pheno$Label <- as.factor(sub_pheno$Label)  # Ensure it's a factor

# Retrieve the feature data for gene symbols
feature_data <- fData(geo_data[[1]]) # empty

# Retrieve the expression data
expression_data <- exprs(geo_data[[1]]) # empty

# Read Raw count data as expression file from NCBI_GEO 
expression_data <- read.delim("GSE243217_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE) # does not have all the samples

# Retrieve the raw data
getGEOSuppFiles("GSE243217")

# Extract the raw data
untar("GSE243217/GSE243217_RAW.tar", exdir = "GSE243217_RAW")
system("gunzip GSE243217_RAW/*.gz")

# Load required libraries
library(dplyr)
library(tidyr)

# Define the directory containing the raw count files
data_dir <- "GSE243217_RAW/"

# List all text files in the directory
files <- list.files(data_dir, pattern = "\\.txt$", full.names = TRUE)

# Read and merge all count files
count_list <- lapply(files, function(f) {
  read.table(f, header = TRUE, row.names = 1, sep = "\t")
})

# Merge all individual sample files into one matrix
# Read all count files into a list, ensuring they have the same gene order
count_list <- lapply(files, function(f) {
  df <- fread(f, header = TRUE, sep = "\t")
  setnames(df, old = names(df)[2], new = gsub("_raw_count.txt", "", basename(f)))  # Rename sample column
  df
})


# Merge all files using 'Reduce' while preserving row names
# Merge all count tables while preserving gene symbols
merged_counts <- Reduce(function(x, y) merge(x, y, by = "gene_symbol", all = TRUE), count_list)

# Convert merged_counts to a standard data.frame
merged_counts_df <- as.data.frame(merged_counts)

# Set gene symbols as row names
rownames(merged_counts_df) <- merged_counts_df$gene_symbol

# Remove the redundant 'gene_symbol' column
merged_counts_df <- dplyr::select(merged_counts_df, -gene_symbol)

# Remove everything after the first underscore (_) in column names
colnames(merged_counts_df) <- gsub("_.*", "", colnames(merged_counts_df))
all(colnames(merged_counts_df) %in% rownames(sub_pheno))
all(rownames(sub_pheno) %in% colnames(merged_counts_df))

# dataset has decimal values, should round the values
merged_counts_df <- round(merged_counts_df)  # Convert to integer values

# NOrmalizing by Deseq2
dds <- DESeqDataSetFromMatrix(countData = merged_counts_df,
                              colData = sub_pheno,
                              design = ~ Label)

# Keep genes with at least 10 counts in total
dds <- dds[ rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds) 

# Normalize counts using variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)  

# Convert to a data frame
normalized_data <- as.data.frame(assay(vsd))
normalized_data$gene_symbol <- rownames(normalized_data)
write.csv(normalized_data, "normalized_with_symbolsGSE243217.csv", row.names = FALSE)

################################################################################
################################################################################
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


present_genes <- genes_of_interest[genes_of_interest %in% normalized_data$gene_symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% normalized_data$gene_symbol] # IFNB1 is missing

# save the present and missing genes for one sample
write.csv(present_genes, "present_genesGSE243217.csv", row.names = TRUE, quote = FALSE)
write.csv(missing_genes, "missing_genesGSE243217.csv", row.names = TRUE, quote = FALSE)

# Filter our genes of interest
fileterd_genes <- normalized_data %>%
  filter(gene_symbol %in% genes_of_interest)

# Transpose the data to make gene symbols as header
transposed_data <- fileterd_genes %>%
  pivot_longer(
    cols = -gene_symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = gene_symbol,
    values_from = Value
  )
# Merge to phenodata to identify sepsis patients
merged_tarnsposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Filter sepsis patients and healthy controls from Covid-19
sepsis_data <- merged_tarnsposed %>%
  filter(Label != "COVID_19")

# Rename the target columns value
sepsis_data <- sepsis_data %>%
  mutate(Label = case_when(
    Label == "healthy_donor" ~ "healthy",
    TRUE ~ Label
  ))

write.csv(sepsis_data, "sepsis_dataGSE243217.csv", row.names = FALSE)

################################################################################
################################# 60 Genes #####################################
################################################################################
# search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A",
  "LY96", "CXCR4", "TRIL", "TLR5", "CXCL13"
)


present_genes <- genes_of_interest[genes_of_interest %in% normalized_data$gene_symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% normalized_data$gene_symbol] # IFNB1 is missing

# save the present and missing genes for one sample
write.csv(present_genes, "present_genesGSE243217.csv", row.names = TRUE, quote = FALSE)
write.csv(missing_genes, "missing_genesGSE243217.csv", row.names = TRUE, quote = FALSE)

# Filter our genes of interest
fileterd_genes <- normalized_data %>%
  filter(gene_symbol %in% genes_of_interest)

# Transpose the data to make gene symbols as header
transposed_data <- fileterd_genes %>%
  pivot_longer(
    cols = -gene_symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = gene_symbol,
    values_from = Value
  )
# Merge to phenodata to identify sepsis patients
merged_tarnsposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Filter sepsis patients and healthy controls from Covid-19
sepsis_data <- merged_tarnsposed %>%
  filter(Label != "COVID_19")

# Rename the target columns value
sepsis_data <- sepsis_data %>%
  mutate(Label = case_when(
    Label == "healthy_donor" ~ "healthy",
    TRUE ~ Label
  ))

write.csv(sepsis_data, "sepsis_60_dataGSE243217.csv", row.names = FALSE)

################################################################################
############################ selected Genes ####################################
################################################################################
# search for genes of interest
genes_of_interest <- c(
  "CALCA", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNB1", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A",
  "LY96", "CXCR4", "TLR5"
)



present_genes <- genes_of_interest[genes_of_interest %in% normalized_data$gene_symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% normalized_data$gene_symbol] # IFNB1 is missing

# save the present and missing genes for one sample
write.csv(present_genes, "present_genesGSE243217.csv", row.names = TRUE, quote = FALSE)
write.csv(missing_genes, "missing_genesGSE243217.csv", row.names = TRUE, quote = FALSE)

# Filter our genes of interest
fileterd_genes <- normalized_data %>%
  filter(gene_symbol %in% genes_of_interest)

# Transpose the data to make gene symbols as header
transposed_data <- fileterd_genes %>%
  pivot_longer(
    cols = -gene_symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = gene_symbol,
    values_from = Value
  )
# Merge to phenodata to identify sepsis patients
merged_tarnsposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Filter sepsis patients and healthy controls from Covid-19
sepsis_selected <- merged_tarnsposed %>%
  filter(Label != "COVID_19")

# Rename the target columns value
sepsis_selected <- sepsis_selected %>%
  mutate(Label = case_when(
    Label == "healthy_donor" ~ "healthy",
    TRUE ~ Label
  ))

write.csv(sepsis_selected, "sepsis_selected_dataGSE243217.csv", row.names = FALSE)

##################################################################### last cluster
normalized_data <- read.csv("normalized_with_symbolsGSE243217.csv")

last_cluster <- c("CD177", "S100A8", "S100A9", "GATA3", "S100A12", "ARG1","ITGAM",
                  "SOCS3", "C3AR1", "IL1R2", "IL10", "FCGR1A", "MMP8", "MAPK14")


present_genes <- last_cluster[last_cluster %in% normalized_data$gene_symbol]
missing_genes <- last_cluster[!last_cluster %in% normalized_data$gene_symbol] # IFNB1 is missing

# save the present and missing genes for one sample
write.csv(present_genes, "present_genesGSE243217.csv", row.names = TRUE, quote = FALSE)
write.csv(missing_genes, "missing_genesGSE243217.csv", row.names = TRUE, quote = FALSE)

# Filter our genes of interest
fileterd_genes <- normalized_data %>%
  filter(gene_symbol %in% last_cluster)

# Transpose the data to make gene symbols as header
transposed_data <- fileterd_genes %>%
  pivot_longer(
    cols = -gene_symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = gene_symbol,
    values_from = Value
  )
# Merge to phenodata to identify sepsis patients
merged_tarnsposed <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Filter sepsis patients and healthy controls from Covid-19
sepsis_selected <- merged_tarnsposed %>%
  filter(Label != "COVID_19")

# Rename the target columns value
sepsis_selected <- sepsis_selected %>%
  mutate(Label = case_when(
    Label == "healthy_donor" ~ "healthy",
    TRUE ~ Label
  ))

write.csv(sepsis_selected, "sepsis_lasttop_GSE243217.csv", row.names = FALSE)
