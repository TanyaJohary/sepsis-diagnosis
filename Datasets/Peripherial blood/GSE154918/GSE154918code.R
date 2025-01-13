
# Retrieve the geo data
geo_data <- getGEO("GSE154918", GSEMatrix = TRUE)

# Retrieve the phenotype data
pheno_data <- pData(geo_data[[1]])
colnames(pheno_data)
sub_pheno <- pheno_data[, c("title", "geo_accession", "characteristics_ch1", "status:ch1")]
write.csv(sub_pheno, "sub_phenoGSE154918.csv", row.names = TRUE)

# Retrieve the feature data of genes and probes
feature_data <- fData(geo_data[[1]])
colnames(feature_data) # was empty

# Retrieve expression data of GSE154918
expression_data <- exprs(geo_data[[1]]) # also empty expression data

# Download the expression file from NCBI_GEO and read the file
expression <- read.table("GSE154918_Schughart_Sepsis_200320.txt", header = TRUE, sep = " ", stringsAsFactors = FALSE, fill = TRUE) # file include the gene symbol
write.csv(expression, "cleaned_expression_dataGSE154918.csv", row.names = FALSE)

# search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)


present_genes <- genes_of_interest[genes_of_interest %in% expression$gene_symbol]
write.csv(present_genes, "present_genesGSE154918.csv", row.names = TRUE)

missing_genes <- genes_of_interest[!genes_of_interest %in% expression$gene_symbol]
write.csv(missing_genes, "missing_genesGSE154918.csv", row.names = TRUE)

# Reload the dataset and sub_phenotype data to filter the genes of interest expression data
expression <- read.csv("cleaned_expression_dataGSE154918.csv")
sub_pheno<- read.csv("sub_phenoGSE154918.csv")

# Remove the extra columns that we do not need anymore
expression <- dplyr::select(expression, -c(Row.names, ENSEMBL_gene_ID, ensembl_transcript_id, entrezgene_id, uniprotswissprot, description, refseq_mrna, transcript_length))

# filter our genes of interest from expression data
filtered_genes <- expression %>%
  filter( gene_symbol %in% genes_of_interest)

# transpose the data to make gene symbol as columns name
transposed_data <- filtered_genes %>%
  pivot_longer(
    cols = -gene_symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = gene_symbol,
    values_from = Value
  )

# dataset include expression changes for heathy control (Hlty), uncomplicated infection (Inf1_P),
# sepsis (Seps_P), septic shock (Shock_P), follow-up of sepsis (Seps_FU), follow-up of septic shock (Shock_FU) groups.
# should exclude the Inf1_P, Seps_FU and Shock_FU group from dataset

# delete extra columns in sub pheno data to merge with expression data
sub_pheno <- dplyr::select(sub_pheno, -c(X, geo_accession, characteristics_ch1))

# Merge the transposed data with phenotype data to identify the sepsis patients
merged_transposed_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "title", all.x = TRUE)

###### Filter the sepsis patients septic shock patients and healthy controls#######
filtered_sepsis_shock <- merged_transposed_data %>%
  filter(status.ch1 %in% c("Hlty", "Seps_P", "Shock_P"))

# clarify the label values
filtered_sepsis_shock <- filtered_sepsis_shock %>%
  mutate(status.ch1 = case_when(
    grepl("Hlty", status.ch1) ~ "Healthy control",
    grepl("Seps_P", status.ch1) ~ "Sepsis",
    grepl("Shock_P", status.ch1) ~ "Septic shock"
  ))
# Define the label column
sepsis_shock_labeled <- filtered_sepsis_shock %>%
  rename_with(~ "Label", matches("status.ch1"))
# Save the dataset contains the sepsis, septic shock patients and healthy controls
write.csv(sepsis_shock_labeled, "sepsis_shock_labeledGSE154918.csv", row.names = FALSE)

###### Filter only the sepsis patients  and healthy controls#######
filtered_sepsis_only <- merged_transposed_data %>%
  filter(status.ch1 %in% c("Hlty", "Seps_P"))

# clarify the label values
filtered_sepsis_only <- filtered_sepsis_only %>%
  mutate(status.ch1 = case_when(
    grepl("Hlty", status.ch1) ~ "Healthy control",
    grepl("Seps_P", status.ch1) ~ "Sepsis"
  ))

# Define the label column
sepsis_only_labeled <- filtered_sepsis_only %>%
  rename_with(~ "Label", matches("status.ch1"))

# Save the dataset contains the sepsis, septic shock patients and healthy controls
write.csv(sepsis_only_labeled, "sepsis_only_labeledGSE154918.csv", row.names = FALSE)

# Make sure it has been correctly saved
x <- read.csv("sepsis_only_labeledGSE154918.csv")
