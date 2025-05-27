
# retrieve the geo data as alist
geo_data <- getGEO("GSE65682", GSEMatrix = TRUE)

# retrieve the expression data
expression_data <- exprs(geo_data[[1]])

# retrieve the phenotype data
pheno_data <- pData(geo_data[[1]])

pheno_data$title <- gsub("\\[.*\\]", "", pheno_data$title) # Remove the numbers inside square brackets from the 'title' column
pheno_data$title <- trimws(pheno_data$title) # Trim any extra spaces that may be left behind

sub_pheno <- pheno_data[, c(1, 2, 12)]
sub_pheno$characteristics_ch1.2 <- gsub("pneumonia diagnoses:", "", sub_pheno$characteristics_ch1.2)
sub_pheno <- sub_pheno %>% 
  rename_with(~"Label", matches("characteristics"))

# Ensure NAs are properly recognized
sub_pheno$Label <- as.character(sub_pheno$Label)  # Convert factor to character
sub_pheno$title <- as.character(sub_pheno$title)  # Convert title column to character

# Clean up whitespace and hidden characters
sub_pheno$Label <- trimws(sub_pheno$Label)            # Remove leading/trailing spaces
sub_pheno$Label <- gsub("[\r\n]", "", sub_pheno$Label)  # Remove newline characters

# Now replace any "NA" string with actual NA
sub_pheno$Label[sub_pheno$Label == "NA"] <- NA
table(sub_pheno$Label, useNA = "always") #see "NA" disappear and <NA> increase.

sub_pheno$Label[which(sub_pheno$title == "healthy subject" & is.na(sub_pheno$Label))] <- "healthy"

sub_pheno <- sub_pheno %>% 
  mutate(Label = case_when(
    Label == "cap" ~ "Sepsis",
    Label == "hap" ~ "Sepsis",
    TRUE ~ Label
  ))

sub_pheno <- sub_pheno[, -1] # remove the title column
write.csv(sub_pheno, "sub_phenoGSE65682.csv", row.names = FALSE)

# retrieve the feature data
feature_data <- fData(geo_data[[1]])
sub_feature <- feature_data[, c(1, 15)]

# merge the expression data and feature data to map the genes symbols
expression_with_symbols <- merge(expression_data, sub_feature, by.x = "row.names", by.y = "ID", all.x = TRUE)

# our genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "IL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

# search for our genes of interest
present_genes <- genes_of_interest[genes_of_interest %in% expression_with_symbols$`Gene Symbol`]
missing_genes <- genes_of_interest[!genes_of_interest %in% expression_with_symbols$`Gene Symbol`]

# save the results
write.csv(present_genes, "present_genesGSE65682.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE65682.csv", row.names = FALSE)

# filter our genes of interest
filtered_genes <- expression_with_symbols %>% 
  filter(`Gene Symbol` %in% present_genes)

filtered_genes <- filtered_genes[, -1]
filtered_genes$`Gene Symbol` <- make.unique(as.character(filtered_genes$`Gene Symbol`))

# Transpose the data
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = `Gene Symbol`,
    values_from = Value
  )

# Merge to sub pheno data to identify the sepsis patients and healthy ones
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# filter our sepsis and healthy groups
sepsis_data <- merged_data %>% 
  filter(Label %in% c("Sepsis", "healthy"))
write.csv(sepsis_data, "sepsis_dataGSE65682.csv", row.names = FALSE)

############################################################## top genes

top_genes <- c("CD177", "ITGAM", "ARG1", "C3AR1", "FCGR1A",
               "S100A8", "S100A9", "S100A12", "HMGB1",
               "GATA3", "SOCS3", "MAPK14", "IL10", "BCL2",
               "IL1R2", "MYD88", "HIF1A")
 
present_top <- top_genes[top_genes %in% expression_with_symbols$`Gene Symbol`]
missing_top <- top_genes[!top_genes %in% expression_with_symbols$`Gene Symbol`]

# filter our genes of interest
filtered_genes <- expression_with_symbols %>% 
  filter(`Gene Symbol` %in% present_top)

filtered_genes <- filtered_genes[, -1]
filtered_genes$`Gene Symbol` <- make.unique(as.character(filtered_genes$`Gene Symbol`))

# Transpose the data
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = `Gene Symbol`,
    values_from = Value
  )

# Merge to sub pheno data to identify the sepsis patients and healthy ones
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# filter our sepsis and healthy groups
sepsis_data <- merged_data %>% 
  filter(Label %in% c("Sepsis", "healthy"))
write.csv(sepsis_data, "sepsis_topGSE65682.csv", row.names = FALSE)

############################################################# top cluster 1

cluster1 <- c("CD177", "ITGAM", "ARG1", "C3AR1", "FCGR1A")

# filter our genes of interest
filtered_genes <- expression_with_symbols %>% 
  filter(`Gene Symbol` %in% cluster1)

filtered_genes <- filtered_genes[, -1]
filtered_genes$`Gene Symbol` <- make.unique(as.character(filtered_genes$`Gene Symbol`))

# Transpose the data
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = `Gene Symbol`,
    values_from = Value
  )

# Merge to sub pheno data to identify the sepsis patients and healthy ones
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# filter our sepsis and healthy groups
sepsis_data <- merged_data %>% 
  filter(Label %in% c("Sepsis", "healthy"))
write.csv(sepsis_data, "sepsis_cluster1_GSE65682.csv", row.names = FALSE)

############################################################ cluster2
cluster2 <- c("S100A8", "S100A9", "S100A12", "HMGB1")

# filter our genes of interest
filtered_genes <- expression_with_symbols %>% 
  filter(`Gene Symbol` %in% cluster2)

filtered_genes <- filtered_genes[, -1]
filtered_genes$`Gene Symbol` <- make.unique(as.character(filtered_genes$`Gene Symbol`))

# Transpose the data
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = `Gene Symbol`,
    values_from = Value
  )

# Merge to sub pheno data to identify the sepsis patients and healthy ones
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# filter our sepsis and healthy groups
sepsis_data <- merged_data %>% 
  filter(Label %in% c("Sepsis", "healthy"))
write.csv(sepsis_data, "sepsis_cluster2_GSE65682.csv", row.names = FALSE)

############################################################# cluster 3

cluster3 <- c("GATA3", "SOCS3", "MAPK14", "IL10", "BCL2",
              "IL1R2", "MYD88", "HIF1A")

# filter our genes of interest
filtered_genes <- expression_with_symbols %>% 
  filter(`Gene Symbol` %in% cluster3)

filtered_genes <- filtered_genes[, -1]
filtered_genes$`Gene Symbol` <- make.unique(as.character(filtered_genes$`Gene Symbol`))

# Transpose the data
transposed_data <- filtered_genes %>% 
  pivot_longer(
    cols = -`Gene Symbol`,
    names_to = "Sample",
    values_to = "Value"
  ) %>% 
  pivot_wider(
    names_from = `Gene Symbol`,
    values_from = Value
  )

# Merge to sub pheno data to identify the sepsis patients and healthy ones
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# filter our sepsis and healthy groups
sepsis_data <- merged_data %>% 
  filter(Label %in% c("Sepsis", "healthy"))
write.csv(sepsis_data, "sepsis_cluster3_GSE65682.csv", row.names = FALSE)
