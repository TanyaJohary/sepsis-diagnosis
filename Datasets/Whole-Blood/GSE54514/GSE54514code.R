
####################prepare data for Random-forest##########################
############################################################################

getGEOSuppFiles("GSE54514")  # Retrieve the raw data for normaliztaion

untar("GSE54514/GSE54514_RAW.tar", exdir= 'data/')

raw_data <- ReadAffy(celfile.path = "/home/rstudio/sepsis/data/") # Raw data is not available


################################################################################
#########################prepare dataset for RF#################################
################################################################################

# Retrieve the datasets list
geo_data <- getGEO("GSE54514", GSEMatrix = TRUE)

# Retrieve the expression data
expression_data <- exprs(geo_data[[1]])
seriesmatrix <- read.delim("GSE54514_series_matrix.txt", header = TRUE, stringsAsFactors = FALSE, skip = 71)
non_norm <- read.delim("GSE54514_non-normalized.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

###### this dataset is already normalized based on related publication

# Retrieve the phenotype data
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c(2,43)]
write.csv(pheno_data, "pheno_data.csv", row.names = FALSE)


# Retrieve the annotation data for gene symbol
featuere_data <- fData(geo_data[[1]])
sub_feature <- featuere_data[, c(1, 14)]

# Merge the expression data with gene symbol by probe ids
merged_data <- merge(expression_data, sub_feature, by.x = "row.names", by.y = "ID", all.x = TRUE)
write.csv(merged_data, "expression_with_symbolGSE54514.csv", row.names = FALSE)

# saerch for genes of interest based on gene symbol
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELA2", "IL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

present_genes <- genes_of_interest[genes_of_interest %in% merged_data$Symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_data$Symbol]

# search for aliases genes
alias_genes <- c(
  "IFNA1", "IFNA2", "IFNB1",
  "ELANE", "HLE", "PMN-E", "ELA2", "HNE", "NE",
  "MDNCF", "NAP-1", "GCP-1", "SCYB8", "LYNAP", "IL-8", "IL8",
  "CGRP-Alpha", "Calcitonin", "CALC1", "CGRP", "PCT", "CGRP1",
  "PTX1", "C-Reactive Protein", "Pentraxin 1",
  "BPIFD2", "Lipopolysaccharide Binding Protein",
  "PTX3", "TSG-14", "TNFAIP5", "TSG14"
)

present_alias_genes <- alias_genes[alias_genes %in% merged_data$Symbol]

write.csv(present_genes, "present_genesGSE54514.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE54514.csv", row.names = FALSE)

# Filter genes of interest from dataset
filtered_genes <- merged_data %>%
  filter(Symbol %in% genes_of_interest)

# Remove the probe id column
filtered_genes <- filtered_genes[, -1]

# Because some genes have duplicated columns we have to make unique each value of column
filtered_genes$Symbol <- make.unique(as.character(filtered_genes$Symbol))

# Transpose the dataset
transposed_data <- filtered_genes %>%
  pivot_longer(
    cols = - Symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Symbol,
    values_from = Value
  )

# REtrieve the pheno data
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c(2,43)]

# Merged phenotype data to expression data to identify sepsis and healthy cases
transposed_merged <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Define the target column
transposed_merged <- transposed_merged %>%
  rename_with(~"Label", matches("disease"))

# Make a diagnostic dataset- sepsis patients all days, D1-D5
sepsis_diagnosis_data <- transposed_merged %>%
  mutate(Label = case_when(
    Label == "sepsis nonsurvivor" ~ "Sepsis",
    Label == "sepsis survivor" ~ "Sepsis",
    TRUE ~ Label
  ))

write.csv(sepsis_diagnosis_data, "sepsis_diagnosis_dataGSE54514.csv", row.names = FALSE)

# Make prognosis dataset- sepsis patients all days, D1-D5
sepsis_prognosis_data <- transposed_merged %>%
  mutate(Label = case_when(
    Label == "sepsis nonsurvivor" ~ "Nonsurvivor",
    Label == "sepsis survivor" ~ "Survivor",
    TRUE ~ Label
  ))

write.csv(sepsis_prognosis_data, "sepsis_prognosis_dataGSE54514.csv", row.names = FALSE)

######################################## Day1 ##################################

# REtrieve the pheno data
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c(2,45)]

# Merged phenotype data to expression data to identify sepsis and healthy cases
transposed_merged <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Define the target column
transposed_merged <- transposed_merged %>%
  rename_with(~"Label", matches("group"))

# Filter sepsis patients for D1 and healthy one 
sepsis_D1 <- transposed_merged %>%
  filter(Label %in% c("HC_D1", "HC_D5", "NS_D1", "S_D1")) %>%
  mutate(Label = case_when(
    Label %in% c("HC_D1", "HC_D5") ~ "Control",
    Label %in% c("NS_D1", "S_D1") ~ "Sepsis"
  ))

write.csv(sepsis_D1, "sepsis_D1_GSE54514.csv", row.names = FALSE)

############################## Day2 ############################################

# Filter sepsis patients for D2 and healthy one 
sepsis_D2 <- transposed_merged %>%
  filter(Label %in% c("HC_D1", "HC_D5", "NS_D2", "S_D2")) %>%
  mutate(Label = case_when(
    Label %in% c("HC_D1", "HC_D5") ~ "Control",
    Label %in% c("NS_D2", "S_D2") ~ "Sepsis"
  ))

write.csv(sepsis_D2, "sepsis_D2_GSE54514.csv", row.names = FALSE)

################################# Day3 #########################################
# Filter sepsis patients for D3 and healthy one 
sepsis_D3 <- transposed_merged %>%
  filter(Label %in% c("HC_D1", "HC_D5", "NS_D3", "S_D3")) %>%
  mutate(Label = case_when(
    Label %in% c("HC_D1", "HC_D5") ~ "Control",
    Label %in% c("NS_D3", "S_D3") ~ "Sepsis"
  ))

write.csv(sepsis_D3, "sepsis_D3_GSE54514.csv", row.names = FALSE)

################################ Day4 ##########################################
# Filter sepsis patients for D4 and healthy one 
sepsis_D4 <- transposed_merged %>%
  filter(Label %in% c("HC_D5", "NS_D4", "S_D4")) %>%
  mutate(Label = case_when(
    Label %in% c("HC_D5") ~ "Control",
    Label %in% c("NS_D4", "S_D4") ~ "Sepsis"
  ))

write.csv(sepsis_D4, "sepsis_D4_GSE54514.csv", row.names = FALSE)

############################# Day5 #############################################

sepsis_D5 <- transposed_merged %>%
  filter(Label %in% c("HC_D5", "NS_D5", "S_D5")) %>%
  mutate(Label = case_when(
    Label %in% c("HC_D5") ~ "Control",
    Label %in% c("NS_D5", "S_D5") ~ "Sepsis"
  ))

write.csv(sepsis_D5, "sepsis_D5_GSE54514.csv", row.names = FALSE)
