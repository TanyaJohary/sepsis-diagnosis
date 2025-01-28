
# Retrieve the gse as a list
geo_data <- getGEO("GSE137342", GSEMatrix = TRUE)

expression_data10558 <- exprs(geo_data[[1]])

pheno_data10558 <- pData(geo_data[[1]])
sub_pheno10588 <- pheno_data10558[, c(2,8)]

feature_data10588 <- fData(geo_data[[1]])
sub_feature10588 <- feature_data10588[, c(1,13)]

# Merge expression data and sub feature to map the gene symbols
merged_data10558 <- merge(expression_data10558, sub_feature10588, by.x = "row.names", by.y = "ID", all.x = TRUE)

# Search for genes of interest 60
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "IL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A",
  "LY96", "CXCR4", "TRIL", "TLR5", "CXCL13"
)

present_genes <- genes_of_interest[genes_of_interest %in% merged_data10558$Symbol] 
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_data10558$Symbol] # PTX3 and MMP8

# Search for aliases
aliases_gene <- c(
  "TNFAIP5", "TSG-14", "TSG14", "Pentraxin 3",
  "CLG1", "PMNL-CL", "MMP-8", "HNC"
)

present_aliases <- aliases_gene[aliases_gene %in% merged_data10558$Symbol] # they are not present, even by aliases


# Save the present and missing genes
write.csv(present_genes, "present_genesGSE137342-10558.csv", row.names = TRUE)
write.csv(missing_genes, "missing_genesGSE137342-10558.csv", row.names = TRUE)

# Filter our genes of interest
filtered_genes10558 <- merged_data10558 %>%
  filter(Symbol %in% genes_of_interest)

# Remove the probbe id column
filtered_genes10558 <- filtered_genes10558[,-1]

# make unique the symbol column to preserve all probes foe each gene
filtered_genes10558$Symbol <- make.unique(as.character(filtered_genes10558$Symbol))

# Transpose the data frame to marge with phenotype data
transposed_10558 <- filtered_genes10558 %>%
  pivot_longer(
    cols = -Symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = Symbol,
    values_from = Value
  )

# Merge the transposed data to sub pheno data to identift the sepsis patients
trans_merged10558 <- merge(transposed_10558, sub_pheno10588, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

# Define the target column
trans_merged10558 <- trans_merged10558 %>%
  rename_with(~"Label", matches("source"))


################################################################################
##### Filter all sepsis patient include severe and septic shock- Day1 ##########
################################################################################

# Filter healthy ones and sepsis patients for day1 of sampling
sepsis_D1_10558 <- trans_merged10558 %>%
  filter((grepl("Day 1", Label) | grepl("Healthy", Label)) & 
           !grepl("Subject 65", Label) & 
           !grepl("Subject 66", Label))

# Simplify the target column values
sepsis_allD1_10558 <- sepsis_D1_10558 %>%
  mutate(Label = case_when(
    grepl("Healthy", Label) ~ "Healthy",
    grepl("Day 1", Label) ~ "Sepsis",
    TRUE ~ Label
  ))
write.csv(sepsis_allD1_10558, "sepsis_allD1_GSE137342_10558.csv", row.names = FALSE)

##################### Filter only sepsis patient exclude septic shock Day1
sepsis_onlyD1_10558 <- sepsis_D1_10558 %>%
  filter(!grepl("Septic Shock", Label))

# Simplify the target column values
sepsis_onlyD1_10558 <- sepsis_onlyD1_10558 %>%
  mutate(Label = case_when(
    grepl("Healthy", Label) ~ "Healthy",
    grepl("Day 1", Label) ~ "Sepsis",
    TRUE ~ Label
  ))
write.csv(sepsis_onlyD1_10558, "sepsis_onlyD1_GSE137342_10558.csv", row.names = FALSE)

##################### Filter only septic shock for Day1

septic_shockD1_10558 <- sepsis_D1_10558 %>%
  filter(!grepl("Sepsis", Label))

# Simplify the target column values
septic_shockD1_10558 <- septic_shockD1_10558 %>%
  mutate(Label = case_when(
    grepl("Healthy", Label) ~ "Healthy",
    grepl("Day 1", Label) ~ "Shock",
    TRUE ~ Label
  ))
write.csv(septic_shockD1_10558, "septic_shockD1_GSE137342_10558.csv", row.names = FALSE)

################################################################################
##################################### Day2 #####################################
################################################################################

# Filter healthy ones and sepsis patients for day2 of sampling
sepsis_D2_10558 <- trans_merged10558 %>%
  filter((grepl("Day 2", Label) | grepl("Healthy", Label)) & 
           !grepl("Subject 65", Label) & 
           !grepl("Subject 66", Label))

# Simplify the target column values
sepsis_allD2_10558 <- sepsis_D2_10558 %>%
  mutate(Label = case_when(
    grepl("Healthy", Label) ~ "Healthy",
    grepl("Day 2", Label) ~ "Sepsis",
    TRUE ~ Label
  ))
write.csv(sepsis_allD2_10558, "sepsis_allD2_GSE137342_10558.csv", row.names = FALSE)

##################### Filter only sepsis patient exclude septic shock Day1
sepsis_onlyD2_10558 <- sepsis_D2_10558 %>%
  filter(!grepl("Septic Shock", Label))

# Simplify the target column values
sepsis_onlyD2_10558 <- sepsis_onlyD2_10558 %>%
  mutate(Label = case_when(
    grepl("Healthy", Label) ~ "Healthy",
    grepl("Day 2", Label) ~ "Sepsis",
    TRUE ~ Label
  ))
write.csv(sepsis_onlyD2_10558, "sepsis_onlyD2_GSE137342_10558.csv", row.names = FALSE)

##################### Filter only septic shock for Day1

septic_shockD2_10558 <- sepsis_D2_10558 %>%
  filter(!grepl("Sepsis", Label))

# Simplify the target column values
septic_shockD2_10558 <- septic_shockD2_10558 %>%
  mutate(Label = case_when(
    grepl("Healthy", Label) ~ "Healthy",
    grepl("Day 2", Label) ~ "Shock",
    TRUE ~ Label
  ))
write.csv(septic_shockD2_10558, "septic_shockD2_GSE137342_10558.csv", row.names = FALSE)

