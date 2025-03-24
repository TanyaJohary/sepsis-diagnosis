
####################prepare data for Random-forest##########################
############################################################################

getGEOSuppFiles("GSE54514")  # Retrieve the raw data for normaliztaion

untar("GSE54514/GSE54514_RAW.tar", exdir= 'data/')

raw_data <- ReadAffy(celfile.path = "/home/rstudio/sepsis/data/") # Raw data is not available


################################################################################
#########################prepare dataset for RF#################################
################################################################################

# Retrieve the geo data list
geo_data <- getGEO("GSE54514", GSEMatrix = TRUE)

# Retrieve the expression data
expression_data <- exprs(geo_data[[1]]) # this data is already normalized based on related publication

# Retrieve the phenotype data
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c(2, 12 , 43)]
sub_pheno <- sub_pheno %>% 
  rename_with(~"Label", matches("disease"))
write.csv(sub_pheno, "sub_phenoGSE54514.csv", row.names = FALSE)


pca_data <- prcomp(t(expression_data), scale. = TRUE)
pca_df <- as.data.frame(pca_data$x)
pca_df$Sample <- rownames(pca_df)
pca_df$batch <- pheno_data$`severity (apacheii):ch1`

ggplot(pca_df, aes(x = PC1, y = PC2, colour = batch))+
  geom_point(size = 3)+
  labs(title = "PCA Plot: icu_acquired_infection_paired:ch1",
       x = "PC1 (First Principal Component)",
       y = "PC2 (Second Principal Component)") +
  theme_minimal() # no batch effect


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
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_data$Symbol] # "CRP""LBP" "IFNA1""IFNA2" "IFNB1" "PTX3"

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


# Merged phenotype data to expression data to identify sepsis and healthy cases
transposed_merged <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)


# Make a diagnostic dataset- sepsis patients all days, D1-D5
sepsis_data <- transposed_merged %>%
  mutate(Label = case_when(
    Label == "sepsis nonsurvivor" ~ "Sepsis",
    Label == "sepsis survivor" ~ "Sepsis",
    TRUE ~ Label
  ))

######################################## Day1 ##################################

# Filter sepsis patients for D1 and healthy one 
sepsis_D1 <- sepsis_data %>%
  filter(Label == "healthy" |
           (Label == "Sepsis" & characteristics_ch1.2 %in% c("group_day: S_D1", "group_day: NS_D1")))
sepsis_D1 <- sepsis_D1[, -76]

write.csv(sepsis_D1, "sepsis_D1_GSE54514.csv", row.names = FALSE)

############################## Day2 ############################################

# Filter sepsis patients for D2 and healthy one 
sepsis_D2 <- sepsis_data %>%
  filter(Label == "healthy" |
           (Label == "Sepsis" & characteristics_ch1.2 %in% c("group_day: S_D2", "group_day: NS_D2")))
sepsis_D2 <- sepsis_D2[, -76]

write.csv(sepsis_D2, "sepsis_D2_GSE54514.csv", row.names = FALSE)

################################# Day3 #########################################
# Filter sepsis patients for D3 and healthy one 
sepsis_D3 <- sepsis_data %>%
  filter(Label == "healthy" |
           (Label == "Sepsis" & characteristics_ch1.2 %in% c("group_day: S_D3", "group_day: NS_D3")))
sepsis_D3 <- sepsis_D3[, -76]

write.csv(sepsis_D3, "sepsis_D3_GSE54514.csv", row.names = FALSE)

################################ Day4 ##########################################
# Filter sepsis patients for D4 and healthy one 
sepsis_D4 <- sepsis_data %>%
  filter(Label == "healthy" |
           (Label == "Sepsis" & characteristics_ch1.2 %in% c("group_day: S_D4", "group_day: NS_D4")))
sepsis_D4 <- sepsis_D4[, -76]

write.csv(sepsis_D4, "sepsis_D4_GSE54514.csv", row.names = FALSE)

############################# Day5 #############################################

sepsis_D5 <- sepsis_data %>%
  filter(Label == "healthy" |
           (Label == "Sepsis" & characteristics_ch1.2 %in% c("group_day: S_D5", "group_day: NS_D5")))
sepsis_D5 <- sepsis_D5[, -76]

write.csv(sepsis_D5, "sepsis_D5_GSE54514.csv", row.names = FALSE)

################################################################################
############################################################### all top genes D1

top_genes <- c("CD177", "ITGAM", "ARG1", "C3AR1", "FCGR1A",
               "S100A8", "S100A9", "S100A12", "HMGB1",
               "GATA3", "SOCS3", "MAPK14", "IL10", "BCL2",
               "IL1R2", "MYD88", "HIF1A")

present_genes <- top_genes[top_genes %in% merged_data$Symbol]
missing_genes <- top_genes[!top_genes %in% merged_data$Symbol] 

# Filter genes of interest from dataset
filtered_genes <- merged_data %>%
  filter(Symbol %in% present_genes)

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


# Merged phenotype data to expression data to identify sepsis and healthy cases
transposed_merged <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)


# Make a diagnostic dataset- sepsis patients all days, D1-D5
sepsis_data <- transposed_merged %>%
  mutate(Label = case_when(
    Label == "sepsis nonsurvivor" ~ "Sepsis",
    Label == "sepsis survivor" ~ "Sepsis",
    TRUE ~ Label
  ))

# Filter sepsis patients for D1 and healthy one 
sepsis_D1 <- sepsis_data %>%
  filter(Label == "healthy" |
           (Label == "Sepsis" & characteristics_ch1.2 %in% c("group_day: S_D1", "group_day: NS_D1")))
sepsis_D1 <- sepsis_D1[, -30]

write.csv(sepsis_D1, "sepsis_alltopD1_GSE54514.csv", row.names = FALSE)

################################################################################
############################################################### top cluster1  D1

cluster1 <- c("CD177", "ITGAM", "C3AR1", "FCGR1A")

present_genes <- cluster1[cluster1 %in% merged_data$Symbol]
missing_genes <- cluster1[!cluster1 %in% merged_data$Symbol] 

# Filter genes of interest from dataset
filtered_genes <- merged_data %>%
  filter(Symbol %in% present_genes)

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


# Merged phenotype data to expression data to identify sepsis and healthy cases
transposed_merged <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)


# Make a diagnostic dataset- sepsis patients all days, D1-D5
sepsis_data <- transposed_merged %>%
  mutate(Label = case_when(
    Label == "sepsis nonsurvivor" ~ "Sepsis",
    Label == "sepsis survivor" ~ "Sepsis",
    TRUE ~ Label
  ))

# Filter sepsis patients for D1 and healthy one 
sepsis_D1 <- sepsis_data %>%
  filter(Label == "healthy" |
           (Label == "Sepsis" & characteristics_ch1.2 %in% c("group_day: S_D1", "group_day: NS_D1")))
sepsis_D1 <- sepsis_D1[, -6]

write.csv(sepsis_D1, "sepsis_cluster1D1_GSE54514.csv", row.names = FALSE)

################################################################################
############################################################### top cluster2  D1

cluster2 <- c("S100A8", "S100A9", "S100A12", "HMGB1")

present_genes <- cluster2[cluster2 %in% merged_data$Symbol]
missing_genes <- cluster2[!cluster2 %in% merged_data$Symbol] 

# Filter genes of interest from dataset
filtered_genes <- merged_data %>%
  filter(Symbol %in% present_genes)

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


# Merged phenotype data to expression data to identify sepsis and healthy cases
transposed_merged <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)


# Make a diagnostic dataset- sepsis patients all days, D1-D5
sepsis_data <- transposed_merged %>%
  mutate(Label = case_when(
    Label == "sepsis nonsurvivor" ~ "Sepsis",
    Label == "sepsis survivor" ~ "Sepsis",
    TRUE ~ Label
  ))

# Filter sepsis patients for D1 and healthy one 
sepsis_D1 <- sepsis_data %>%
  filter(Label == "healthy" |
           (Label == "Sepsis" & characteristics_ch1.2 %in% c("group_day: S_D1", "group_day: NS_D1")))
sepsis_D1 <- sepsis_D1[, -7]

write.csv(sepsis_D1, "sepsis_cluster2D1_GSE54514.csv", row.names = FALSE)

################################################################################
############################################################### top cluster3  D1

cluster3 <- c("GATA3", "SOCS3", "MAPK14", "IL10", "BCL2",
              "IL1R2", "MYD88", "HIF1A", "ARG1")

present_genes <- cluster3[cluster3 %in% merged_data$Symbol]
missing_genes <- cluster3[!cluster3 %in% merged_data$Symbol] 

# Filter genes of interest from dataset
filtered_genes <- merged_data %>%
  filter(Symbol %in% present_genes)

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


# Merged phenotype data to expression data to identify sepsis and healthy cases
transposed_merged <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)


# Make a diagnostic dataset- sepsis patients all days, D1-D5
sepsis_data <- transposed_merged %>%
  mutate(Label = case_when(
    Label == "sepsis nonsurvivor" ~ "Sepsis",
    Label == "sepsis survivor" ~ "Sepsis",
    TRUE ~ Label
  ))

# Filter sepsis patients for D1 and healthy one 
sepsis_D1 <- sepsis_data %>%
  filter(Label == "healthy" |
           (Label == "Sepsis" & characteristics_ch1.2 %in% c("group_day: S_D1", "group_day: NS_D1")))
sepsis_D1 <- sepsis_D1[, -21]

write.csv(sepsis_D1, "sepsis_cluster3D1_GSE54514.csv", row.names = FALSE)
