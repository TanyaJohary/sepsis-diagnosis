# GSE69528

geo_data <- getGEO("GSE69528", GSEMatrix = TRUE)

pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c(2, 34)]

expression_data <- exprs(geo_data[[1]])

feature_data <- fData(geo_data[[1]])
sub_feature <- feature_data[, c(1, 13)]

pca_data <- prcomp(t(expression_data), scale. = TRUE)

# Create a PCA dataframe with metadata
pca_df <- data.frame(
  PC1 = pca_data$x[,1],
  PC2 = pca_data$x[,2],
  Batch = pheno_data$`pathogens:ch1`  
)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA of Expression Data Before Batch Correction") +
  labs(color = "pathogens") # no batch effect except for biological diffrences

expression_with_symbols <- merge(expression_data, sub_feature, by.x = "row.names", by.y = "ID", all.x = TRUE)

# Search for our clusters

cluster1 <- c("GATA3", "SOCS3", "MAPK14", "IL10", "BCL2",
              "IL1R2", "ARG1", "MMP9")

present_genes <- cluster1[cluster1 %in% expression_with_symbols$Symbol] 
missing_genes <- cluster1[!cluster1 %in% expression_with_symbols$Symbol]

filtered_genes <- expression_with_symbols %>% 
  filter(Symbol %in% present_genes)

filtered_genes$Symbol <- make.unique(filtered_genes$Symbol)
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

merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

sepsis_data <- merged_data %>% 
  filter(`study group:ch1` != "Uninfected type 2 diabetes mellitus")

sepsis_data <- sepsis_data %>% 
  rename_with(~"Label", matches("study group"))

cluster_1 <- sepsis_data %>% 
  mutate(Label = case_when(
    grepl("Uninfected healthy", Label) ~ "Healthy",
    grepl("Septicemic melioidosis", Label) ~ "Sepsis",
    grepl("Other sepsis", Label) ~ "Sepsis"
  ))

write.csv(cluster_1, "cluster1_GSE69528.csv", row.names = FALSE)

############################################################## cluster2
cluster2 <- c("FCGR1A", "CD177", "ITGAM", "C3AR1")

present_genes <- cluster2[cluster2 %in% expression_with_symbols$Symbol] 
missing_genes <- cluster2[!cluster2 %in% expression_with_symbols$Symbol]

filtered_genes <- expression_with_symbols %>% 
  filter(Symbol %in% present_genes)

filtered_genes$Symbol <- make.unique(filtered_genes$Symbol)
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

merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

sepsis_data <- merged_data %>% 
  filter(`study group:ch1` != "Uninfected type 2 diabetes mellitus")

sepsis_data <- sepsis_data %>% 
  rename_with(~"Label", matches("study group"))

cluster_2 <- sepsis_data %>% 
  mutate(Label = case_when(
    grepl("Uninfected healthy", Label) ~ "Healthy",
    grepl("Septicemic melioidosis", Label) ~ "Sepsis",
    grepl("Other sepsis", Label) ~ "Sepsis"
  ))

write.csv(cluster_2, "cluster2_GSE69528.csv", row.names = FALSE)


############################################################## cluster3
cluster3 <- c("S100A8", "S100A9", "S100A12")

present_genes <- cluster3[cluster3 %in% expression_with_symbols$Symbol] 
missing_genes <- cluster3[!cluster3 %in% expression_with_symbols$Symbol]

filtered_genes <- expression_with_symbols %>% 
  filter(Symbol %in% present_genes)

filtered_genes$Symbol <- make.unique(filtered_genes$Symbol)
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

merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

sepsis_data <- merged_data %>% 
  filter(`study group:ch1` != "Uninfected type 2 diabetes mellitus")

sepsis_data <- sepsis_data %>% 
  rename_with(~"Label", matches("study group"))

cluster_3 <- sepsis_data %>% 
  mutate(Label = case_when(
    grepl("Uninfected healthy", Label) ~ "Healthy",
    grepl("Septicemic melioidosis", Label) ~ "Sepsis",
    grepl("Other sepsis", Label) ~ "Sepsis"
  ))

write.csv(cluster_3, "cluster3_GSE69528.csv", row.names = FALSE)
