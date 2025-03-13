
# Retrieve the geo data
geo_data <- getGEO("GSE154918", GSEMatrix = TRUE)

# Retrieve the phenotype data
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c("geo_accession", "status:ch1")]
sub_pheno <- sub_pheno %>% 
  rename_with(~"Label", matches("status"))
sub_pheno$Label <- as.factor(sub_pheno$Label)
write.csv(sub_pheno, "sub_phenoGSE154918.csv", row.names = FALSE)


# Download the raw expression file from NCBI_GEO and read the file
raw_data<- read.delim("GSE154918_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(raw_data) <- raw_data$GeneID
raw_data <- raw_data[, -1]

all(colnames(raw_data) %in% rownames(sub_pheno))
all(rownames(sub_pheno) %in% colnames(raw_data))

# Create a Deseq2 objects
dds <- DESeqDataSetFromMatrix(countData = raw_data,
                              colData = sub_pheno,
                              design = ~ Label)

# Normalize the data
dds <- DESeq(dds)

dds <- dds[rowSums(counts(dds)) > 10, ] 

# variance stabilization
vsd <- vst(dds, blind = FALSE)
norm_data <- assay(vsd)
norm_data <- as.data.frame(norm_data)

zero_var_genes <- apply(norm_data, 1, function(x) var(x) == 0)
sum(zero_var_genes)  # Number of genes with zero variance

# Create PCA objects
pca_data <- prcomp(t(norm_data), scale. = TRUE)

#extract pca
pca_df <- as.data.frame(pca_data$x)
pca_df$Sample <- rownames(pca_df)
pca_df$batch <- pheno_data$`Sex:ch1`

ggplot(pca_df, aes(x = PC1, y = PC2, colour = batch ))+
  geom_point(size = 3)+
  labs(title = "PCA Plot: icu_acquired_infection_paired:ch1",
       x = "PC1 (First Principal Component)",
       y = "PC2 (Second Principal Component)") +
  theme_minimal() # there is no technical and oother batch effect

# Read the annotation data
annot_data <- read.delim("Human.GRCh38.p13.annot.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sub_annot <- annot_data[, c(1,2)]

# Merge normalized data to annotation data to map the genes symbols
norm_with_symbols <- merge(norm_data, sub_annot, by.x = "row.names", by.y = "GeneID", all.x = TRUE)
write.csv(norm_with_symbols, "normalized_with_symbolsGSE154918.csv", row.names = FALSE)

# search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)


present_genes <- genes_of_interest[genes_of_interest %in% norm_with_symbols$Symbol]
write.csv(present_genes, "present_genesGSE154918.csv", row.names = TRUE)

missing_genes <- genes_of_interest[!genes_of_interest %in% norm_with_symbols$Symbol]
write.csv(missing_genes, "missing_genesGSE154918.csv", row.names = TRUE)


# Remove the probe id column
norm_with_symbols <- norm_with_symbols[, -1]

# filter our genes of interest from expression data
filtered_genes <- norm_with_symbols %>%
  filter( Symbol %in% present_genes)

# transpose the data to make gene symbol as columns name
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

# dataset include expression changes for heathy control (Hlty), uncomplicated infection (Inf1_P),
# sepsis (Seps_P), septic shock (Shock_P), follow-up of sepsis (Seps_FU), follow-up of septic shock (Shock_FU) groups.
# should exclude the Inf1_P, Seps_FU and Shock_FU group from dataset

# Merge the transposed data with phenotype data to identify the sepsis patients
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

###### Filter the sepsis patients septic shock patients and healthy controls#######
filtered_sepsis <- merged_data %>%
  filter(Label %in% c("Hlty", "Seps_P", "Shock_P")) %>% 
  droplevels()

# clarify the label values
sepsis_data <- filtered_sepsis %>%
  mutate(Label = case_when(
    grepl("Hlty", Label) ~ "Healthy",
    grepl("Seps_P", Label) ~ "Sepsis",
    grepl("Shock_P", Label) ~ "Sepsis"
  ))

# Save the dataset contains the sepsis, septic shock patients and healthy controls
write.csv(sepsis_data, "sepsis_dataGSE154918.csv", row.names = FALSE)

############################################################# Selected genes

selected_genes <- c( "TREM1", "PLAUR", "IL6", "HMGB1", "ELANE", "CXCL8", "IL10", "IL1B",
                                       "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
                                       "C5AR1", "LCN2", "CXCL10", "IFNG", "CCL25", "CX3CR1", "P2RX7",
                                       "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88",
                                       "NLRP3", "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1",
                                       "IL1R2", "CD177", "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3",
                                       "GATA3", "CCR2", "HIF1A", "LY96", "CXCR4", "CXCL12" )


present_genes <- selected_genes[selected_genes %in% norm_with_symbols$Symbol]
missing_genes <- selected_genes[!selected_genes %in% norm_with_symbols$Symbol]


# filter our selected from expression data
filtered_genes <- norm_with_symbols %>%
  filter( Symbol %in% present_genes)

# transpose the data to make gene symbol as columns name
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

# dataset include expression changes for heathy control (Hlty), uncomplicated infection (Inf1_P),
# sepsis (Seps_P), septic shock (Shock_P), follow-up of sepsis (Seps_FU), follow-up of septic shock (Shock_FU) groups.
# should exclude the Inf1_P, Seps_FU and Shock_FU group from dataset

# Merge the transposed data with phenotype data to identify the sepsis patients
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

###### Filter the sepsis patients septic shock patients and healthy controls#######
filtered_sepsis <- merged_data %>%
  filter(Label %in% c("Hlty", "Seps_P", "Shock_P")) %>% 
  droplevels()

# clarify the label values
sepsis_selected <- filtered_sepsis %>%
  mutate(Label = case_when(
    grepl("Hlty", Label) ~ "Healthy",
    grepl("Seps_P", Label) ~ "Sepsis",
    grepl("Shock_P", Label) ~ "Sepsis"
  ))

# Save the dataset contains the sepsis, septic shock patients and healthy controls
write.csv(sepsis_selected, "sepsis_selectedGSE154918.csv", row.names = FALSE)

################################################################## all top genes

expression_with_symbols <- read.csv("normalized_with_symbolsGSE154918.csv")

top_genes <- c("CD177", "ITGAM", "ARG1", "C3AR1", "FCGR1A",
               "S100A8", "S100A9", "S100A12", "HMGB1",
               "GATA3", "SOCS3", "MAPK14", "IL10", "BCL2",
               "IL1R2", "MYD88", "HIF1A")

present_genes <- top_genes[top_genes %in% expression_with_symbols$Symbol]
missing_genes <- top_genes[!top_genes %in% expression_with_symbols$Symbol]

# filter our selected from expression data
filtered_genes <- expression_with_symbols %>%
  filter( Symbol %in% present_genes)
filtered_genes <- filtered_genes[, -1]

# transpose the data to make gene symbol as columns name
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

# dataset include expression changes for heathy control (Hlty), uncomplicated infection (Inf1_P),
# sepsis (Seps_P), septic shock (Shock_P), follow-up of sepsis (Seps_FU), follow-up of septic shock (Shock_FU) groups.
# should exclude the Inf1_P, Seps_FU and Shock_FU group from dataset

# Merge the transposed data with phenotype data to identify the sepsis patients
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

###### Filter the sepsis patients septic shock patients and healthy controls#######
filtered_sepsis <- merged_data %>%
  filter(Label %in% c("Hlty", "Seps_P", "Shock_P")) %>% 
  droplevels()

# clarify the label values
sepsis_selected <- filtered_sepsis %>%
  mutate(Label = case_when(
    grepl("Hlty", Label) ~ "Healthy",
    grepl("Seps_P", Label) ~ "Sepsis",
    grepl("Shock_P", Label) ~ "Sepsis"
  ))

# Save the dataset contains the sepsis, septic shock patients and healthy controls
write.csv(sepsis_selected, "sepsis_alltop_GSE154918.csv", row.names = FALSE)

############################################################# top cluster 1

cluster1 <- c("CD177", "ITGAM", "C3AR1", "FCGR1A")

present_genes <- cluster1[cluster1 %in% expression_with_symbols$Symbol]
missing_genes <- cluster1[!cluster1 %in% expression_with_symbols$Symbol]

# filter our selected from expression data
filtered_genes <- expression_with_symbols %>%
  filter( Symbol %in% present_genes)
filtered_genes <- filtered_genes[, -1]

# transpose the data to make gene symbol as columns name
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

# dataset include expression changes for heathy control (Hlty), uncomplicated infection (Inf1_P),
# sepsis (Seps_P), septic shock (Shock_P), follow-up of sepsis (Seps_FU), follow-up of septic shock (Shock_FU) groups.
# should exclude the Inf1_P, Seps_FU and Shock_FU group from dataset

# Merge the transposed data with phenotype data to identify the sepsis patients
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

###### Filter the sepsis patients septic shock patients and healthy controls#######
filtered_sepsis <- merged_data %>%
  filter(Label %in% c("Hlty", "Seps_P", "Shock_P")) %>% 
  droplevels()

# clarify the label values
sepsis_selected <- filtered_sepsis %>%
  mutate(Label = case_when(
    grepl("Hlty", Label) ~ "Healthy",
    grepl("Seps_P", Label) ~ "Sepsis",
    grepl("Shock_P", Label) ~ "Sepsis"
  ))

# Save the dataset contains the sepsis, septic shock patients and healthy controls
write.csv(sepsis_selected, "sepsis_cluster1_GSE154918.csv", row.names = FALSE)

################################################################# cluster2
cluster2 <- c("S100A8", "S100A9", "S100A12", "HMGB1")

present_genes <- cluster2[cluster2 %in% expression_with_symbols$Symbol]
missing_genes <- cluster2[!cluster2 %in% expression_with_symbols$Symbol]

# filter our selected from expression data
filtered_genes <- expression_with_symbols %>%
  filter( Symbol %in% present_genes)
filtered_genes <- filtered_genes[, -1]

# transpose the data to make gene symbol as columns name
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

# dataset include expression changes for heathy control (Hlty), uncomplicated infection (Inf1_P),
# sepsis (Seps_P), septic shock (Shock_P), follow-up of sepsis (Seps_FU), follow-up of septic shock (Shock_FU) groups.
# should exclude the Inf1_P, Seps_FU and Shock_FU group from dataset

# Merge the transposed data with phenotype data to identify the sepsis patients
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

###### Filter the sepsis patients septic shock patients and healthy controls#######
filtered_sepsis <- merged_data %>%
  filter(Label %in% c("Hlty", "Seps_P", "Shock_P")) %>% 
  droplevels()

# clarify the label values
sepsis_selected <- filtered_sepsis %>%
  mutate(Label = case_when(
    grepl("Hlty", Label) ~ "Healthy",
    grepl("Seps_P", Label) ~ "Sepsis",
    grepl("Shock_P", Label) ~ "Sepsis"
  ))

# Save the dataset contains the sepsis, septic shock patients and healthy controls
write.csv(sepsis_selected, "sepsis_cluster2_GSE154918.csv", row.names = FALSE)

################################################################ cluster 3

cluster3 <- c("GATA3", "SOCS3", "MAPK14", "IL10", "BCL2",
              "IL1R2", "MYD88", "HIF1A", "ARG1")

present_genes <- cluster3[cluster3 %in% expression_with_symbols$Symbol]
missing_genes <- cluster3[!cluster3 %in% expression_with_symbols$Symbol]

# filter our selected from expression data
filtered_genes <- expression_with_symbols %>%
  filter( Symbol %in% present_genes)
filtered_genes <- filtered_genes[, -1]

# transpose the data to make gene symbol as columns name
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

# dataset include expression changes for heathy control (Hlty), uncomplicated infection (Inf1_P),
# sepsis (Seps_P), septic shock (Shock_P), follow-up of sepsis (Seps_FU), follow-up of septic shock (Shock_FU) groups.
# should exclude the Inf1_P, Seps_FU and Shock_FU group from dataset

# Merge the transposed data with phenotype data to identify the sepsis patients
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

###### Filter the sepsis patients septic shock patients and healthy controls#######
filtered_sepsis <- merged_data %>%
  filter(Label %in% c("Hlty", "Seps_P", "Shock_P")) %>% 
  droplevels()

# clarify the label values
sepsis_selected <- filtered_sepsis %>%
  mutate(Label = case_when(
    grepl("Hlty", Label) ~ "Healthy",
    grepl("Seps_P", Label) ~ "Sepsis",
    grepl("Shock_P", Label) ~ "Sepsis"
  ))

# Save the dataset contains the sepsis, septic shock patients and healthy controls
write.csv(sepsis_selected, "sepsis_cluster3_GSE154918.csv", row.names = FALSE)

################################################################################
########################################################### last top genes

last_cluster <- c("CD177", "S100A8", "S100A9", "GATA3", "S100A12", "ARG1","ITGAM",
                  "SOCS3", "C3AR1", "IL1R2", "IL10", "FCGR1A", "MMP8", "MAPK14")

present_genes <- last_cluster[last_cluster %in% expression_with_symbols$Symbol]
missing_genes <- last_cluster[!last_cluster %in% expression_with_symbols$Symbol]

# filter our selected from expression data
filtered_genes <- expression_with_symbols %>%
  filter( Symbol %in% present_genes)
filtered_genes <- filtered_genes[, -1]

# transpose the data to make gene symbol as columns name
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

# dataset include expression changes for heathy control (Hlty), uncomplicated infection (Inf1_P),
# sepsis (Seps_P), septic shock (Shock_P), follow-up of sepsis (Seps_FU), follow-up of septic shock (Shock_FU) groups.
# should exclude the Inf1_P, Seps_FU and Shock_FU group from dataset

# Merge the transposed data with phenotype data to identify the sepsis patients
merged_data <- merge(transposed_data, sub_pheno, by.x = "Sample", by.y = "geo_accession", all.x = TRUE)

###### Filter the sepsis patients septic shock patients and healthy controls#######
filtered_sepsis <- merged_data %>%
  filter(Label %in% c("Hlty", "Seps_P", "Shock_P")) %>% 
  droplevels()

# clarify the label values
sepsis_selected <- filtered_sepsis %>%
  mutate(Label = case_when(
    grepl("Hlty", Label) ~ "Healthy",
    grepl("Seps_P", Label) ~ "Sepsis",
    grepl("Shock_P", Label) ~ "Sepsis"
  ))

# Save the dataset contains the sepsis, septic shock patients and healthy controls
write.csv(sepsis_selected, "sepsis_lasttop_GSE154918.csv", row.names = FALSE)
