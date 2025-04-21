
################################################################################
#############################prepare data#######################################
################################################################################
library(GEOquery)
library(dplyr)
geo_data <- getGEO("GSE63311", GSEMatrix = TRUE)

expression_data <- exprs(geo_data[[1]]) # empty
expression_data <- read.csv("GSE63311_MappedCounts.csv")
colnames(expression_data)[1] <- "GeneID"
row.names(expression_data) <- expression_data$GeneID
expression_data <- expression_data[, -1]

# Extract and view phenotype data
pheno_data <- pData(geo_data[[1]])
sub_pheno <- pheno_data[, c("description", "group:ch1")]
rownames(sub_pheno) <- sub_pheno$description

sub_pheno <- sub_pheno %>%
  rename_with(~"Label", matches("group"))

sub_pheno <- sub_pheno %>%
  mutate(Label = case_when(
    Label == "Control T0" ~ "Control",
    Label == "No Sepsis T0" ~ "NO_sepsis",
    Label == "Sepsis T0" ~ "Sepsis",
    TRUE ~ Label
  ))


# check if the colnames in expression data is the same as rownames in sub_pheno data
all(colnames(expression_data) %in% rownames(sub_pheno))
setdiff(colnames(expression_data), rownames(sub_pheno)) # "S936742" this was not match, should remove it
setdiff(rownames(sub_pheno), colnames(expression_data)) 
expression_data <- expression_data[, colnames(expression_data) != "S936742"] # remove extra column
all(colnames(expression_data) == rownames(sub_pheno)) # not in a order
sub_pheno <- sub_pheno[colnames(expression_data),]

# normalizing the dataset
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("DESeq2")

library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = expression_data,
  colData = sub_pheno,
  design = ~ Label
)

# Filter lowly expressed genes
dds <- dds[rowSums(counts(dds)) > 10, ] 

# Estimate Size Factors for Normalization
dds <- DESeq(dds)

# Apply Variance Stabilizing Transformation (VST)
vsd <- vst(dds, blind = FALSE)

# Get normalized counts
normalized_data <- assay(vsd)

# check for batch effect
pca_data <- prcomp(t(normalized_data), scale. = TRUE)

# Extract PCA results
pca_df <- as.data.frame(pca_data$x)
pca_df$Sample <- rownames(pca_df)
pca_df$batch <- pheno_data$`culture positive:ch1`

# visualize the pca results
ggplot(pca_df, aes(x = PC1, y = PC2, colour = batch))+
  geom_point(size = 3)+
  labs(title = "PCA Plot: icu_acquired_infection_paired:ch1",
       x = "PC1 (First Principal Component)",
       y = "PC2 (Second Principal Component)") +
  theme_minimal() # no batch effect


# connect to ensembl database
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(normalized_data),
  mart = ensembl
)

# Merge gene symbols to vst normalized data
merged_vst_symbol <- merge(normalized_data, gene_mapping, by.x = "row.names", by.y = "ensembl_gene_id", all.x = TRUE)

# search for genes of interest
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

present_genes <- genes_of_interest[genes_of_interest %in% merged_vst_symbol$hgnc_symbol]
missing_genes <- genes_of_interest[!genes_of_interest %in% merged_vst_symbol$hgnc_symbol]

genes_ids <- c("ENSG00000110680", "ENSG00000129988", "ENSG00000131142")
present_ids <- genes_ids[genes_ids %in% merged_vst_symbol$Row.names] ## no result, genes CALCA,LBP and CCL25 are miissing in dataset

write.csv(present_genes, "present_genesGSE63311.csv", row.names = FALSE)
write.csv(missing_genes, "missing_genesGSE63311.csv", row.names = FALSE)

# Filter our genes of interest from dataset
filtered_vst_genes <- merged_vst_symbol %>%
  filter(hgnc_symbol %in% genes_of_interest)

# Remove the gene ensembl id
filtered_vst_genes <- filtered_vst_genes[,-1]

# Transpose the daatset
transposed_vst_genes <- filtered_vst_genes %>%
  pivot_longer(
    cols = -hgnc_symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = hgnc_symbol,
    values_from = Value
  )

# Merge the transposed data to pheno data to identify sepsis patients
merged_transposed_vst <- merge(transposed_vst_genes, sub_pheno, by.x = "Sample", by.y = "description", all.x = TRUE)

# Filter just sepsis and control ones
filtered_sepsis_vst <- merged_transposed_vst %>%
  filter(Label %in% c("Sepsis", "Control"))
write.csv(filtered_sepsis_vst, "filtered_sepsis_vstGSE63311.csv", row.names = FALSE)

############################################################# All top genes

top_genes <- c(
  "S100A8", "S100A9", "S100A12",
  "IL10", "MAPK14", "SOCS3", "GATA3", "ARG1", "IL1R2", "MMP9", "BCL2",
  "FCGR1A", "ITGAM", "C3AR1", "CD177")

# Filter our genes of interest from dataset
filtered_genes <- merged_vst_symbol %>%
  filter(hgnc_symbol %in% top_genes)

# Remove the gene ensembl id
filtered_genes <- filtered_genes[,-1]

# Transpose the daatset
transposed_genes <- filtered_genes %>%
  pivot_longer(
    cols = -hgnc_symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = hgnc_symbol,
    values_from = Value
  )

# Merge the transposed data to pheno data to identify sepsis patients
merged_transposed <- merge(transposed_genes, sub_pheno, by.x = "Sample", by.y = "description", all.x = TRUE)

# Filter just sepsis and control ones
filtered_sepsis <- merged_transposed %>%
  filter(Label %in% c("Sepsis", "Control"))
write.csv(filtered_sepsis, "topgenes_GSE63311.csv", row.names = FALSE)

############################################################# top cluster 1

cluster1 <- c("IL10", "MAPK14", "SOCS3", "GATA3", "ARG1", "IL1R2", "MMP9", "BCL2")

# Filter our genes of interest from dataset
filtered_genes <- merged_vst_symbol %>%
  filter(hgnc_symbol %in% cluster1)

# Remove the gene ensembl id
filtered_genes <- filtered_genes[,-1]

# Transpose the daatset
transposed_genes <- filtered_genes %>%
  pivot_longer(
    cols = -hgnc_symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = hgnc_symbol,
    values_from = Value
  )

# Merge the transposed data to pheno data to identify sepsis patients
merged_transposed <- merge(transposed_genes, sub_pheno, by.x = "Sample", by.y = "description", all.x = TRUE)

# Filter just sepsis and control ones
filtered_sepsis <- merged_transposed %>%
  filter(Label %in% c("Sepsis", "Control"))
write.csv(filtered_sepsis, "cluster1_GSE63311.csv", row.names = FALSE)

################################################################################
################################################################# cluster2
cluster2 <- c("ITGAM", "CD177", "C3AR1", "FCGR1A")

# Filter our genes of interest from dataset
filtered_genes <- merged_vst_symbol %>%
  filter(hgnc_symbol %in% cluster2)

# Remove the gene ensembl id
filtered_genes <- filtered_genes[,-1]

# Transpose the daatset
transposed_genes <- filtered_genes %>%
  pivot_longer(
    cols = -hgnc_symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = hgnc_symbol,
    values_from = Value
  )

# Merge the transposed data to pheno data to identify sepsis patients
merged_transposed <- merge(transposed_genes, sub_pheno, by.x = "Sample", by.y = "description", all.x = TRUE)

# Filter just sepsis and control ones
filtered_sepsis <- merged_transposed %>%
  filter(Label %in% c("Sepsis", "Control"))
write.csv(filtered_sepsis, "cluster2_GSE63311.csv", row.names = FALSE)

################################################################ cluster 3

cluster3 <- c("S100A8", "S100A9", "S100A12")

# Filter our genes of interest from dataset
filtered_genes <- merged_vst_symbol %>%
  filter(hgnc_symbol %in% cluster3)

# Remove the gene ensembl id
filtered_genes <- filtered_genes[,-1]

# Transpose the daatset
transposed_genes <- filtered_genes %>%
  pivot_longer(
    cols = -hgnc_symbol,
    names_to = "Sample",
    values_to = "Value"
  ) %>%
  pivot_wider(
    names_from = hgnc_symbol,
    values_from = Value
  )

# Merge the transposed data to pheno data to identify sepsis patients
merged_transposed <- merge(transposed_genes, sub_pheno, by.x = "Sample", by.y = "description", all.x = TRUE)

# Filter just sepsis and control ones
filtered_sepsis <- merged_transposed %>%
  filter(Label %in% c("Sepsis", "Control"))
write.csv(filtered_sepsis, "cluster3_GSE63311.csv", row.names = FALSE)
