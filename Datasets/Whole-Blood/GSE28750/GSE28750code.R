
# read the series matrix file
seriesmatrix <- read.delim("GSE28750_series_matrix.txt", header = TRUE, stringsAsFactors = FALSE, skip = 76)
write.csv(seriesmatrix, "seriesmatrixGSE28750.csv", row.names = TRUE)

geo_data <- getGEO("GSE28750", GSEMatrix = TRUE)

# Phenotypic Data
pheno_data <- pData(geo_data[[1]])
write.csv(pheno_data, "pheno_dataGSE28750.csv", row.names = TRUE)

# feature data
feature_data <- fData(geo_data[[1]])
write.csv(feature_data, "feature_dataGSE28750.csv", row.names = TRUE)

# extract the gene symbols that we need
sub_feature <- feature_data[, c("ID", "Gene Title", "Gene Symbol")]
write.csv(sub_feature, "sub_featureGSE28750.csv", row.names = TRUE)

# merge the expression data with gene symbols
expression_with_symbols <- merge(seriesmatrix, sub_feature, by.x = "ID_REF", by.y = "ID", all.x = TRUE )
write.csv(expression_with_symbols, "expression_with_symbolGSE28750.csv", row.names = TRUE)

# genes of interest names
genes_of_interest <- c(
  "CALCA", "CRP", "TREM1", "PLAUR", "LBP", "IL6", "HMGB1", "ELANE", "CXCL8",
  "IL10", "IL1B", "S100A8", "S100A9", "S100A12", "CD14", "FCGR1A", "ITGAM", "C3AR1",
  "C5AR1", "LCN2", "CXCL10", "IFNG", "IFNA1", "IFNA2", "IFNB1", "CCL19", "CCL25", "CX3CR1",
  "P2RX7", "PTX3", "TNFSF10", "MMP8", "MMP9", "HLADRA", "HLA-DRA", "TNF", "MYD88", "NLRP3",
  "TLR2", "TLR4", "NOTCH1", "BCL2", "PDCD1", "CCL2", "ARG1", "IL1R2", "CD177",
  "OLFM4", "MAPK14", "VCAM1", "ICAM1", "SOCS3", "GATA3", "CCR7", "CCR2", "HIF1A"
)

# find the present and missing genes
present_genes <- genes_of_interest[genes_of_interest %in% expression_with_symbols$Gene.Symbol]
write.csv(present_genes, "present_genesGSE28750.csv", row.names = TRUE)

missing_genes <- genes_of_interest[!genes_of_interest %in% expression_with_symbols$Gene.Symbol]
write.csv(missing_genes, "missing_genesGSE28750.csv", row.names = TRUE)


# search for aliases genes
alias_genes <- c(
  "IFNA1", "IFNA2", "IFNB1",
  "ELANE", "HLE", "PMN-E", "ELA2", "HNE", "NE",
  "MDNCF", "NAP-1", "GCP-1", "SCYB8", "LYNAP", "IL-8", "IL8"
)

# reload the expression with symbol and search for aliases
expression_with_symbols <- read.csv("expression_with_symbolGSE28750.csv", header = TRUE)
present_alias_genes <- alias_genes[alias_genes %in% expression_with_symbols$Gene.Symbol]
missing_alias_genes <- alias_genes[!alias_genes %in% expression_with_symbols$Gene.Symbol]

###########################################################################
####################### prepare data for Random forest#####################
###########################################################################

pheno_data <- read.csv("pheno_dataGSE28750.csv")

# Trim data 
data <- read.csv("expression_with_symbolGSE28750.csv")
data <- dplyr::select(data, -c("X", "ID_REF", "Gene.Title"))

# Remove rows where all values are NA
data <- data[rowSums(is.na(data)) < ncol(data), ]
write.csv(data, "expression_with_symbolGSE28750.csv", row.names = FALSE)

# filter out our genes of interest from expression data
 filtered_genes <- data %>%
   filter(Gene.Symbol %in% genes_of_interest)
 write.csv(filtered_genes, "filtered_genes_with_symbolGSE28750.csv", row.names = FALSE)
 
