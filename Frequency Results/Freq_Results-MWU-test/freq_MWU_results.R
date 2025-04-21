###############################################################################
#              TOP-25% BY ADJUSTED P-VALUE (LOWEST P)                  #
###############################################################################

library(dplyr)
library(tidyr)

# 0)  Settings -----------------------------------------------------------
data_folder <- "~/sepsis/MWU-tops/"
file_pattern <- "*.csv"

#  alias map for gene synonyms
alias_map <- c(
  "IL8"     = "CXCL8",
  "ELA2"    = "ELANE",
  "HLA.DRA" = "HLA-DRA"
)

# Possible column names for gene, adjusted p-value
gene_col_candidates <- c("Gene", "Feature", "GeneID")
pval_col_candidates <- c("adj_p", "padj", "adj_p_value", "p_val_adj", "adjusted_p")

# majority fraction for "majority vote" approach
majority_fraction <- 0.5

###############################################################################
# 1) List Files
###############################################################################
file_paths <- list.files(data_folder, pattern=file_pattern, full.names=TRUE)
cat("Found", length(file_paths), "files.\n")
if (length(file_paths) == 0) stop("No files found.")


###############################################################################
# 2) Read & Standardize Function
###############################################################################
#  - Renames gene column -> "GeneID"
#  - Renames p-value column -> "adj_p"
#  - Applies alias_map if needed
read_and_standardize_p <- function(file_path) {
  df <- read.csv(file_path, stringsAsFactors=FALSE)
  
  # Identify gene column
  actual_gene_col <- intersect(gene_col_candidates, colnames(df))
  if (length(actual_gene_col) == 0) {
    stop(paste("No gene column found in:", file_path))
  }
  if (length(actual_gene_col) > 1) {
    stop(paste("Multiple possible gene columns found in:", file_path, 
               "\nCandidates:", paste(actual_gene_col, collapse=", ")))
  }
  
  if (actual_gene_col != "GeneID") {
    colnames(df)[colnames(df) == actual_gene_col] <- "GeneID"
  }
  
  # Identify adjusted p-value column
  actual_pval_col <- intersect(pval_col_candidates, colnames(df))
  if (length(actual_pval_col) == 0) {
    stop(paste("No recognized p-value column found in:", file_path))
  }
  if (length(actual_pval_col) > 1) {
    stop(paste("Multiple possible p-value columns found in:", file_path,
               "\nCandidates:", paste(actual_pval_col, collapse=", ")))
  }
  
  # Rename p-value column -> "adj_p"
  if (actual_pval_col != "adj_p") {
    colnames(df)[colnames(df) == actual_pval_col] <- "adj_p"
  }
  
  # Apply small alias map if needed
  df$GeneID <- ifelse(df$GeneID %in% names(alias_map),
                      alias_map[df$GeneID],
                      df$GeneID)
  
  return(df)
}


###############################################################################
# METHOD A: Per-dataset top 25% (lowest p-values), then Combine
###############################################################################
# A.1) Extract top genes from each dataset
get_top25_by_p <- function(f) {
  df <- read_and_standardize_p(f)
  
  # We'll find the 25th percentile of p-values
  # (lowest 25% => p < p_25 cutoff)
  p_25 <- quantile(df$adj_p, 0.25, na.rm=TRUE)
  
  # Keep genes below that cutoff
  # i.e., genes with p < p_25
  top_df <- df %>%
    filter(adj_p < p_25)
  
  return(top_df$GeneID)
}

all_top_sets <- lapply(file_paths, get_top25_by_p)
names(all_top_sets) <- basename(file_paths)

# Intersection: genes in top-25% of *every* dataset
intersection_genes <- Reduce(intersect, all_top_sets)

# Union: genes in top-25% of *any* dataset
union_genes <- Reduce(union, all_top_sets)

# Frequency table: how many datasets each gene appears in
all_genes_df <- bind_rows(lapply(names(all_top_sets), function(ds) {
  data.frame(GeneID = all_top_sets[[ds]],
             Dataset= ds,
             stringsAsFactors = FALSE)
}))

gene_freq <- all_genes_df %>%
  group_by(GeneID) %>%
  summarize(Count = n_distinct(Dataset)) %>%
  arrange(desc(Count))

# Majority vote: genes in at least X of the Y datasets
num_files <- length(file_paths)
threshold <- ceiling(num_files * majority_fraction)

majority_genes <- gene_freq %>%
  filter(Count >= threshold)

write.table(intersection_genes, file = output_file_intersection_A,
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(union_genes,        file = output_file_union_A,
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.csv(gene_freq,           file = output_file_freq_A,
          row.names=FALSE)
write.csv(majority_genes,      file = output_file_majority_A,
          row.names=FALSE)



cat("\n--- METHOD A: Top 25% by Adjusted p-value (lowest p) ---\n")
cat("Intersection Genes:", length(intersection_genes), "\n")
cat("Union Genes:", length(union_genes), "\n")
cat("Majority Genes (≥", threshold, "):", nrow(majority_genes), "\n")



