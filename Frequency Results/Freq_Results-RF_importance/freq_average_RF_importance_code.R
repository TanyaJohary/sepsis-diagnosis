###############################################################################
#                  UNIFIED SCRIPT FOR TOP-25% GENE SELECTION IN RF_RESULTS BETWEEN DATASETS                  #
###############################################################################
# This script demonstrates two methods:
#  - Method A: Top 25% within each dataset, then combine sets
#  - Method B: Merge all datasets, compute meta-importance, then top 25%
#
# It also:
#  - Handles varied column names for gene (Gene/Feature/GeneID) 
#  - Handles varied column names for importance (Mean_Importance, MeanImportance, avg_importance).
#  - Applies a small alias map to unify known synonyms (IL8 -> CXCL8, etc.).
###############################################################################

# Load necessary packages
library(dplyr)
library(purrr)
library(tidyr)

###############################################################################
# STEP 0: Configurations
###############################################################################
# 1) Path to folder containing our CSV files
data_folder <- "~/sepsis/average-feature-importance/" 

# 2) Pattern to match our CSV files
file_pattern <- "*.csv"

# 3) Possible column names for Gene ID
gene_col_candidates <- c("Gene", "Feature", "GeneID")

# 4) Possible column names for importance
importance_col_candidates <- c("Mean_Importance", "MeanImportance", "avg_importance")

# 5) Number of datasets to consider as "majority" (≥ fraction of total)
majority_fraction <- 0.5  # e.g., 0.5 means “in at least half of the datasets”

# 6) Filenames for output 
output_file_intersection_A <- "intersection_genes_methodA.txt"
output_file_union_A        <- "union_genes_methodA.txt"
output_file_freq_A         <- "gene_frequency_methodA.csv"
output_file_majority_A     <- "majority_genes_methodA.csv"

output_file_merged_B       <- "merged_importance_methodB.csv"
output_file_top25_B        <- "top25_genes_methodB.csv"

# 7) Small manual alias map for known synonyms
alias_map <- c(
  "IL8"     = "CXCL8",   # IL8 gets replaced by CXCL8
  "ELA2"    = "ELANE",   # ELA2 -> ELANE
  "HLA.DRA" = "HLA-DRA"  # HLA.DRA -> HLA-DRA
)

###############################################################################
# STEP 1: Discover and List All CSV Files
###############################################################################
file_paths <- list.files(
  path       = data_folder,
  pattern    = file_pattern,
  full.names = TRUE
)

cat("Found", length(file_paths), "files.\n")
if (length(file_paths) == 0) {
  stop("No files found. Check our 'data_folder' and 'file_pattern' settings.")
}

###############################################################################
# STEP 2: Function to Read & Standardize a Single CSV
###############################################################################
#  - Renames the gene column to "GeneID"
#  - Renames the importance column to "avg_importance"
#  - Replaces known aliases (IL8 -> CXCL8, etc.) via alias_map
read_and_standardize <- function(file_path,
                                 gene_candidates = gene_col_candidates,
                                 imp_candidates  = importance_col_candidates) {
  
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # 1) Identify the gene column
  actual_gene_col <- intersect(gene_candidates, colnames(df))
  if (length(actual_gene_col) == 0) {
    stop(paste("No gene column found in:", file_path))
  }
  if (length(actual_gene_col) > 1) {
    stop(paste("Multiple possible gene columns found in:", file_path, 
               "\nCandidates:", paste(actual_gene_col, collapse=", ")))
  }
  
  # Rename it to "GeneID"
  if (actual_gene_col != "GeneID") {
    colnames(df)[colnames(df) == actual_gene_col] <- "GeneID"
  }
  
  # 2) Identify the importance column
  actual_imp_col <- intersect(imp_candidates, colnames(df))
  if (length(actual_imp_col) == 0) {
    stop(paste("No importance column found in:", file_path))
  }
  if (length(actual_imp_col) > 1) {
    stop(paste("Multiple possible importance columns found in:", file_path,
               "\nCandidates:", paste(actual_imp_col, collapse=", ")))
  }
  
  # Rename it to "avg_importance"
  if (actual_imp_col != "avg_importance") {
    colnames(df)[colnames(df) == actual_imp_col] <- "avg_importance"
  }
  
  # 3) Apply our manual alias map to unify known synonyms
  #    Only change the GeneID if it matches a key in alias_map
  df$GeneID <- ifelse(df$GeneID %in% names(alias_map),
                      alias_map[df$GeneID],
                      df$GeneID)
  
  # Return the standardized data frame
  return(df)
}

###############################################################################
# METHOD A: Top 25% in Each Dataset Separately, Then Combine
###############################################################################

# A.1) Function to get top-25% genes for a single file
get_top25_genes_per_dataset <- function(f) {
  df <- read_and_standardize(file_path = f,
                             gene_candidates = gene_col_candidates,
                             imp_candidates  = importance_col_candidates)
  
  # Calculate 75th percentile of the "avg_importance"
  cutoff_75 <- quantile(df$avg_importance, 0.75, na.rm = TRUE)
  
  # Filter to genes above that cutoff
  top_df <- df %>%
    filter(avg_importance > cutoff_75)
  
  # Return just the gene IDs
  return(top_df$GeneID)
}

# A.2) Apply to all files
all_top_sets <- lapply(file_paths, get_top25_genes_per_dataset)
names(all_top_sets) <- basename(file_paths)

# A.3) Intersection & Union
intersection_genes_A <- Reduce(intersect, all_top_sets)
union_genes_A        <- Reduce(union, all_top_sets)

# A.4) Frequency Table (How many datasets does each gene appear in?)
all_genes_df_A <- bind_rows(
  lapply(names(all_top_sets), function(ds_name) {
    data.frame(
      GeneID  = all_top_sets[[ds_name]],
      Dataset = ds_name,
      stringsAsFactors = FALSE
    )
  })
)

gene_freq_A <- all_genes_df_A %>%
  group_by(GeneID) %>%
  summarize(Count = n_distinct(Dataset)) %>%
  arrange(desc(Count))

# A.5) Majority Vote (≥ fraction of datasets)
num_files <- length(file_paths)
threshold_majority <- ceiling(num_files * majority_fraction)

majority_genes_A <- gene_freq_A %>%
  filter(Count >= threshold_majority)

###############################################################################
# STEP 3: Write Out Method A Results 
###############################################################################
write.table(intersection_genes_A, file = output_file_intersection_A,
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(union_genes_A,        file = output_file_union_A,
            quote=FALSE, row.names=FALSE, col.names=FALSE)
write.csv(gene_freq_A,           file = output_file_freq_A,
          row.names=FALSE)
write.csv(majority_genes_A,      file = output_file_majority_A,
          row.names=FALSE)

cat("\n--- METHOD A RESULTS ---\n")
cat("Intersection Genes:", length(intersection_genes_A), "\n")
cat("Union Genes:", length(union_genes_A), "\n")
cat("Majority Genes (≥", threshold_majority, "):", nrow(majority_genes_A), "\n")

###############################################################################
# METHOD B: Merge All Datasets, Compute Meta-Importance, Then Top 25% (for testing and comparing)
###############################################################################
# B.1) Create a function to read & standardize, but keep only (GeneID, avg_importance)
#      Then rename avg_importance to a dataset-specific column for merging
read_for_merge <- function(f) {
  df <- read_and_standardize(file_path = f,
                             gene_candidates = gene_col_candidates,
                             imp_candidates  = importance_col_candidates)
  
  ds_name <- gsub("\\.csv$", "", basename(f))          # remove .csv
  import_col <- paste0("Importance_", ds_name)         # e.g., "Importance_dataset1"
  
  df_out <- df %>%
    dplyr::select(GeneID, avg_importance) %>%
    rename(!!import_col := avg_importance)
  
  return(df_out)
}

# B.2) Read & merge them by GeneID
df_list_B <- lapply(file_paths, read_for_merge)
merged_data_B <- Reduce(function(x, y) full_join(x, y, by = "GeneID"), df_list_B)

# B.3) Compute "meta-importance" (e.g., rowMeans across all Importance_* columns)
imp_cols_B <- setdiff(colnames(merged_data_B), "GeneID")
merged_data_B$MetaImportance <- rowMeans(merged_data_B[, imp_cols_B], na.rm = TRUE)

# B.4) Top 25% by MetaImportance
cutoff_75_B <- quantile(merged_data_B$MetaImportance, 0.75, na.rm=TRUE)
top25_B <- merged_data_B %>%
  filter(MetaImportance > cutoff_75_B)

# B.5) Write Out Method B Results 
write.csv(merged_data_B, output_file_merged_B, row.names=FALSE)
write.csv(top25_B,       output_file_top25_B, row.names=FALSE)

cat("\n--- METHOD B RESULTS ---\n")
cat("Total Genes in Merged Table:", nrow(merged_data_B), "\n")
cat("Top 25% Genes by MetaImportance:", nrow(top25_B), "\n")

###############################################################################
#                           
###############################################################################
