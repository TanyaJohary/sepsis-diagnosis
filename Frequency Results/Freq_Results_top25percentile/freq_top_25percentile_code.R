# Load needed packages
library(dplyr)

# 1. List our CSV files
#    Adjust the folder path and pattern as needed
file_paths <- list.files(
  path = "~/sepsis/top25percentile/",  
  pattern = "*.csv", 
  full.names = TRUE
)

# Sanity check: How many files do we have?
length(file_paths)

# 2. Read gene sets from each file
#    Assume each file has a column named 'GeneID' or similar
all_gene_sets <- list()

for (f in file_paths) {
  df <- read.csv(f)  # Adjust sep = " " if needed
  # Extract the gene column as characters
  genes <- as.character(df$Gene)  
  
  # Use file name (without path) as the list element name
  dataset_name <- basename(f)
  all_gene_sets[[dataset_name]] <- genes
}

# 3. Compute intersection: Genes that appear in *all* top-25% lists
intersection_genes <- Reduce(intersect, all_gene_sets)

# 4. Compute union: Genes that appear in *any* top-25% list
union_genes <- Reduce(union, all_gene_sets)

# 5. Calculate frequency for each gene across all files
#    (How many files does each gene appear in?)
#    First convert each gene set to a small data frame with columns: GeneID, Dataset
all_genes_df <- bind_rows(lapply(names(all_gene_sets), function(ds) {
  data.frame(GeneID = all_gene_sets[[ds]], 
             Dataset = ds, 
             stringsAsFactors = FALSE)
}))

# Count how many distinct files each gene appears in
gene_freq <- all_genes_df %>%
  group_by(GeneID) %>%
  summarize(Count = n_distinct(Dataset)) %>%
  arrange(desc(Count))

# 6. Majority Vote: Genes that appear in >= 50% of files

num_files <- length(file_paths)
threshold <- ceiling(num_files * 0.5)  # for 50%

majority_genes <- gene_freq %>%
  filter(Count >= threshold) %>%
  arrange(desc(Count))

# 7. Write our results 
write.table(intersection_genes, 
            file = "intersection_genes.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

write.table(union_genes, 
            file = "union_genes.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

write.csv(gene_freq, 
          "gene_frequency_table.csv", 
          row.names = FALSE)

write.csv(majority_genes, 
          "majority_genes.csv", 
          row.names = FALSE)

# 8. Quick Summaries
cat("Number of files:", num_files, "\n")
cat("Intersection Genes:", length(intersection_genes), "\n")
cat("Union Genes:", length(union_genes), "\n")
cat("Majority Genes (Count >= ", threshold, "):", nrow(majority_genes), "\n")

# The objects 'intersection_genes', 'union_genes', 'gene_freq', and 'majority_genes'
# are now in our environment for further analysis or inspection.
