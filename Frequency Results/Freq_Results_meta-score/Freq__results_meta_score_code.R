library(dplyr)
library(purrr)

#  each CSV has columns: Gene, adjusted_p, EffectSize, predicted_importance
file_paths <- list.files("~/sepsis/effect_size/", pattern="*.csv", full.names=TRUE)

# Function to get top 25% of meta_score from a single file
get_top25_metascore <- function(f) {
  df <- read.csv(f, header=TRUE, sep=" ", stringsAsFactors=FALSE)
  
  # rename columns if needed
  # if we have "GeneID" instead of "Gene", or "adjP" instead of "adjusted_p", adapt accordingly
  
  df <- df %>%
    mutate(
      neg_log_p   = -log10(adjusted_p),
      abs_effect  = abs(EffectSize),
      meta_score  = neg_log_p * abs_effect * predicted_importance
    )
  
  score_75 <- quantile(df$meta_score, 0.75, na.rm=TRUE)
  
  df_top <- df %>%
    filter(meta_score > score_75)
  
  return(df_top$Gene)
}

# Apply to all files
all_top_sets <- lapply(file_paths, get_top25_metascore)
names(all_top_sets) <- basename(file_paths)

# Intersection, union, majority vote
intersection_genes <- Reduce(intersect, all_top_sets)
union_genes <- Reduce(union, all_top_sets)

# Frequency table for majority vote
library(tidyr)
all_genes_df <- bind_rows(
  lapply(names(all_top_sets), function(ds) {
    data.frame(
      Gene   = all_top_sets[[ds]],
      Source = ds,
      stringsAsFactors = FALSE
    )
  })
)

gene_freq <- all_genes_df %>%
  group_by(Gene) %>%
  summarize(Count = n_distinct(Source)) %>%
  arrange(desc(Count))

#  majority fraction = 0.5
num_files <- length(file_paths)
threshold <- ceiling(num_files * 0.5)

majority_genes <- gene_freq %>%
  filter(Count >= threshold)

# Write results
write.table(intersection_genes, "intersection_meta_score.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(union_genes,        "union_meta_score.txt",        quote=FALSE, row.names=FALSE, col.names=FALSE)
write.csv(gene_freq,            "gene_freq_meta_score.csv",     row.names=FALSE)
write.csv(majority_genes,       "majority_genes_meta_score.csv", row.names=FALSE)
