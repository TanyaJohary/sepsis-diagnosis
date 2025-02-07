################################################################################
############### Perform the Mann-Whitney U Test for All Genes ##################
################################################################################

# Load necessary libraries
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(coin, ggplot2, dplyr, tidyr, readr)

# Load the dataset
data <- read.csv("T2_sepsis_labeledGSE69063.csv") %>% select(-Sample)

# Ensure 'Label' is a factor
data$Label <- as.factor(data$Label)

# List of all genes (all columns except "Label")
gene_columns <- setdiff(names(data), "Label")

cat("Detected gene columns:\n")
print(gene_columns)

# Create a list to store results (avoiding slow rbind operations)
mann_whitney_list <- list()

# Perform the Mann-Whitney U test for each gene
for (gene in gene_columns) {
  cat("Testing gene:", gene, "\n")
  
  # Ensure the gene column is numeric
  gene_data <- as.numeric(data[[gene]])
  
  # Perform the Wilcoxon rank-sum test
  test_result <- wilcox.test(gene_data ~ data$Label, exact = FALSE)
  
  # Store results in a named list
  mann_whitney_list[[gene]] <- list(
    Gene = gene,
    W = unname(test_result$statistic),  # Extracts W-value properly
    p_value = test_result$p.value
  )
}

# Convert list to dataframe
mann_whitney_results <- bind_rows(mann_whitney_list)

# Debugging: Print column names to check correctness
print(colnames(mann_whitney_results))

# Ensure correct data types and order by adjusted p-value
mann_whitney_results <- mann_whitney_results %>%
  mutate(W = as.numeric(W), p_value = as.numeric(p_value)) %>%
  mutate(adjusted_p = p.adjust(p_value, method = "BH")) %>%
  arrange(adjusted_p)

# Save results
write.csv(mann_whitney_results, "mann_whitney_results_adjustedGSE60424.csv", row.names = FALSE)


# Debug: Print summary of significant genes
cat("Number of genes with raw p-value < 0.05:", sum(mann_whitney_results$p_value < 0.05), "\n")
cat("Number of genes with adjusted p-value < 0.05:", sum(mann_whitney_results$adjusted_p < 0.05), "\n")

# Separate significant and non-significant genes
significant_genes <- filter(mann_whitney_results, adjusted_p < 0.05)
non_significant_genes <- filter(mann_whitney_results, adjusted_p >= 0.05)

# Save significant and non-significant genes
write_csv(significant_genes, "significant_genesGSE69063_T2.csv")
write_csv(non_significant_genes, "non_significant_genesGSE69063_T2.csv")

# Print significant genes summary
cat("Significant genes (ordered from most to least significant):\n")
print(significant_genes)

################################################################################
################### Plot All Significant Genes Ordered by Significance #########
################################################################################
# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)



if (nrow(significant_genes) > 0) {
  
  # Select only significant genes
  ordered_genes <- significant_genes$Gene
  
  # Reshape the data into long format for plotting
  data_long <- data %>%
    select(Label, all_of(ordered_genes)) %>%
    pivot_longer(cols = -Label, names_to = "Gene", values_to = "Expression")
  
  # Convert Gene to factor for ordered plotting
  data_long$Gene <- factor(data_long$Gene, levels = ordered_genes)
  
  # Improved boxplot with better color, grid, and styling
  p <- ggplot(data_long, aes(x = Label, y = Expression, fill = Label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Reduce opacity and remove outliers
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) + # Add jittered points for visibility
    facet_wrap(~ Gene, scales = "free_y") +
    scale_fill_manual(values = c("sepsis" = "orange", "healthy control" = "grey")) +  # Custom colors
    labs(
      title = "Expression of Significant Genes by Sepsis Status\n(Ordered by Significance)",
      subtitle = "Mann-Whitney U test (BH-adjusted p-value < 0.05)",
      x = "Sepsis Status", 
      y = "Expression Level"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(color = "gray90", linetype = "dashed"), # Light grid
      strip.text = element_text(size = 8, face = "bold", color = "darkgrey"), # Better facet labels
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Display and save the plot
  print(p)
  ggsave("significant_genes_boxplot_orderedGSE60424.png", plot = p, width = 14, height = 10, dpi = 300)
  
} else {
  cat("No significant genes found (adjusted p-value < 0.05).\n")
}
