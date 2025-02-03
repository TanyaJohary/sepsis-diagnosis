################################################################################
############### Perform the Mann-Whitney U Test for One Gene ###################
################################################################################

# Load the dataset
data <- read.csv("sepsis_dataGSE243217.csv")
data <- dplyr::select(data, -Sample)

# Ensure the label column is a factor
data$Label <- as.factor(data$Label)

# This tests whether GeneX expression is significantly different between "Sepsis" and "Healthy" groups.
wilcox_test(TLR4 ~ Label, data = data, distribution = "exact")

################################################################################
############### Perform the Mann-Whitney U Test for all Genes ###################
################################################################################
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(coin, ggplot2, dplyr, tidyr)


# List of all genes (all columns except "Label")
gene_columns <- colnames(data)[colnames(data) != "Label"]

# Print gene columns for debugging
cat("Detected gene columns:\n")
print(gene_columns)

# Create an empty dataframe to store the Mann-Whitney test results
mann_whitney_results <- data.frame(
  Gene = character(),
  W = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each gene and perform the Wilcoxon test
for (gene in gene_columns) {
  cat("Testing gene:", gene, "\n")  # Debug print
  
  # Ensure the gene column is numeric (convert if needed)
  gene_data <- as.numeric(data[[gene]])
  
  # Perform the test comparing expression by Label
  test_result <- wilcox.test(gene_data ~ data$Label, exact = FALSE)
  
  # Append the results to the dataframe
  mann_whitney_results <- rbind(mann_whitney_results, 
                                data.frame(Gene = gene, 
                                           W = test_result$statistic, 
                                           p_value = test_result$p.value))
}

# Adjust p-values for multiple testing using the Benjamini-Hochberg (BH) method
mann_whitney_results$adjusted_p <- p.adjust(mann_whitney_results$p_value, method = "BH")

# Save the full results to a CSV file
write.csv(mann_whitney_results, "mann_whitney_results_adjustedGSE243217.csv", row.names = FALSE)

# Debug: Print counts of significant genes using raw and adjusted thresholds
cat("Number of genes with raw p_value < 0.05:", sum(mann_whitney_results$p_value < 0.05), "\n")
cat("Number of genes with adjusted p < 0.05:", sum(mann_whitney_results$adjusted_p < 0.05), "\n")

# Filter significant genes (adjusted p-value < 0.05) and order from most to least significant
significant_genes <- mann_whitney_results %>% 
  filter(adjusted_p < 0.05) %>% 
  arrange(adjusted_p)

# Save the ordered significant genes to a CSV file
write.csv(significant_genes, "significant_genesGSE243217.csv", row.names = FALSE)

# Print the ordered significant genes
cat("Significant genes (ordered from most to least significant):\n")
print(significant_genes)

# Extract non-significant genes (adjusted p-value >= 0.05) and order them
non_significant_genes <- mann_whitney_results %>% 
  filter(adjusted_p >= 0.05) %>% 
  arrange(adjusted_p)

# Save the non-significant genes to a CSV file
write.csv(non_significant_genes, "non_significant_genesGSE243217.csv", row.names = FALSE)

# Print the non-significant genes
cat("Non-significant genes (adjusted p-value >= 0.05):\n")
print(non_significant_genes)

################################################################################
################### Plot All Significant Genes Ordered by Significance #########
################################################################################

if (nrow(significant_genes) > 0) {
  
  # Order genes by adjusted p-value (lowest = most significant first)
  ordered_genes <- significant_genes %>% 
    arrange(adjusted_p) %>% 
    pull(Gene)
  
  # Subset the data to include only the "Label" column and the significant gene columns
  data_subset <- data %>% select(Label, all_of(ordered_genes))
  
  # Reshape the data to long format: one row per gene measurement
  data_long <- data_subset %>%
    pivot_longer(cols = -Label, names_to = "Gene", values_to = "Expression")
  
  # Convert 'Gene' into a factor with levels ordered by significance (most to least significant)
  data_long$Gene <- factor(data_long$Gene, levels = ordered_genes)
  
  # Create a facetted boxplot for the significant genes
  p <- ggplot(data_long, aes(x = Label, y = Expression, fill = Label)) +
    geom_boxplot() +
    facet_wrap(~ Gene, scales = "free_y") +
    labs(title = "Expression of Significant Genes by Sepsis Status\n(Ordered by Significance)",
         x = "Sepsis Status", 
         y = "Expression Level") +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.spacing = unit(1.5, "lines")  # Adjust spacing between facets here
    )
  
  # Display the plot
  print(p)
  
  # Save the plot as a PNG file
  ggsave("significant_genes_boxplot_orderedGSE243217.png", plot = p, width = 12, height = 8)
  
} else {
  cat("No significant genes found (adjusted p-value < 0.05).\n")
}
