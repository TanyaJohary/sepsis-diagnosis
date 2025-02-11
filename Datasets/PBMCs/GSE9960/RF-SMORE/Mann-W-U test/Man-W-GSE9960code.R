###############################################################################
# Mann-Whitney U Test for All Genes with Effect Size, Direction, and Enhanced Plot
###############################################################################

# Load necessary libraries
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  coin, effsize, ggplot2, dplyr, tidyr, readr
)

# Load the dataset
data <- read.csv("sepsis_dataGSE9960.csv")
data <- dplyr::select(data, -Sample)  # remove 'Sample' column if present

# Ensure 'Label' is a factor
data$Label <- as.factor(data$Label)

# List of all genes (all columns except "Label")
gene_columns <- setdiff(names(data), "Label")

cat("Detected gene columns:\n")
print(gene_columns)

# Create a list to store results (efficiently)
mann_whitney_list <- list()

###############################################################################
# Perform the Mann-Whitney U test for each gene, compute effect size + direction
###############################################################################
for (gene in gene_columns) {
  # Convert the gene column to numeric just to be sure
  gene_data <- as.numeric(data[[gene]])
  
  # Split into two groups (assuming 'Sepsis' vs. 'Control')
  sepsis_vals  <- gene_data[data$Label == "Sepsis"]
  control_vals <- gene_data[data$Label == "Control"]
  
  # Wilcoxon rank-sum test (Mann-Whitney U)
  test_result <- wilcox.test(sepsis_vals, control_vals, exact = FALSE)
  
  # Cliff’s delta effect size
  # (If either group is empty, cliff.delta() will fail, so double-check lengths)
  effect_size <- NA
  if (length(sepsis_vals) > 0 && length(control_vals) > 0) {
    cliff_res <- cliff.delta(sepsis_vals, control_vals)
    effect_size <- cliff_res$estimate  # numeric effect size
  }
  
  # Calculate direction of change
  mean_sepsis  <- mean(sepsis_vals, na.rm = TRUE)
  mean_control <- mean(control_vals, na.rm = TRUE)
  direction    <- ifelse(mean_sepsis > mean_control, "Up_in_Sepsis", "Down_in_Sepsis")
  
  # Store in list
  mann_whitney_list[[gene]] <- list(
    Gene        = gene,
    Mean_Sepsis = mean_sepsis,
    Mean_Control = mean_control,
    Direction   = direction,
    W           = unname(test_result$statistic),  # W-value
    p_value     = test_result$p.value,
    effect_size = effect_size
  )
}

###############################################################################
# Convert list to a data frame & adjust p-values
###############################################################################
mann_whitney_results <- dplyr::bind_rows(mann_whitney_list) %>%
  mutate(
    # Convert to numeric just in case
    W          = as.numeric(W),
    p_value    = as.numeric(p_value),
    effect_size = as.numeric(effect_size),
    adjusted_p = p.adjust(p_value, method = "BH")  # BH-corrected p-values
  ) %>%
  arrange(adjusted_p)  # order by significance

# Save full results
write.csv(mann_whitney_results, "mann_whitney_results_adjusted_GSE9960.csv", row.names = FALSE)

# Debug: Summaries of significance
cat("Number of genes with raw p-value < 0.05:",
    sum(mann_whitney_results$p_value < 0.05), "\n")
cat("Number of genes with adjusted p-value < 0.05:",
    sum(mann_whitney_results$adjusted_p < 0.05), "\n")

# Separate significant and non-significant
significant_genes <- dplyr::filter(mann_whitney_results, adjusted_p < 0.05)
non_significant_genes <- dplyr::filter(mann_whitney_results, adjusted_p >= 0.05)

write_csv(significant_genes,"significant_genes_GSE9960.csv")
write_csv(non_significant_genes, "non_significant_genes_GSE9960.csv")

cat("Significant genes (ordered by adjusted p-value):\n")
print(significant_genes)

###############################################################################
# Plot Top 25 Significant Genes, Annotating Adjusted p-values
###############################################################################
if (nrow(significant_genes) > 0) {
  
  # Select only the top 25 most significant genes
  top_25_genes <- significant_genes %>%
    slice_min(order_by = adjusted_p, n = 25, with_ties = FALSE)
  
  # Create a new column that combines gene name + scientific p-value
  top_25_genes <- top_25_genes %>%
    mutate(Gene_label = paste0(
      Gene,
      "\n(adj p=", formatC(adjusted_p, format = "e", digits = 2), ")"
    ))
  
  # We'll keep track of these new labels in the same order
  ordered_labels <- top_25_genes$Gene_label
  
  # Reshape the data into long format for plotting, but only for these top 25 genes
  data_long <- data %>%
    dplyr::select(Label, dplyr::all_of(top_25_genes$Gene)) %>%
    tidyr::pivot_longer(cols = -Label, names_to = "Gene", values_to = "Expression")
  
  # Create a lookup from original gene to the annotated label
  label_map <- setNames(top_25_genes$Gene_label, top_25_genes$Gene)
  
  # Add a new column "Gene_label" to data_long
  data_long$Gene_label <- label_map[data_long$Gene]
  
  # Convert to factor with the desired order
  data_long$Gene_label <- factor(data_long$Gene_label, levels = ordered_labels)
  
  # Improved boxplot with better color, grid, and styling
  p <- ggplot(data_long, aes(x = Label, y = Expression, fill = Label)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +       # remove outliers from box
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +      # add jittered points
    facet_wrap(~ Gene_label, scales = "free_y", ncol = 5) + # 5 columns of facets
    scale_fill_manual(values = c("Sepsis" = "#79FFBC", "Control" = "grey")) +
    labs(
      title = "Expression of Top 25 Significant Genes by Sepsis Status",
      subtitle = "Mann-Whitney U test (BH-adjusted p-value < 0.05)",
      x = "Sepsis Status", 
      y = "Expression Level"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
      strip.text = element_text(size = 8, face = "bold", color = "darkgrey"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Display and save the plot
  print(p)
  ggsave("top_25_significant_genes_boxplot_GSE9960.png",
         plot = p, width = 14, height = 10, dpi = 300)
  
} else {
  cat("No significant genes found (adjusted p-value < 0.05).\n")
}

