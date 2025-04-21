# Read the final merged ranking file that contains the effect size for each gene
final_merged <- read.csv("~/sepsis/effect_size/final_merged_with_effect_sizes_GSE54514.csv", sep = " ")

# Compute the median and top 25% (75th percentile) of predicted importance
median_pred_imp <- median(final_merged$predicted_importance, na.rm=TRUE)
top25_pred_imp <- quantile(final_merged$predicted_importance, 0.75, na.rm=TRUE)

print(median_pred_imp)  # Above this = above-average importance # 0.3259046
print(top25_pred_imp)   # Above this = top 25% most important genes # 0.4594151   

# Filter genes in the top 25% most important
top_genes <- final_merged %>%
  filter(predicted_importance > top25_pred_imp)

# Save the results into a new CSV file
write.csv(top_genes, "top_25_percent_genes_GSE54514.csv", row.names = FALSE, quote = FALSE)

# View the filtered results
View(top_genes)
