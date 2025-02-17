###############################################################################
# COMPLETE SCRIPT: Mann–Whitney vs. Random Forest + Labeling "Significant & High Importance" Genes
###############################################################################

# 1. Load required libraries
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, readr, ggplot2, corrr, rcompanion, ggrepel, tidyr, ggVennDiagram, boot)

###############################################################################
# 2. Read Input Data
###############################################################################
mann_whitney_results <- read_csv("mann_whitney_results_adjusted_D2_GSE54514.csv")
rf_importance        <- read_csv("average_gini_importance_D2_GSE54514.csv")

data <- read.csv("sepsis_D2_GSE54514.csv")
data$Label <- factor(data$Label, levels = c("Control", "Sepsis"))

###############################################################################
# 3. Prepare Mann-Whitney Ranking
###############################################################################
mann_whitney_ranking <- mann_whitney_results %>%
  as.data.frame() %>%
  dplyr::select(Gene, adjusted_p) %>%
  dplyr::arrange(adjusted_p) %>%
  dplyr::mutate(Rank_Stat = dplyr::row_number())  # 1 = most significant

write_csv(mann_whitney_ranking, "mann-w-ranking_D2_GSE154918.csv")

###############################################################################
# 4. Prepare Random Forest Ranking
###############################################################################
rf_ranking <- rf_importance %>%
  as.data.frame() %>%
  dplyr::select(Feature, Mean_Importance) %>%
  dplyr::arrange(desc(Mean_Importance)) %>%
  dplyr::mutate(Rank_RF = dplyr::row_number()) %>%
  dplyr::rename(Gene = Feature)

write_csv(rf_ranking, "rf-ranking-D2_GSE154918.csv")
###############################################################################
# 5. Merge Rankings & Compute Rank Difference
###############################################################################
merged_ranking <- dplyr::inner_join(mann_whitney_ranking, rf_ranking, by = "Gene") %>%
  dplyr::mutate(
    Rank_Diff = abs(Rank_Stat - Rank_RF)
  )

write_csv(merged_ranking, "merged_rankings_comparison_D2_GSE154918.csv")

###############################################################################
# 6. Correlations and Formal Tests
###############################################################################
spearman_corr <- cor(merged_ranking$Rank_Stat, merged_ranking$Rank_RF, method = "spearman")
kendall_corr  <- cor(merged_ranking$Rank_Stat, merged_ranking$Rank_RF, method = "kendall")

cat("\n===== Basic Correlation Estimates =====\n")
cat("Spearman's Rank Correlation: ", spearman_corr, "\n")
cat("Kendall's Tau Correlation:   ", kendall_corr, "\n")

# Formal tests
spearman_test <- cor.test(~ Rank_Stat + Rank_RF, data = merged_ranking, method = "spearman")
kendall_test  <- cor.test(~ Rank_Stat + Rank_RF, data = merged_ranking, method = "kendall")

cat("\n===== Spearman Correlation Test =====\n")
print(spearman_test)
cat("\n===== Kendall's Tau Correlation Test =====\n")
print(kendall_test)

# 6.1 Extract & Save Correlation Test Results to CSV
spearman_df <- data.frame(
  Test       = "Spearman",
  Estimate   = unname(spearman_test$estimate),
  Statistic  = unname(spearman_test$statistic),
  P_Value    = spearman_test$p.value,
  Alternative_Hypothesis = spearman_test$alternative
)

kendall_df <- data.frame(
  Test       = "Kendall",
  Estimate   = unname(kendall_test$estimate),
  Statistic  = unname(kendall_test$statistic),
  P_Value    = kendall_test$p.value,
  Alternative_Hypothesis = kendall_test$alternative
)

cor_test_results <- rbind(spearman_df, kendall_df)
write_csv(cor_test_results, "correlation_test_results_D2_GSE154918.csv")

###############################################################################
# 7. Plot: Mann–Whitney vs. RF Ranks (color by Rank Diff)
###############################################################################
p_comparison <- ggplot(merged_ranking, aes(x = Rank_Stat, y = Rank_RF, color = Rank_Diff)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", color = "#FF4A47") +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black", size = 1) +
  scale_color_gradient(low = "#249699", high = "#FFC868") +
  labs(
    title = "Comparison of Feature Rankings",
    subtitle = paste0("Spearman: ", round(spearman_corr, 3), 
                      " | Kendall's Tau: ", round(kendall_corr, 3)),
    x = "Mann–Whitney Rank (Statistical Test)",
    y = "Random Forest Rank (Feature Importance)",
    color = "Rank Difference"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "darkgray")
  )

ggsave("comparison_ranks_D2_GSE154918.png", p_comparison, width = 6, height = 5, dpi = 300)
print(p_comparison)

###############################################################################
# 8. Histogram of Rank Differences
###############################################################################
p_hist_diff <- ggplot(merged_ranking, aes(x = Rank_Diff)) +
  geom_histogram(bins = 30, fill = "#69b3a2", color = "white", alpha = 0.8) +
  labs(
    title = "Distribution of Rank Differences",
    x = "Absolute Difference in Ranks (|Rank_Stat - Rank_RF|)",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)

ggsave("hist_rank_diff_D2_GSE154918.png", p_hist_diff, width = 6, height = 4, dpi = 300)
print(p_hist_diff)

###############################################################################
# 9. Bar Plot of Top N Genes with Largest Rank Discrepancies
###############################################################################
N <- 10
top_discrepancies <- merged_ranking %>%
  arrange(desc(Rank_Diff)) %>%
  slice_head(n = N)

p_bar_diff <- ggplot(top_discrepancies, aes(x = reorder(Gene, Rank_Diff), y = Rank_Diff)) +
  geom_col(fill = "#ff7f50", alpha = 0.8) +
  coord_flip() +
  labs(
    title = paste0("Top ", N, " Genes with Largest Ranking Discrepancies"),
    x = "Gene",
    y = "Absolute Rank Difference"
  ) +
  theme_minimal(base_size = 14)

ggsave("bar_top_discrepancies_D2_GSE154918.png", p_bar_diff, width = 6, height = 4, dpi = 300)
print(p_bar_diff)

###############################################################################
# 10. Histogram of Mann–Whitney Adjusted p-values
###############################################################################
p_hist_pval <- ggplot(mann_whitney_ranking, aes(x = adjusted_p)) +
  geom_histogram(bins = 30, fill = "#008b8b", color = "white", alpha = 0.8) +
  labs(
    title = "Distribution of Mann–Whitney Adjusted p-values",
    x = "Adjusted p-value",
    y = "Count"
  ) +
  theme_minimal(base_size = 14)

ggsave("hist_mannwhitney_pvals_D2_GSE154918.png", p_hist_pval, width = 6, height = 4, dpi = 300)
print(p_hist_pval)

###############################################################################
# 12. Label Genes That Are Highly Significant & Above the Red Trend Line
###############################################################################
# (a) Fit a linear model to Mean_Importance ~ log10(adjusted_p)
model_lm <- lm(Mean_Importance ~ log10(adjusted_p), data = merged_ranking)

# (b) Create a new column with predicted values from this model
merged_ranking <- merged_ranking %>%
  mutate(predicted_importance = predict(model_lm, newdata = merged_ranking))

# (c) Define our criteria. Example: 
#     - "Highly significant" means adjusted_p < 1e-4
#     - "Above the line" means Mean_Importance > predicted_importance
highlight_genes <- merged_ranking %>%
  filter(adjusted_p < 1e-4, Mean_Importance > predicted_importance)

# (d) Plot with gene labels for the highlight set
p_scatter_pval_importance_labeled <- ggplot(merged_ranking, aes(x = adjusted_p, y = Mean_Importance)) +
  geom_point(alpha = 0.8, color = "#847FE5") +
  geom_smooth(method = "lm", se = FALSE, color = "#FF4867") +
  scale_x_log10() +
  labs(
    title = "Mann–Whitney p-values vs. Random Forest Importance (Labeled Outliers)",
    x = "Adjusted p-value (log scale)",
    y = "Mean Importance (RF)"
  ) +
  theme_minimal(base_size = 14) +
  geom_text_repel(
    data = highlight_genes,
    aes(label = Gene),
    color = "grey",
    fontface = "bold",
    size = 3,
    max.overlaps = 15
  )

ggsave("scatter_pval_importance_labeled_D2_GSE154918.png", p_scatter_pval_importance_labeled, width = 7, height = 5, dpi = 300)
print(p_scatter_pval_importance_labeled)

###############################################################################
# 13. Effect Size Calculations (Using rcompanion::wilcoxonR)
###############################################################################
# 1) Pivot from wide -> long to get columns: [Sample, Label, Gene, Expression].

long_data <- data %>%
  pivot_longer(
    cols = -c(Sample, Label),  # exclude the ID columns, IFNA1 has same value in all samples
    names_to = "Gene",
    values_to = "Expression"
  )

# 2) For each Gene, compute Wilcoxon effect size using wilcoxonR()
#    Ensure Label has exactly 2 levels (e.g., c("Healthy","Sepsis")).
effect_sizes <- long_data %>%
  group_by(Gene) %>%
  summarize(
    # wilcoxonR() returns a list with r (effect size) and confidence interval if ci=TRUE
    res = list(wilcoxonR(Expression, Label, ci=TRUE))
  ) %>%
  rowwise() %>%
  mutate(
    eff       = res$r,      # rank-based effect size
    conf.low  = res$ci[1],  # lower CI
    conf.high = res$ci[2]   # upper CI
  ) %>%
  select(-res)

# 3) Inspect effect size results
print(effect_sizes)

###############################################################################
# 14. Merge Effect Sizes with Merged Rankings
###############################################################################
# we want to incorporate the effect sizes into our main merged_ranking table:

# (a) Rename columns for clarity
effect_sizes <- effect_sizes %>%
  rename(
    EffectSize = eff
    # LowerCI   = conf.low,
    # UpperCI   = conf.high
  )

# (b) Join to merged_ranking
final_merged <- merged_ranking %>%
  left_join(effect_sizes, by = "Gene")

# (c) Save final table
write_csv(final_merged, "final_merged_with_effect_sizes_D2_GSE154918.csv")

# (d) Quick check
colnames(final_merged)
head(final_merged)

###############################################################################
# (2) TOP-N OVERLAP / VENN DIAGRAM
###############################################################################
# Compare top 10 from Mann–Whitney vs. top 10 from RF
top_10_mw <- mann_whitney_ranking %>%
  slice_head(n = 10) %>%
  pull(Gene)

top_10_rf <- rf_ranking %>%
  slice_head(n = 10) %>%
  pull(Gene)

venn_sets <- list(
  "Top 10 MW" = top_10_mw,
  "Top 10 RF" = top_10_rf
)

p_venn <- ggVennDiagram(venn_sets, label_alpha = 0) +
  scale_fill_gradient(low = "#FEEBD0", high = "#F0A789") +
  labs(title = "Overlap of Top 10 Genes: MW vs. RF")

ggsave("venn_top10_mw_rf_D2_GSE154918.png", p_venn, width = 5, height = 4, dpi = 300)
print(p_venn)


# 3) Find overlap (intersection) of these two vectors
overlap_genes <- intersect(top_10_mw, top_10_rf)

# 4) Print them
cat("Overlap Genes:\n")
print(overlap_genes)

# 5)  Save to a CSV
write.csv(data.frame(OverlapGenes = overlap_genes),
          "top10_overlap_genes_D2_GSE154918.csv",
          row.names = FALSE)

###############################################################################
# (7) BOOTSTRAPPED SPEARMAN CORRELATION
###############################################################################
# Example: Bootstrapping the correlation between Rank_Stat and Rank_RF
# using the 'boot' package. We'll create a function returning spearman corr.

spearman_func <- function(data, indices) {
  d <- data[indices, ]
  cor(d$Rank_Stat, d$Rank_RF, method = "spearman")
}

set.seed(123)  # for reproducibility
boot_res <- boot(
  data = merged_ranking,
  statistic = spearman_func,
  R = 1000  # number of bootstrap samples
)

# Basic summary
cat("\nBootstrap Mean of Spearman:", mean(boot_res$t), "\n")
cat("Original Spearman correlation:", boot_res$t0, "\n")

# Confidence interval (percentile method)
boot_ci <- boot.ci(boot_res, type = "perc")
cat("Bootstrap 95% CI for Spearman correlation (percentile):\n")
print(boot_ci$percent[4:5])  # lower, upper bounds

write.csv(boot_res$t, "boot_spearman_samples_D2_GSE54514.csv", row.names=FALSE)

###############################################################################
# End of Script
###############################################################################
