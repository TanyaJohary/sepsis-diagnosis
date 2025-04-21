#  we already have:
# rf_gene_freq: data.frame with GeneID, Count_RF
rf_gene_freq <- read.csv("gene_frequency_methodA.csv")

# mwu_gene_freq: data.frame with GeneID, Count_MW
mw_gene_freq <- read.csv("gene_frequency_MWU_methodA.csv")

# combine the frequency n both methods
combined_freq <- full_join(rf_gene_freq, mw_gene_freq, by="GeneID") %>%
  rename(
    Count_RF = Count.x,
    Count_MW = Count.y
  ) %>% 
  mutate(
    Count_RF = ifelse(is.na(Count_RF), 0, Count_RF),
    Count_MW = ifelse(is.na(Count_MW), 0, Count_MW),
    Combined_Freq = Count_RF + Count_MW
  ) %>%
  arrange(desc(Combined_Freq))

# Keep genes with Combined_Freq >= 5 (or top 50 genes)
final_genes <- combined_freq %>% filter(Combined_Freq >= 5)
write.csv(final_genes,"joint_freq.csv")
