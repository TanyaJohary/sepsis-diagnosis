# Create a list of character vectors, one per dataset (or line) overlapped betweem Random forest and Mann-W U test:
top_overlap <- list(
  # 1st line (10 genes)
  c("FCGR1A","S100A9","S100A12","S100A8","GATA3","IL10","CD177","MYD88","TNFSF10","C3AR1"),
  
  # 2nd line (9 genes)
  c("CD177","GATA3","S100A12","S100A8","ARG1","IL1R2","MAPK14","BCL2","CX3CR1"),
  
  # 3rd line (10 genes)
  c("ARG1","GATA3","IL10","MMP9","S100A8","CD177","S100A12","SOCS3","TLR2","HLA-DRA"),
  
  # 4th line (8 genes)
  c("CD177","GATA3","HLA-DRA","IL10","MAPK14","MMP9","S100A12","S100A8"),
  
  # 5th line (9 genes)
  c("FCGR1A","ARG1","CD177","S100A8","S100A9","S100A12","SOCS3","C3AR1","MAPK14"),
  
  # 6th line (8 genes)
  c("MAPK14","S100A8","S100A12","IL1R2","SOCS3","ARG1","CD177","S100A9"),
  
  # 7th line (7 genes)
  c("CD177","GATA3","ITGAM","MAPK14","MMP8","S100A12","S100A9"),
  
  # 8th line (8 genes)
  c("CD177","ITGAM","MMP8","MMP9","MYD88","S100A12","S100A8","S100A9"),
  
  # 9th line (7 genes)
  c("ARG1","CD177","GATA3","ITGAM","MAPK14","MMP9","S100A12"),
  
  # 10th line (9 genes)
  c("C3AR1","CCR7","CD177","ELANE","MMP8","S100A12","S100A8","S100A9")
)

# Now 'top_overlap' is a list with 10 elements (one per dataset/line).
# Combine all genes into one character vector
all_genes <- unlist(top_overlap)

# Count how many times each gene appears
freq_table <- table(all_genes)

# Sort in decreasing order 
freq_table <- sort(freq_table, decreasing = TRUE)

freq_df <- as.data.frame(freq_table)
colnames(freq_df) <- c("Gene", "Count")
freq_df

write.csv(freq_df, "top_overlap_frequency.csv", row.names=FALSE)

