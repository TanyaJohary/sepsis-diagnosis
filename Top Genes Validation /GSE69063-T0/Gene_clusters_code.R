library(dplyr)


data <- read.csv("T0_sepsis_GSE69063.csv")

# Define genes of interest
genes_of_interest <- c(
  "CD177","FCGR1A","C3AR1","ITGAM","S100A12","S100A8","S100A9",
  "ARG1","SOCS3","MAPK14","BCL2","IL10","IL1R2","GATA3","MMP9"
)

# Build a regex pattern like: "^(CD177|FCGR1A|C3AR1|...)"
pattern <- paste0("^(", paste(genes_of_interest, collapse = "|"), ")")

topgenes <- data %>%
  dplyr::select(
    Sample, 
    matches(pattern), # any column whose name starts with one of the listed genes
    Label
  )

write.csv(topgenes, "topgenes_T0_GSE69063.csv",row.names = FALSE)

################################################################### cluster1
# Define genes of interest
genes_of_interest <- c(
  "ARG1", "SOCS3", "MAPK14", "BCL2", "IL10", "IL1R2", "GATA3", "MMP9"
)

# Build a regex pattern like: "^(CD177|FCGR1A|C3AR1|...)"
pattern <- paste0("^(", paste(genes_of_interest, collapse = "|"), ")")

cluster1 <- data %>%
  dplyr::select(
    Sample, 
    matches(pattern), # any column whose name starts with one of the listed genes
    Label
  )

write.csv(cluster1, "cluster1_T0_GSE69063.csv",row.names = FALSE)
##################################################################### cluster2

# Define genes of interest
genes_of_interest <- c("CD177", "FCGR1A", "C3AR1", "ITGAM") # FCGR1A missing

# Build a regex pattern like: "^(CD177|FCGR1A|C3AR1|...)"
pattern <- paste0("^(", paste(genes_of_interest, collapse = "|"), ")")

cluster2 <- data %>%
  dplyr::select(
    Sample, 
    matches(pattern), # any column whose name starts with one of the listed genes
    Label
  )

write.csv(cluster2, "cluster2_T0_GSE69063.csv",row.names = FALSE)
##################################################################### cluster3

# Define genes of interest
genes_of_interest <- c("S100A12", "S100A8", "S100A9")

# Build a regex pattern like: "^(CD177|FCGR1A|C3AR1|...)"
pattern <- paste0("^(", paste(genes_of_interest, collapse = "|"), ")")

cluster3 <- data %>%
  dplyr::select(
    Sample, 
    matches(pattern), # any column whose name starts with one of the listed genes
    Label
  )

write.csv(cluster3, "cluster3_T0_GSE69063.csv",row.names = FALSE)

