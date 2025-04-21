data <- read.csv("sepsis_dataGSE185263.csv")

topgenes <- data %>% 
  dplyr::select("Sample", "CD177", "FCGR1A", "C3AR1", "ITGAM", "S100A12", "S100A9", "S100A8",
                "ARG1", "SOCS3", "MAPK14", "BCL2", "IL10", "IL1R2", "GATA3", "MMP9", "Label")
write.csv(topgenes, "topgenesGSE185263.csv", row.names = FALSE)


cluster1 <- data %>% 
  dplyr::select("Sample", "S100A12", "S100A9", "S100A8", "Label")
write.csv(cluster1, "cluster1_GSE185263.csv", row.names = FALSE)


cluster2 <- data %>% 
  dplyr::select("Sample", "CD177", "FCGR1A", "C3AR1", "ITGAM", "Label")
write.csv(cluster2, "cluster2_GSE185263.csv", row.names = FALSE)

cluster3 <- data %>% 
  dplyr::select("Sample", "ARG1", "SOCS3", "MAPK14", "BCL2", "IL10", "IL1R2", "GATA3", "MMP9", "Label")
write.csv(cluster3, "cluster3_GSE185263.csv", row.names = FALSE)
