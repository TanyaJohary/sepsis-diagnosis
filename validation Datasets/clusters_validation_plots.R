# Install and load required package
install.packages("fmsb")
library(fmsb)

# Create the data frame
data <- data.frame(
  MCC    = c(1.0, 0.0, 0.804, 0.723, 0.849, 0.872),
  F1     = c(1.0, 0.0, 0.900, 0.887, 0.943, 0.943),
  AUROC  = c(1.0, 0.0, 0.969, 0.944, 0.967, 0.987),
  TPR    = c(1.0, 0.0, 0.877, 0.853, 0.925, 0.918),
  TNR    = c(1.0, 0.0, 0.921, 0.874, 0.935, 0.961),
  PPV    = c(1.0, 0.0, 0.953, 0.946, 0.972, 0.981),
  NPV    = c(1.0, 0.0, 0.870, 0.871, 0.872, 0.866)
)

# Add row names
rownames(data) <- c("Max", "Min", "Cluster 1", "Cluster 2", "Cluster 3", "Full Panel")

# Set plot colors
plot_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")

# Open PNG device to save the plot
png("radar_plot_cluster_performance.png", width = 1800, height = 1800, res = 300)

# Plot the radar chart
radarchart(data,
           axistype = 1,
           pcol = plot_colors,
           plty = 1,
           plwd = 3,
           cglcol = "grey60",
           cglty = 1,
           cglwd = 0.8,
           axislabcol = "black",
           caxislabels = seq(0, 1, 0.2),
           vlcex = 1.1,
           title = "Diagnostic Performance of Gene Clusters")

# Add legend
legend("topright", legend = rownames(data)[3:6],
       col = plot_colors, lty = 1, lwd = 3, bty = "n", cex = 0.8)

# Close PNG device
dev.off()


###############################################################################
############################ cluster 1 results ################################
###############################################################################

# Create the performance matrix for Cluster 1
cluster1_data <- data.frame(
  MCC   = c(0.8586, 0.6818, 0.8773),
  F1    = c(0.9603, 0.7724, 0.9687),
  AUROC = c(0.9896, 0.9315, 0.9871),
  TPR   = c(0.9437, 0.7225, 0.9657),
  TNR   = c(0.936,  0.922,  0.905),
  PPV   = c(0.9808, 0.9016, 0.976),
  NPV   = c(0.8606, 0.8284, 0.92)
)

# Add the required max and min rows
cluster1_data <- rbind(
  rep(1.0, 7),   # Max
  rep(0.0, 7),   # Min
  cluster1_data
)

# Add row names
rownames(cluster1_data) <- c("Max", "Min", "GSE69528", "GSE60424", "GSE63311")

# Colors for each dataset
colors <- c("#1b9e77", "#d95f02", "#7570b3")

# Save as high-resolution PNG
png("cluster1_radar_plot.png", width = 2000, height = 1500, res = 300)

# Plot radar chart
radarchart(cluster1_data,
           axistype = 1,
           pcol = colors,
           plty = 1,
           plwd = 3,
           cglcol = "grey60",
           cglty = 1,
           cglwd = 0.8,
           axislabcol = "black",
           caxislabels = seq(0, 1, 0.2),
           vlcex = 1.1,
           title = "Cluster 1 – Diagnostic Performance Across Validation Datasets")

# Add legend
legend("topright", legend = rownames(cluster1_data)[3:5],
       col = colors, lty = 1, lwd = 3, bty = "n", cex = 0.8)

# Close PNG device
dev.off()

###############################################################################
####################### cluster 2 results #####################################
###############################################################################


# Cluster 2 performance data
cluster2_data <- data.frame(
  MCC   = c(0.8354, 0.745,  0.5901),
  F1    = c(0.953,  0.8254, 0.8821),
  AUROC = c(0.992,  0.9475, 0.8928),
  TPR   = c(0.9318, 0.775,  0.8527),
  TNR   = c(0.93,   0.938,  0.7551),
  PPV   = c(0.9786, 0.9265, 0.9333),
  NPV   = c(0.8349, 0.8559, 0.6537)
)

# Add max and min reference rows
cluster2_data <- rbind(
  rep(1.0, 7),   # Max
  rep(0.0, 7),   # Min
  cluster2_data
)

# Label rows
rownames(cluster2_data) <- c("Max", "Min", "GSE69528", "GSE60424", "GSE63311")

# Define colors for lines
colors <- c("#1b9e77", "#d95f02", "#7570b3")

# Save as high-res PNG
png("cluster2_radar_plot.png", width = 2000, height = 1500, res = 300)

# Plot radar chart
radarchart(cluster2_data,
           axistype = 1,
           pcol = colors,
           plty = 1,
           plwd = 3,
           cglcol = "grey60",
           cglty = 1,
           cglwd = 0.8,
           axislabcol = "black",
           caxislabels = seq(0, 1, 0.2),
           vlcex = 1.1,
           title = "Cluster 2 – Diagnostic Performance Across Validation Datasets")

# Add legend
legend("topright", legend = rownames(cluster2_data)[3:5],
       col = colors, lty = 1, lwd = 3, bty = "n", cex = 0.8)

# Close device
dev.off()

###############################################################################
########################### cluster 3 results #################################
###############################################################################

# Cluster 3 performance data
cluster3_data <- data.frame(
  MCC   = c(0.8672, 0.8934, 0.7865),
  F1    = c(0.9649, 0.9302, 0.934),
  AUROC = c(0.9796, 0.973,  0.9485),
  TPR   = c(0.9543, 0.9175, 0.9028),
  TNR   = c(0.926,  0.964,  0.915),
  PPV   = c(0.9776, 0.961,  0.977),
  NPV   = c(0.8797, 0.9458, 0.7896)
)

# Add max and min reference rows
cluster3_data <- rbind(
  rep(1.0, 7),   # Max
  rep(0.0, 7),   # Min
  cluster3_data
)

# Label the rows
rownames(cluster3_data) <- c("Max", "Min", "GSE69528", "GSE60424", "GSE63311")

# Set distinct colors
colors <- c("#1b9e77", "#d95f02", "#7570b3")

# Save high-resolution PNG
png("cluster3_radar_plot.png", width = 2000, height = 1500, res = 300)

# Generate radar plot
radarchart(cluster3_data,
           axistype = 1,
           pcol = colors,
           plty = 1,
           plwd = 3,
           cglcol = "grey60",
           cglty = 1,
           cglwd = 0.8,
           axislabcol = "black",
           caxislabels = seq(0, 1, 0.2),
           vlcex = 1.1,
           title = "Cluster 3 – Diagnostic Performance Across Validation Datasets")

# Add legend
legend("topright", legend = rownames(cluster3_data)[3:5],
       col = colors, lty = 1, lwd = 3, bty = "n", cex = 0.8)

# Close PNG device
dev.off()
