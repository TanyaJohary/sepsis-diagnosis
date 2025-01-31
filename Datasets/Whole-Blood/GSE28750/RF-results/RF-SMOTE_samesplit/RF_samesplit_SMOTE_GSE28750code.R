
################################################################################
############################# SMOTE RESAMPLING #################################
################################################################################

################################################################################
# Improved Code: Random Forest with SMOTE, Noise Analysis, Feature Removal, and Robust Evaluation
################################################################################

# Load necessary libraries
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(randomForest, caret, pROC, dplyr, ggplot2, parallel, doParallel, smotefamily, ROSE)

# Set up parallel computing
num_cores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(num_cores)
registerDoParallel(cl)

################################################################################
# Load and Prepare Data
################################################################################
data <- read.csv("sepsis_data_GSE28750.csv")
data <- dplyr::select(data, -Sample)
data <- data %>%
  mutate(Label = as.factor(Label))  # Convert the target variable to a factor

###############################################################################
# 1) Pre-compute the partitions for all repeats
###############################################################################

n_repeats <- 100
# Create a list of partitions so that each iteration i uses the same split
all_partitions <- vector("list", n_repeats)
for (i in 1:n_repeats) {
  all_partitions[[i]] <- createDataPartition(data$Label, p = 0.8, list = FALSE)
}



###############################################################################
# 2) Baseline Random Forest with SMOTE, using the same precomputed splits
###############################################################################
results <- data.frame(
  MCC = numeric(n_repeats),
  F1  = numeric(n_repeats),
  AUC = numeric(n_repeats),
  TPR = numeric(n_repeats),
  TNR = numeric(n_repeats),
  PPV = numeric(n_repeats),
  NPV = numeric(n_repeats)
)

for (i in 1:n_repeats) {
  train_index <- all_partitions[[i]]
  
  # Split
  train_data <- data[train_index, ]
  test_data  <- data[-train_index, ]
  
  # Apply SMOTE
  smote_train_data <- SMOTE(train_data[, -ncol(train_data)], train_data$Label, 
                            K = 5, dup_size = 2)$data
  colnames(smote_train_data)[ncol(smote_train_data)] <- "Label"
  smote_train_data$Label <- as.factor(smote_train_data$Label)
  
  # Train Random Forest
  rf_model <- randomForest(Label ~ ., data = smote_train_data,
                           ntree = 100,
                           mtry  = sqrt(ncol(smote_train_data) - 1))
  
  # Predictions
  probs <- predict(rf_model, newdata = test_data, type = "prob")[, 2]
  predictions <- ifelse(probs > 0.5, "SEPSIS", "HEALTHY")
  
  # Skip if not enough classes present
  if (length(unique(test_data$Label)) < 2 || length(unique(predictions)) < 2) {
    cat("Iteration", i, "skipped due to class imbalance in test set or predictions.\n")
    next
  }
  
  # Metrics
  confusion <- confusionMatrix(factor(predictions, levels = levels(test_data$Label)), 
                               test_data$Label, positive = "SEPSIS")
  TP <- confusion$table[2, 2]
  FP <- confusion$table[1, 2]
  FN <- confusion$table[2, 1]
  TN <- confusion$table[1, 1]
  
  AUC_val <- ifelse(length(unique(probs)) > 1,
                    pROC::auc(as.numeric(test_data$Label), probs), NA)
  
  MCC_val <- ifelse(((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)) > 0,
                    ((TP * TN) - (FP * FN)) / 
                      sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)),
                    NA)
  
  results[i, ] <- c(MCC_val,
                    confusion$byClass["F1"],
                    AUC_val,
                    confusion$byClass["Sensitivity"],
                    confusion$byClass["Specificity"],
                    confusion$byClass["Pos Pred Value"],
                    confusion$byClass["Neg Pred Value"])
}

# Save detailed results
write.csv(results, "repeated_splits_SMOTE_GSE28750.csv", row.names = FALSE)

# Average baseline metrics
average_metrics <- colMeans(results, na.rm = TRUE)
cat("Baseline (All Features) Average Metrics with SMOTE:\n")
print(average_metrics)


for (i in seq_along(average_metrics)) {
  cat(names(average_metrics)[i], ":", average_metrics[i], "\n")
}

# Create a data frame for the average metrics
average_metrics_df <- data.frame(
  Metric = names(average_metrics),
  Value = average_metrics
)

# Save the data frame to a CSV file
write.csv(average_metrics_df, "average_metrics_SMOTE_GSE28750.csv", row.names = FALSE)

# Print confirmation
cat("Average metrics saved to 'average_metrics_SMOTE_GSE28750.csv'\n")

###############################################################################
# 3) Feature Removal Analysis, re-using the SAME partitions
###############################################################################
features <- setdiff(names(data), "Label")
feature_removal_results <- data.frame(
  Feature = character(),
  MCC = numeric(),
  F1  = numeric(),
  AUC = numeric(),
  TPR = numeric(),
  TNR = numeric(),
  PPV = numeric(),
  NPV = numeric(),
  stringsAsFactors = FALSE
)

for (feature_to_remove in features) {
  cat("Removing feature:", feature_to_remove, "\n")
  
  # Remove 1 feature from the dataset
  modified_data <- data[, !names(data) %in% feature_to_remove]
  
  # Container for repeated splits
  metrics <- data.frame(
    MCC = numeric(n_repeats),
    F1  = numeric(n_repeats),
    AUC = numeric(n_repeats),
    TPR = numeric(n_repeats),
    TNR = numeric(n_repeats),
    PPV = numeric(n_repeats),
    NPV = numeric(n_repeats)
  )
  
  for (i in 1:n_repeats) {
    train_index <- all_partitions[[i]]
    
    train_data <- modified_data[train_index, ]
    test_data  <- modified_data[-train_index, ]
    
    # SMOTE on the new training set
    smote_train_data <- SMOTE(train_data[, -ncol(train_data)], train_data$Label,
                              K = 5, dup_size = 2)$data
    colnames(smote_train_data)[ncol(smote_train_data)] <- "Label"
    smote_train_data$Label <- as.factor(smote_train_data$Label)
    
    # Train RF
    rf_model <- randomForest(Label ~ ., data = smote_train_data,
                             ntree = 100,
                             mtry  = sqrt(ncol(smote_train_data) - 1))
    
    probs <- predict(rf_model, newdata = test_data, type = "prob")[, 2]
    predictions <- ifelse(probs > 0.5, "SEPSIS", "HEALTHY")
    
    if (length(unique(test_data$Label)) < 2 || length(unique(predictions)) < 2) {
      next
    }
    
    confusion <- confusionMatrix(factor(predictions, levels = levels(test_data$Label)), 
                                 test_data$Label, positive = "SEPSIS")
    TP <- confusion$table[2, 2]
    FP <- confusion$table[1, 2]
    FN <- confusion$table[2, 1]
    TN <- confusion$table[1, 1]
    
    AUC_val <- ifelse(length(unique(probs)) > 1,
                      pROC::auc(as.numeric(test_data$Label), probs), NA)
    
    MCC_val <- ifelse(((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)) > 0,
                      ((TP * TN) - (FP * FN)) / 
                        sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)),
                      NA)
    
    metrics[i, ] <- c(MCC_val,
                      confusion$byClass["F1"],
                      AUC_val,
                      confusion$byClass["Sensitivity"],
                      confusion$byClass["Specificity"],
                      confusion$byClass["Pos Pred Value"],
                      confusion$byClass["Neg Pred Value"])
  }
  
  avg_metrics <- colMeans(metrics, na.rm = TRUE)
  feature_removal_results <- rbind(feature_removal_results, 
                                   data.frame(Feature = feature_to_remove, 
                                              MCC     = avg_metrics[1],
                                              F1      = avg_metrics[2],
                                              AUC     = avg_metrics[3],
                                              TPR     = avg_metrics[4],
                                              TNR     = avg_metrics[5],
                                              PPV     = avg_metrics[6],
                                              NPV     = avg_metrics[7],
                                              stringsAsFactors = FALSE))
}


# Save feature removal results
write.csv(feature_removal_results, "feature_removal_results_SMOTE_GSE28750.csv", row.names = FALSE)

###############################################################################
# 4) Compare "All Features" vs. "Removed One Feature" Using the Same Splits
###############################################################################

# Plot MCC impact from feature removal
baseline_mcc <- average_metrics["MCC"]
cat("\nBaseline MCC (using all features):", baseline_mcc, "\n\n")

# 2) Combine baseline and feature-removal MCCs into one data frame
mcc_comparison_df <- rbind(
  data.frame(Feature = "All Features (Baseline)", MCC = baseline_mcc),
  feature_removal_results[, c("Feature", "MCC")]
)

# 3) Reorder factors so that the largest MCC is displayed at the bottom of the plot
#    (because geom_bar + coord_flip will place the first level at the bottom).
mcc_comparison_df$Feature <- reorder(mcc_comparison_df$Feature, -mcc_comparison_df$MCC)

# 4) Create a bar plot of MCC for baseline vs. each feature removed
ggplot(mcc_comparison_df, aes(x = Feature, y = MCC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal(base_size = 6) +  # Set a larger "base" text size
  theme(
    # customize each text element 
    plot.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 5),
    # legend.text = element_text(size = 12), etc.
  ) +
  labs(
    title = "MCC for All Features (Baseline) vs. Removal of Each Feature",
    x = "Feature Removed",
    y = "MCC"
  )
ggsave("feature_removal_mcc_plot_SMOTE_GSE28750.png", width = 10, height = 12, dpi = 300)

################################################################################
# 4) Sanity Check: Adding Noise + Using Same Splits
################################################################################

noise_levels <- c(0, 10, 20, 30, 40, 50)
sanity_check_results <- data.frame(
  Noise_Level = integer(),
  MCC = numeric(),
  F1  = numeric(),
  AUC = numeric(),
  TPR = numeric(),
  TNR = numeric(),
  PPV = numeric(),
  NPV = numeric()
)

for (noise_level in noise_levels) {
  cat("Adding noise at level:", noise_level, "%\n")
  
  # Create a noisy copy of the data
  noisy_data <- data
  if (noise_level > 0) {
    num_noisy_rows <- floor(nrow(data) * noise_level / 100)
    noisy_rows <- sample(seq_len(nrow(data)), num_noisy_rows)
    
    # Simple approach: multiply the chosen rows' numeric features by ~1 ± small factor
    # Adjust as needed for your data type
    numeric_cols <- which(sapply(noisy_data, is.numeric))
    for (r in noisy_rows) {
      for (col in numeric_cols) {
        noisy_data[r, col] <- noisy_data[r, col] * runif(1, min = 0.8, max = 1.2)
      }
    }
  }
  
  # We'll run the same n_repeats using the *same partitions*
  metrics <- data.frame(
    MCC = numeric(n_repeats),
    F1  = numeric(n_repeats),
    AUC = numeric(n_repeats),
    TPR = numeric(n_repeats),
    TNR = numeric(n_repeats),
    PPV = numeric(n_repeats),
    NPV = numeric(n_repeats)
  )
  
  for (i in 1:n_repeats) {
    train_index <- all_partitions[[i]]
    train_data  <- noisy_data[train_index, ]
    test_data   <- noisy_data[-train_index, ]
    
    # SMOTE on the noisy training
    smote_train_data <- SMOTE(train_data[, -ncol(train_data)], train_data$Label, 
                              K = 5, dup_size = 2)$data
    colnames(smote_train_data)[ncol(smote_train_data)] <- "Label"
    smote_train_data$Label <- as.factor(smote_train_data$Label)
    
    rf_model <- randomForest(Label ~ ., data = smote_train_data, ntree = 100,
                             mtry = sqrt(ncol(smote_train_data) - 1))
    
    probs <- predict(rf_model, newdata = test_data, type = "prob")[, 2]
    predictions <- ifelse(probs > 0.5, "SEPSIS", "HEALTHY")
    
    if (length(unique(test_data$Label)) < 2 || length(unique(predictions)) < 2) {
      next
    }
    
    confusion <- confusionMatrix(factor(predictions, levels = levels(test_data$Label)),
                                 test_data$Label, positive = "SEPSIS")
    TP <- confusion$table[2, 2]
    FP <- confusion$table[1, 2]
    FN <- confusion$table[2, 1]
    TN <- confusion$table[1, 1]
    
    AUC_val <- ifelse(length(unique(probs)) > 1,
                      auc(as.numeric(test_data$Label), probs), NA)
    
    MCC_val <- ifelse(((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)) > 0,
                      ((TP * TN) - (FP * FN)) / 
                        sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)),
                      NA)
    
    metrics[i, ] <- c(MCC_val,
                      confusion$byClass["F1"],
                      AUC_val,
                      confusion$byClass["Sensitivity"],
                      confusion$byClass["Specificity"],
                      confusion$byClass["Pos Pred Value"],
                      confusion$byClass["Neg Pred Value"])
  }
  
  avg_metrics <- colMeans(metrics, na.rm = TRUE)
  sanity_check_results <- rbind(
    sanity_check_results,
    data.frame(
      Noise_Level = noise_level,
      MCC = avg_metrics[1],
      F1  = avg_metrics[2],
      AUC = avg_metrics[3],
      TPR = avg_metrics[4],
      TNR = avg_metrics[5],
      PPV = avg_metrics[6],
      NPV = avg_metrics[7]
    )
  )
}

write.csv(sanity_check_results, "sanity_check_results_SMOTE_GSE28750.csv", row.names = FALSE)

# Plot
noise_plot <- ggplot(sanity_check_results, aes(x = Noise_Level)) +
  geom_line(aes(y = MCC, color = "MCC"), size = 1.2) +
  geom_line(aes(y = F1, color = "F1"), size = 1.2) +
  geom_line(aes(y = AUC, color = "AUC"), size = 1.2) +
  labs(
    title = "Impact of Noise on Model Performance (Same Splits)",
    x = "Noise Level (%)",
    y = "Metric Value",
    color = "Metrics"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("MCC" = "blue", "F1" = "green", "AUC" = "red"))

print(noise_plot)
ggsave("impact_of_noise_plot_SMOTE_GSE28750.png", plot = noise_plot, width = 10, height = 6, dpi = 300)

# Stop parallel cluster
stopCluster(cl)

################################################################################
# End
################################################################################