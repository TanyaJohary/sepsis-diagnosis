################################################################################
# Improved Code: Random Forest with Repeated Split, Feature Removal, 
# and Sanity Check (Same Splits)
################################################################################

# Load necessary libraries
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(randomForest, caret, pROC, dplyr, ggplot2, parallel, doParallel)

# Set up parallel computing
num_cores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(num_cores)
registerDoParallel(cl)

################################################################################
# 1) Load and Prepare Data
################################################################################
data <- read.csv("sepsis_dataGSE243217.csv")
data <- dplyr::select(data, -Sample)  # Remove column 'Sample' if it exists
data$Label <- as.factor(data$Label)   # Ensure Label is a factor

n_repeats <- 100

#set a seed for reproducible splits
#set.seed(123)

# Precompute the same stratified train/test partitions for every iteration
all_partitions <- vector("list", n_repeats)
for (i in 1:n_repeats) {
  all_partitions[[i]] <- createDataPartition(data$Label, p = 0.6, list = FALSE)
}

################################################################################
# 2) Baseline Random Forest (Using the Same Splits)
################################################################################
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
  train_data  <- data[train_index, ]
  test_data   <- data[-train_index, ]
  
  # Train Random Forest model
  rf_model <- randomForest(Label ~ ., data = train_data,
                           ntree = 100,
                           mtry = sqrt(ncol(train_data) - 1))
  
  # Predict on test data
  probs <- predict(rf_model, newdata = test_data, type = "prob")[, 2]
  predictions <- ifelse(probs > 0.5, "Sepsis", "healthy")
  
  # Skip iteration if only one class present in test data or predictions
  if (length(unique(test_data$Label)) < 2 || length(unique(predictions)) < 2) {
    cat("Baseline iteration", i, "skipped due to class imbalance in test set or predictions.\n")
    next
  }
  
  # Calculate metrics
  confusion <- confusionMatrix(factor(predictions, levels = levels(test_data$Label)),
                               test_data$Label, positive = "Sepsis")
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
write.csv(results, "repeated_splits_metrics_GSE243217.csv", row.names = FALSE)

# Average metrics
average_metrics <- colMeans(results, na.rm = TRUE)
cat("Average Metrics (All Features):\n")
print(average_metrics)

# Create a data frame for the average metrics
average_metrics_df <- data.frame(
  Metric = names(average_metrics),
  Value  = average_metrics
)
write.csv(average_metrics_df, "average_metrics_GSE243217.csv", row.names = FALSE)
cat("Average metrics saved to 'average_metrics_GSE243217.csv'\n")

################################################################################
# 3) Feature Removal Analysis (Using the Same Splits)
################################################################################
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
  
  # Create a dataset with the chosen feature removed
  modified_data <- data[, !names(data) %in% feature_to_remove]
  
  # Prepare container for repeated results
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
    train_data  <- modified_data[train_index, ]
    test_data   <- modified_data[-train_index, ]
    
    # Train RF on the modified dataset
    rf_model <- randomForest(Label ~ ., data = train_data,
                             ntree = 100,
                             mtry = sqrt(ncol(train_data) - 1))
    
    probs <- predict(rf_model, newdata = test_data, type = "prob")[, 2]
    predictions <- ifelse(probs > 0.5, "Sepsis", "healthy")
    
    if (length(unique(test_data$Label)) < 2 || length(unique(predictions)) < 2) {
      next
    }
    
    confusion <- confusionMatrix(factor(predictions, levels = levels(test_data$Label)),
                                 test_data$Label, positive = "Sepsis")
    TP <- confusion$table[2, 2]
    FP <- confusion$table[1, 2]
    FN <- confusion$table[2, 1]
    TN <- confusion$table[1, 1]
    
    AUC_val <- ifelse(length(unique(probs)) > 1,
                      pROC::auc(as.numeric(test_data$Label), probs),
                      NA)
    
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
                                   data.frame(
                                     Feature = feature_to_remove,
                                     MCC     = avg_metrics[1],
                                     F1      = avg_metrics[2],
                                     AUC     = avg_metrics[3],
                                     TPR     = avg_metrics[4],
                                     TNR     = avg_metrics[5],
                                     PPV     = avg_metrics[6],
                                     NPV     = avg_metrics[7],
                                     stringsAsFactors = FALSE
                                   ))
}

# Save feature removal results
write.csv(feature_removal_results, "feature_removal_results_GSE243217.csv", row.names = FALSE)

# Plot MCC impact from feature removal
baseline_mcc <- average_metrics["MCC"]
cat("\nBaseline MCC (using all features):", baseline_mcc, "\n\n")

# Combine baseline and feature-removal MCCs into one data frame
mcc_comparison_df <- rbind(
  data.frame(Feature = "All Features (Baseline)", MCC = baseline_mcc),
  feature_removal_results[, c("Feature", "MCC")]
)

# Reorder so that the highest MCC is at the bottom (when we do coord_flip())
mcc_comparison_df$Feature <- reorder(mcc_comparison_df$Feature, -mcc_comparison_df$MCC)

# Create a bar plot
ggplot(mcc_comparison_df, aes(x = Feature, y = MCC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal(base_size = 6) +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text  = element_text(size = 5)
  ) +
  labs(
    title = "MCC for All Features (Baseline) vs. Removal of Each Feature",
    x = "Feature Removed",
    y = "MCC"
  )
ggsave("feature_removal_mcc_plot_GSE243217.png", width = 10, height = 12, dpi = 300)

################################################################################
# 4) Sanity Check: Adding Noise to Data (Using the Same Splits)
################################################################################
noise_levels <- c(0, 10, 20, 30, 40, 50)  # Noise percentages
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
  cat("\nAdding noise at level:", noise_level, "%\n")
  
  # Make a copy of the original data
  noisy_data <- data
  
  # Inject noise
  if (noise_level > 0) {
    num_noisy_rows <- floor(nrow(data) * noise_level / 100)
    noisy_rows <- sample(seq_len(nrow(data)), num_noisy_rows)
    
    # Example: multiply chosen rows by random factor [0.8, 1.2] in numeric columns
    numeric_cols <- which(sapply(noisy_data, is.numeric))
    for (r in noisy_rows) {
      for (col in numeric_cols) {
        noisy_data[r, col] <- noisy_data[r, col] * runif(1, min = 0.8, max = 1.2)
      }
    }
  }
  
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
    
    # Train Random Forest on the noisy data
    rf_model <- randomForest(Label ~ ., data = train_data,
                             ntree = 100,
                             mtry = sqrt(ncol(train_data) - 1))
    
    probs <- predict(rf_model, newdata = test_data, type = "prob")[, 2]
    predictions <- ifelse(probs > 0.5, "Sepsis", "healthy")
    
    if (length(unique(test_data$Label)) < 2 || length(unique(predictions)) < 2) {
      next
    }
    
    confusion <- confusionMatrix(factor(predictions, levels = levels(test_data$Label)),
                                 test_data$Label, positive = "Sepsis")
    TP <- confusion$table[2, 2]
    FP <- confusion$table[1, 2]
    FN <- confusion$table[2, 1]
    TN <- confusion$table[1, 1]
    
    AUC_val <- ifelse(length(unique(probs)) > 1,
                      pROC::auc(as.numeric(test_data$Label), probs),
                      NA)
    
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
  sanity_check_results <- rbind(sanity_check_results,
                                data.frame(
                                  Noise_Level = noise_level,
                                  MCC = avg_metrics[1],
                                  F1  = avg_metrics[2],
                                  AUC = avg_metrics[3],
                                  TPR = avg_metrics[4],
                                  TNR = avg_metrics[5],
                                  PPV = avg_metrics[6],
                                  NPV = avg_metrics[7]
                                ))
}

# Save sanity check results
write.csv(sanity_check_results, "sanity_check_results_GSE243217.csv", row.names = FALSE)

# Plot impact of noise on metrics
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
ggsave("impact_of_noise_on_model_performance_GSE243217.png",
       plot = noise_plot, width = 10, height = 6, dpi = 300)

# Stop parallel cluster
stopCluster(cl)
