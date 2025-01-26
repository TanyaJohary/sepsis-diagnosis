################################################################################
# Improved Code: Random Forest with Repeated Split, Feature Removal, and Sanity Check
################################################################################

# Load necessary libraries
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(randomForest, caret, pROC, dplyr, ggplot2, parallel, doParallel)

# Set up parallel computing
num_cores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(num_cores)
registerDoParallel(cl)

################################################################################
# Load and Prepare Data
################################################################################
data <- read.csv("sepsis_dataGSE100159.csv")
data <- dplyr::select(data, -Sample)
data <- data %>%
  mutate(Label = as.factor(Label))  # Convert the target variable to a factor

################################################################################
# Repeated Train-Test Splits with Random Forest
################################################################################
n_repeats <- 100
results <- data.frame(
  MCC = numeric(n_repeats),
  F1 = numeric(n_repeats),
  AUC = numeric(n_repeats),
  TPR = numeric(n_repeats),
  TNR = numeric(n_repeats),
  PPV = numeric(n_repeats),
  NPV = numeric(n_repeats)
)

for (i in 1:n_repeats) {
  # Stratified 80/20 split
  train_index <- createDataPartition(data$Label, p = 0.8, list = FALSE)
  train_data <- data[train_index, ]
  test_data <- data[-train_index, ]
  
  # Train Random Forest model
  rf_model <- randomForest(Label ~ ., data = train_data, ntree = 100, mtry = sqrt(ncol(train_data) - 1))
  
  # Predict on test data
  probs <- predict(rf_model, newdata = test_data, type = "prob")[, 2]
  predictions <- ifelse(probs > 0.5, "Sepsis", "Control")
  
  # Skip iteration if only one class present in test data or predictions
  if (length(unique(test_data$Label)) < 2 || length(unique(predictions)) < 2) {
    cat("Iteration", i, "skipped due to class imbalance in test set or predictions.\n")
    next
  }
  
  # Calculate metrics
  confusion <- confusionMatrix(factor(predictions, levels = levels(test_data$Label)), 
                               test_data$Label, positive = "Sepsis")
  TP <- confusion$table[2, 2]
  FP <- confusion$table[1, 2]
  FN <- confusion$table[2, 1]
  TN <- confusion$table[1, 1]
  AUC <- ifelse(length(unique(probs)) > 1, auc(as.numeric(test_data$Label), probs), NA)
  MCC <- ifelse((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) > 0, 
                ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)), NA)
  results[i, ] <- c(MCC, confusion$byClass["F1"], AUC, confusion$byClass["Sensitivity"], 
                    confusion$byClass["Specificity"], confusion$byClass["Pos Pred Value"], 
                    confusion$byClass["Neg Pred Value"])
}

# Save detailed results
write.csv(results, "repeated_splits_metrics100159.csv", row.names = FALSE)
# Average metrics
average_metrics <- colMeans(results, na.rm = TRUE)
cat("Average Metrics:\n")
for (i in seq_along(average_metrics)) {
  cat(names(average_metrics)[i], ":", average_metrics[i], "\n")
}

# Create a data frame for the average metrics
average_metrics_df <- data.frame(
  Metric = names(average_metrics),
  Value = average_metrics
)

# Save the data frame to a CSV file
write.csv(average_metrics_df, "average_metrics100159.csv", row.names = FALSE)

# Print confirmation
cat("Average metrics saved to 'average_metrics100159.csv'\n")
###############################################################################
# Feature Removal Analysis
################################################################################
features <- setdiff(names(data), "Label")
feature_removal_results <- data.frame(
  Feature = character(),
  MCC = numeric(),
  F1 = numeric(),
  AUC = numeric(),
  TPR = numeric(),
  TNR = numeric(),
  PPV = numeric(),
  NPV = numeric(),
  stringsAsFactors = FALSE
)

for (feature_to_remove in features) {
  cat("Removing feature:", feature_to_remove, "\n")
  
  # Create a temporary dataset with the feature removed
  modified_data <- data[, !names(data) %in% feature_to_remove]
  
  metrics <- data.frame(
    MCC = numeric(n_repeats),
    F1 = numeric(n_repeats),
    AUC = numeric(n_repeats),
    TPR = numeric(n_repeats),
    TNR = numeric(n_repeats),
    PPV = numeric(n_repeats),
    NPV = numeric(n_repeats)
  )
  
  for (i in 1:n_repeats) {
    train_index <- createDataPartition(modified_data$Label, p = 0.8, list = FALSE)
    train_data <- modified_data[train_index, ]
    test_data <- modified_data[-train_index, ]
    
    rf_model <- randomForest(Label ~ ., data = train_data, ntree = 100, mtry = sqrt(ncol(train_data) - 1))
    probs <- predict(rf_model, newdata = test_data, type = "prob")[, 2]
    predictions <- ifelse(probs > 0.5, "Sepsis", "Control")
    
    if (length(unique(test_data$Label)) < 2 || length(unique(predictions)) < 2) {
      next
    }
    
    confusion <- confusionMatrix(factor(predictions, levels = levels(test_data$Label)), 
                                 test_data$Label, positive = "Sepsis")
    TP <- confusion$table[2, 2]
    FP <- confusion$table[1, 2]
    FN <- confusion$table[2, 1]
    TN <- confusion$table[1, 1]
    AUC <- ifelse(length(unique(probs)) > 1, auc(as.numeric(test_data$Label), probs), NA)
    MCC <- ifelse((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) > 0, 
                  ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)), NA)
    metrics[i, ] <- c(MCC, confusion$byClass["F1"], AUC, confusion$byClass["Sensitivity"], 
                      confusion$byClass["Specificity"], confusion$byClass["Pos Pred Value"], 
                      confusion$byClass["Neg Pred Value"])
  }
  
  avg_metrics <- colMeans(metrics, na.rm = TRUE)
  feature_removal_results <- rbind(feature_removal_results, 
                                   data.frame(Feature = feature_to_remove, 
                                              MCC = avg_metrics[1],
                                              F1 = avg_metrics[2],
                                              AUC = avg_metrics[3],
                                              TPR = avg_metrics[4],
                                              TNR = avg_metrics[5],
                                              PPV = avg_metrics[6],
                                              NPV = avg_metrics[7]))
}

# Save feature removal results
write.csv(feature_removal_results, "feature_removal_results100159.csv", row.names = FALSE)

# Plot MCC impact from feature removal
ggplot(feature_removal_results, aes(x = reorder(Feature, -MCC), y = MCC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Impact of Feature Removal on MCC", x = "Feature", y = "MCC") +
  theme_minimal() +
  theme(axis.text = element_text(size = 7), axis.title = element_text(size = 14))

ggsave("feature_removal_mcc_plot100159.png", width = 10, height = 10, dpi = 300)

################################################################################
# Sanity Check: Adding Noise to Data
################################################################################
noise_levels <- seq(0, 50, by = 5)  # Noise percentages in smaller increments
n_iterations <- 100  # Number of iterations to average results per noise level

sanity_check_results <- data.frame(
  Noise_Level = integer(),
  MCC = numeric(),
  F1 = numeric(),
  AUC = numeric(),
  TPR = numeric(),
  TNR = numeric(),
  PPV = numeric(),
  NPV = numeric()
)

for (noise_level in noise_levels) {
  cat("Adding noise at level:", noise_level, "%\n")
  
  # Initialize storage for metrics at the current noise level
  metrics_list <- data.frame(
    MCC = numeric(n_iterations),
    F1 = numeric(n_iterations),
    AUC = numeric(n_iterations),
    TPR = numeric(n_iterations),
    TNR = numeric(n_iterations),
    PPV = numeric(n_iterations),
    NPV = numeric(n_iterations)
  )
  
  for (iteration in 1:n_iterations) {
    set.seed(123 + iteration)  # Ensure reproducibility with variation across iterations
    
    # Add normalized noise to the dataset
    noisy_data <- data
    if (noise_level > 0) {
      num_noisy_rows <- floor(nrow(data) * noise_level / 100)
      noisy_rows <- sample(1:nrow(data), num_noisy_rows)
      
      # Normalize features before adding noise
      noisy_features <- noisy_data[noisy_rows, -ncol(noisy_data)]
      noisy_features <- scale(noisy_features) + 
        rnorm(n = num_noisy_rows * (ncol(data) - 1), mean = 0, sd = 1)
      
      # Denormalize features back to original scale
      noisy_features <- noisy_features * apply(data[, -ncol(data)], 2, sd) +
        apply(data[, -ncol(data)], 2, mean)
      noisy_data[noisy_rows, -ncol(noisy_data)] <- noisy_features
    }
    
    # Stratified train-test split
    train_index <- createDataPartition(noisy_data$Label, p = 0.8, list = FALSE)
    train_data <- noisy_data[train_index, ]
    test_data <- noisy_data[-train_index, ]
    
    # Apply SMOTE to balance the training dataset
    smote_recipe <- recipe(Label ~ ., data = train_data) %>%
      step_smote(Label) %>%
      prep()
    train_data_smote <- bake(smote_recipe, new_data = NULL)
    
    # Train Random Forest model on SMOTE-balanced dataset
    rf_model <- randomForest(Label ~ ., data = train_data_smote, ntree = 100, mtry = sqrt(ncol(train_data_smote) - 1))
    
    # Predict on test data
    probs <- predict(rf_model, newdata = test_data, type = "prob")[, 2]
    predictions <- ifelse(probs > 0.5, "Sepsis", "Control")
    
    # Skip iteration if only one class present in test data or predictions
    if (length(unique(test_data$Label)) < 2 || length(unique(predictions)) < 2) {
      next
    }
    
    # Calculate metrics
    confusion <- confusionMatrix(factor(predictions, levels = levels(test_data$Label)), 
                                 test_data$Label, positive = "Sepsis")
    TP <- confusion$table[2, 2]
    FP <- confusion$table[1, 2]
    FN <- confusion$table[2, 1]
    TN <- confusion$table[1, 1]
    AUC <- ifelse(length(unique(probs)) > 1, auc(as.numeric(test_data$Label), probs), NA)
    MCC <- ifelse((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) > 0, 
                  ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)), NA)
    metrics_list[iteration, ] <- c(
      MCC, confusion$byClass["F1"], AUC, confusion$byClass["Sensitivity"],
      confusion$byClass["Specificity"], confusion$byClass["Pos Pred Value"], 
      confusion$byClass["Neg Pred Value"]
    )
  }
  
  # Average metrics across iterations for the current noise level
  avg_metrics <- colMeans(metrics_list, na.rm = TRUE)
  sanity_check_results <- rbind(sanity_check_results, 
                                data.frame(Noise_Level = noise_level, 
                                           MCC = avg_metrics[1],
                                           F1 = avg_metrics[2],
                                           AUC = avg_metrics[3],
                                           TPR = avg_metrics[4],
                                           TNR = avg_metrics[5],
                                           PPV = avg_metrics[6],
                                           NPV = avg_metrics[7]))
}

# Save sanity check results
write.csv(sanity_check_results, "sanity_check_results100159.csv", row.names = FALSE)

# Plot impact of noise on metrics
noise_plot <- ggplot(sanity_check_results, aes(x = Noise_Level)) +
  geom_line(aes(y = MCC, color = "MCC"), size = 1.2) +
  geom_line(aes(y = F1, color = "F1 Score"), size = 1.2) +
  geom_line(aes(y = AUC, color = "AUC"), size = 1.2) +
  labs(title = "Impact of Noise on Model Performance", x = "Noise Level (%)", y = "Metric Value", color = "Metrics") +
  theme_minimal() +
  scale_color_manual(values = c("MCC" = "blue", "F1 Score" = "green", "AUC" = "red"))

print(noise_plot)
ggsave("impact_of_noise_on_model_performance100159.png", plot = noise_plot, width = 10, height = 6, dpi = 300)


# Stop parallel cluster
stopCluster(cl)
