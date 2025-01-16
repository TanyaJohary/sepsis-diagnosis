# Load necessary libraries
library(randomForest)
library(caret)
library(pROC)  
library(dplyr)

# Set seed for reproducibility
set.seed(123)

# Load and prepare data
data <- read.csv("T2_sepsis_labeledGSE69063.csv")
data <- dplyr::select(data, -Sample)
data <- data %>% mutate(Label = as.factor(Label))  # Ensure Label is a factor
data$Label <- factor(data$Label, levels = c("healthy control", "sepsis"))  # Explicitly set levels

# Define the number of repetitions
n_repeats <- 100

# Initialize a data frame to store metrics for each split
results <- data.frame(
  MCC = numeric(n_repeats),
  F1 = numeric(n_repeats),
  AUC = numeric(n_repeats),
  TPR = numeric(n_repeats),
  TNR = numeric(n_repeats),
  PPV = numeric(n_repeats),
  NPV = numeric(n_repeats)
)

# Perform 100 random train-test splits
for (i in 1:n_repeats) {
  # Random 80/20 split
  train_index <- createDataPartition(data$Label, p = 0.8, list = FALSE)
  train_data <- data[train_index, ]
  test_data <- data[-train_index, ]
  
  # Ensure Label levels are in the correct order
  test_data$Label <- factor(test_data$Label, levels = c("healthy control", "sepsis"))
  
  # Train Random Forest model
  rf_model <- randomForest(Label ~ ., data = train_data, ntree = 100, mtry = sqrt(ncol(train_data) - 1))
  
  # Predict on test data
  predictions <- predict(rf_model, newdata = test_data)
  predictions <- factor(predictions, levels = levels(test_data$Label))  # Match levels with test_data$Label
  
  # Predict probabilities for AUC
  probs <- predict(rf_model, newdata = test_data, type = "prob")[, 2]  # Probabilities for "sepsis"
  
  # Confusion matrix
  confusion <- confusionMatrix(predictions, test_data$Label, positive = "sepsis")
  TP <- confusion$table[2, 2]
  FP <- confusion$table[1, 2]
  FN <- confusion$table[2, 1]
  TN <- confusion$table[1, 1]
  
  # Calculate metrics
  MCC <- ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  F1 <- confusion$byClass["F1"]
  
  # Check AUC validity
  if (length(unique(probs)) > 1) {
    AUC <- auc(response = test_data$Label, predictor = probs, levels = c("healthy control", "sepsis"), direction = "<")
  } else {
    AUC <- NA  # Set to NA if probs lacks variability
  }
  
  # Debugging: Print current iteration details
  cat("Iteration:", i, "AUC:", AUC, "Labels:", test_data$Label, "Probs:", probs, "\n")
  
  TPR <- confusion$byClass["Sensitivity"]  # True Positive Rate
  TNR <- confusion$byClass["Specificity"]  # True Negative Rate
  PPV <- confusion$byClass["Pos Pred Value"]  # Positive Predictive Value
  NPV <- confusion$byClass["Neg Pred Value"]  # Negative Predictive Value
  
  # Store metrics in the results data frame
  results[i, ] <- c(MCC, F1, AUC, TPR, TNR, PPV, NPV)
}

# Calculate average metrics over all splits
average_metrics <- colMeans(results, na.rm = TRUE)

# Save detailed results to a CSV file
write.csv(results, "repeated_splits_metrics.csv", row.names = FALSE)

# Save average metrics to a text file
sink("average_metrics.txt")
cat("Average Metrics Over 100 Splits:\n")
print(average_metrics)
sink()

# Print average metrics to console
print(average_metrics)


#################################################################
#################leave-one-feature-out testing###################
#################################################################


# Load necessary libraries
library(randomForest)
library(caret)
library(pROC)  
library(dplyr)
library(ggplot2)

# Set seed for reproducibility
set.seed(123)

# Load and prepare our dataset
data <- read.csv("T2_sepsis_labeledGSE69063.csv")
data <- dplyr::select(data, -Sample)
data <- data %>% mutate(Label = as.factor(Label))  # Ensure Label is a factor

# Explicitly set "Sepsis" as the positive class
data$Label <- factor(data$Label, levels = c("sepsis", "healthy control"))

# Verify levels
print(levels(data$Label))  

# Get feature names (excluding the target column)
features <- setdiff(names(data), "Label")

# Initialize a data frame to store average metrics for each removed feature
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

# Loop over each feature to remove it
for (feature_to_remove in features) {
  cat("Removing feature:", feature_to_remove, "\n")  # Progress indicator
  
  # Remove the feature from the dataset
  modified_data <- data %>% select(-all_of(feature_to_remove))
  
  # Initialize metrics storage for this feature removal
  metrics <- data.frame(
    MCC = numeric(100),
    F1 = numeric(100),
    AUC = numeric(100),
    TPR = numeric(100),
    TNR = numeric(100),
    PPV = numeric(100),
    NPV = numeric(100)
  )
  
  # Perform 100 train/test splits
  for (i in 1:100) {
    cat("Iteration:", i, "for feature:", feature_to_remove, "\n")  # Progress indicator
    
    # Random 80/20 split
    train_index <- createDataPartition(modified_data$Label, p = 0.8, list = FALSE)
    train_data <- modified_data[train_index, ]
    test_data <- modified_data[-train_index, ]
    
    # Ensure test_data levels match original data levels
    test_data$Label <- factor(test_data$Label, levels = levels(data$Label))
    
    # Train Random Forest model
    rf_model <- randomForest(Label ~ ., data = train_data, ntree = 100, mtry = sqrt(ncol(train_data) - 1))  # Reduced ntree for speed
    
    # Predict on test data
    predictions <- predict(rf_model, newdata = test_data)
    probs <- predict(rf_model, newdata = test_data, type = "prob")[, 2]  # Sepsis probabilities
    
    # Confusion matrix: specify positive class as "Sepsis"
    confusion <- confusionMatrix(predictions, test_data$Label, positive = "sepsis")
    
    # AUC calculation: explicitly set levels and direction
    if (length(unique(probs)) > 1) {
      roc_obj <- roc(test_data$Label, probs, levels = c("sepsis", "healthy control"), direction = "<")
      auc_val <- auc(roc_obj)
    } else {
      auc_val <- NA  # Handle cases where AUC cannot be calculated
    }
    
    # Extract metrics from confusion matrix and calculations
    TP <- confusion$table[2, 2]
    FP <- confusion$table[1, 2]
    FN <- confusion$table[2, 1]
    TN <- confusion$table[1, 1]
    
    MCC <- ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    F1 <- confusion$byClass["F1"]
    TPR <- confusion$byClass["Sensitivity"]  # True Positive Rate
    TNR <- confusion$byClass["Specificity"]  # True Negative Rate
    PPV <- confusion$byClass["Pos Pred Value"]  # Positive Predictive Value
    NPV <- confusion$byClass["Neg Pred Value"]  # Negative Predictive Value
    
    # Store metrics for this iteration
    metrics[i, ] <- c(MCC, F1, auc_val, TPR, TNR, PPV, NPV)
  }
  
  # Calculate average metrics for this feature removal
  avg_metrics <- colMeans(metrics, na.rm = TRUE)
  
  # Store results
  feature_removal_results <- rbind(
    feature_removal_results,
    data.frame(
      Feature = feature_to_remove,  # Ensure the feature name is correctly recorded
      MCC = avg_metrics[1],
      F1 = avg_metrics[2],
      AUC = avg_metrics[3],
      TPR = avg_metrics[4],
      TNR = avg_metrics[5],
      PPV = avg_metrics[6],
      NPV = avg_metrics[7],
      stringsAsFactors = FALSE
    )
  )
  
  # Save intermediate results after each feature
  write.csv(feature_removal_results, "feature_removal_results.csv", row.names = FALSE)
}

# Verify the structure of the final results
str(feature_removal_results)

# Plot
ggplot(feature_removal_results, aes(x = reorder(Feature, -MCC), y = MCC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Impact of Gene Removal on MCC", x = "Gene", y = "MCC") +
  theme_minimal()
