
################################################################################
############################# SMOTE RESAMPLING #################################
################################################################################

################################################################################
# Improved Code: Random Forest with SMOTE, Noise Analysis, Feature Removal, and Robust Evaluation
################################################################################

###############################################################################
# Load and Prepare Data
###############################################################################
# Load necessary libraries
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(randomForest, caret, pROC, dplyr, ggplot2, parallel, doParallel, smotefamily, ROSE, viridis)

# Set up parallel computing
num_cores <- detectCores() - 1  # Leave one core free
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Read in your data
data <- read.csv("sepsis_allD2_GSE137342_10558.csv")

# Remove the Sample column and ensure Label is a factor
data <- dplyr::select(data, -Sample)
data <- data %>% mutate(Label = as.factor(Label))

# Precompute partitions for all repeats
n_repeats <- 100
all_partitions <- vector("list", n_repeats)
for (i in 1:n_repeats) {
  all_partitions[[i]] <- createDataPartition(data$Label, p = 0.8, list = FALSE)
}

# Storage for feature importance
all_feature_importance <- data.frame()
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
  
  # Split into train and test
  train_data <- data[train_index, ]
  test_data  <- data[-train_index, ]
  
  
  # Calculate the SMOTE oversampling ratio (avoid negative, force integer)
  dup_size_value <- max(
    as.integer(ceiling(
      (sum(train_data$Label == "Sepsis") / sum(train_data$Label == "Healthy")) - 1
    )),
    0
  )
  
  # Apply SMOTE (from smotefamily)
  smote_output <- SMOTE(
    X        = train_data[, -ncol(train_data)],  # predictors only
    target   = train_data$Label,                 # factor target
    K        = 5,
    dup_size = dup_size_value
  )
  
  smote_train_data <- smote_output$data
  
  # Rename the last column to "Label" and ensure factor levels match
  colnames(smote_train_data)[ncol(smote_train_data)] <- "Label"
  smote_train_data$Label <- factor(smote_train_data$Label, levels = levels(train_data$Label))
  
  # Verify the class distribution
  print(table(train_data$Label))
  print(table(smote_train_data$Label))
  
  # Train Random Forest on SMOTE-adjusted training set
  rf_model <- randomForest(Label ~ ., 
                           data  = smote_train_data,
                           ntree = 100,
                           mtry  = floor(sqrt(ncol(smote_train_data) - 1)) # typical heuristic
  )
  
  # -----------------------------
  # 🔹 Extract Feature Importance
  # -----------------------------
  feature_importance <- as.data.frame(importance(rf_model))
  feature_importance$Feature <- rownames(feature_importance)
  feature_importance$Iteration <- i  # Track which iteration this belongs to
  
  # Store feature importance for later analysis
  all_feature_importance <- rbind(all_feature_importance, feature_importance)
  
  # Predict probabilities and classes on the test set
  probs <- predict(rf_model, newdata = test_data, type = "prob")[, "Sepsis"]
  predictions <- ifelse(probs > 0.5, "Sepsis", "Healthy")
  
  # Skip iteration if classes are missing in either the test set or predictions
  if (length(unique(test_data$Label)) < 2 || length(unique(predictions)) < 2) {
    cat("Iteration", i, "skipped due to class imbalance in test set or predictions.\n")
    next
  }
  
  # Confusion matrix (caret)
  confusion <- confusionMatrix(
    data      = factor(predictions, levels = levels(test_data$Label)),
    reference = test_data$Label, 
    positive  = "Sepsis"
  )
  
  # Extract confusion matrix terms
  TP <- confusion$table[2, 2]
  FP <- confusion$table[1, 2]
  FN <- confusion$table[2, 1]
  TN <- confusion$table[1, 1]
  
  # Compute AUC (only if probs have > 1 unique value)
  AUC_val <- ifelse(length(unique(probs)) > 1,
                    pROC::auc(as.numeric(test_data$Label), probs), 
                    NA)
  
  # Compute MCC (only if denominator is nonzero)
  denom <- (TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)
  MCC_val <- ifelse(denom > 0,
                    ((TP * TN) - (FP * FN)) / sqrt(denom),
                    NA)
  
  # Save results
  results[i, ] <- c(
    MCC_val,
    confusion$byClass["F1"],
    AUC_val,
    confusion$byClass["Sensitivity"],  # TPR
    confusion$byClass["Specificity"],  # TNR
    confusion$byClass["Pos Pred Value"],
    confusion$byClass["Neg Pred Value"]
  )
}



###############################################################################
# Save and Summarize Results
###############################################################################
# Write the per-iteration metrics
write.csv(results, "repeated_splits_SMOTE_allD2_GSE137342.csv", row.names = FALSE)

# Average baseline metrics
average_metrics <- colMeans(results, na.rm = TRUE)
cat("Baseline (All Features) Average Metrics with SMOTE:\n")
print(average_metrics)

# Print each metric
for (i in seq_along(average_metrics)) {
  cat(names(average_metrics)[i], ":", average_metrics[i], "\n")
}

# Create and save a data frame for the average metrics
average_metrics_df <- data.frame(
  Metric = names(average_metrics),
  Value  = average_metrics
)
write.csv(average_metrics_df, "average_metrics_SMOTE_allD2_GSE137342.csv", row.names = FALSE)
cat("Average metrics saved to 'average_metrics_SMOTE_allD2_GSE137342.csv'\n")

# Save feature importance data
write.csv(all_feature_importance, "feature_importance_SMOTE_allD2_GSE137342.csv", row.names = FALSE)
cat("Feature importance results saved to 'feature_importance_SMOTE_allD2_GSE137342.csv'\n")

# Compute average importance across all iterations
avg_importance <- all_feature_importance %>%
  group_by(Feature) %>%
  summarize(Mean_Importance = mean(MeanDecreaseGini, na.rm = TRUE)) %>%
  arrange(desc(Mean_Importance))

# Save average feature importance
write.csv(avg_importance, "average_feature_importance_SMOTE_allD2_GSE137342.csv", row.names = FALSE)
cat("Average feature importance saved to 'average_feature_importance_SMOTE_allD2_GSE137342.csv'\n")

# ----------------------------------------
# 📊 Enhanced Feature Importance Plot
# ----------------------------------------

# Select top 20 most important features for better visualization
top_n_features <- 50
top_features <- avg_importance %>% top_n(top_n_features, Mean_Importance)

# Create the plot
feature_plot <- ggplot(top_features, aes(x = reorder(Feature, Mean_Importance), y = Mean_Importance, fill = Mean_Importance)) +
  geom_bar(stat = "identity", width = 0.7, color = "black") +  # Add black borders for better contrast
  coord_flip() +  # Flip axes for readability
  scale_fill_viridis(option = "magma", direction = -1, begin = 0.2, end = 0.9) +  # Better color scale
  theme_light(base_size = 14) +  # Modern, clean theme
  labs(
    title = paste("Top", top_n_features, "Most Important Features (Random Forest + SMOTE)"),
    subtitle = "Feature importance based on Mean Decrease in Gini",
    x = "Feature",
    y = "Importance Score",
    fill = "Mean Decrease in Gini"
  ) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Center title
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),  # Subtitle styling
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  ) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 10))  # Improve color legend layout

# Save and display the plot
print(feature_plot)
ggsave("top_feature_importance_plot_SMOTE_allD2_GSE137342.png", plot = feature_plot, width = 12, height = 7, dpi = 300)


###############################################################################
# 3) Feature Removal Analysis, re-using the SAME partitions
###############################################################################
###############################################################################
# Feature Removal with SMOTE
###############################################################################
#   1) 'data' loaded with columns for features + "Label"
#   2) 'all_partitions' already created (100 repeated splits)

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
  
  # Remove the one specified feature
  modified_data <- data[, !names(data) %in% feature_to_remove]
  
  # Prepare a container for repeated-split metrics
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
    
    # Split into train and test
    train_data <- modified_data[train_index, ]
    test_data  <- modified_data[-train_index, ]
    
    # SMOTE Oversampling
    # Ensure sepsis is indeed the minority class and duplication is integer
    dup_size_value <- max(
      as.integer(ceiling(
        (sum(train_data$Label == "Healthy") / sum(train_data$Label == "Sepsis")) - 1
      )),
      0
    )
    
    # Apply SMOTE (from smotefamily)
    smote_output <- SMOTE(
      X        = train_data[, -ncol(train_data)],  # all predictors
      target   = train_data$Label,                 # factor target
      K        = 5,
      dup_size = dup_size_value
    )
    
    smote_train_data <- smote_output$data
    
    # Rename last column to "Label" and ensure factor levels
    colnames(smote_train_data)[ncol(smote_train_data)] <- "Label"
    smote_train_data$Label <- factor(smote_train_data$Label,
                                     levels = levels(train_data$Label))
    
    # Train a Random Forest using the SMOTE-augmented training data
    rf_model <- randomForest(
      Label ~ ., 
      data  = smote_train_data,
      ntree = 100,
      mtry  = floor(sqrt(ncol(smote_train_data) - 1))  # heuristic
    )
    
    # Predict probabilities and classes on the test data
    probs <- predict(rf_model, newdata = test_data, type = "prob")[, "Sepsis"]
    predictions <- ifelse(probs > 0.5, "Sepsis", "Healthy")
    
    # If either the test set or predictions have only one class, skip
    if (length(unique(test_data$Label)) < 2 || length(unique(predictions)) < 2) {
      next
    }
    
    # Evaluate performance
    confusion <- confusionMatrix(
      data      = factor(predictions, levels = levels(test_data$Label)),
      reference = test_data$Label, 
      positive  = "Sepsis"
    )
    
    # Extract confusion matrix elements
    TP <- confusion$table[2, 2]
    FP <- confusion$table[1, 2]
    FN <- confusion$table[2, 1]
    TN <- confusion$table[1, 1]
    
    # AUC (only if probs have > 1 unique value)
    AUC_val <- ifelse(length(unique(probs)) > 1,
                      pROC::auc(as.numeric(test_data$Label), probs),
                      NA)
    
    # MCC (only if denominator > 0)
    denom <- (TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)
    MCC_val <- ifelse(denom > 0,
                      ((TP * TN) - (FP * FN)) / sqrt(denom),
                      NA)
    
    metrics[i, ] <- c(
      MCC_val,
      confusion$byClass["F1"],
      AUC_val,
      confusion$byClass["Sensitivity"],
      confusion$byClass["Specificity"],
      confusion$byClass["Pos Pred Value"],
      confusion$byClass["Neg Pred Value"]
    )
  }
  
  # Compute average metrics across repeats for this feature
  avg_metrics <- colMeans(metrics, na.rm = TRUE)
  
  feature_removal_results <- rbind(
    feature_removal_results,
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
    )
  )
}

# Save results to CSV
write.csv(feature_removal_results, "feature_removal_results_SMOTE_allD2_GSE137342.csv", row.names = FALSE)

###############################################################################
# 4) Compare "All Features" vs. "Removed One Feature" Using the Same Splits
###############################################################################
# we have a variable 'average_metrics' from our baseline (all features) scenario:
baseline_mcc <- average_metrics["MCC"]
cat("\nBaseline MCC (using all features):", baseline_mcc, "\n\n")

# 2) Combine baseline and feature-removal MCCs into one data frame
mcc_comparison_df <- rbind(
  data.frame(Feature = "All Features (Baseline)", MCC = baseline_mcc),
  feature_removal_results[, c("Feature", "MCC")]
)

# 3) Reorder factors so that the largest MCC is displayed at the bottom of the plot
mcc_comparison_df$Feature <- reorder(mcc_comparison_df$Feature, -mcc_comparison_df$MCC)

# 4) Create a bar plot of MCC for baseline vs. each feature removed
mcc_plot <- ggplot(mcc_comparison_df, aes(x = Feature, y = MCC, fill = MCC)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +  # Wider bars for clarity
  coord_flip() +  # Flip axes for readability
  scale_fill_gradient(low = "#E1735A", high = "#138A78") +  # Better color contrast
  theme_light(base_size = 14) +  # Clean theme with readable text
  labs(
    title = "Impact of Feature Removal on MCC Score",
    subtitle = "Baseline vs. Individual Feature Removal",
    x = "Feature Removed",
    y = "MCC Score",
    fill = "MCC"
  ) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),  # Center title
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),  # Subtitle styling
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 7),
    legend.position = "right",
  ) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 10))  # Improve legend layout


ggsave("feature_removal_mcc_plot_SMOTE_allD2_GSE137342.png", width = 10, height = 12, dpi = 300)


###############################################################################
# 4) Sanity Check: Adding Noise + Using Same Splits (with SMOTE Fix)
###############################################################################

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
    
    # Simple approach: multiply numeric features by ~1 ± small factor
    numeric_cols <- which(sapply(noisy_data, is.numeric))
    for (r in noisy_rows) {
      for (col in numeric_cols) {
        noisy_data[r, col] <- noisy_data[r, col] * runif(1, min = 0.8, max = 1.2)
      }
    }
  }
  
  # We'll run the same n_repeats using the same partitions
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
    
    # SMOTE on the noisy training (ensure integer oversampling)
    dup_size_value <- max(
      as.integer(ceiling(
        (sum(train_data$Label == "healthy") / sum(train_data$Label == "Sepsis")) - 1
      )),
      0
    )
    
    smote_output <- SMOTE(
      X        = train_data[, -ncol(train_data)],  # all predictors except Label
      target   = train_data$Label,                 # factor target
      K        = 5,
      dup_size = dup_size_value
    )
    
    smote_train_data <- smote_output$data
    
    # Rename the last column to "Label" and align factor levels
    colnames(smote_train_data)[ncol(smote_train_data)] <- "Label"
    smote_train_data$Label <- factor(smote_train_data$Label,
                                     levels = levels(train_data$Label))
    
    # Train RF on SMOTE-adjusted data
    rf_model <- randomForest(
      Label ~ ., 
      data  = smote_train_data,
      ntree = 100,
      mtry  = floor(sqrt(ncol(smote_train_data) - 1))
    )
    
    # Predict and evaluate
    probs <- predict(rf_model, newdata = test_data, type = "prob")[, "Sepsis"]
    predictions <- ifelse(probs > 0.5, "Sepsis", "Healthy")
    
    # Skip if any class is missing in test labels or predictions
    if (length(unique(test_data$Label)) < 2 || length(unique(predictions)) < 2) {
      next
    }
    
    confusion <- confusionMatrix(
      data      = factor(predictions, levels = levels(test_data$Label)),
      reference = test_data$Label, 
      positive  = "Sepsis"
    )
    
    TP <- confusion$table[2, 2]
    FP <- confusion$table[1, 2]
    FN <- confusion$table[2, 1]
    TN <- confusion$table[1, 1]
    
    # AUC (using pROC)
    AUC_val <- ifelse(
      length(unique(probs)) > 1,
      pROC::auc(as.numeric(test_data$Label), probs),
      NA
    )
    
    # MCC
    denom <- (TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)
    MCC_val <- ifelse(
      denom > 0,
      ((TP * TN) - (FP * FN)) / sqrt(denom),
      NA
    )
    
    # Store metrics for this iteration
    metrics[i, ] <- c(
      MCC_val,
      confusion$byClass["F1"],
      AUC_val,
      confusion$byClass["Sensitivity"],
      confusion$byClass["Specificity"],
      confusion$byClass["Pos Pred Value"],
      confusion$byClass["Neg Pred Value"]
    )
  }
  
  # Average metrics across repeats
  avg_metrics <- colMeans(metrics, na.rm = TRUE)
  sanity_check_results <- rbind(
    sanity_check_results,
    data.frame(
      Noise_Level = noise_level,
      MCC         = avg_metrics[1],
      F1          = avg_metrics[2],
      AUC         = avg_metrics[3],
      TPR         = avg_metrics[4],
      TNR         = avg_metrics[5],
      PPV         = avg_metrics[6],
      NPV         = avg_metrics[7]
    )
  )
}

# Save results
write.csv(sanity_check_results, "sanity_check_results_SMOTE_allD2_GSE137342.csv", row.names = FALSE)

# Plot impact of noise on metrics
noise_plot <- ggplot(sanity_check_results, aes(x = Noise_Level)) +
  # Add line plots with distinct markers for better visibility
  geom_line(aes(y = MCC, color = "MCC"), size = 1.2) +
  geom_point(aes(y = MCC, color = "MCC"), size = 3, shape = 16) + # Add circular points
  
  geom_line(aes(y = F1, color = "F1"), size = 1.2) +
  geom_point(aes(y = F1, color = "F1"), size = 3, shape = 17) + # Add triangular points
  
  geom_line(aes(y = AUC, color = "AUC"), size = 1.2) +
  geom_point(aes(y = AUC, color = "AUC"), size = 3, shape = 15) + # Add square points
  
  # Customize labels and titles
  labs(
    title = "Impact of Noise on Model Performance",
    subtitle = "Model robustness evaluated at different noise levels",
    x = "Noise Level (%)",
    y = "Performance Metric",
    color = "Metrics"
  ) +
  
  # Improve theme for clarity
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold")
  ) +
  
  # Use a better color scheme (Blue-Green-Yellow)
  scale_color_manual(values = c("MCC" = "#098689",  # Blue
                                "F1" = "#F75F5C",   # Green
                                "AUC" = "#E69F00")) # Orange

# Display & Save the Plot
print(noise_plot)
ggsave("noise_plot_allD2_GSE137342.png", plot = noise_plot, width = 10, height = 6, dpi = 300)

# Stop parallel cluster
stopCluster(cl)

################################################################################
# End
################################################################################
