###############################################################################
# Load Libraries
###############################################################################
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  randomForest, caret, pROC, dplyr, ggplot2,
  parallel, doParallel, smotefamily, ROSE, viridis, foreach
)

###############################################################################
# 0) Parallel Setup & Reproducibility
###############################################################################
num_cores <- detectCores() - 1  # leave one core free
cl <- makeCluster(num_cores)
registerDoParallel(cl)

#set.seed(123)  # partial reproducibility

###############################################################################
# 1) Load Data and Create Stratified Partitions
###############################################################################
data <- read.csv("cluster3_T0_GSE69063.csv")

# Remove 'Sample' column and ensure 'Label' is factor
data <- dplyr::select(data, -Sample)
data$Label <- factor(data$Label, levels = c("healthy control","sepsis"))  # consistent factor levels

n_repeats <- 100

# Create repeated stratified partitions (80% train, 20% test)
all_partitions <- vector("list", n_repeats)
for (i in seq_len(n_repeats)) {
  all_partitions[[i]] <- createDataPartition(data$Label, p = 0.8, list = FALSE)
}

###############################################################################
# 2) Baseline Random Forest + SMOTE
#    (Store iteration messages & the full confusionMatrix for iteration #1)
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

# Lists for confusion & distributions
confusion_list  <- vector("list", n_repeats)
train_dist_list <- vector("list", n_repeats)
test_dist_list  <- vector("list", n_repeats)
smote_dist_list <- vector("list", n_repeats)

# For feature importance across iterations
all_feature_importance <- data.frame()

# Track how many iterations were skipped
skipped_count <- 0

# We'll store the iteration messages (printed after the loop)
iteration_messages <- character(n_repeats)

# We'll also store the entire confusionMatrix for iteration #1
sample_conf_mat <- NULL

###############################################################################
# 2a) Parallel loop
###############################################################################
results_list <- foreach(i = seq_len(n_repeats), 
                        .packages = c("randomForest","smotefamily","caret","pROC","dplyr")) %dopar% {
                          # Build a message for this iteration (not printed yet)
                          iteration_msg <- paste0(
                            "Iteration ", i, ":\n",
                            "Setting levels: control = 1, case = 2\n",
                            "Setting direction: controls < cases\n\n"
                          )
                          
                          train_index <- all_partitions[[i]]
                          train_data  <- data[train_index, ]
                          test_data   <- data[-train_index, ]
                          
                          # Factor levels
                          train_data$Label <- factor(train_data$Label, levels = c("healthy control","sepsis"))
                          test_data$Label  <- factor(test_data$Label, levels = c("healthy control","sepsis"))
                          
                          # Capture distributions
                          train_dist <- table(train_data$Label)
                          test_dist  <- table(test_data$Label)
                          
                          # If training set has only one class => skip
                          if (length(unique(train_data$Label)) < 2) {
                            return(list(
                              metrics            = c(NA, NA, NA, NA, NA, NA, NA),
                              feature_importance = NULL,
                              skipped            = TRUE,
                              confusion          = NULL,
                              train_dist         = train_dist,
                              test_dist          = test_dist,
                              smote_dist         = NULL,
                              iteration_message  = iteration_msg,
                              full_conf_obj      = NULL
                            ))
                          }
                          
                          n_Sepsis  <- sum(train_data$Label == "sepsis")
                          n_Healthy <- sum(train_data$Label == "healthy control")
                          
                          if (n_Sepsis == 0 || n_Healthy == 0) {
                            return(list(
                              metrics            = c(NA, NA, NA, NA, NA, NA, NA),
                              feature_importance = NULL,
                              skipped            = TRUE,
                              confusion          = NULL,
                              train_dist         = train_dist,
                              test_dist          = test_dist,
                              smote_dist         = NULL,
                              iteration_message  = iteration_msg,
                              full_conf_obj      = NULL
                            ))
                          }
                          
                          # Compute SMOTE ratio
                          dup_size_value <- max(
                            as.integer(ceiling((n_Healthy / n_Sepsis) - 1)),
                            0
                          )
                          
                          # SMOTE
                          smote_output <- SMOTE(
                            X        = train_data[, -ncol(train_data)],
                            target   = train_data$Label,
                            K        = 5,
                            dup_size = dup_size_value
                          )
                          
                          smote_train_data <- smote_output$data
                          colnames(smote_train_data)[ncol(smote_train_data)] <- "Label"
                          smote_train_data$Label <- factor(smote_train_data$Label, levels = c("healthy control","sepsis"))
                          
                          # Capture SMOTE distribution
                          smote_dist <- table(smote_train_data$Label)
                          
                          # Train Random Forest
                          rf_model <- randomForest(
                            Label ~ .,
                            data  = smote_train_data,
                            ntree = 100,
                            mtry  = floor(sqrt(ncol(smote_train_data) - 1))
                          )
                          
                          # Feature importance
                          feature_imp <- as.data.frame(importance(rf_model))
                          feature_imp$Feature <- rownames(feature_imp)
                          feature_imp$Iteration <- i
                          
                          # Predict
                          probs <- predict(rf_model, newdata = test_data, type = "prob")[, "sepsis"]
                          preds <- ifelse(probs > 0.5, "sepsis", "healthy control")
                          
                          if (length(unique(test_data$Label)) < 2 || length(unique(preds)) < 2) {
                            return(list(
                              metrics            = c(NA, NA, NA, NA, NA, NA, NA),
                              feature_importance = feature_imp,
                              skipped            = TRUE,
                              confusion          = NULL,
                              train_dist         = train_dist,
                              test_dist          = test_dist,
                              smote_dist         = smote_dist,
                              iteration_message  = iteration_msg,
                              full_conf_obj      = NULL
                            ))
                          }
                          
                          # Build confusionMatrix
                          confusion <- confusionMatrix(
                            data      = factor(preds, levels = c("healthy control","sepsis")),
                            reference = test_data$Label,
                            positive  = "sepsis"
                          )
                          cm_table <- confusion$table
                          
                          # We only store the entire confusionMatrix for iteration #1
                          if (i == 1) {
                            full_conf_obj <- confusion
                          } else {
                            full_conf_obj <- NULL
                          }
                          
                          # Compute metrics
                          TP <- cm_table[2,2]
                          FP <- cm_table[1,2]
                          FN <- cm_table[2,1]
                          TN <- cm_table[1,1]
                          
                          AUC_val <- if (length(unique(probs)) > 1) {
                            pROC::auc(as.numeric(test_data$Label), probs)
                          } else {
                            NA
                          }
                          
                          denom <- (TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)
                          MCC_val <- if (denom > 0) {
                            ((TP * TN) - (FP * FN)) / sqrt(denom)
                          } else {
                            NA
                          }
                          
                          metrics_vector <- c(
                            MCC_val,
                            confusion$byClass["F1"],
                            AUC_val,
                            confusion$byClass["Sensitivity"],
                            confusion$byClass["Specificity"],
                            confusion$byClass["Pos Pred Value"],
                            confusion$byClass["Neg Pred Value"]
                          )
                          
                          list(
                            metrics            = metrics_vector,
                            feature_importance = feature_imp,
                            skipped            = FALSE,
                            confusion          = cm_table,
                            train_dist         = train_dist,
                            test_dist          = test_dist,
                            smote_dist         = smote_dist,
                            iteration_message  = iteration_msg,
                            full_conf_obj      = full_conf_obj
                          )
                        }

###############################################################################
# 2b) Post-processing the parallel results
###############################################################################
for (i in seq_along(results_list)) {
  out <- results_list[[i]]
  
  # Store iteration messages
  iteration_messages[i] <- out$iteration_message
  
  if (out$skipped) {
    skipped_count <- skipped_count + 1
    results[i, ] <- rep(NA, 7)
  } else {
    results[i, ] <- out$metrics
    all_feature_importance <- rbind(all_feature_importance, out$feature_importance)
  }
  
  confusion_list[[i]]  <- out$confusion
  train_dist_list[[i]] <- out$train_dist
  test_dist_list[[i]]  <- out$test_dist
  smote_dist_list[[i]] <- out$smote_dist
  
  # If iteration #1 wasn't skipped, we store the entire confusionMatrix object
  # for later inspection
  if (!is.null(out$full_conf_obj)) {
    # This is only non-NULL if i == 1 and was not skipped
    sample_conf_mat <- out$full_conf_obj
  }
}

# Print iteration messages in sequence
cat("\n==== Iteration Messages ====\n")
for (i in seq_along(iteration_messages)) {
  cat(iteration_messages[i])
}

cat("Number of skipped iterations:", skipped_count, "out of", n_repeats, "\n")

# If iteration #1 was not skipped, sample_conf_mat holds the entire confusionMatrix
if (!is.null(sample_conf_mat)) {
  cat("\n==== Confusion Matrix for iteration #1 (Full caret object) ====\n")
  print(sample_conf_mat)          # prints the table + overall
  cat("\nOverall Statistics:\n")
  print(sample_conf_mat$overall)  # includes Accuracy, Kappa, etc.
  cat("\nBy-Class Statistics:\n")
  print(sample_conf_mat$byClass)  # includes Sensitivity, Specificity, etc.
  
} else {
  cat("\nNo full confusionMatrix was stored (iteration #1 may have been skipped).\n")
}

###############################################################################
# Save Baseline Results
###############################################################################
write.csv(results, "repeated_splits_SMOTE_T0_cluster3_GSE69063.csv", row.names = FALSE)

avg_metrics <- colMeans(results, na.rm = TRUE)
sd_metrics  <- apply(results, 2, sd, na.rm = TRUE)

cat("\nBaseline (All Features) Average Metrics with SMOTE:\n")
print(avg_metrics)

cat("\nStandard Deviations of Metrics:\n")
print(sd_metrics)

average_metrics_df <- data.frame(
  Metric = colnames(results),
  Mean   = avg_metrics,
  SD     = sd_metrics
)
write.csv(average_metrics_df, "average_metrics_SMOTE_T0_cluster3_GSE69063.csv", row.names = FALSE)

# Save feature importance data
write.csv(all_feature_importance, "feature_importance_SMOTE_T0_cluster3_GSE69063.csv", row.names = FALSE)

# Compute average importance
avg_importance <- all_feature_importance %>%
  group_by(Feature) %>%
  summarize(
    Mean_Importance = mean(MeanDecreaseGini, na.rm = TRUE),
    SD_Importance   = sd(MeanDecreaseGini,  na.rm = TRUE)
  ) %>%
  arrange(desc(Mean_Importance))

write.csv(avg_importance, "average_feature_importance_SMOTE_T0_cluster3_GSE69063.csv", row.names = FALSE)

###############################################################################
# Plot Feature Importance (with error bars)
###############################################################################
top_n_features <- 50
top_features <- avg_importance %>%
  top_n(top_n_features, Mean_Importance)

feature_plot <- ggplot(
  top_features,
  aes(x = reorder(Feature, Mean_Importance), y = Mean_Importance)
) +
  geom_bar(stat = "identity", aes(fill = Mean_Importance), width = 0.7, color = "black") +
  geom_errorbar(aes(
    ymin = Mean_Importance - SD_Importance,
    ymax = Mean_Importance + SD_Importance
  ), width = 0.3) +
  coord_flip() +
  scale_fill_viridis(option = "magma", direction = -1, begin = 0.2, end = 0.9) +
  theme_light(base_size = 14) +
  labs(
    title    = paste("Top", top_n_features, "Features (MeanDecreaseGini ± 1 SD)"),
    subtitle = "Random Forest + SMOTE",
    x        = "Feature",
    y        = "Importance",
    fill     = "Importance"
  ) +
  theme(
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),
    axis.title.x  = element_text(size = 14, face = "bold"),
    axis.title.y  = element_text(size = 14, face = "bold"),
    axis.text.x   = element_text(size = 12),
    axis.text.y   = element_text(size = 7),
    legend.position = "right"
  ) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 10))

ggsave("top_feature_importance_plot_SMOTE_T0_cluster3_GSE69063.png", 
       plot = feature_plot, width = 12, height = 7, dpi = 300)

###############################################################################
# 3) Feature Removal Analysis (Re-using the SAME Partitions)
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

for (f_idx in seq_along(features)) {
  feature_to_remove <- features[f_idx]
  cat("Removing feature:", feature_to_remove, "\n")
  
  modified_data <- data[, !names(data) %in% feature_to_remove]
  
  fr_metrics <- matrix(NA, nrow = n_repeats, ncol = 7)
  colnames(fr_metrics) <- c("MCC","F1","AUC","TPR","TNR","PPV","NPV")
  
  removal_list <- foreach(i = seq_len(n_repeats), .packages = c("randomForest","caret","smotefamily","pROC")) %dopar% {
    train_index <- all_partitions[[i]]
    train_data  <- modified_data[train_index, ]
    test_data   <- modified_data[-train_index, ]
    
    train_data$Label <- factor(train_data$Label, levels = c("healthy control","sepsis"))
    test_data$Label  <- factor(test_data$Label, levels = c("healthy control","sepsis"))
    
    if (length(unique(train_data$Label)) < 2) {
      return(rep(NA,7))
    }
    
    n_Sepsis  <- sum(train_data$Label == "sepsis")
    n_Healthy <- sum(train_data$Label == "healthy control")
    
    if (n_Sepsis == 0 || n_Healthy == 0) {
      return(rep(NA,7))
    }
    
    dup_size_value <- max(
      as.integer(ceiling((n_Healthy / n_Sepsis) - 1)),
      0
    )
    
    smote_output <- SMOTE(
      X        = train_data[, -ncol(train_data)],
      target   = train_data$Label,
      K        = 5,
      dup_size = dup_size_value
    )
    
    smote_train_data <- smote_output$data
    colnames(smote_train_data)[ncol(smote_train_data)] <- "Label"
    smote_train_data$Label <- factor(smote_train_data$Label, levels = c("healthy control","sepsis"))
    
    rf_model <- randomForest(
      Label ~ .,
      data  = smote_train_data,
      ntree = 100,
      mtry  = floor(sqrt(ncol(smote_train_data) - 1))
    )
    
    probs <- predict(rf_model, newdata = test_data, type = "prob")[, "sepsis"]
    preds <- ifelse(probs > 0.5, "sepsis", "healthy control")
    
    if (length(unique(test_data$Label)) < 2 || length(unique(preds)) < 2) {
      return(rep(NA,7))
    }
    
    confusion <- confusionMatrix(
      data      = factor(preds, levels = c("healthy control","sepsis")),
      reference = test_data$Label,
      positive  = "sepsis"
    )
    
    cm_table <- confusion$table
    TP <- cm_table[2,2]
    FP <- cm_table[1,2]
    FN <- cm_table[2,1]
    TN <- cm_table[1,1]
    
    AUC_val <- if (length(unique(probs)) > 1) {
      pROC::auc(as.numeric(test_data$Label), probs)
    } else {
      NA
    }
    
    denom  <- (TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)
    MCC_val <- if (denom > 0) {
      ((TP * TN) - (FP * FN)) / sqrt(denom)
    } else {
      NA
    }
    
    c(
      MCC_val,
      confusion$byClass["F1"],
      AUC_val,
      confusion$byClass["Sensitivity"],
      confusion$byClass["Specificity"],
      confusion$byClass["Pos Pred Value"],
      confusion$byClass["Neg Pred Value"]
    )
  }
  
  for (i in seq_len(n_repeats)) {
    fr_metrics[i,] <- removal_list[[i]]
  }
  
  avg_metrics <- colMeans(fr_metrics, na.rm = TRUE)
  
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

write.csv(feature_removal_results, "feature_removal_results_SMOTE_T0_cluster3_GSE69063.csv", row.names = FALSE)

# Compare "All Features" vs. "Removed Feature"
baseline_mcc <- avg_metrics["MCC"]  # from the baseline 'avg_metrics'
cat("\nBaseline MCC (using all features):", baseline_mcc, "\n\n")

mcc_comparison_df <- rbind(
  data.frame(Feature = "All Features (Baseline)", MCC = baseline_mcc),
  feature_removal_results[, c("Feature","MCC")]
)

mcc_comparison_df$Feature <- reorder(mcc_comparison_df$Feature, -mcc_comparison_df$MCC)

mcc_plot <- ggplot(mcc_comparison_df, aes(x = Feature, y = MCC, fill = MCC)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  coord_flip() +
  scale_fill_gradient(low = "#E1735A", high = "#138A78") +
  theme_light(base_size = 14) +
  labs(
    title    = "Impact of Feature Removal on MCC Score",
    subtitle = "Baseline vs. Individual Feature Removal",
    x        = "Feature Removed",
    y        = "MCC Score",
    fill     = "MCC"
  ) +
  theme(
    plot.title    = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    axis.title.x  = element_text(size = 14, face = "bold"),
    axis.title.y  = element_text(size = 8,  face = "bold"),
    axis.text.x   = element_text(size = 12),
    axis.text.y   = element_text(size = 7),
    legend.position = "right"
  ) +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 10))

ggsave("feature_removal_mcc_plot_SMOTE_T0_cluster3_GSE69063.png", 
       plot = mcc_plot, width = 10, height = 12, dpi = 300)

###############################################################################
# 4) Sanity Check: Adding Noise + Using Same Partitions
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
  cat("Noise level:", noise_level, "%\n")
  
  # Make a noisy copy
  noisy_data <- data
  
  if (noise_level > 0) {
    set.seed(999 + noise_level)  # optional seed for noise
    num_noisy_rows <- floor(nrow(data) * noise_level / 100)
    noisy_rows <- sample(seq_len(nrow(data)), num_noisy_rows)
    
    numeric_cols <- which(sapply(noisy_data, is.numeric))
    for (r in noisy_rows) {
      for (col in numeric_cols) {
        noisy_data[r, col] <- noisy_data[r, col] * runif(1, min = 0.8, max = 1.2)
      }
    }
  }
  
  noise_metrics <- matrix(NA, nrow = n_repeats, ncol = 7)
  colnames(noise_metrics) <- c("MCC","F1","AUC","TPR","TNR","PPV","NPV")
  
  noise_results_list <- foreach(i = seq_len(n_repeats), 
                                .packages = c("randomForest","caret","smotefamily","pROC")) %dopar% {
                                  train_index <- all_partitions[[i]]
                                  train_data  <- noisy_data[train_index, ]
                                  test_data   <- noisy_data[-train_index, ]
                                  
                                  train_data$Label <- factor(train_data$Label, levels = c("healthy control","sepsis"))
                                  test_data$Label  <- factor(test_data$Label, levels = c("healthy control","sepsis"))
                                  
                                  n_Sepsis  <- sum(train_data$Label == "sepsis")
                                  n_Healthy <- sum(train_data$Label == "healthy control")
                                  
                                  if (n_Sepsis == 0 || n_Healthy == 0) {
                                    return(rep(NA,7))
                                  }
                                  
                                  dup_size_value <- max(
                                    as.integer(ceiling((n_Healthy / n_Sepsis) - 1)),
                                    0
                                  )
                                  
                                  smote_output <- SMOTE(
                                    X        = train_data[, -ncol(train_data)],
                                    target   = train_data$Label,
                                    K        = 5,
                                    dup_size = dup_size_value
                                  )
                                  
                                  smote_train_data <- smote_output$data
                                  colnames(smote_train_data)[ncol(smote_train_data)] <- "Label"
                                  smote_train_data$Label <- factor(smote_train_data$Label, levels = c("healthy control","sepsis"))
                                  
                                  rf_model <- randomForest(
                                    Label ~ .,
                                    data  = smote_train_data,
                                    ntree = 100,
                                    mtry  = floor(sqrt(ncol(smote_train_data) - 1))
                                  )
                                  
                                  probs <- predict(rf_model, newdata = test_data, type = "prob")[, "sepsis"]
                                  preds <- ifelse(probs > 0.5, "sepsis", "healthy control")
                                  
                                  if (length(unique(test_data$Label)) < 2 || length(unique(preds)) < 2) {
                                    return(rep(NA,7))
                                  }
                                  
                                  confusion <- confusionMatrix(
                                    data      = factor(preds, levels = c("healthy control","sepsis")),
                                    reference = test_data$Label,
                                    positive  = "sepsis"
                                  )
                                  
                                  cm_table <- confusion$table
                                  TP <- cm_table[2,2]
                                  FP <- cm_table[1,2]
                                  FN <- cm_table[2,1]
                                  TN <- cm_table[1,1]
                                  
                                  AUC_val <- if (length(unique(probs)) > 1) {
                                    pROC::auc(as.numeric(test_data$Label), probs)
                                  } else {
                                    NA
                                  }
                                  
                                  denom  <- (TP + FP)*(TP + FN)*(TN + FP)*(TN + FN)
                                  MCC_val <- if (denom > 0) {
                                    ((TP * TN) - (FP * FN)) / sqrt(denom)
                                  } else {
                                    NA
                                  }
                                  
                                  c(
                                    MCC_val,
                                    confusion$byClass["F1"],
                                    AUC_val,
                                    confusion$byClass["Sensitivity"],
                                    confusion$byClass["Specificity"],
                                    confusion$byClass["Pos Pred Value"],
                                    confusion$byClass["Neg Pred Value"]
                                  )
                                }
  
  for (i in seq_len(n_repeats)) {
    noise_metrics[i,] <- noise_results_list[[i]]
  }
  
  avg_nm <- colMeans(noise_metrics, na.rm = TRUE)
  
  sanity_check_results <- rbind(
    sanity_check_results,
    data.frame(
      Noise_Level = noise_level,
      MCC         = avg_nm[1],
      F1          = avg_nm[2],
      AUC         = avg_nm[3],
      TPR         = avg_nm[4],
      TNR         = avg_nm[5],
      PPV         = avg_nm[6],
      NPV         = avg_nm[7]
    )
  )
}

write.csv(sanity_check_results, "sanity_check_results_SMOTE_T0_cluster3_GSE69063.csv", row.names = FALSE)

# Plot noise vs. metrics
noise_plot <- ggplot(sanity_check_results, aes(x = Noise_Level)) +
  geom_line(aes(y = MCC, color = "MCC"), size = 1.2) +
  geom_point(aes(y = MCC, color = "MCC"), size = 3, shape = 16) +
  
  geom_line(aes(y = F1,  color = "F1"),  size = 1.2) +
  geom_point(aes(y = F1, color = "F1"),  size = 3, shape = 17) +
  
  geom_line(aes(y = AUC, color = "AUC"), size = 1.2) +
  geom_point(aes(y = AUC, color = "AUC"), size = 3, shape = 15) +
  
  labs(
    title    = "Impact of Noise on Model Performance",
    subtitle = "Model robustness at different noise levels",
    x        = "Noise Level (%)",
    y        = "Performance Metric",
    color    = "Metrics"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),
    axis.title.x  = element_text(size = 14, face = "bold"),
    axis.title.y  = element_text(size = 14, face = "bold"),
    axis.text.x   = element_text(size = 12),
    axis.text.y   = element_text(size = 12),
    legend.position = "top",
    legend.title    = element_text(size = 12, face = "bold")
  ) +
  scale_color_manual(values = c(
    "MCC" = "#098689", 
    "F1"  = "#F75F5C", 
    "AUC" = "#E69F00"
  ))

ggsave("noise_plot_SMOTE_T0_cluster3_GSE69063.png", 
       plot = noise_plot, width = 10, height = 6, dpi = 300)

###############################################################################
# Shut down the parallel cluster
###############################################################################
stopCluster(cl)

cat("\nAll done.\n")

