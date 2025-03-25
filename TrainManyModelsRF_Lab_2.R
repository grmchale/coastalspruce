library(dplyr)
library(caret)
library(ranger)
library(randomForest)

# ---------------------------
# 1. Prepare Data
# ---------------------------
# Define predictor variables (all columns except the dynamic response variable)
predictor_vars <- names(df_rf)[!(names(df_rf) %in% dvarb)]
predictor_vars <- names(df_rf)[2:486]
print(predictor_vars)

# Define response variable (e.g., "TWD")
response_var <- dvarb

# (Optional) RF model name (for later use)
RF_NAME <- "5nm_VIstats"

# Select data for modeling: keep only the response and predictor variables and remove rows with missing values.
model_data <- df_rf %>%
  select(all_of(c(response_var, predictor_vars))) %>%
  drop_na()

# ---------------------------
# 2. Hyperparameter Tuning using caret
# ---------------------------
# Set up repeated cross-validation (here: 5-fold repeated 5 times)
set.seed(1234)  # For reproducibility
train_control <- trainControl(method = "repeatedcv",
                              number = 5,
                              repeats = 5,
                              savePredictions = TRUE)

# Use caret to tune the model using the ranger engine. We focus on maximizing Rsquared.
rf_tuned <- train(
  TWD ~ ., 
  data = model_data,
  method = "ranger",
  tuneLength = 5,      # Tests a range of hyperparameters
  num.trees = 1000,
  importance = "impurity",
  trControl = train_control
)

# Display best tuned hyperparameters and results
print(rf_tuned$bestTune)
print(rf_tuned$results)

# Optionally, display the results for the best model configuration:
best_results <- rf_tuned$results %>%
  filter(mtry == rf_tuned$bestTune$mtry,
         min.node.size == rf_tuned$bestTune$min.node.size,
         splitrule == rf_tuned$bestTune$splitrule)
print(best_results)

# ---------------------------
# 3. Ensemble Runs to Average Performance and Variable Importance
# ---------------------------
# We now run 50 independent randomForest models using the best mtry from caret.
# Note: We assume the response variable is "TWD" as used in the formula.
set.seed(123)  # For reproducibility of the ensemble runs
n_runs <- 50

# Define the predictor names for importance matrix (exclude the response variable)
predictors <- setdiff(names(model_data), response_var)

# Initialize a data frame to store evaluation metrics for each run
rf_metrics <- data.frame(Rsquared = numeric(n_runs),
                         RMSE = numeric(n_runs),
                         MAE = numeric(n_runs))

# Initialize a matrix to store variable importance (permutation importance)
imp_matrix <- matrix(NA, nrow = length(predictors), ncol = n_runs,
                     dimnames = list(predictors, NULL))

for (i in 1:n_runs) {
  rf <- randomForest(
    TWD ~ ., 
    data = model_data,
    ntree = 500, 
    mtry = rf_tuned$bestTune$mtry, 
    importance = TRUE
  )
  
  # Extract OOB (out-of-bag) predictions
  preds <- rf$predicted
  actuals <- model_data$TWD
  
  # Compute evaluation metrics
  rmse_val <- sqrt(mean((preds - actuals)^2))
  mae_val <- mean(abs(preds - actuals))
  # For randomForest, OOB R² is stored in the last element of the rsq vector.
  rsq_val <- rf$rsq[length(rf$rsq)]
  
  rf_metrics[i, ] <- c(rsq_val, rmse_val, mae_val)
  
  # Extract permutation-based variable importance (MeanDecreaseAccuracy)
  imp <- importance(rf, type = 1)[, "%IncMSE"]
  imp_matrix[names(imp), i] <- imp
}

# Average the evaluation metrics across the 50 runs
avg_metrics <- colMeans(rf_metrics)
print("Average Evaluation Metrics from 50 RF runs:")
print(avg_metrics)

# Optionally, also view the standard deviation for insight into variability:
sd_metrics <- apply(rf_metrics, 2, sd)
print("Standard Deviation of Evaluation Metrics:")
print(sd_metrics)

# Average variable importance across the 50 runs
avg_imp <- rowMeans(imp_matrix, na.rm = TRUE)
avg_imp_sorted <- sort(avg_imp, decreasing = TRUE)
print("Average Variable Importance (MeanDecreaseAccuracy):")
print(avg_imp_sorted)
