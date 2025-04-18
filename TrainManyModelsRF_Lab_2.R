library(boot)
library(randomForest)

#----------------DROP NA COLUMNS IF NEED BE (e.g. ZENTROPY)---------------------
# First, drop rows with NA in the dependent variable (dvarb)
df_rf <- df_rf %>% 
  drop_na(!!sym(dvarb))
# Then, filter out predictors that contain any NA values
clean_predictors <- predictors[!sapply(df_rf[, predictors, drop = FALSE], function(x) any(is.na(x)))]
# Finally, select only the dependent variable and the cleaned predictors
df_rf <- df_rf %>% 
  select(!!sym(dvarb), all_of(clean_predictors))

#-----------------CHOOSE DEPENDENT + PREDICTOR VARIABLES-------------------------
# Dependent variable (can be dynamically changed)
dvarb <- "zg_fraction"  
# Predictor columns (dynamic)
predictors <- colnames(df_rf)[3:544]

#-----------------PCA FOR SPECTRAL LIBRARY-------------------------------------

# Define dependent variable and original predictor columns
original_predictors <- colnames(df_rf)[3:544]

# Run PCA on the original predictors (assuming df_rf has no missing values in these columns)
pca_res <- prcomp(df_rf[, original_predictors], center = TRUE, scale. = TRUE)

# Look at the explained variance (proportion and cumulative)
explained <- summary(pca_res)$importance[2, ]       # Proportion of variance explained
cum_explained <- summary(pca_res)$importance[3, ]     # Cumulative variance explained
print(cum_explained)

# Decide on number of PCs:
# Here we use PC1 and PC2, and include PC3 if it adds at least 5% extra variance (i.e., if cum_explained[3] - cum_explained[2] >= 0.05).
n_pcs <- ifelse((cum_explained[3] - cum_explained[2]) >= 0.05, 3, 2)
cat("Using", n_pcs, "principal components.\n")

# Extract the chosen PC scores and rename them
pc_scores <- as.data.frame(pca_res$x[, 1:n_pcs])
colnames(pc_scores) <- paste0("PC", 1:n_pcs)

# Add the PC scores to the original data frame
df_rf <- cbind(df_rf, pc_scores)

# Update predictors to be the selected PCs
predictors <- paste0("PC", 1:n_pcs)
print(predictors)
#---------------------BOOTSTRAP SETUP------------------------------
# Number of bootstrap replicates
n_boot <- 50  # adjust as needed

# Storage for results
results <- data.frame(
  Iteration = 1:n_boot,
  RMSE_train = numeric(n_boot),
  RMSE_test = numeric(n_boot),
  Rsq_train = numeric(n_boot),
  Rsq_test = numeric(n_boot)
)
#--------------------BOOTSTRAP WITH REGRESSION - PCA INPUT------------------------
set.seed(123)

for(i in 1:n_boot) {
  
  # Create an 80/20 split without replacement so that the test set is guaranteed to have data
  all_indices <- seq_len(nrow(df_rf))
  train_indices <- sample(all_indices, size = floor(0.7 * nrow(df_rf)), replace = FALSE)
  test_indices <- setdiff(all_indices, train_indices)
  
  train_data <- df_rf[train_indices, c(dvarb, predictors)]
  test_data  <- df_rf[test_indices, c(dvarb, predictors)]
  
  # Fit a linear regression model
  formula_reg <- as.formula(paste(dvarb, "~", paste(predictors, collapse = "+")))
  lm_model <- lm(formula_reg, data = train_data)
  
  # Make predictions on both train and test data
  pred_train <- predict(lm_model, newdata = train_data)
  pred_test  <- predict(lm_model, newdata = test_data)
  
  # Get actual values
  actual_train <- train_data[[dvarb]]
  actual_test  <- test_data[[dvarb]]
  
  # Calculate RMSE for train and test
  rmse_train <- sqrt(mean((actual_train - pred_train)^2))
  rmse_test  <- sqrt(mean((actual_test - pred_test)^2))
  
  # Calculate R-squared (using squared correlation)
  rsq_train <- cor(actual_train, pred_train)^2
  rsq_test  <- cor(actual_test, pred_test)^2
  
  # Store the performance metrics
  results$RMSE_train[i] <- rmse_train
  results$RMSE_test[i]  <- rmse_test
  results$Rsq_train[i]  <- rsq_train
  results$Rsq_test[i]   <- rsq_test
}

# View the first few rows of results
head(results)

#---------------------BOOTSTRAP WITH RANDOM FOREST------------------------------
set.seed(123)

# Bootstrap iterations
for(i in 1:n_boot){
  
  # Create 80/20 random split
  sample_indices <- sample(nrow(df_rf), size = 0.7 * nrow(df_rf), replace = TRUE)
  
  train_data <- df_rf[sample_indices, c(dvarb, predictors)]
  test_data  <- df_rf[-unique(sample_indices), c(dvarb, predictors)] # ensure test is different
  
  # Fit Random Forest model dynamically
  formula_rf <- as.formula(paste(dvarb, "~", paste(predictors, collapse = "+")))
  
  rf_model <- randomForest(
    formula_rf, 
    data = train_data,
    ntree = 200
  )
  
  # Predictions
  pred_train <- predict(rf_model, train_data)
  pred_test <- predict(rf_model, test_data)
  
  # Calculate performance metrics
  actual_train <- train_data[[dvarb]]
  actual_test <- test_data[[dvarb]]
  
  # RMSE and R-squared
  rmse_train <- sqrt(mean((actual_train - pred_train)^2))
  rmse_test  <- sqrt(mean((actual_test - pred_test)^2))
  
  rsq_train <- cor(actual_train, pred_train)^2
  rsq_test  <- cor(actual_test, pred_test)^2
  
  # Store results
  results$RMSE_train[i] <- rmse_train
  results$RMSE_test[i] <- rmse_test
  results$Rsq_train[i] <- rsq_train
  results$Rsq_test[i] <- rsq_test
}

# View summarized results
head(results)

# ----------Summarize variability ACROSS BOOTSTRAPS------------------------------
summary_results <- data.frame(
  Metric = c("RMSE_train", "RMSE_test", "Rsq_train", "Rsq_test"),
  Mean = c(mean(results$RMSE_train), mean(results$RMSE_test),
           mean(results$Rsq_train), mean(results$Rsq_test)),
  SD = c(sd(results$RMSE_train), sd(results$RMSE_test),
         sd(results$Rsq_train), sd(results$Rsq_test)),
  CI_lower = c(
    quantile(results$RMSE_train, 0.025),
    quantile(results$RMSE_test, 0.025),
    quantile(results$Rsq_train, 0.025),
    quantile(results$Rsq_test, 0.025)
  ),
  CI_upper = c(
    quantile(results$RMSE_train, 0.975),
    quantile(results$RMSE_test, 0.975),
    quantile(results$Rsq_train, 0.975),
    quantile(results$Rsq_test, 0.975)
  )
)

# Print summary
summary_results

################################################################################

# Function to bootstrap correlation for a given predictor
boot_corr_fn <- function(data, indices, predictor){
  sampled_data <- data[indices, ]
  return(cor(sampled_data$TWD, sampled_data[[predictor]], use="complete.obs"))
}

# Prepare storage for results
cor_results <- data.frame(
  Predictor = predictor_names,
  Correlation = NA,
  CI_Lower = NA,
  CI_Upper = NA
)

# Run bootstrapping for correlations (can take some minutes)
set.seed(123)

for(i in seq_along(predictor_names)){
  predictor <- predictor_names[i]
  boot_res <- boot(df_rf, boot_corr_fn, R=1000, predictor=predictor)
  
  cor_results$Correlation[i] <- mean(boot_res$t)
  ci <- boot.ci(boot_res, type="perc")
  cor_results$CI_Lower[i] <- ci$percent[4]  # 2.5% lower bound
  cor_results$CI_Upper[i] <- ci$percent[5]  # 97.5% upper bound
}

# Display the first few correlation results
head(cor_results)

########################################################################################################
# ---------------------------RANDOM FOREST WITH REPEATS-----------------------------------------------
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
