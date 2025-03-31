##################################################################################
#-------CREATE CORRELATION MATRIX BETWEEN VARIABLES - ZG_DAYS VS. SPECTRAL, ETC.-----

# Define the columns of interest (columns 3 to 2830)
all_vars <- names(df_rf)[3:2830]
# Exclude 'zg_fraction' if it's present in those columns
other_vars <- setdiff(all_vars, "zg_fraction")

# Compute correlations between 'zg_fraction' and the selected variables
correlations <- sapply(other_vars, function(x) {
  cor(df_rf$zg_fraction, df_rf[[x]], use = "complete.obs")
})
cor_df <- data.frame(variable = other_vars, correlation = correlations)

# Sort the dataframe by the absolute value of the correlation (largest first)
cor_df <- cor_df[order(-abs(cor_df$correlation)), ]
head(cor_df)  # view the top correlations

write.csv(cor_df,"./R_outputs/modelling/correlations/zg_days_VIstatsMERGED.csv")

# Optional: function to filter correlations based on a threshold
filter_correlations <- function(cor_df, threshold, metric = "correlation") {
  if (metric == "correlation") {
    # Filter by absolute correlation value
    return(cor_df[abs(cor_df$correlation) > threshold, ])
  } else if (metric == "rsquared") {
    # Filter by R-squared value (correlation squared)
    return(cor_df[(cor_df$correlation)^2 > threshold, ])
  } else {
    stop("Invalid metric. Choose 'correlation' or 'rsquared'.")
  }
}

# Example: Filter for variables with an absolute correlation greater than 0.3
filtered_df <- filter_correlations(cor_df, threshold = 0.3, metric = "correlation")
head(filtered_df)

#####################################################################################################
#------------------TESTING OF MODELS USING PCA INPUTS------------------------------------
# Load required packages
library(caret)
library(randomForest)
library(e1071)     # for SVM
# xgboost is supported by caret (xgbTree), so no explicit library(xgboost) is needed here

# Define dynamic dependent variable and predictors (assuming PCA has been run and PC scores added)
dvarb <- "zg_fraction"
# For example, using only PC1 and PC2 (or include "PC3" if desired)
predictors <- c("PC1", "PC2")  

# Build the model formula dynamically
model_formula <- as.formula(paste(dvarb, "~", paste(predictors, collapse = "+")))
cat("Model formula is:", deparse(model_formula), "\n")

# Define training control (10-fold cross-validation)
set.seed(123)
train_control <- trainControl(method = "cv", number = 5)

# Train a linear regression model
model_lm <- train(model_formula, 
                  data = df_rf, 
                  method = "lm", 
                  trControl = train_control)
# Train a random forest model (even with 2 predictors, we can still try it)
model_rf <- train(model_formula, 
                  data = df_rf, 
                  method = "rf", 
                  trControl = train_control,
                  ntree = 200)
# Train a support vector machine model with a radial kernel
model_svm <- train(model_formula, 
                   data = df_rf, 
                   method = "svmRadial", 
                   trControl = train_control)

# Train an XGBoost model (if available)
model_xgb <- train(model_formula, 
                   data = df_rf, 
                   method = "xgbTree", 
                   trControl = train_control)
# Compare model performance using resampling results
models_list <- list(LinearRegression = model_lm, 
                    RandomForest = model_rf, 
                    SVM = model_svm, 
                    XGBoost = model_xgb)
# Gather resample results
results_resamples <- resamples(models_list)
summary(results_resamples)

# Optionally, visualize performance differences
bwplot(results_resamples)
