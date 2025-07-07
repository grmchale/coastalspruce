##################################################################################
#-------CREATE CORRELATION MATRIX BETWEEN SEVERAL VARIABLES----------------------
# Ensure the output directory exists
dir.create("./R_outputs/modelling/correlations/", recursive = TRUE, showWarnings = FALSE)

# Define the column indices clearly
spectral_cols <- 2:1961       # Columns of spectral variables
zg_cols <- 1962:1967          # Columns of zg_fraction variants

# Extract variable names
spectral_vars <- names(df_rf)[spectral_cols]
zg_vars <- names(df_rf)[zg_cols]

# Initialize an empty correlation dataframe (1960 rows, 6 columns)
cor_matrix <- data.frame(variable = spectral_vars)

# Calculate Spearman correlations and store results
for (zg_var in zg_vars) {
  cor_matrix[[zg_var]] <- sapply(spectral_vars, function(spec_var) {
    cor(df_rf[[zg_var]], df_rf[[spec_var]], method = "spearman", use = "complete.obs")
  })
}

# Export the combined correlation matrix
write.csv(cor_matrix, "./R_outputs/modelling/correlations/zg_fractions_all_VIstats.csv", row.names = FALSE)

# Additionally, create and export separate sorted dataframes for each zg_fraction variant
for (zg_var in zg_vars) {
  cor_df_individual <- data.frame(
    variable = spectral_vars,
    correlation = cor_matrix[[zg_var]]
  )
  
  # Sort by absolute correlation (strongest first)
  cor_df_individual <- cor_df_individual[order(-abs(cor_df_individual$correlation)), ]
  
  # Export each dataframe
  filename <- paste0("./R_outputs/modelling/correlations/", zg_var, "_correlations.csv")
  write.csv(cor_df_individual, filename, row.names = FALSE)
}

##################################################################################
#-------CREATE CORRELATION MATRIX BETWEEN 2 VARIABLES - ZG_DAYS VS. SPECTRAL, ETC.-----

# Define the columns of interest (columns 3 to 2830)
all_vars <- names(df_rf)[3:2830]
# Exclude 'zg_fraction' if it's present in those columns
other_vars <- setdiff(all_vars, "zg_fraction")

# Compute correlations between 'zg_fraction' and the selected variables
correlations <- sapply(other_vars, function(x) {
  cor(df_rf$zg_fraction, df_rf[[x]], method = "spearman", use = "complete.obs")
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
cor04_df <- filter_correlations(cor_df, threshold = 0.4, metric = "correlation")
head(filtered_df)

#####################################################################################################
#----------------- REMOVE CORRELATED PREDICTORS USING SPEARMAN MATRIX----------------
# Load libraries
library(corrplot)
library(caret)

# PARAMETERS
top_n <- 100        # number of top correlated variables to select
cor_threshold <- 0.9  # intercorrelation threshold for removal

# STEP 1: Filter top N variables by absolute Spearman correlation
top_vars_df <- cor_df[order(-abs(cor_df$correlation)), ][1:top_n, ]
top_vars <- top_vars_df$variable

# STEP 2: Extract subset of data
top_data <- df_rf[, top_vars]

# OPTION A — MANUAL SPEARMAN CORRELATION MATRIX
# Compute and visualize Spearman correlation matrix
spearman_matrix <- cor(top_data, method = "spearman", use = "complete.obs")

# Visual inspection
corrplot(spearman_matrix, method = "color", type = "lower", tl.cex = 0.6,
         order = "hclust", addrect = 3)

library(htmlwidgets)

# Create a plot and export as an interactive HTML
html_file <- "R_outputs/modelling/spearman_corrplot_interactive.html"

# Use htmlwidgets for zoomable image (not full interactivity, but high-res)
svg("R_outputs/modelling/zoomed_corrplot.svg", width = 12, height = 12)
corrplot(spearman_matrix, method = "color", type = "lower", tl.cex = 0.7,
         order = "hclust", addrect = 3)
dev.off()


# Extract highly correlated variable pairs
high_corr_pairs <- which(abs(spearman_matrix) > cor_threshold & abs(spearman_matrix) < 1, arr.ind = TRUE)
manual_corr_info <- data.frame(
  var1 = rownames(spearman_matrix)[high_corr_pairs[, 1]],
  var2 = colnames(spearman_matrix)[high_corr_pairs[, 2]],
  correlation = spearman_matrix[high_corr_pairs]
)

# OPTIONAL: Print or inspect top pairs
head(manual_corr_info)

# STEP 3 (OPTION B) — AUTOMATED REMOVAL WITH CARET
remove_vars <- findCorrelation(spearman_matrix, cutoff = cor_threshold, names = TRUE)
reduced_vars <- setdiff(top_vars, remove_vars)

# Output: list of final variables
cat("Final variable count:", length(reduced_vars), "\n")
print(reduced_vars[1:10])  # show first few

# Save results
write.csv(reduced_vars, "./R_outputs/modelling/selected_vars_reduced.csv", row.names = FALSE)
write.csv(manual_corr_info, "./R_outputs/modelling/intercorrelated_pairs_manual.csv", row.names = FALSE)

###################################################################################################
#REMOVE CORRELATED PREDICTORS USING *magic* EMPHASIZING HIGHLY CORRELATED VARIABLES TO ZG_FRACTION
library(caret)
library(dplyr)
library(stringr)

# Parameters
top_n <- 300               # Number of top variables based on correlation
cor_threshold <- 0.9        # Threshold to define intercorrelation groups

# Step 1: Select top N variables by absolute Spearman correlation
top_vars_df <- cor_df %>%
  arrange(desc(abs(correlation))) %>%
  slice(1:top_n)

top_vars <- top_vars_df$variable

# Step 2: Subset original data frame for selected variables
top_data <- df_rf[, top_vars]

# Step 3: Compute Spearman correlation matrix among predictors
spearman_matrix <- cor(top_data, method = "spearman", use = "complete.obs")

# Step 4: Identify groups of highly correlated variables
high_corr <- findCorrelation(spearman_matrix, cutoff = cor_threshold, names = TRUE, exact = TRUE)

# Function to prioritize Median/Mean/SD over Q25/Q75
priority_score <- function(var_name) {
  if (str_detect(var_name, "Median")) return(3)
  if (str_detect(var_name, "Mean")) return(2)
  if (str_detect(var_name, "SD")) return(2)
  if (str_detect(var_name, "Q25|Q75")) return(1)
  return(0)
}

# Step 5: Automated selection strategy:
# - For each correlated set, keep the variable with highest absolute correlation to zg_fraction
# - If ties or close values, prioritize Median/Mean/SD over Q25/Q75

# Initialize variables to keep track of removal
vars_to_remove <- c()

# Loop through each high-correlation group
for (var in high_corr) {
  # Find variables highly correlated with current variable
  correlated_group <- colnames(spearman_matrix)[abs(spearman_matrix[var, ]) >= cor_threshold]
  
  # Ensure group has multiple variables
  if (length(correlated_group) <= 1) next
  
  # Get correlations to zg_fraction for this group
  group_cor_df <- top_vars_df %>%
    filter(variable %in% correlated_group) %>%
    mutate(priority = sapply(variable, priority_score)) %>%
    arrange(desc(abs(correlation)), desc(priority))
  
  # Select the best variable to KEEP (highest abs(correlation), higher priority)
  var_to_keep <- group_cor_df$variable[1]
  
  # Mark the others for removal
  vars_to_remove <- c(vars_to_remove, setdiff(correlated_group, var_to_keep))
}

# Final selected variables after removing redundancies
final_selected_vars <- setdiff(top_vars, vars_to_remove)

# Display results clearly
cat("Original top variables:", length(top_vars), "\n")
cat("Variables removed due to redundancy:", length(vars_to_remove), "\n")
cat("Final variable count:", length(final_selected_vars), "\n")

# Show first few final variables
head(final_selected_vars)
print(final_selected_vars)

# (Optional) Save results
write.csv(final_selected_vars, "./R_outputs/modelling/final_selected_vars.csv", row.names = FALSE)

############################################################################################
#--------------QUICK RF BABYYYYYYY--------------------------
library(caret)
library(randomForest)

# Data setup
model_data <- df_rf[, c("zg_fraction", final_selected_vars)]

# RF tuning grid
tune_grid <- expand.grid(
  mtry = c(2, 3, 4, 5),
  splitrule = "variance",
  min.node.size = c(3, 5, 7, 10)
)

# Train control with cross-validation
control <- trainControl(method = "repeatedcv",
                        number = 5,
                        repeats = 3)

# Train the RF model
rf_model <- train(zg_fraction ~ .,
                  data = model_data,
                  method = "ranger",
                  tuneGrid = tune_grid,
                  trControl = control,
                  importance = 'permutation')

# Check results
print(rf_model)
varImp(rf_model)

# Create output directory if it doesn't exist
output_dir <- "G:/HyperspectralUAV/R_outputs/modelling/rf/ALLvariables_zg_fraction_rf"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
# Save the RF model object
saveRDS(rf_model, file = file.path(output_dir, "zg_fraction_rf_model.rds"))
# Save the variable importance table
importance_vals <- varImp(rf_model)$importance
write.csv(importance_vals, file = file.path(output_dir, "zg_fraction_variable_importance.csv"), row.names = TRUE)
# Save model performance metrics (e.g., RMSE, Rsquared)
performance_metrics <- rf_model$results
write.csv(performance_metrics, file = file.path(output_dir, "zg_fraction_model_performance.csv"), row.names = FALSE)

##################################################################################################
#------------------TESTING OF MODELS USING PCA INPUTS (DEPRECIATED)----------------------------------
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
