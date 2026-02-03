# Set working directory to coastalspruce repo

## Read in spectra + dendro data
# These are the resampled spectral libraries, VIs, and dendro metrics for each spruce crown
dendro_spectra <- read.csv(
  file = "./data/dendro_spectra_joined.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

##############################################################################
################# PLSR MODELS TO PREDICT DENDRO ATTRIBUTES ####################
#install.packages("pls")
library(pls)

# Select spectral predictors - 1 nm spectral library
X <- dendro_spectra[, paste0("1nm_", 398:999)]

# remove rows with missing data (if any)
complete_idx <- complete.cases(X, dendro_spectra$TWD_drone, dendro_spectra$cshrink)

X <- X[complete_idx, ]

### 1) PLSR FOR TWD_DRONE - TWD ON DAY OF DRONE FLIGHT ###
Y_twd <- dendro_spectra$TWD_drone[complete_idx]

pls_twd <- plsr(
  Y_twd ~ as.matrix(X),
  ncomp = 7, # sq rt of n (53)
  scale = TRUE,
  validation = "CV"
)

# Choose optimal number of components
validationplot(pls_twd, val.type = "RMSEP")
# Summary
summary(pls_twd)
# R² vs components
plot(R2(pls_twd), legendpos = "topright")
# Predictions
twd_pred <- predict(pls_twd, ncomp = which.min(RMSEP(pls_twd)$val[1, , ]))

### 2) PLSR FOR CSHRINK - CONSECUTIVE DAYS PRE-FLIGHT OF NO NET GROWTH ###
Y_cshrink <- dendro_spectra$cshrink[complete_idx]

pls_cshrink <- plsr(
  Y_cshrink ~ as.matrix(X),
  ncomp = 7, # sq rt of n (53)
  scale = TRUE,
  validation = "CV"
)

# Choose optimal number of components
validationplot(pls_cshrink, val.type = "RMSEP")
# Summary
summary(pls_cshrink)
# R² vs components
plot(R2(pls_cshrink), legendpos = "topright")
# Predictions
cshrink_pred <- predict(pls_cshrink, ncomp = which.min(RMSEP(pls_cshrink)$val[1, , ]))

######## 3) LOOPING ACROSS 4 SPECTRAL RESAMPLES WITH PLS ##########
library(pls)

# Function to fit PLSR and extract performance metrics
fit_plsr_metrics <- function(X, Y, ncomp = 7, scale = TRUE) {
  
  # Fit PLSR model with cross-validation
  pls_model <- plsr(
    Y ~ as.matrix(X),
    ncomp = ncomp,
    scale = scale,
    validation = "CV"
  )
  
  # Find optimal number of components (minimum RMSEP)
  rmsep_vals <- RMSEP(pls_model)$val[1, , ]
  opt_ncomp_index <- which.min(rmsep_vals)
  opt_ncomp <- opt_ncomp_index - 1  # Convert index to actual component number
  
  # Extract metrics
  rmsep_cv <- min(rmsep_vals)  # CV RMSEP at optimal ncomp
  r2_cv <- R2(pls_model)$val[1, , opt_ncomp_index]  # CV R²
  
  # Return metrics as a list
  return(list(
    opt_ncomp = opt_ncomp,
    RMSEP_CV = rmsep_cv,
    R2_CV = r2_cv,
    model = pls_model
  ))
}

# Define spectral libraries (wavelength resolutions)
spectral_resolutions <- c("1nm", "5nm", "10nm", "15nm")

# Define response variables
response_vars <- list(
  TWD_drone = dendro_spectra$TWD_drone,
  cshrink = dendro_spectra$cshrink
)

# Initialize results storage
results <- data.frame(
  Library = character(),
  Predictor = character(),
  nComp_optimal = integer(),
  RMSEP_CV = numeric(),
  R2_CV = numeric(),
  stringsAsFactors = FALSE
)

# Loop over spectral libraries and response variables
for (resolution in spectral_resolutions) {
  
  # Automatically detect all columns for this resolution
  pattern <- paste0("^", resolution, "_")
  spec_cols <- grep(pattern, names(dendro_spectra), value = TRUE)
  
  if (length(spec_cols) == 0) {
    warning(paste("No columns found for resolution:", resolution))
    next
  }
  
  X <- dendro_spectra[, spec_cols]
  
  # Loop over response variables
  for (response_name in names(response_vars)) {
    
    Y <- response_vars[[response_name]]
    
    # Remove rows with missing data
    complete_idx <- complete.cases(X, Y)
    X_clean <- X[complete_idx, ]
    Y_clean <- Y[complete_idx]
    
    # Fit model and get metrics
    metrics <- fit_plsr_metrics(X_clean, Y_clean, ncomp = 7)
    
    # Store results
    results <- rbind(results, data.frame(
      Library = resolution,
      Predictor = response_name,
      nComp_optimal = metrics$opt_ncomp,
      RMSEP_CV = metrics$RMSEP_CV,
      R2_CV = metrics$R2_CV
    ))
    
    # Print progress
    cat(sprintf("Completed: %s | %s (n columns = %d)\n", 
                resolution, response_name, length(spec_cols)))
  }
}

# Display results table
print(results)

### 4) PLSR GRAPH FOR CSHRINK ###
# Optimal number of components
n_opt <- which.min(RMSEP(pls_cshrink)$val[1, , ])

# Cross-validated predictions
cshrink_cv_pred <- pls_cshrink$validation$pred[, , n_opt]

# Observed values
cshrink_obs <- Y_cshrink

# CV R²
r2_cv <- cor(cshrink_obs, cshrink_cv_pred)^2

# RMSEP
rmsep_cv <- RMSEP(pls_cshrink)$val[1, n_opt + 1, 1]

# Plot
plot(
  cshrink_obs,
  cshrink_cv_pred,
  pch = 16,
  xlab = "Observed cshrink",
  ylab = "PLSR-predicted cshrink (CV)",
  main = paste0(
    "PLSR Prediction of cshrink (", n_opt, " components)\n",
    "CV R² = ", round(r2_cv, 2),
    ", RMSEP = ", round(rmsep_cv, 2),
    ", n = ", length(cshrink_obs)
  )
)
abline(0, 1, lwd = 2, lty = 2)
abline(lm(cshrink_cv_pred ~ cshrink_obs), col = "black")

#############################################################################
################ RF MODELS TO PREDICT DENDRO ATTRIBUTES ######################
# Filter out stats first before RF'ing
dendro_clean <- dendro_spectra[
  , !grepl("(Q75|Q25|Mean|SD)$", names(dendro_spectra))
]

library(randomForest)
library(caret)

######### 1) RF BOOSTRAP LOOP FOR TWD_DRONE #################
library(randomForest)
library(caret)  # for some helper functions if needed

### CUSTOMIZABLE PARAMETERS ###

# Response variable (cshrink or TWD_drone)
response_var <- "cshrink"

# Spectral library resolutions to test
resolutions <- c("1nm", "5nm", "10nm", "15nm")

# Random Forest parameters
ntree <- 500  # number of trees (same for all models)
# mtry will be auto-calculated as sqrt(p) for each resolution

# Bootstrap parameters
n_bootstrap <- 100  # number of bootstrap iterations
train_prop <- 0.8   # proportion for training (80/20 split)

# Data frame
data <- dendro_clean

# Set seed for reproducibility
set.seed(123)

#### COLUMN EXTRACTION SETUP ###

# Function to extract predictor columns for a given resolution
get_predictors <- function(resolution, data) {
  
  # Pattern for spectral wavelengths: "1nm_", "5nm_", etc.
  spec_pattern <- paste0("^", resolution, "_")
  spec_cols <- grep(spec_pattern, names(data), value = TRUE)
  
  # Pattern for vegetation indices: "index_1nm_Median", "index_5nm_Median", etc.
  index_pattern <- paste0("_", resolution, "_Median$")
  index_cols <- grep(index_pattern, names(data), value = TRUE)
  
  # Combine both
  predictor_cols <- c(spec_cols, index_cols)
  
  # Print summary
  cat(sprintf("%s: %d spectral + %d indices = %d total predictors\n", 
              resolution, length(spec_cols), length(index_cols), length(predictor_cols)))
  
  return(predictor_cols)
}

### BOOTSTRAP RANDOM FOREST FUNCTION ###

bootstrap_rf <- function(X, Y, n_bootstrap, train_prop, ntree) {
  
  n <- nrow(X)
  n_train <- floor(n * train_prop)
  
  # Storage for bootstrap results
  rmse_vec <- numeric(n_bootstrap)
  r2_vec <- numeric(n_bootstrap)
  
  # Run bootstrap iterations
  for (i in 1:n_bootstrap) {
    
    # Random 80/20 split
    train_idx <- sample(1:n, n_train, replace = FALSE)
    test_idx <- setdiff(1:n, train_idx)
    
    X_train <- X[train_idx, ]
    Y_train <- Y[train_idx]
    X_test <- X[test_idx, ]
    Y_test <- Y[test_idx]
    
    # Fit Random Forest
    rf_model <- randomForest(
      x = X_train,
      y = Y_train,
      ntree = ntree,
      mtry = floor(sqrt(ncol(X_train))),
      importance = FALSE  # Skip importance for bootstrap iterations (saves time)
    )
    
    # Predict on test set
    Y_pred <- predict(rf_model, X_test)
    
    # Calculate metrics
    rmse_vec[i] <- sqrt(mean((Y_test - Y_pred)^2))
    r2_vec[i] <- 1 - sum((Y_test - Y_pred)^2) / sum((Y_test - mean(Y_test))^2)
    
    # Progress indicator
    if (i %% 20 == 0) {
      cat(sprintf("  Bootstrap iteration %d/%d complete\n", i, n_bootstrap))
    }
  }
  
  # Return summary statistics
  return(list(
    rmse_mean = mean(rmse_vec),
    rmse_sd = sd(rmse_vec),
    r2_mean = mean(r2_vec),
    r2_sd = sd(r2_vec),
    rmse_all = rmse_vec,
    r2_all = r2_vec
  ))
}

#MAIN LOOP: FIT RF FOR EACH RESOLUTION

# Storage for results
results_summary <- data.frame(
  Resolution = character(),
  n_predictors = integer(),
  RMSE_mean = numeric(),
  RMSE_sd = numeric(),
  R2_mean = numeric(),
  R2_sd = numeric(),
  OOB_error = numeric(),
  stringsAsFactors = FALSE
)

# Storage for variable importance (top 20 for each model)
var_importance_list <- list()

# Storage for full models
models_list <- list()

# Loop over resolutions
for (res in resolutions) {
  
  cat(sprintf("\n========== Processing %s ==========\n", res))
  
  # Extract predictors
  predictor_cols <- get_predictors(res, data)
  cat(sprintf("Number of predictors: %d\n", length(predictor_cols)))
  
  # Prepare data
  X <- data[, predictor_cols]
  Y <- data[[response_var]]
  
  # Remove rows with missing data
  complete_idx <- complete.cases(X, Y)
  X_clean <- X[complete_idx, ]
  Y_clean <- Y[complete_idx]
  
  cat(sprintf("Sample size after removing NAs: %d\n", nrow(X_clean)))
  
  # Run bootstrap
  cat("Running bootstrap iterations...\n")
  bootstrap_results <- bootstrap_rf(X_clean, Y_clean, n_bootstrap, train_prop, ntree)
  
  # Fit final model on full data for OOB error and variable importance
  cat("Fitting final model on full data for OOB and variable importance...\n")
  final_model <- randomForest(
    x = X_clean,
    y = Y_clean,
    ntree = ntree,
    mtry = floor(sqrt(ncol(X_clean))),
    importance = TRUE
  )
  
  # Extract OOB error (MSE) and convert to RMSE
  oob_mse <- final_model$mse[ntree]
  oob_rmse <- sqrt(oob_mse)
  
  # Get variable importance (top 20)
  var_imp <- importance(final_model)
  var_imp_sorted <- var_imp[order(var_imp[, 1], decreasing = TRUE), , drop = FALSE]
  top20_imp <- head(var_imp_sorted, 20)
  var_importance_list[[res]] <- top20_imp
  
  # Store model
  models_list[[res]] <- final_model
  
  # Store results
  results_summary <- rbind(results_summary, data.frame(
    Resolution = res,
    n_predictors = length(predictor_cols),
    RMSE_mean = bootstrap_results$rmse_mean,
    RMSE_sd = bootstrap_results$rmse_sd,
    R2_mean = bootstrap_results$r2_mean,
    R2_sd = bootstrap_results$r2_sd,
    OOB_RMSE = oob_rmse
  ))
  
  cat(sprintf("✓ %s complete!\n", res))
}

# DISPLAY RESULTS

cat("\n========== RESULTS SUMMARY ==========\n")
print(results_summary)

# Display top 10 most important variables for each resolution
cat("\n========== TOP 10 VARIABLE IMPORTANCE ==========\n")
for (res in resolutions) {
  cat(sprintf("\n%s:\n", res))
  print(head(var_importance_list[[res]], 10))
}
# SAVE RESULTS WITH RESPONSE VARIABLE NAME

# Create dynamically named results dataframe
results_df_name <- paste0("results_", response_var)
assign(results_df_name, results_summary)

cat(sprintf("\n✓ Results saved as: %s\n", results_df_name))


#### 2) RF FOR TWD_DRONE ######
# Create predictor matrix
X <- dendro_clean[
                  which(names(dendro_clean) == "1nm_398"):
                    which(names(dendro_clean) == "Vogelmann4_15nm_Median")]

X <- as.data.frame(X)

# Response
Y_twd <- dendro_clean$TWD_drone

# TWD_DRONE RF MODEL
set.seed(123)

rf_twd <- randomForest(
  x = X_twd_sel,
  y = Y_twd,
  ntree = 1000,
  importance = TRUE
)

# Predictions
twd_pred <- predict(rf_twd, X_twd_sel)

# Performance metrics
twd_rmse <- RMSE(twd_pred, Y_twd)
twd_r2   <- R2(twd_pred, Y_twd)

twd_rmse
twd_r2

### 3) VSURF + RF FOR CSHRINK ###
# Response
Y_cshrink <- dendro_clean$cshrink

# VSURF variable selection
set.seed(123)

vsurf_cshrink <- VSURF(
  x = X,
  y = Y_cshrink,
  ntree = 1000,
  parallel = TRUE
)

# Use prediction step variables
vars_cshrink <- vsurf_cshrink$varselect.pred

# Reduced predictor set 
X_cshrink_sel <- X[, vars_cshrink, drop = FALSE]

# Random forest model
set.seed(123)

rf_cshrink <- randomForest(
  x = X_cshrink_sel,
  y = Y_cshrink,
  ntree = 1000,
  importance = TRUE
)

# Predictions
cshrink_pred <- predict(rf_cshrink, X_cshrink_sel)

# Performance metrics
cshrink_rmse <- RMSE(cshrink_pred, Y_cshrink)
cshrink_r2   <- R2(cshrink_pred, Y_cshrink)

cshrink_rmse
cshrink_r2