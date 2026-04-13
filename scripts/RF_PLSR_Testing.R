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

#### 1) Try smoothing and cutting off spectra north of 850 nm ####

#### 1a) Smoothing with Savitsky-Golay (SG) ####
#install.packages("signal")
library(signal)
library(ggplot2)
library(tidyr)
library(dplyr)
## Editable parameters
# Choose spectral library: "1nm" or "5nm"
spectral_res <- "5nm"

# SG paramters
sg_window <- 11       # Must be odd (11 is a suitable starting size)
sg_poly   <- 3        # Polynomial order (must be < sg_window, 3 is a suitable starting size)

# Rows to plot (pick 1 or 2 integers from 1–53)
plot_rows <- c(1)

## Setup for SG
# Build spectral column names based on chosen resolution
if (spectral_res == "1nm") {
  spec_cols <- paste0("1nm_", 398:999)
} else if (spectral_res == "5nm") {
  spec_cols <- paste0("5nm_", seq(398, 999, by = 5))
} else {
  stop("spectral_res must be '1nm' or '5nm'")
}
# Extract spectral matrix
spectra_matrix <- as.matrix(dendro_spectra[, spec_cols])

## Smoothing with SG!
# Apply SG filter row-wise
sg_filter    <- sgolay(p = sg_poly, n = sg_window)
smoothed_matrix <- t(apply(spectra_matrix, 1, function(row) sgolayfilt(row, sg_filter)))

# Name smoothed columns
smoothed_cols <- paste0(spec_cols, "_sg")
colnames(smoothed_matrix) <- smoothed_cols

# Drop any previously appended SG columns before re-appending
sg_cols_to_remove <- grepl("_sg$", colnames(dendro_spectra))
if (any(sg_cols_to_remove)) {
  dendro_spectra <- dendro_spectra[, !sg_cols_to_remove]
}
# Append smoothed columns to original dataframe
dendro_spectra <- cbind(dendro_spectra, as.data.frame(smoothed_matrix))

## Plot smoothing results
# Extract wavelengths from column names
wavelengths <- as.numeric(
  gsub(paste0(spectral_res, "_"), "", spec_cols)
)

# Build long-format data for selected rows
plot_data <- lapply(plot_rows, function(i) {
  data.frame(
    wavelength = wavelengths,
    original   = spectra_matrix[i, ],
    smoothed   = smoothed_matrix[i, ],
    sample     = paste0("Row ", i)
  )
}) |> bind_rows()

plot_long <- plot_data |>
  pivot_longer(cols = c(original, smoothed),
               names_to  = "type",
               values_to = "reflectance")

sg_plot <- ggplot(plot_long, aes(x = wavelength, y = reflectance,
                                 color = type, linetype = type)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ sample, ncol = 1) +
  scale_color_manual(values = c(original = "steelblue", smoothed = "firebrick"),
                     labels  = c(original = "Original", smoothed = "SG Smoothed")) +
  scale_linetype_manual(values = c(original = "dashed", smoothed = "solid"),
                        labels  = c(original = "Original", smoothed = "SG Smoothed")) +
  labs(
    title    = paste0("Savitzky-Golay Smoothing (", spectral_res, " library)"),
    subtitle = paste0("Window = ", sg_window, ", Polynomial order = ", sg_poly),
    x        = "Wavelength (nm)",
    y        = "Reflectance",
    color    = NULL, linetype = NULL
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

sg_plot

ggsave(
  filename = "outputs/dendro_savitzkygolay_11_3.png",
  plot     = sg_plot,
  width    = 8,
  height   = 6,
  dpi      = 300
)

### 1c) Data prep!!!! ###
library(pls)

# Select spectral predictors - 1 nm spectral library
X <- dendro_spectra[, paste0("1nm_", 398:999)]

# remove rows with missing data (if any)
complete_idx <- complete.cases(X, dendro_spectra$TWD_drone, dendro_spectra$cshrink)

X <- X[complete_idx, ]

#### 2a) PLSR FOR TWD_DRONE - TWD ON DAY OF DRONE FLIGHT (TEST) ####
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

#### 2b) PLSR FOR CSHRINK - CONSECUTIVE DAYS PRE-FLIGHT OF NO NET GROWTH (TEST) ####
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
# Export
write.csv(results, file = "outputs/dendro_PLSR_CV_specres.csv", row.names = FALSE)

###### 4) PLSRs WITH BEST PERFORMING SPEC LIB #######
library(pls)
library(caret)
library(ggplot2)

# ---- SETTINGS + DATA PREP ----
TARGET_VAR     <- "cshrink"   # "TWD_drone" or "cshrink" or other
RUN_PERMTEST   <- TRUE         # TRUE to run permutation test, FALSE to skip
N_PERMS        <- 999           # number of permutations, 999 is good start, 10,000 is standard for publications (if RUN_PERMTEST = TRUE)
NCOMP_MAX      <- 7             # max components (sqrt of n ~ 7 for n=53)
N_FOLDS        <- 10            # k for k-fold CV
N_REPEATS_KFOLD <- 100           # number of CV repeats
N_REPEATS_PERM <- 10            # number of CV repeats for the permutation
SEED           <- 42            # 42 is the default!

## Data prep
# Select which spectral libray to go with
#X <- dendro_spectra[, paste0("1nm_", 398:999)] #(1 nm, not smoothed)
#X <- dendro_spectra[, paste0("5nm_", seq(398, 999, by = 5))] #(5 nm, not smoothed)
X <- dendro_spectra[, paste0("5nm_", seq(398, 898, by = 5), "_sg")] #(5nm, smoothed, removes > 900nm)

# Clean up selected library
complete_idx <- complete.cases(X, dendro_spectra[[TARGET_VAR]])
X <- X[complete_idx, ]
#wavelengths <- as.numeric(gsub("5nm_", "", colnames(X))) # not smoothed columns
wavelengths <- as.numeric(gsub("5nm_|_sg", "", colnames(X))) # smoothed columns
Y <- dendro_spectra[[TARGET_VAR]][complete_idx]
X_mat <- as.matrix(X)

cat("Modelling:", TARGET_VAR, "\n")
cat("n =", length(Y), "\n")

# ---- REPEATED K-FOLD CV ----
set.seed(SEED)
ctrl <- trainControl(
  method = "repeatedcv",
  number = N_FOLDS,
  repeats = N_REPEATS_KFOLD,
  savePredictions = "final"
)

pls_model <- train(
  x = X_mat,
  y = Y,
  method = "pls",
  tuneGrid = data.frame(ncomp = 1:NCOMP_MAX),
  trControl = ctrl,
  preProcess = "scale"
)

# Select optimal components
n_opt <- pls_model$bestTune$ncomp
cat("Optimal ncomp:", n_opt, "\n")

# Cross validation metrics
cv_preds <- pls_model$pred[pls_model$pred$ncomp == n_opt, ]

rmsep_cv <- sqrt(mean((cv_preds$obs - cv_preds$pred)^2))
r2_cv    <- cor(cv_preds$obs, cv_preds$pred)^2
bias_cv  <- mean(cv_preds$pred - cv_preds$obs)
rpd_cv   <- sd(Y) / rmsep_cv

cat("\n--- CV Metrics ---\n")
cat("RMSEP:", round(rmsep_cv, 4), "\n")
cat("R²:   ", round(r2_cv, 4), "\n")
cat("Bias: ", round(bias_cv, 4), "\n")
cat("RPD:  ", round(rpd_cv, 4), "\n")

# REFIT FINAL MODEL ON ALL DATA
pls_final <- plsr(
  Y ~ X_mat,
  ncomp = n_opt,
  scale = TRUE,
  validation = "none"
)

# Variable importance scores
vip_func <- function(model, ncomp) {
  W <- model$loading.weights[, 1:ncomp, drop = FALSE]
  Q <- matrix(model$Yloadings[1:ncomp], nrow = ncomp)  # force matrix for single response
  T <- model$scores[, 1:ncomp, drop = FALSE]
  
  SS <- as.vector(Q^2) * colSums(T^2)
  W_norm <- W / sqrt(colSums(W^2))
  
  vip <- sqrt(nrow(W) * rowSums(W_norm^2 %*% diag(SS)) / sum(SS))
  return(vip)
}

vip_scores <- vip_func(pls_final, ncomp = n_opt)
length(vip_scores)  # should be 121 (for full 5nm library)

# (Optional) permutation test
if (RUN_PERMTEST) {
  
  # Separate lightweight ctrl for permutations
  ctrl_perm <- trainControl(
    method = "repeatedcv",
    number = N_FOLDS,
    repeats = N_REPEATS_PERM,              
    savePredictions = "final"
  )
  
  cat("\nRunning permutation test (", N_PERMS, "permutations)...\n")
  set.seed(SEED)
  perm_r2 <- numeric(N_PERMS)
  
  for (i in 1:N_PERMS) {
    Y_perm <- sample(Y)
    pls_perm <- train(
      x = X_mat,
      y = Y_perm,
      method = "pls",
      tuneGrid = data.frame(ncomp = n_opt),
      trControl = ctrl_perm,     # <-- uses lighter ctrl, not the main one
      preProcess = "scale"
    )
    preds_perm <- pls_perm$pred
    perm_r2[i] <- cor(preds_perm$obs, preds_perm$pred)^2
  }
  
  p_val <- mean(perm_r2 >= r2_cv)
  cat("Permutation test p-value:", p_val, "\n")
  
  # Permutation null distribution plot
  hist(perm_r2,
       main = paste("Permutation Test -", TARGET_VAR),
       xlab = "Permuted R²",
       breaks = 30,
       col = "lightgrey")
  abline(v = r2_cv, col = "red", lwd = 2)
  legend("topright", legend = paste("Observed R² =", round(r2_cv, 3)),
         col = "red", lwd = 2)
} else {
  p_val <- NA
}

# ---- PLOTS ----

# Observed vs predicted
p_obs <- ggplot(cv_preds, aes(x = obs, y = pred)) +
  geom_abline(slope = 1, intercept = 0,
              color = "red", linetype = "dashed", linewidth = 0.9) +
  geom_point(size = 2.5, alpha = 0.8, shape = 21, fill = NA, color = "black") +
  annotate("text",
           x     = min(cv_preds$obs),
           y     = max(cv_preds$pred),
           label = paste0("R² = ", round(r2_cv, 3),
                          "\nRMSEP = ", round(rmsep_cv, 3),
                          "\nRPD = ", round(rpd_cv, 3)),
           hjust = 0, vjust = 1, size = 3.8) +
  labs(title = paste(TARGET_VAR, "- Observed vs Predicted (Repeated k-fold CV)"),
       x = paste("Observed", TARGET_VAR),
       y = paste("Predicted", TARGET_VAR)) +
  theme_bw(base_size = 14)
print(p_obs)

# VIP scores
vip_df <- data.frame(Wavelength = wavelengths, VIP = vip_scores)

p_vip <- ggplot(vip_df, aes(x = Wavelength, y = VIP)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 1, color = "red",
             linetype = "dashed", linewidth = 0.9) +
  labs(title = paste("VIP Scores -", TARGET_VAR),
       x = "Wavelength (nm)",
       y = "VIP Score") +
  theme_bw(base_size = 14)
print(p_vip)

# Permutation test
if (RUN_PERMTEST) {
  perm_df <- data.frame(perm_r2 = perm_r2)
  
  p_perm <- ggplot(perm_df, aes(x = perm_r2)) +
    geom_histogram(bins = 30, fill = "lightgrey", color = "white") +
    geom_vline(xintercept = r2_cv, color = "red",
               linetype = "solid", linewidth = 0.9) +
    annotate("text",
             x     = r2_cv,
             y     = Inf,
             label = paste0("Observed R² = ", round(r2_cv, 3)),
             hjust = 1.2, vjust = 4, size = 3.8, color = "red") +
    labs(title = paste("Permutation Test -", TARGET_VAR),
         x = "Permuted R²",
         y = "Frequency") +
    theme_bw(base_size = 14)
  print(p_perm)
}

# ---- RESULTS TABLE ----
results_row <- data.frame(
  Library        = "5nm",
  Predictor      = TARGET_VAR,
  nComp_optimal  = n_opt,
  RMSEP_CV       = rmsep_cv,
  R2_CV          = r2_cv,
  Bias           = bias_cv,
  RPD            = rpd_cv,
  Perm_p         = p_val
)
print(results_row)

# ---- EXPORT ----
# Model metrics
write.csv(results_row,
          file = paste0("outputs/dendro_PLSR_results_", TARGET_VAR, ".csv"),
          row.names = FALSE)

# Observed vs predicted values
write.csv(data.frame(Observed = cv_preds$obs, Predicted = cv_preds$pred),
          file = paste0("outputs/dendro_PLSR_obs_vs_pred_", TARGET_VAR, ".csv"),
          row.names = FALSE)

# VIP scores
write.csv(vip_df,
          file = paste0("outputs/dendro_PLSR_VIP_", TARGET_VAR, ".csv"),
          row.names = FALSE)

# Plots
ggsave(paste0("outputs/dendro_PLSR_obs_vs_pred_", TARGET_VAR, ".png"),
       plot = p_obs, width = 7, height = 7, dpi = 300)

ggsave(paste0("outputs/dendro_PLSR_VIP_", TARGET_VAR, ".png"),
       plot = p_vip, width = 10, height = 5, dpi = 300)

if (RUN_PERMTEST) {
  ggsave(paste0("outputs/dendro_PLSR_permtest_", TARGET_VAR, ".png"),
         plot = p_perm, width = 7, height = 6, dpi = 300)
}

#############################################################################
################ RF MODELS TO PREDICT DENDRO ATTRIBUTES ######################
# Filter out stats first before RF'ing

## READ IN SPECTRAL + STRUCTURAL DATA + CLEAN UP
lidar_metrics <- readRDS("./data/lidar_metrics.rds")
chrono <- readRDS("./data/chrono_VIstats_metrics.rds")
# Dendrometer metrics
dendro_spectra <- read.csv(
  file = "./data/dendro_spectra_joined.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

library(dplyr)
# Remove extraneous stat columns from dendro_spectra
dendro_filt <- dendro_spectra %>%
  select(-matches("(SD|Q25|Q75|Mean)$"))

# Pull age, DBH, BAI_2024 from raw chrono + calculate age_2
chrono_dendro <- chrono %>%
  select(TreeID, age, DBH, BAI_2024) %>%
  mutate(age_2 = 2024 - age)

# Build RF dataframe
model_df <- dendro_filt %>%
  left_join(chrono_dendro,    by = "TreeID") %>%
  left_join(lidar_metrics,  by = "TreeID")

########## 1) RANDOM FOREST WITH STRUCTURAL + SPECTRAL INDICES ##############
#install.packages("ranger")
library(ranger)
library(caret)
library(ggplot2)
library(dplyr)

# ---- SETTINGS ----
TARGET_VAR      <- "TWD_drone"   # "TWD_drone" or "cshrink"
RUN_PERMTEST    <- FALSE        # TRUE to run permutation test, FALSE to skip
N_PERMS         <- 999          # number of permutations (if RUN_PERMTEST = TRUE)
N_FOLDS         <- 10           # k for k-fold CV
N_REPEATS_KFOLD <- 100          # number of CV repeats
N_REPEATS_PERM  <- 10           # number of CV repeats for permutation
NUM_TREES       <- 500          # number of trees in the forest
MTRY            <- NULL         # predictors sampled per split (NULL = default: p/3)
TOP_N           <- 20           # number of top variables to show in importance plot
SEED            <- 42

# ---- DATA PREP ----

# Select predictor columns
rf_predictors <- model_df %>%
  select(rugosity:Area_m2, DBH, ARI1_5nm_Median:Vogelmann4_5nm_Median) %>%
  select(-any_of("zentropy")) # too many N/As

X_rf <- rf_predictors
Y_rf <- model_df[[TARGET_VAR]]

# Remove incomplete cases
complete_idx <- complete.cases(X_rf, Y_rf)
X_rf <- X_rf[complete_idx, ]
Y_rf <- Y_rf[complete_idx]

cat("Modelling:", TARGET_VAR, "\n")
cat("n =", length(Y_rf), "\n")
cat("Predictors:", ncol(X_rf), "\n")

# ---- REPEATED K-FOLD CV ----
set.seed(SEED)

mtry_val <- if (is.null(MTRY)) floor(ncol(X_rf) / 3) else MTRY

ctrl <- trainControl(
  method          = "repeatedcv",
  number          = N_FOLDS,
  repeats         = N_REPEATS_KFOLD,
  savePredictions = "final"
)

rf_model <- train(
  x         = X_rf,
  y         = Y_rf,
  method    = "ranger",
  num.trees = NUM_TREES,
  tuneGrid  = data.frame(
    mtry              = mtry_val,
    splitrule         = "variance",
    min.node.size     = 5
  ),
  trControl  = ctrl,
  importance = "impurity"
)

# CV predictions
cv_preds <- rf_model$pred

rmse_cv <- sqrt(mean((cv_preds$obs - cv_preds$pred)^2))
r2_cv   <- cor(cv_preds$obs, cv_preds$pred)^2
bias_cv <- mean(cv_preds$pred - cv_preds$obs)
rpd_cv  <- sd(Y_rf) / rmse_cv

cat("\n--- CV Metrics ---\n")
cat("RMSE:", round(rmse_cv, 4), "\n")
cat("R²:  ", round(r2_cv,   4), "\n")
cat("Bias:", round(bias_cv,  4), "\n")
cat("RPD: ", round(rpd_cv,   4), "\n")

# ---- REFIT FINAL MODEL ON ALL DATA ----
set.seed(SEED)
rf_final <- ranger(
  formula       = as.formula(paste(TARGET_VAR, "~ .")),
  data          = cbind(X_rf, setNames(data.frame(Y_rf), TARGET_VAR)),
  num.trees     = NUM_TREES,
  mtry          = mtry_val,
  importance    = "impurity",
  min.node.size = 5
)

# Variable importance
importance_df <- data.frame(
  Variable   = names(rf_final$variable.importance),
  Importance = rf_final$variable.importance
) %>%
  arrange(desc(Importance)) %>%
  slice_head(n = TOP_N)

# ---- OPTIONAL PERMUTATION TEST ----
if (RUN_PERMTEST) {
  
  ctrl_perm <- trainControl(
    method          = "repeatedcv",
    number          = N_FOLDS,
    repeats         = N_REPEATS_PERM,
    savePredictions = "final"
  )
  
  cat("\nRunning permutation test (", N_PERMS, "permutations)...\n")
  set.seed(SEED)
  perm_r2 <- numeric(N_PERMS)
  
  for (i in 1:N_PERMS) {
    Y_perm   <- sample(Y_rf)
    rf_perm  <- train(
      x         = X_rf,
      y         = Y_perm,
      method    = "ranger",
      num.trees = NUM_TREES,
      tuneGrid  = data.frame(
        mtry          = mtry_val,
        splitrule     = "variance",
        min.node.size = 5
      ),
      trControl  = ctrl_perm,
      importance = "none"   # skip importance during perms for speed
    )
    preds_perm  <- rf_perm$pred
    perm_r2[i] <- cor(preds_perm$obs, preds_perm$pred)^2
  }
  
  p_val <- mean(perm_r2 >= r2_cv)
  cat("Permutation test p-value:", p_val, "\n")
  
} else {
  p_val <- NA
}

# ---- PLOTS ----

# Observed vs predicted
p_obs <- ggplot(cv_preds, aes(x = obs, y = pred)) +
  geom_abline(slope = 1, intercept = 0,
              color = "red", linetype = "dashed", linewidth = 0.9) +
  geom_point(size = 2.5, alpha = 0.8, shape = 21,
             fill = NA, color = "black") +
  annotate("text",
           x     = min(cv_preds$obs),
           y     = max(cv_preds$pred),
           label = paste0("R² = ",   round(r2_cv,   3),
                          "\nRMSE = ", round(rmse_cv, 3),
                          "\nRPD = ",  round(rpd_cv,  3)),
           hjust = 0, vjust = 1, size = 3.8) +
  labs(title = paste(TARGET_VAR, "- Observed vs Predicted (Repeated k-fold CV)"),
       x     = paste("Observed", TARGET_VAR),
       y     = paste("Predicted", TARGET_VAR)) +
  theme_bw(base_size = 14)
print(p_obs)

# Variable importance bar chart
p_imp <- ggplot(importance_df,
                aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title    = paste("Top", TOP_N, "Variable Importance -", TARGET_VAR),
       x        = NULL,
       y        = "Importance (Impurity)") +
  theme_bw(base_size = 14)
print(p_imp)

# Permutation test plot
if (RUN_PERMTEST) {
  perm_df <- data.frame(perm_r2 = perm_r2)
  
  p_perm <- ggplot(perm_df, aes(x = perm_r2)) +
    geom_histogram(bins = 30, fill = "lightgrey", color = "white") +
    geom_vline(xintercept = r2_cv, color = "red",
               linetype = "solid", linewidth = 0.9) +
    annotate("text",
             x     = r2_cv,
             y     = Inf,
             label = paste0("Observed R² = ", round(r2_cv, 3)),
             hjust = 1.2, vjust = 4, size = 3.8, color = "red") +
    labs(title = paste("Permutation Test -", TARGET_VAR),
         x     = "Permuted R²",
         y     = "Frequency") +
    theme_bw(base_size = 14)
  print(p_perm)
}

# ---- RESULTS TABLE ----
results_row <- data.frame(
  Model      = "RandomForest",
  Predictor  = TARGET_VAR,
  RMSE_CV    = rmse_cv,
  R2_CV      = r2_cv,
  Bias       = bias_cv,
  RPD        = rpd_cv,
  Perm_p     = p_val
)
print(results_row)

# ---- EXPORT ----
write.csv(results_row,
          file      = paste0("outputs/dendro_RF_results_", TARGET_VAR, ".csv"),
          row.names = FALSE)

write.csv(data.frame(Observed = cv_preds$obs, Predicted = cv_preds$pred),
          file      = paste0("outputs/dendro_RF_obs_vs_pred_", TARGET_VAR, ".csv"),
          row.names = FALSE)

write.csv(importance_df,
          file      = paste0("outputs/dendro_RF_importance_", TARGET_VAR, ".csv"),
          row.names = FALSE)

ggsave(paste0("outputs/dendro_RF_obs_vs_pred_", TARGET_VAR, ".png"),
       plot = p_obs, width = 7, height = 7, dpi = 300)

ggsave(paste0("outputs/dendro_RF_importance_", TARGET_VAR, ".png"),
       plot = p_imp, width = 8, height = 6, dpi = 300)

if (RUN_PERMTEST) {
  ggsave(paste0("outputs/dendro_RF_permtest_", TARGET_VAR, ".png"),
         plot = p_perm, width = 7, height = 6, dpi = 300)
}

######### 1b) RF BOOSTRAP LOOP FOR TWD_DRONE (TEST) #################
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