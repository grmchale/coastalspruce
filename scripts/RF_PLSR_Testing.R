# Set working directory to coastalspruce repo

## Read in spectra + dendro data
# These are the resampled spectral libraries, VIs, and dendro metrics for each spruce crown
dendro_spectra <- read.csv(
  file = "./data/dendro_spectra_joined.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

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
  ncomp = 20,
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
  ncomp = 20,
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

### PLSR GRAPH FOR CSHRINK ###
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

################ RF MODELS TO PREDICT DENDRO ATTRIBUTES ######################
# Filter out stats first before RF'ing
dendro_clean <- dendro_spectra[
  , !grepl("(Q75|Q25|Mean|SD)$", names(dendro_spectra))
]
#install.packages("VSURF")
library(VSURF)
library(randomForest)
library(caret)
# Create predictor matrix
X <- dendro_clean[
                  which(names(dendro_clean) == "1nm_398"):
                    which(names(dendro_clean) == "Vogelmann4_15nm_Median")]

X <- as.data.frame(X)

#### 1) VSURF + RF FOR TWD_DRONE ###
# Response
Y_twd <- dendro_clean$TWD_drone

# VSURF VARIABLE SELECTION
set.seed(123)
vsurf_twd <- VSURF(
  x = X,
  y = Y_twd,
  ntree = 1000,
  parallel = TRUE
)

# Use prediction step variables
vars_twd <- vsurf_twd$varselect.pred

# Reduced predictor set
X_twd_sel <- X[, vars_twd, drop = FALSE]

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

### 2) VSURF + RF FOR CSHRINK ###
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