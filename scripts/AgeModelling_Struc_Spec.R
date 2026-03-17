########## DATA READ IN ########

crown_rugosity <-read.csv("./LiDAR/Normalized/LiDAR_crowns/crown_rugosity.csv") # Rugosity of crowns
crown_metrics <-read.csv("./LiDAR/Normalized/LiDAR_crowns/crown_metrics.csv") # Point cloud metrics
crown_areas <- read.csv("./LiDAR/Normalized/LiDAR_crowns/crown_areas.csv") # Area of each crown

lidar_metrics <- crown_rugosity |>
  left_join(crown_metrics, by = "TreeID") |>
  left_join(crown_areas, by = "TreeID")

# Write to RDS
saveRDS(lidar_metrics, "./data/lidar_metrics.rds")

# Read back in
lidar_metrics <- readRDS("./data/lidar_metrics.rds")

################## RANDOM FOREST TO PREDICT AGE #################################
install.packages("dplyr")
library(dplyr)
chrono <- readRDS("./data/chrono_VIstats_metrics.rds")
chrono_filt <- chrono %>%
  select(-matches("(SD|Q25|Q75|Mean)$"))
chrono_filt <- chrono_filt %>%
  left_join(
    dplyr::select(crown_rugosity, TreeID, rugosity),
    by = "TreeID"
  )
# safest: fully qualify dplyr verbs to avoid select() conflicts
chrono_rf <- dplyr::left_join(
  chrono_filt,
  crown_metrics,                 # do NOT drop TreeID here
  by = "TreeID"
)

# Random Forest for AGE (robust, repeated)
# Inputs: chrono_rf, spectral cols 7:790 + lidar cols "rugosity":"zpcum9"
# Output: ranked predictors by mean %IncMSE + stability


library(randomForest)

## 1) Subset to rows with AGE
df <- chrono_rf
if (!"age" %in% names(df)) stop("Column 'age' not found in chrono_rf.")
df <- df[!is.na(df$age), , drop = FALSE]   # expect ~216 rows

## 2) Build predictor matrix: spectral 7:790 and lidar rugosity:zpcum9
# spectral
pmin <- 7
pmax <- min(790, ncol(df))
if (pmin > ncol(df)) stop("Data has fewer than 7 columns; check chrono_rf.")
spec_idx <- pmin:pmax

# lidar range by names (inclusive)
lidar_start <- match("rugosity", names(df))
lidar_end   <- match("zpcum9",  names(df))
if (is.na(lidar_start) || is.na(lidar_end)) {
  stop("Could not find 'rugosity' and/or 'zpcum9' in chrono_rf column names.")
}
lidar_idx <- seq(lidar_start, lidar_end)

# combine and keep unique (in case of overlap)
pred_idx <- unique(c(spec_idx, lidar_idx))

X_raw <- df[, pred_idx, drop = FALSE]
y     <- df[["age"]]

## 3) Keep numeric predictors only
is_num <- vapply(X_raw, is.numeric, logical(1))
X_raw  <- X_raw[, is_num, drop = FALSE]

## 4) Handle missing values and constants
# drop columns with >20% NA
na_frac <- vapply(X_raw, function(z) mean(is.na(z)), numeric(1))
X_raw   <- X_raw[, na_frac <= 0.20, drop = FALSE]

# median-impute remaining NA per column
for (j in seq_len(ncol(X_raw))) {
  z <- X_raw[[j]]
  if (anyNA(z)) {
    med <- median(z, na.rm = TRUE)
    z[is.na(z)] <- med
    X_raw[[j]] <- z
  }
}

# drop zero-variance columns
sd_vec <- vapply(X_raw, function(z) sd(z, na.rm = TRUE), numeric(1))
X_raw  <- X_raw[, sd_vec > 0, drop = FALSE]

# Optional: standardize (RF doesn't require, harmless)
X <- as.data.frame(scale(X_raw))

## 5) Repeated RF for robust importance
set.seed(41)
B <- 200                         # number of runs; raise to 500 if you want extra stability
p <- ncol(X)
mtry_val <- max(1, floor(sqrt(p)))
ntree_val <- 1000

# storage
imp_mat  <- matrix(NA_real_, nrow = p, ncol = B, dimnames = list(colnames(X), NULL))
oob_mse  <- numeric(B)

# track top-K frequency (stability)
K <- 20
top_hits <- integer(p); names(top_hits) <- colnames(X)

for (b in seq_len(B)) {
  set.seed(41 + b)
  rf_fit <- randomForest(
    x = X, y = y,
    importance = TRUE,
    mtry = mtry_val,
    ntree = ntree_val,
    nodesize = 5
  )
  # store OOB MSE for this run
  oob_mse[b] <- rf_fit$mse[ntree_val]
  
  # permutation importance
  imp <- importance(rf_fit, type = 1)[, 1]   # %IncMSE; named by predictor
  imp[is.na(imp)] <- 0
  imp_mat[, b] <- imp[colnames(X)]          # align by name
  
  # stability: count appearance in top-K for this run
  ord <- order(imp, decreasing = TRUE)
  top_vars <- names(imp)[ord][seq_len(min(K, length(ord)))]
  top_hits[top_vars] <- top_hits[top_vars] + 1L
}

## 6) Aggregate importance + stability
mean_incMSE <- rowMeans(imp_mat, na.rm = TRUE)
sd_incMSE   <- apply(imp_mat, 1, sd, na.rm = TRUE)
rank_rf     <- rank(-mean_incMSE, ties.method = "average")
topK_freq   <- top_hits / B

rf_age_summary <- data.frame(
  Predictor   = names(mean_incMSE),
  MeanIncMSE  = as.numeric(mean_incMSE),
  SDIncMSE    = as.numeric(sd_incMSE),
  RF_Rank     = as.numeric(rank_rf),
  TopK_Freq   = as.numeric(topK_freq),
  row.names   = NULL
)

# sort by rank (lower is better)
rf_age_summary <- rf_age_summary[order(rf_age_summary$RF_Rank), ]

## 7) Report: overall robustness + strongest predictors
cat(sprintf("\nRepeated RF runs: B = %d, ntree = %d, mtry = %d, predictors p = %d, rows n = %d\n",
            B, ntree_val, mtry_val, p, nrow(X)))
cat(sprintf("OOB RMSE (mean ± sd): %.2f ± %.2f\n\n",
            sqrt(mean(oob_mse)), sd(sqrt(oob_mse))))

# Show the top 25 predictors (adjust as you like)
head(rf_age_summary, 25)

# assuming you already have: X, y, B, mtry_val, ntree_val, nodesize=5
set.seed(41)
oob_rmse <- numeric(B)
oob_mae  <- numeric(B)
oob_r2   <- numeric(B)

for (b in seq_len(B)) {
  set.seed(41 + b)
  rf_fit <- randomForest(
    x = X, y = y,
    importance = TRUE,
    mtry = mtry_val,
    ntree = ntree_val,
    nodesize = 5
  )
  
  # OOB predictions are in rf_fit$predicted for regression
  pred_oob <- rf_fit$predicted
  resid    <- y - pred_oob
  
  oob_rmse[b] <- sqrt(mean(resid^2, na.rm = TRUE))
  oob_mae[b]  <- mean(abs(resid), na.rm = TRUE)
  oob_r2[b]   <- 1 - mean(resid^2, na.rm = TRUE) / var(y, na.rm = TRUE)
}

cat(sprintf("OOB RMSE (mean ± sd): %.2f ± %.2f\n", mean(oob_rmse), sd(oob_rmse)))
cat(sprintf("OOB  MAE (mean ± sd): %.2f ± %.2f\n", mean(oob_mae),  sd(oob_mae)))
cat(sprintf("OOB   R² (mean ± sd): %.3f ± %.3f\n\n", mean(oob_r2),  sd(oob_r2)))

# Write top predictors to disk
write.csv(rf_age_summary, "./R_outputs/modelling/rf/age_prediction_rf/age_strongestpredictors.csv")

################### LINEAR MIXED EFFECTS MODEL FOR AGE ##########################
library(lme4)
library(dplyr)

# 1) Filter usable rows
df_lme <- chrono_rf %>%
  filter(!is.na(age),
         !is.na(zmax),
         !is.na(Datt3_1nm_Median),
         !is.na(rugosity),
         !is.na(Site))

# 2) Scale predictors
df_lme <- df_lme %>%
  mutate(
    zmax_sc      = scale(zmax),
    Datt3_sc     = scale(Datt3_1nm_Median),
    rugosity_sc  = scale(rugosity)
  )

# 3) Fit mixed model (random intercept by Site)
mod_age <- lmer(
  age ~ zmax_sc + Datt3_sc + rugosity_sc + (1 | Site),
  data = df_lme,
  REML = FALSE
)

# 4) Summarize model
summary(mod_age)

# Test model
# --- 1. Marginal vs Conditional R² ---
# Given: mod_age is your lmer model:
# mod_age <- lmer(age ~ zmax_sc + Datt3_sc + rugosity_sc + (1 | Site), data = df_lme, REML = FALSE)

# Fixed-effects variance (variance of predictions with random effects turned off)
yhat_fix <- predict(mod_age, re.form = NA)   # fixed-part only
var_fix  <- var(yhat_fix, na.rm = TRUE)

# Random-effects variance components (sum all G matrix variances)
vc <- as.data.frame(VarCorr(mod_age))
var_rand <- sum(vc$vcov, na.rm = TRUE)

# Residual (sigma^2)
var_res <- sigma(mod_age)^2

# R2 calculations
R2_marginal   <- var_fix / (var_fix + var_rand + var_res)
R2_conditional<- (var_fix + var_rand) / (var_fix + var_rand + var_res)

cat(sprintf("Marginal R² (fixed only):     %.3f\n", R2_marginal))
cat(sprintf("Conditional R² (fixed+random): %.3f\n", R2_conditional))


# --- 2. Diagnostics ---
par(mfrow = c(1, 2))
plot(fitted(mod_age), resid(mod_age),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2, col = "red")

qqnorm(resid(mod_age), main = "Normal Q-Q")
qqline(resid(mod_age), col = "red")
par(mfrow = c(1, 1))

# --- 3. Leave-one-site-out CV ---
sites <- unique(df_lme$Site)
cv_results <- data.frame(Site = sites, RMSE = NA, R2 = NA)

for (s in sites) {
  train <- df_lme[df_lme$Site != s, ]
  test  <- df_lme[df_lme$Site == s, ]
  
  fit <- lmer(age ~ zmax_sc + Datt3_sc + rugosity_sc + (1 | Site),
              data = train, REML = FALSE)
  
  preds <- predict(fit, newdata = test, allow.new.levels = TRUE)
  rmse  <- sqrt(mean((test$age - preds)^2))
  r2    <- cor(test$age, preds)^2
  
  cv_results[cv_results$Site == s, c("RMSE","R2")] <- c(rmse, r2)
}

print(cv_results)

############## AGE vs ZMAX LINEAR REGRESSION ##########
# Keep only rows with both age and zmax
df_lin <- chrono_rf[!is.na(chrono_rf$age) & !is.na(chrono_rf$zmax), ]

# Fit linear regression
lm_age_zmax <- lm(age ~ zmax, data = df_lin)

# Model summary
summary(lm_age_zmax)

# Quick plot with regression line
plot(df_lin$zmax, df_lin$age,
     xlab = "zmax", ylab = "Age",
     main = "Linear regression: Age ~ zmax",
     pch = 16, col = "darkblue")
abline(lm_age_zmax, col = "red", lwd = 2)