########## DATA READ IN + CLEANING ########

#crown_rugosity <-read.csv("./LiDAR/Normalized/LiDAR_crowns/crown_rugosity.csv") # Rugosity of crowns
#crown_metrics <-read.csv("./LiDAR/Normalized/LiDAR_crowns/crown_metrics.csv") # Point cloud metrics
#crown_areas <- read.csv("./LiDAR/Normalized/LiDAR_crowns/crown_areas.csv") # Area of each crown

#lidar_metrics <- crown_rugosity |>
#  left_join(crown_metrics, by = "TreeID") |>
#  left_join(crown_areas, by = "TreeID")

# Write to RDS
#saveRDS(lidar_metrics, "./data/lidar_metrics.rds")

# Read back in
lidar_metrics <- readRDS("./data/lidar_metrics.rds")
chrono <- readRDS("./data/chrono_VIstats_metrics.rds")
#install.packages("dplyr")
library(dplyr)
# Filter out extraneous stats
chrono_filt <- chrono %>%
  select(-matches("(SD|Q25|Q75|Mean)$"))
# Join spectra + struc metrics together
struc_spec <- chrono_filt |> inner_join(lidar_metrics, by = "TreeID")
# Create new field representing tree age in years (rather than germination year)
struc_spec <- struc_spec |> mutate(age_2 = 2024 - age)
# Remove crowns with no age
struc_spec <- struc_spec |> filter(!is.na(age))

################## RANDOM FOREST TO PREDICT AGE #################################

# Random Forest for AGE (robust, repeated)
# Inputs: struc_spec, spectral cols 7:790 + lidar cols "rugosity":"zpcum9"
# Output: ranked predictors by mean %IncMSE + stability

library(randomForest)


## 1) Build predictor matrix: spectral 7:790 and lidar rugosity:zpcum9
# spectral
df <- struc_spec
pmin <- 7
pmax <- min(398, ncol(df))
if (pmin > ncol(df)) stop("Data has fewer than 7 columns; check chrono_rf.")
spec_idx <- pmin:pmax

# lidar range by names (inclusive)
lidar_start <- match("rugosity", names(df))
lidar_end   <- match("Area_m2",  names(df))
if (is.na(lidar_start) || is.na(lidar_end)) {
  stop("Could not find 'rugosity' and/or 'Area_m2' in chrono_rf column names.")
}
lidar_idx <- seq(lidar_start, lidar_end)

# combine and keep unique (in case of overlap)
pred_idx <- unique(c(spec_idx, lidar_idx))

X_raw <- df[, pred_idx, drop = FALSE]
y     <- df[["age_2"]]

## 2) Keep numeric predictors only
is_num <- vapply(X_raw, is.numeric, logical(1))
X_raw  <- X_raw[, is_num, drop = FALSE]

## 3) Handle missing values and constants
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

## 4) Repeated RF for robust importance
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

## 5) Aggregate importance + stability
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

## 6) Report: overall robustness + strongest predictors
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
df_lin <- struc_spec[!is.na(struc_spec$age_2) & !is.na(struc_spec$zmax), ]

# Fit linear regression
lm_age_zmax <- lm(age ~ zmax, data = df_lin)

# Model summary
summary(lm_age_zmax)

# Quick plot with regression line
plot(df_lin$zmax, df_lin$age_2,
     xlab = "zmax", ylab = "Age",
     main = "Linear regression: Age ~ zmax",
     pch = 16, col = "darkblue")
abline(lm_age_zmax, col = "red", lwd = 2)

########## QUANTILE REGRESSIONS ##################
#install.packages("quantreg")
#install.packages("quantregForest")
#install.packages("ggplot2")
library(quantreg)

#### (TEST) Single predictor (one veg index) ###
qr_single <- rq(age_2 ~ EGFR_5nm_Median, data = struc_spec, tau = c(0.10, 0.25, 0.50, 0.75, 0.90))
summary(qr_single)

# Pseudo R squared
rq.fit <- rq(age_2 ~ EGFR_5nm_Median, data = struc_spec, tau = 0.50)
rq.null <- rq(age_2 ~ 1, data = struc_spec, tau = 0.50)

# Pseudo R² per quantile
1 - (summary(rq.fit)$rho / summary(rq.null)$rho)

taus <- c(0.10, 0.25, 0.50, 0.75, 0.90)

pseudo_r2 <- sapply(taus, function(tau) {
  fit  <- rq(age_2 ~ EGFR_5nm_Median, data = struc_spec, tau = tau)
  null <- rq(age_2 ~ 1,               data = struc_spec, tau = tau)
  1 - (summary(fit)$rho / summary(null)$rho)
})

data.frame(tau = taus, pseudo_R2 = round(pseudo_r2, 3))

#### (TEST) Multiple predictors (spectral + structural) ###
qr_multi <- rq(age_2 ~ EGFR_5nm_Median + zmax + zq95 + zsd + rugosity,
               data = struc_spec,
               tau = c(0.25, 0.50, 0.75, 0.90))
summary(qr_multi)

#### (TEST) Quantile regression forest ###
library(quantregForest)
library(ggplot2)

# Set up your predictor matrix — use your top predictors from the RF importance output
predictors <- c("zmax", "zq95", "zsd", "zmean", "rugosity",
                "Datt3_1nm_Median", "D1_5nm_Median", "Boochs_1nm_Median", "EGFN_5nm_Median")
X_qrf <- struc_spec[, predictors]
y_qrf <- struc_spec$age_2
# Fit quantile regression forest
qrf_model <- quantregForest(x = X_qrf, y = y_qrf, ntree = 1000, keep.inbag = TRUE)
# Predict quantiles for each observation
qrf_preds <- predict(qrf_model, newdata = X_qrf, what = c(0.10, 0.50, 0.90))
colnames(qrf_preds) <- c("q10", "q50", "q90")

# ── 1. OOB-based RMSE on the median (q50) prediction ──────────────────────────
# This is comparable to your earlier RF RMSE
oob_preds_median <- predict(qrf_model, what = 0.50)  # uses OOB by default
rmse_qrf <- sqrt(mean((y_qrf - oob_preds_median)^2, na.rm = TRUE))
mae_qrf  <- mean(abs(y_qrf - oob_preds_median), na.rm = TRUE)
cat(sprintf("QRF Median RMSE: %.2f\n", rmse_qrf))
cat(sprintf("QRF Median MAE:  %.2f\n", mae_qrf))

# ── 2. Pinball loss (quantile-specific accuracy) ───────────────────────────────
# This is the proper loss function for quantile predictions.
# For a given tau, it penalizes over- and under-prediction asymmetrically.
pinball <- function(y, q_hat, tau) {
  mean(ifelse(y >= q_hat, tau * (y - q_hat), (tau - 1) * (y - q_hat)))
}

oob_q10 <- predict(qrf_model, what = 0.10)
oob_q50 <- predict(qrf_model, what = 0.50)
oob_q90 <- predict(qrf_model, what = 0.90)

cat(sprintf("Pinball loss q10: %.3f\n", pinball(y_qrf, oob_q10, 0.10)))
cat(sprintf("Pinball loss q50: %.3f\n", pinball(y_qrf, oob_q50, 0.50)))
cat(sprintf("Pinball loss q90: %.3f\n", pinball(y_qrf, oob_q90, 0.90)))

# ── 3. Coverage: do 80% of obs fall inside the 10th–90th interval? ────────────
# Ideally ~80% of actual ages should fall between q10 and q90 predictions
coverage_80 <- mean(y_qrf >= oob_q10 & y_qrf <= oob_q90, na.rm = TRUE)
cat(sprintf("80%% interval coverage: %.1f%%\n", coverage_80 * 100))
# If this is >> 80%, your intervals are too wide (overconfident uncertainty)
# If this is << 80%, your intervals are too narrow (underconfident)

# ── 4. Interval width ─────────────────────────────────────────────────────────
# How wide are the prediction intervals on average?
interval_width <- oob_q90 - oob_q10
cat(sprintf("Mean 80%% interval width: %.1f years\n", mean(interval_width)))

# ── 5. Observed vs predicted plot (median) ────────────────────────────────────
diag_df <- data.frame(
  observed  = y_qrf,
  predicted = as.numeric(oob_q50),
  q10       = as.numeric(oob_q10),
  q90       = as.numeric(oob_q90)
)

ggplot(diag_df, aes(x = observed, y = predicted)) +
  geom_errorbar(aes(ymin = q10, ymax = q90), alpha = 0.3, color = "steelblue") +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Observed Age", y = "Predicted Age (median + 80% interval)",
       title = "QRF: Observed vs Predicted") +
  theme_bw()

# ── 6. Residual plot ──────────────────────────────────────────────────────────
diag_df$residual <- diag_df$observed - diag_df$predicted

ggplot(diag_df, aes(x = predicted, y = residual)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "loess", se = FALSE, color = "steelblue") +
  labs(x = "Predicted Age", y = "Residual (Observed - Predicted)",
       title = "QRF Residuals") +
  theme_bw()

#### Continuous binning model + model metrics + plot ###

library(quantreg)
library(ggplot2)
# Binning approach analog — age predicts EGFR
qr_continuous <- rq(EGFR_5nm_Median ~ age_2, 
                    data = struc_spec,
                    tau  = c(0.10, 0.25, 0.50, 0.75, 0.90))

taus <- c(0.10, 0.25, 0.50, 0.75, 0.90)

qr_fit <- rq(EGFR_5nm_Median ~ age_2, data = struc_spec, tau = taus)

# Build prediction dataframe across the full age range
age_seq  <- data.frame(age_2 = seq(min(struc_spec$age_2), 
                                   max(struc_spec$age_2), 
                                   length.out = 200))

# Get predicted EGFR at each tau across the age sequence
pred_list <- lapply(taus, function(tau) {
  m <- rq(EGFR_5nm_Median ~ age_2, data = struc_spec, tau = tau)
  data.frame(age_2 = age_seq$age_2,
             EGFR_pred = predict(m, newdata = age_seq),
             tau = as.factor(tau))
})

pred_df <- do.call(rbind, pred_list)

## Extract model metrics from each tau
taus <- c(0.10, 0.25, 0.50, 0.75, 0.90)

# Fit each tau separately so we can extract metrics
models_by_tau <- lapply(taus, function(tau) {
  rq(EGFR_5nm_Median ~ age_2, data = struc_spec, tau = tau)
})
names(models_by_tau) <- paste0("tau_", taus)

# Null model per tau (intercept only)
null_by_tau <- lapply(taus, function(tau) {
  rq(EGFR_5nm_Median ~ 1, data = struc_spec, tau = tau)
})

# Extract AIC and R1 for each tau
metrics_df <- data.frame(
  tau      = taus,
  AIC      = sapply(models_by_tau, AIC),
  R1       = mapply(function(m, n) 1 - m$rho / n$rho, models_by_tau, null_by_tau),
  coef_int = sapply(models_by_tau, function(m) coef(m)[1]),  # intercept
  coef_age = sapply(models_by_tau, function(m) coef(m)[2])   # slope on age_2
)

metrics_df <- metrics_df |>
  mutate(across(c(AIC, R1, coef_int, coef_age), ~ round(.x, 3)))

print(metrics_df)
write.csv(metrics_df, "outputs/EGFR_age_metricsbytau.csv", row.names = FALSE)

# Plot
p <- ggplot() +
  geom_point(data = struc_spec, aes(x = age_2, y = EGFR_5nm_Median),
             alpha = 0.4, size = 1.5) +
  geom_line(data = pred_df, aes(x = age_2, y = EGFR_pred, 
                                color = tau, linetype = tau), 
            linewidth = 0.9) +
  scale_color_manual(values = c("0.1"  = "#2166ac",
                                "0.25" = "#74add1",
                                "0.5"  = "#000000",
                                "0.75" = "#f46d43",
                                "0.9"  = "#d73027"),
                     name = "Quantile") +
  scale_linetype_manual(values = c("0.1"  = "dashed",
                                   "0.25" = "dotdash",
                                   "0.5"  = "solid",
                                   "0.75" = "dotdash",
                                   "0.9"  = "dashed"),
                        name = "Quantile") +
  labs(x     = "Tree Age",
       y     = "EGFR 5nm Median",
       title = "Quantile Regression: Tree Age vs. EGFR") +
  theme_bw(base_size = 14)

ggsave("outputs/EGFR_age_quantreg.png", plot = p,
       width = 8, height = 6, units = "in", dpi = 300)

#### LICHEN 90TH QUANTILE APPROACH + TESTS MULTIPLE MODELS ####
library(quantreg)
library(tidyverse)

# ── 1. Define candidate models 
XVAR <- c(
  "EGFR_5nm_Median",
  "EGFR_5nm_Median + zq95",
  "EGFR_5nm_Median + zq95 + zsd + rugosity",
  "EGFR_5nm_Median + zq95 + zsd",
  "zmax + zsd + rugosity",
  "zq95",
  "zmax",
  "rugosity",
  "Datt3_1nm_Median",
  "Datt3_1nm_Median + zq95 + zsd"
)

XVAR_names <- c(
  "EGFR_only",
  "EGFR_zq95",
  "EGFR_full_structural",
  "EGFR_zq95_zsd",
  "structural_only",
  "zq95_only",
  "zmax_only",
  "rugosity_only",
  "Datt3_only",
  "Datt3_zq95_zsd"
)

# ── 2. Fit all models at tau = 0.90 ───────────────────────────────────────────
formulas <- lapply(paste("age_2 ~", XVAR), as.formula)
models   <- lapply(formulas, function(f) rq(f, tau = 0.90, data = struc_spec))
names(models) <- XVAR_names

# ── 3. Null model for pseudo-R² ───────────────────────────────────────────────
null_model <- rq(age_2 ~ 1, tau = 0.90, data = struc_spec)

# ── 4. Pseudo-R² (R1) and AIC for each model ──────────────────────────────────
pseudo_r2 <- sapply(models, function(m) 1 - m$rho / null_model$rho)
aic_vals  <- sapply(models, AIC)

model_comparison <- data.frame(
  Model     = XVAR_names,
  Formula   = XVAR,
  AIC       = round(aic_vals, 2),
  PseudoR2  = round(pseudo_r2, 3)
) |> arrange(AIC)

print(model_comparison)
write.csv(model_comparison, "outputs/age_quantreg_modelcomparison.csv", row.names = FALSE)

# ── 5. Bootstrapped confidence intervals for best model ───────────────────────
# Replace "EGFR_full_structural" with whichever model wins above
best_model_name <- model_comparison$Model[1]
best_model      <- models[[best_model_name]]

pred_ci <- predict(best_model, struc_spec,
                   interval  = "confidence",
                   level     = 0.95,
                   type      = "percentile",
                   se        = "boot",
                   R         = 1000)  # bootstrap resamples

pred_df_2 <- cbind(struc_spec, as.data.frame(pred_ci))
# pred_df now has columns: fit, lower, higher

# ── 6. Figure: age ~ EGFR with 90th quantile fit + CI ─────────────────────────
tiff(paste0("outputs/EGFR_age_90th_quantile_", format(Sys.time(), "%y%d%m"), ".tiff"),
     width = 1200, height = 1200, units = "px", pointsize = 24)

plot(age_2 ~ EGFR_5nm_Median, data = struc_spec,
     xlab = "EGFR 5nm Median",
     ylab = "Tree Age",
     cex.lab = 1.7, cex.axis = 1.7)

# Sort by predictor for clean line plotting
pred_sorted <- pred_df_2[order(pred_df_2$EGFR_5nm_Median), ]
points(pred_sorted$EGFR_5nm_Median, pred_sorted$fit,    col = "blue",  pch = 19)
points(pred_sorted$EGFR_5nm_Median, pred_sorted$lower,  col = "black", pch = 4, cex = 0.25)
points(pred_sorted$EGFR_5nm_Median, pred_sorted$higher, col = "black", pch = 4, cex = 0.25)

# Add pseudo-R² and model formula to plot
best_r2 <- pseudo_r2[best_model_name]
text(max(struc_spec$EGFR_5nm_Median) * 0.6,
     max(struc_spec$age_2, na.rm = TRUE) * 0.95,
     labels = paste0("R1 = ", round(best_r2, 3)), cex = 1.4)
text(max(struc_spec$EGFR_5nm_Median) * 0.6,
     max(struc_spec$age_2, na.rm = TRUE) * 0.88,
     labels = model_comparison$Formula[1], cex = 0.9)

legend("topright",
       legend = c("Raw Data", "90th Quantile Fit", "95% Bootstrap CI"),
       col    = c("black", "blue", "black"),
       pch    = c(1, 19, 4))

dev.off()

