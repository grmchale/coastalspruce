## BASAL AREA INCREMENT DATA OF RED SPRUCE CROWNS ##
###### DATA READ IN ##########
library(dplyr)
lidar_metrics <- readRDS("./data/lidar_metrics.rds")
chrono <- readRDS("./data/chrono_VIstats_metrics.rds")
# Filter out extraneous stats
chrono_filt <- chrono %>%
  select(-matches("(SD|Q25|Q75|Mean)$"))
# Join spectra + struc metrics together
struc_spec <- chrono_filt |> inner_join(lidar_metrics, by = "TreeID")
# Remove crowns with no BAI 2024
struc_spec <- struc_spec |> filter(!is.na(BAI_2024))

########## RANDOM FOREST FOR PREDICTING BAI ###########

# Random Forest for BAI 2024
# Inputs: struc_spec, spectral cols 7:398 + lidar cols "rugosity":"zpcum9"
# Output: ranked predictors by mean %IncMSE + stability

library(randomForest)

## 1) Build predictor matrix: spectral 7:398 and lidar rugosity:zpcum9
# spectral
df <- struc_spec
pmin <- 7
pmax <- min(398, ncol(df))
if (pmin > ncol(df)) stop("Data has fewer than 7 columns; check chrono_rf.")
spec_idx <- pmin:pmax

# lidar range by names (inclusive): DBH to Area_m2
lidar_start <- match("DBH", names(df))
lidar_end   <- match("Area_m2", names(df))

if (is.na(lidar_start) || is.na(lidar_end)) {
  stop("Could not find 'DBH' and/or 'Area_m2' in struc_spec column names.")
}
# full range
lidar_idx <- seq(lidar_start, lidar_end)

# exclude problematic columns
exclude_cols <- c("age", "zentropy")
exclude_idx  <- match(exclude_cols, names(df))

# remove any matches (and handle NA safely)
lidar_idx <- setdiff(lidar_idx, exclude_idx[!is.na(exclude_idx)])

# combine and keep unique (in case of overlap)
pred_idx <- unique(c(spec_idx, lidar_idx))

X_raw <- df[, pred_idx, drop = FALSE]
y     <- df[["BAI_2024"]]

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
#X <- as.data.frame(scale(X_raw))
X <- X_raw

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

valid <- !is.na(y)
X_raw <- X_raw[valid, ]
y <- y[valid]

## 5) Aggregate importance + stability
mean_incMSE <- rowMeans(imp_mat, na.rm = TRUE)
sd_incMSE   <- apply(imp_mat, 1, sd, na.rm = TRUE)
rank_rf     <- rank(-mean_incMSE, ties.method = "average")
topK_freq   <- top_hits / B

rf_BAI_summary <- data.frame(
  Predictor   = names(mean_incMSE),
  MeanIncMSE  = as.numeric(mean_incMSE),
  SDIncMSE    = as.numeric(sd_incMSE),
  RF_Rank     = as.numeric(rank_rf),
  TopK_Freq   = as.numeric(topK_freq),
  row.names   = NULL
)

# sort by rank (lower is better)
rf_BAI_summary <- rf_BAI_summary[order(rf_BAI_summary$RF_Rank), ]

## 6) Report: overall robustness + strongest predictors
cat(sprintf("\nRepeated RF runs: B = %d, ntree = %d, mtry = %d, predictors p = %d, rows n = %d\n",
            B, ntree_val, mtry_val, p, nrow(X)))
cat(sprintf("OOB RMSE (mean ± sd): %.2f ± %.2f\n\n",
            sqrt(mean(oob_mse)), sd(sqrt(oob_mse))))

# Show the top 25 predictors (adjust as you like)
head(rf_BAI_summary, 25)

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

####### CREATE CORRELATION MATRIX  - BAI VS. SPECTRAL PREDICTORS #########

## Data read in
library(dplyr)
lidar_metrics <- readRDS("./data/lidar_metrics.rds")
chrono <- readRDS("./data/chrono_VIstats_metrics.rds")
# Filter out extraneous stats
chrono_filt <- chrono %>%
  select(-matches("(SD|Q25|Q75|Mean)$"))
# Join spectra + struc metrics together
struc_spec <- chrono_filt |> inner_join(lidar_metrics, by = "TreeID")
# Remove crowns with no BAI 2024
struc_spec <- struc_spec |> filter(!is.na(BAI_2024))

## Spearman + Pearson correlation matrices - just medians
library(stringr)
library(tidyr) 
library(tibble)
# X: spectral + structural
X <- struc_spec %>%
  select(ARI1_1nm_Median:WBI_15nm_Median,  # spectral
         DBH:Area_m2) %>%                   # structural
  select(where(is.numeric))

# Y: dendro
Y <- struc_spec %>% select(BAI_2024)  # keep as data frame column for cor()

# Spearman correlation matrix (rows = X vars, cols = Y var)
cor_spearman <- cor(as.matrix(X), as.matrix(Y),
                    method = "spearman",
                    use = "pairwise.complete.obs")

cor_spearman_df <- as.data.frame(cor_spearman) %>%
  rownames_to_column("Variable") %>%
  rename(Spearman_r = BAI_2024) %>%
  arrange(desc(abs(Spearman_r)))

## Pearson correlation matrix
cor_pearson <- cor(as.matrix(X), as.matrix(Y),
                   method = "pearson",
                   use = "pairwise.complete.obs")

cor_pearson_df <- as.data.frame(cor_pearson) %>%
  rownames_to_column("Variable") %>%
  rename(Pearson_r = BAI_2024) %>%
  arrange(desc(abs(Pearson_r)))

## Combined table
cor_combined <- left_join(cor_spearman_df, cor_pearson_df, by = "Variable")

print(cor_combined)

## Export results to CSV =
out_dir <- "./outputs"
write.csv(cor_combined, file.path(out_dir, "BAI2024_spearpears_corr.csv"))

# REPLE_5nm_Median is best across Spearman + Pearson correlations

##### Consensus ranking across Spearman, Pearson, and RF #####

# 1. Rank by absolute correlation value (rank 1 = strongest)
spearman_ranked <- cor_spearman_df %>%
  mutate(Spearman_Rank = rank(-abs(Spearman_r))) %>%
  select(Variable, Spearman_r, Spearman_Rank)

pearson_ranked <- cor_pearson_df %>%
  mutate(Pearson_Rank = rank(-abs(Pearson_r))) %>%
  select(Variable, Pearson_r, Pearson_Rank)

# 2. Prep RF df to match naming
rf_ranked <- rf_BAI_summary %>%
  rename(Variable = Predictor) %>%
  select(Variable, RF_Rank)

# 3. Join all three
consensus_df <- spearman_ranked %>%
  left_join(pearson_ranked, by = "Variable") %>%
  left_join(rf_ranked,     by = "Variable") %>%
  # 4. Average rank across all three (lower = better)
  mutate(Mean_Rank = rowMeans(cbind(Spearman_Rank, Pearson_Rank, RF_Rank),
                              na.rm = TRUE)) %>%
  arrange(Mean_Rank)

print(consensus_df)

write.csv(consensus_df, file.path(out_dir, "BAI2024_consensus_ranking.csv"), row.names = FALSE)


######## LME MODELS FOR BAI 2024, TESTS MULTIPLE ##########

### Create model list ####
# a priori selection based on correlations and random forest above
library(lme4)
library(MuMIn)

# --- 1. Candidate predictors ---
candidate_vars <- c("Maccioni_5nm_Median", "zpcum9", "zpcum8",
                    "Area_m2", "REPLE_5nm_Median", "PRI_5nm_Median",
                    "PRInorm_5nm_Median")

# --- 2. Collinearity check
cor_candidates <- cor(struc_spec[, candidate_vars],
                      use = "pairwise.complete.obs",
                      method = "pearson")

# Tidy the correlation matrix into a flagged table, threshold (Pearson |r| > 0.70)
cor_cand_df <- as.data.frame(as.table(cor_candidates)) %>%
  rename(Var1 = Var1, Var2 = Var2, Pearson_r = Freq) %>%
  filter(Var1 != Var2,                      # remove self-correlations
         as.character(Var1) < as.character(Var2)) %>%  # remove duplicates
  mutate(Collinear = abs(Pearson_r) > 0.70) %>%
  arrange(desc(abs(Pearson_r)))

cat("=== Collinearity Table (all pairs) ===\n")
print(cor_cand_df)

cat("\n=== Flagged pairs (|r| > 0.70) ===\n")
print(filter(cor_cand_df, Collinear))

out_dir <- "./outputs"
write.csv(cor_cand_df, file.path(out_dir, "BAI2024_predictor_collinearity.csv"), row.names = FALSE)

# --- 3. Model setup ---
# Prep data
model_data <- struc_spec %>%
  select(all_of(candidate_vars), BAI_2024, Site) %>%
  filter(complete.cases(.))

# Structural vars
S  <- c("zpcum9", "zpcum8", "Area_m2")
# Spectral vars
Sp <- c("Maccioni_5nm_Median", "REPLE_5nm_Median",
        "PRI_5nm_Median", "PRInorm_5nm_Median")

# Known collinear pairs to exclude
excluded_pairs <- list(
  c("PRI_5nm_Median",    "PRInorm_5nm_Median"),
  c("zpcum8",            "zpcum9")
)

# Helper: checks if a model contains a known collinear pair
is_collinear <- function(vars, excluded_pairs) {
  for (pair in excluded_pairs) {
    if (all(pair %in% vars)) return(TRUE)
  }
  return(FALSE)
}

# --- 4. Build model formula list ---
model_formulas <- list()

# (a) Null model
model_formulas[["null"]] <- BAI_2024 ~ 1 + (1 | Site)

# (b) Single predictor models (all 7)
for (v in candidate_vars) {
  name <- paste0("single_", v)
  model_formulas[[name]] <- as.formula(paste("BAI_2024 ~", v, "+ (1 | Site)"))
}

# (c) Structural-only combinations (2-var, no collinear pairs)
struct_combos <- combn(S, 2, simplify = FALSE)
for (combo in struct_combos) {
  if (!is_collinear(combo, excluded_pairs)) {
    name <- paste0("struct_", paste(combo, collapse = "_"))
    model_formulas[[name]] <- as.formula(
      paste("BAI_2024 ~", paste(combo, collapse = " + "), "+ (1 | Site)")
    )
  }
}
# Full structural (3-var)
if (!is_collinear(S, excluded_pairs)) {
  model_formulas[["struct_all"]] <- as.formula(
    paste("BAI_2024 ~", paste(S, collapse = " + "), "+ (1 | Site)")
  )
}

# (d) Spectral-only combinations (2-var, no collinear pairs)
spec_combos <- combn(Sp, 2, simplify = FALSE)
for (combo in spec_combos) {
  if (!is_collinear(combo, excluded_pairs)) {
    name <- paste0("spec_", paste(combo, collapse = "_"))
    model_formulas[[name]] <- as.formula(
      paste("BAI_2024 ~", paste(combo, collapse = " + "), "+ (1 | Site)")
    )
  }
}

# (e) Mixed structural + spectral (1 structural + 1 spectral, no collinear pairs)
for (sv in S) {
  for (spv in Sp) {
    combo <- c(sv, spv)
    if (!is_collinear(combo, excluded_pairs)) {
      name <- paste0("mixed_", sv, "_", spv)
      model_formulas[[name]] <- as.formula(
        paste("BAI_2024 ~", paste(combo, collapse = " + "), "+ (1 | Site)")
      )
    }
  }
}

# (f) 2 structural + 1 spectral models ---
safe_struct_pairs <- list(
  c("zpcum8", "Area_m2"),
  c("zpcum9", "Area_m2")
)

for (pair in safe_struct_pairs) {
  for (spv in Sp) {
    combo <- c(pair, spv)
    name <- paste0("mixed2s_", pair[1], "_", pair[2], "_", spv)
    model_formulas[[name]] <- as.formula(
      paste("BAI_2024 ~", paste(combo, collapse = " + "), "+ (1 | Site)")
    )
  }
}
cat("\n=== Candidate models (", length(model_formulas), "total) ===\n")
for (nm in names(model_formulas)) {
  cat(nm, ":", deparse(model_formulas[[nm]]), "\n")
}

#### Test model list ####
library(lme4)
library(MuMIn)  # AICc, r.squaredGLMM

# --- Fit all models ---
fitted_models <- list()
for (nm in names(model_formulas)) {
  fitted_models[[nm]] <- tryCatch(
    lmer(model_formulas[[nm]], data = model_data, REML = FALSE),
    error   = function(e) { cat("ERROR in", nm, ":", e$message, "\n"); NULL },
    warning = function(w) { cat("WARNING in", nm, ":", w$message, "\n")
      suppressWarnings(lmer(model_formulas[[nm]], data = model_data, REML = FALSE)) }
  )
}

# --- Extract metrics ---
model_summary <- do.call(rbind, lapply(names(fitted_models), function(nm) {
  m <- fitted_models[[nm]]
  # add this line right after:  m <- fitted_models[[nm]]
  singular <- isSingular(m)
  if (is.null(m)) {
    return(data.frame(Model = nm, K = NA, AICc = NA, Delta_AICc = NA,
                      AICc_Weight = NA, Marginal_R2 = NA,
                      Conditional_R2 = NA, RMSE = NA))
  }
  
  # K = fixed effects + random intercept variance + residual variance
  k       <- length(fixef(m)) + 2
  
  # R-squared (marginal + conditional)
  r2      <- tryCatch(r.squaredGLMM(m), error = function(e) c(NA, NA))
  
  # RMSE
  rmse    <- sqrt(mean(residuals(m)^2))
  
  data.frame(
    Model          = nm,
    K              = k,
    AICc           = AICc(m),
    Delta_AICc     = NA,
    AICc_Weight    = NA,
    Marginal_R2    = round(r2[1], 3),
    Conditional_R2 = round(r2[2], 3),
    RMSE           = round(rmse, 2),
    Singular       = singular       # <-- add this
  )
}))

# --- Compute Delta_AICc and AICc weights ---
model_summary <- model_summary %>%
  filter(!is.na(AICc)) %>%
  mutate(
    Delta_AICc  = round(AICc - min(AICc), 2),
    AICc_Weight = round(exp(-0.5 * Delta_AICc) / sum(exp(-0.5 * Delta_AICc)), 3)
  ) %>%
  arrange(Delta_AICc)

print(model_summary)

write.csv(model_summary, file.path(out_dir, "BAI2024_LMEmodel_comparison.csv"), row.names = FALSE)

##### Move forward with the best model ####
top_model_lm <- lm(BAI_2024 ~ zpcum8 + Area_m2 + REPLE_5nm_Median, data = model_data)  # adjust predictors
summary(top_model_lm)

