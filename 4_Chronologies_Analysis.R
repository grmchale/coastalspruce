####################################################################################
################ COMPARE CHRONOLOGY METRICS AND SPECTRA ! ##########################
setwd("G:/HyperspectralUAV")
###################### READ IN DATAFRAMES ################################
library(dplyr)

chronologies_VIstats_combined <- readRDS("./R_outputs/speclib_chronologies/veg_indices/chronologies_VIstats_combined.rds")
treeage <- read.csv("G:/Colbys_Plots/spectra_analyses/treeage.csv") # Recruitment age for crown
BAI_wide <- read.csv("G:/Colbys_Plots/spectra_analyses/BAI_wide.csv") # BAI from last 20 years + avg. BAI from last few
cindex <- read.csv("G:/Colbys_Plots/spectra_analyses/hegyibytree.csv") # Competition index for each tree
chronology_crown_attr <-read.csv("G:/Colbys_Plots/spectra_analyses/chronology_crowns.csv") # Site, Species, Tag columns
crown_rugosity <-read.csv("G:/LiDAR/Normalized/LiDAR_crowns/crown_rugosity.csv") # Rugosity of crowns
crown_metrics <-read.csv("G:/LiDAR/Normalized/LiDAR_crowns/crown_metrics.csv") # Point cloud metrics
crown_areas <- read.csv("G:/LiDAR/Normalized/LiDAR_crowns/crown_areas.csv") # Area of each crown
crown_latlong <- read.csv("G:/LiDAR/Latitude/chroncrown_latlon.csv")

###################### PREP DATAFRAME ################################
# Merge duplicate file in spectra df
merge_pair <- function(df, base_id, suffix) {
  ids <- c(base_id, paste0(base_id, suffix))
  if (!all(ids %in% df$Sample)) {
    warning("Skipping merge: not all IDs found: ", paste(ids, collapse = " & "))
    return(df)
  }
  num_cols <- names(df)[sapply(df, is.numeric)]
  merged_vals <- df %>%
    filter(Sample %in% ids) %>%
    summarise(across(all_of(num_cols), ~ mean(.x, na.rm = TRUE)))
  merged_row <- tibble(Sample = base_id) %>% bind_cols(merged_vals)
  df %>%
    filter(!Sample %in% ids) %>%
    bind_rows(merged_row) %>%
    select(names(df))  # preserve original column order
}
chronologies_VIstats_clean <- chronologies_VIstats_combined %>%
  merge_pair("GW_PIRU_14.4_9.5_D", "_6") %>%
  merge_pair("GW_PIRU_12_37.5_D", "_7")

# Join attributes to VI dataframe
chronologies_VIstats_clean <- chronologies_VIstats_clean %>% rename(TreeID = Sample)
crown_keep <- chronology_crown_attr %>%
  select(TreeID, CC, Site, Tag, Dndrmtr) %>%
  distinct(TreeID, .keep_all = TRUE)   # handles the 2 duplicate GW_* rows
chronologies_VIstats_joined <- chronologies_VIstats_clean %>%
  left_join(crown_keep, by = "TreeID") %>%
  relocate(CC, Site, Tag, Dndrmtr, .after = TreeID)

# Create join field for chronology metrics
chronologies_VIstats_joined <- chronologies_VIstats_joined %>%
  mutate(ChronID = paste0(Site, Tag)) %>%   # e.g., "CC" + "67" -> "CC67"
  relocate(ChronID, .after = Tag)

# Join chronology metrics from 3 dfs - keep only needed columns
BAI_keep <- BAI_wide %>%
  select(TreeID, all_of(names(BAI_wide)[3:ncol(BAI_wide)])) %>%
  distinct(TreeID, .keep_all = TRUE)
cindex_keep <- cindex %>%
  select(TreeID, cindex, DBH) %>%
  distinct(TreeID, .keep_all = TRUE)
treeage_keep <- treeage %>%
  select(SeriesID, age) %>%
  distinct(SeriesID, .keep_all = TRUE)

chrono_final <- chronologies_VIstats_joined %>%
  left_join(BAI_keep,   by = c("ChronID" = "TreeID")) %>%
  left_join(cindex_keep, by = c("ChronID" = "TreeID")) %>%
  left_join(treeage_keep, by = c("ChronID" = "SeriesID"))

write.csv(chrono_final,"./R_outputs/speclib_chronologies/veg_indices/chrono_VIstats_metrics.csv")
saveRDS(chrono_final,"./R_outputs/speclib_chronologies/veg_indices/chrono_VIstats_metrics.rds")

###################### ANALYSIS: CORRELATION MATRIX #############################
chrono <- readRDS("./R_outputs/speclib_chronologies/veg_indices/chrono_VIstats_metrics.rds")

library(dplyr)

##### Compute raw correlation matrix #######
# Split spectral vs dendro variables
spectral_vars <- chrono[, 7:1966]
dendro_vars   <- chrono[, 1967:2004]

# Compute Spearman correlations (pairwise complete obs)
cor_matrix <- cor(
  spectral_vars,
  dendro_vars,
  method = "spearman",
  use = "pairwise.complete.obs"
)

# Convert to dataframe for easy viewing
cor_df <- as.data.frame(cor_matrix)

# Optional: add a column for spectral variable names
cor_df <- tibble::rownames_to_column(cor_df, var = "Spectral_Variable")

# Largest correlation for every dendrochronology variable
# assumes cor_df has Spectral_Variable as first column
largest_corr <- lapply(names(cor_df)[-1], function(col) {
  tmp <- cor_df %>%
    select(Spectral_Variable, all_of(col)) %>%
    arrange(desc(abs(.data[[col]]))) %>%
    slice(1)
  tmp$Dendro_Var <- col
  tmp
}) %>%
  bind_rows()

# Top 5 correlations for every dendrochronology variable
top5_corr <- lapply(names(cor_df)[-1], function(col) {
  tmp <- cor_df %>%
    select(Spectral_Variable, all_of(col)) %>%
    arrange(desc(abs(.data[[col]]))) %>%
    slice(1:5)
  tmp$Dendro_Var <- col
  tmp
}) %>%
  bind_rows()

####### Compute correlation matrix without SD, Q25, Q75 #######
chrono_filt <- chrono %>%
  select(-matches("(SD|Q25|Q75)$"))

# --- Split spectral vs dendrochronology variables ---
spectral_vars_f <- chrono_filt[, 7:790]
dendro_vars_f   <- chrono_filt[, 791:828]

# --- Compute Spearman correlations ---
cor_matrix_filt <- cor(
  spectral_vars_f,
  dendro_vars_f,
  method = "spearman",
  use = "pairwise.complete.obs"
)

# --- Convert to dataframe for viewing ---
cor_df_clean <- as.data.frame(cor_matrix_filt)
cor_df_clean <- tibble::rownames_to_column(cor_df_clean, var = "Spectral_Variable")

largest_corr_clean <- lapply(names(cor_df_clean)[-1], function(col) {
  tmp <- cor_df_clean %>%
    select(Spectral_Variable, all_of(col)) %>%
    arrange(desc(abs(.data[[col]]))) %>%
    slice(1)
  tmp$Dendro_Var <- col
  tmp
}) %>%
  bind_rows()

top5_corr_clean <- lapply(names(cor_df_clean)[-1], function(col) {
  tmp <- cor_df_clean %>%
    select(Spectral_Variable, all_of(col)) %>%
    arrange(desc(abs(.data[[col]]))) %>%
    slice(1:5)
  tmp$Dendro_Var <- col
  tmp
}) %>%
  bind_rows()

####### Compute correlation with just dendrometer trees ########
library(dplyr)
library(tibble)

# --- Filter to crowns with BAI_2024 ---
df_bai <- chrono_filt %>% filter(!is.na(BAI_2024))

# --- Split spectral vs dendrochronology variables ---
spectral_vars_bai <- df_bai[, 7:790]
dendro_vars_bai   <- df_bai[, 791:829]

# --- Compute Spearman correlations ---
cor_matrix_bai <- cor(
  spectral_vars_bai,
  dendro_vars_bai,
  method = "spearman",
  use = "pairwise.complete.obs"
)

# --- Convert to dataframe for viewing ---
cor_df_bai <- as.data.frame(cor_matrix_bai)
cor_df_bai <- rownames_to_column(cor_df_bai, var = "Spectral_Variable")

# columns to iterate (all dendro vars)
cols <- names(cor_df_bai)[-1]  # drop "Spectral_Variable"

# Largest correlation per dendro var -
largest_corr_bai <- do.call(rbind, lapply(cols, function(col) {
  ord <- order(abs(cor_df_bai[[col]]), decreasing = TRUE)
  i   <- ord[1]
  data.frame(
    Spectral_Variable = cor_df_bai$Spectral_Variable[i],
    Correlation       = cor_df_bai[[col]][i],
    Dendro_Var        = col,
    stringsAsFactors  = FALSE
  )
}))
row.names(largest_corr_bai) <- NULL
# optional: sort by |Correlation|
largest_corr_bai <- largest_corr_bai[order(-abs(largest_corr_bai$Correlation)), ]

# -Top 5 correlations per dendro var --
top5_corr_bai <- do.call(rbind, lapply(cols, function(col) {
  ord <- order(abs(cor_df_bai[[col]]), decreasing = TRUE)
  k   <- seq_len(min(5, length(ord)))
  idx <- ord[k]
  data.frame(
    Spectral_Variable = cor_df_bai$Spectral_Variable[idx],
    Correlation       = cor_df_bai[[col]][idx],
    Dendro_Var        = rep(col, length(idx)),
    stringsAsFactors  = FALSE
  )
}))
row.names(top5_corr_bai) <- NULL
# optional: sort nicely
top5_corr_bai <- top5_corr_bai[order(top5_corr_bai$Dendro_Var, -abs(top5_corr_bai$Correlation)), ]



####### Compute correlation matrix without I or S trees #######
# 1. Filter out I and S crown classes
chrono_noIS <- chrono_filt %>%
  dplyr::filter(!(CC %in% c("I", "S")))
# 2. Spearman correlation matrix
cor_df_noIS <- cor(
  chrono_noIS[, 7:790],   # spectral variables
  chrono_noIS[, 791:828],# dendrochronology variables
  method = "spearman",
  use = "pairwise.complete.obs"
)

cor_df_noIS <- as.data.frame(cor_df_noIS)

# 3. Largest correlation (per dendro var)
largest_corr_noIS <- apply(cor_df_noIS, 2, function(col) {
  col[which.max(abs(col))]
})

largest_corr_noIS <- data.frame(
  DendroVar = names(largest_corr_noIS),
  CorrValue = as.numeric(largest_corr_noIS),
  SpectralVar = apply(cor_df_noIS, 2, function(col) {
    rownames(cor_df_noIS)[which.max(abs(col))]
  })
)

# 4. Top 5 correlations (per dendro var)
top5_corr_noIS <- lapply(names(cor_df_noIS), function(var) {
  vals <- cor_df_noIS[[var]]
  names(vals) <- rownames(cor_df_noIS)
  vals <- sort(vals, decreasing = TRUE, na.last = NA)
  
  # sort by absolute correlation
  vals <- sort(vals, decreasing = TRUE, na.last = NA, 
               method = "quick", index.return = FALSE)
  top5 <- head(vals[order(abs(vals), decreasing = TRUE)], 5)
  
  data.frame(
    DendroVar = var,
    SpectralVar = names(top5),
    CorrValue = as.numeric(top5),
    AbsCorr = abs(as.numeric(top5))
  )
})

top5_corr_noIS <- do.call(rbind, top5_corr_noIS)

####### Compute correlation matrix, just medians no S trees (FINAL) #######
library(dplyr)
chrono <- readRDS("./R_outputs/speclib_chronologies/veg_indices/chrono_VIstats_metrics.rds")
chrono_final <- chrono %>%
  # remove columns ending in SD, Q25, Q75, Mean
  select(-matches("(SD|Q25|Q75|Mean)$")) %>%
  # filter out S crown classes
  filter(!CC %in% c("S")) %>%
  # join with crown_metrics and crown_rugosity by TreeID
  left_join(crown_metrics,  by = "TreeID") %>%
  left_join(crown_rugosity, by = "TreeID")

## Spearman correlation matrix
# X: spectral + structural
X <- chrono_final %>%
  select(ARI1_1nm_Median:WBI_15nm_Median,  # spectral
         zmax:rugosity) %>%                 # structural
  select(where(is.numeric))

# Y: dendro
Y <- chrono_final %>%
  select(BAI_2000:age) %>%
  select(where(is.numeric))

# Spearman correlation matrix (rows = X vars, cols = Y vars)
cor_final <- cor(as.matrix(X), as.matrix(Y),
                    method = "spearman",
                    use = "pairwise.complete.obs")
cor_df_final <- as.data.frame(cor_final)

library(stringr)
library(tidyr) 

# Convert correlation matrix to just indices
cor_long <- cor_df_final %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Index") %>%
  pivot_longer(-Index, names_to = "Dendro", values_to = "Correlation")

cor_vi <- cor_long %>%
  filter(str_ends(Index, "Median"))

# Parse index name + resample (split at underscores) ---
cor_vi <- cor_vi %>%
  mutate(
    Index_clean = str_remove(Index, "_[0-9]+nm_Median$"),  # e.g. "SR8"
    Resample    = str_extract(Index, "[0-9]+nm"),          # e.g. "15nm"
    AbsCorr     = abs(Correlation)
  )

# --- Top 20 for BAI_2024 ---
top20_BAI_2024 <- cor_vi %>%
  filter(Dendro == "BAI_2024") %>%
  arrange(desc(AbsCorr)) %>%
  slice_head(n = 20) %>%
  select(Index_clean, Resample, Correlation, AbsCorr)

# --- Top 20 for age ---
top20_age <- cor_vi %>%
  filter(Dendro == "age") %>%
  arrange(desc(AbsCorr)) %>%
  slice_head(n = 20) %>%
  select(Index_clean, Resample, Correlation, AbsCorr)

install.packages("DT")
library(DT)

# Interactive, copy-friendly tables in Viewer
library(dplyr)
library(DT)

# Round correlations to 3 decimal places
top20_BAI_2024_round <- top20_BAI_2024 %>%
  mutate(
    Correlation = round(Correlation, 3),
    AbsCorr     = round(AbsCorr, 3)
  )

top20_age_round <- top20_age %>%
  mutate(
    Correlation = round(Correlation, 3),
    AbsCorr     = round(AbsCorr, 3)
  )

# Send to Viewer
datatable(top20_BAI_2024_round,
          caption = "Top 20 vegetation indices correlated with BAI_2024",
          options = list(pageLength = 20, autoWidth = TRUE))

datatable(top20_age_round,
          caption = "Top 20 vegetation indices correlated with Age",
          options = list(pageLength = 20, autoWidth = TRUE))

# = EXPORT RESULTS TO CSV =
out_dir <- "G:/HyperspectralUAV/R_outputs/modelling/correlations/dendrochronology_correlations"
write.csv(cor_df_final,
          file.path(out_dir, "cor_df_final.csv"),
          row.names = TRUE)
write.csv(top20_age_round,
          file.path(out_dir, "top20_age_round.csv"),
          row.names = FALSE)
write.csv(top20_BAI_2024_round,
          file.path(out_dir, "top20_BAI_2024_round.csv"),
          row.names = FALSE)

####### Use chrono_age to generate scatterplots (spec_var on y-axis) #######
library(dplyr)
library(ggplot2)

# = USER SETTINGS =
spec_var        <- "REPLE_5nm_Median"  # <- spectral column in chrono_age (Y axis)
dendro_var      <- "BAI_2024"          # <- dendro column in chrono_age (X axis)
reg_type        <- "linear"           # <- "linear", "quadratic", or "exponential"
show_regression <- TRUE               # <- TRUE to draw regression line + R² label on plot
point_alpha     <- 0.7
point_size      <- 3.0

# Manual plot labels
plot_title <- "Red Edge Position Linear Estimation (REPLE) vs 2024 BAI"
x_label    <- "BAI in 2024"
y_label    <- "REPLE (5nm median)"

# = PREP =
stopifnot(spec_var %in% names(chrono_age),
          dendro_var %in% names(chrono_age))

df_xy <- chrono_age %>%
  dplyr::select(Site,
                x = dplyr::all_of(dendro_var),   # dendro on x
                y = dplyr::all_of(spec_var)) %>% # spectral on y
  tidyr::drop_na()

# = FIT + METRICS =
r2_text <- ""
fit_line_df <- NULL
metrics_note <- NULL

if (reg_type == "linear") {
  fit <- lm(y ~ x, data = df_xy)
  sm  <- summary(fit)
  r2  <- sm$r.squared
  adj <- sm$adj.r.squared
  slope <- coef(fit)[["x"]]
  intercept <- coef(fit)[["(Intercept)"]]
  p_value <- sm$coefficients["x", "Pr(>|t|)"]
  sigma <- sm$sigma
  fstat <- unname(sm$fstatistic["value"])
  df_num <- unname(sm$fstatistic["numdf"])
  df_den <- unname(sm$fstatistic["dendf"])
  
  r2_text <- sprintf("R² = %.3f", r2)
  metrics_note <- sprintf(
    "\n--- Linear model: %s ~ %s ---\nIntercept: %.4f\nSlope: %.4f\nR²: %.4f\nAdj R²: %.4f\nResidual SE: %.4f\nF-statistic: %.3f on %d and %d df\np-value (slope): %.4g\n",
    spec_var, dendro_var, intercept, slope, r2, adj, sigma, fstat, df_num, df_den, p_value
  )
  
} else if (reg_type == "quadratic") {
  fit <- lm(y ~ poly(x, 2, raw = TRUE), data = df_xy)
  sm  <- summary(fit)
  r2  <- sm$r.squared
  adj <- sm$adj.r.squared
  beta0 <- coef(fit)[["(Intercept)"]]
  beta1 <- coef(fit)[["poly(x, 2, raw = TRUE)1"]]
  beta2 <- coef(fit)[["poly(x, 2, raw = TRUE)2"]]
  sigma <- sm$sigma
  fstat <- unname(sm$fstatistic["value"])
  df_num <- unname(sm$fstatistic["numdf"])
  df_den <- unname(sm$fstatistic["dendf"])
  
  r2_text <- sprintf("R² = %.3f", r2)
  metrics_note <- sprintf(
    "\n--- Quadratic model: %s ~ %s + %s^2 ---\nIntercept: %.4f\nx coef: %.4f\nx^2 coef: %.4f\nR²: %.4f\nAdj R²: %.4f\nResidual SE: %.4f\nF-statistic: %.3f on %d and %d df\n",
    spec_var, dendro_var, dendro_var, beta0, beta1, beta2, r2, adj, sigma, fstat, df_num, df_den
  )
  
} else if (reg_type == "exponential") {
  df_fit <- dplyr::filter(df_xy, y > 0)
  if (nrow(df_fit) < 2) stop("Not enough positive y values for exponential fit.")
  fit <- lm(log(y) ~ x, data = df_fit)
  sm  <- summary(fit)
  preds <- exp(predict(fit, newdata = df_fit))
  r2 <- 1 - sum((df_fit$y - preds)^2, na.rm = TRUE) /
    sum((df_fit$y - mean(df_fit$y, na.rm = TRUE))^2, na.rm = TRUE)
  a <- coef(fit)[["(Intercept)"]]
  b <- coef(fit)[["x"]]
  A <- exp(a)
  p_value <- sm$coefficients["x", "Pr(>|t|)"]
  sigma_log <- sm$sigma
  xgrid <- data.frame(x = seq(min(df_fit$x, na.rm = TRUE),
                              max(df_fit$x, na.rm = TRUE),
                              length.out = 200))
  fit_line_df <- transform(xgrid, y = exp(predict(fit, newdata = xgrid)))
  
  r2_text <- sprintf("R² = %.3f", r2)
  metrics_note <- sprintf(
    "\n--- Exponential model: %s = A * exp(b * %s) ---\nA: %.6f\nb: %.6f\nR²: %.4f\nResidual SE (log-space): %.4f\np-value (b): %.4g\n",
    spec_var, dendro_var, A, b, r2, sigma_log, p_value
  )
  
} else {
  stop("reg_type must be one of: 'linear', 'quadratic', 'exponential'")
}

# Print metrics to console
cat(metrics_note)

# = PLOT =
p <- ggplot(df_xy, aes(x = x, y = y, color = Site)) +
  geom_point(alpha = point_alpha, size = point_size) +
  labs(
    title = plot_title,
    x = x_label,
    y = y_label,
    color = "Site"
  ) +
  theme_minimal(base_size = 14)

if (show_regression) {
  if (reg_type == "linear") {
    p <- p + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black")
  } else if (reg_type == "quadratic") {
    p <- p + geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "black")
  } else if (reg_type == "exponential") {
    p <- p + geom_line(data = fit_line_df, aes(x = x, y = y), inherit.aes = FALSE, color = "black")
  }
  p <- p + annotate("text", x = Inf, y = Inf, label = r2_text,
                    hjust = 1.02, vjust = 1.5, size = 4.2)
}

p

### Quick QQ plots
# QQ plots for three variables
par(mfrow = c(1, 3))  # put 3 plots side by side

# BAI_2024
qqnorm(chrono_age$BAI_2024, main = "QQ Plot: BAI_2024")
qqline(chrono_age$BAI_2024, col = "red")

# PRInorm_10nm_Median
qqnorm(chrono_age$PRInorm_10nm_Median, main = "QQ Plot: PRInorm (10nm)")
qqline(chrono_age$PRInorm_10nm_Median, col = "red")

# REPLE_5nm_Median
qqnorm(chrono_age$REPLE_5nm_Median, main = "QQ Plot: REPLE (5nm)")
qqline(chrono_age$REPLE_5nm_Median, col = "red")

par(mfrow = c(1, 1))  # reset layout



####
# Quick scatterplot test
library(dplyr)
library(ggplot2)

# Summarize means per site
site_means <- chrono_filt %>%
  filter(!is.na(BAI_2024), !is.na(PRInorm_10nm_Mean)) %>%
  group_by(Site) %>%
  summarise(
    PRI_mean = mean(PRInorm_10nm_Mean, na.rm = TRUE),
    BAI_mean = mean(BAI_2024, na.rm = TRUE),
    .groups = "drop"
  )

print(site_means)

# Scatterplot of site-level means
ggplot(site_means, aes(x = PRI_mean, y = BAI_mean, label = Site)) +
  geom_point(size = 3, color = "blue") +
  geom_text(vjust = -1, size = 3) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Mean PRI (10nm)",
    y = "Mean BAI_2024",
    title = "Site-level means of PRI vs BAI_2024"
  )

################ AGE vs. BEST INDEX: BIN INTO AGE GROUPS ################
library(dplyr)
chrono <- readRDS("./R_outputs/speclib_chronologies/veg_indices/chrono_VIstats_metrics.rds")
chrono_final <- chrono %>%
  # remove columns ending in SD, Q25, Q75, Mean
  select(-matches("(SD|Q25|Q75|Mean)$")) %>%
  # filter out S crown classes
  filter(!CC %in% c("S")) #%>%
  # join with crown_metrics and crown_rugosity by TreeID
  #left_join(crown_metrics,  by = "TreeID") %>%
  #left_join(crown_rugosity, by = "TreeID")

chrono_age <- chrono_final %>%
  filter(!is.na(age)) %>%
  mutate(age_new = 2024 - age)

# Bin ages into 5-year intervals and summarize
age_bins_5 <- chrono_age %>%
  mutate(age_bin = cut(
    age_new,
    breaks = seq(floor(min(age_new, na.rm = TRUE)),
                 ceiling(max(age_new, na.rm = TRUE)),
                 by = 5),
    right = TRUE,
    include.lowest = TRUE
  )) %>%
  group_by(age_bin) %>%
  summarise(
    n          = n(),
    EGFR_mean  = mean(EGFR_5nm_Median, na.rm = TRUE),
    EGFR_median= median(EGFR_5nm_Median, na.rm = TRUE),
    EGFR_sd    = sd(EGFR_5nm_Median, na.rm = TRUE),
    EGFN_mean  = mean(EGFN_5nm_Median, na.rm = TRUE),
    EGFN_median= median(EGFN_5nm_Median, na.rm = TRUE),
    EGFN_sd    = sd(EGFN_5nm_Median, na.rm = TRUE)
  )
# QQ plot for EGFR mean values across bins
qqnorm(age_bins_5$EGFR_mean, main = "QQ Plot of EGFR bin means")
qqline(age_bins_5$EGFR_mean, col = "red")

# QQ plot for raw EGFR values
qqnorm(chrono_age$EGFR_5nm_Median, main = "QQ Plot of EGFR_5nm_Median")
qqline(chrono_age$EGFR_5nm_Median, col = "red")

# QQ plot for raw EGFN values
qqnorm(chrono_age$EGFN_5nm_Median, main = "QQ Plot of EGFN_5nm_Median")
qqline(chrono_age$EGFN_5nm_Median, col = "blue")

# Filter rows with age_new between 146 and 171
age_rows_146_171 <- chrono_age %>%
  filter(age_new >= 146 & age_new <= 171) %>%
  select(TreeID, CC, Site, age_new)



library(dplyr)
library(ggplot2)
library(stringr)

### Scatterplot for binned age vs. specta
value_col <- "EGFR_mean"   # e.g., "EGFR_mean", "EGFR_median", "EGFN_mean", "EGFN_median"
sd_col    <- "EGFR_sd"     # matching SD column, e.g., "EGFR_sd" or "EGFN_sd"
y_lab     <- "EGFR (mean ± SD)"
point_alpha <- 0.9
point_size  <- 3
errorbar_width <- 0.0

# Helper to parse an age_bin label like "[146,151]" or "(146,151]"
parse_bin_bounds <- function(lbl) {
  core <- gsub("\\[|\\]|\\(|\\)", "", lbl)
  parts <- strsplit(core, ",")[[1]]
  as.numeric(parts)
}

# Build plotting df from age_bins_5
plot_df <- age_bins_5 %>%
  mutate(
    age_bin_chr = as.character(age_bin),
    bounds      = lapply(age_bin_chr, parse_bin_bounds),
    bin_start   = vapply(bounds, function(z) z[1], numeric(1)),
    bin_end     = vapply(bounds,   function(z) z[2], numeric(1)),
    bin_mid     = (bin_start + bin_end) / 2,
    # compact label like "146–151"
    bin_lab     = sprintf("%g–%g", bin_start, bin_end)
  ) %>%
  arrange(bin_start)

# Ticks: ensure 1:1 breaks/labels
# Explicit namespace call to avoid conflicts
tick_df <- plot_df %>%
  dplyr::select(bin_mid, bin_lab) %>%
  dplyr::distinct() %>%
  dplyr::arrange(bin_mid)

# ---- Linear model on bin midpoints ----
mod_df <- plot_df %>%
  select(bin_mid, value = all_of(value_col)) %>%
  filter(is.finite(bin_mid), is.finite(value))

model <- lm(value ~ bin_mid, data = mod_df)
sm    <- summary(model)

# Extract model metrics
intercept <- unname(coef(model)[1])
slope     <- unname(coef(model)[2])
p_value   <- sm$coefficients["bin_mid", "Pr(>|t|)"]
r2        <- sm$r.squared
adj_r2    <- sm$adj.r.squared
df_res    <- sm$df[2]
sigma     <- sm$sigma
fstat     <- unname(sm$fstatistic["value"])
df_num    <- unname(sm$fstatistic["numdf"])
df_den    <- unname(sm$fstatistic["dendf"])

# Print metrics to console
cat("\n--- Linear model: ", value_col, " ~ bin_mid ---\n", sep = "")
cat(sprintf("Intercept: %.4f\n", intercept))
cat(sprintf("Slope:     %.4f\n", slope))
cat(sprintf("R²:        %.4f\n", r2))
cat(sprintf("Adj R²:    %.4f\n", adj_r2))
cat(sprintf("Residual SE (sigma): %.4f on %d df\n", sigma, df_res))
cat(sprintf("F-statistic: %.3f on %d and %d df\n", fstat, df_num, df_den))
cat(sprintf("p-value (slope): %.4g\n\n", p_value))

# --- X breaks every 15 years based on bin midpoints ---
xmin <- floor(min(plot_df$bin_mid, na.rm = TRUE) / 15) * 15
xmax <- ceiling(max(plot_df$bin_mid, na.rm = TRUE) / 15) * 15
x_breaks <- seq(xmin, xmax, by = 15)

# Top-of-plot y position for the R² label
y_top <- max(plot_df[[value_col]] + plot_df[[sd_col]], na.rm = TRUE)

ggplot(plot_df, aes(x = bin_mid, y = .data[[value_col]])) +
  geom_point(size = point_size, alpha = point_alpha) +
  geom_errorbar(
    aes(
      ymin = .data[[value_col]] - .data[[sd_col]],
      ymax = .data[[value_col]] + .data[[sd_col]]
    ),
    width = errorbar_width,
    na.rm = TRUE
  ) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_x_continuous(
    name = "Age (binned in 5 yr intervals)",
    breaks = x_breaks,
    labels = x_breaks,
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    y = y_lab,
    title = paste0("Enhanced Green Fraction - Ratio (EGFR) vs. Binned Age")
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  annotate(
    "text",
    x = xmin, y = y_top, hjust = 0, vjust = 1.2,
    label = sprintf("R² = %.3f", r2)
  )

# ---- Plot ----
ggplot(plot_df, aes(x = bin_mid, y = .data[[value_col]])) +
  geom_point(size = point_size, alpha = point_alpha) +
  geom_errorbar(
    aes(
      ymin = .data[[value_col]] - .data[[sd_col]],
      ymax = .data[[value_col]] + .data[[sd_col]]
    ),
    width = errorbar_width,
    na.rm = TRUE
  ) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_x_continuous(
    name   = "Age (binned in 5 yr intervals)",
    breaks = tick_df$bin_mid,
    labels = tick_df$bin_lab,
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    y = y_lab,
    title = paste0("EGFR Mean Age per Bin (±SD, linear regression)")
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




################## ANALYSIS: LINEAR MIXED EFFECTS MODEL ############################
library(dplyr)
library(lme4)
chrono <- readRDS("./R_outputs/speclib_chronologies/veg_indices/chrono_VIstats_metrics.rds")
chrono_filt <- chrono %>%
  select(-matches("(SD|Q25|Q75)$"))
chrono_filt <- chrono_filt %>%
  left_join(
    dplyr::select(crown_rugosity, TreeID, rugosity),
    by = "TreeID"
  )

# make a working dataset
df_mod <- chrono_filt %>%
  filter(!is.na(BAI_2024)) %>%
  mutate(
    rugosity_sc = scale(rugosity),
    age_sc = scale(age),
    cindex_sc = scale(cindex),
    PRI_sc = scale(PRInorm_10nm_Median)
  )
# Model A: age + rugosity
m_age_rug <- lm(BAI_2024 ~ rugosity_sc + PRI_sc, data = df_mod)

# Model B: age + PRI
m_age_pri <- lm(BAI_2024 ~ age_sc + PRI_sc, data = df_mod)

# Model C: cindex + PRI
m_cindex_pri <- lm(BAI_2024 ~ cindex_sc + PRI_sc, data = df_mod)

# compare summaries
summary(m_age_rug)
summary(m_age_pri)
summary(m_cindex_pri)

# compare AIC
AIC(m_age_rug, m_age_pri, m_cindex_pri)

# combined fixed-effects model
m_combo <- lm(BAI_2024 ~ PRI_sc + age_sc + cindex_sc, data = df_mod)
summary(m_combo)

# pri only
m_pri_only <- lm(BAI_2024 ~ PRI_sc, data = df_mod)
summary(m_pri_only)
AIC(m_pri_only, m_combo)

# combined fixed-effects with random site
m_combo_RE <- lmer(BAI_2024 ~ PRI_sc + age_sc + cindex_sc + (1 | Site),
                   data = df_mod, REML = FALSE)
summary(m_combo_RE)

# combined fixed effects with fixed site
m_combo_FE <- lm(BAI_2024 ~ PRI_sc + age_sc + cindex_sc + Site, data = df_mod)
summary(m_combo_FE)

# Test cross validation
sites <- unique(df_mod$Site)
cv_results <- data.frame(Site=character(), Model=character(),
                         RMSE=double(), R2=double(), stringsAsFactors=FALSE)

for (s in sites) {
  train <- df_mod %>% filter(Site != s)
  test  <- df_mod %>% filter(Site == s)
  
  # PRI-only
  fit1 <- lm(BAI_2024 ~ PRI_sc, data=train)
  pred1 <- predict(fit1, newdata=test)
  rmse1 <- sqrt(mean((test$BAI_2024 - pred1)^2))
  r21   <- cor(test$BAI_2024, pred1)^2
  cv_results <- rbind(cv_results, data.frame(Site=s, Model="PRI_only", RMSE=rmse1, R2=r21))
  
  # PRI + age + cindex
  fit2 <- lm(BAI_2024 ~ PRI_sc + age_sc + cindex_sc, data=train)
  pred2 <- predict(fit2, newdata=test)
  rmse2 <- sqrt(mean((test$BAI_2024 - pred2)^2))
  r22   <- cor(test$BAI_2024, pred2)^2
  cv_results <- rbind(cv_results, data.frame(Site=s, Model="PRI+age+cindex", RMSE=rmse2, R2=r22))
}

cv_results

################## NON-PARAMETRIC MODELS #################################
library(dplyr)
chrono <- readRDS("./R_outputs/speclib_chronologies/veg_indices/chrono_VIstats_metrics.rds")
chrono_filt <- chrono %>%
  select(-matches("(SD|Q25|Q75)$"))
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

# Attempt at LMEM
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

#AGE VS ZMAX
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




################# OUTLIER ANALYSIS - BAI 2024 #########################
library(dplyr)
chrono <- readRDS("./R_outputs/speclib_chronologies/veg_indices/chrono_VIstats_metrics.rds")
chrono_filt <- chrono %>%
  select(-matches("(SD|Q25|Q75)$"))
chrono_bai24 <- chrono_filt %>%
  filter(!is.na(BAI_2024))
# --- 1. Calculate IQR-based thresholds ---
Q1 <- quantile(chrono_bai24$PRInorm_10nm_Median, 0.25, na.rm = TRUE)
Q3 <- quantile(chrono_bai24$PRInorm_10nm_Median, 0.75, na.rm = TRUE)
IQR_val <- Q3 - Q1

lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

# --- 2. Flag outliers ---
outliers <- chrono_bai24 %>%
  filter(PRInorm_10nm_Median < lower_bound |
           PRInorm_10nm_Median > upper_bound)

# --- 3. Preview them ---
outliers


# --- 1) IQR thresholds (same rule you used) ---
Q1 <- quantile(chrono_bai24$PRInorm_10nm_Median, 0.25, na.rm = TRUE)
Q3 <- quantile(chrono_bai24$PRInorm_10nm_Median, 0.75, na.rm = TRUE)
IQR_val <- Q3 - Q1
lower_bound <- Q1 - 1.5 * IQR_val
upper_bound <- Q3 + 1.5 * IQR_val

# --- 2) Remove outliers & NAs in either axis variable ---
chrono_bai24_clean <- chrono_bai24 %>%
  filter(!is.na(PRInorm_10nm_Median),
         !is.na(BAI_2024),
         PRInorm_10nm_Median >= lower_bound,
         PRInorm_10nm_Median <= upper_bound)

# --- 3) Fit linear model & get R² on cleaned data ---
mod <- lm(BAI_2024 ~ PRInorm_10nm_Median, data = chrono_bai24_clean)
r2  <- summary(mod)$r.squared

# Where to place the R² label (top-left of the panel)
x_pos <- min(chrono_bai24_clean$PRInorm_10nm_Median, na.rm = TRUE)
y_pos <- max(chrono_bai24_clean$BAI_2024, na.rm = TRUE)

# --- 4) Plot ---
p <- ggplot(chrono_bai24_clean,
            aes(x = PRInorm_10nm_Median, y = BAI_2024, label = ChronID)) +
  geom_point(size = 2, alpha = 0.75, color = "blue") +
  # If you have many labels, swap the next line for geom_text_repel(size = 3)
  geom_text(vjust = -0.8, size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  annotate("text", x = x_pos, y = y_pos,
           label = sprintf("R² = %.3f", r2),
           hjust = 0, vjust = 1, size = 5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "PRInorm_10nm_Median vs. BAI_2024 (Outliers Removed)",
    x = "PRInorm_10nm_Median",
    y = "BAI_2024"
  )

p

##### Look at scatterplots of just dendrometer trees ####
# = USER SETTINGS =
spec_var  <- "Vogelmann4_1nm_Median"      # <- name of a spectral column in chrono_filt
dendro_var<- "BAI_2024"               # <- name of a dendro column in chrono_filt
reg_type  <- "linear"   # <- "linear", "quadratic", or "exponential"
point_alpha <- 0.7
point_size  <- 1.8

# = PREP =
stopifnot(spec_var %in% names(chrono_bai24),
          dendro_var %in% names(chrono_bai24))

df_xy <- chrono_bai24 %>%
  select(Site, x = all_of(spec_var), y = all_of(dendro_var)) %>%
  tidyr::drop_na()

# = FIT + R² =
r2_text <- ""
fit_line_df <- NULL

if (reg_type == "linear") {
  fit <- lm(y ~ x, data = df_xy)
  r2  <- summary(fit)$r.squared
  r2_text <- sprintf("R² = %.3f (linear)", r2)
  
} else if (reg_type == "quadratic") {
  fit <- lm(y ~ poly(x, 2, raw = TRUE), data = df_xy)
  r2  <- summary(fit)$r.squared
  r2_text <- sprintf("R² = %.3f (quadratic)", r2)
  
} else if (reg_type == "exponential") {
  df_fit <- df_xy %>% filter(y > 0)
  if (nrow(df_fit) < 2) stop("Not enough positive y values for exponential fit.")
  fit <- lm(log(y) ~ x, data = df_fit)
  preds <- exp(predict(fit, newdata = df_fit))
  r2 <- 1 - sum((df_fit$y - preds)^2, na.rm = TRUE) / sum((df_fit$y - mean(df_fit$y, na.rm = TRUE))^2, na.rm = TRUE)
  r2_text <- sprintf("R² = %.3f (exponential)", r2)
  
  xgrid <- data.frame(x = seq(min(df_fit$x, na.rm = TRUE),
                              max(df_fit$x, na.rm = TRUE),
                              length.out = 200))
  fit_line_df <- transform(xgrid, y = exp(predict(fit, newdata = xgrid)))
} else {
  stop("reg_type must be one of: 'linear', 'quadratic', 'exponential'")
}

# = PLOT =
p <- ggplot(df_xy, aes(x = x, y = y, color = Site)) +
  geom_point(alpha = point_alpha, size = point_size) +
  labs(title = paste(spec_var, "vs", dendro_var),
       x = spec_var, y = dendro_var, color = "Site") +
  theme_minimal(base_size = 14)

if (reg_type == "linear") {
  p <- p + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black")
} else if (reg_type == "quadratic") {
  p <- p + geom_smooth(method = "lm", formula = y ~ poly(x, 2, raw = TRUE), se = FALSE, color = "black")
} else if (reg_type == "exponential") {
  p <- p + geom_line(data = fit_line_df, aes(x = x, y = y), inherit.aes = FALSE, color = "black")
}

# Add R² annotation
p <- p + annotate("text", x = Inf, y = Inf, label = r2_text,
                  hjust = 1.02, vjust = 1.5, size = 4.2)

p

## Visusalize using age or DBH
library(dplyr)
library(ggplot2)

# = USER SETTINGS =
spec_var   <- "REPLE_5nm_Median"  # x-axis spectral column in chrono_bai24
dendro_var <- "BAI_2024"               # y-axis dendro column in chrono_bai24
size_by    <- "DBH"                    # choose "Age" or "DBH" (or any numeric column)
point_alpha <- 0.7
size_range  <- c(1.8, 6)               # visual size range of points

# = PREP =
req_cols <- c("Site", spec_var, dendro_var, size_by)
stopifnot(all(req_cols %in% names(chrono_bai24)))

df_xy <- chrono_bai24 %>%
  select(Site,
         x   = all_of(spec_var),
         y   = all_of(dendro_var),
         size_var = all_of(size_by)) %>%
  tidyr::drop_na()

# = PLOT =
ggplot(df_xy, aes(x = x, y = y, color = Site, size = size_var)) +
  geom_point(alpha = point_alpha) +
  scale_size_continuous(name = size_by, range = size_range) +
  labs(
    title = paste(spec_var, "vs", dendro_var),
    x = spec_var, y = dendro_var, color = "Site"
  ) +
  theme_minimal(base_size = 14)

# QQ plot for REPLE_5nm_Median
qqnorm(chrono_bai24$REPLE_5nm_Median, main = "QQ Plot: REPLE_5nm_Median")
qqline(chrono_bai24$REPLE_5nm_Median, col = "red")

# QQ plot for BAI_2024
qqnorm(chrono_bai24$BAI_2024, main = "QQ Plot: BAI_2024")
qqline(chrono_bai24$BAI_2024, col = "red")

par(mfrow = c(1, 1))  # reset layout









############## MULTISPECTRAL ANALYSIS: DO RELATIONSHIPS HOLD? ##################
rds_path <- "G:/HyperspectralUAV/R_outputs/speclib_chronologies/dataframes/chronologies_speclib_bytree_1nm.rds"
# read the spectral library
chronologies_speclib_1nm <- readRDS(rds_path)
## Create new multispectral bands
# INPUT
df_in   <- chronologies_speclib_1nm
nm_min  <- 398
nm_max  <- 999

# --- Detect wavelength columns (either "1nm_###" or just "###") ---
nms <- names(df_in)
pat_pref <- grepl("^1nm_\\d+$", nms)
pat_bare <- grepl("^\\d+$",      nms)

if (any(pat_pref)) {
  nm_cols <- nms[pat_pref]
  nm_vals <- as.integer(sub("^1nm_", "", nm_cols))
} else if (any(pat_bare)) {
  nm_cols <- nms[pat_bare]
  nm_vals <- as.integer(nm_cols)
} else {
  stop("No wavelength columns detected (expected names like '1nm_398'..'1nm_999' or '398'..'999').")
}

# keep only those within desired range and sort
keep <- nm_vals >= nm_min & nm_vals <= nm_max
nm_cols <- nm_cols[keep]
nm_vals <- nm_vals[keep]
ord <- order(nm_vals)
nm_cols <- nm_cols[ord]
nm_vals <- nm_vals[ord]

# map wavelength -> column name
wl2col <- setNames(nm_cols, nm_vals)
avail_nm <- as.integer(names(wl2col))

# --- Sentinel-2A bands (center nm, bandwidth nm) ---
s2A <- data.frame(
  band   = c("B1","B2","B3","B4","B5","B6","B7","B8","B8A","B9"),
  center = c(442.7,492.4,559.8,664.6,704.1,740.5,782.8,832.8,864.7,945.1),
  bw     = c(21,    66,   36,   31,   15,   15,   20,   106,  21,    20),
  stringsAsFactors = FALSE
)

mk_win <- function(center, bw) {
  c(lo = ceiling(center - bw/2), hi = floor(center + bw/2))
}

# compute integer windows per band (clipped to available nm range)
wins <- lapply(seq_len(nrow(s2A)), function(i) {
  w <- mk_win(s2A$center[i], s2A$bw[i])
  w["lo"] <- max(nm_min, w["lo"])
  w["hi"] <- min(nm_max, w["hi"])
  unname(w)
})
names(wins) <- s2A$band

# --- Compute row-wise MEAN for each band window ---
df_out <- df_in
for (b in s2A$band) {
  lo <- wins[[b]][1]; hi <- wins[[b]][2]
  if (!is.finite(lo) || !is.finite(hi) || lo > hi) {
    df_out[[b]] <- NA_real_
    next
  }
  wl_seq <- lo:hi
  wl_seq <- wl_seq[wl_seq %in% avail_nm]
  if (length(wl_seq) == 0) {
    df_out[[b]] <- NA_real_
    next
  }
  sel_cols <- wl2col[as.character(wl_seq)]
  m <- as.matrix(df_in[, sel_cols, drop = FALSE])
  df_out[[b]] <- rowMeans(m, na.rm = TRUE)
}

# RESULT
chronologies_speclib_sentinel <- df_out

# Optional quick checks:
names(chronologies_speclib_sentinel)[(ncol(df_in)+1):ncol(chronologies_speclib_sentinel)]
data.frame(band = names(wins), lo = sapply(wins, `[`, 1), hi = sapply(wins, `[`, 2))

### Create indices using Sentinel bands

# inputs
df <- chronologies_speclib_sentinel
L_SAVI <- 0.5   # common SAVI soil-adjustment (0 to 1); change if you prefer

# helpers
getv <- function(d, nm) if (nm %in% names(d)) d[[nm]] else rep(NA_real_, nrow(d))
safe_div <- function(num, den) { out <- num/den; out[!is.finite(out)] <- NA_real_; out }

# shorthand reflectances (Sentinel-2A band means)
R2  <- getv(df, "B2")   # ~Blue 490
R3  <- getv(df, "B3")   # ~Green 560
R4  <- getv(df, "B4")   # ~Red 665
R5  <- getv(df, "B5")   # ~RE 705
R6  <- getv(df, "B6")   # ~RE 740
R7  <- getv(df, "B7")   # ~RE 783
R8  <- getv(df, "B8")   # ~NIR 833
R8A <- getv(df, "B8A")  # if you later want narrow-NIR variants

# --- indices you requested ---

# NDVI
df$NDVI <- safe_div(R8 - R4, R8 + R4)
# NDRE (use 705 nm red-edge)
df$NDRE <- safe_div(R8 - R5, R8 + R5)
# Red edge CI (CIre = NIR/RE - 1) using 705 nm
df$CIre <- safe_div(R8, R5) - 1
# EVI (standard coefficients)
df$EVI <- 2.5 * safe_div(R8 - R4, R8 + 6*R4 - 7.5*R2 + 1)
# SAVI (with L parameter)
df$SAVI <- (1 + L_SAVI) * safe_div(R8 - R4, R8 + R4 + L_SAVI)
# PSND proxy: (NIR - Blue) / (NIR + Blue); Blue≈B2 (490 nm)
df$PSND_proxy <- safe_div(R8 - R2, R8 + R2)
# PARS proxy: original ~ R746/R513; approximate with B6 (~740) / B3 (~560)
df$PARS_proxy <- safe_div(R6, R3)
# PRI proxy: per your instruction, just use green reflectance
df$PRI_green_proxy <- R3
# REPLE (REP, linear-extrapolation style; 700≈B5, 740≈B6, 780≈B7)
Rre <- (R4 + R7) / 2
df$REPLE <- 700 + 40 * safe_div((Rre - R5), (R6 - R5))
# Also include the standalone ]Red-edge reflectances]
df$RE705  <- R5       # ~705 nm

# result
chronologies_speclib_sentinel_vis <- df

#### Prep correlation matrix df
chrono_multi<- chrono %>%
  # remove columns ending in SD, Q25, Q75, Mean
  select(-matches("(SD|Q25|Q75|Mean)$")) %>%
  # filter out S crown classes
  filter(!CC %in% c("S")) %>%
  # join with crown_metrics and crown_rugosity by TreeID
  left_join(crown_metrics,  by = "TreeID") %>%
  left_join(crown_rugosity, by = "TreeID") %>%
  left_join(crown_areas, by = "TreeID")

chrono_multi <- chrono_multi %>%
  left_join(
    chronologies_speclib_sentinel_vis %>%
      select(Sample, B1:RE705),
    by = c("TreeID" = "Sample")
  )
# output directory
out_dir <- "G:/HyperspectralUAV/R_outputs/modelling/correlations/dendrochronology_correlations"

# file paths
csv_path <- file.path(out_dir, "chrono_multi.csv")
rds_path <- file.path(out_dir, "chrono_multi.rds")
# write files
write.csv(chrono_multi, csv_path, row.names = FALSE)
saveRDS(chrono_multi, rds_path)

##### Do correlation matrix (multi) #####
# output directory
out_dir <- "G:/HyperspectralUAV/R_outputs/modelling/correlations/dendrochronology_correlations"

# --- code to read the .rds back in: hyperspectral, multispectral, structural, dendro all in one
chrono_multi <- readRDS(out_dir, "chrono_multi.rds")

## Spearman correlation matrix
# X: spectral + structural
X <- chrono_multi %>%
  select(B1:RE705) %>%
  select(where(is.numeric))

# Y: dendro
Y <- chrono_multi %>%
  select(BAI_2000:age) %>%
  select(where(is.numeric))

# Spearman correlation matrix (rows = X vars, cols = Y vars)
cor_multi <- cor(as.matrix(X), as.matrix(Y),
                 method = "spearman",
                 use = "pairwise.complete.obs")
cor_df_multi <- as.data.frame(cor_multi)

# Convert correlation matrix into long-form dataframe
cor_long <- cor_multi %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Spectral") %>%
  tidyr::pivot_longer(-Spectral, names_to = "Dendro", values_to = "Correlation") %>%
  mutate(AbsCorr = abs(Correlation))

# --- Top 10 for BAI_2024 ---
top10_BAI_2024 <- cor_long %>%
  filter(Dendro == "BAI_2024") %>%
  arrange(desc(AbsCorr)) %>%
  slice_head(n = 12) %>%
  mutate(
    Correlation = round(Correlation, 3),
    AbsCorr     = round(AbsCorr, 3)
  )

# --- Top 10 for age ---
top10_age <- cor_long %>%
  filter(Dendro == "age") %>%
  arrange(desc(AbsCorr)) %>%
  slice_head(n = 11) %>%
  mutate(
    Correlation = round(Correlation, 3),
    AbsCorr     = round(AbsCorr, 3)
  )

library(DT)
# Send to Viewer
datatable(top10_BAI_2024,
          caption = "Top 10 multispectral vegetation indices correlated with BAI_2024",
          options = list(pageLength = 20, autoWidth = TRUE))

datatable(top10_age,
          caption = "Top 10 multispectral vegetation indices correlated with Age",
          options = list(pageLength = 20, autoWidth = TRUE))

# Export tables
out_dir <- "G:/HyperspectralUAV/R_outputs/modelling/correlations/dendrochronology_correlations"
write.csv(top10_BAI_2024, file.path(out_dir, "top10_BAI_2024.csv"), row.names = FALSE)
write.csv(top10_age,     file.path(out_dir, "top10_age.csv"),     row.names = FALSE)
write.csv(cor_df_multi,     file.path(out_dir, "cor_df_multi.csv"),     row.names = FALSE)

# Optional RDS export
saveRDS(top10_BAI_2024, file.path(out_dir, "top10_BAI_2024.rds"))
saveRDS(top10_age,     file.path(out_dir, "top10_age.rds"))
saveRDS(cor_df_multi,     file.path(out_dir, "cor_df_multi.rds"))

#################### PCA + ORDINATION ANALYSIS #############################
library(dplyr)
### Prep df ####
merge_pair <- function(df, base_id, suffix) {
  ids <- c(base_id, paste0(base_id, suffix))
  if (!all(ids %in% df$Sample)) {
    warning("Skipping merge: not all IDs found: ", paste(ids, collapse = " & "))
    return(df)
  }
  num_cols <- names(df)[sapply(df, is.numeric)]
  merged_vals <- df %>%
    filter(Sample %in% ids) %>%
    summarise(across(all_of(num_cols), ~ mean(.x, na.rm = TRUE)))
  merged_row <- tibble(Sample = base_id) %>% bind_cols(merged_vals)
  df %>%
    filter(!Sample %in% ids) %>%
    bind_rows(merged_row) %>%
    select(names(df))  # preserve original column order
}
# --- 1. Clean up duplicates in chronologies_speclib_1nm ---
chronologies_speclib_1nm_clean <- chronologies_speclib_1nm %>%
  merge_pair("GW_PIRU_14.4_9.5_D", "_6") %>%
  merge_pair("GW_PIRU_12_37.5_D", "_7")

# --- 2. Join onto chrono_multi (214 rows is the base) ---
chrono_joined <- chrono_multi %>%
  left_join(chronologies_speclib_1nm_clean,
            by = c("TreeID" = "Sample"))

# --- Set output directory ---
out_dir <- "G:/HyperspectralUAV/R_outputs/modelling/correlations/dendrochronology_correlations"

# --- Save ---
write.csv(chrono_joined,
          file.path(out_dir, "chrono_multi_speclib.csv"),
          row.names = FALSE)

saveRDS(chrono_joined,
        file.path(out_dir, "chrono_multi_speclib.rds"))

### Get cooking ####
out_dir <- "G:/HyperspectralUAV/R_outputs/modelling/correlations/dendrochronology_correlations"
# --- Read back in --- add lat/long!
chrono_multi_speclib <- readRDS(file.path(out_dir, "chrono_multi_speclib.rds"))
chrono_multi_speclib <- chrono_multi_speclib %>%
  left_join(crown_latlong %>% select(TreeID, Lat, Long),
            by = "TreeID")

# check
names(chrono_multi_speclib)[1:20]   # confirm Lat/Long are present

library(dplyr)
library(vegan)

### PCA on only numeric spectral columns (398–999) ####
suppressPackageStartupMessages({
  library(dplyr)
  library(vegan)
  library(tibble)
  library(ggplot2)
})

# --- 1) PCA on raw spectra (398–999) ---
X <- chrono_multi_speclib %>% select(`398`:`999`)
pca_spectra <- rda(scale(X))

# Variance explained (first 10 PCs)
print(summary(pca_spectra)$cont$importance[2, 1:10])

# --- 2) Eigenvector (loading) curves for PC1 & PC2 ---
loadings_df <- as.data.frame(scores(pca_spectra, choices = 1:2, display = "species")) %>%
  rownames_to_column("Wavelength") %>%
  mutate(Wavelength = as.numeric(Wavelength))

ggplot(loadings_df, aes(x = Wavelength)) +
  geom_line(aes(y = PC1, color = "PC1")) +
  geom_line(aes(y = PC2, color = "PC2")) +
  labs(y = "Loading (eigenvector coefficient)",
       title = "PCA Loadings by Wavelength (398–999 nm)") +
  scale_color_manual(values = c("PC1" = "blue", "PC2" = "red"),
                     name = NULL) +   # removes "colour" as title
  theme_minimal()

# --- 3) envfit: DBH, structure, and coordinates (NAs OK per-variable) ---
Y <- chrono_multi_speclib %>%
  mutate(across(c(DBH, rugosity, Area_m2, zq95, Lat, Long), as.numeric)) %>%
  select(DBH, rugosity, Area_m2, zq95, Lat, Long)

set.seed(1)
fit <- envfit(pca_spectra, Y, permutations = 999)

# --- 4) PCA plot: crowns colored by Site; all arrows (non-sig red, sig black) ---
site_fac <- factor(chrono_multi_speclib$Site)

plot(pca_spectra, display = "sites", type = "n",
     main = "PCA of 1nm Spectral Library (398–999 nm) with Variable Overlays")
points(pca_spectra, display = "sites",
       col = site_fac, pch = 19, cex = 1.2)
legend("topright", legend = levels(site_fac),
       col = seq_along(levels(site_fac)), pch = 19, cex = 1.5, bty = "n")

# Draw ALL variables in red first (non-sig color)
plot(fit, p.max = 1, col = "red")
# Overplot significant ones in black
plot(fit, p.max = 0.05, col = "black")

# --- 5) Tidy envfit results table ---
envfit_df <- as.data.frame(fit$vectors$arrows * sqrt(fit$vectors$r))
envfit_df$Variable <- rownames(envfit_df)
envfit_df$r2   <- fit$vectors$r
envfit_df$pval <- fit$vectors$pvals
envfit_df <- envfit_df %>% select(Variable, PC1 = V1, PC2 = V2, r2, pval)

#--- 6) Export results ---
out_dir <- "G:/HyperspectralUAV/R_outputs/pca/chronologies_pca"

# --- Save ---
write.csv(envfit_df,
          file.path(out_dir, "speclib_pca_envfit.csv"),
          row.names = TRUE)
write.csv(envfit_stats,
          file.path(out_dir, "speclib_pca_envfit_stats.csv"),
          row.names = TRUE)
# save the full PCA object (best for reloading in R) ---
saveRDS(pca_spectra,
        file.path(out_dir, "speclib_pca_object.rds"))

# export variance explained ---
var_explained <- summary(pca_spectra)$cont$importance[2, ]
write.csv(data.frame(PC = names(var_explained),
                     Variance_Explained = var_explained),
          file.path(out_dir, "speclib_pca_variance.csv"),
          row.names = FALSE)
# export PCA site scores (crowns in ordination space)
site_scores <- as.data.frame(scores(pca_spectra, display = "sites"))
site_scores$TreeID <- chrono_multi_speclib$TreeID
site_scores$Site   <- chrono_multi_speclib$Site
write.csv(site_scores,
          file.path(out_dir, "speclib_pca_sites.csv"),
          row.names = FALSE)
# export PCA loadings (wavelength contributions)
loadings <- scores(pca_spectra, display = "species")
loadings_df <- as.data.frame(loadings) %>%
  rownames_to_column("Wavelength")
write.csv(loadings_df,
          file.path(out_dir, "speclib_pca_loadings.csv"),
          row.names = FALSE)

##### PCA on only vegetation indices ARI1 through WBI #####
suppressPackageStartupMessages({
  library(dplyr)
  library(vegan)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})

# --- 1) PCA on vegetation indices (scaled) ---
X <- chrono_multi_speclib %>% select(ARI1_1nm_Median:WBI_15nm_Median)
pca_vi <- rda(scale(X))
#princomp()

# Variance explained (first 10 PCs)
print(summary(pca_vi)$cont$importance[2, 1:10])

# --- 2) Loadings for PC1 & PC2 (indices) ---
vi_loadings <- as.data.frame(scores(pca_vi, choices = 1:2, display = "species")) %>%
  rownames_to_column("Index")

# make a tidy long df and keep top contributors for readability
vi_long <- vi_loadings %>%
  tidyr::pivot_longer(cols = c(PC1, PC2), names_to = "PC", values_to = "Loading") %>%
  group_by(PC) %>%
  slice_max(order_by = abs(Loading), n = 20) %>%   # top 20 by abs loading per PC
  ungroup() %>%
  arrange(PC, Loading)

ggplot(vi_long, aes(x = Loading, y = reorder(Index, Loading), color = PC)) +
  geom_point() +
  facet_wrap(~ PC, scales = "free_y") +
  labs(title = "Top index loadings for PC1 and PC2",
       x = "Loading (eigenvector coefficient)", y = NULL) +
  scale_color_manual(values = c("PC1" = "blue", "PC2" = "red"), name = NULL) +
  theme_minimal(base_size = 12)

# --- 3) envfit with structure & coordinates (NAs OK per-variable) ---
Y <- chrono_multi_speclib %>%
  mutate(across(c(DBH, rugosity, Area_m2, zq95, Lat, Long), as.numeric)) %>%
  select(DBH, rugosity, Area_m2, zq95, Lat, Long)

set.seed(1)
fit_vi <- envfit(pca_vi, Y, permutations = 999)

# --- 4) PCA plot: crowns colored by Site; arrows (non-sig red, sig black) ---
site_fac <- factor(chrono_multi_speclib$Site)

plot(pca_vi, display = "sites", type = "n",
     main = "PCA of Vegetation Indices with Variable Overlays")
points(pca_vi, display = "sites",
       col = site_fac, pch = 19, cex = 1.2)
legend("topright", legend = levels(site_fac),
       col = seq_along(levels(site_fac)), pch = 19, cex = 1.2, bty = "n")

# Draw ALL variables in red first (non-sig color), then overplot sig in black
plot(fit_vi, p.max = 1, col = "red")
plot(fit_vi, p.max = 0.05, col = "black")

# --- 5) Tidy envfit results table ---
envfit_df <- as.data.frame(fit_vi$vectors$arrows * sqrt(fit_vi$vectors$r))
envfit_df$Variable <- rownames(envfit_df)
envfit_df$r2   <- fit_vi$vectors$r
envfit_df$pval <- fit_vi$vectors$pvals
#envfit_df <- envfit_df %>% select(Variable, PC1 = V1, PC2 = V2, r2, pval)

# add a simple stats table too (handy for quick viewing/export)
envfit_stats <- data.frame(
  Variable = rownames(fit_vi$vectors$arrows),
  r2   = fit_vi$vectors$r,
  pval = fit_vi$vectors$pvals,
  stringsAsFactors = FALSE
)

# --- 6) Export results ---
out_dir <- "G:/HyperspectralUAV/R_outputs/pca/chronologies_pca"

# Save envfit tables
write.csv(envfit_df,
          file.path(out_dir, "vi_pca_envfit.csv"),
          row.names = TRUE)
write.csv(envfit_stats,
          file.path(out_dir, "vi_pca_envfit_stats.csv"),
          row.names = TRUE)

# Save the full PCA object (reloadable in R)
saveRDS(pca_vi,
        file.path(out_dir, "vi_pca_object.rds"))

# Export variance explained
var_explained_vi <- summary(pca_vi)$cont$importance[2, ]
write.csv(data.frame(PC = names(var_explained_vi),
                     Variance_Explained = var_explained_vi),
          file.path(out_dir, "vi_pca_variance.csv"),
          row.names = FALSE)

# Export PCA site scores (crowns in ordination space)
vi_site_scores <- as.data.frame(scores(pca_vi, display = "sites"))
vi_site_scores$TreeID <- chrono_multi_speclib$TreeID
vi_site_scores$Site   <- chrono_multi_speclib$Site
write.csv(vi_site_scores,
          file.path(out_dir, "vi_pca_sites.csv"),
          row.names = FALSE)

# Export PCA loadings (index contributions)
write.csv(vi_loadings,
          file.path(out_dir, "vi_pca_loadings.csv"),
          row.names = FALSE)

##### PCA on vegetation indices ARI1 through WBI - using princomp() #######
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(vegan)
  library(ggplot2)
})

# --- 1) PCA on vegetation indices (scaled) using prcomp ---
X <- chrono_multi_speclib %>% select(ARI1_1nm_Median:WBI_15nm_Median)
pca_vi <- prcomp(X, scale. = TRUE, center = TRUE)

# Variance explained (first 10 PCs)
var_explained <- (pca_vi$sdev^2) / sum(pca_vi$sdev^2)
print(var_explained[1:10])

# --- 2) Loadings for PC1 & PC2 (indices) ---
vi_loadings <- as.data.frame(pca_vi$rotation[, 1:2]) %>%
  rownames_to_column("Index") %>%
  rename(PC1 = PC1, PC2 = PC2)

# tidy for plotting top contributors
vi_long <- vi_loadings %>%
  pivot_longer(cols = c(PC1, PC2), names_to = "PC", values_to = "Loading") %>%
  group_by(PC) %>%
  slice_max(order_by = abs(Loading), n = 20) %>%
  ungroup() %>%
  arrange(PC, Loading)

ggplot(vi_long, aes(x = Loading, y = reorder(Index, Loading), color = PC)) +
  geom_point() +
  facet_wrap(~ PC, scales = "free_y") +
  labs(title = "Top index loadings for PC1 and PC2",
       x = "Loading (eigenvector coefficient)", y = NULL) +
  scale_color_manual(values = c("PC1" = "blue", "PC2" = "red"), name = NULL) +
  theme_minimal(base_size = 12)

# --- 3) envfit with structure & coordinates ---
Y <- chrono_multi_speclib %>%
  mutate(across(c(DBH, rugosity, Area_m2, zq95), as.numeric)) %>% # Lat, Long), as.numeric)) %>%
  select(DBH, rugosity, Area_m2, zq95) #Lat, Long)

# envfit requires a vegan-style ordination
# wrap PCA scores into an ordiplot object
pca_sites <- pca_vi$x[, 1:2]   # first 2 PCs
colnames(pca_sites) <- c("PC1", "PC2")
ord <- vegan::ordiplot(pca_sites, choices = c(1, 2), plot = FALSE)

set.seed(1)
fit_vi <- envfit(ord, Y, permutations = 999)

# --- 4) PCA plot: crowns colored by Site + envfit arrows ---
site_fac <- factor(chrono_multi_speclib$Site)

plot(ord, type = "n",
     main = "PCA of Vegetation Indices with Variable Overlays (prcomp)")
points(ord, col = site_fac, pch = 19, cex = 1.2)
legend("topright", legend = levels(site_fac),
       col = seq_along(levels(site_fac)), pch = 19, cex = 1.2, bty = "n")

plot(fit_vi, col = "red", p.max = 1)      # all
plot(fit_vi, col = "black", p.max = 0.05) # sig only

# --- 5) Tidy envfit results table ---
envfit_df <- as.data.frame(fit_vi$vectors$arrows * sqrt(fit_vi$vectors$r))
envfit_df$Variable <- rownames(envfit_df)
envfit_df$r2   <- fit_vi$vectors$r
envfit_df$pval <- fit_vi$vectors$pvals
envfit_df <- envfit_df %>% select(Variable, PC1 = V1, PC2 = V2, r2, pval)

# also keep simple stats table
envfit_stats <- data.frame(
  Variable = rownames(fit_vi$vectors$arrows),
  r2   = fit_vi$vectors$r,
  pval = fit_vi$vectors$pvals,
  stringsAsFactors = FALSE
)

#### PCA without highly correlated indices (test) ######
library(caret)   # for findCorrelation()
library(dplyr)
library(vegan)

# --- 1. Select VI columns ---
X <- chrono_multi_speclib %>% select(ARI1_1nm_Median:WBI_15nm_Median)

# --- 2. Remove highly correlated indices (r > 0.95, adjust if needed) ---
cor_mat <- cor(X, use = "pairwise.complete.obs")
high_cor <- findCorrelation(cor_mat, cutoff = 0.95)   # indices of cols to drop
X_pruned <- X[, -high_cor]

cat("Original:", ncol(X), "indices\n")
cat("After pruning:", ncol(X_pruned), "indices\n")

# --- 3. PCA on pruned indices (scaled) ---
pca_vi_pruned <- rda(scale(X_pruned))

# Variance explained (first 10 PCs)
print(summary(pca_vi_pruned)$cont$importance[2, 1:10])

# --- 4. Loadings for inspection ---
vi_loadings_pruned <- as.data.frame(scores(pca_vi_pruned, choices = 1:2, display = "species")) %>%
  tibble::rownames_to_column("Index")

head(vi_loadings_pruned)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)    # <-- add this
  library(vegan)
  library(ggplot2)
})

# --- 1) envfit on the pruned PCA ---
Y <- chrono_multi_speclib %>%
  mutate(across(c(DBH, rugosity, Area_m2, zq95, Lat, Long), as.numeric)) %>%
  select(DBH, rugosity, Area_m2, zq95, Lat, Long)

set.seed(1)
fit_vi_pruned <- envfit(pca_vi_pruned, Y, permutations = 999)

# --- 2) Dot plot of top loadings for PC1 & PC2 (pruned indices) ---
vi_loadings_pruned <- as.data.frame(scores(pca_vi_pruned, choices = 1:2, display = "species")) %>%
  tibble::rownames_to_column("Index")

vi_long_pruned <- vi_loadings_pruned %>%
  pivot_longer(cols = c(PC1, PC2), names_to = "PC", values_to = "Loading") %>%
  group_by(PC) %>%
  slice_max(order_by = abs(Loading), n = 20) %>%
  ungroup() %>%
  arrange(PC, Loading)

ggplot(vi_long_pruned, aes(x = Loading, y = reorder(Index, Loading), color = PC)) +
  geom_point() +
  facet_wrap(~ PC, scales = "free_y") +
  labs(title = "Top index loadings for PC1 and PC2 (pruned)",
       x = "Loading (eigenvector coefficient)", y = NULL) +
  scale_color_manual(values = c("PC1" = "blue", "PC2" = "red"), name = NULL) +
  theme_minimal(base_size = 12)

# --- 3) PCA plot (sites colored) with ALL arrows (non-sig red, sig black) ---
site_fac <- factor(chrono_multi_speclib$Site)

plot(pca_vi_pruned, display = "sites", type = "n",
     main = "PCA of Pruned Vegetation Indices with Variable Overlays")
points(pca_vi_pruned, display = "sites",
       col = site_fac, pch = 19, cex = 1.1)
legend("topright", legend = levels(site_fac),
       col = seq_along(levels(site_fac)), pch = 19, cex = 1.2, bty = "n")

# draw all variables, then overplot significant ones
plot(fit_vi_pruned, p.max = 1, col = "red",  labels = TRUE)
plot(fit_vi_pruned, p.max = 0.05, col = "black", labels = TRUE)

# scale envfit vectors to PCA biplot space
vec  <- scores(fit_vi_pruned, display = "vectors")
sc   <- ordiArrowMul(vec)
ends <- vec * sc

# significance-based colors
pvals <- fit_vi_pruned$vectors$pvals
cols  <- ifelse(!is.na(pvals) & pvals <= 0.05, "black", "red")

# draw arrows for all variables
arrows(x0 = 0, y0 = 0, x1 = ends[,1], y1 = ends[,2],
       length = 0.06, col = cols, lwd = 1.1)

# add labels slightly beyond arrow tips so they don't sit on the heads
text(ends[,1] * 1.07, ends[,2] * 1.07,
     labels = rownames(ends), col = cols, cex = 0.9)

# --- 4) (optional) tidy envfit results table for quick inspection ---
envfit_df_pruned <- as.data.frame(fit_vi_pruned$vectors$arrows * sqrt(fit_vi_pruned$vectors$r))
envfit_df_pruned$Variable <- rownames(envfit_df_pruned)
envfit_df_pruned$r2   <- fit_vi_pruned$vectors$r
envfit_df_pruned$pval <- fit_vi_pruned$vectors$pvals
envfit_df_pruned <- envfit_df_pruned %>% select(Variable, PC1 = V1, PC2 = V2, r2, pval)

envfit_df_pruned
