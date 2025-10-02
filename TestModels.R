###################### PREP LINES ############################
#setwd("G:/HyperspectralUAV")
library(dplyr)
library(readr)

# Helper to remove 'X' that only precedes numerics or 'nm' tokens, e.g. "X1nm_X748" -> "1nm_748"
remove_x_numeric_nm <- function(colnames_vec) {
  # Remove 'X' at start or after '_' when followed by digits or by digits+'nm'
  gsub("(^|_)X(?=(\\d|\\d*nm))", "\\1", colnames_vec, perl = TRUE)
}

# --- Read CSVs without "X" name-munging ---
# readr::read_csv never forces syntactic names, so numbers are kept as-is
VIstats_speclib <- read_csv(
  "./R_outputs/speclib_dendrometers/dendrometer_canopy_bytreeid-speclibs_MERGED.csv",
  show_col_types = FALSE
)
dendro_metrics <- read_csv(
  "G:/Dendrometers/dendro_metrics.csv",
  show_col_types = FALSE
)
# --- Clean any legacy leading 'X' in column names (e.g., X1nm_X748) ---
names(VIstats_speclib) <- remove_x_numeric_nm(names(VIstats_speclib))

# If a stray column literally named "X" or "...1" exists from prior writes, drop it
VIstats_speclib <- VIstats_speclib %>% select(-any_of("X"))
VIstats_speclib <- VIstats_speclib %>% select(-any_of("...1"))

# --- Join: keep only the 53 Trees that are present in both (inner join) ---
# (If you prefer to keep all 59 from VIstats with NAs for unmatched, switch to left_join)
dendro_spectra <- VIstats_speclib %>%
  inner_join(dendro_metrics, by = "TreeID")

# Write joined dendrometer–spectra dataframe to CSV
write.csv(dendro_spectra,
          "G:/Dendrometers/dendro_spectra_joined.csv",
          row.names = FALSE)
############# READ IN FOR CODE BELOW ########################
library(dplyr)
library(readr)
dendro_spectra <- read.csv("G:/Dendrometers/dendro_spectra_joined.csv",
                           check.names = FALSE)

##################################################################################
#-------REMOVE OTHER STATS FROM DF BESIDES MEAN AND MEDIAN-----------------
dendro_spectra_clean <- dendro_spectra %>%
  select(-matches("(Q25|Q75|SD)$"))
# Remove raw wavelength columns (those starting with 1nm, 5nm, 10nm, or 15nm)
dendro_spectra_clean <- dendro_spectra_clean %>%
  select(-matches("^(1nm|5nm|10nm|15nm)"))

# Write to CSV
write_csv(dendro_spectra_clean, "G:/Dendrometers/dendro_spectra_clean.csv")
# Read back in (column names preserved, no munging)
dendro_spectra_clean <- read_csv("G:/Dendrometers/dendro_spectra_clean.csv",
                                 show_col_types = FALSE)
################## SAMPLE SCATTERPLOTS BETWEEN ZG_DAY_# vs. INDEX #####################
# USER SETTINGS
chosen_vi   <- "Boochs_1nm_Median"    # Vegetation index column name
chosen_zg   <- "zg_day_7"             # zg_day column name
overlay     <- "none"                 # Options: "loess", "gam", "binned", "robust", "none"
start_x_at  <- 0                      # Where to start x-axis (e.g., 4 to match earlier plots)
show_pearson <- TRUE                  # Add Pearson r label?
pearson_label_text <- "Spearman r = -0.76"  # Custom text for subtitle

# Custom titles/labels
plot_title  <- "Boochs Index vs Shrinkage Days (7-day window)"
x_axis_lab  <- "Total Shrinkage Days (7-day pre-flight window)"
y_axis_lab  <- "Boochs Index (1nm Median )"

# ==
# FUNCTION
# ==
make_vi_zg_plot <- function(
    data,
    chosen_vi,
    chosen_zg,
    overlay = c("loess", "gam", "binned", "robust", "none"),
    start_x_at = 4,
    show_pearson = FALSE,
    pearson_label_text = "r = ...",
    plot_title = NULL,
    x_axis_lab = NULL,
    y_axis_lab = NULL
) {
  overlay <- match.arg(overlay)
  
  if (!all(c(chosen_vi, chosen_zg) %in% names(data))) {
    stop("Chosen columns not found. Check 'chosen_vi' and/or 'chosen_zg'.")
  }
  
  plot_df <- data %>%
    select(all_of(c(chosen_vi, chosen_zg))) %>%
    rename(VI = !!sym(chosen_vi), zg = !!sym(chosen_zg)) %>%
    mutate(zg = suppressWarnings(as.numeric(zg))) %>%
    tidyr::drop_na(VI, zg)
  
  p <- ggplot(plot_df, aes(x = zg, y = VI)) +
    geom_point(alpha = 0.6, size = 2,
               position = position_jitter(width = 0.15, height = 0)) +
    scale_x_continuous(limits = c(start_x_at, NA)) +
    labs(
      title = plot_title %||% paste0(chosen_vi, " vs ", chosen_zg),
      subtitle = if (show_pearson) pearson_label_text else NULL,
      x = x_axis_lab %||% chosen_zg,
      y = y_axis_lab %||% chosen_vi
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none")
  
  if (overlay == "loess") {
    p <- p + geom_smooth(method = "loess", se = TRUE,
                         color = "black", linetype = "dashed", linewidth = 0.8)
  } else if (overlay == "gam") {
    p <- p + stat_smooth(method = "gam", formula = y ~ s(x, k = 6),
                         se = TRUE, color = "black", linetype = "dashed", linewidth = 0.8)
  } else if (overlay == "binned") {
    summ <- plot_df %>%
      group_by(zg) %>%
      summarise(
        med = median(VI, na.rm = TRUE),
        q1  = quantile(VI, 0.25, na.rm = TRUE),
        q3  = quantile(VI, 0.75, na.rm = TRUE),
        .groups = "drop"
      )
    p <- p +
      geom_line(data = summ, aes(x = zg, y = med), color = "black") +
      geom_pointrange(data = summ, aes(x = zg, y = med, ymin = q1, ymax = q3),
                      color = "black", size = 0.3)
  } else if (overlay == "robust") {
    p <- p + stat_smooth(method = MASS::rlm, se = TRUE,
                         color = "black", linetype = "dashed", linewidth = 0.8)
  }
  
  return(p)
}

# ===
# EXAMPLE CALL
# ===
make_vi_zg_plot(
  data = dendro_spectra_clean,
  chosen_vi = chosen_vi,
  chosen_zg = chosen_zg,
  overlay = overlay,
  start_x_at = start_x_at,
  show_pearson = show_pearson,
  pearson_label_text = pearson_label_text,
  plot_title = plot_title,
  x_axis_lab = x_axis_lab,
  y_axis_lab = y_axis_lab
)

############## GENERATE CORRELATION MATRIX + SUMMARY TABLE FOR ZG_DAY COLUMNS ##########
# --- Packages ---
library(dplyr)
library(stringr)
library(readr)
library(tibble)
library(purrr)

# INPUT: dendro_spectra_clean (53 rows)

# 1) Identify columns (only up to 94)
zg_day_cols <- names(dendro_spectra_clean) %>%
  keep(~ str_detect(., "^zg_day_\\d+$")) %>%                 # only zg_day columns
  keep(~ as.integer(str_remove(., "^zg_day_")) <= 94) %>%    # restrict to <= 94
  sort() 

# Spectral *index* features: keep mean/median summaries only (not raw wavelengths)
# Adjust the pattern if your naming differs.
vi_feature_cols <- names(dendro_spectra_clean) %>%
  keep(~ str_detect(., "(?i)_(Mean|Median)$")) %>%
  setdiff(zg_day_cols)

# 2) Build correlation matrix (rows = index features, cols = zg_day lags)
#    Each cell = Spearman rho
cor_mat <- sapply(zg_day_cols, function(zg) {
  sapply(vi_feature_cols, function(vi) {
    cor(dendro_spectra_clean[[zg]], dendro_spectra_clean[[vi]],
        method = "spearman", use = "complete.obs")
  })
}, simplify = "matrix")

# Make rows = indices, cols = lags
cor_mat <- t(cor_mat)             # currently zg as rows, vi as cols
cor_mat <- t(cor_mat)             # flip back: vi as rows, zg as cols (double t ensures)
rownames(cor_mat) <- vi_feature_cols
colnames(cor_mat) <- zg_day_cols

# Also as a data frame with an explicit feature column
zg_vi_corr_matrix <- as_tibble(cor_mat, rownames = "spectral_feature")

# 3) Summary table: for each index, grab the strongest |rho| across all lags
summary_list <- apply(cor_mat, 1, function(rhos_by_lag) {
  idx_max <- which.max(abs(rhos_by_lag))
  best_rho <- rhos_by_lag[idx_max]
  best_lag_name <- names(rhos_by_lag)[idx_max]
  best_lag_days <- suppressWarnings(as.integer(sub("^zg_day_", "", best_lag_name)))
  
  list(
    spectral_feature = NA_character_,   # filled later
    max_abs_rho      = abs(best_rho),
    rho_signed       = best_rho,
    sign             = ifelse(best_rho >= 0, "+", "–"),
    best_zg_var      = best_lag_name,
    lag_days_back    = best_lag_days
  )
})
zg_vi_corr_summary <- bind_rows(summary_list) |>
  mutate(spectral_feature = vi_feature_cols, .before = 1) |>
  arrange(desc(max_abs_rho))

# 4) (Optional) Write outputs for reuse
dir.create("./R_outputs/modelling/correlations", recursive = TRUE, showWarnings = FALSE)
write_csv(zg_vi_corr_matrix,
          "./R_outputs/modelling/correlations/zg_vi_corr_matrix.csv")
write_csv(zg_vi_corr_summary,
          "./R_outputs/modelling/correlations/zg_vi_corr_summary.csv")

# 5) Quick peek
print(zg_vi_corr_summary, n = 10)

################## FIND MAX CORREALTION FOR EVERY ZG_DAY (UP TO 94) ###################
# 1) Identify lag columns (zg_day_1 .. zg_day_94), sorted numerically
lag_cols <- names(zg_vi_corr_matrix) %>%
  keep(~ str_detect(., "^zg_day_\\d+$")) %>%
  keep(~ as.integer(str_remove(., "^zg_day_")) <= 94) %>%
  (\(x) x[order(as.integer(str_remove(x, "^zg_day_")))])()

# 2) For each lag, find the index with the max |rho|
zg_lag_winners <- map_dfr(lag_cols, function(lag_col) {
  rhos <- zg_vi_corr_matrix[[lag_col]]
  i_max <- which.max(abs(rhos))
  best_rho <- rhos[i_max]
  
  tibble(
    zg_day        = lag_col,
    lag_days_back = as.integer(str_remove(lag_col, "^zg_day_")),
    max_abs_rho   = abs(best_rho),
    rho_signed    = best_rho,
    sign          = ifelse(best_rho >= 0, "+", "–"),
    best_index    = zg_vi_corr_matrix$spectral_feature[i_max]
  )
}) %>%
  arrange(lag_days_back)

# (Optional) Save it
write_csv(zg_lag_winners, "./R_outputs/modelling/correlations/zg_lag_winners.csv")

#--PLOT THIS HOE OUT--
library(ggplot2)
library(dplyr)
library(stringr)
library(readr)

zg_lag_winners <- read_csv("./R_outputs/modelling/correlations/zg_lag_winners.csv",
                           show_col_types = FALSE)

# Start from zg_lag_winners you already created
plot_data <- zg_lag_winners %>%
  mutate(index_family = str_remove(best_index, "_.*$"))

# Explicit colors for selected families
custom_colors <- c(
  "Boochs"     = "#4C78A8",  # medium blue
  "CARI"       = "#F58518",  # warm orange
  "PRI"        = "#55A868",  # light/soft green
  "DPI"        = "#E45756",  # coral red
  "Datt3"      = "#B279A2",  # muted purple
  "TCARIOSAVI" = "#72B7B2"   # seafoam teal
)

p <- ggplot(plot_data, aes(x = lag_days_back,
                           y = max_abs_rho,
                           fill = index_family)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = custom_colors, drop = TRUE) +  # drop unused from legend
  scale_x_continuous(limits = c(4, NA)) +
  labs(
    x = "Days prior to UAV flight (N-day shrinkage window)",
    y = "Absolute Value of Max Spearman ρ",
    fill = "Index",
    title = "Strongest spectral–shrinkage correlation by lag"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

# Optional LOESS line
p + geom_smooth(aes(group = 1), method = "loess",
                color = "black", se = FALSE, linetype = "dashed")

################## TRACK BEST PERFORMING INDICES OVER TIME ###################
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readr)

# --- Inputs assumed: zg_vi_corr_matrix (with 'spectral_feature' + zg_day_1...zg_day_94) ---

# 1) Pick the exact variants you want to plot
target_features <- c(
  "Boochs_1nm_Median",
  "PRI_15nm_Median",
  "CARI_1nm_Mean",
  "Datt3_5nm_Mean",
  "DPI_1nm_Median",
  "TCARIOSAVI_1nm_Median"
)

# 2) Build index-family palette (your colors)
family_colors <- c(
  "Boochs"     = "#4C78A8",  # medium blue
  "CARI"       = "#F58518",  # warm orange
  "PRI"        = "#55A868",  # light/soft green
  "DPI"        = "#E45756",  # coral red
  "Datt3"      = "#B279A2",  # muted purple
  "TCARIOSAVI" = "#72B7B2"   # seafoam teal
)

# 3) Reshape: keep only your targets, pivot lags long, compute |rho|, add lag integer + family
plot_lines <- zg_vi_corr_matrix %>%
  filter(spectral_feature %in% target_features) %>%
  pivot_longer(
    cols = starts_with("zg_day_"),
    names_to = "zg_day",
    values_to = "rho"
  ) %>%
  mutate(
    lag_days_back = as.integer(str_remove(zg_day, "^zg_day_")),
    abs_rho       = abs(rho),
    index_family  = str_remove(spectral_feature, "_.*$")
  ) %>%
  filter(lag_days_back >= 4, lag_days_back <= 94)

# 4) Plot: one line per *variant*, colored by *family*; start x at 4
p_lines <- ggplot(
  plot_lines,
  aes(x = lag_days_back, y = abs_rho,
      group = spectral_feature,
      color = index_family)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.6) +
  scale_color_manual(values = family_colors, drop = TRUE) +
  scale_x_continuous(limits = c(4, 94)) +
  labs(
    title = "Absolute correlation profiles for select indices across pre-flight shrinkage windows",
    x = "Days prior to UAV flight (N-day shrinkage window)",
    y = "Absolute Value of Spearman ρ",
    color = "Index"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

p_lines



################# TEST INDICES AGAINST NEW DENDRO VARIABLES #########################
library(dplyr)
library(readr)

dendro_spectra_clean <- read_csv("G:/Dendrometers/dendro_spectra_clean.csv",
                                 show_col_types = FALSE)

# --- Define columns ---
spectral_cols <- 2:777
dendro_metrics <- c("cTWD", "cTWD_drone", "TWD_drone", "cshrink")

# --- Extract variable names ---
spectral_vars <- names(dendro_spectra_clean)[spectral_cols]

# --- Build correlation matrix for all 4 at once (like before) ---
dm_vi_corr_matrix <- sapply(dendro_metrics, function(dm_var) {
  sapply(spectral_vars, function(spec_var) {
    cor(dendro_spectra_clean[[dm_var]], dendro_spectra_clean[[spec_var]],
        method = "spearman", use = "complete.obs")
  })
})

dm_vi_corr_matrix <- as_tibble(dm_vi_corr_matrix, rownames = "spectral_feature")
write_csv(dm_vi_corr_matrix, "./R_outputs/modelling/correlations/dm_vi_corr_matrix.csv")

# --- Split into 4 separate correlation data frames ---
corr_cTWD <- dm_vi_corr_matrix %>%
  select(spectral_feature, cTWD) %>%
  mutate(abs_cTWD = abs(cTWD))

corr_cTWD_drone <- dm_vi_corr_matrix %>%
  select(spectral_feature, cTWD_drone) %>%
  mutate(abs_cTWD_drone = abs(cTWD_drone))

corr_TWD_drone <- dm_vi_corr_matrix %>%
  select(spectral_feature, TWD_drone) %>%
  mutate(abs_TWD_drone = abs(TWD_drone))

corr_cshrink <- dm_vi_corr_matrix %>%
  select(spectral_feature, cshrink) %>%
  mutate(abs_cshrink = abs(cshrink))

# --- (Optional) Save them ---
#write_csv(dm_vi_corr_matrix, "./R_outputs/modelling/correlations/dm_vi_corr_matrix.csv")
write_csv(corr_cTWD,       "./R_outputs/modelling/correlations/corr_cTWD.csv")
write_csv(corr_cTWD_drone, "./R_outputs/modelling/correlations/corr_cTWD_drone.csv")
write_csv(corr_TWD_drone,  "./R_outputs/modelling/correlations/corr_TWD_drone.csv")
write_csv(corr_cshrink,    "./R_outputs/modelling/correlations/corr_cshrink.csv")

# Same dataset, but removing MEAN index values
# --- Create new df removing spectral features ending with "Mean" ---
dm_vi_corr_matrix_noMean <- dm_vi_corr_matrix %>%
  filter(!str_detect(spectral_feature, "Mean$"))

# --- Split into 4 separate correlation data frames (noMean) ---
corr_cTWD_noMean <- dm_vi_corr_matrix_noMean %>%
  select(spectral_feature, cTWD) %>%
  mutate(abs_cTWD = abs(cTWD))

corr_cTWD_drone_noMean <- dm_vi_corr_matrix_noMean %>%
  select(spectral_feature, cTWD_drone) %>%
  mutate(abs_cTWD_drone = abs(cTWD_drone))

corr_TWD_drone_noMean <- dm_vi_corr_matrix_noMean %>%
  select(spectral_feature, TWD_drone) %>%
  mutate(abs_TWD_drone = abs(TWD_drone))

corr_cshrink_noMean <- dm_vi_corr_matrix_noMean %>%
  select(spectral_feature, cshrink) %>%
  mutate(abs_cshrink = abs(cshrink))

# --- (Optional) Save them ---
write_csv(dm_vi_corr_matrix_noMean,   "./R_outputs/modelling/correlations/dm_vi_corr_matrix_noMean.csv")
write_csv(corr_cTWD_noMean,           "./R_outputs/modelling/correlations/corr_cTWD_noMean.csv")
write_csv(corr_cTWD_drone_noMean,     "./R_outputs/modelling/correlations/corr_cTWD_drone_noMean.csv")
write_csv(corr_TWD_drone_noMean,      "./R_outputs/modelling/correlations/corr_TWD_drone_noMean.csv")
write_csv(corr_cshrink_noMean,        "./R_outputs/modelling/correlations/corr_cshrink_noMean.csv")

###### Create tables ###########
# ==
# USER SETTINGS
# ==
top_n <- 10
out_docx <- "G:/Dendrometers/top_indices_all_metrics.docx"

library(dplyr)
library(stringr)
library(flextable)
library(officer)

# ==
# HELPERS
# ==
# Build a top-N flextable for a single corr df, with 3-dec rounding and no leading zero
make_top_corr_flextable <- function(corr_df, top_n = 5) {
  stopifnot(is.data.frame(corr_df), "spectral_feature" %in% names(corr_df))
  abs_cols   <- grep("^abs_", names(corr_df), value = TRUE)
  metric_col <- setdiff(names(corr_df), c("spectral_feature", abs_cols))
  if (length(metric_col) != 1) stop("Could not uniquely identify the metric column.")
  metric_col <- metric_col[[1]]
  
  top_tbl <- corr_df %>%
    transmute(
      Index        = str_remove(spectral_feature, "_.*$"),
      Resample_nm  = str_match(spectral_feature, "^[^_]+_([^_]+)")[, 2],
      Correlation  = .data[[metric_col]],
      Abs_Corr     = abs(Correlation)
    ) %>%
    arrange(desc(Abs_Corr)) %>%
    slice_head(n = top_n) %>%
    # format to 3 decimals, drop leading zero (e.g., 0.538 -> .538)
    mutate(
      Correlation = sub("^0\\.", ".", sprintf("%.3f", Correlation)),
      Abs_Corr    = sub("^0\\.", ".", sprintf("%.3f", Abs_Corr))
    )
  
  caption_txt <- sprintf("Top %d spectral indices for %s (Spearman ρ)", top_n, metric_col)
  
  ft <- flextable(top_tbl)
  ft <- set_caption(ft, caption = caption_txt)
  ft <- theme_box(ft)         # clear borders around all cells
  ft <- align(ft, align = "center", part = "all")
  ft <- autofit(ft)
  ft
}

# ==
# BUILD TABLES
# =
ft_cshrink     <- make_top_corr_flextable(corr_cshrink_noMean,     top_n)
ft_cTWD        <- make_top_corr_flextable(corr_cTWD_noMean,        top_n)
ft_cTWD_drone  <- make_top_corr_flextable(corr_cTWD_drone_noMean,  top_n)
ft_TWD_drone   <- make_top_corr_flextable(corr_TWD_drone_noMean,   top_n)

# (Optional) view in RStudio Viewer; you can also copy-paste to Google Docs directly
ft_cshrink
ft_cTWD
ft_cTWD_drone
ft_TWD_drone

# =
# EXPORT ALL TABLES INTO ONE DOCX
# =
doc <- read_docx()

doc <- body_add_par(doc, "Top Indices by Correlation — cshrink", style = "heading 1")
doc <- body_add_flextable(doc, ft_cshrink)
doc <- body_add_par(doc, "")  # spacer

doc <- body_add_par(doc, "Top Indices by Correlation — cTWD", style = "heading 1")
doc <- body_add_flextable(doc, ft_cTWD)
doc <- body_add_par(doc, "")

doc <- body_add_par(doc, "Top Indices by Correlation — cTWD_drone", style = "heading 1")
doc <- body_add_flextable(doc, ft_cTWD_drone)
doc <- body_add_par(doc, "")

doc <- body_add_par(doc, "Top Indices by Correlation — TWD_drone", style = "heading 1")
doc <- body_add_flextable(doc, ft_TWD_drone)

print(doc, target = out_docx)

##### Test Normality + Scatterplot ####
dendro_spectra_clean <- read_csv("G:/Dendrometers/dendro_spectra_clean.csv",
                                 show_col_types = FALSE)
# ==
# USER SETTINGS (Diagnostics), test linearity

chosen_x <- "Datt3_5nm_Median"   # any column in dendro_spectra_clean
chosen_y <- "cshrink"            # any column in dendro_spectra_clean

# ==
# CODE (Diagnostics)
# ==
library(dplyr)
library(broom)

# optional tests if the package exists
have_lmtest <- requireNamespace("lmtest", quietly = TRUE)

# Build clean analysis df
df_xy <- dendro_spectra_clean %>%
  select(all_of(c(chosen_x, chosen_y))) %>%
  rename(x = !!chosen_x, y = !!chosen_y) %>%
  mutate(
    x = suppressWarnings(as.numeric(x)),
    y = suppressWarnings(as.numeric(y))
  ) %>%
  tidyr::drop_na()

# Quick prechecks
issues <- character(0)
if (!is.numeric(df_xy$x) || !is.numeric(df_xy$y)) issues <- c(issues, "Non-numeric values detected.")
if (nrow(df_xy) < 10)                               issues <- c(issues, "Fewer than 10 complete cases.")
if (dplyr::n_distinct(df_xy$x) < 8)                issues <- c(issues, "Too few unique x values (< 8).")
if (var(df_xy$x) == 0 || var(df_xy$y) == 0)        issues <- c(issues, "Zero variance in x or y.")

lm_fit <- NULL
reset_p <- NA_real_
bp_p <- NA_real_
shapiro_p <- NA_real_
r2_val <- NA_real_

cat("---- Linear model quick checks ----\n")
if (length(issues)) {
  cat("Precheck issues: ", paste(issues, collapse = " | "), "\n", sep = "")
} else {
  lm_fit <- lm(y ~ x, data = df_xy)
  gl <- glance(lm_fit)
  r2_val <- gl$r.squared
  
  if (have_lmtest) {
    reset_p  <- tryCatch(lmtest::resettest(lm_fit)$p.value, error = function(e) NA_real_)
    bp_p     <- tryCatch(lmtest::bptest(lm_fit)$p.value,     error = function(e) NA_real_)
  }
  resids <- resid(lm_fit)
  if (length(resids) >= 3 && length(resids) <= 5000) {
    shapiro_p <- tryCatch(shapiro.test(resids)$p.value, error = function(e) NA_real_)
  }
  
  if (!is.na(reset_p))   cat(sprintf("RESET nonlinearity p = %.3f\n", reset_p))
  if (!is.na(bp_p))      cat(sprintf("Breusch–Pagan p = %.3f\n", bp_p))
  if (!is.na(shapiro_p)) cat(sprintf("Shapiro–Wilk (residuals) p = %.3f\n", shapiro_p))
  if (is.finite(r2_val)) cat(sprintf("Model R² = %.3f\n", r2_val))
  if (!is.na(reset_p) && reset_p < 0.05) cat("⚠️  RESET suggests nonlinearity (consider LOESS/GAM/binned summaries).\n")
  if (!is.na(bp_p) && bp_p < 0.05)       cat("⚠️  Heteroskedasticity detected (consider robust SE or transform).\n")
}
cat("----------------------------------\n")

# Scatterplot

# USER SETTINGS (Plot)

library(ggplot2)
show_r2  <- TRUE     # show R^2 in subtitle (uses r2_val from diagnostics if available)
add_lm   <- TRUE     # add linear regression line if diagnostics passed
add_loess_reference <- FALSE  # add dashed LOESS for visual nonlinearity check

plot_title <- "Datt3 vs. cshrink"         # NULL -> auto "chosen_y vs chosen_x"
x_lab <- "Datt3"                           # NULL -> chosen_x
y_lab <- NULL                              # NULL -> chosen_y

point_alpha <- 0.8
point_size  <- 3

library(ggplot2)

# Ensure df_xy exists from Chunk A
stopifnot(exists("df_xy"))

# Compose labels safely without %||%
final_title <- if (is.null(plot_title)) paste0(chosen_y, " vs ", chosen_x) else plot_title
final_xlab  <- if (is.null(x_lab)) chosen_x else x_lab
final_ylab  <- if (is.null(y_lab)) chosen_y else y_lab

subtitle_txt <- NULL
if (show_r2 && exists("r2_val") && is.finite(r2_val)) {
  subtitle_txt <- sprintf("R² = %.3f", r2_val)
  subtitle_txt <- sub("^R² = 0\\.", "R² = .", subtitle_txt)  # optional: drop leading zero
}

p <- ggplot(df_xy, aes(x = x, y = y)) +
  geom_point(alpha = point_alpha, size = point_size) +
  labs(
    title = final_title,
    subtitle = subtitle_txt,
    x = final_xlab,
    y = final_ylab
  ) +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, NA))   # zooms without filtering data   # lock y-axis lower bound at 0

if (add_loess_reference) {
  p <- p + geom_smooth(method = "loess", se = TRUE,
                       color = "gray40", linetype = "dashed", linewidth = 0.8)
}

if (add_lm && exists("lm_fit") && !is.null(lm_fit) && length(issues) == 0) {
  p <- p + geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.9)
}

p

lm_quad <- lm(cTWD_drone ~ poly(Datt3_5nm_Median, 2), data = dendro_spectra_clean)
summary(lm_quad)

##### Scatterplot and Trend Lines ######
dendro_spectra_clean <- read_csv("G:/Dendrometers/dendro_spectra_clean.csv",
                                 show_col_types = FALSE)
# 
# USER SETTINGS
# 
chosen_x <- "Vogelmann3_5nm_Median"
chosen_y <- "TWD_drone"

point_alpha <- 0.7
point_size  <- 2.5

# ===
# PREP
# ===
library(dplyr)
library(ggplot2)
library(broom)
library(mgcv)

df_xy <- dendro_spectra_clean %>%
  select(all_of(c(chosen_x, chosen_y))) %>%
  rename(x = !!chosen_x, y = !!chosen_y) %>%
  tidyr::drop_na() %>%
  mutate(x = as.numeric(x),
         y = as.numeric(y)) %>%
  filter(y > 0)   # needed for log/exponential models

# ===
# MODELS
# ===

# Linear
fit_lin <- lm(y ~ x, data = df_xy)

# Quadratic polynomial
fit_quad <- lm(y ~ poly(x, 2, raw = TRUE), data = df_xy)

# Exponential decay: y = a * exp(-b*x)
start_a <- max(df_xy$y, na.rm = TRUE)
fit_exp <- tryCatch(
  nls(y ~ a * exp(-b * x),
      data = df_xy,
      start = list(a = start_a, b = 1)),
  error = function(e) NULL
)

# GAM smoother
fit_gam <- gam(y ~ s(x, k = 5), data = df_xy)

# ==
# R² CALCULATIONS
# ===
rsq <- list()

rsq$linear   <- glance(fit_lin)$r.squared
rsq$quadratic <- glance(fit_quad)$r.squared

if (!is.null(fit_exp)) {
  pred_exp <- predict(fit_exp)
  ss_res <- sum((df_xy$y - pred_exp)^2)
  ss_tot <- sum((df_xy$y - mean(df_xy$y))^2)
  rsq$exponential <- 1 - ss_res/ss_tot
} else {
  rsq$exponential <- NA
}

pred_gam <- predict(fit_gam)
ss_res <- sum((df_xy$y - pred_gam)^2)
ss_tot <- sum((df_xy$y - mean(df_xy$y))^2)
rsq$gam <- 1 - ss_res/ss_tot

print("R² values:")
print(rsq)

# ===
# USER SETTING: which model to plot
# choices: "linear", "quadratic", "exponential", "gam"
# ===
model_choice <- "exponential"

# ==
# PLOT (single selected model + R² on plot)
# ===
# base scatter
p <- ggplot(df_xy, aes(x = x, y = y)) +
  geom_point(alpha = point_alpha, size = point_size) +
  theme_minimal(base_size = 14) +
  coord_cartesian(ylim = c(0, NA))

# choose layer + R²
rsq_selected <- NA_real_

if (model_choice == "linear") {
  p <- p + geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 1,
                       formula = y ~ x)
  rsq_selected <- summary(fit_lin)$r.squared
  
} else if (model_choice == "quadratic") {
  p <- p + geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1,
                       formula = y ~ poly(x, 2, raw = TRUE))
  rsq_selected <- summary(fit_quad)$r.squared
  
} else if (model_choice == "gam") {
  p <- p + geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = FALSE,
                       color = "darkgreen", linewidth = 1)
  pred_gam <- predict(fit_gam)
  ss_res <- sum((df_xy$y - pred_gam)^2)
  ss_tot <- sum((df_xy$y - mean(df_xy$y))^2)
  rsq_selected <- 1 - ss_res/ss_tot
  
} else if (model_choice == "exponential") {
  if (!is.null(fit_exp)) {
    newdat <- data.frame(x = seq(min(df_xy$x), max(df_xy$x), length.out = 200))
    newdat$y <- predict(fit_exp, newdata = newdat)
    p <- p + geom_line(data = newdat, aes(x = x, y = y),
                       color = "purple", linewidth = 1)
    pred_exp <- predict(fit_exp)
    ss_res <- sum((df_xy$y - pred_exp)^2)
    ss_tot <- sum((df_xy$y - mean(df_xy$y))^2)
    rsq_selected <- 1 - ss_res/ss_tot
  } else {
    message("Exponential model failed to fit; falling back to linear.")
    p <- p + geom_smooth(method = "lm", se = FALSE, color = "blue", linewidth = 1,
                         formula = y ~ x)
    rsq_selected <- summary(fit_lin)$r.squared
    model_choice <- "linear"
  }
}

# labels (with R²)
title_txt <- paste(chosen_x, "vs", chosen_y)
subtitle_txt <- if (is.finite(rsq_selected)) sub("^0\\.", ".", sprintf("R² = %.3f", rsq_selected)) else NULL

p <- p +
  labs(
    title = "Vogelmann3 vs. Tree Water Deficit (TWD) on Flight Day",
    subtitle = subtitle_txt,
    x = "Vogelmann3",
    y = "TWD"
  )

p

################### TEST DENDRO METRICS IN MULTISPECTRAL SPACE - SENTINEL2A ###################
# Start here
dendro_spectra <- read.csv("G:/Dendrometers/dendro_spectra_joined.csv",
                           check.names = FALSE)
# Create Sentinel-2A bands
# user choices
s2_version <- "A"     # "A", "B", or "avg"
input_df   <- dendro_spectra
nm_prefix  <- "1nm_"
nm_min     <- 398
nm_max     <- 999

# band definitions (nm)
s2A <- data.frame(
  band = c("B1","B2","B3","B4","B5","B6","B7","B8","B8A","B9"),
  center = c(442.7,492.4,559.8,664.6,704.1,740.5,782.8,832.8,864.7,945.1),
  bw     = c(21,   66,   36,   31,   15,   15,   20,   106,  21,   20),
  stringsAsFactors = FALSE
)
s2B <- data.frame(
  band = c("B1","B2","B3","B4","B5","B6","B7","B8","B8A","B9"),
  center = c(442.2,492.1,559.0,664.9,703.8,739.1,779.7,832.9,864.0,943.2),
  bw     = c(21,   66,   36,   31,   16,   15,   20,   106,  22,   21),
  stringsAsFactors = FALSE
)

# helper to make integer window
mk_win <- function(center, bw) {
  c(lo = ceiling(center - bw/2), hi = floor(center + bw/2))
}

# build windows per band
bands <- s2A$band
wins <- lapply(bands, function(b){
  Aw <- mk_win(s2A[s2A$band==b, "center"], s2A[s2A$band==b, "bw"])
  Bw <- mk_win(s2B[s2B$band==b, "center"], s2B[s2B$band==b, "bw"])
  if (s2_version == "A") {
    w <- Aw
  } else if (s2_version == "B") {
    w <- Bw
  } else {
    w <- c(lo = min(Aw["lo"], Bw["lo"]), hi = max(Aw["hi"], Bw["hi"]))  # union
  }
  # clip to available domain
  w["lo"] <- max(nm_min, w["lo"])
  w["hi"] <- min(nm_max, w["hi"])
  unname(w)  # returns numeric length-2
})
names(wins) <- bands

# identify available 1-nm columns
all_cols <- colnames(input_df)
is_1nm   <- grepl(paste0("^", nm_prefix, "\\d+$"), all_cols)
nm_cols  <- all_cols[is_1nm]
nm_vals  <- as.integer(sub(paste0("^", nm_prefix), "", nm_cols))
wl2col   <- setNames(nm_cols, nm_vals)
avail_nm <- sort(as.integer(names(wl2col)))

# compute and append stats
out_df <- input_df
for (b in bands) {
  lo <- wins[[b]][1]; hi <- wins[[b]][2]
  if (!is.finite(lo) || !is.finite(hi) || lo > hi) {
    out_df[[paste0(b, "_mean")]]   <- NA_real_
    out_df[[paste0(b, "_median")]] <- NA_real_
    next
  }
  wl_seq <- lo:hi
  wl_seq <- wl_seq[wl_seq %in% avail_nm]
  if (length(wl_seq) == 0) {
    out_df[[paste0(b, "_mean")]]   <- NA_real_
    out_df[[paste0(b, "_median")]] <- NA_real_
    next
  }
  sel_cols <- wl2col[as.character(wl_seq)]
  m <- as.matrix(out_df[, sel_cols, drop = FALSE])
  out_df[[paste0(b, "_mean")]]   <- rowMeans(m, na.rm = TRUE)
  out_df[[paste0(b, "_median")]] <- apply(m, 1, function(x){
    x <- x[!is.na(x)]
    if (length(x)==0) NA_real_ else median(x)
  })
}

dendro_spectra_sentinel <- out_df

# optional: quick check of windows actually used
data.frame(
   band = bands,
   lo = sapply(wins, `[`, 1),
   hi = sapply(wins, `[`, 2)
 )

# Generate indices for Sentinel

# user inputs
df <- dendro_spectra_sentinel   # input with B#_mean / B#_median columns

# helpers
safe_div <- function(num, den) {
  out <- num / den
  out[!is.finite(out)] <- NA_real_
  out
}
getv <- function(d, nm) {
  if (nm %in% names(d)) d[[nm]] else rep(NA_real_, nrow(d))
}

add_indices_for <- function(d, suffix = c("mean","median")) {
  suffix <- match.arg(suffix)
  
  # map sentinel bands to reflectance shortcuts for this suffix
  R550 <- getv(d, paste0("B3_",  suffix))   # ~550 nm (Green)
  R670 <- getv(d, paste0("B4_",  suffix))   # ~665 nm (Red)
  R700 <- getv(d, paste0("B5_",  suffix))   # ~705 nm (Red-edge)
  R740 <- getv(d, paste0("B6_",  suffix))   # ~740 nm (Red-edge)
  R780 <- getv(d, paste0("B7_",  suffix))   # ~783 nm (Red-edge)
  R800 <- getv(d, paste0("B8_",  suffix))   # ~833 nm (NIR broad)
  # B8A available if you want NDVI_8A later: getv(d, paste0("B8A_", suffix))
  
  # NDVI
  ndvi <- safe_div(R800 - R670, R800 + R670)
  
  # CARI (broadband form; 700≈B5)
  a <- safe_div(R700 - R550, 150)
  b <- R550 - a * 550
  cari <- safe_div(R700, R670) * abs(a * 670 + b - R670) * sqrt(a^2 + 1)
  
  # TCARI2/OSAVI2 (705/750 variant using 740 as proxy for 750)
  tcari2  <- 3 * ((R740 - R700) - 0.2 * (R740 - R550)) * safe_div(R740, R700)
  osavi2  <- (1 + 0.16) * safe_div((R740 - R700), (R740 + R700 + 0.16))
  tcari_osavi <- safe_div(tcari2, osavi2)
  
  # REP (Guyot–Baret style; 700≈B5, 740≈B6, 780≈B7)
  Rre <- (R670 + R780) / 2
  rep_val <- 700 + 40 * safe_div((Rre - R700), (R740 - R700))
  
  # Boochs (approx 703): use R705 ≈ B5
  boochs <- R700
  
  # Datt3 proxy (derivative ratio not possible): use red-edge ratio 740/705 = B6/B5
  datt3 <- safe_div(R740, R700)
  
  # PRI proxy (green only, per your instruction)
  pri_green <- R550
  
  # Vogelmann3 proxy (derivative ratio not possible): use 740/705 = B6/B5
  vog3_proxy <- safe_div(R740, R700)
  
  # attach with appropriate names
  d[[paste0("NDVI_",            suffix)]] <- ndvi
  d[[paste0("CARI_",            suffix)]] <- cari
  d[[paste0("TCARI_OSAVI_",     suffix)]] <- tcari_osavi
  d[[paste0("REP_",             suffix)]] <- rep_val
  d[[paste0("Boochs_",          suffix)]] <- boochs
  d[[paste0("Datt3_",           suffix)]] <- datt3
  d[[paste0("PRI_green_",       suffix)]] <- pri_green
  d[[paste0("VOG3_proxy_",      suffix)]] <- vog3_proxy
  
  d
}

# apply for both mean and median, keeping all existing columns
df <- add_indices_for(df, "mean")
df <- add_indices_for(df, "median")

# final output
dendro_spectra_s2_indices <- df

# Correlation matrix - Sentinel bands/indices (rows) vs dendro metrics (columns)

# input dataframe
df <- dendro_spectra_s2_indices

# define dendro metric columns
dendro_cols <- c("cTWD", "cTWD_drone", "cshrink", "TWD_drone",
                 names(df)[2778:2804])

# define Sentinel band/index columns
spectral_cols <- names(df)[2886:2921]

# extract sub-dataframes
dendro_mat   <- df[, dendro_cols, drop = FALSE]
spectral_mat <- df[, spectral_cols, drop = FALSE]

# compute Spearman correlation matrix (rows = spectral, cols = dendro)
sentinel_mat <- cor(spectral_mat, dendro_mat,
                    method = "spearman", use = "pairwise.complete.obs")

# view result
dim(sentinel_mat)         # [# spectral vars x # dendro metrics]
head(sentinel_mat[, 1:5]) # peek at first 5 dendro metrics

# final output object
spearman_matrix <- sentinel_mat

# Export
# keep only mean-based Sentinel rows (drop medians)
mean_rows <- !grepl("_median$", rownames(spearman_matrix), ignore.case = TRUE)
spearman_mean_only <- spearman_matrix[mean_rows, , drop = FALSE]

# export to CSV
out_dir  <- "G:/HyperspectralUAV/R_outputs/modelling/correlations"
out_file <- file.path(out_dir, "spearman_sentinel_means_vs_dendro.csv")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(spearman_mean_only, out_file, row.names = TRUE)

#-------CREATE CORRELATION MATRIX BETWEEN SEVERAL VARIABLES (OLD)----------------------
# Ensure the output directory exists
dir.create("./R_outputs/modelling/correlations/", recursive = TRUE, showWarnings = FALSE)

# Define the column indices clearly
spectral_cols <- 2:1602       # Columns of spectral variables
dendro_cols <- 1962:1967          # Columns of zg_fraction variants

# Extract variable names
spectral_vars <- names(dendro_spectra)[spectral_cols]
zg_vars <- names(dendro_spectra)[zg_cols]

# Initialize an empty correlation dataframe (1960 rows, 6 columns)
cor_matrix <- data.frame(variable = spectral_vars)

# Calculate Spearman correlations and store results
for (zg_var in zg_vars) {
  cor_matrix[[zg_var]] <- sapply(spectral_vars, function(spec_var) {
    cor(dendro_spectra[[zg_var]], dendro_spectra[[spec_var]], method = "spearman", use = "complete.obs")
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
################# REMOVE CORRELATED PREDICTORS USING *magic* EMPHASIZING HIGHLY CORRELATED VARIABLES TO ZG_FRACTION
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
