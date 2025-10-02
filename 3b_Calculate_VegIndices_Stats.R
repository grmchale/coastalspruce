library(dplyr)

### CALCULATE SD, MEAN, AND QUANTILES FOR VEGETATION INDICES ###

# base path and target nm spacings
base_path <- "./R_outputs/speclib_chronologies/veg_indices/"
NMs <- c("1nm", "5nm", "10nm", "15nm")

for (NM in NMs) {
  # READ IN VEGETATION INDICES (per-pixel)
  VIs_df <- readRDS(paste0(base_path, "chronologies_VIs_", NM, ".rds"))
  
  # Assuming columns 3:end are vegetation indices
  VI_cols <- names(VIs_df)[3:ncol(VIs_df)]
  
  # CALCULATE SD, MEAN, AND QUANTILES FOR VEGETATION INDICES W/ ARIs
  VI_stats_df <- VIs_df %>%
    group_by(Sample) %>%
    summarize(
      across(
        all_of(VI_cols),
        list(
          Mean   = ~ mean(.x, na.rm = TRUE),
          SD     = ~ sd(.x, na.rm = TRUE),
          Q25    = ~ quantile(.x, 0.25, na.rm = TRUE),
          Median = ~ quantile(.x, 0.50, na.rm = TRUE),
          Q75    = ~ quantile(.x, 0.75, na.rm = TRUE)
        ),
        .names = paste0("{.col}_", NM, "_{.fn}")
      ),
      .groups = "drop"
    )
  
  # create objects named chronologies_VIstats_#nm in the environment
  assign(paste0("chronologies_VIstats_", NM), VI_stats_df)
}

# (optional) quick confirmation
cat("Created: ", paste0("chronologies_VIstats_", NMs, collapse = ", "), "\n")

# -----------Merge 1, 5, 10, 15 nm spectral libraries OR stats(?)---------------
library(dplyr)
library(purrr)

NM_values  <- c("1nm", "5nm", "10nm", "15nm")
df_names   <- paste0("chronologies_VIstats_", NM_values)

# grab the existing data frames from the environment
VIstats_list <- mget(df_names, inherits = TRUE)

# inner-join on Sample across all four
chronologies_VIstats_combined <- reduce(VIstats_list, ~ inner_join(.x, .y, by = "Sample"))

# quick peek
dim(chronologies_VIstats_combined)
head(names(chronologies_VIstats_combined), 20)

# Save the combined dataframe explicitly as a single RDS file
saveRDS(chronologies_VIstats_combined, "./R_outputs/speclib_chronologies/veg_indices/chronologies_VIstats_combined.rds")
