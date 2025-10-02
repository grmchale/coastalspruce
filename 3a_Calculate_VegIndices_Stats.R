library(dplyr)

###CALCULATE SD, MEAN, AND QUANTILES FOR VEGETATION INDICES###

# READ IN VEGETATION INDICES----------------------------------------------------
#Resampled vegetation indices nm spacing (1,5,10,15)
NM <- "5nm"
#Read in vegetation indices by pixel
VIs_df <- readRDS(paste0("./R_outputs/speclib_chronologies/veg_indices/chronologies_VIs_", NM, ".rds"))

# ---------CALCULATE SD, MEAN, AND QUANTILES FOR VEGETATION INDICES W/ ARIs ------------------

# Assuming columns 12:end are vegetation indices
VI_cols <- names(VIs_df)[12:ncol(VIs_df)]
print(VI_cols)

# Group by TreeID and calculate mean, SD, and quantiles for each vegetation index explicitly
VI_stats_df <- VIs_df %>%
  group_by(TreeID) %>%
  summarize(across(
    all_of(VI_cols),
    list(
      Mean = ~mean(.x, na.rm = TRUE),
      SD = ~sd(.x, na.rm = TRUE),
      Q25 = ~quantile(.x, 0.25, na.rm = TRUE),
      Median = ~quantile(.x, 0.50, na.rm = TRUE),
      Q75 = ~quantile(.x, 0.75, na.rm = TRUE)
    ),
    .names = paste0("{.col}_", NM, "_{.fn}")
  )) %>%
  ungroup()

#Save to .csv
write.csv(
  VI_stats_df,
  paste0("./R_outputs/speclib_dendrometers/veg_indices_stats/dendrometer_VIstats_", NM, ".csv"),
  row.names = FALSE
)

# -----------Merge 1, 5, 10, 15 nm spectral libraries OR stats(?)---------------
library(purrr)
NM_values <- c("1nm", "5nm", "10nm", "15nm")

# Read and store each csv file in a list explicitly
VIstats_list <- lapply(NM_values, function(nm){
  read.csv(paste0("./R_outputs/speclib_dendrometers/veg_indices_stats/dendrometer_VIstats_", nm, ".csv"))
})
# Explicitly join them into one dataframe using TreeID for safety
VIstats_combined <- reduce(VIstats_list, ~inner_join(.x, .y, by = "TreeID"))

# Save the combined dataframe explicitly as a single CSV
write.csv(VIstats_combined, "./R_outputs/speclib_dendrometers/veg_indices_stats/dendrometer_VIstats_MERGED.csv", row.names = FALSE)
