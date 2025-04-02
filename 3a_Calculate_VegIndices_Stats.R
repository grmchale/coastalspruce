library(dplyr)

###CALCULATE SD, MEAN, AND QUANTILES FOR VEGETATION INDICES###

# READ IN VEGETATION INDICES----------------------------------------------------
#Resampled vegetation indices nm spacing (1,5,10,15)
NM <- "5nm"
#Read in vegetation indices by pixel
VIs_df <- read.csv(paste0("./R_outputs/speclib_dendrometers/veg_indices/dendrometer_VIs_", NM, ".csv"))

# ----------------ADD IN ARI1, ARI2, and WBI FIRST------------------------------
# Read in spectral library
speclib <- read.csv(paste0("./R_outputs/speclib_dendrometers/dendrometer_canopy_speclib_", NM, ".csv"))

# Create ARI1, ARI2, and WBI
if (NM == "1nm") {
  speclib <- speclib %>%
    mutate(
      ARI1 = (1 / `X550`) - (1 / `X700`),
      ARI2 = `X800` * ((1 / `X550`) - (1 / `X700`)),
      WBI  = (`X970` / `X900`)
    )
} else if (NM == "15nm") {
  speclib <- speclib %>%
    mutate(
      ARI1 = (1 / `X548`) - (1 / `X698`),
      ARI2 = `X803` * ((1 / `X548`) - (1 / `X698`)),  # Explicitly use X803 for 15nm
      WBI  = (`X968` / `X893`)                         # Use X893 for 15nm
    )
} else {
  speclib <- speclib %>%
    mutate(
      ARI1 = (1 / `X548`) - (1 / `X698`),
      ARI2 = `X798` * ((1 / `X548`) - (1 / `X698`)),
      WBI  = (`X968` / `X898`)
    )
}

# Directly pass calculated columns into VIs_df
VIs_df <- VIs_df %>%
  mutate(
    ARI1 = speclib$ARI1,
    ARI2 = speclib$ARI2,
    WBI  = speclib$WBI
  ) %>%
  relocate(ARI1, ARI2, .after = 11) %>%
  relocate(WBI, .after = last_col())

#----------

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
