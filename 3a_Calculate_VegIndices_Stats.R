library(dplyr)

###CALCULATE SD, MEAN, AND QUANTILES FOR VEGETATION INDICES###

# READ IN VEGETATION INDICES-------------------------------------------------------
#Resampled vegetation indices nm spacing (1,5,10,15)
NM <- "5nm"
#Read in vegetation indices by pixel
VIs_df <- read.csv(paste0("./R_outputs/speclib_dendrometers/veg_indices/dendrometer_VIs_", NM, ".csv"))

# ADD IN ARI1 AND ARI2 FIRST------------------------------------------------------
# Read in spectral library
speclib <- read.csv(paste0("./R_outputs/speclib_dendrometers/dendrometer_canopy_speclib_", NM, ".csv"))
# Calculate ARI1 and ARI2 using values from columns - 1 NM
speclib <- speclib %>%
  mutate(
    ARI1 = (1 / `X550`) - (1 / `X700`),
    ARI2 = `X800` * ((1 / `X550`) - (1 / `X700`))
  )
# Calculate ARI1 and ARI2 using values from columns - 5, 10, 15 NM
speclib <- speclib %>%
  mutate(
    ARI1 = (1 / `X548`) - (1 / `X698`),
    ARI2 = `X798` * ((1 / `X548`) - (1 / `X698`)) #Use 'X803 for 15nm
  )
# Directly pass calculated columns into VIs_df after column 11
VIs_df <- VIs_df %>%
  mutate(ARI1 = speclib$ARI1, ARI2 = speclib$ARI2) %>%  
  relocate(ARI1, ARI2, .after = 11)

# CALCULATE SD, MEAN, AND QUANTILES FOR VEGETATION INDICES W/ ARIs ------------------
# Assuming columns 12:end are vegetation indices
VI_cols <- names(VIs_df)[12:ncol(VIs_df)]

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

#Merge 1, 5, 10, 15 nm spectral libraries -------------------------------------------------------
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
