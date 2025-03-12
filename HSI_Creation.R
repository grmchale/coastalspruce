#Read in canopy VIs
canopy_VIs<-read.csv("./R_outputs/speclib_dendrometers/veg_indices/dendrometer_VIs_5nm_bytreeid.csv")

# Scale the indices and combine into a singular stress index (HSI)
indices <- c("ARI1", "ARI2", "CRI1", "CRI2", "PSRI")

# Scaling indices to a standard range (0-1)
canopy_VIs_scaled <- canopy_VIs
canopy_VIs_scaled[indices] <- lapply(canopy_VIs[indices], function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
})

# Compute HSI as equally-weighted average
canopy_VIs_scaled$HSI <- rowMeans(canopy_VIs_scaled[indices], na.rm = TRUE)

# Check resulting dataframe
head(canopy_VIs_scaled)
