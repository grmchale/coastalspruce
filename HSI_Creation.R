library(dplyr)
library(ggplot2)
###Read in canopy VIs and dendro attributes###
canopy_VIs<-read.csv("./R_outputs/speclib_dendrometers/veg_indices/dendrometer_VIs_15nm_bytreeid.csv")

# Scale the indices and combine into a singular stress index (HSI)
indices <- c("ARI1", "ARI2", "CRI1", "CRI2", "PSRI", "PRI")

# Scaling indices to a standard range (0-1)
canopy_VIs_scaled <- canopy_VIs
canopy_VIs_scaled[indices] <- lapply(canopy_VIs[indices], function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
})

# Compute HSI as equally-weighted average
canopy_VIs_scaled$HSI <- rowMeans(canopy_VIs_scaled[indices], na.rm = TRUE)

# Check resulting dataframe
head(canopy_VIs_scaled)

#-----------BOX PLOT FOR NEW CUSTOM HSI---------------------------

# BOX PLOT BY TREEID!!
# Aggregate the pixel-level data by TreeID using the median value for each index
agg_data <- canopy_VIs_scaled %>%
  group_by(TreeID) %>%
  summarise(across(c(HSI), median, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Site = substr(TreeID, 1, 2))  # Extract site code from TreeID

# Pivot the aggregated data to long format for plotting
long_agg_data <- agg_data %>%
  pivot_longer(cols = c(HSI),
               names_to = "Index",
               values_to = "Value")

# Create the boxplot for each index across sites using the aggregated data
ggplot(long_agg_data, aes(x = Site, y = Value, fill = Site)) +
  geom_boxplot() +
  facet_wrap(~Index, scales = "free_y") +
  labs(title = "Combined Tree Canopy Stress Indices (15 nm Resample) - Tree Median Value",
       x = "Site",
       y = "Index Value") +
  theme_minimal()

#-------------SCATTERPLOTS OF VIS VS DENDRO ATTRIBUTES-------------------------

###Join dendro attributes to VIs###
# Perform the joins using 'Dndrmtr' as the key
canopy_VIs_joined <- canopy_VIs_scaled %>%
  left_join(cumulative_zg %>% select(Dndrmtr, zg_days), by = "Dndrmtr") %>%
  left_join(drone_dendro %>% select(Dndrmtr, TWD), by = "Dndrmtr")

###Quick scatter plot and simple models for VI vs dendro attributes###
VI <- "CRI1"
# Scatterplot: VI vs zg_days
ggplot(canopy_VIs_joined, aes(x = zg_days, y = .data[[VI]])) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  geom_smooth(method = "lm", formula = y ~ log(x), color = "red") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "green") +
  labs(title = paste(VI, "vs. zg_days (Linear, Log, Polynomial)"), x = "zg_days", y = VI)

# Linear model for VI vs zg_days
lm_zg_linear <- lm(reformulate("zg_days", response = VI), data = canopy_VIs_joined)
summary(lm_zg_linear)

# Logarithmic model for VI vs zg_days
lm_zg_log <- lm(as.formula(paste(VI, "~ log(zg_days)")), data = canopy_VIs_joined)
summary(lm_zg_log)

# Polynomial (quadratic) model for VI vs zg_days
lm_zg_poly <- lm(as.formula(paste(VI, "~ poly(zg_days, 2)")), data = canopy_VIs_joined)
summary(lm_zg_poly)

# Scatterplot: VI vs TWD
ggplot(canopy_VIs_joined, aes(x = TWD, y = .data[[VI]])) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_smooth(method = "lm", formula = y ~ log(x), color = "blue") +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), color = "green") +
  labs(title = paste(VI, "vs. TWD (Linear, Log, Polynomial)"), x = "TWD", y = VI)

# Linear model for VI vs TWD
lm_twd_linear <- lm(reformulate("TWD", response = VI), data = canopy_VIs_joined)
summary(lm_twd_linear)

# Logarithmic model for VI vs TWD
lm_twd_log <- lm(as.formula(paste(VI, "~ log(TWD)")), data = canopy_VIs_joined)
summary(lm_twd_log)

# Polynomial (quadratic) model for VI vs TWD
lm_twd_poly <- lm(as.formula(paste(VI, "~ poly(TWD, 2)")), data = canopy_VIs_joined)
summary(lm_twd_poly)

