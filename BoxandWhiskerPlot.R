#Code for creating box and whisker plots with spectral indicies across sites

#Read in df or .csv
MASTER_TreeID_speclib

library(dplyr)
library(tidyr)
library(ggplot2)
library(dplyr)

#Adding ARI1 and ARI2 to the df
MASTER_TreeID_speclib_clean <- MASTER_TreeID_speclib_clean %>%
  mutate(
    ARI1_raw = (1 / X550) - (1 / X700),
    ARI2_raw = X800 * ((1 / X550) - (1 / X700))
  )

# Compute the min and max for ARI1_raw and ARI2_raw
ARI1_min <- min(MASTER_TreeID_speclib_clean$ARI1_raw, na.rm = TRUE)
ARI1_max <- max(MASTER_TreeID_speclib_clean$ARI1_raw, na.rm = TRUE)

ARI2_min <- min(MASTER_TreeID_speclib_clean$ARI2_raw, na.rm = TRUE)
ARI2_max <- max(MASTER_TreeID_speclib_clean$ARI2_raw, na.rm = TRUE)

# Apply min-max scaling to get ARI1 and ARI2 into [0,1] range
MASTER_TreeID_speclib_clean <- MASTER_TreeID_speclib_clean %>%
  mutate(
    ARI1 = (ARI1_raw - ARI1_min) / (ARI1_max - ARI1_min),
    ARI2 = (ARI2_raw - ARI2_min) / (ARI2_max - ARI2_min)
  )


# Clean the data to remove outliers
MASTER_TreeID_speclib_clean <- MASTER_TreeID_speclib %>%
  filter(CRI1 > -1e10)  # Adjust threshold as needed

# Create the boxplot - STRESS INDICES

# Extract site code and pivot indices to long format
long_data <- MASTER_TreeID_speclib_clean %>%
  mutate(Site = substr(TreeID, 1, 2)) %>% 
  # Pivot CRI1, CRI2, PSRI into long format
  pivot_longer(cols = c(CRI1, CRI2, PSRI, ARI1, ARI2), 
               names_to = "Index", 
               values_to = "Value")

# Now 'long_data' has columns: TreeID, Site, Index (CRI1/CRI2/PSRI), and Value

# This createS 3 separate panels, one for each index, and inside each panel
# you'll see boxplots for CE, ET, and RI.

ggplot(long_data, aes(x = Site, y = Value, fill = Site)) +
  geom_boxplot() +
  facet_wrap(~Index, scales = "free_y") + 
  # 'scales = "free_y"' is optional. If you think all indices have similar scales,
  # you can leave it as scales = "fixed".
  labs(title = "Comparison of Stress Indicies (CRI1, CRI2, PSRI, ARI1, ARI2) Across Sites",
       x = "Site",
       y = "Index Value") +
  theme_minimal()

# If you'd prefer to see all 9 boxes in a single panel (3 sites × 3 indices),
# you can rearrange the aes mapping. For example, put Index on the x-axis and color by Site:

 #ggplot(long_data, aes(x = Index, y = Value, fill = Site)) +
  # geom_boxplot(position = position_dodge(width = 0.8)) +
  # labs(title = "Comparison of CRI1, CRI2, and PSRI across Sites",
  #      x = "Index",
  #      y = "Index Value") +
  # theme_minimal()

# In this alternative layout, you'll get 3 groups on the x-axis (CRI1, CRI2, PSRI),
# and each group will have 3 boxplots (one for each site) side-by-side.

#Create the boxplot - NDVI INDICES
# Extract site code and pivot indices to long format
ndvi_data <- MASTER_TreeID_speclib_clean %>%
  mutate(Site = substr(TreeID, 1, 2)) %>% 
  # Pivot CRI1, CRI2, PSRI into long format
  pivot_longer(cols = c(NDVI, mNDVI, GreenNDVI), 
               names_to = "Index", 
               values_to = "Value")

# Now 'ndvi_data' has columns: TreeID, Site, Index (NDVI/mNDVI/GreenNDVI), and Value

# This createS 3 separate panels, one for each index, and inside each panel
# you'll see boxplots for CE, ET, and RI.

ggplot(ndvi_data, aes(x = Site, y = Value, fill = Site)) +
  geom_boxplot() +
  facet_wrap(~Index, scales = "free_y") + 
  # 'scales = "free_y"' is optional. If you think all indices have similar scales,
  # you can leave it as scales = "fixed".
  labs(title = "Comparison of NDVI Indicies (NDVI, mNDVI, and GreenNDVI) Across Sites",
       x = "Site",
       y = "Index Value") +
  theme_minimal()
