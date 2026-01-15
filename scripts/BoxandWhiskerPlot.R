#Code for creating box and whisker plots with spectral indicies across sites

library(dplyr)
library(tidyr)
library(ggplot2)
library(dplyr)

# Read in VI indices by pixels
VIs <- read.csv("./R_outputs/speclib_dendrometers/veg_indices/dendrometer_VIs_1nm.csv")
VIs_clean <- read.csv("./R_outputs/speclib_dendrometers/veg_indices/dendrometer_VIs_15nm_bytreeid.csv")
# Read in spectral library
speclib <- read.csv(paste0("./R_outputs/speclib_dendrometers/dendrometer_canopy_speclib_1nm.csv"))

# ------------------------------------------------------------------------------
# ADD IN ARI1 AND ARI2!

# Clean the data to remove outliers
#VIs_clean <- VIs %>%
 # filter(CRI1 > -1e10)  # Adjust threshold as needed

speclib <- speclib %>%
  mutate(
    ARI1 = (1 / `X550`) - (1 / `X700`),
    ARI2 = `X800` * ((1 / `X550`) - (1 / `X700`)) #Use 'X803 for 15nm
  )
# Directly pass calculated columns into VIs_df after column 11
VIs_clean <- VIs %>%
  mutate(ARI1 = speclib$ARI1, ARI2 = speclib$ARI2) %>%  
  relocate(ARI1, ARI2, .after = 11)

# ------------------------------------------------------------------------------
# BOX PLOT BY PIXEL!!
# Extract site code from TreeID and pivot indices to long format
long_data <- VIs_clean %>%
  mutate(Site = substr(TreeID, 1, 2)) %>% 
  pivot_longer(cols = c(CRI1, CRI2, PSRI, ARI1, ARI2, NDVI), 
               names_to = "Index", 
               values_to = "Value")

# Create the boxplot for each index across sites
ggplot(long_data, aes(x = Site, y = Value, fill = Site)) +
  geom_boxplot() +
  facet_wrap(~Index, scales = "free_y") + 
  labs(title = "Dendrometer Tree Canopies Stress Indices (5 nm Resample) - All Pixels",
       x = "Site",
       y = "Index Value") +
  theme_minimal()

# ------------------------------------------------------------------------------
# BOX PLOT BY TREEID!!
# Aggregate the pixel-level data by TreeID using the median value for each index
agg_data <- VIs_clean %>%
  group_by(TreeID) %>%
  summarise(across(c(CRI1, CRI2, PSRI, ARI1, ARI2, PRI), median, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Site = substr(TreeID, 1, 2))  # Extract site code from TreeID

# Pivot the aggregated data to long format for plotting
long_agg_data <- agg_data %>%
  pivot_longer(cols = c(CRI1, CRI2, PSRI, ARI1, ARI2, PRI),
               names_to = "Index",
               values_to = "Value")

# Create the boxplot for each index across sites using the aggregated data
ggplot(long_agg_data, aes(x = Site, y = Value, fill = Site)) +
  geom_boxplot() +
  facet_wrap(~Index, scales = "free_y") +
  labs(title = "Dendrometer Tree Canopies Stress Indices (15 nm Resample) - Tree Median Value",
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
