# Load necessary packages
library(dplyr)
library(vegan)
library(ggplot2)

####################### ADONIS perMANOVA ############################
# Randomly sample 200 rows per site
set.seed(123) # for reproducibility
VIs_df_sampled <- VIs_df %>%
  group_by(Sample) %>%
  slice_sample(n = 200) %>%
  ungroup()

# Extract the vegetation indices you want to test
indices_matrix <- VIs_df_sampled %>%
  select(ARI1, ARI2, CRI1, CRI2, PSRI, PRI)

# Run perMANOVA (adonis2)
adonis_results <- adonis2(indices_matrix ~ Sample, 
                          data = VIs_df_sampled, 
                          method = "euclidean", 
                          permutations = 500)

# Check results
print(adonis_results)

####################### PAIRWISE PERMANOVA ############################
library(pairwiseAdonis)

# Perform pairwise perMANOVA comparisons
pairwise_results <- pairwise.adonis(indices_matrix, 
                                    factors = VIs_df_sampled$Sample, 
                                    sim.method = "euclidean", 
                                    p.adjust.m = "fdr") # adjust for multiple comparisons

# View pairwise comparison results
pairwise_results
#Write to .csv
write.csv(adonis_results, "./R_outputs/modelling/perMANOVA/adonispma_results_amoebas_v1.csv")
write.csv(pairwise_results,"./R_outputs/modelling/perMANOVA/pairwise_results_amoebas_v1.csv")

################# Non-metric Multidimensional Scaling (NMDS) ################
# NMDS using Euclidean distance
set.seed(42) # for reproducibility
nmds <- metaMDS(VIs_df_sampled[, c("ARI1", "ARI2", "CRI1", "CRI2", "PSRI", "PRI")],
                distance = "euclidean", 
                k = 2, 
                trymax = 100)

# Extract site scores
nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$Site <- VIs_df_sampled$Sample

# Plot using ggplot2
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Site)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  theme(legend.position = "right") +
  ggtitle("NMDS: Vegetation Index Differences Between Sites")
# Plot using ggplot2 with axis limits to remove outliers
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Site)) +
  geom_point(alpha = 0.5) +
  coord_cartesian(xlim = c(-50, 100), ylim = c(-100, 100)) +  # zooms in without dropping data
  theme_minimal() +
  theme(legend.position = "right") +
  ggtitle("NMDS: Vegetation Index Differences Between Sites")

############### BOX PLOT CREATION ##################################

###Read in canopy VIs and dendro attributes###
canopy_VIs<-VIs_df
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

##### Box plot ####
library(tidyr)

# Define custom site order
site_order <- c("GI", "CE", "HI", "BI_North", "BI_South", "CC_North", "CC_South",
                "ET", "FP", "WP", "GW", "RI_IP", "RI_GB")

# Random sample of 1000 pixels per site
set.seed(124)
sampled_data <- canopy_VIs_scaled %>%
  group_by(Sample) %>%
  slice_sample(n = 1000) %>%
  ungroup() %>%
  mutate(Sample = factor(Sample, levels = site_order))  # apply site order once

### 1. PRI Boxplot
ggplot(sampled_data, aes(x = Sample, y = PRI, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4) +
  coord_cartesian(ylim = c(0.45, 0.6)) +
  labs(title = "PRI Distributions Across Coastal Spruce Sites",
       x = "Site",
       y = "PRI Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

### 2. PSRI Boxplot
ggplot(sampled_data, aes(x = Sample, y = PSRI, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4) +
  coord_cartesian(ylim = c(0.275, .425)) +
  labs(title = "PSRI Distributions Across Coastal Spruce Sites",
       x = "Site",
       y = "PSRI Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

### 3. CRI1 Boxplot
ggplot(sampled_data, aes(x = Sample, y = CRI1, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4) +
  coord_cartesian(ylim = c(0.085, .12)) +
  labs(title = "CRI1 Distributions Across Coastal Spruce Sites",
       x = "Site",
       y = "CRI1 Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

install.packages("patchwork")
library(patchwork)

# Create shared theme elements
shared_theme <- theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.text = element_blank())  # removes facet-style titles

# PRI
p1 <- ggplot(sampled_data, aes(x = Sample, y = PRI, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4) +
  coord_cartesian(ylim = c(0.5, 0.6)) +
  ylab("PRI Value") +
  xlab("") +
  shared_theme

# PSRI
p2 <- ggplot(sampled_data, aes(x = Sample, y = PSRI, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4) +
  coord_cartesian(ylim = c(0.2, 0.4)) +
  ylab("PSRI Value") +
  xlab("") +
  shared_theme

# CRI1
p3 <- ggplot(sampled_data, aes(x = Sample, y = CRI1, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4) +
  coord_cartesian(ylim = c(0, .2)) +
  ylab("CRI1 Value") +
  xlab("Site") +  # only this one gets the x-axis label
  shared_theme

# Combine vertically
p1 / p2 / p3

