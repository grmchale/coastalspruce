#setwd("paintrock_leaf_spectra")
source("Functions/writeSLI_source.R")
source("Functions/lecospectR.R")
require(Polychrome)
require(vegan)
require(glue)

paintrock_spectra_df<-read.csv("./R_outputs/paintrock_spectra_clean.csv")

#Summarize reflectance by speccies and rearrange columns such that the first 
#is wavelength and every subsequent column is a sample, in this case median reflectance by species.

#spruce_spectra_df<-MASTER_Site_speclib %>% 
#    group_by(Site) %>%
#    dplyr::select(`X398`:`X999`) %>% #dim
#    pivot_longer(cols = `X398`:`X999`,  names_to  = "Wavelength", values_to = "Reflectance") %>% #dim
#    mutate(Wavelength = gsub("X","",Wavelength)) %>%
#    group_by(taxon_code,Wavelength) %>%  
#    dplyr::summarise(Reflectance = median(Reflectance))%>%
#    mutate(Wavelength = as.numeric(Wavelength)) %>%
#    as.data.frame() %>% 
#    pivot_wider(names_from = taxon_code, values_from = Reflectance) %>%
#    mutate(Wavelength = as.numeric(Wavelength)) %>%
#    dplyr::arrange(Wavelength)

#writeSLI(spruce_spectra_df,"output/paintrock_spectra.sli", wavl.units = "Nanometers")

library(dplyr)
library(tidyr)
library(ggplot2)

#SPECTRAL PROFILES USING AVERAGE REFLECTANCE VALUES

# Assuming MASTER_Site_speclib has
# - A "Site" column with 3 sites.
# - Columns named X398 through X999 for reflectance values.

site_spectra <- MASTER_Site_speclib %>%
  dplyr::select(-X) %>%
  pivot_longer(cols = starts_with("X"),    # All columns that start with "X"
               names_to = "Wavelength", 
               values_to = "Reflectance") %>%
  mutate(Wavelength = as.numeric(gsub("X", "", Wavelength))) %>%
  arrange(Wavelength)

# Now plot the spectral profiles for each site in one figure
ggplot(site_spectra, aes(x = Wavelength, y = Reflectance, color = Site)) +
  geom_line(size=1.5) +
  labs(title = "Spectral Profiles by Site",
       x = "Wavelength (nm)",
       y = "Reflectance") +
  theme_minimal()

#SPECTRAL PROFILES USING MEDIAN REFLECTANCE VALUES

# Prepare the spectral data
site_spectra <- MASTER_TreeID_speclib_2 %>%
  mutate(Site = substr(TreeID, 1, 2)) %>%               # Extract Site code
  pivot_longer(
    cols = starts_with("X"),                            # Select reflectance columns
    names_to = "Wavelength", 
    values_to = "Reflectance"
  ) %>%
  mutate(
    Wavelength = as.numeric(gsub("X", "", Wavelength)), # Convert Wavelength to numeric
    Reflectance = as.numeric(Reflectance)              # Ensure Reflectance is numeric
  ) %>%
  group_by(Site, Wavelength) %>%
  summarise(
    Reflectance = median(Reflectance, na.rm = TRUE)    # Aggregate by median
  ) %>%
  ungroup() %>%
  arrange(Wavelength)

# Verify Reflectance range
print(range(site_spectra$Reflectance))  # Should be between 0 and 1

# Define custom colors for sites
custom_colors <- c("CE" = "red", "ET" = "green", "RI" = "blue")

# Create the spectral profile plot
spectral_plot <- ggplot(site_spectra, aes(x = Wavelength, y = Reflectance, color = Site)) +
  geom_line(size = 1.5) +                                      # Thicker lines
  scale_color_manual(values = custom_colors) +                 # Assign custom colors
  labs(
    title = "Spectral Profiles by Site - Median Values",
    x = "Wavelength (nm)",
    y = "Reflectance"
  ) +
  theme_minimal() +                                            # Minimal theme
  ylim(0, 0.3)                                                   # Fix y-axis range

# Display the plot
print(spectral_plot)

# Save the plot to a PNG file (optional)
ggsave("Spectral_Profiles_by_Site.png", plot = spectral_plot, width = 10, height = 6, dpi = 300)

###ALL OF THE TREE SPECTRAL PROFILES: EACH LINE IS A DIFFERENT TREE###

# Prepare the spectral data
spectral_data <- Canopy_image_spectra_bytreeid %>%
  # Extract the first two characters of TreeID to get the Site code
  mutate(Site = substr(TreeID, 1, 2)) %>%
  
  # Transform the data from wide to long format using the wavelength columns (398 to 999)
  pivot_longer(
    cols = all_of(as.character(seq(398, 998, by = 15))),
    #cols = all_of(as.character(398:999)),   # Now the columns are named "398" to "999"
    names_to = "Wavelength",                # New column for wavelength names
    values_to = "Reflectance"               # New column for reflectance values
  ) %>%
  
  # Clean and convert the data
  mutate(
    Wavelength = as.numeric(Wavelength),    # Convert wavelength strings to numeric
    Reflectance = as.numeric(Reflectance)   # Ensure Reflectance is numeric
  ) %>%
  
  # Arrange the data by TreeID and Wavelength for plotting
  arrange(TreeID, Wavelength)

# Verify the range of Reflectance
print(range(spectral_data$Reflectance))  # Should output values between 0 and 1

# Define custom colors for the sites
custom_colors <- c(
  "CE" = "red", 
  "CC" = "orange", 
  "GI" = "blue", 
  "FP" = "green", 
  "RI" = "purple", 
  "HI" = "brown"
)

# Create the spectral profile plot
spectral_plot_vr <- ggplot(spectral_data, aes(x = Wavelength, y = Reflectance, group = TreeID, color = Site)) +
  geom_line(size = 1) +                                     # Thicker lines; adjust 'size' as needed
  scale_color_manual(values = custom_colors) +              # Assign custom colors
  labs(
    title = "Reflectance of Individual Trees by Site, 15 nm Smoothed - Visble Range",
    x = "Wavelength (nm)",
    y = "Reflectance"
  ) +
  theme_minimal() +                                         # Minimal theme for clarity
  ylim(.15, .35) +                                             # Fix y-axis range (adjust as needed)
  coord_cartesian(xlim = c(750, 850))                         # Limit x-axis to wavelengths 398-700

# Display the plot
print(spectral_plot_vr)

# Save the plot to a PNG file (optional)
ggsave("Spectral_Profiles_by_Site.png", plot = spectral_plot, width = 10, height = 6, dpi = 300)



