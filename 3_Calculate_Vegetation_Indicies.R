source("Functions/lecospectR.R")
#source("./Scripts/Utilities/install_dependencies.R")

library(spectrolab)
library(terra)
library(raster)
library(rgdal)
library(hsdar)
library(rgeos)
library(sf)
library(doFuture)
source("Functions/spectral_operations.R")

###CALCULATE VEGETATION INDICIES FOR LOADED IN SPECTRAL LIBRARY###
#RDS file for using when returning to code
plots_image_spectra<-readRDS("./R_outputs/speclib_dendrometers/dataframes/dendrometer_canopy_speclib_15nm.rds")

#Cast the spectral library to a data frame, select which version of resampled spectral library (5, 10, 15 nm)
plots_image_spectra_df<- speclib_to_df(plots_image_spectra)
plots_image_spectra_df<- plots_image_spectra_df %>% select(-c(3:10)) #remove other attributes besides wavelengths

#Check unique entries
plots_image_spectra_df %>% group_by(TreeID) %>% tally() %>% 
  filter(n>1) #%>% ungroup() %>% select(Flight, Field, Plot) %>% unique %>% print(n=200)
getAnywhere("get_required_veg_indices")

#Calculate vegetation indices for the pixels
plots_image_spectra_VIs<-get_vegetation_indices(plots_image_spectra_df, NULL)
#plots_image_spectra_VIs<-calc_veg_index(plots_image_spectra_df, NA)

#Adds categorical variables back in
plots_image_spectra_VIs<-cbind(as.data.frame(plots_image_spectra)[,1:10],plots_image_spectra_VIs) 

##Create median VI value per TreeID
median_VIs <- plots_image_spectra_VIs %>%
  #select(-106) %>%
  group_by(TreeID) %>%
  summarise(across(10:104, ~ median(.x, na.rm = TRUE)))
# Extract a unique row per TreeID from the original dataframe
plot_info <- plots_image_spectra_VIs %>%
  select(TreeID, sample_name, Site, Species, X, Y, CC, Dndrmtr, DBH, Tag) %>%
  distinct(TreeID, .keep_all = TRUE)
# Join the median_VIs (which has one row per TreeID) with the unique plot_info
median_VIs_final <- median_VIs %>%
  left_join(plot_info, by = "TreeID") %>%
  select(TreeID, sample_name, Site, Species, X, Y, CC, Dndrmtr, DBH, Tag, everything())
# Read in the other dataframe and calculate ARI1 and ARI2 using values from columns "550" and "700"
#speclib_bytreeid <- readRDS("./R_outputs/speclib_dendrometers/dataframes/dendrometer_canopy_bytreeid.rds") %>%
  #mutate(
    #ARI1 = (1 / `550`) - (1 / `700`),
    #ARI2 = `800` * ((1 / `550`) - (1 / `700`))
 # )
#Use this one for 5nm, 10m, and 15nm:
speclib_bytreeid <- readRDS("./R_outputs/speclib_dendrometers/dataframes/dendrometer_canopy_bytreeid_15nm.rds") %>%
  mutate(
    ARI1 = (1 / `548`) - (1 / `698`),
    ARI2 = `803` * ((1 / `548`) - (1 / `698`)) #Use '798' for 5nm and 10nm
  )
# Join ARI1 and ARI2 from speclib_bytreeid into median_VIs_final (keyed by TreeID)
median_VIs_final <- median_VIs_final %>%
  left_join(speclib_bytreeid %>% select(TreeID, ARI1, ARI2), by = "TreeID") %>%
  # Move ARI1 and ARI2 to immediately after column 10
  select(1:10, ARI1, ARI2, everything())

#Export VIs to .csv
write.csv(plots_image_spectra_VIs,"./R_outputs/speclib_dendrometers/veg_indices/dendrometer_VIs_15nm.csv")
write.csv(median_VIs_final,"./R_outputs/speclib_dendrometers/veg_indices/dendrometer_VIs_15nm_bytreeid.csv")

###CREATE SEPERATE .CSVs PARSING OUT SEPERATE PORTIONS OF THE STATS###
library(dplyr)
library(tidyr)

pivot_wide <- plots_image_spectra_VIs_stats_df_final %>%
  # Select only the columns we want to keep
  select(Site, Species, X, Y, CC, UID, VegIndices, Mean) %>%
  
  # Pivot the data so that each VegIndices becomes a column
  pivot_wider(
    names_from = VegIndices,
    values_from = Mean
  )
#Join DBH, Dendrometer, Notes
pivot_wide_2 <- pivot_wide %>%
  left_join(
    Overstory %>% select(TreeID, DBH, Dendrometer, Notes),
    by = c("UID" = "TreeID")
  )
#Add in SD just in case
pivot_wide_SD <- plots_image_spectra_VIs_stats_df_final %>%
  # Keep the key identifying columns plus VegIndices, Mean, and SD
  select(Site, Species, X, Y, CC, UID, VegIndices, Mean, SD) %>% 
  pivot_wider(
    names_from = VegIndices, 
    values_from = c(Mean, SD), 
    names_sep = ""  # This will create columns like Boochs2Mean, Boochs2SD, etc.
  )
#Join DBH, Dendrometer, Notes to SD
pivot_wide_SD <- pivot_wide_SD %>%
  left_join(
    Overstory %>% select(TreeID, DBH, Dendrometer, Notes),
    by = c("UID" = "TreeID")
  )

#Write both of these puppies to a .csv
write.csv(pivot_wide_2, "./R_outputs/speclib/VIs_stats_1nm_means.csv")
write.csv(pivot_wide_SD, "./R_outputs/speclib/VIs_stats_1nm_means_SD.csv")

####Bin MASTER csv by TreeID/Site (only bands 398-999 here)####
library(dplyr)

MASTER_joined<-readRDS("G:/HyperspectralUAV/R_outputs/speclib/MASTER_joined.rds")

MASTER_Site_speclib <- MASTER_joined %>%
  group_by(Site) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

MASTER_TreeID_speclib_2<-MASTER_TreeID_speclib[, c(1,5:606)]
write.csv(MASTER_TreeID_speclib_2, "./R_outputs/speclib/TreeID_speclib_summary.csv" )
