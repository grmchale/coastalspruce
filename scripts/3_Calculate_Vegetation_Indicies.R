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
library(torch)
source("./Functions/spectral_operations.R")

#Resampled vegetation indices nm spacing (1,5,10,15)
NM <- "15nm"

#RDS file for using when returning to code
#speclib<-readRDS(paste0("./R_outputs/speclib_chronologies/dataframes/chronologies_speclib_", NM, ".rds"))
speclib_df<-Canopy_image_spectra_15nm_df

#Cast the spectral library to a data frame
#speclib_df<- speclib_df %>% select(-c(3:10)) #remove other attributes besides wavelengths

#Check unique entries
speclib_df %>% group_by(Sample) %>% tally() %>% 
  filter(n>1) #%>% ungroup() %>% select(Flight, Field, Plot) %>% unique %>% print(n=200)
getAnywhere("get_required_veg_indices")

############ CALCULATE VEGETATION INDICIES FOR LOADED IN SPECTRAL LIBRARY ##########

#Calculate vegetation indices for the pixels
VIs_df<-get_vegetation_indices(speclib_df, NULL)
#VIs_df<-calc_veg_index(speclib_df, NA)

# Create ARI1, ARI2, and WBI
if (NM == "1nm") {
  speclib_df <- speclib_df %>%
    mutate(
      ARI1 = (1 / `550`) - (1 / `700`),
      ARI2 = `800` * ((1 / `550`) - (1 / `700`)),
      WBI  = (`970` / `900`)
    )
} else if (NM == "15nm") {
  speclib_df <- speclib_df %>%
    mutate(
      ARI1 = (1 / `548`) - (1 / `698`),
      ARI2 = `803` * ((1 / `548`) - (1 / `698`)),  # Explicitly use X803 for 15nm
      WBI  = (`968` / `893`)                         # Use X893 for 15nm
    )
} else {
  speclib_df <- speclib_df %>%
    mutate(
      ARI1 = (1 / `548`) - (1 / `698`),
      ARI2 = `798` * ((1 / `548`) - (1 / `698`)),
      WBI  = (`968` / `898`)
    )
}

# Directly pass calculated columns into VIs_df
VIs_df <- VIs_df %>%
  mutate(
    ARI1 = speclib_df$ARI1,
    ARI2 = speclib_df$ARI2,
    WBI  = speclib_df$WBI
  ) %>%
  relocate(ARI1, ARI2, .before = 1) %>%
  relocate(WBI, .after = last_col())
#Adds categorical variables back in
VIs_df<-cbind(as.data.frame(speclib_df)[,1:3],VIs_df) 

############ CALCULATE MEDIAN VI VALUE FOR EACH TREEID ###########################
library(dplyr)
##Create median VI value per TreeID
median_VIs <- VIs_df %>%
  group_by(Sample) %>%
  summarise(across(ARI1:WBI, ~ median(.x, na.rm = TRUE)), .groups = "drop")
# Extract a unique row per TreeID from the original dataframe
plot_info <- VIs_df %>%
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
#speclib_bytreeid <- readRDS("./R_outputs/speclib_dendrometers/dataframes/dendrometer_canopy_bytreeid_15nm.rds") %>%
#  mutate(
#    ARI1 = (1 / `548`) - (1 / `698`),
#    ARI2 = `803` * ((1 / `548`) - (1 / `698`)) #Use '798' for 5nm and 10nm
#  )
# Join ARI1 and ARI2 from speclib_bytreeid into median_VIs_final (keyed by TreeID)
median_VIs_final <- median_VIs_final %>%
  left_join(speclib_bytreeid %>% select(TreeID, ARI1, ARI2), by = "TreeID") %>%
  # Move ARI1 and ARI2 to immediately after column 10
  select(1:10, ARI1, ARI2, everything())

# CSVs
write.csv(
  VIs_df,
  paste0("./R_outputs/speclib_amoebas_final/veg_indices/amoebas_VIs_", NM, ".csv"),
  row.names = FALSE
)
write.csv(
  median_VIs,
  paste0("./R_outputs/speclib_chronologies/veg_indices/chronologies_VIs_", NM, "_bytreeid.csv"),
  row.names = FALSE
)
# RDS
saveRDS(
  VIs_df,
  paste0("./R_outputs/speclib_amoebas_final/veg_indices/amoebas_VIs_", NM, ".rds")
)
saveRDS(
  median_VIs,
  paste0("./R_outputs/speclib_chronologies/veg_indices/chronologies_VIs_", NM, "_bytreeid.rds")
)

