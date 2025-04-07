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
NM <- "1nm"

#RDS file for using when returning to code
speclib<-readRDS(paste0("./R_outputs/speclib_amoebas/dataframes/amoebas_speclib_", NM, ".rds"))
speclib<-Canopy_image_spectra_df

#Cast the spectral library to a data frame
speclib_df<- speclib
#speclib_df<- speclib_df %>% select(-c(3:10)) #remove other attributes besides wavelengths

#Check unique entries
speclib_df %>% group_by(Sample) %>% tally() %>% 
  filter(n>1) #%>% ungroup() %>% select(Flight, Field, Plot) %>% unique %>% print(n=200)
getAnywhere("get_required_veg_indices")

###################### OPTIONAL: FILTER OUT BRIGHTEST PIXELS #######################################

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
      ARI1 = (1 / `X548`) - (1 / `X698`),
      ARI2 = `X803` * ((1 / `X548`) - (1 / `X698`)),  # Explicitly use X803 for 15nm
      WBI  = (`X968` / `X893`)                         # Use X893 for 15nm
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
VIs_df<-cbind(as.data.frame(speclib_df)[,1:2],VIs_df) 

############ CALCULATE MEDIAN VI VALUE FOR EACH TREEID ###########################
##Create median VI value per TreeID
median_VIs <- VIs_df %>%
  #select(-106) %>%
  group_by(TreeID) %>%
  summarise(across(10:104, ~ median(.x, na.rm = TRUE)))
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
write.csv(VIs_df,"./R_outputs/speclib_amoebas/veg_indices/amoeba_VIs_1nm.csv")
write.csv(median_VIs_final,"./R_outputs/speclib_dendrometers/veg_indices/dendrometer_VIs_15nm_bytreeid.csv")
