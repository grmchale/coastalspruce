source("Functions/lecospectR.R")
library(spectrolab)
library(terra)
library(raster)
library(rgdal)
library(hsdar)
library(rgeos)
library(sf)
library(doFuture)

# Set directory (where the cropped hyperspectral .dat files are stored)
path <- "./R_outputs/canopy_spectra_amoebas/"

# List all .dat files of FULL canopies
allfiles <- list.files(path)
# Update pattern to match .dat files (not .ENVI)
imgs <- subset(allfiles, grepl("\\.dat$", allfiles))
print(imgs)

#### PARSE LOOP ####
Canopy_labeled <- lapply(1:length(imgs), function(x){ 
  # Bring in image
  tst <- terra::rast(paste0(path, imgs[x]))
  tst_names <- names(tst)
  
  # Create data frame with bands and values
  df <- as.data.frame(tst)
  
  # Extract and apply band names using your custom function
  new_names <- extract_bands(df)
  names(df) <- new_names
  
  # Filter bands and convert to a spectral library object
  df <- filter_bands(df)
  df <- df_to_speclib(df, type = "spectrolab")
  df <- spectrolab::resample(df, new_bands = seq(398, 999, 1), fwhm = 1)
  
  # Use the .dat file name (without the "_Hyp_PIRU" portion and extension) as metadata.
  file_name <- gsub("_Hyp_PIRU", "", imgs[x])
  file_name <- gsub("\\.dat$", "", file_name)
  
  # Correct metadata assignment (using `nrow(df)` instead of `length(df)`)
  meta(df) <- data.frame(Sample = rep(file_name, times = nrow(df)))
  
  return(df)
})


Canopy_image_spectra <- Reduce(spectrolab::combine, Canopy_labeled)
Canopy_image_spectra
dim(as.data.frame(Canopy_image_spectra))
#as.data.frame(Canopy_image_spectra) %>% dplyr::select(-names_drop)# %>%
  #dplyr::filter(Site == "HB3") %>% 
  #dplyr::filter(Canopy_Type == "val") %>% dim
#Write full spectra
names_drop<-colnames(as.data.frame(Canopy_image_spectra)[,7:25])
names_drop

#Resample to 5 nm, 10 nm, and 15 nm
Canopy_image_spectra_5nm<-spectrolab::resample(Canopy_image_spectra, new_bands = seq(398, 999, 5), fwhm = 5)
Canopy_image_spectra_10nm<-spectrolab::resample(Canopy_image_spectra, new_bands = seq(398, 999, 10), fwhm = 10)
Canopy_image_spectra_15nm<-spectrolab::resample(Canopy_image_spectra, new_bands = seq(398, 999, 15), fwhm = 15)

#Convert spectral objects to dataframes
Canopy_image_spectra_df <- as.data.frame(Canopy_image_spectra)
Canopy_image_spectra_5nm_df <- as.data.frame(Canopy_image_spectra_5nm)
Canopy_image_spectra_10nm_df <- as.data.frame(Canopy_image_spectra_10nm)
Canopy_image_spectra_15nm_df <- as.data.frame(Canopy_image_spectra_15nm)

###FILTER OUT EXTRANEOUS COLUMNS AND NO DATA VALUES!
#No data are pixels filtered out due to being part of the understory (no data value of '1.0000000467011e-34')
# Create a cleaned data frame by dropping columns and filtering out rows
Canopy_image_spectra_clean <- as.data.frame(Canopy_image_spectra_10nm_df) %>%
  dplyr::select(-names_drop) %>%
  dplyr::filter(rowSums(across(everything(), ~ . <= 1.000000e-30)) <= ncol(.)/2)

###CREATE KEY FIELD
Canopy_image_spectra_clean <- Canopy_image_spectra_clean %>%
  # Create the new TreeID field by concatenating the values with underscores
  mutate(TreeID = paste(Site, Species, X, Y, CC, sep = "_")) %>%
  # Move TreeID to be the second column (after the first column)
  relocate(TreeID, .after = 1)

###JOIN ATTRIBUTES USING KEY FIELD
dendro_attributes<- read.csv("./R_outputs/speclib_dendrometers/dendro_attributetable_hyp.csv")
# Select only needed fields from dendro_attributes
dendro_subset <- dendro_attributes %>%
  select(TreeID, Dndrmtr, DBH, Tag)
#Join subset with Canopy_image_spectra_clean using the TreeID key.
Canopy_image_spectra_clean <- Canopy_image_spectra_clean %>%
  left_join(dendro_subset, by = "TreeID") %>%
  relocate(Dndrmtr, DBH, Tag, .after = 7)
# Count rows with missing values in one of the join fields, e.g., Tag
missing_matches <- sum(is.na(Canopy_image_spectra_clean$Tag))
cat("Rows with missing Dndrmtr:", missing_matches, "\n")

###PRODUCE PIXEL TOTAL FOR EACH SITE!
Canopy_image_spectra_grouped <- as.data.frame(Canopy_image_spectra_clean) %>%
  dplyr::group_by(Site) %>%
  dplyr::tally()
# View and export the grouped data
print(Canopy_image_spectra_grouped)
write.csv(Canopy_image_spectra_grouped, "./R_outputs/speclib_dendrometers/dendrometer_canopy_sitepixeltotal.csv")

###PRODUCE MEDIAN REFLECTANCE VALUE FOR EACH SITE (398-999)!
# Create a vector of the wavelength column names as characters
wavelengths <- as.character(398:999)
# Group the cleaned data by 'Site', computing the median for the wavelength columns
# and taking the first value for the non-numeric columns.
Canopy_image_spectra_bysite <- Canopy_image_spectra_clean %>%
  group_by(Site) %>%
  summarize(
    across(all_of(wavelengths), ~ median(.x, na.rm = TRUE)),
    CC = first(CC),
    Species = first(Species),
    .groups = "drop"
  )
#View and export by site spectral library
print(Canopy_image_spectra_bysite)
write.csv(Canopy_image_spectra_bysite, "./R_outputs/speclib_dendrometers/dendrometer_canopy_bysite.csv")

#IGNORE BELOW LINE WHEN RERUNNING THIS CODE!
#Canopy_image_spectra_clean <- readRDS("./R_outputs/speclib_dendrometers/dataframes/dendrometer_canopy_speclib_15nm.rds")
###PRODUCE MEDIAN REFLECTANCE VALUE FOR EACH TREEID (398-999)!
#wavelengths <- as.character(seq(398, 998, by = 15)) #For 15 nm resampled
Canopy_image_spectra_bytreeid <- Canopy_image_spectra_clean %>%
  group_by(TreeID) %>%
  summarize(
    across(all_of(wavelengths), ~ median(.x, na.rm = TRUE)),
    CC = first(CC),
    Species = first(Species),
    .groups = "drop"
  )
#View and export by tree spectral library
#print(Canopy_image_spectra_bytreeid)
write.csv(Canopy_image_spectra_bytreeid, "./R_outputs/speclib_dendrometers/dendrometer_canopy_bytreeid_15nm.csv")
saveRDS(Canopy_image_spectra_bytreeid, "./R_outputs/speclib_dendrometers/dataframes/dendrometer_canopy_bytreeid_15nm.rds")

###WRITING SPECTRAL LIBRARIES TO .CSVs
#Write resampled spectral libraries to csvs
write.csv(Canopy_image_spectra_clean, "./R_outputs/speclib_dendrometers/dendrometer_canopy_speclib_15nm.csv")
#write.csv(Canopy_image_spectra_5nm, "./R_outputs/speclib_dendrometers/dendrometer_canopy_speclib_5nm.csv")
#write.csv(Canopy_image_spectra_10nm, "./R_outputs/speclib_dendrometers/dendrometer_canopy_speclib_10nm.csv")
#write.csv(Canopy_image_spectra_15nm, "./R_outputs/speclib_dendrometers/dendrometer_canopy_speclib_15nm.csv")


#write.csv(as.data.frame(Canopy_image_spectra_5nm) %>% dplyr::select(-names_drop), "./R_outputs/speclib/spruce_canopy_speclib_5nm.csv")
#write.csv(as.data.frame(Canopy_image_spectra_10nm) %>% dplyr::select(-names_drop), "./R_outputs/speclib/spruce_canopy_speclib_10nm.csv")
#write.csv(as.data.frame(Canopy_image_spectra_15nm) %>% dplyr::select(-names_drop), "./R_outputs/speclib/spruce_canopy_speclib_15nm.csv")

#Pull in master/joined .csv

#Use these lines to store variables for later use - saves files as dataframes
#canopy_image_csv<-read.csv("./R_outputs/speclib_dendrometers/dendrometer_canopy_speclib_1nm.csv")
saveRDS(Canopy_image_spectra_clean,"./R_outputs/speclib_dendrometers/dataframes/dendrometer_canopy_speclib_10nm.rds")

#Test the canopy spectra!
#tst<-rast("output/canopy_spectra/HB3_R19_val.ENVI")
#tst
#plot(tst["483.389 nm"])
#as.data.frame(Canopy_image_spectra) %>% 
  #dplyr::select(-names_drop) %>% 
  #group_by(Site, TreeID, Canopy_Type) %>% 
  #tally %>% print(n=300)