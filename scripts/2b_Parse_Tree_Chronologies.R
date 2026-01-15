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
path <- "./R_outputs/canopy_spectra_chronologies/"

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
  
  # Use the .dat file name (without the extension) as metadata.
  file_name <- gsub("\\.dat$", "", imgs[x])
  
  # Correct metadata assignment (using `nrow(df)` instead of `length(df)`)
  meta(df) <- data.frame(Sample = rep(file_name, times = nrow(df)))
  
  return(df)
})

# Feed spectral object to canopy_image_spectra
Canopy_image_spectra <- Reduce(spectrolab::combine, Canopy_labeled)
Canopy_image_spectra
dim(as.data.frame(Canopy_image_spectra))

################################### RESAMPLE AND SAVE DATAFRAMES ####################################
#Resample to 5 nm, 10 nm, and 15 nm
Canopy_image_spectra_5nm<-spectrolab::resample(Canopy_image_spectra, new_bands = seq(398, 999, 5), fwhm = 5)
Canopy_image_spectra_10nm<-spectrolab::resample(Canopy_image_spectra, new_bands = seq(398, 999, 10), fwhm = 10)
Canopy_image_spectra_15nm<-spectrolab::resample(Canopy_image_spectra, new_bands = seq(398, 999, 15), fwhm = 15)

#Convert spectral objects to dataframes
#Canopy_image_spectra_df <- as.data.frame(Canopy_image_spectra)
Canopy_image_spectra_df <- as.data.frame(Canopy_image_spectra)
Canopy_image_spectra_5nm_df <- as.data.frame(Canopy_image_spectra_5nm)
Canopy_image_spectra_10nm_df <- as.data.frame(Canopy_image_spectra_10nm)
Canopy_image_spectra_15nm_df <- as.data.frame(Canopy_image_spectra_15nm)

# Rename columns in the dataframe
library(dplyr)
library(stringr)

# Check the resulting values:
unique(Canopy_image_spectra_df$Sample)

### WRITING SPECTRAL LIBRARIES TO RDS FILES
saveRDS(Canopy_image_spectra_df,      "./R_outputs/speclib_chronologies/dataframes/chronologies_speclib_1nm.rds")
saveRDS(Canopy_image_spectra_5nm_df,  "./R_outputs/speclib_chronologies/dataframes/chronologies_speclib_5nm.rds")
saveRDS(Canopy_image_spectra_10nm_df, "./R_outputs/speclib_chronologies/dataframes/chronologies_speclib_10nm.rds")
saveRDS(Canopy_image_spectra_15nm_df, "./R_outputs/speclib_chronologies/dataframes/chronologies_speclib_15nm.rds")

### WRITING SPECTRAL LIBRARIES TO .CSVs
write.csv(Canopy_image_spectra,      "./R_outputs/speclib_chronologies/chronologies_speclib_1nm.csv")
write.csv(Canopy_image_spectra_5nm,  "./R_outputs/speclib_chronologies/chronologies_speclib_5nm.csv")
write.csv(Canopy_image_spectra_10nm, "./R_outputs/speclib_chronologies/chronologies_speclib_10nm.csv")
write.csv(Canopy_image_spectra_15nm, "./R_outputs/speclib_chronologies/chronologies_speclib_15nm.csv")

######## -------------ADDITIONAL PROCESSING -----------------------############

###FILTER OUT EXTRANEOUS COLUMNS AND NO DATA VALUES!
#No data are pixels filtered out due to being part of the understory (no data value of '1.0000000467011e-34')
# Create a cleaned data frame by dropping columns and filtering out rows
#Canopy_image_spectra_clean <- as.data.frame(Canopy_image_spectra_df) %>%
 # dplyr::filter(rowSums(across(everything(), ~ . <= 1.000000e-30)) <= ncol(.)/2)

###PRODUCE PIXEL TOTAL FOR EACH TREEID!
Canopy_image_spectra_grouped <- as.data.frame(Canopy_image_spectra_df) %>%
  dplyr::group_by(Sample) %>%
  dplyr::tally()
# View and export the grouped data
print(Canopy_image_spectra_grouped)
write.csv(Canopy_image_spectra_grouped, "./R_outputs/speclib_chronologies/chronologies_speclib_pixeltotal.csv")

###PRODUCE MEDIAN REFLECTANCE VALUE FOR EACH SITE (398-999)!
# Create a vector of the wavelength column names as characters
wavelengths <- as.character(398:999)
# Group the cleaned data by 'Site', computing the median for the wavelength columns
# and taking the first value for the non-numeric columns.
Canopy_image_spectra_bytree <- Canopy_image_spectra_df %>%
  group_by(Sample) %>%
  summarize(
    across(all_of(wavelengths), ~ median(.x, na.rm = TRUE)),
    #CC = first(CC),
    #Species = first(Species),
    .groups = "drop"
  )
#View and export by site spectral library
print(Canopy_image_spectra_bytree)
write.csv(Canopy_image_spectra_bytree, "./R_outputs/speclib_chronologies/chronologies_speclib_bytree_1nm.csv")
saveRDS(Canopy_image_spectra_bytree,"./R_outputs/speclib_chronologies/dataframes/chronologies_speclib_bytree_1nm.rds")
