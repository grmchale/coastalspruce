source("Functions/lecospectR.R")
library(spectrolab)
library(terra)
library(raster)
library(rgdal)
library(hsdar)
library(rgeos)
library(sf)
library(doFuture)
library(tools)

# Set directory (where the cropped hyperspectral .dat files are stored)
path <- "G:/LiD-Hyp/_filtered_hyp_amoebas"

# List all .dat files of FULL canopies
imgs <- list.files(path, pattern = "\\.dat$", full.names = FALSE)
print(imgs)

#### PARSE LOOP ####
Canopy_labeled <- lapply(seq_along(imgs), function(x){ 
  # Bring in image
  tst <- terra::rast(file.path(path, imgs[x]))
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
  
  # Metadata parsing
  file_stem <- tools::file_path_sans_ext(basename(imgs[x]))           # ex) "GI_North_1m_aggregated"
  file_name <- sub("_1m_aggregated$", "", file_stem, ignore.case = TRUE) # ex) -> "GI_North"
  
  # Assign metadata (one row per spectrum)
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

# Add Amoeba_ID (first two characters of Sample) to each dataframe
Canopy_image_spectra_df$Amoeba_ID      <- substr(Canopy_image_spectra_df$Sample, 1, 2)
Canopy_image_spectra_5nm_df$Amoeba_ID  <- substr(Canopy_image_spectra_5nm_df$Sample, 1, 2)
Canopy_image_spectra_10nm_df$Amoeba_ID <- substr(Canopy_image_spectra_10nm_df$Sample, 1, 2)
Canopy_image_spectra_15nm_df$Amoeba_ID <- substr(Canopy_image_spectra_15nm_df$Sample, 1, 2)
reorder_cols <- function(df) {
  cols <- names(df)
  # put Sample and Amoeba_ID first, in that order, then everything else
  new_order <- c("Sample", "Amoeba_ID", setdiff(cols, c("Sample", "Amoeba_ID")))
  df[, new_order]
}

Canopy_image_spectra_df      <- reorder_cols(Canopy_image_spectra_df)
Canopy_image_spectra_5nm_df  <- reorder_cols(Canopy_image_spectra_5nm_df)
Canopy_image_spectra_10nm_df <- reorder_cols(Canopy_image_spectra_10nm_df)
Canopy_image_spectra_15nm_df <- reorder_cols(Canopy_image_spectra_15nm_df)


# Rename columns in the dataframe
library(dplyr)
library(stringr)

#Canopy_image_spectra_15nm_df <- Canopy_image_spectra_15nm_df %>%
#  mutate(Sample = case_when(
#    Sample %in% c("FP_West_aggregated", "FP_East_aggregated") ~ "FP",
#    Sample %in% c("GI_North_aggregated", "GI_South_aggregated") ~ "GI",
#    Sample %in% c("HI_East_aggregated", "HI_West_aggregated") ~ "HI",
#    Sample == "BI_Patch_1_aggregated" ~ "BI_South",
#    Sample %in% c("GW_East_aggregated", "GW_West_aggregated") ~ "GW",
#    Sample == "RI_GB_South_aggregated" ~ "RI_GB",
#    Sample == "RI_IP_Plot_aggregated" ~ "RI_IP",
#    TRUE ~ str_remove(Sample, "_aggregated$")
#  ))

# Check the resulting values:
unique(Canopy_image_spectra_df$Sample)

###WRITING SPECTRAL LIBRARIES TO RDS FILES
saveRDS(Canopy_image_spectra_df,"./R_outputs/speclib_amoebas_final/dataframes/amoebas_speclib_1nm.rds")
saveRDS(Canopy_image_spectra_5nm_df,"./R_outputs/speclib_amoebas_final/dataframes/amoebas_speclib_5nm.rds")
saveRDS(Canopy_image_spectra_10nm_df,"./R_outputs/speclib_amoebas_final/dataframes/amoebas_speclib_10nm.rds")
saveRDS(Canopy_image_spectra_15nm_df,"./R_outputs/speclib_amoebas_final/dataframes/amoebas_speclib_15nm.rds")

###WRITING SPECTRAL LIBRARIES TO .CSVs
write.csv(Canopy_image_spectra_clean, "./R_outputs/speclib_dendrometers/dendrometer_canopy_speclib_15nm.csv")
#write.csv(Canopy_image_spectra_5nm, "./R_outputs/speclib_dendrometers/dendrometer_canopy_speclib_5nm.csv")
#write.csv(Canopy_image_spectra_10nm, "./R_outputs/speclib_dendrometers/dendrometer_canopy_speclib_10nm.csv")
#write.csv(Canopy_image_spectra_15nm, "./R_outputs/speclib_dendrometers/dendrometer_canopy_speclib_15nm.csv")

######## -------------ADDITIONAL PROCESSING -----------------------############

###FILTER OUT EXTRANEOUS COLUMNS AND NO DATA VALUES!
#No data are pixels filtered out due to being part of the understory (no data value of '1.0000000467011e-34')
# Create a cleaned data frame by dropping columns and filtering out rows
Canopy_image_spectra_clean <- as.data.frame(Canopy_image_spectra_df) %>%
  dplyr::filter(rowSums(across(everything(), ~ . <= 1.000000e-30)) <= ncol(.)/2)

###PRODUCE PIXEL TOTAL FOR EACH SITE!
Canopy_image_spectra_grouped <- as.data.frame(Canopy_image_spectra_clean) %>%
  dplyr::group_by(Amoeba_ID) %>%
  dplyr::tally()
# View and export the grouped data
print(Canopy_image_spectra_grouped)
write.csv(Canopy_image_spectra_grouped, "./R_outputs/speclib_amoebas_final/amoeba_speclib_pixeltotal.csv")

###PRODUCE MEDIAN REFLECTANCE VALUE FOR EACH SITE (398-999)!
# Create a vector of the wavelength column names as characters
wavelengths <- as.character(398:999)
# Group the cleaned data by 'Site', computing the median for the wavelength columns
# and taking the first value for the non-numeric columns.
Canopy_image_spectra_bysite <- Canopy_image_spectra_clean %>%
  group_by(Amoeba_ID) %>%
  summarize(
    across(all_of(wavelengths), ~ median(.x, na.rm = TRUE)),
    #CC = first(CC),
    #Species = first(Species),
    .groups = "drop"
  )
#View and export by site spectral library
print(Canopy_image_spectra_bysite)
write.csv(Canopy_image_spectra_bysite, "./R_outputs/speclib_amoebas_final/amoeba_speclib_bysite_1nm.csv")
saveRDS(Canopy_image_spectra_bysite,"./R_outputs/speclib_amoebas_final/dataframes/amoebas_speclib_bysite_1nm.rds")
