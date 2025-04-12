#This script crops and masks an image by a shapefile. The output is a raster of each tree canopy for that image.

source("Functions/lecospectR.R")

#Load in packages
library(spectrolab)
library(terra)
library(raster)
library(rgdal)
library(hsdar)
library(rgeos)
library(sf)
library(doFuture)

###Load in ALL lecospec packages??###
packages <- c(
  "gdalUtils", "rgeos", "sf", "sp", "mapview", "maptools",
  "doParallel", "randomForest", "ranger", "missForest", "caret",
  "e1071", "raster", "spectrolab", "tidyverse", "useful", "SpaDES",
  "SpaDES.tools", "Polychrome", "gplots", "rasterVis", "RColorBrewer",
  "naniar", "doSNOW", "snow", "rjson", "arrow", "RhpcBLASctl", "plotly",
  "readr", "terra", "pivottabler", "webshot", "xgboost", "permute",
  "pls", "landscapemetrics", "landscapetools", "stars", "maptiles"
)
packages <- unique(packages)
for(pkg in packages) {
  library(pkg, character.only = TRUE)
}

#Define file paths - NO .DATs ON THE END HERE?
GI_image_path = "G:/HyperspectralUAV/Hyperspectral_Clips/Understory_filtered/GI_Colby_filtered.dat"
CE_image_path = "G:/HyperspectralUAV/Hyperspectral_Clips/Understory_filtered/CE_Colby_filtered_1.dat"
HI_image_path = "G:/HyperspectralUAV/Hyperspectral_Clips/Understory_filtered/HI_Colby_filtered.dat"
CC_image_path = "G:/HyperspectralUAV/Hyperspectral_Clips/Understory_filtered/CC_Colby_filtered_1.dat"
FP_image_path = "G:/HyperspectralUAV/Hyperspectral_Clips/Understory_filtered/FP_Colby_east_filtered_2.dat"
FP2_image_path = "G:/HyperspectralUAV/Hyperspectral_Clips/Understory_filtered/FP_Colby_west_filtered_2.dat"
RI_image_path = "G:/HyperspectralUAV/Hyperspectral_Clips/Understory_filtered/RI-IP_Colby_filtered.dat"

#Load hyperspectral clips to an object
spruce_imgs<-c(
  #CC_image_path,
  #CE_image_path,
  #FP_image_path,
  #FP2_image_path,
  #GI_image_path,
  HI_image_path)#,
  #RI_image_path)

# Define the directory containing the filtered hyperspectral .dat files
filtered_img_dir <- "G:/LiD-Hyp/filtered_hyp"
# Gather all .dat files 
spruce_imgs <- list.files(filtered_img_dir, pattern = "\\.dat$", full.names = TRUE)
print(spruce_imgs)

canopies_path = "G:/LiD-Hyp/Hyperspectral_PIRU_shapefiles/"
canopies_sites<-list.files(canopies_path)[grepl("*.shp$",list.files(canopies_path))]
canopies<-lapply(1:length(canopies_sites), function(x) {terra::vect(paste(canopies_path, canopies_sites[x], sep=""))})

canopies
#, add=TRUE)

#FULL CANOPIES#
lapply(1:length(spruce_imgs), function(x) {
  tst_img <- terra::rast(spruce_imgs[x])
  tst_names <- names(tst_img)
  tst_quads <- canopies[[x]]
  # Remove any file extension from the shapefile name
  shp_name <- tools::file_path_sans_ext(names(canopies)[x])
  lapply(1:length(tst_quads), function(i) {
    tst_crop <- terra::crop(tst_img, tst_quads[i])
    tst_mask <- terra::mask(tst_crop, tst_quads[i])
    bandnames(tst_mask) <- tst_names
    # Use .dat extension for ENVI files
    writeRaster(tst_mask, paste0("./R_outputs/canopy_spectra_amoebas/", 
                                 shp_name, ".dat"), overwrite = TRUE, filetype = "ENVI")
    rm(tst_crop)
    rm(tst_mask)
  })
  rm(tst_img)
  rm(tst_quads)
  gc()
})

tst_canopy_files<-list.files("./R_outputs/canopy_spectra_amoebas/")
tst_canopy<-terra::rast(paste("./R_outputs/canopy_spectra_amoebas/", tst_canopy_files[1], sep=""))
tst_canopy_files
########################################################################################################################
##################### SPRUCE AMOEBA WORKFLOW (CANOPIES ABOVE) #####################
# Define the directory containing the filtered hyperspectral .dat files
filtered_img_dir <- "G:/LiD-Hyp/filtered_hyp"
# Gather all .dat files 
spruce_imgs <- list.files(filtered_img_dir, pattern = "\\.dat$", full.names = TRUE)
print(spruce_imgs)

# Define the directory containing the shapefiles
canopies_path <- "G:/LiD-Hyp/Hyperspectral_PIRU_shapefiles/"
# Gather all shapefiles (only those ending with .shp)
canopies_sites <- list.files(canopies_path, pattern = "\\.shp$", full.names = FALSE)
# Read in the shapefiles and assign names (without extension)
canopies <- lapply(canopies_sites, function(f) {
  terra::vect(file.path(canopies_path, f))
})
names(canopies) <- tools::file_path_sans_ext(canopies_sites)

print(canopies)

# Loop over each hyperspectral image and crop/mask it using the corresponding shapefile
lapply(seq_along(spruce_imgs), function(x) {
  # Read hyperspectral image and get band names
  tst_img <- terra::rast(spruce_imgs[x])
  tst_names <- names(tst_img)
  # Get the corresponding canopy polygons
  tst_quads <- canopies[[x]]
  # Get the shapefile name (without extension)
  shp_name <- names(canopies)[x]
  
  # Loop over each polygon in the shapefile (if multiple polygons exist)
  lapply(seq_len(nrow(tst_quads)), function(i) {
    tst_crop <- terra::crop(tst_img, tst_quads[i,])
    tst_mask <- terra::mask(tst_crop, tst_quads[i,])
    bandnames(tst_mask) <- tst_names
    # Write the output using the shapefile name (without any additional suffix)
    writeRaster(tst_mask, 
                file.path("./R_outputs/canopy_spectra_amoebas", paste0(shp_name, ".dat")),
                overwrite = TRUE, filetype = "ENVI")
    rm(tst_crop, tst_mask)
  })
  rm(tst_img, tst_quads)
  gc()
})

# Verify the outputs
tst_canopy_files <- list.files("./R_outputs/canopy_spectra_amoebas/")
print(tst_canopy_files)
tst_canopy <- terra::rast(file.path("./R_outputs/canopy_spectra_amoebas", tst_canopy_files[1]))
print(tst_canopy_files)
