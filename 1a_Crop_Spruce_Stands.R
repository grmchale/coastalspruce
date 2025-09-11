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

###########################################################################################################
########################### CROP SPRUCE CANOPIES (CHRONOLOGY WORKFLOW) ####################################

# =========================
# USER-DEFINED PATHS & SETTINGS
# =========================
filtered_img_dir <- "G:/LiD-Hyp/_filtered_hyp_EVI02"
polygon_shapefile <- "G:/HyperspectralUAV/Hyperspectral_Vectors/chronology_crowns/Chronology_Crowns1.shp"
polygon_id_field <- "HypCLIP"
output_dir <- "./R_outputs/canopy_spectra_chronologies/"
n_preview <- 5  # Number of previews before batch processing

# =========================
# LOAD FILES
# =========================
# List of available .dat raster files
raster_files <- list.files(filtered_img_dir, pattern = "\\.dat$", full.names = TRUE)
raster_names <- tools::file_path_sans_ext(basename(raster_files))

# Load the full canopy shapefile (multiple polygons)
canopies <- terra::vect(polygon_shapefile)

# =========================
# FILTER POLYGONS TO THOSE MATCHING AVAILABLE RASTERS
# =========================
canopies$HypCLIP_base <- tools::file_path_sans_ext(basename(canopies[[polygon_id_field]]))
valid_polygons <- canopies[canopies$HypCLIP_base %in% raster_names, ]
split_polys <- split(valid_polygons, valid_polygons$HypCLIP_base)

cat("Number of raster files with matching polygons:", length(split_polys), "\n")

# =========================
# PREVIEW FIRST n_preview CROPPED POLYGONS
# =========================
preview_names <- names(split_polys)[1:min(n_preview, length(split_polys))]

for (name in preview_names) {
  raster_path <- file.path(filtered_img_dir, paste0(name, ".dat"))
  r <- terra::rast(raster_path)
  polys <- split_polys[[name]]
  for (i in 1:length(polys)) {
    crop_i <- terra::crop(r, polys[i])
    mask_i <- terra::mask(crop_i, polys[i])
    plotRGB(mask_i, r = 133, g = 84, b = 41, stretch = "lin", main = paste(name, "-", i))
    Sys.sleep(1)  # Pause briefly between previews
  }
}

# Ask user whether to continue
cat("\nPreview complete. Proceed with batch processing? [y/n]: ")
response <- tolower(scan(what = character(), nmax = 1, quiet = TRUE))
if (response != "y") {
  stop("Batch processing aborted by user.")
}

# =========================
# MAIN CROPPING LOOP
# =========================
lapply(names(split_polys), function(name) {
  raster_path <- file.path(filtered_img_dir, paste0(name, ".dat"))
  r <- terra::rast(raster_path)
  original_names <- names(r)
  polys <- split_polys[[name]]
  
  lapply(1:length(polys), function(i) {
    clipped <- terra::mask(terra::crop(r, polys[i]), polys[i])
    bandnames(clipped) <- original_names
    
    # Get TreeID and sanitize
    tree_id <- as.character(polys[i][["TreeID"]])
    tree_id <- gsub("[^A-Za-z0-9_\\-]", "_", tree_id)
    out_name <- file.path(output_dir, paste0(tree_id, ".dat"))
    
    writeRaster(clipped, out_name, overwrite = TRUE, filetype = "ENVI")
    rm(clipped)
  })
  
  rm(r)
  gc()
})

###########################################################################################################
########################### CROP SPRUCE CANOPIES (DENDROMETER WORKFLOW) ###################################

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
