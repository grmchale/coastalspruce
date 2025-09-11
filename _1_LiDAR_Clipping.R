library(lidR)
library(sf)
library(dplyr)
library(tools)

############################################################################################
############### CLIP LIDAR AMOEBAS! ############################

# Define directories
shapefile_path <- "G:/LiDAR/LiDAR_Shapefiles/Amoebas/lidar_amoebas.shp"
las_dir <- "G:/LiDAR/Normalized/Clipped_to_PlotsPatches"
output_dir <- "G:/LiDAR/Normalized/LiDAR_Amoebas"

# Read amoeba shapefile
amoebas <- st_read(shapefile_path)

# Optional: filter for a single ID if needed (comment this line out to process all)
# amoebas <- filter(amoebas, ID == "BI")

# Loop over each row of the shapefile
for (i in 1:nrow(amoebas)) {
  row <- amoebas[i, ]
  las_filename <- as.character(row$LiDAR_clip)
  output_id <- as.character(row$ID)
  
  # Full path to input LAS file
  las_path <- file.path(las_dir, las_filename)
  
  # Check if the LAS file exists
  if (!file.exists(las_path)) {
    warning(paste("LAS file not found:", las_path, "- Skipping."))
    next
  }
  
  # Read LAS file
  las <- readLAS(las_path)
  if (is.empty(las)) {
    warning(paste("LAS file", las_path, "is empty - Skipping."))
    next
  }
  
  # Ensure CRS match
  if (st_crs(row) != st_crs(las)) {
    projection(las) <- st_crs(row)$wkt
  }
  
  # Clip to the amoeba polygon
  las_clip <- clip_roi(las, row)
  
  if (is.empty(las_clip)) {
    warning(paste("No points in clipped output for", output_id, "- Skipping."))
    next
  }
  
  # Define output path
  output_path <- file.path(output_dir, paste0(output_id, ".las"))
  
  # Save the clipped LAS
  writeLAS(las_clip, output_path)
  message(paste("Saved:", output_path))
}

##########################################################################################
####################### CLIP LIDAR TO CROWNS ############################
