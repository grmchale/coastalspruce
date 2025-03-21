# ---------------------- CONVERT SHAPEFILES TO NAD 1983 UTM 19N ---------------------------
library(sf)

# Define input and output directories
input_folder <- "G:/LiDAR/Dendrometer_Crowns_LiDAR"
output_folder <- "G:/LiDAR/Dendrometer_Crowns_LiDAR_NAD83"

# Create the output directory if it does not exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# List all shapefiles in the input folder
shapefiles <- list.files(input_folder, pattern = "\\.shp$", full.names = TRUE)

# Define source and target CRS using EPSG codes.
# (Replace these codes if your specific CRS definitions differ.)
source_crs <- 6348   # NAD 1983 (2011) UTM Zone 19N (verify this code)
target_crs <- 26919  # NAD 1983 UTM Zone 19N

# Loop through each shapefile, transform and write to output folder
for (shp in shapefiles) {
  # Read the shapefile
  layer <- st_read(shp, quiet = TRUE)
  
  # Ensure the source CRS is assigned (if not already set in the file)
  st_crs(layer) <- source_crs
  
  # Transform the layer to the target CRS
  layer_trans <- st_transform(layer, crs = target_crs)
  
  # Build the output file path (using the same base name)
  output_file <- file.path(output_folder, basename(shp))
  
  # Write the transformed shapefile, preserving all attributes
  st_write(layer_trans, output_file, delete_layer = TRUE, quiet = TRUE)
  
  cat("Reprojected and saved:", basename(shp), "\n")
}
# ----------------------------------------------------------------------------------
# # # ELEVATION, SLOPE, ASPECT, AND DISTANCE TO OCEAN FOR EACH SHAPEFILE (BY TREEID) # # #
# Load required libraries
library(sf)
library(raster)
library(dplyr)

# ------------------ DISTANCE TO OCEAN ------------------------------------------

# MERGE ALL CROWN SHAPEFILES INTO A SINGLE SHAPEFILE!! #

# Define file paths
input_dir <- "G:/LiDAR/Dendrometer_Crowns_LiDAR_NAD83"
output_dir <- "G:/LiDAR/Dendrometer_Crowns_LiDAR_Merged"
output_file <- file.path(output_dir, "Dendrometer_Crowns_LiDAR_Merged.shp")

# Create the output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# List all shapefiles in the input directory
shapefile_list <- list.files(input_dir, pattern = "\\.shp$", full.names = TRUE)

# Read each shapefile into a list of sf objects
sf_list <- lapply(shapefile_list, st_read, quiet = TRUE)

# Merge all shapefiles using bind_rows (which handles different columns)
merged_sf <- bind_rows(sf_list)

# Check the number of features (should be 59)
cat("Total number of features:", nrow(merged_sf), "\n")

# Write the merged shapefile
st_write(merged_sf, output_file, delete_dsn = TRUE)

cat("Merged shapefile written to:", output_file, "\n")

# USE THE NEAR TOOL OVER IN ARCGIS PRO! COASTALSPRUCEPRELIMS.APRX----------

# Read IN the merged shapefile
merged_sf <- st_read("G:/LiDAR/Dendrometer_Crowns_LiDAR_Merged/Dendrometer_Crowns_LiDAR_Merged.shp", quiet = TRUE)
# Extract the attribute table (drop the geometry)
ocean_df <- st_drop_geometry(merged_sf)

# -----------------------ELEVATION, SLOPE, ASPECT-------------------------------------

# ----------- STEP 1: Define File Paths and List Files -----------

# Folder containing reprojected shapefiles
shapefile_dir <- "G:/LiDAR/Dendrometer_Crowns_LiDAR_NAD83"
# Folder containing corresponding DEM TIFs (make sure they are in the same order)
dem_dir <- "G:/Misc/DEMs/_DENDROMETER_SITES"
# Ocean coastline shapefile path

# List and sort shapefiles and DEM files
shapefile_list <- sort(list.files(shapefile_dir, pattern = "\\.shp$", full.names = TRUE))
dem_list <- sort(list.files(dem_dir, pattern = "\\.tif$", full.names = TRUE))


# ----------- STEP 2: Loop Through Shapefile/DEM Pairs and Calculate Metrics -----------

# Initialize an empty list to store results for each shapefile
esa_list <- list()

for(i in seq_along(shapefile_list)) {
  cat("\nProcessing file:", basename(shapefile_list[i]), "\n")
  
  # Read current shapefile and its corresponding DEM raster
  shp <- st_read(shapefile_list[i], quiet = TRUE)
  dem <- raster(dem_list[i])
  
  # --------------------- DEM Analysis ---------------------
  
  # Calculate slope and aspect from DEM (in degrees)
  slope_raster <- terrain(dem, opt = "slope", unit = "degrees")
  aspect_raster <- terrain(dem, opt = "aspect", unit = "degrees")
  
  # Compute centroids of each feature (using the center of geometry)
  shp_centroids <- st_centroid(shp)
  
  # Convert centroids to Spatial objects for raster extraction
  centroids_sp <- as(shp_centroids, "Spatial")
  
  # Extract elevation, slope, and aspect at each centroid location
  elev_vals <- extract(dem, centroids_sp)
  slope_vals <- extract(slope_raster, centroids_sp)
  aspect_vals <- extract(aspect_raster, centroids_sp)
  
  # --------------------- Build a Results Data Frame ---------------------
  
  # Assumes the shapefile has an attribute "TreeID" for each feature.
  df <- data.frame(
    TreeID = shp$TreeID,
    Elevation = elev_vals,
    Slope = slope_vals,      
    Aspect = aspect_vals    
  )
  
  # Optionally add a column to identify the site or source file
  df$Site <- basename(shapefile_list[i])
  
  # Append the current results to the list
  esa_list[[i]] <- df
}

# ----------- STEP 3: Combine Results -----------

# Combine all the per-file data frames into one final data frame
esa_results <- do.call(rbind, esa_list)

# View final results
print(esa_results)

#----------------- JOIN ELEVATION, SLOPE, ASPECT, AND DISTANCE TO COAST IN SINGLE DATAFRAME ------------------------#

# Join the dataframes using "TreeID"
geospatial_df <- inner_join(ocean_df, esa_results, by = "TreeID") %>%
  rename(Distance_to_Ocean = NEAR_DIST) %>%
  select(TreeID, Elevation, Slope, Aspect, Distance_to_Ocean)
head(geospatial_df)

# Save it!
write.csv(geospatial_df, "G:/LiDAR/Geospatial_Variables/elv_slope_asp_ocean_dendrometers.csv" )



