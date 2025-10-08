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
library(lidR)
library(sf)
library(dplyr)

# Define file paths
las_dir <- "G:/LiDAR/Normalized/Clipped_to_PlotsPatches"
crown_shp_path <- "G:/LiDAR/Chronology_Crowns_LiDAR/Chronology_Crowns_LiDAR1.shp"
output_dir <- "G:/LiDAR/Normalized/LiDAR_crowns"

# Read in crown polygons
crowns <- st_read(crown_shp_path)

# List all LAS files in directory
las_files <- list.files(las_dir, pattern = "\\.las$", full.names = TRUE)

# Loop over each LAS file
for (las_path in las_files) {
  # Read LAS file
  las <- readLAS(las_path)
  if (is.empty(las)) {
    warning(paste("LAS file", basename(las_path), "is empty - Skipping."))
    next
  }
  
  # Ensure CRS match
  if (st_crs(crowns) != st_crs(las)) {
    projection(las) <- st_crs(crowns)$wkt
  }
  
  # Find crowns that intersect with this LAS file’s extent
  las_bbox <- st_as_sfc(st_bbox(las))
  overlapping_crowns <- crowns[st_intersects(crowns, las_bbox, sparse = FALSE), ]
  
  if (nrow(overlapping_crowns) == 0) {
    message(paste("No overlapping crowns found for", basename(las_path)))
    next
  }
  
  # Loop through overlapping crowns and clip
  for (i in 1:nrow(overlapping_crowns)) {
    crown <- overlapping_crowns[i, ]
    tree_id <- as.character(crown$TreeID)
    
    # Clip LAS to crown polygon
    las_clip <- clip_roi(las, crown)
    
    if (is.empty(las_clip)) {
      warning(paste("No points found in TreeID", tree_id, "- Skipping."))
      next
    }
    
    # Define output path
    output_path <- file.path(output_dir, paste0(tree_id, ".las"))
    
    # Save clipped LAS
    writeLAS(las_clip, output_path)
    message(paste("Saved:", output_path))
  }
}

##################### CALCULATE RUGOSITY FOR CROWNS ############################
library(lidR)
library(raster)
library(tools)

# Define paths
las_dir <- "G:/LiDAR/Normalized/LiDAR_crowns"
dsm_dir <- file.path(las_dir, "DSMs")
dir.create(dsm_dir, showWarnings = FALSE)  # Create DSM folder if it doesn't exist

# List .las files
las_files <- list.files(las_dir, pattern = "\\.las$", full.names = TRUE)

# Initialize results dataframe
rugosity_df <- data.frame(TreeID = character(), rugosity = numeric(), stringsAsFactors = FALSE)

# Loop over each .las file (each crown)
for (las_path in las_files) {
  tree_id <- file_path_sans_ext(basename(las_path))
  
  # Read LAS
  las <- readLAS(las_path)
  if (is.empty(las)) {
    warning(paste("LAS empty for TreeID", tree_id, "- Skipping."))
    next
  }
  
  # Generate DSM (resolution = 0.1 m to preserve detail)
  dsm <- grid_canopy(las, res = 0.1, algorithm = p2r(subcircle = 0.15))
  
  # Save DSM raster
  dsm_path <- file.path(dsm_dir, paste0(tree_id, "_DSM.tif"))
  writeRaster(dsm, dsm_path, format = "GTiff", overwrite = TRUE)
  
  # Calculate vertical rugosity (standard deviation of DSM values)
  dsm_vals <- values(dsm)
  dsm_vals <- dsm_vals[!is.na(dsm_vals)]  # Remove NA
  if (length(dsm_vals) == 0) {
    warning(paste("DSM has no valid values for TreeID", tree_id))
    next
  }
  
  rugosity_val <- sd(dsm_vals)
  
  # Append to result dataframe
  rugosity_df <- rbind(rugosity_df, data.frame(TreeID = tree_id, rugosity = rugosity_val))
  
  message(paste("Processed TreeID:", tree_id, "Rugosity:", round(rugosity_val, 3)))
}

# View or export the result
# View(rugosity_df)
write.csv(rugosity_df, "G:/LiDAR/Normalized/LiDAR_crowns/crown_rugosity.csv", row.names = FALSE)

########################### AMOEBAS RUGOSITY CALCULATION ###############################
library(lidR)
library(raster)
library(tools)

# Define directories
las_dir <- "G:/LiDAR/Normalized/LiDAR_Amoebas"
dsm_dir <- file.path(las_dir, "DSMs")
dir.create(dsm_dir, showWarnings = FALSE)

# List LAS files
las_files <- list.files(las_dir, pattern = "\\.las$", full.names = TRUE)

# Initialize output dataframe
rugosity_df <- data.frame(ID = character(),
                          vertical_rugosity = numeric(),
                          horizontal_rugosity = numeric(),
                          stringsAsFactors = FALSE)

# Loop over each amoeba .las file
for (las_path in las_files) {
  amoeba_id <- file_path_sans_ext(basename(las_path))
  
  # Read LAS file
  las <- readLAS(las_path)
  if (is.empty(las)) {
    warning(paste("Empty LAS for", amoeba_id, "- Skipping."))
    next
  }
  
  # Generate DSM
  dsm <- grid_canopy(las, res = 0.25, algorithm = p2r(subcircle = 0.3))
  
  # Save DSM
  dsm_path <- file.path(dsm_dir, paste0(amoeba_id, "_DSM.tif"))
  writeRaster(dsm, dsm_path, format = "GTiff", overwrite = TRUE)
  
  # Get DSM values
  dsm_vals <- values(dsm)
  dsm_vals <- dsm_vals[!is.na(dsm_vals)]
  if (length(dsm_vals) == 0) {
    warning(paste("DSM has no valid values for", amoeba_id, "- Skipping."))
    next
  }
  
  # Vertical rugosity = SD of height
  vert_rug <- sd(dsm_vals)
  
  # Horizontal rugosity = SD of slope (in degrees)
  slope <- terrain(dsm, opt = "slope", unit = "degrees")
  slope_vals <- values(slope)
  slope_vals <- slope_vals[!is.na(slope_vals)]
  horiz_rug <- sd(slope_vals)
  
  # Append to output
  rugosity_df <- rbind(rugosity_df, data.frame(
    ID = amoeba_id,
    vertical_rugosity = vert_rug,
    horizontal_rugosity = horiz_rug
  ))
  
  message(paste("Processed", amoeba_id, 
                "| Vertical rugosity:", round(vert_rug, 3), 
                "| Horizontal rugosity:", round(horiz_rug, 3)))
}

# View or export
# View(rugosity_df)
write.csv(rugosity_df, file.path(las_dir, "amoeba_rugosity.csv"), row.names = FALSE)

################# CALCULATE POINT CLOUD METRICS FOR CROWNS #####################
library(lidR)
library(tools)
library(dplyr)

# Define directories
las_dir <- "G:/LiDAR/Normalized/LiDAR_crowns"
output_csv <- file.path(las_dir, "crown_metrics.csv")

# List all LAS files
las_files <- list.files(las_dir, pattern = "\\.las$", full.names = TRUE)

# Initialize output dataframe
all_metrics <- data.frame()

# Loop through each LAS file
for (las_path in las_files) {
  tree_id <- file_path_sans_ext(basename(las_path))
  
  # Read LAS
  las <- readLAS(las_path)
  if (is.empty(las)) {
    warning(paste("LAS file", tree_id, "is empty — skipping."))
    next
  }
  
  # Compute point cloud metrics using lidR's cloud_metrics()
  metrics <- cloud_metrics(las, .stdmetrics_z)
  
  # Add TreeID
  metrics$TreeID <- tree_id
  
  # Bind to master dataframe
  all_metrics <- bind_rows(all_metrics, metrics)
  
  message(paste("Processed:", tree_id))
}

# Reorder columns: TreeID first
all_metrics <- all_metrics %>%
  select(TreeID, everything())

# Write to CSV
write.csv(all_metrics, output_csv, row.names = FALSE)
message(paste("Saved crown metrics to:", output_csv))

################# CALCULATE AREA FOR EACH CROWN ###############################
library(sf)
library(dplyr)

# Define file paths
crown_shp_path <- "G:/LiDAR/Chronology_Crowns_LiDAR/Chronology_Crowns_LiDAR1.shp"
output_csv <- "G:/LiDAR/Normalized/LiDAR_crowns/crown_areas.csv"

# Read shapefile
crowns <- st_read(crown_shp_path)

# Ensure projected CRS for area calculation
if (st_is_longlat(crowns)) {
  stop("CRS is not projected. Please reproject to a CRS in meters before calculating area.")
}

# Calculate area and drop geometry before saving
crown_areas <- crowns %>%
  mutate(Area_m2 = as.numeric(st_area(.))) %>%
  st_drop_geometry() %>%
  select(TreeID, Area_m2)

# Write clean CSV
write.csv(crown_areas, output_csv, row.names = FALSE)
message(paste("Saved crown areas to:", output_csv))

################## LAT/LONG FOR EVERY CROWN ###########################
library(sf)
library(dplyr)

# Define paths
crown_shp_path <- "G:/LiDAR/Chronology_Crowns_LiDAR/Chronology_Crowns_LiDAR1.shp"
output_csv <- "G:/LiDAR/Latitude/chroncrown_latlon.csv"

# Read crowns
crowns <- st_read(crown_shp_path)

# Transform to WGS84 (EPSG:4326) for lat/lon
crowns_wgs <- st_transform(crowns, crs = 4326)

# Get centroids
centroids <- st_centroid(crowns_wgs)

# Extract coordinates
coords <- st_coordinates(centroids)

# Build output table
crown_latlon <- crowns_wgs %>%
  st_drop_geometry() %>%
  mutate(Long = coords[, 1],
         Lat  = coords[, 2]) %>%
  select(TreeID, Lat, Long)

# Write to CSV
write.csv(crown_latlon, output_csv, row.names = FALSE)
message(paste("Saved crown lat/lon to:", output_csv))

# Unload sf package
detach("package:sf", unload = TRUE)

################# LiDAR METRICS FOR AMOEBAS ############################
library(lidR)
library(tools)
library(dplyr)

# Define paths
las_dir <- "G:/LiDAR/Normalized/LiDAR_Amoebas"
output_csv <- file.path(las_dir, "amoeba_metrics.csv")

# List LAS files
las_files <- list.files(las_dir, pattern = "\\.las$", full.names = TRUE)

# Initialize results
all_metrics <- data.frame()

# Loop through each .las file
for (las_path in las_files) {
  amoeba_id <- file_path_sans_ext(basename(las_path))
  
  # Read LAS file
  las <- readLAS(las_path)
  if (is.empty(las)) {
    warning(paste("LAS is empty for:", amoeba_id, "- Skipping."))
    next
  }
  
  # Compute standard metrics
  metrics <- cloud_metrics(las, .stdmetrics_z)
  metrics$ID <- amoeba_id
  
  # Combine into results table
  all_metrics <- bind_rows(all_metrics, metrics)
  
  message(paste("Processed:", amoeba_id))
}

# Reorder columns: ID first
all_metrics <- all_metrics %>%
  select(ID, everything())

# Write to CSV
write.csv(all_metrics, output_csv, row.names = FALSE)
message(paste("Saved amoeba metrics to:", output_csv))

################# AMOEBA METRICS:LAT,LONG,ELV,ASPECT,SLOPE,AREA ##################
library(sf)
library(raster)
library(dplyr)
library(tools)

# --- 1. Define input/output paths ---
shapefile_path <- "G:/LiDAR/LiDAR_Shapefiles/Amoebas/lidar_amoebas.shp"
dem_dir <- "G:/Misc/DEMs/ALL"
output_csv <- "G:/LiDAR/Normalized/LiDAR_Amoebas/amoeba_latlon_topo.csv"

# --- 2. Read shapefile ---
amoebas <- st_read(shapefile_path)

# --- 3. Calculate area in square meters ---
amoebas$Area_m2 <- as.numeric(st_area(amoebas))  # Use original UTM CRS

# --- 4. Calculate centroids and transform to WGS84 for lat/lon ---
centroids <- st_centroid(amoebas)
centroids_wgs <- st_transform(centroids, crs = 4326)
coords <- st_coordinates(centroids_wgs)
amoebas$Long <- coords[, 1]
amoebas$Lat <- coords[, 2]

# --- 5. Define DEM mapping to IDs ---
dem_map <- list(
  ET = "EdgarTennisDEM.tif",
  FP = "ForbesPondDEM.tif",
  GI = "GerrishIsland_DEM.tif",
  GW = "GreatWassDEM.tif",
  HI = "HogIslandDEM.tif",
  RI = "RoqueIslandDEM.tif",
  WP = "WillardPointDEM.tif",
  CE = "CapeElizabethDEM.tif",
  CC = "CrockettCoveDEM.tif",
  BI = "CrockettCoveDEM.tif"
)

# --- 6. Initialize columns for elevation, slope, aspect ---
amoebas$Elevation <- NA
amoebas$Slope_deg <- NA
amoebas$Aspect_deg <- NA

# --- 7. Loop through rows and extract DEM-based metrics ---
for (i in 1:nrow(amoebas)) {
  row <- amoebas[i, ]
  id <- as.character(row$ID)
  
  # Lookup DEM filename
  dem_file <- dem_map[[id]]
  if (is.null(dem_file)) {
    warning(paste("No DEM mapped for ID:", id, "- Skipping."))
    next
  }
  dem_path <- file.path(dem_dir, dem_file)
  if (!file.exists(dem_path)) {
    warning(paste("DEM file not found:", dem_path))
    next
  }
  
  # Read DEM and compute slope/aspect
  dem <- raster(dem_path)
  slope <- terrain(dem, opt = "slope", unit = "degrees")
  aspect <- terrain(dem, opt = "aspect", unit = "degrees")
  
  # Get centroid in DEM's CRS
  centroid_utm <- st_coordinates(st_centroid(row))
  xy <- data.frame(x = centroid_utm[1], y = centroid_utm[2])
  
  # Extract values
  amoebas$Elevation[i]   <- extract(dem,    xy)
  amoebas$Slope_deg[i]   <- extract(slope,  xy)
  amoebas$Aspect_deg[i]  <- extract(aspect, xy)
}

# --- 8. Final selection and export ---
amoeba_summary <- st_drop_geometry(amoebas)[, c("ID", "Area_m2", "Lat", "Long", "Elevation", "Slope_deg", "Aspect_deg")]

write.csv(amoeba_summary, output_csv, row.names = FALSE)
message(paste("✅ Saved summary to:", output_csv))






