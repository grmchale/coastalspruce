############# AGGREGATE UAV PIXELS TO COARSER SPATIAL RESOLUTION ###############
library(terra)

# Paths
input_dir <- "G:/LiD-Hyp/hyp_files"
output_dir <- "G:/HyperspectralUAV/R_outputs/canopy_spectra_amoebas/aggregated_amoebas/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# List .dat files
amoebas <- list.files(input_dir, pattern = "\\.dat$", full.names = TRUE)

# Set desired resolution and aggregation factor
target_res <- 1  # desired resolution in meters
agg_factor <- round(target_res / 0.1)

# --- Test function on one .dat file ---
test_resample_aggregate <- function(test_file, agg_factor, target_res) {
  cat("Testing on file:", basename(test_file), "\n")
  
  # Read raster
  rast <- rast(test_file)
  NAflag(rast) <- NaN
  
  # Aggregate
  aggregated <- aggregate(rast, fact = agg_factor, fun = mean, na.rm = TRUE)
  plotRGB(aggregated, r=10, g=20, b=30, scale=1, stretch="lin", main="Aggregated Result (bands 10,20,30)")
  
  # Resample
  #template <- rast(extent=ext(rast), resolution=target_res, crs=crs(rast))
  #resampled <- resample(rast, template, method="bilinear")
  #plotRGB(resampled, r=10, g=20, b=30, scale=1, stretch="lin", main="Resampled Result (bands 10,20,30)")
}

# Run the test function on first .dat file only:
test_resample_aggregate(amoebas[1], agg_factor, target_res)

#########--- Loop through all files (once satisfied with the test) ---##########
for (file in amoebas) {
  cat("Processing:", basename(file), "\n")
  
  # Read raster
  rast <- rast(file)
  NAflag(rast) <- NaN
  
  # Aggregate
  aggregated <- aggregate(rast, fact = agg_factor, fun = mean, na.rm = TRUE)
  agg_outfile <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(file)), "_aggregated.dat"))
  writeRaster(aggregated, agg_outfile, overwrite = TRUE, filetype = "ENVI")
  
  # Resample
  #template <- rast(extent=ext(rast), resolution=target_res, crs=crs(rast))
  #resampled <- resample(rast, template, method = "bilinear")
  #res_outfile <- file.path(output_dir, paste0(tools::file_path_sans_ext(basename(file)), "_resampled.dat"))
  #writeRaster(resampled, res_outfile, overwrite = TRUE, filetype = "ENVI")
}

####### COUNT NUMBER OF PIXELS IN EACH FILE ############################
# Directory with aggregated .dat files
output_dir <- "G:/HyperspectralUAV/R_outputs/canopy_spectra_amoebas/aggregated_amoebas/"

# List aggregated .dat files
agg_files <- list.files(output_dir, pattern = "\\.dat$", full.names = TRUE)

# Function to count pixels (non-NA cells)
count_pixels <- function(raster_file) {
  rast <- rast(raster_file)
  
  # Count non-NA cells in the first band (assumes all bands have same NA pattern)
  pixel_count <- global(!is.na(rast[[1]]), fun = "sum", na.rm = TRUE)[1,1]
  
  return(pixel_count)
}

# Create summary dataframe
pixel_counts <- data.frame(
  File = basename(agg_files),
  Pixel_Count = sapply(agg_files, count_pixels)
)

print(pixel_counts)
