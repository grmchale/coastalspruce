############# AGGREGATE UAV PIXELS TO 1 m RESOLUTION (Using Aggregate) ###############
library(terra)

# Define input and output directories
input_dir <- "./R_outputs/canopy_spectra_amoebas"
output_dir <- "G:/LiD-Hyp/_aggregated_hyp_1m"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# List all .dat files
amoebas <- list.files(input_dir, pattern = "\\.dat$", full.names = TRUE)
print(amoebas)

# Set aggregation factor (assuming original resolution ~0.1 m, adjust if needed)
agg_factor <- 10  # aggregate factor of 10 for 1 m resolution

# Test aggregation on the first raster to visually check results
rast_orig_test <- rast(amoebas[1])
NAflag(rast_orig_test) <- NaN

# Aggregate the first raster
aggregated_test <- aggregate(rast_orig_test, fact = agg_factor, fun = mean, na.rm = TRUE)

# Plot RGB composite to check aggregation
plotRGB(aggregated_test, r=10, g=20, b=30, scale=1, stretch="lin", main="Aggregated Test Raster")

# Loop through each file and aggregate
for (file in amoebas) {
  cat("Aggregating:", basename(file), "\n")
  
  # Read raster
  rast_orig <- rast(file)
  NAflag(rast_orig) <- NaN
  
  # Aggregate raster
  aggregated_rast <- aggregate(rast_orig, fact = agg_factor, fun = mean, na.rm = TRUE)
  
  # Define output filename and write raster
  aggregated_outfile <- file.path(
    output_dir,
    paste0(tools::file_path_sans_ext(basename(file)), "_1m_aggregated.dat")
  )
  
  writeRaster(aggregated_rast, aggregated_outfile, overwrite = TRUE, filetype = "ENVI")
}

cat("Aggregation complete!", "\n")

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
