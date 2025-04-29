#setwd("G:/LiD-Hyp")
source("./Functions/lecospectR.R")
library(terra)
library(grid)
library(gridExtra)

# Set thresholds (adjustable)
evi_threshold <- 0.2
ndvi_threshold <- 0.6

# Band specifications for indices and visualization
# For EVI
nir_band_evi <- 219
red_band_evi <- 149
blue_band_evi <- 44

# For NDVI
nir_band_ndvi <- 219
red_band_ndvi <- 155

# For RGB visualization
rgb_red_band <- 133
rgb_green_band <- 84
rgb_blue_band <- 41

# Specify directory containing ENVI .dat files
dir_envi <- "G:/LiD-Hyp/_aggregated_hyp"

# Create output directories (with checks to prevent errors on rerun)
filtered_evi_dir <- "G:/LiD-Hyp/filtered_evi"
filtered_evi_ndvi_dir <- "G:/LiD-Hyp/filtered_evi_ndvi"
masks_dir <- "G:/LiD-Hyp/filtered_masks"

# Safely create directories if they don't exist
for (dir_path in c(filtered_evi_dir, filtered_evi_ndvi_dir, masks_dir)) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    cat("Created directory:", dir_path, "\n")
  } else {
    cat("Directory already exists:", dir_path, "\n")
  }
}

# Function to create RGB visualization of specified files
create_rgb_grid <- function(file_indices, directory, 
                            red_band = rgb_red_band, 
                            green_band = rgb_green_band, 
                            blue_band = rgb_blue_band,
                            title = "RGB Visualization") {
  
  # List all .dat files in the directory
  all_files <- list.files(directory, pattern = "\\.dat$", full.names = TRUE)
  
  if (length(all_files) == 0) {
    stop("No .dat files found in directory:", directory)
  }
  
  # Filter to selected files
  selected_files <- all_files[file_indices]
  valid_files <- selected_files[file.exists(selected_files)]
  
  if (length(valid_files) == 0) {
    stop("None of the specified files were found in directory:", directory)
  }
  
  # Create a 2x2 grid for plotting
  par(mfrow = c(2, 2), mar = c(2, 2, 3, 1))
  
  # Process each file
  for (i in seq_along(valid_files)) {
    file <- valid_files[i]
    cat("Visualizing file:", basename(file), "\n")
    
    # Read the image
    img <- rast(file)
    
    # Create RGB composite
    if (nlyr(img) >= max(red_band, green_band, blue_band)) {
      r <- img[[red_band]]
      g <- img[[green_band]]
      b <- img[[blue_band]]
      
      rgb_img <- c(r, g, b)
      names(rgb_img) <- c("red", "green", "blue")
      
      # Plot with a title showing the file number
      plotRGB(rgb_img, r=1, g=2, b=3, stretch="lin", 
              main=paste("File", file_indices[i]), axes=FALSE)
    } else {
      # If the file doesn't have enough bands, display a message
      plot(1, 1, type="n", axes=FALSE, xlab="", ylab="")
      text(1, 1, paste("File", file_indices[i], "\nInsufficient bands"), cex=1.2)
    }
  }
  
  # Add overall title
  mtext(title, side = 3, line = -1.5, outer = TRUE, cex = 1.5)
  
  # Reset plotting parameters
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}

# Step 1: Filter by EVI
filter_by_evi <- function(file, threshold = evi_threshold,
                          nir = nir_band_evi, 
                          red = red_band_evi, 
                          blue = blue_band_evi,
                          inspect = FALSE) {
  
  # Read the hyperspectral image
  img <- rast(file)
  file_name <- tools::file_path_sans_ext(basename(file))
  cat("\nProcessing EVI filter for:", file_name, "\n")
  
  # Calculate EVI
  nir_band <- img[[nir]]
  red_band <- img[[red]]
  blue_band <- img[[blue]]
  
  cat("  - Calculating EVI using bands NIR:", nir, "Red:", red, "Blue:", blue, "\n")
  evi <- 2.5 * (nir_band - red_band) / (nir_band + 6 * red_band - 7.5 * blue_band + 1)
  
  # Create mask based on EVI threshold
  cat("  - Applying EVI threshold of", threshold, "\n")
  evi_mask <- ifel(evi < threshold, NA, evi)
  binary_evi_mask <- ifel(evi < threshold, NA, 1)
  
  # Apply mask to image
  filtered_img <- img * binary_evi_mask
  
  # Save EVI mask
  mask_file <- file.path(masks_dir, paste0("evi_mask_", file_name, ".dat"))
  writeRaster(evi_mask, mask_file, overwrite = TRUE, filetype = "ENVI")
  
  # Save EVI filtered image
  output_file <- file.path(filtered_evi_dir, paste0("evi_", file_name, ".dat"))
  writeRaster(filtered_img, output_file, overwrite = TRUE, filetype = "ENVI")
  
  # Display the mask if requested
  if (inspect) {
    plot(evi_mask, main = paste("EVI Mask for", file_name))
  }
  
  return(list(
    filtered_img = filtered_img,
    mask = evi_mask,
    output_file = output_file
  ))
}

# Step 2: Filter EVI-filtered images by NDVI
filter_by_ndvi <- function(evi_filtered_file, threshold = ndvi_threshold,
                           nir = nir_band_ndvi, 
                           red = red_band_ndvi,
                           inspect = FALSE) {
  
  # Read the EVI-filtered image
  img <- rast(evi_filtered_file)
  file_name <- tools::file_path_sans_ext(basename(evi_filtered_file))
  # Extract the original file name without the "evi_" prefix
  original_name <- sub("^evi_", "", file_name)
  
  cat("\nProcessing NDVI filter for:", file_name, "\n")
  
  # Calculate NDVI
  nir_band <- img[[nir]]
  red_band <- img[[red]]
  
  cat("  - Calculating NDVI using bands NIR:", nir, "Red:", red, "\n")
  ndvi <- (nir_band - red_band) / (nir_band + red_band)
  
  # Create mask based on NDVI threshold
  cat("  - Applying NDVI threshold of", threshold, "\n")
  ndvi_mask <- ifel(ndvi < threshold, NA, ndvi)
  binary_ndvi_mask <- ifel(ndvi < threshold, NA, 1)
  
  # Apply mask to image
  filtered_img <- img * binary_ndvi_mask
  
  # Save NDVI mask
  mask_file <- file.path(masks_dir, paste0("ndvi_mask_", original_name, ".dat"))
  writeRaster(ndvi_mask, mask_file, overwrite = TRUE, filetype = "ENVI")
  
  # Save NDVI filtered image
  output_file <- file.path(filtered_evi_ndvi_dir, paste0("ndvi_evi_", original_name, ".dat"))
  writeRaster(filtered_img, output_file, overwrite = TRUE, filetype = "ENVI")
  
  # Display the mask if requested
  if (inspect) {
    plot(ndvi_mask, main = paste("NDVI Mask for", file_name))
  }
  
  return(list(
    filtered_img = filtered_img,
    mask = ndvi_mask,
    output_file = output_file
  ))
}

# List .dat files in the directory
files <- list.files(dir_envi, pattern = "\\.dat$", full.names = TRUE)
if (length(files) == 0) {
  stop("No .dat files found in directory:", dir_envi)
} else {
  cat("Found", length(files), "hyperspectral .dat files\n")
}

# Files to visualize (indices 1, 5, 11, 15)
vis_indices <- c(1, 5, 11, 15)
# Make sure we don't exceed the available files
vis_indices <- vis_indices[vis_indices <= length(files)]

# Test the EVI filtering on the visualization files
evi_results <- list()
for (idx in vis_indices) {
  if (idx <= length(files)) {
    evi_result <- filter_by_evi(files[idx], threshold = evi_threshold, inspect = TRUE)
    evi_results[[length(evi_results) + 1]] <- evi_result
  }
}

# Visualize the EVI filtered files
cat("\nVisualizing EVI filtered files (2x2 grid)...\n")
create_rgb_grid(vis_indices, filtered_evi_dir, 
                title = paste("EVI Filtered (threshold =", evi_threshold, ")"))

# Test the NDVI filtering on the EVI-filtered visualization files
ndvi_results <- list()
evi_filtered_files <- list.files(filtered_evi_dir, pattern = "\\.dat$", full.names = TRUE)
for (idx in vis_indices) {
  if (idx <= length(evi_filtered_files)) {
    ndvi_result <- filter_by_ndvi(evi_filtered_files[idx], threshold = ndvi_threshold, inspect = TRUE)
    ndvi_results[[length(ndvi_results) + 1]] <- ndvi_result
  }
}

# Visualize the NDVI filtered files
cat("\nVisualizing NDVI filtered files (2x2 grid)...\n")
create_rgb_grid(vis_indices, filtered_evi_ndvi_dir, 
                title = paste("EVI+NDVI Filtered (thresholds =", evi_threshold, "+", ndvi_threshold, ")"))

# Prompt user to continue with full processing
cat("\n-----------------------------------------\n")
cat("Test visualizations complete.\n")
cat("EVI threshold:", evi_threshold, "  NDVI threshold:", ndvi_threshold, "\n")
cat("Would you like to process all files? (y/n): ")
answer <- readline()

if (tolower(substr(answer, 1, 1)) == "y") {
  cat("\nProcessing all files...\n")
  
  # Step 1: Process all files with EVI filter
  for (file in files) {
    evi_result <- filter_by_evi(file, threshold = evi_threshold)
  }
  
  # Step 2: Process all EVI-filtered files with NDVI filter
  evi_filtered_files <- list.files(filtered_evi_dir, pattern = "\\.dat$", full.names = TRUE)
  for (file in evi_filtered_files) {
    ndvi_result <- filter_by_ndvi(file, threshold = ndvi_threshold)
  }
  
  cat("\nAll files have been processed.\n")
  cat("EVI filtered files are in:", filtered_evi_dir, "\n")
  cat("NDVI+EVI filtered files are in:", filtered_evi_ndvi_dir, "\n")
  cat("Masks are in:", masks_dir, "\n")
} else {
  cat("\nFull processing cancelled. Only test files were processed.\n")
}

# Function to run full process without prompt (for automation)
run_full_process <- function(evi_thresh = evi_threshold, ndvi_thresh = ndvi_threshold) {
  cat("\nRunning full automated process with EVI threshold:", evi_thresh, 
      "and NDVI threshold:", ndvi_thresh, "\n")
  
  # Process all files with EVI filter
  for (file in files) {
    filter_by_evi(file, threshold = evi_thresh, inspect = FALSE)
  }
  
  # Process all EVI-filtered files with NDVI filter
  evi_filtered_files <- list.files(filtered_evi_dir, pattern = "\\.dat$", full.names = TRUE)
  for (file in evi_filtered_files) {
    filter_by_ndvi(file, threshold = ndvi_thresh, inspect = FALSE)
  }
  
  cat("\nAll files have been processed.\n")
}

# Example: To run full process without prompt, uncomment and modify:
# run_full_process(evi_thresh = 0.25, ndvi_thresh = 0.7)