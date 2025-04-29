#setwd("G:/LiD-Hyp")
source("./Functions/lecospectR.R")
library(terra)

threshvalue <- 0.12
Vindex <- "VIAVG" # Options: "FCI1", "EVI", "VIAVG", "MSR705", or "NDVI"

# For MSR705, we need two threshold values (min and max)
threshmin <- 3  # Minimum threshold for range-based indices (like MSR705)
threshmax <- 8  # Maximum threshold for range-based indices (like MSR705)

# Define a test function for one ENVI hyperspectral file
test_filter_file <- function(file, threshold = threshvalue, 
                             threshmin = threshmin, threshmax = threshmax,
                             index = Vindex, 
                             red_band = 144, red_edge_band = 179,  # new red edge band for FCI1
                             nir_band = 219, blue_band = 25, 
                             band750 = 192, band705 = 168, band445 = 27, 
                             nir_avg_bands = 190:241, 
                             inspect = FALSE) {
  # Read in the ENVI image (assumes .dat file with accompanying header)
  img <- rast(file)
  cat("Processing file:", file, "\n")
  cat("Image has", nlyr(img), "bands.\n")
  cat("Using index:", index, "\n")
  
  # Calculate vegetation index based on the selected type
  if (toupper(index) == "FCI1") {
    # FCI1: red x red edge
    red <- img[[red_band]]
    red_edge <- img[[red_edge_band]]
    veg_index <- red * red_edge
    cat("Processing FCI1 index - masking values ABOVE threshold", threshold, "\n")
  } else if (toupper(index) == "EVI") {
    if (is.null(blue_band)) {
      stop("For EVI calculation, please specify the blue_band parameter.")
    }
    red <- img[[red_band]]
    nir <- img[[nir_band]]
    blue <- img[[blue_band]]
    veg_index <- 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
    cat("Processing EVI index - masking values BELOW threshold", threshold, "\n")
  } else if (toupper(index) == "MSR705") {
    # Ensure the required band parameters are provided
    if (is.null(band750) || is.null(band705) || is.null(band445)) {
      stop("For MSR705 calculation, please specify band750, band705, and band445 parameters.")
    }
    r750 <- img[[band750]]
    r705 <- img[[band705]]
    r445 <- img[[band445]]
    veg_index <- (r750 - r445) / (r705 - r445)
    cat("Processing MSR705 index - masking values OUTSIDE range [", threshmin, ",", threshmax, "]\n")
  } else if (toupper(index) == "VIAVG" || toupper(index) == "IRAVG") {
    if (is.null(nir_avg_bands)) {
      stop("For VIAVG/IRAVG calculation, please specify the nir_avg_bands parameter as a vector of band indices.")
    }
    # Calculate the average of all specified near-infrared bands
    veg_index <- app(img[[nir_avg_bands]], fun = mean)
    cat("Processing VIAVG/IRAVG index - masking values BELOW threshold", threshold, "\n")
  } else if (toupper(index) == "NDVI") {
    # NDVI calculation: (NIR - Red) / (NIR + Red)
    if (is.null(red_band) || is.null(nir_band)) {
      stop("For NDVI calculation, please specify both red_band and nir_band parameters.")
    }
    red <- img[[red_band]]
    nir <- img[[nir_band]]
    veg_index <- (nir - red) / (nir + red)
    cat("Processing NDVI index - masking values BELOW threshold", threshold, "\n")
  } else {
    stop("Index must be one of 'FCI1', 'EVI', 'MSR705', 'VIAVG', or 'NDVI'.")
  }
  
  # Create a VI mask based on the index type and thresholds
  if (toupper(index) == "FCI1") {
    # For FCI1: Mask out values ABOVE the threshold
    veg_index_mask <- ifel(veg_index > threshold, NA, veg_index)
    binary_mask <- ifel(veg_index > threshold, NA, 1)
  } else if (toupper(index) == "MSR705") {
    # For MSR705: Mask out values OUTSIDE the range [threshmin, threshmax]
    veg_index_mask <- ifel(veg_index < threshmin | veg_index > threshmax, NA, veg_index)
    binary_mask <- ifel(veg_index < threshmin | veg_index > threshmax, NA, 1)
  } else {
    # For EVI, VIAVG/IRAVG, NDVI: Mask out values BELOW the threshold
    veg_index_mask <- ifel(veg_index < threshold, NA, veg_index)
    binary_mask <- ifel(veg_index < threshold, NA, 1)
  }
  
  # Option to inspect the VI mask interactively
  if (inspect) {
    plot(veg_index_mask, main = paste("Masked", index, "Image"))
  }
  
  # Write the VI mask raster to disk in the filtered_masks folder
  mask_out_dir <- "G:/LiD-Hyp/filtered_masks"
  if (!dir.exists(mask_out_dir)) dir.create(mask_out_dir, recursive = TRUE)
  mask_file <- file.path(mask_out_dir, paste0("mask_", tools::file_path_sans_ext(basename(file)), ".dat"))
  writeRaster(veg_index_mask, mask_file, overwrite = TRUE, filetype = "ENVI")
  
  # Apply the binary mask to the original hyperspectral image to filter out undesired pixels
  filtered_img <- img * binary_mask
  
  # Write the filtered hyperspectral image to disk in the filtered_hyp folder
  hyp_out_dir <- "G:/LiD-Hyp/filtered_hyp"
  if (!dir.exists(hyp_out_dir)) dir.create(hyp_out_dir, recursive = TRUE)
  hyp_file <- file.path(hyp_out_dir, paste0("filtered_", tools::file_path_sans_ext(basename(file)), ".dat"))
  writeRaster(filtered_img, hyp_file, overwrite = TRUE, filetype = "ENVI")
  
  # Return both outputs for further inspection if needed
  return(list(veg_index_mask = veg_index_mask, filtered_img = filtered_img))
}

# Specify directory containing ENVI .dat files
dir_envi <- "G:/LiD-Hyp/_aggregated_hyp"

# List .dat files in the directory (adjust the pattern if necessary)
files <- list.files(dir_envi, pattern = "\\.dat$", full.names = TRUE)
print(files)

# Test the filtering function on the first file
if (length(files) > 0) {
  # You can modify these values for testing
  test_result <- test_filter_file(
    files[1], 
    threshold = threshvalue,  # Main threshold value for most indices
    threshmin = threshmin,    # Min threshold for range-based indices
    threshmax = threshmax,    # Max threshold for range-based indices
    index = Vindex,           # Index type: "FCI1", "EVI", "VIAVG", "MSR705", or "NDVI"
    red_band = 144, 
    red_edge_band = 179, 
    nir_band = 219, 
    blue_band = 25, 
    band750 = 192, 
    band705 = 168, 
    band445 = 27, 
    nir_avg_bands = 190:241, 
    inspect = TRUE
  )
  
} else {
  cat("No .dat files found in the directory:", dir_envi, "\n")
}

# Once satisfied with the test, loop over all files:
for (file in files) {
  test_filter_file(
    file, 
    threshold = threshvalue,
    threshmin = threshmin,
    threshmax = threshmax,
    index = Vindex, 
    red_band = 144, 
    red_edge_band = 179,
    nir_band = 219, 
    blue_band = 25,
    band750 = 192, 
    band705 = 168, 
    band445 = 27,
    nir_avg_bands = 190:241
  )
}

####################### DISPLAY MASKED FILES ###################################
# Load necessary libraries
library(grid)
library(gridExtra)

# Define the function to display filtered images in a grid
display_filtered_images <- function(directory = "G:/LiD-Hyp/filtered_hyp", 
                                    red_band = 133, 
                                    green_band = 84,  # As specified for true color
                                    blue_band = 41,
                                    max_images = 9,   # Maximum number of images to display
                                    stretch = TRUE,   # Apply stretch for better visualization
                                    rows = 3,         # Number of rows in the grid
                                    cols = 3) {       # Number of columns in the grid
  
  # List all .dat files in the filtered directory
  files <- list.files(directory, pattern = "\\.dat$", full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No .dat files found in the directory:", directory)
  }
  # Limit the number of files if needed
  if (length(files) > max_images) {
    cat("More than", max_images, "files found. Displaying only the first", max_images, "files.\n")
    files <- files[1:max_images]
  }
  # Function to create an RGB plot from a hyperspectral image
  create_rgb_plot <- function(file) {
    # Read the image
    img <- rast(file)
    
    # Check if the image has enough bands
    if (nlyr(img) < max(red_band, green_band, blue_band)) {
      warning("Image ", basename(file), " doesn't have enough bands for RGB display.")
      return(NULL)
    }
    # Extract RGB bands
    r <- img[[red_band]]
    g <- img[[green_band]]
    b <- img[[blue_band]]
    
    # Create RGB composite
    rgb_img <- c(r, g, b)
    names(rgb_img) <- c("red", "green", "blue")
    
    # Return the RGB composite
    return(rgb_img)
  }
  # Create a list to store all RGB images
  rgb_images <- list()
  file_names <- character(length(files))
  
  # Process each file
  for (i in seq_along(files)) {
    cat("Processing file", i, "of", length(files), ":", basename(files[i]), "\n")
    rgb_img <- create_rgb_plot(files[i])
    if (!is.null(rgb_img)) {
      rgb_images[[i]] <- rgb_img
      file_names[i] <- tools::file_path_sans_ext(basename(files[i]))
    }
  }
  # Filter out any NULL entries
  valid_indices <- which(!sapply(rgb_images, is.null))
  rgb_images <- rgb_images[valid_indices]
  file_names <- file_names[valid_indices]
  
  if (length(rgb_images) == 0) {
    stop("No valid RGB composites could be created.")
  }
  
  # Set up the plotting window
  current_rows <- min(rows, ceiling(length(rgb_images) / cols))
  current_cols <- min(cols, length(rgb_images))
  
  # Create a multi-panel plot
  par(mfrow = c(current_rows, current_cols))
  
  # Plot each RGB image with stretch
  for (i in seq_along(rgb_images)) {
    if (stretch) {
      # Apply a quantile stretch for better visualization
      plotRGB(rgb_images[[i]], r = 1, g = 2, b = 3, stretch = "lin", 
              main = file_names[i], axes = FALSE)
    } else {
      plotRGB(rgb_images[[i]], r = 1, g = 2, b = 3, 
              main = file_names[i], axes = FALSE)
    }
  }
  
  # Reset the plotting window
  par(mfrow = c(1, 1))
  
  cat("Displayed", length(rgb_images), "filtered images in true color RGB.\n")
}

# Example usage:
# Set the directory containing filtered .dat files
filtered_dir <- "G:/LiD-Hyp/filtered_hyp"

# Display the filtered images in a 3x3 grid
display_filtered_images(
  directory = filtered_dir,
  red_band = 133,   # Red band
  green_band = 84,  # Green band for true color
  blue_band = 41,   # Blue band
  max_images = 17,   # Adjust as needed
  stretch = TRUE,   # Apply linear stretch for better visualization
  rows = 6,
  cols = 3
)

# To save the grid as an image, you can wrap the function call in:
# png("filtered_images_grid.png", width = 1200, height = 1200, res = 150)
# display_filtered_images(...)
# dev.off()