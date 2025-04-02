#setwd("G:/LiD-Hyp")
source("./Functions/lecospectR.R")
library(terra)

threshvalue <- 0.2
Vindex <- "EVI" # Options: "NDVI", "EVI", or "MSR705"

# Define a test function for one ENVI hyperspectral file with options for NDVI, EVI, or MSR705
test_filter_file <- function(file, threshold = threshvalue, index = Vindex, 
                             red_band = 155, nir_band = 219, blue_band = 25, 
                             band750 = 192, band705 = 168, band445 = 27, nir_avg_bands = 190:241, 
                             inspect = FALSE) {
  # Read in the ENVI image (assumes .dat file with accompanying header)
  img <- rast(file)
  cat("Processing file:", file, "\n")
  cat("Image has", nlyr(img), "bands.\n")
  
  # Calculate vegetation index based on the selected type
  if (toupper(index) == "NDVI") {
    red <- img[[red_band]]
    nir <- img[[nir_band]]
    veg_index <- (nir - red) / (nir + red)
  } else if (toupper(index) == "EVI") {
    if (is.null(blue_band)) {
      stop("For EVI calculation, please specify the blue_band parameter.")
    }
    red <- img[[red_band]]
    nir <- img[[nir_band]]
    blue <- img[[blue_band]]
    veg_index <- 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
  } else if (toupper(index) == "MSR705") {
    # Ensure the required band parameters are provided
    if (is.null(band750) || is.null(band705) || is.null(band445)) {
      stop("For mSR705 calculation, please specify band750, band705, and band445 parameters.")
    }
    r750 <- img[[band750]]
    r705 <- img[[band705]]
    r445 <- img[[band445]]
    veg_index <- (r750 - r445) / (r705 - r445)
  } else if (toupper(index) == "IRAVG") {
    if (is.null(nir_avg_bands)) {
      stop("For IRavg calculation, please specify the nir_avg_bands parameter as a vector of band indices.")
    }
    # Calculate the average of all specified near-infrared bands
    veg_index <- app(img[[nir_avg_bands]], fun = mean)
  } else {
    stop("Index must be one of 'NDVI', 'EVI', 'MSR705', or 'IRavg'.")
  }
  
  # Create a VI mask: pixels meeting the threshold retain the veg_index; others become NA
  veg_index_mask <- ifel(veg_index < threshold, NA, veg_index)
  
  # Option to inspect the VI mask interactively
  if (inspect) {
    plot(veg_index_mask, main = paste("Masked", index, "Image"))
  }
  
  # Write the VI mask raster to disk in the filtered_masks folder
  mask_out_dir <- "G:/LiD-Hyp/filtered_masks"
  if (!dir.exists(mask_out_dir)) dir.create(mask_out_dir, recursive = TRUE)
  mask_file <- file.path(mask_out_dir, paste0("mask_", tools::file_path_sans_ext(basename(file)), ".dat"))
  writeRaster(veg_index_mask, mask_file, overwrite = TRUE, filetype = "ENVI")
  
  # Create a binary mask: pixels meeting the threshold become 1, others NA
  binary_mask <- ifel(veg_index < threshold, NA, 1)
  
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
dir_envi <- "G:/LiD-Hyp/hyp_files"

# List .dat files in the directory (adjust the pattern if necessary)
files <- list.files(dir_envi, pattern = "\\.dat$", full.names = TRUE)
print(files)

# Test the filtering function on the first file
if (length(files) > 0) {
  # For NDVI filtering (adjust parameters as needed)
  test_result <- test_filter_file(files[1], threshold = threshvalue, index = Vindex, 
                                  red_band = 155, nir_band = 219, blue_band = 25, 
                                  band750 = 192, band705 = 168, band445 = 27, nir_avg_bands = 190:241, 
                                  inspect = TRUE)
  
} else {
  cat("No .dat files found in the directory:", dir_envi, "\n")
}

# Once satisfied with the test, loop over all files:
for (file in files) {
  test_filter_file(file, threshold = threshvalue, index = Vindex, 
                   red_band = 155, nir_band = 219, blue_band = 25,
                   band750 = 192, band705 = 168, band445 = 27)
}

# Display all of the files:
dev.new()
par(mfrow = c(1, 1))

hyp_out_dir <- "G:/LiD-Hyp/filtered_hyp"
filtered_files <- list.files(hyp_out_dir, pattern = "\\.dat$", full.names = TRUE)

cat("Displaying filtered hyperspectral images from", hyp_out_dir, "\n")

for (f in filtered_files) {
  cat("Displaying:", f, "\n")
  img_filtered <- rast(f)
  plot(img_filtered, main = paste("Filtered Image:", basename(f)))
  invisible(readline(prompt = "Press [enter] to continue to the next image"))
}

