# --- Load dependencies ---
source("./Functions/lecospectR.R")
library(terra)

# --- Customize inputs here ---
# Index and threshold settings
Vindex <- "EVI"       # Options: "FCI1", "EVI", "VIAVG", "MSR705", "NDVI"
threshvalue <- 0.2    # For NDVI, EVI, VIAVG
threshmin   <- 3      # For MSR705
threshmax   <- 8      # For MSR705

# Band indices
bands <- list(
  red           = 144,
  red_edge      = 179,
  nir           = 219,
  blue          = 25,
  band750       = 192,
  band705       = 168,
  band445       = 27,
  nir_avg_bands = 190:241
)

# File paths
dir_envi     <- "G:/LiD-Hyp/hyp_files"
out_mask_dir <- "G:/LiD-Hyp/vi_masks/EVI02"
out_hyp_dir  <- "G:/LiD-Hyp/_filtered_hyp_EVI02"

# --- Ensure output directories exist ---
dir.create(out_mask_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_hyp_dir, recursive = TRUE, showWarnings = FALSE)

# --- Main processing function ---
test_filter_file <- function(file, index, threshold, threshmin, threshmax, bands, inspect = FALSE) {
  img <- rast(file)
  cat("Processing:", basename(file), "| Bands:", nlyr(img), "| Index:", index, "\n")
  
  # Compute vegetation index
  vi <- switch(toupper(index),
               "FCI1" = img[[bands$red]] * img[[bands$red_edge]],
               "EVI"  = {
                 red <- img[[bands$red]]
                 nir <- img[[bands$nir]]
                 blue <- img[[bands$blue]]
                 2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)
               },
               "MSR705" = {
                 r750 <- img[[bands$band750]]
                 r705 <- img[[bands$band705]]
                 r445 <- img[[bands$band445]]
                 (r750 - r445) / (r705 - r445)
               },
               "VIAVG" = app(img[[bands$nir_avg_bands]], fun = mean),
               "IRAVG" = app(img[[bands$nir_avg_bands]], fun = mean),
               "NDVI"  = {
                 red <- img[[bands$red]]
                 nir <- img[[bands$nir]]
                 (nir - red) / (nir + red)
               },
               stop("Unsupported vegetation index.")
  )
  
  # Masking logic
  mask_expr <- switch(toupper(index),
                      "FCI1"   = vi > threshold,
                      "MSR705" = vi < threshmin | vi > threshmax,
                      vi < threshold  # Default for other indices
  )
  
  veg_index_mask <- ifel(mask_expr, NA, vi)
  binary_mask    <- ifel(mask_expr, NA, 1)
  
  if (inspect) plot(veg_index_mask, main = paste("Masked", index))
  
  # Output paths
  mask_out <- file.path(out_mask_dir, paste0("mask_", tools::file_path_sans_ext(basename(file)), ".dat"))
  hyp_out  <- file.path(out_hyp_dir, basename(file))  # SAME filename as input
  
  writeRaster(veg_index_mask, mask_out, overwrite = TRUE, filetype = "ENVI")
  writeRaster(img * binary_mask, hyp_out, overwrite = TRUE, filetype = "ENVI")
  
  invisible(list(veg_index_mask = veg_index_mask, filtered_img = img * binary_mask))
}

# --- Batch processing ---
files <- list.files(dir_envi, pattern = "\\.dat$", full.names = TRUE)
if (length(files) == 0) stop("No .dat files found in: ", dir_envi)

# Inspect xth file interactively
test_filter_file(
  file      = files[8],
  index     = Vindex,
  threshold = threshvalue,
  threshmin = threshmin,
  threshmax = threshmax,
  bands     = bands,
  inspect   = TRUE
)

# Process all files
lapply(files, function(f) {
  test_filter_file(
    file      = f,
    index     = Vindex,
    threshold = threshvalue,
    threshmin = threshmin,
    threshmax = threshmax,
    bands     = bands
  )
})

######################### DISPLAY MASKED IMAGES ################################
library(terra)

# --- Display function for filtered hyperspectral images in true color ---
display_filtered_images <- function(
    directory = "G:/LiD-Hyp/_filtered_hyp_EVI02",   # Folder with filtered .dat files
    red_band = 133, green_band = 84, blue_band = 41,  # True color bands
    max_images = 9, rows = 3, cols = 3,
    stretch = TRUE
) {
  files <- list.files(directory, pattern = "\\.dat$", full.names = TRUE)
  if (length(files) == 0) stop("No .dat files found in: ", directory)
  
  # Limit to max_images and initialize
  files <- head(files, max_images)
  
  # Load and prepare RGB composites
  get_rgb <- function(file) {
    img <- rast(file)
    if (nlyr(img) < max(red_band, green_band, blue_band)) {
      warning("Skipping", basename(file), "- not enough bands.")
      return(NULL)
    }
    rgb <- c(img[[red_band]], img[[green_band]], img[[blue_band]])
    names(rgb) <- c("red", "green", "blue")
    rgb
  }
  
  rgb_images <- lapply(files, get_rgb)
  valid_idx <- which(!sapply(rgb_images, is.null))
  rgb_images <- rgb_images[valid_idx]
  file_names <- tools::file_path_sans_ext(basename(files[valid_idx]))
  
  if (length(rgb_images) == 0) stop("No valid RGB composites to display.")
  
  # Set up plotting layout
  par(mfrow = c(min(rows, ceiling(length(rgb_images)/cols)), min(cols, length(rgb_images))))
  for (i in seq_along(rgb_images)) {
    plotRGB(rgb_images[[i]], r = 1, g = 2, b = 3, stretch = if (stretch) "lin" else "none",
            main = file_names[i], axes = FALSE)
  }
  par(mfrow = c(1, 1))  # Reset layout
  cat("Displayed", length(rgb_images), "true color RGB images.\n")
}

# Run the function and display the images
display_filtered_images(
  directory = "G:/LiD-Hyp/_filtered_hyp_EVI02",
  red_band = 133,
  green_band = 84,
  blue_band = 41,
  max_images = 17,  # or however many you want to preview
  stretch = TRUE,
  rows = 6,
  cols = 3
)

