#This script crops and masks an image by a shapefile. The output is a raster of each tree canopy for that image.
#setwd()
source("Functions/lecospectR.R")

# Load packages
library(spectrolab)
library(terra)
library(raster)
library(rgdal)
library(hsdar)
library(rgeos)
library(sf)
library(doFuture)

############################# CHRONOLOGY WORKFLOW PREP ########################################
suppressPackageStartupMessages({
  library(terra)
})

polygon_shapefile <- "G:/HyperspectralUAV/Hyperspectral_Vectors/chronology_crowns/Chronology_Crowns1.shp"
out_dir           <- "G:/LiD-Hyp/chronology_shapefiles_split"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

polys <- terra::vect(polygon_shapefile)
fld   <- "HypCLIP"
stopifnot(fld %in% names(polys))

# Get attribute values as text (avoid factor/integer codes)
# Just the attribute table (no geometry)
attrs <- as.data.frame(polys)
vals_raw <- as.character(attrs[[fld]])

# ---- build grouping key (plain ASCII only) ----
strip_dat   <- function(x) tools::file_path_sans_ext(x)
sanitize_fn <- function(x) {
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- sub("^_+", "", x)
  substr(x, 1, 200)
}

keys <- sanitize_fn(strip_dat(vals_raw))

# Attach normalized key with a simple name to avoid parser quirks
polys$split_key <- keys

# Unique groups (exclude NA/empty)
groups <- sort(unique(keys[!is.na(keys) & nzchar(keys)]))
if (!length(groups)) stop("No non-empty HypCLIP groups found.")

cat("Found ", length(groups), " unique HypCLIP groups.\n", sep = "")

# ---- write one shapefile per group ----
counts <- integer(length(groups))
names(counts) <- groups

for (g in groups) {
  sel <- polys[polys[["split_key"]] == g, ]   # preserves ALL attributes
  counts[g] <- nrow(sel)
  
  out_shp <- file.path(out_dir, paste0(g, "_shapefile.shp"))
  
  terra::writeVector(sel, out_shp, filetype = "ESRI Shapefile", overwrite = TRUE)
  
  cat(sprintf("Wrote %-40s (%3d features)\n", basename(out_shp), nrow(sel)))
}

cat("\nSummary (features per output):\n")
print(sort(counts, decreasing = TRUE))
cat("\nDone. Shapefiles written to:\n", out_dir, "\n")


###########################################################################################################
########################### CROP SPRUCE CANOPIES (CHRONOLOGY WORKFLOW) ####################################

# =========================
# PACKAGES
# =========================
suppressPackageStartupMessages({
  library(terra)   # rast(), vect(), crop(), mask(), writeRaster(), writeVector()
})

# =========================
# USER PATHS
# =========================
shp_dir         <- "G:/LiD-Hyp/chronology_shapefiles_split"
raster_dir      <- "G:/LiD-Hyp/_filtered_hyp_EVI02"
out_dir         <- "./R_outputs/canopy_spectra_chronologies/"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================
# HELPERS
# =========================
# derive base key from a split shapefile name, e.g. "CC_North_shapefile.shp" -> "CC_North"
key_from_shp <- function(shp_path) {
  base <- tools::file_path_sans_ext(basename(shp_path))
  sub("_shapefile$", "", base, perl = TRUE)
}

# make safe filenames from TreeID
safe_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9_\\-\\.]+", "_", x)  # keep letters/digits/_/-/. (allow dots in TreeID)
  x <- sub("^_+", "", x)
  substr(x, 1, 200)                           # conservative for Windows
}

# =========================
# LIST INPUTS
# =========================
# all split shapefiles
shps <- list.files(shp_dir, pattern = "\\.shp$", full.names = TRUE)
if (!length(shps)) stop("No .shp files found in: ", shp_dir)

# all available .dat rasters (ENVI)
dat_paths <- list.files(raster_dir, pattern = "\\.dat$", ignore.case = TRUE, full.names = TRUE)
if (!length(dat_paths)) stop("No .dat rasters found in: ", raster_dir)

dat_keys <- tools::file_path_sans_ext(basename(dat_paths))
# map normalized keys to raster paths
raster_map <- stats::setNames(dat_paths, dat_keys)

# =========================
# MAIN LOOP
# =========================
total_written <- 0L

for (shp in shps) {
  key <- key_from_shp(shp)                  # e.g., "CC_North"
  if (!(key %in% names(raster_map))) {
    message("Skipping: no matching .dat for shapefile key '", key, "'")
    next
  }
  
  # read raster and polygons
  r_path <- raster_map[[key]]
  r <- try(terra::rast(r_path), silent = TRUE)
  if (inherits(r, "try-error")) {
    warning("Failed to read raster: ", r_path)
    next
  }
  poly <- try(terra::vect(shp), silent = TRUE)
  if (inherits(poly, "try-error")) {
    warning("Failed to read shapefile: ", shp)
    next
  }
  
  # ensure CRS alignment (reproject polygons if needed)
  if (!is.null(crs(r)) && !is.null(crs(poly)) && !terra::same.crs(r, poly)) {
    poly <- terra::project(poly, r)
  }
  
  lyr_names <- names(r)   # keep original band names
  
  # QC: need TreeID field
  if (!("TreeID" %in% names(poly))) {
    warning("Shapefile lacks 'TreeID' field: ", basename(shp), " — skipping.")
    next
  }
  
  nfeat <- nrow(poly)
  cat(sprintf("Processing %-30s  ->  %s  (%d polygons)\n",
              basename(shp), basename(r_path), nfeat))
  
  # iterate polygons
  for (i in seq_len(nfeat)) {
    tree_id_raw <- poly[i, ][["TreeID"]]
    tree_id     <- safe_name(tree_id_raw)
    if (!nzchar(tree_id) || is.na(tree_id)) tree_id <- sprintf("%s_poly%03d", key, i)
    
    out_path <- file.path(out_dir, paste0(tree_id, ".dat"))
    # avoid accidental overwrite if duplicate TreeIDs
    if (file.exists(out_path)) {
      out_path <- file.path(out_dir, paste0(tree_id, "_", i, ".dat"))
    }
    
    # crop + mask
    # (crop first for speed, then mask to polygon boundary)
    clipped <- try({
      r_crop <- terra::crop(r, poly[i, ])
      terra::mask(r_crop, poly[i, ])
    }, silent = TRUE)
    
    if (inherits(clipped, "try-error")) {
      warning("Clip failed for TreeID=", tree_id, " (", basename(shp), " idx=", i, ")")
      next
    }
    
    # restore band names
    names(clipped) <- lyr_names
    
    # write ENVI .dat + .hdr
    terra::writeRaster(clipped, out_path, overwrite = TRUE, filetype = "ENVI")
    total_written <- total_written + 1L
    
    # tidy
    rm(clipped); gc(verbose = FALSE)
  }
  
  rm(r, poly); gc(verbose = FALSE)
}

cat("\nDone. Wrote ", total_written, " clipped rasters to:\n", normalizePath(out_dir), "\n", sep = "")


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
#### Create .shp files for clipping ####
# Packages
suppressPackageStartupMessages({
  library(sf)       # works fine with 2020-era R
})

# Paths (use forward slashes on Windows in R) 
in_shp  <- "G:/HyperspectralUAV/Hyperspectral_Vectors/hyperspectral_amoebas_shp.shp"
out_dir <- "G:/HyperspectralUAV/Hyperspectral_Vectors/amoebas"

# Create output folder if needed
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Read shapefile
amoebas <- st_read(in_shp, quiet = TRUE)

# Check required field -
field_name <- "HYP_Clip"
stopifnot(field_name %in% names(amoebas))

# Build file-safe names from HYP_Clip (strip .dat and sanitize) 
strip_dat <- function(x) sub("\\.dat$", "", x, ignore.case = TRUE)
sanitize  <- function(x) {
  x2 <- gsub("[^A-Za-z0-9_\\-]+", "_", x)        # keep letters, numbers, _, -
  x2 <- sub("^_+", "", x2)                       # drop leading underscores
  x2 <- substr(x2, 1, 200)                       # shapefile base name length safety
  ifelse(nchar(x2) == 0, "feature", x2)
}

base_names <- sanitize(strip_dat(amoebas[[field_name]]))

# Ensure uniqueness (append index to duplicates)
make_unique <- function(v) {
  ave(v, v, FUN = function(x) {
    if (length(x) == 1) return(x)
    paste0(x, "_", seq_along(x))
  })
}
uniq_names <- make_unique(base_names)

# Write each feature as its own shapefile
for (i in seq_len(nrow(amoebas))) {
  feat   <- amoebas[i, , drop = FALSE]
  fname  <- paste0(uniq_names[i], ".shp")
  out_shp <- file.path(out_dir, fname)
  
  # st_write with delete_dsn=TRUE overwrites existing files with same name
  st_write(feat, dsn = out_shp, driver = "ESRI Shapefile",
           delete_dsn = TRUE, quiet = TRUE)
  message("Wrote: ", out_shp)
}

message("Done. Wrote ", nrow(amoebas), " shapefiles to: ", out_dir)

### Do the clipping ###

suppressPackageStartupMessages({
  library(terra)
  library(tools)
})

# ============================
# Paths
# ============================
filtered_img_dir <- "G:/LiD-Hyp/hyp_files"  # .dat files live here
canopies_path    <- "G:/HyperspectralUAV/Hyperspectral_Vectors/amoebas"  # split shapefiles
out_dir          <- "./R_outputs/canopy_spectra_amoebas"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================
# Gather inputs (preserve your style)
# ============================
# Gather all .dat files 
spruce_imgs <- list.files(filtered_img_dir, pattern = "\\.dat$", full.names = TRUE)
print(spruce_imgs)

# Gather all shapefiles (only those ending with .shp)
canopies_sites <- list.files(canopies_path, pattern = "\\.shp$", full.names = TRUE)

# Read in the shapefiles and assign names (without extension)
canopies <- lapply(canopies_sites, function(f) terra::vect(f))
names(canopies) <- tools::file_path_sans_ext(basename(canopies_sites))

print(canopies)

processed <- character(); skipped <- character()

# ============================
# Loop: image-by-image, mask with matching shapefile
# ============================
invisible(lapply(seq_along(spruce_imgs), function(x) {
  # Read hyperspectral image and get band names
  tst_img   <- terra::rast(spruce_imgs[x])
  tst_names <- names(tst_img)
  
  # Match shapefile by stem (image base name without extension)
  img_stem  <- tools::file_path_sans_ext(basename(spruce_imgs[x]))
  if (!img_stem %in% names(canopies)) {
    message("[SKIP] No matching shapefile for: ", img_stem)
    skipped <<- c(skipped, img_stem)
    rm(tst_img); gc()
    return(NULL)
  }
  
  # Get the corresponding canopy polygon (each file should have one feature)
  tst_quads <- canopies[[img_stem]]
  shp_name  <- img_stem
  
  # Reproject polygons to raster CRS if needed
  if (!terra::same.crs(tst_img, tst_quads)) {
    tst_quads <- terra::project(tst_quads, terra::crs(tst_img))
  }
  
  # If multiple polygons ended up inside a file, process each; otherwise this runs once
  lapply(seq_len(nrow(tst_quads)), function(i) {
    # Safe crop/mask: skip on any error (covers no-overlap cases)
    tst_crop <- try(terra::crop(tst_img, tst_quads[i, ]), silent = TRUE)
    if (inherits(tst_crop, "try-error")) {
      message("  [SKIP] Crop failed for ", img_stem, " (feature ", i, ")")
      skipped <<- c(skipped, paste0(img_stem, "_crop_", i))
      return(NULL)
    }
    tst_mask <- try(terra::mask(tst_crop, tst_quads[i, ]), silent = TRUE)
    if (inherits(tst_mask, "try-error")) {
      message("  [SKIP] Mask failed for ", img_stem, " (feature ", i, ")")
      skipped <<- c(skipped, paste0(img_stem, "_mask_", i))
      return(NULL)
    }
    
    # Set band names the terra way
    names(tst_mask) <- tst_names
    
    # Write the output using the image stem (preserves 1:1 pairing)
    out_path <- file.path(out_dir, paste0(shp_name, ".dat"))
    terra::writeRaster(tst_mask, out_path, overwrite = TRUE, filetype = "ENVI")
    message("  [WROTE] ", out_path)
    processed <<- c(processed, shp_name)
    
    rm(tst_crop, tst_mask); gc()
    NULL
  })
  
  rm(tst_img, tst_quads); gc()
  NULL
}))

# Optional: summary
cat("\nSummary\n-------\nProcessed unique outputs:", length(unique(processed)), "\n")
if (length(skipped)) {
  cat("Skipped entries:", length(skipped), "\n")
  print(unique(skipped))
}

library(terra)

# Directory with your clipped outputs
out_dir <- "./R_outputs/canopy_spectra_amoebas"
out_dats <- list.files(out_dir, pattern = "\\.dat$", full.names = TRUE)

# Function to find nearest band to wavelength
nearest_band <- function(r, target_nm) {
  wv <- suppressWarnings(as.numeric(gsub("[^0-9\\.]", "", names(r))))
  if (all(is.na(wv))) stop("Band names are not wavelengths")
  which.min(abs(wv - target_nm))
}

# Set up grid: square-ish layout
n <- length(out_dats)
ncol <- ceiling(sqrt(n))
nrow <- ceiling(n / ncol)

par(mfrow = c(nrow, ncol), mar = c(1,1,2,1))  # small margins

for (f in out_dats) {
  r <- rast(f)
  # Find RGB bands
  b_red   <- nearest_band(r, 660)
  b_green <- nearest_band(r, 560)
  b_blue  <- nearest_band(r, 470)
  
  plotRGB(r, r = b_red, g = b_green, b = b_blue,
          stretch = "lin",
          axes = FALSE, main = tools::file_path_sans_ext(basename(f)))
}



# Define the directory containing the filtered hyperspectral .dat files
filtered_img_dir <- "G:/LiD-Hyp/hyp_files"
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


