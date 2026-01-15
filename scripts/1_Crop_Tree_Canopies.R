#This script crops and masks an image by a shapefile. The output is a raster of each tree canopy for that image.

#Set up lecospec functions
source("./Scripts/Utilities/install_dependencies.R")
source("Functions/lecospectR.R")

#View warnings
warnings()

#Install packages that failed in first run (spectrolab, hsdar, and doFuture)
install.packages("devtools")
library(devtools)

#rgdal, get the fricken hell in here
install.packages("rgdal")
install_github("cran/rgdal")
library(rgdal)

# Install spectrolab and hsdar from GitHub
install_github("meireles/spectrolab")
install_github("cran/hsdar")
install("hsdar")
install.packages("hsdar")
install_github("HenrikBengtsson/doFuture")
install_github("cran/rgeos")

#Load in packages
library(spectrolab)
library(terra)
library(raster)
library(rgdal)
library(hsdar)
library(rgeos)
library(sf)
library(doFuture)

###Load in ALL lecospec packages??###
packages <- c(
  "gdalUtils", "rgeos", "sf", "sp", "mapview", "maptools",
  "doParallel", "randomForest", "ranger", "missForest", "caret",
  "e1071", "raster", "spectrolab", "tidyverse", "useful", "SpaDES",
  "SpaDES.tools", "Polychrome", "gplots", "rasterVis", "RColorBrewer",
  "naniar", "doSNOW", "snow", "rjson", "arrow", "RhpcBLASctl", "plotly",
  "readr", "terra", "pivottabler", "webshot", "xgboost", "permute",
  "pls", "landscapemetrics", "landscapetools", "stars", "maptiles"
)
packages <- unique(packages)
for(pkg in packages) {
  library(pkg, character.only = TRUE)
}

#Define file paths - NO .DATs ON THE END HERE?
#GI_image_path = "G:/HyperspectralUAV/Hyperspectral_Clips/Understory_filtered/GI_Colby_filtered.dat"
#CE_image_path = "G:/HyperspectralUAV/Hyperspectral_Clips/Understory_filtered/CE_Colby_filtered_1.dat"
HI_image_path = "G:/HyperspectralUAV/Hyperspectral_Clips/Understory_filtered/HI_Colby_filtered.dat"
#CC_image_path = "G:/HyperspectralUAV/Hyperspectral_Clips/Understory_filtered/CC_Colby_filtered_1.dat"
#FP_image_path = "G:/HyperspectralUAV/Hyperspectral_Clips/Understory_filtered/FP_Colby_east_filtered_2.dat"
#FP2_image_path = "G:/HyperspectralUAV/Hyperspectral_Clips/Understory_filtered/FP_Colby_west_filtered_2.dat"
#RI_image_path = "G:/HyperspectralUAV/Hyperspectral_Clips/Understory_filtered/RI-IP_Colby_filtered.dat"

#Load hyperspectral clips to an object
spruce_imgs<-c(
  #CC_image_path,
  #CE_image_path,
  #FP_image_path,
  #FP2_image_path,
  #GI_image_path,
  HI_image_path)#,
  #RI_image_path)

canopies_path = "G:/HyperspectralUAV/Hyperspectral_Vectors/hi_crowns/"
canopies_sites<-list.files(canopies_path)[grepl("*.shp$",list.files(canopies_path))]
canopies<-lapply(1:length(canopies_sites), function(x) {terra::vect(paste(canopies_path, canopies_sites[x], sep=""))})

canopies
#, add=TRUE)

#FULL CANOPIES#
lapply(1:length(spruce_imgs),  
       function(x) {
#x=1
         tst_img <- terra::rast(spruce_imgs[x])
         tst_names<-names(tst_img)
         tst_quads<-canopies[[x]]
         #metadata(tst_mask)<-tst_quads$CLASS_NAME
         lapply(1:length(tst_quads), function (i) {
         #i=1   
         tst_crop <- terra::crop(tst_img, tst_quads[i])
           tst_mask <- terra::mask(tst_crop, tst_quads[i])
           bandnames(tst_mask)<-tst_names
           writeRaster(tst_mask, paste("./R_outputs/canopy_spectra_dendrometers/", tst_quads[i]$TreeID, ".ENVI",sep=""), overwrite = TRUE)
           rm(tst_crop)
           rm(tst_mask)})
         rm(tst_img)
         rm(tst_quads)
         gc()})

tst_canopy_files<-list.files("R_outputs/canopy_spectra_dendrometers/")
tst_canopy<-terra::rast(paste("R_outputs/canopy_spectra_dendrometers/", tst_canopy_files[1], sep=""))
tst_canopy_files
