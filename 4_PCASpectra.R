#setwd("paintrock_leaf_spectra")
source("Functions/writeSLI_source.R")
source("Functions/lecospectR.R")
#install.packages("vegan")  # Run this if you haven't installed vegan yet
library(vegan)
require(Polychrome)
require(vegan)
require(glue)
library(dplyr)

#START HERE! DEFINE img_mat!!!
img_mat<-read.csv("G:/HyperspectralUAV/R_outputs/speclib_dendrometers/dendrometer_canopy_bytreeid.csv")

#Store df of canopy spectra by TreeID
img_df <- readRDS("G:/HyperspectralUAV/R_outputs/speclib_dendrometers/dataframes/dendrometer_canopy_bytreeid.rds")
spectra_treeid<-img_df

# Convert matrix to data frame
#img_df <- as.data.frame(img_mat)

# Remove unwanted columns and rename columns as needed
img_df <- img_df %>% 
  select(-X, -CC, -Species, -TreeID)
#Remove X from columns
img_df <- img_df %>%
  rename_with(.fn = function(names) {
    sapply(names, function(name) {
      if (grepl("^X[0-9]+$", name)) {
        num <- as.numeric(sub("^X", "", name))
        if (num >= 398 & num <= 999) {
          return(sub("^X", "", name))
        }
      }
      name
    })
  })
#Convert columns to numeric
img_df <- apply(img_df, 2, function(x) as.numeric(as.character(x)))
img_df <- img_df[,2:603]
# Convert back to matrix if needed
img_mat <- as.matrix(img_df)

#saveRDS(img_df,"G:/HyperspectralUAV/R_outputs/speclib_dendrometers/dataframes/dendrometer_canopy_bytreeid.rds")
head(img_mat)

#ASK PETER ABOUT MULTIVARIATE ANALYSIS HERE!!!

#Multivariate analysis of PFT groups 
TreeID_adonis<-adonis2(img_mat~as.factor(MASTER_TreeID_speclib$Site), method="euclidean", permutations=100)
#tree_genus_adonis<-adonis2(img_mat~as.factor(tree_spectra$taxon_genus), method="euclidean", permutations=100)
TreeID_adonis
#
#> tree_adonis
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 100
#
#adonis2(formula = img_mat ~ as.factor(tree_spectra$taxon_code), permutations = 100, method = "euclidean")
#                                    Df SumOfSqs     R2      F   Pr(>F)   
#as.factor(tree_spectra$taxon_code)  26  8886255 0.5734 13.338 0.009901 **
#Residual                           258  6611102 0.4266                   
#Total                              284 15497357 1.0000                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#> tree_genus_adonis
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 100
#
#adonis2(formula = img_mat ~ as.factor(tree_spectra$taxon_genus), permutations = 100, method = "euclidean")
#                                     Df SumOfSqs      R2      F   Pr(>F)   
#as.factor(tree_spectra$taxon_genus)  18  8268942 0.53357 16.905 0.009901 **
#Residual                            266  7228416 0.46643                   
#Total                               284 15497357 1.00000                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

colnames(img_mat)

#img_mat<-as.numeric(img_mat[1:nrow(img_mat),])

###BUILD PCA WITH AND W/O SQRT TRANSFORM###
img_pca<-prcomp(img_mat, center = TRUE, scale =TRUE) #, center=FALSE, scale=FALSE)
img_pca_pr<-prcomp(img_mat)#, center=FALSE, scale=FALSE)
# Extract the first two characters from TreeID to represent the site
site_codes <- substring(spectra_treeid$TreeID, 1, 2)
# Determine unique site codes (e.g., CE, ET, RI, CC, GI, HI)
unique_sites <- unique(site_codes)
# Assign a distinct color to each site code
my_colors <- c("red", "green", "blue", "black", "purple", "orange")[1:length(unique_sites)]
# Match each row's site code to the appropriate color
point_colors <- my_colors[match(site_codes, unique_sites)]
### Assign plotting symbols (optional, can be changed as desired)
#my_symbols <- 1:length(unique_sites)
#point_symbols <- my_symbols[match(site_codes, unique_sites)]

# Calculate PCA scores for the first two components
pca_scores <- vegan::scores(img_pca)[, 2:3]
# Plot to the RStudio Plots window
plot(pca_scores, 
     col = point_colors,
     pch = 16,
     main = "Reflectance by Site, Dendrometer Crowns - PCA Axes 1 & 2")
legend("topright", 
       legend = unique_sites, 
       col = my_colors, 
       pch = 16, 
       cex = 1)
# Export!
#dev.copy(jpeg, "./R_outputs/figures/dendrometers/PCA_dendrometers_1_2_v1.jpg")
#dev.off()


###EXAMINING THE PCA COMPONENTS!!!!

#Set up PCA path
pca_path <- "./R_outputs/pca/dendrometers_pca"
## Variance explained by each principal component
summary(img_pca)
summary(img_pca_pr)

##Eigenvectors
img_pca$rotation

pca_wavelengths <- (as.data.frame(img_pca$rotation))
write.csv(pca_wavelengths, file.path(pca_path,"pca_eigenvectors_wavelengths.csv"), row.names = TRUE)

##Eigenvalues
eigenvalues <- (as.data.frame(img_pca$sdev)^2)
colnames(eigenvalues)[1:2] <- c("PCA Band", "Eigenvalue")
write.csv(eigenvalues, file.path(pca_path,"pca_eigenvalues.csv"), row.names = TRUE)

##The principal component scores for each observation in img_mat
head(img_pca$x)
pca_scores <- (as.data.frame(img_pca$x))
pca_scores <- cbind(TreeID = spectra_treeid[, 1], pca_scores)
pca_scores$Site <- substr(spectra_treeid$TreeID, 1, 2)
pca_scores <- pca_scores[, c("TreeID", "Site", setdiff(names(pca_scores), c("TreeID", "Site")))]

write.csv(pca_scores, file.path(pca_path,"pca_scores.csv"), row.names = TRUE)

##PCA plot Axes 2 vs 3
#jpeg("output/PCA_trees_Axes23.jpg")
#plot(vegan::scores(img_pca)[,2:3], col=tree_spectra$Color, pch=tree_spectra$ColorNum)
#title(main="PCA of PFT Reflectance")
#legend(x = -55, y =-10, legend=unique(tree_spectra$taxon_genus), lty=1,  pch = unique(tree_spectra$ColorNum), col=unique(tree_spectra$Color), cex=1)
#dev.off()
