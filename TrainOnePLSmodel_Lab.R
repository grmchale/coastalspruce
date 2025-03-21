# Load required libraries
library(dplyr)
library(caret)

# Read in the spectral library and dendro attributes
spec_library <-read.csv("./R_outputs/speclib_dendrometers/dendrometer_canopy_bytreeid-speclibs_MERGED.csv") #All VI stats, all speclibs for each resample
drone_dendro<-read.csv("G:/Dendrometers/drone_dendro.csv")

# Join only the TWD column from dendro attributes to the spectral library
spec_library <- spec_library %>%
  left_join(select(drone_dendro, TreeID, TWD), by = "TreeID") %>%
  filter(!is.na(TWD)) %>%               # Remove rows without TWD
  mutate(TWD = as.numeric(TWD))         # Ensure TWD is numeric

# Identify spectral band columns
band_names <- colnames(spec_library[, 3:827])
# Spectral bands - 1 and 5nm for 400 to 800 nm, 10 and 15 nm for 800 to 1000 nm
band_names <- c(
  colnames(spec_library)[3:405],
  colnames(spec_library)[605:686],
  colnames(spec_library)[766:786],
  colnames(spec_library)[814:827]
)
print(band_names)

# Create a cleaned dataframe with only TWD and spectral bands
spec_library_clean <- spec_library %>%
  select(TWD, all_of(band_names))

# Set seed for reproducibility
set.seed(1235)

# Option 1: Use LOOCV for such a small dataset
#ctrl <- trainControl(method = "LOOCV")

# Option 2: Alternatively, use 5-fold CV (uncomment the following line to use it)
ctrl <- trainControl(method = "cv", number = 5)

# Reduce the number of PLS components to test given the small sample size
ncomp <- 10

# Fit the Partial Least Squares (PLS) regression model to predict TWD
plsFit <- train(
  TWD ~ .,
  data = spec_library_clean,
  method = "pls",
  maxit = 10000,
  trControl = ctrl,
  tuneLength = ncomp
)

# Print model summary and performance metrics
print(plsFit)

# If you prefer an independent test set, you could do an 80/20 split:
inTrain <- createDataPartition(y = spec_library_clean$TWD, p = 0.8, list = FALSE)
training <- spec_library_clean[inTrain, ]
testing <- spec_library_clean[-inTrain, ]

# Refit model on training data using the same control settings (or adjust as needed)
plsFit2 <- train(
  TWD ~ .,
  data = training,
  method = "pls",
  maxit = 10000,
  trControl = ctrl,
  tuneLength = ncomp
)

# Predict TWD on the testing set
plsFit_pred <- predict(plsFit2, newdata = testing)

# Plot predicted vs. observed TWD
plot(plsFit_pred, testing$TWD,
     xlab = "Predicted TWD",
     ylab = "Observed TWD",
     main = "PLS Regression: Predicted vs. Observed TWD")
abline(0, 1, col = "red", lwd = 2)

# Evaluate model performance with R-squared and RMSE
r2_val <- R2(plsFit_pred, testing$TWD)
rmse_val <- RMSE(plsFit_pred, testing$TWD)

print(paste("R-squared:", round(r2_val, 3)))
print(paste("RMSE:", round(rmse_val, 3)))

# Get variable importance for the fitted model
importance <- varImp(plsFit, scale = FALSE)
print(importance)

# Plot only the top 20 variables
plot(importance, top = 20, main = "Top 20 Variable Importance for PLS Model")













