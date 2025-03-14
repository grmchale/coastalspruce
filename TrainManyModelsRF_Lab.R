source("Functions/lecospectR.R")
library(tidyverse)
library(ranger)
library(caret)
library(yardstick)   # model metrics (R2, RMSE)

#Read in data
speclib_rf<-read.csv("./R_outputs/speclib_dendrometers/dendrometer_canopy_bytreeid.csv") #Note 1nm, 5nm, 10nm here
canopy_VIs<-read.csv("./R_outputs/speclib_dendrometers/veg_indices/dendrometer_VIs_1nm_bytreeid.csv") #Note 1nm, 5nm, 10nm here
drone_dendro<-read.csv("G:/Dendrometers/drone_dendro.csv")
cumulative_zg<-read.csv("G:/Dendrometers/cumulative_zg.csv")

#Make a list of  band names for data filtering later
#band_names<-colnames(speclib_rf[,30:ncol(speclib_rf)])
#VI_names<-colnames(canopy_VIs[,6:ncol(canopy_VIs)])

##Decide which RF dataframe to use, set dataframe name to match workflow beneath
#Joined TWD, spectral library & canopy VIs
df_rf <- speclib_rf %>% 
  inner_join(canopy_VIs, by = "TreeID", keep = FALSE) %>%  
  inner_join(drone_dendro %>% select(Dndrmtr, TWD), by = "Dndrmtr") %>%  
  drop_na(TWD) %>%
  select(-matches("X.x|CC.x|Species.x|X.1|sample_name|Site|Species.y|X.y|Y|CC.y|Dndrmtr|DBH|Tag"))

#Joined TWD & canopy VIs
df_rf <- canopy_VIs %>%
  inner_join(drone_dendro %>% select(Dndrmtr, TWD), by = "Dndrmtr") %>%
  drop_na(TWD)  # remove rows with no TWD value

# Identify predictor columns
predictor_vars <- names(canopy_VIs)[12:108] #For just vegetation indicies & TWD
predictor_vars <- names(df_rf)[2:640] #For vegetation indicies, spectral library, and TWD, adjust as needed
print(predictor_vars)
# Define response variable
response_var <- "TWD"
# Select data for modeling
model_data <- df_rf %>%
  select(all_of(c(response_var, predictor_vars))) %>%
  drop_na()  # Remove rows with missing values

#Split into training/testing sets
set.seed(1234) # For reproducibility
trainIndex <- createDataPartition(model_data[[response_var]], p = 0.7, list = FALSE)
train_data <- model_data[trainIndex, ]
test_data  <- model_data[-trainIndex, ]

# Train the Random Forest regression using ranger
rf_model <- ranger(
  formula = TWD ~ .,  # formula shortcut for all predictors
  data = train_data,
  num.trees = 1000,    # 500 is reasonable?
  importance = 'impurity'  # variable importance
)
print(rf_model)

# Predict on test data & evaluate performance
rf_preds <- predict(rf_model, test_data)$predictions
# Combine predictions with observed
results <- test_data %>%
  select(TWD) %>%
  mutate(Predicted = rf_preds)
# Evaluate model
metrics(results, truth = TWD, estimate = Predicted)

# Investigate feature importance
importance_vals <- data.frame(
  Vegetation_Index = names(rf_model$variable.importance),
  Importance = rf_model$variable.importance
) %>%
  arrange(desc(Importance))
# Display the top 10 vegetation indices
head(importance_vals, 10)
importance_vals %>% 
  slice_max(order_by = Importance, n = 15) %>%
  ggplot(aes(reorder(Vegetation_Index, Importance), Importance)) +
  geom_col(fill = "forestgreen") +
  coord_flip() +
  labs(x = "Vegetation Index", y = "Importance", title = "Top 15 Vegetation Indices Predicting TWD")

# Save results
write.csv(results, "RandomForest_TWD_predictions.csv", row.names = FALSE)
saveRDS(rf_model, file = "RandomForest_TWD_Model.rds")









