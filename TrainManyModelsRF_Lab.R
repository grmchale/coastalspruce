#source("Functions/lecospectR.R")
library(tidyverse)
library(ranger)
library(caret)
library(yardstick)   # model metrics (R2, RMSE)

## Read in data
#Spectral data
speclib_rf<-read.csv("./R_outputs/speclib_dendrometers/dendrometer_canopy_bytreeid_1nm.csv") #Note 1nm, 5nm, 10nm here
canopy_VIs<-read.csv("./R_outputs/speclib_dendrometers/veg_indices/dendrometer_VIs_1nm_bytreeid.csv") #Note 1nm, 5nm, 10nm here
VIstats_merged <-read.csv("./R_outputs/speclib_dendrometers/veg_indices_stats/dendrometer_VIstats_MERGED.csv") #All VI stats for each nm resample
VIstats_speclib <-read.csv("./R_outputs/speclib_dendrometers/dendrometer_canopy_bytreeid-speclibs_MERGED.csv") #All VI stats, all speclibs for each resample
#Dendrometer data
drone_dendro<-read.csv("G:/Dendrometers/drone_dendro.csv")
cumulative_zg<-read.csv("G:/Dendrometers/cumulative_zg.csv")

## Decide which RF dataframe to use, set dataframe name to match workflow beneath, your onramp!
#Joined TWD & canopy VIs
df_rf <- canopy_VIs %>%
  inner_join(drone_dendro %>% select(Dndrmtr, TWD), by = "Dndrmtr") %>%
  drop_na(TWD)  # remove rows with no TWD value
#Joined TWD, spectral library & canopy VIs
df_rf <- speclib_rf %>% 
  inner_join(canopy_VIs, by = "TreeID", keep = FALSE) %>%  
  inner_join(drone_dendro %>% select(Dndrmtr, TWD), by = "Dndrmtr") %>%  
  drop_na(TWD) %>%
  select(-matches("X.x|CC.x|Species.x|X.1|sample_name|Site|Species.y|X.y|Y|CC.y|Dndrmtr|DBH|Tag"))
#Joined TWD, VI stats (SD, Mean, Median, Quantiles) for 1,5,10,15 nm resamples
df_rf <- VIstats_merged %>%
  inner_join(drone_dendro %>% select(TreeID, TWD), by = "TreeID") %>%
  drop_na(TWD)
#Joined TWD, VI stats (SD, Mean, Median, Quantiles) for 1,5,10,15 nm resamples + spectral library (at 1,5,10,15nm)
df_rf <- VIstats_speclib %>%
  inner_join(drone_dendro %>% select(TreeID, TWD), by = "TreeID") %>%
  drop_na(TWD) %>%
  select(-X) %>%  # Remove the column named "X" if it exists
  rename_with(~ gsub("^X", "", .), starts_with("X"))  # Remove "X" prefix from column names
# ------------------------------------------------------------------------------------------------

#Define dependent and independent variables & model name
# Identify predictor columns
predictor_vars <- names(canopy_VIs)[12:108] #For just vegetation indicies & TWD
predictor_vars <- names(df_rf)[2:207] #For vegetation indicies, spectral library, and TWD, adjust as needed
predictor_vars <- names(df_rf)[2:1941] #For VI stats, 1,5,10,15 nm
predictor_vars <- names(df_rf)[2:2766] #For VI stats, 1,5,10,15 nm, spectral libraries (1,5,10,15nm)
print(predictor_vars)
# Define response variable
response_var <- "TWD"
#RF model name (for naming later)
RF_NAME <- "VIstats_speclib-1-5-10-15nm"

# ------------------------------------------------------------------------------------------------

# Select data for modeling
model_data <- df_rf %>%
  select(all_of(c(response_var, predictor_vars))) %>%
  drop_na()  # Remove rows with missing values

#Split into training/testing sets
set.seed(1234) # For reproducibility, each time the code is run, you get the same results if the parameters are the same!
train_control <- trainControl(method = "cv", number = 5) #cv = cross validation, then 4 training subsets, 1 testing (number=5)

# Train RF model with cross-validation using caret
rf_tuned <- train(
  TWD ~ ., data = model_data,
  method = "ranger",
  tuneLength = 5,      # quickly tests a range of hyperparameters
  num.trees = 1000,
  importance = "impurity",
  trControl = train_control
)
rf_tuned$bestTune #Displays the best hyperparameters selected by caret
rf_tuned$results  #Outputs the evaluation metrics (RMSE, R², MAE) for each combination of hyperparameters

# Summary results for the best tuned model
rf_tuned$results %>%
  filter(mtry == rf_tuned$bestTune$mtry,
         min.node.size == rf_tuned$bestTune$min.node.size,
         splitrule == rf_tuned$bestTune$splitrule)

## Feature Importance from cross-validated model
importance_vals <- data.frame(
  Vegetation_Index = names(rf_tuned$finalModel$variable.importance),
  Importance = rf_tuned$finalModel$variable.importance
) %>%
  arrange(desc(Importance))

# View top 10 vegetation indices
head(importance_vals, 10)
importance_vals %>% 
  slice_max(order_by = Importance, n = 15) %>%
  ggplot(aes(reorder(Vegetation_Index, Importance), Importance)) +
  geom_col(fill = "forestgreen") +
  coord_flip() +
  labs(x = "Variable", y = "Importance", title = paste("Top 15 Variables Predicting TWD -", RF_NAME))

## -------------------------------------------------------------

## Save cross-validated model, importance, and importance graph
# Define file paths clearly
plot_filepath <- paste0("./R_outputs/modelling/rf/importance values_graphs/", RF_NAME, "_importance.png")
model_filepath <- paste0("./R_outputs/modelling/rf/rfmodel_dataframes/", RF_NAME, ".rds")
importance_filepath <- paste0("./R_outputs/modelling/rf/importance values_csvs/", RF_NAME, "_imp.csv")

# Save plot, RDS model, and importance CSV
ggsave(plot = last_plot(), filename = plot_filepath, width = 8, height = 6, dpi = 300)
saveRDS(rf_tuned, file = model_filepath)
write.csv(importance_vals, importance_filepath, row.names = FALSE)

##################################################################################

