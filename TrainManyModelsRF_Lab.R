#source("Functions/lecospectR.R")
#setwd("G:/HyperspectralUAV")
library(tidyverse)
library(ranger)
library(caret)
library(dplyr)
library(yardstick)   # model metrics (R2, RMSE)

## Read in data ---------------------------------------------------------------------------------------------
#Spectral data
speclib_rf<-read.csv("./R_outputs/speclib_dendrometers/dendrometer_canopy_bytreeid_1nm.csv") #Note 1nm, 5nm, 10nm here
canopy_VIs<-read.csv("./R_outputs/speclib_dendrometers/veg_indices/dendrometer_VIs_5nm_bytreeid.csv") #Note 1nm, 5nm, 10nm here
canopy_VI_stats<-read.csv("./R_outputs/speclib_dendrometers/veg_indices_stats/dendrometer_VIstats_5nm.csv") #Note wavelength
VIstats_merged <-read.csv("./R_outputs/speclib_dendrometers/veg_indices_stats/dendrometer_VIstats_MERGED.csv") #All VI stats for each nm resample
VIstats_speclib <-read.csv("./R_outputs/speclib_dendrometers/dendrometer_canopy_bytreeid-speclibs_MERGED.csv") #All VI stats, all speclibs for each resample
#LiDAR data
LiDAR_metrics <- read.csv("G:/LiDAR/Metrics/Dendro_Crown_Metrics.csv")
#Geospatial data
easo <- read.csv("G:/LiDAR/Geospatial_Variables/elv_slope_asp_ocean_dendrometers.csv") #elevation, slope, aspect, distance to coast
latitude_df <- read.csv("G:/LiDAR/Latitude/Dendrometer_Crowns_Latitude.csv")
#In-situ data
DBH <- read.csv("./R_outputs/speclib_dendrometers/veg_indices/dendrometer_VIs_1nm_bytreeid.csv") %>%
  select(TreeID, DBH)
# ...age <- read.csv TO BE CONTINUED...
#Dendrometer + HOBO data
drone_dendro<-read.csv("G:/Dendrometers/drone_dendro.csv")
cumulative_zg<-read.csv("G:/Dendrometers/cumulative_zg.csv")
uptoflight_zg<-read.csv("G:/Dendrometers/cumulative_zg_uptoflight.csv")
uptoflight_intervals_zg<-read.csv("G:/Dendrometers/uptoflight_intervals_zg.csv")
dendro_metrics<-read.csv("G:/Dendrometers/dendro_metrics.csv")

# ------------------------------------------------------------------------------------------------

## Decide which RF dataframe to use, set dataframe
# Define dependent/response df and variable
dvarb <- "zg_fraction"
df_dep <- "uptoflight_intervals_zg"

#Joined dvarb & canopy VIs
df_rf <- canopy_VIs %>%
  inner_join(get(df_dep) %>% select(TreeID, !!dvarb), by = "TreeID") %>%
  drop_na(!!sym(dvarb))  # remove rows with no dvarb value
# Joined dvarbs & canopy VI stats (all nms)
df_rf <- VIstats_merged %>%
  inner_join(
    get(df_dep) %>% select(TreeID, zg_fraction, zg_fraction_7, zg_fraction_15, 
                           zg_fraction_30, zg_fraction_60, zg_fraction_90), 
    by = "TreeID"
  ) #%>%
  #drop_na(all_of(dvarb))
#Joined dvarb, spectral library & canopy VIs
df_rf <- speclib_rf %>% 
  #inner_join(canopy_VIs, by = "TreeID", keep = FALSE) %>%  
  inner_join(get(df_dep) %>% select(TreeID, !!dvarb), by = "TreeID") %>%  
  drop_na(!!sym(dvarb)) %>%
  select(-matches("X.x|CC.x|Species.x|X.1|sample_name|Site|Species.y|X.y|Y|CC.y|Dndrmtr|DBH|Tag"))
#Joined dvarb, VI stats (SD, Mean, Median, Quantiles) for 1,5,10,15 nm resamples
df_rf <- canopy_VI_stats %>%
  inner_join(get(df_dep) %>% select(TreeID, !!dvarb), by = "TreeID") %>%
  drop_na(!!sym(dvarb)) %>%
  
#Joined dvarb, VI stats (SD, Mean, Median, Quantiles) for 1,5,10,15 nm resamples + spectral library (at 1,5,10,15nm)
df_rf <- VIstats_speclib %>%
  inner_join(get(df_dep) %>% select(TreeID, !!dvarb), by = "TreeID") %>%
  drop_na(!!sym(dvarb)) %>%
  select(-X) %>%  # Remove the column named "X" if it exists
  rename_with(~ gsub("^X", "", .), starts_with("X"))  # Remove "X" prefix from column names

#Joined dvarb and LiDAR metrics
df_rf <- LiDAR_metrics %>% 
  inner_join(get(df_dep) %>% select(TreeID, !!dvarb), by = "TreeID") %>%
  drop_na(!!sym(dvarb))
#Joined dvarb, LiDAR & 5 nm VI stats
df_rf <- LiDAR_metrics %>%
  inner_join(canopy_VI_stats, by = "TreeID", keep = FALSE) %>%  
  inner_join(get(df_dep) %>% select(TreeID, !!dvarb), by = "TreeID") %>%  
  drop_na(!!sym(dvarb)) %>%
  select(-TreeID)
#Joined dvarb, LiDAR, merged VI stats + speclib, easo,
df_rf <- VIstats_speclib %>%
  inner_join(LiDAR_metrics, by = "TreeID") %>%
  inner_join(easo, by = "TreeID") %>%
  inner_join(latitude_df, by = "TreeID") %>%
  inner_join(DBH, by = "TreeID") %>%
  inner_join(get(df_dep) %>% select(TreeID, !!dvarb), by = "TreeID") %>%
  drop_na(!!sym(dvarb))
# ------------------------------------------------------------------------------------------------

#Define dependent and independent variables & model name
# Identify predictor columns
predictor_vars <- names(canopy_VIs)[12:108] #For just vegetation indicies & TWD
predictor_vars <- names(df_rf)[2:207] #For vegetation indicies, spectral library, and TWD, adjust as needed
predictor_vars <- names(df_rf)[2:1941] #For VI stats, 1,5,10,15 nm
predictor_vars <- names(df_rf)[2:2766] #For VI stats, 1,5,10,15 nm, spectral libraries (1,5,10,15nm)
predictor_vars <- names(df_rf)[1:56] #For only LiDAR variables
predictor_vars <- names(df_rf)[!(names(df_rf) %in% dvarb)] #Dynamic read in????

print(predictor_vars)
# Define response variable
response_var <- dvarb
#RF model name (for naming later)
RF_NAME <- "LiDAR_metrics-5nm_VIstats"

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

# ---------------------------------------------------------------------------------------------------

# Feature Importance from cross-validated model
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

