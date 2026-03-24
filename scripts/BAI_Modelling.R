## WORK WITH BASAL AREA INCREMENT DATA OF RED SPRUCE CROWNS ##
###### DATA READ IN ##########
library(dplyr)
lidar_metrics <- readRDS("./data/lidar_metrics.rds")
chrono <- readRDS("./data/chrono_VIstats_metrics.rds")
# Filter out extraneous stats
chrono_filt <- chrono %>%
  select(-matches("(SD|Q25|Q75|Mean)$"))
# Join spectra + struc metrics together
struc_spec <- chrono_filt |> inner_join(lidar_metrics, by = "TreeID")
# Create new field representing tree age in years (rather than germination year)
struc_spec <- struc_spec |> mutate(age_2 = 2024 - age)
# Remove crowns with no age
struc_spec <- struc_spec |> filter(!is.na(age))