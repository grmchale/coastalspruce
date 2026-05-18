# Read back in
lidar_metrics <- readRDS("./data/lidar_metrics.rds")
chrono <- readRDS("./data/chrono_VIstats_metrics.rds")
dendro_spectra <- read.csv("./data/dendro_spectra_joined.csv")
#install.packages("dplyr")
library(dplyr)
# Filter out extraneous stats
chrono_filt <- chrono %>%
  select(-matches("(SD|Q25|Q75|Mean)$"))
# Join spectra + struc metrics together
struc_spec <- chrono_filt |> inner_join(lidar_metrics, by = "TreeID")
# Create new field representing tree age in years (rather than germination year)
struc_spec <- struc_spec |> mutate(age_2 = 2024 - age)


######### FIND NUMBER OF RED SPRUCE CROWNS PER SITE WITH EACH METRIC #########
library(dplyr)

# a) Sample sizes for Dendrometers
dndrmtr_by_site <- dendro_spectra %>%
  mutate(Site = substr(TreeID, 1, 2)) %>%
  group_by(Site) %>%
  summarise(n_Dndrmtr = n())

dndrmtr_total <- nrow(dendro_spectra)

# b) Sample sizes for BAI_2024
bai_by_site <- struc_spec %>%
  group_by(Site) %>%
  summarise(n_BAI_2024 = sum(!is.na(BAI_2024)))

bai_total <- sum(!is.na(struc_spec$BAI_2024))

# c) Sample sizes for age
age_by_site <- struc_spec %>%
  group_by(Site) %>%
  summarise(n_age = sum(!is.na(age)))

age_total <- sum(!is.na(struc_spec$age))

# View per-site results
dndrmtr_by_site
bai_by_site
age_by_site

# View totals
cat("Total rows with Dndrmtr:", dndrmtr_total, "\n")
cat("Total rows with BAI_2024:", bai_total, "\n")
cat("Total rows with age:", age_total, "\n")

#### COMBINE TABLES INTO ONE! ####

# Get dendrometer counts per site from dendro_spectra
dndrmtr_by_site <- dendro_spectra %>%
  mutate(Site = substr(TreeID, 1, 2)) %>%
  group_by(Site) %>%
  summarise(n_Dndrmtr = n())

# Get BAI and age counts per site from struc_spec, then join
sample_sizes <- struc_spec %>%
  group_by(Site) %>%
  summarise(
    n_BAI_2024 = sum(!is.na(BAI_2024)),
    n_age = sum(!is.na(age))
  ) %>%
  left_join(dndrmtr_by_site, by = "Site")
sample_sizes %>%
  bind_rows(
    summarise(.,
              Site = "Total",
              n_BAI_2024 = sum(n_BAI_2024, na.rm = TRUE),
              n_age = sum(n_age, na.rm = TRUE),
              n_Dndrmtr = sum(n_Dndrmtr, na.rm = TRUE)
    )
  )

sample_sizes
