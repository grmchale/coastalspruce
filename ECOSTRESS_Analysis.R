###Use "Point" AppEEARS extraction .csv file###

###CREATE YEAR COLUMN FROM DATE - ESI###
# Read in the CSV (modify file path as needed)
ECOSTRESS_ESI <- read.csv("./4thRun-ESIL4-2018-22-ECO4ESIPTJPL-001-results.csv", stringsAsFactors = FALSE)
ECOSTRESS_ESI$Year <- substr(ECOSTRESS_ESI$Date, 1, 4)
#New column names
names(ECOSTRESS_ESI)[names(ECOSTRESS_ESI) == "ECO4ESIPTJPL_001_Evaporative_Stress_Index_PT_JPL_ESIavg"] <- "ESI"
names(ECOSTRESS_ESI)[names(ECOSTRESS_ESI) == "ECO4ESIPTJPL_001_Evaporative_Stress_Index_PT_JPL_PET"] <- "PET"
head(ECOSTRESS_ESI)

###DISPLAY CATEGORY BY YEAR - ESI###
library(dplyr)
library(ggplot2)
# 1. Summarize data to get mean ESI for each (Year, Category) pair
ECOSTRESS_ESI_summary <- ECOSTRESS_ESI %>%
  group_by(Year, Category) %>%
  summarize(median_ESI = median(ESI, na.rm = TRUE))
# 2. Create a grouped bar chart
ggplot(ECOSTRESS_ESI_summary, aes(x = Year, y = median_ESI, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  labs(title = "Median Evaporative Stress Index by Year and Forest Type - All Sites",
       x = "Year",
       y = "Average ESI") +
  theme_minimal()

###CREATE YEAR COLUMN FROM DATE - Evapotranspiration###
# Read in the CSV (modify file path as needed)
ECOSTRESS_Evap <- read.csv("./allpoints4thrun_may-sept2018-22_L3_ET_v1/4thRun-ETL3-2018-22-ECO3ETPTJPL-001-results.csv", stringsAsFactors = FALSE)
ECOSTRESS_Evap$Year <- substr(ECOSTRESS_Evap$Date, 1, 4)
#New column names + drop columns
colnames(ECOSTRESS_Evap) <- sub("^ECO3ETPTJPL_001_EVAPOTRANSPIRATION_PT_JPL_", "", colnames(ECOSTRESS_Evap))
ECOSTRESS_Evap <- ECOSTRESS_Evap[, !names(ECOSTRESS_Evap) %in% c("ETinstUncertainty_bitmask", "ETinstUncertainty_ContinuousLayer", "ETinstUncertainty_ContinuousLayer_Description")]
head(ECOSTRESS_Evap)

###DISPLAY CATEGORY BY YEAR - Evapotranspiration###
library(dplyr)
library(ggplot2)
# 1. Summarize data to get median ET for each (Year, Category) pair
ECOSTRESS_Evap_summary <- ECOSTRESS_Evap %>%
  group_by(Year, Category) %>%
  summarize(median_Evapotranspiration_daily = median(ETdaily, na.rm = TRUE))
# 2. Create a grouped bar chart
ggplot(ECOSTRESS_Evap_summary, aes(x = Year, y = median_Evapotranspiration_daily, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7)) +
  labs(title = "Median ET Daily by Year and Forest Type - All Sites",
       x = "Year",
       y = "Median ET Daily (W/m^2)") +
  theme_minimal()


