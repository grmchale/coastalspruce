###### TESTING SUB DAILY TREE WATER DEFECIT (TWD) ######
## Uses sample trees from CC + RI

# Read in the CSV files from subdaily_TWD folder
df_99 <- read.csv("data/subdaily_TWD/data_92243799_2024_10_18_0.csv", sep = ";")
df_05 <- read.csv("data/subdaily_TWD/data_92243705_2024_10_25_0.csv",  sep = ";")

# Drop the last 3 columns from each df
df_99 <- df_99[, 1:(ncol(df_99) - 3)]
df_05 <- df_05[, 1:(ncol(df_05) - 3)]

### Work with df_05 (Roque Island, RI) ###
# Filter rows where the second column starts with "2024.08.30"
df_05_drone <- df_05[grepl("^2024\\.08\\.30", df_05[, 2]), ]
# Extract the time portion from the second column into a new "time" column
df_05$time <- sub(".*\\s", "", df_05[, 2])

# Create plot
library(ggplot2)

ggplot(df_05_drone, aes(x = time, y = X24.25)) +
  geom_line(group = 1) +
  labs(x = "Time", y = "X24.25") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(ggplot2)
library(dplyr)

# Filter to the date range and convert datetime
df_05_summer <- df_05 %>%
  filter(as.POSIXct(.[, 2], format = "%Y.%m.%d %H:%M") >= as.POSIXct("2024.05.18", format = "%Y.%m.%d") &
           as.POSIXct(.[, 2], format = "%Y.%m.%d %H:%M") <= as.POSIXct("2024.10.25", format = "%Y.%m.%d")) %>%
  mutate(datetime = as.POSIXct(.[, 2], format = "%Y.%m.%d %H:%M"))

# Plot
ggplot(df_05_summer, aes(x = datetime, y = X1)) +
  geom_line() +
  scale_x_datetime(date_breaks = "2 weeks", date_labels = "%b %d") +
  labs(x = "Date", y = "X24.25") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
