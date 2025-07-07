library(dplyr)

cum_zg_suma<-read.csv("G:/Dendrometers/cum_zg_sums.csv")
drone_dendro<-read.csv("G:/Dendrometers/drone_dendro.csv")
uptoflight_zg<-read.csv("G:/Dendrometers/cumulative_zg_uptoflight.csv")

cum_zg_suma <- cum_zg_suma %>%
  left_join(drone_dendro %>% select(Dndrmtr, TreeID), 
            by = c("id" = "Dndrmtr"))

cum_zg_suma <- cum_zg_suma %>%
  left_join(uptoflight_zg %>% select(Dndrmtr, zg_fraction), 
            by = c("id" = "Dndrmtr"))

cum_zg_suma <- cum_zg_suma %>%
  mutate(
    zg_fraction_7 = zg_7 / 7,
    zg_fraction_15 = zg_15 / 15,
    zg_fraction_30 = zg_30 / 30,
    zg_fraction_60 = zg_60 / 60,
    zg_fraction_90 = zg_90 / 90
  )

write.csv(cum_zg_suma,"G:/Dendrometers/uptoflight_intervals_zg.csv")

###################################################################################
## TEST SPECTRA VS ZG FRACTION-------------------------------------------------
# Boochs vs. zg_fraction_7
library(ggplot2)


# Fit linear model first
fit <- lm(zg_fraction_7 ~ Boochs_1nm_Median, data = df_rf)
rsquared <- summary(fit)$r.squared

# Create plot
ggplot(df_rf, aes(x = Boochs_1nm_Median, y = zg_fraction_7)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Scatterplot of Boochs_1nm_Median vs. zg_fraction_7",
    subtitle = paste("R² =", round(rsquared, 3)),
    x = "Boochs_1nm_Median",
    y = "zg_fraction_7"
  ) +
  theme_minimal()

# PRI_15nm_Median vs. zg_fraction_15

# Fit linear model
fit <- lm(zg_fraction_15 ~ PRI_15nm_Median, data = df_rf)
rsquared <- summary(fit)$r.squared

# Create scatterplot with best-fit line
ggplot(df_rf, aes(x = PRI_15nm_Median, y = zg_fraction_15)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Scatterplot of PRI_15nm_Median vs. zg_fraction_15",
    subtitle = paste("R² =", round(rsquared, 3)),
    x = "PRI_15nm_Median",
    y = "zg_fraction_15"
  ) +
  theme_minimal()

# PRI_15nm_Median vs. zg_fraction_30

# Fit linear model
fit <- lm(zg_fraction_30 ~ PRI_15nm_Median, data = df_rf)
rsquared <- summary(fit)$r.squared

# Create scatterplot with best-fit line
ggplot(df_rf, aes(x = PRI_15nm_Median, y = zg_fraction_30)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Scatterplot of PRI_15nm_Median vs. zg_fraction_30",
    subtitle = paste("R² =", round(rsquared, 3)),
    x = "PRI_15nm_Median",
    y = "zg_fraction_30"
  ) +
  theme_minimal()

