library(readr)
library(dplyr)

# --- Define file path ---
file_path <- "G:/HyperspectralUAV/R_outputs/speclib_chronologies/veg_indices/chrono_VIstats_metrics.csv"

# --- Read the CSV ---
veg_indices <- read_csv(file_path)

age_counts <- veg_indices %>%
  group_by(Site) %>%
  summarise(
    n_with_age = sum(!is.na(age)),
    n_total = n(),
    pct_with_age = round(100 * n_with_age / n_total, 1)
  ) %>%
  arrange(desc(n_with_age))

# View results
print(age_counts, n = Inf)

total_with_age <- sum(age_counts$n_with_age)
total_with_age

### Filter out S first ###
# --- Filter out rows where CC == "S" ---
veg_indices_filtered <- veg_indices %>%
  filter(CC != "S")

# --- Count rows with non-missing 'age' per Site ---
age_counts <- veg_indices_filtered %>%
  group_by(Site) %>%
  summarise(
    n_with_age = sum(!is.na(age)),
    n_total = n(),
    pct_with_age = round(100 * n_with_age / n_total, 1)
  ) %>%
  arrange(desc(n_with_age))

# --- View per-site summary ---
print(age_counts, n = Inf)

# --- Calculate total counts across all sites ---
summary_totals <- age_counts %>%
  summarise(
    total_with_age = sum(n_with_age),
    total_rows = sum(n_total),
    pct_with_age_overall = round(100 * total_with_age / total_rows, 1)
  )

# --- View overall totals ---
print(summary_totals)

library(dplyr)
library(stringr)

# --- 1. Read CSV ---
df <- read.csv("G:/HyperspectralUAV/R_outputs/speclib_amoebas_final/amoebas_VIs_and_spectra_combined.csv")

# --- 2. Remove leading "X" from wavelength columns (e.g., X398_1nm -> 398_1nm) ---
names(df) <- names(df) %>% 
  str_replace("^X([0-9]{3}_1nm)$", "\\1")

# --- 3. Remove the "_1nm" from any wavelength columns ---
names(df) <- names(df) %>% 
  str_replace("([0-9]{3})_1nm$", "\\1")

# --- 4. Identify wavelength columns (now just 398–999) ---
spec_cols <- grep("^([4-9][0-9]{2}|999)$", names(df), value = TRUE)

# --- 5. Relocate wavelength columns after sample_name ---
df2 <- df %>%
  relocate(all_of(spec_cols), .after = sample_name)

# --- 6. Rename columns ---
df2 <- df2 %>%
  rename(
    Site  = Amoeba_ID,
    Pixel = sample_name
  )

#
library(readr)

df <- read_csv("G:/Colbys_Plots/amoeba_climate_5_8_24.csv")

library(dplyr)

ppt_summary <- df %>%
  group_by(ID) %>%
  summarise(total_ppt = sum(ppt, na.rm = TRUE))

library(dplyr)
library(ggplot2)

# Define your desired ID order
site_order <- c("GI", "CE", "HI", "BI", "CC", "ET", "FP", "WP", "GW", "RI")

ppt_summary <- df %>%
  group_by(ID) %>%
  summarise(total_ppt = sum(ppt, na.rm = TRUE)) %>%
  mutate(ID = factor(ID, levels = site_order))

ggplot(ppt_summary, aes(x = ID, y = total_ppt)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Site",
    y = "Total Precipitation",
    title = "Total PPT by Site"
  ) +
  coord_cartesian(ylim = c(200, max(ppt_summary$total_ppt))) +
  theme_minimal()

######### Display a spectral profile!########
# 1. Read in the CSV ---------------------------------------------

# Use forward slashes so you don't have to escape backslashes
file_path <- "G:/HyperspectralUAV/R_outputs/speclib_chronologies/chronologies_speclib_bytree_1nm.csv"

# check.names = FALSE prevents R from adding "X" in front of numeric column names
chronos <- read.csv(file_path, check.names = FALSE)

# 2. Set up wavelengths and extract the first 10 spectra ----------

# Wavelengths 398–850 (assumes columns are named "398", "399", ..., "850")
wavelengths <- 398:850
spec_cols   <- as.character(wavelengths)

# First 10 trees (rows)
# If the file has fewer than 10 rows this will still work safely
n_to_plot <- min(10, nrow(chronos))
spec_mat  <- as.matrix(chronos[1:n_to_plot, spec_cols])

# 3. Plot setup --------------------------------------------------------------
col_blue   <- "#56B4E9"   # light blue
col_orange <- "#E69F00"   # bright orange

plot(wavelengths,
     spec_mat[1, ],
     type = "n",
     xlab = "Wavelength (nm)",
     ylab = "Reflectance",
     ylim = range(spec_mat, na.rm = TRUE),
     xlim = c(398, 850))   # explicitly stop at 800 nm

# 4. Draw thicker lines ------------------------------------------------------

# Blue (first 5)
for (i in 1:min(5, n_to_plot)) {
  lines(wavelengths, spec_mat[i, ], col = col_blue, lwd = 3)
}

# Orange (next 5)
if (n_to_plot > 5) {
  for (i in 6:n_to_plot) {
    lines(wavelengths, spec_mat[i, ], col = col_orange, lwd = 3)
  }
}

########### CREATE BAR CHART OF CORRELATIONS FOR CSHRINK ##################
# Read in the CSV -------------------------------------------------------------
file_path <- "G:/HyperspectralUAV/R_outputs/modelling/correlations/dendrometer_correlations/corr_cshrink_noMean.csv"
df <- read.csv(file_path, check.names = FALSE)

# Order by abs_shrink from largest to smallest --------------------------------
df_ordered <- df[order(-df$abs_cshrink), ]   # the minus sign sorts descending

# View the first few rows (optional)
head(df_ordered)

library(dplyr)
library(stringr)

# df_ordered already exists from previous step

df_filtered <- df_ordered %>%
  # Extract base index name BEFORE the "_<resolution>nm"
  mutate(base_index = str_extract(spectral_feature, "^[^_]+")) %>%
  
  # Group by base index name
  group_by(base_index) %>%
  
  # Keep the row with the maximum abs_cshrink within each group
  slice_max(abs_cshrink, n = 1, with_ties = FALSE) %>%
  
  ungroup()

library(dplyr)
library(ggplot2)

plot_top_indices <- function(df, n_top = 20, show_labels = TRUE,
                             bar_color = "#006744") {
  
  # Determine how many to plot
  if (is.infinite(n_top)) {
    n_plot <- nrow(df)
  } else {
    n_plot <- min(n_top, nrow(df))
  }
  
  df_plot <- df %>%
    arrange(desc(abs_cshrink)) %>%
    slice_head(n = n_plot)
  
  p <- ggplot(df_plot,
              aes(x = reorder(base_index, -abs_cshrink),
                  y = abs_cshrink)) +
    geom_col(fill = bar_color) +
    labs(x = "Index",
         y = "abs_cshrink",
         title = paste0("Top ", n_plot, " indices by |cshrink|")) +
    theme_minimal(base_size = 12)
  
  # Optionally remove index labels if cluttered
  if (!show_labels) {
    p <- p + theme(axis.text.x = element_blank())
  } else {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  # ---- Y-axis range here ----
  p <- p + coord_cartesian(ylim = c(0.4, 0.7))
  
  print(p)
}


# Examples:
plot_top_indices(df_filtered, n_top = 5, show_labels = TRUE)
plot_top_indices(df_filtered, n_top = 20, show_labels = FALSE)

###### CREATE NEW SCATTERPLOTS FOR AGU ########
# Read the RDS file -----------------------------------------------------------
file_path <- "G:/HyperspectralUAV/R_outputs/speclib_chronologies/veg_indices/chrono_VIstats_metrics.rds"

library(dplyr)
veg_indices <- readRDS(file_path)
veg_indices_clean <- veg_indices %>%
  dplyr::filter(!is.na(BAI_2024))

# Create the plot
library(dplyr)
library(ggplot2)

# USER OPTIONS 
add_lm_line <- FALSE   # TRUE = add overall linear regression line
add_r2      <- FALSE   # TRUE = display R² from that regression
point_size  <- 4      # adjust point size here
# ===

# Keep only the sites you care about, in the order you specified
sites_keep <- c("GI", "CE", "HI", "CC", "FP", "RI")

#df_plot <- veg_indices_clean %>% filter(Site %in% sites_keep) %>% mutate(Site = factor(Site, levels = sites_keep))

df_plot <- veg_indices_clean %>%
  filter(Site %in% sites_keep) %>%
  filter(
    BAI_2024 > (quantile(BAI_2024, 0.25) - 1.5 * IQR(BAI_2024)) &
      BAI_2024 < (quantile(BAI_2024, 0.75) + 1.5 * IQR(BAI_2024)),
    PRInorm_10nm_Median > (quantile(PRInorm_10nm_Median, 0.25) - 1.5 * IQR(PRInorm_10nm_Median)) &
      PRInorm_10nm_Median < (quantile(PRInorm_10nm_Median, 0.75) + 1.5 * IQR(PRInorm_10nm_Median))
  ) %>%
  mutate(Site = factor(Site, levels = sites_keep))


# Custom colors for each site
site_colors <- c(
  "GI" = "#FF7A9C",
  "CE" = "#E99136",
  "HI" = "#C3A719",
  "CC" = "#19C367",
  "FP" = "#19BFE8",
  "RI" = "#FC70DA"
)

# Base plot: points colored by Site
p <- ggplot(df_plot,
            aes(x = BAI_2024,
                y = PRInorm_10nm_Median,
                color = Site)) +
  geom_point(size = point_size) +
  scale_color_manual(values = site_colors,
                     breaks = sites_keep) +
  labs(
    x = "BAI 2024",
    y = "PRInorm (10 nm, Median)",
    color = "Site"
  ) +
  theme_minimal(base_size = 14)

# Add linear regression line (overall, not per-site), if requested
if (add_lm_line) {
  p <- p +
    geom_smooth(
      data = df_plot,
      aes(x = BAI_2024, y = PRInorm_10nm_Median),
      method = "lm",
      se = FALSE,
      inherit.aes = FALSE,
      color = "black",
      linetype = "dashed"
    )
}

# Add R² annotation, if requested
if (add_r2) {
  lm_fit <- lm(PRInorm_10nm_Median ~ BAI_2024, data = df_plot)
  r2_val <- summary(lm_fit)$r.squared
  
  p <- p +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = paste0("R² = ", round(r2_val, 3)),
      hjust = 1.1, vjust = 1.5,
      size = 5
    )
}

# Draw plot
print(p)









