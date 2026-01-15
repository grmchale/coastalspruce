#####################  A M O E B A S ############################
#############ANOVA + PERMANOVA + PCA + NMDS ANALYSIS ############
############################# PREP DATAFRAME ########################
VIs_1 <-readRDS("./R_outputs/speclib_amoebas_final/amoebas_VIs_and_spectra_combined.rds")

library(dplyr)
library(readr)
library(readxl)

# --- 1. Read in base dataframe ---
# Assuming VIs_1 already exists in your environment

# --- 2. Read and prepare supporting dataframes ---

# amoeba_latlon_topo.csv → keep all except ID
latlon <- read_csv("G:/LiDAR/Normalized/LiDAR_Amoebas/amoeba_latlon_topo.csv") %>%
  rename(Amoeba_ID = ID)
# amoeba_metrics.csv → keep all except ID
metrics <- read_csv("G:/LiDAR/Normalized/LiDAR_Amoebas/amoeba_metrics.csv") %>%
  rename(Amoeba_ID = ID)
# amoeba_rugosity.csv → keep all except ID
rugosity <- read_csv("G:/LiDAR/Normalized/LiDAR_Amoebas/amoeba_rugosity.csv") %>%
  rename(Amoeba_ID = ID)
# amoeba_dtc.csv → keep only Dist_to_Coast
dtc <- read_csv("G:/LiDAR/Normalized/LiDAR_Amoebas/amoeba_dtc.csv") %>%
  select(ID, Dist_to_Coast) %>%
  rename(Amoeba_ID = ID)
# full_site_metrics.xlsx → keep only median, depth, cindex
site_metrics <- read_excel("G:/Colbys_Plots/full_site_metrics.xlsx") %>%
  select(Site, median, depth, cindex)

# --- 3. Join everything together ---

VIs_1_joined <- VIs_1 %>%
  left_join(latlon,    by = "Amoeba_ID") %>%
  left_join(metrics,   by = "Amoeba_ID") %>%
  left_join(rugosity,  by = "Amoeba_ID") %>%
  left_join(dtc,       by = "Amoeba_ID") %>%
  left_join(site_metrics, by = c("Amoeba_ID" = "Site"))

library(lubridate)
# --- 1. Rename fields ---
VIs_1_joined <- VIs_1_joined %>%
  rename(
    median_age = median,
    soil_depth = depth
  )

# --- 2. Create DOY mapping for each Amoeba_ID ---
# First, define a named vector of dates
amoeba_dates <- c(
  "CC" = "2024-08-12",
  "ET" = "2024-08-12",
  "BI" = "2024-08-12",
  "HI" = "2024-08-13",
  "CE" = "2024-08-14",
  "GI" = "2024-08-14",
  "FP" = "2024-08-29",
  "WP" = "2024-08-29",
  "GW" = "2024-08-29",
  "RI" = "2024-08-30"
)
# Map the date based on Amoeba_ID
VIs_1_joined <- VIs_1_joined %>%
  mutate(
    Date = as.Date(amoeba_dates[Amoeba_ID]),
    DOY = yday(Date)
  ) %>%
  select(-Date)  # Optionally remove the intermediate Date column

library(readr)
climate_5_8_24 <- read_csv("G:/Colbys_Plots/amoeba_climate_5_8_24.csv")

# Consolidate climate data by site
climate_summary <- climate_5_8_24 %>%
  group_by(ID) %>%
  summarise(
    ppt = sum(ppt),           # Total growing season precipitation
    tmin = mean(tmin),        # Mean minimum temperature
    tmean = mean(tmean),      # Mean temperature
    tmax = mean(tmax),        # Mean maximum temperature
    tdmean = mean(tdmean),    # Mean dewpoint temperature
    vpdmin = mean(vpdmin),    # Mean minimum VPD
    vpdmax = mean(vpdmax)     # Mean maximum VPD
  )
VIs_1_joined_2 <- VIs_1_joined %>%
  left_join(climate_summary, by = c("Amoeba_ID" = "ID"))
# Save both formats
saveRDS(VIs_1_joined_2, "./R_outputs/speclib_amoebas_final/amoebas_VIs_spectra_varbs.rds")
write.csv(VIs_1_joined_2, "./R_outputs/speclib_amoebas_final/amoebas_VIs_spectra_varbs.csv", row.names = FALSE)
######################## READ IN DATAFRAME(S)! ###########################
setwd("G:/HyperspectralUAV/")
# Load necessary packages
library(dplyr)
library(vegan)
library(ggplot2)

# Read .rds
VIs_df <- readRDS("./R_outputs/speclib_amoebas_final/amoebas_VIs_spectra_varbs.rds")

# Create random subsample
set.seed(42)  # for reproducibility
VI_subset <- VIs_df %>%
  group_by(Amoeba_ID) %>%     
  sample_n(size = 400) %>%  # 400 pixels per site
  ungroup()
table(VI_subset$Amoeba_ID)
nrow(VI_subset)

##################### PCA WITH ORDINATION !! ############################
library(vegan)
library(ggplot2)

#### PCA with spectral columns (1nm library), ALL 4,000 ROWS ####
start_col <- which(names(VI_subset) == "398_1nm")
end_col   <- which(names(VI_subset) == "850_1nm")

spec_mat <- VI_subset[, start_col:end_col]

# Standardize (mean=0, sd=1)
spec_scaled <- scale(spec_mat)

# PCA using vegan (for later use in envfit)
pca_res <- rda(spec_scaled)  # same as PCA (prcomp) but integrates with envfit

# Check variance explained
summary(pca_res)

# Extract loadings from RDA (same as PCA)
pca_loadings <- as.data.frame(scores(pca_res, display = "species"))
pca_loadings$Wavelength <- as.numeric(gsub("_1nm", "", rownames(pca_loadings)))

# Plot PC1 and PC2 loadings
ggplot(pca_loadings, aes(x = Wavelength)) +
  geom_line(aes(y = PC1), color = "#1F78B4", size = 0.8) +
  geom_line(aes(y = PC2), color = "#E31A1C", size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCA Loadings by Wavelength (PC1 & PC2)",
    x = "Wavelength (nm)",
    y = "Loading Weight"
  ) +
  scale_x_continuous(breaks = seq(400, 850, 50))

# Test environmental variables
# 1. Define environmental variables
env_vars <- VI_subset %>%
  dplyr::select(
    Area_m2, Lat, Long, Elevation, Slope_deg, Aspect_deg,
    zmean, zsd, zq95, vertical_rugosity, horizontal_rugosity,
    Dist_to_Coast, median_age, soil_depth, cindex, DOY,
    ppt, tmin, tmean, tmax, tdmean, vpdmin, vpdmax
  )
# 2. Fit environmental vectors onto PCA ordination
set.seed(123)  # reproducibility for permutation test
env_fit <- envfit(pca_res, env_vars, permutations = 999)
print(env_fit)
# 3. Plot ordination with fitted environmental vectors
plot(pca_res, display = "sites", type = "points", col = "gray70", main = "PCA + Environmental Vectors")
plot(env_fit, p.max = 0.05, col = "red")  # show only significant arrows (p < 0.05)
# 4. Extract envfit results into a tidy table
env_arrows <- as.data.frame(env_fit$vectors$arrows)
env_table <- env_arrows %>%
  mutate(
    Variable = rownames(env_arrows),
    r2 = env_fit$vectors$r,
    p  = env_fit$vectors$pvals
  ) %>%
  rename(PC1 = 1, PC2 = 2) %>%  # rename the first two columns to PC1/PC2
  select(Variable, r2, p, PC1, PC2) %>%
  arrange(desc(r2))
# 5. View top environmental gradients
print(env_table[1:10, ], row.names = FALSE)  # top 10 strongest relationships
library(vegan)
# 6) Colors for sites (auto palette; tweak if you want)
sites <- factor(VI_subset$Amoeba_ID)
site_cols <- c(
  "BI" = "#8FD17F",
  "CC" = "#4E79A7",  # muted blue
  "CE" = "#F28E2B",  # warm orange
  "ET" = "#E15759",  # soft red
  "FP" = "#76B7B2",  # turquoise
  "GI" = "#387C44",  # green
  "GW" = "#EDC948",  # yellow
  "HI" = "#B07AA1",  # violet
  "RI" = "#FF9DA7",  # light pink
  "WP" = "#9C755F"   # brown-grey
)
# 7) PCA site scores
sc_sites <- scores(pca_res, display = "sites")
# 8) Base PCA frame + colored points
plot(pca_res, display = "sites", type = "n",
     main = "PCA with Top 6 Environmental Variables",
     xlab = "PC1", ylab = "PC2",
     xlim = c(-3, 2))  # <-- manually set PC1 axis range
abline(h = 0, v = 0, lty = 3, col = "grey75")
points(sc_sites[,1], sc_sites[,2],
       pch = 16, cex = 0.7, col = site_cols[as.character(sites)])  # ← smaller points
legend("topleft", legend = levels(sites), pch = 16, pt.cex = 1,
       col = site_cols, bty = "n")

# 9) Top 10 env vectors (by r²)
r2 <- env_fit$vectors$r
keep <- names(sort(r2, decreasing = TRUE))[seq_len(min(6, length(r2)))]
vec_all <- scores(env_fit, display = "vectors")
vec10   <- vec_all[keep, , drop = FALSE]
# scale arrows like plot(envfit)
mul <- vegan:::ordiArrowMul(vec10, sc_sites)

# 10) Black arrows + labels for top 10
mul <- vegan:::ordiArrowMul(vec10, sc_sites) * 9  # try 3–5× to taste

# --- Replot arrows and labels ---
arrows(0, 0, vec10[,1]*mul, vec10[,2]*mul,
       length = 0.08, col = "black", lwd = 2.3)
text(vec10[,1]*mul*1.08, vec10[,2]*mul*1.08,
     labels = rownames(vec10), cex = 1.0, col = "black")

# 11) Export results
# Create output directory if it doesn't exist
out_dir <- "./R_outputs/pca/amoebas_pca"
#Extract PCA results
pca_scores   <- as.data.frame(scores(pca_res, display = "sites"))    # site (sample) scores
pca_loadings <- as.data.frame(scores(pca_res, display = "species"))  # variable loadings (wavelengths)
pca_eig      <- eigenvals(pca_res)                                  # eigenvalues
pca_var      <- data.frame(
  PC = paste0("PC", seq_along(pca_eig)),
  Eigenvalue = pca_eig,
  Proportion = pca_eig / sum(pca_eig),
  Cumulative = cumsum(pca_eig / sum(pca_eig))
)
# Write CSVs
write.csv(pca_scores,   file.path(out_dir, "pca_scores_rda.csv"),   row.names = TRUE)
write.csv(pca_loadings, file.path(out_dir, "pca_loadings_rda.csv"), row.names = TRUE)
write.csv(pca_var,      file.path(out_dir, "pca_variance_rda.csv"), row.names = FALSE)
# Save RDS objects for easy reload
saveRDS(pca_res,        file.path(out_dir, "pca_res_rda.rds"))
saveRDS(pca_scores,     file.path(out_dir, "pca_scores_rda.rds"))
saveRDS(pca_loadings,   file.path(out_dir, "pca_loadings_rda.rds"))
saveRDS(pca_var,        file.path(out_dir, "pca_variance_rda.rds"))
# Save envfit table
write.csv(env_table, file.path(out_dir, "envfit_table.csv"), row.names = FALSE)
saveRDS(env_table,    file.path(out_dir, "envfit_table.rds"))

### PCA with spectral columns (1nm library), COLLAPSED INTO SITES N=10 ####
library(dplyr)
library(vegan)

# --- 1) Identify spectral range (398–999 nm) ---
start_col <- which(names(VI_subset) == "398_1nm")
end_col   <- which(names(VI_subset) == "850_1nm")

stopifnot(length(start_col) == 1L, length(end_col) == 1L, end_col > start_col)

spec_cols <- names(VI_subset)[start_col:end_col]

# --- 2) Collapse to site medians (n = 10 rows) ---
site_spec_med <- VI_subset %>%
  group_by(Amoeba_ID) %>%
  summarize(across(all_of(spec_cols), ~median(.x, na.rm = TRUE)), .groups = "drop")

# --- 3) One row of environmental variables per site ---
env_vars_site <- VI_subset %>%
  group_by(Amoeba_ID) %>%
  summarize(
    Area_m2 = first(Area_m2), Lat = first(Lat), Long = first(Long),
    Elevation = first(Elevation), Slope_deg = first(Slope_deg), Aspect_deg = first(Aspect_deg),
    zmean = first(zmean), zsd = first(zsd), zq95 = first(zq95),
    vertical_rugosity = first(vertical_rugosity), horizontal_rugosity = first(horizontal_rugosity),
    Dist_to_Coast = first(Dist_to_Coast), median_age = first(median_age),
    soil_depth = first(soil_depth), cindex = first(cindex), DOY = first(DOY),
    ppt = first(ppt), tmin = first(tmin), tmean = first(tmean), tmax = first(tmax),
    tdmean = first(tdmean), vpdmin = first(vpdmin), vpdmax = first(vpdmax),
    .groups = "drop"
  )

# --- 4) Matrix for PCA (scale = shape emphasis) ---
spec_mat_site   <- as.matrix(site_spec_med[, spec_cols, drop = FALSE])
spec_scaled_site <- scale(spec_mat_site)

# --- 5) PCA via vegan::rda on site medians ---
pca_site <- rda(spec_scaled_site)
summary(pca_site)  # variance explained etc.

# --- 5a) View the weightings
# Extract variable (wavelength) loadings
pca_loadings_site <- as.data.frame(scores(pca_site, display = "species"))

# Add wavelength column (numeric)
pca_loadings_site$Wavelength <- as.numeric(gsub("_1nm", "", rownames(pca_loadings_site)))

# Plot PC1 & PC2 loadings vs wavelength
ggplot(pca_loadings_site, aes(x = Wavelength)) +
  geom_line(aes(y = PC1, color = "PC1"), size = 0.9) +
  geom_line(aes(y = PC2, color = "PC2"), size = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  scale_color_manual(values = c("PC1" = "#1F78B4", "PC2" = "#E31A1C")) +
  theme_minimal(base_size = 14) +
  labs(title = "PCA Loadings by Wavelength (Site-Median Spectra)",
       x = "Wavelength (nm)",
       y = "Loading Weight",
       color = "Principal Component") +
  scale_x_continuous(breaks = seq(400, 1000, 50))


# --- 6) Envfit at true site scale (valid p-values) ---
env_fit_site <- envfit(pca_site, env_vars_site[, -match("Amoeba_ID", names(env_vars_site))], permutations = 999)
print(env_fit_site)

# --- 7) Ranked envfit table (site level) ---
env_arrows_site <- as.data.frame(env_fit_site$vectors$arrows)
env_table_site <- env_arrows_site %>%
  mutate(
    Variable = rownames(env_arrows_site),
    r2 = env_fit_site$vectors$r,
    p  = env_fit_site$vectors$pvals
  ) %>%
  rename(PC1 = 1, PC2 = 2) %>%
  select(Variable, r2, p, PC1, PC2) %>%
  arrange(desc(r2))

print(env_table_site, row.names = FALSE)

# --- 8) Base plot: sites + black arrows-
# Ensure site color map exists (reuses your earlier palette)
stopifnot(exists("site_cols"))
sites_all <- factor(VI_subset$Amoeba_ID)  # keeps the same levels/order as your big PCA
site_levels_fixed <- levels(sites_all)

# Scores
sc_sites_site <- scores(pca_site, display = "sites")

# Base frame + points (colored by site with your palette)
plot(pca_site, display = "sites", type = "n",
     main = "PCA, Site Median Reflectance (398–850 nm) + Environmental Variables",
     xlab = "PC1", ylab = "PC2")
abline(h = 0, v = 0, lty = 3, col = "grey75")

# Map colors using your palette; preserve factor levels
sites_site <- factor(site_spec_med$Amoeba_ID, levels = site_levels_fixed)
pt_cols_site <- site_cols[as.character(sites_site)]

# Plot points
points(sc_sites_site[,1], sc_sites_site[,2],
       pch = 16, cex = 2.1, col = pt_cols_site)

# Add site labels next to each point
text(sc_sites_site[,1], sc_sites_site[,2],
     labels = sites_site, pos = 3,  # 3 = above point; try 4 for right, 1 for below
     cex = 0.9, font = 2, col = pt_cols_site)

# Top env vectors by r² (to match earlier)
r2s   <- env_fit_site$vectors$r
keep6 <- names(sort(r2s, decreasing = TRUE))[seq_len(min(2, length(r2s)))] # Change to 6?
vec   <- scores(env_fit_site, display = "vectors")[keep6, , drop = FALSE]

# Scale arrows like before (same multiplier and thickness)
mul <- vegan:::ordiArrowMul(vec, sc_sites_site) * 9

# Draw arrows
arrows(0, 0, vec[,1]*mul, vec[,2]*mul,
       length = 0.08, col = "black", lwd = 2.3)

# --- label positions (start from your usual 1.08x)
lab_scale <- 1.08
lab_x <- vec[,1]*mul*lab_scale
lab_y <- vec[,2]*mul*lab_scale

# axis-range–aware offsets (robust to plot scale)
scores_mat <- if (exists("sc_sites_site")) sc_sites_site else sc_sites
dx <- 0.06 * diff(range(scores_mat[,1]))  # ~5% of x-range
dy <- 0.05 * diff(range(scores_mat[,2]))  # ~5% of y-range

# nudge specific labels
if ("Dist_to_Coast" %in% rownames(vec)) {
  lab_x["Dist_to_Coast"] <- lab_x["Dist_to_Coast"] - dx     # left
}
if ("Elevation" %in% rownames(vec)) {
  lab_x["Elevation"] <- lab_x["Elevation"] + 0.08 * diff(range(scores_mat[,1]))  # right
  lab_y["Elevation"] <- lab_y["Elevation"] + 0.01 * diff(range(scores_mat[,2]))  # down
}

# draw labels
text(lab_x, lab_y, labels = rownames(vec), cex = 1.0, col = "black")

# --- 9) (Optional) export site-level results ---
out_dir <- "./R_outputs/pca/amoebas_pca_site"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(site_spec_med, file.path(out_dir, "site_median_spectra.csv"), row.names = FALSE)
write.csv(env_vars_site, file.path(out_dir, "site_env.csv"), row.names = FALSE)
write.csv(as.data.frame(scores(pca_site, display = "sites")), file.path(out_dir, "pca_scores_site.csv"))
write.csv(as.data.frame(scores(pca_site, display = "species")), file.path(out_dir, "pca_loadings_site.csv"))
write.csv(env_table_site, file.path(out_dir, "envfit_table_site.csv"), row.names = FALSE)
saveRDS(pca_site, file.path(out_dir, "pca_site.rds")) 
saveRDS(env_fit_site, file.path(out_dir, "env_fit_site.rds"))


################# ANOVA WITH SELECT INDICES #############################
# PRI_15nm,REPLE_5nm,NDVI_5nm,Vogelmann3_5nm,PRInorm_5nm,Datt3_5nm
library(dplyr)
library(rstatix)

#### PRInorm Welch ANOVA ####
# Test normality
# --- 1. Fit a one-way ANOVA model ---
anova_model <- aov(PRInorm_5nm ~ Amoeba_ID, data = VI_subset)

# --- 2. Check residual normality ---
# Shapiro-Wilk test (sensitive for large samples)
shapiro.test(residuals(anova_model))

# Visual check
par(mfrow = c(1,2))
hist(residuals(anova_model), main = "Residual Histogram", xlab = "Residuals")
qqnorm(residuals(anova_model)); qqline(residuals(anova_model), col = "red")
par(mfrow = c(1,1))

# Bartlett's test (requires approximate normality)
bartlett.test(PRInorm_5nm ~ Amoeba_ID, data = VI_subset)
# Brown-Forsythe version of Levene's test (no car package)
oneway.test(PRInorm_5nm ~ Amoeba_ID, data = VI_subset, var.equal = FALSE)

# Welch ANOVA
welch_result <- oneway.test(PRInorm_5nm ~ Amoeba_ID, data = VI_subset)
print(welch_result)

# Games–Howell post-hoc
gh_result <- VI_subset %>%
  games_howell_test(PRInorm_5nm ~ Amoeba_ID)

# View results
gh_result
pairwise.t.test(VI_subset$PRInorm_5nm, VI_subset$Amoeba_ID,
                p.adjust.method = "bonferroni", pool.sd = FALSE)

library(ggplot2)

means_ci <- VI_subset %>%
  group_by(Amoeba_ID) %>%
  summarize(
    n = n(),
    mean = mean(PRInorm_5nm, na.rm = TRUE),
    se = sd(PRInorm_5nm, na.rm = TRUE)/sqrt(n),
    ci_low  = mean + qt(0.025, df = n - 1) * se,
    ci_high = mean + qt(0.975, df = n - 1) * se,
    .groups = "drop"
  )
ggplot(means_ci, aes(x = reorder(Amoeba_ID, mean), y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
  coord_flip() +
  labs(title = "PRInorm_5nm by Site (mean ± 95% CI)",
       x = "Site", y = "PRInorm_5nm") +
  theme_minimal(base_size = 12)

# Box and whisker plot


# Export results
out_dir <- "./R_outputs/modelling/welchANOVA/PRInorm"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Welch ANOVA 
saveRDS(welch_result, file = file.path(out_dir, "PRInorm_welch_result.rds"))
capture.output(print(welch_result),
               file = file.path(out_dir, "PRInorm_welch_result.txt"))
# Games–Howell post-hoc
gh_path_csv <- file.path(out_dir, "PRInorm_games_howell_results.csv")
gh_path_rds <- file.path(out_dir, "PRInorm_games_howell_results.rds")
write.csv(gh_result, gh_path_csv, row.names = FALSE)
saveRDS(gh_result, gh_path_rds)

# Bonferroni pairwise t-tests
bonf_result <- pairwise.t.test(VI_subset$PRInorm_5nm, VI_subset$Amoeba_ID,
                               p.adjust.method = "bonferroni",
                               pool.sd = FALSE)
write.csv(bonf_result$p.value,
          file = file.path(out_dir, "PRInorm_pairwise_bonferroni_pmatrix.csv"),
          row.names = TRUE)

# Keep a copy of VI_subset site means too (helps with plotting/interpretation)
site_summary <- VI_subset %>%
  group_by(Amoeba_ID) %>%
  summarize(n = n(),
            mean = mean(PRInorm_5nm, na.rm = TRUE),
            sd   = sd(PRInorm_5nm, na.rm = TRUE),
            .groups = "drop")
write.csv(site_summary, file.path(out_dir, "PRInorm_site_summary.csv"), row.names = FALSE)

### Datt3 Welch ANOVA #####
# Test normality
# --- 1. Fit a one-way ANOVA model ---
anova_model_datt <- aov(Datt3_5nm ~ Amoeba_ID, data = VI_subset)

# --- 2. Check residual normality ---
# Shapiro-Wilk test (sensitive for large samples)
shapiro.test(residuals(anova_model_datt))

# Visual check
par(mfrow = c(1,2))
hist(residuals(anova_model_datt), main = "Residual Histogram", xlab = "Residuals")
qqnorm(residuals(anova_model_datt)); qqline(residuals(anova_model_datt), col = "red")
par(mfrow = c(1,1))

# Bartlett's test (requires approximate normality)
bartlett.test(Datt3_5nm ~ Amoeba_ID, data = VI_subset)
# Brown-Forsythe version of Levene's test (no car package)
oneway.test(Datt3_5nm ~ Amoeba_ID, data = VI_subset, var.equal = FALSE)

# Welch ANOVA
welch_result_datt <- oneway.test(Datt3_5nm ~ Amoeba_ID, data = VI_subset)
print(welch_result_datt)

# Games–Howell post-hoc
gh_result_datt <- VI_subset %>%
  games_howell_test(Datt3_5nm ~ Amoeba_ID)

# View results
gh_result_datt
pairwise.t.test(VI_subset$Datt3_5nm, VI_subset$Amoeba_ID,
                p.adjust.method = "bonferroni", pool.sd = FALSE)

library(ggplot2)

means_ci_datt <- VI_subset %>%
  group_by(Amoeba_ID) %>%
  summarize(
    n = n(),
    mean = mean(Datt3_5nm, na.rm = TRUE),
    se = sd(Datt3_5nm, na.rm = TRUE)/sqrt(n),
    ci_low  = mean + qt(0.025, df = n - 1) * se,
    ci_high = mean + qt(0.975, df = n - 1) * se,
    .groups = "drop"
  )
ggplot(means_ci_datt, aes(x = reorder(Amoeba_ID, mean), y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
  coord_flip() +
  labs(title = "Datt3_5nm by Site (mean ± 95% CI)",
       x = "Site", y = "Datt3_5nm") +
  theme_minimal(base_size = 12)

# Export results
out_dir <- "./R_outputs/modelling/welchANOVA/Datt3"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Welch ANOVA 
saveRDS(welch_result_datt, file = file.path(out_dir, "Datt3_welch_result.rds"))
capture.output(print(welch_result_datt),
               file = file.path(out_dir, "Datt3_welch_result.txt"))
# Games–Howell post-hoc
gh_path_csv <- file.path(out_dir, "Datt3_games_howell_results.csv")
gh_path_rds <- file.path(out_dir, "Datt3_games_howell_results.rds")
write.csv(gh_result_datt, gh_path_csv, row.names = FALSE)
saveRDS(gh_result_datt, gh_path_rds)

# Bonferroni pairwise t-tests
bonf_result_datt <- pairwise.t.test(VI_subset$Datt3_5nm, VI_subset$Amoeba_ID,
                               p.adjust.method = "bonferroni",
                               pool.sd = FALSE)
write.csv(bonf_result_datt$p.value,
          file = file.path(out_dir, "Datt3_pairwise_bonferroni_pmatrix.csv"),
          row.names = TRUE)

# Keep a copy of VI_subset site means too (helps with plotting/interpretation)
site_summary_datt <- VI_subset %>%
  group_by(Amoeba_ID) %>%
  summarize(n = n(),
            mean = mean(Datt3_5nm, na.rm = TRUE),
            sd   = sd(Datt3_5nm, na.rm = TRUE),
            .groups = "drop")
write.csv(site_summary_datt, file.path(out_dir, "Datt3_site_summary.csv"), row.names = FALSE)

####
# 🌿 Customizable Boxplot for Vegetation Indices
####
library(ggplot2)
library(scales)   # for hue_pal()

# On-Ramps (edit these) 
index_col  <- "Datt3_5nm"                     # ← choose your index column
plot_title <- "Datt3 Index Across Sites"
x_title    <- "Site"
y_title    <- "Datt3 (5nm)"
save_path  <- "./R_outputs/figures/amoebas/Datt3_5nm_boxplot.png"

# Site order 
site_order <- c("GI", "CE", "HI", "BI", "CC", "ET", "FP", "WP", "GW", "RI")
VI_subset$Amoeba_ID <- factor(VI_subset$Amoeba_ID, levels = site_order)

# Auto-generate rainbow palette 
n_sites <- length(site_order)
site_cols <- hue_pal(h = c(0, 360), c = 100, l = 65)(n_sites)

# Plot 
p_box <- ggplot(VI_subset, aes_string(x = "Amoeba_ID", y = index_col, fill = "Amoeba_ID")) +
  geom_boxplot(width = 0.7, outlier.size = 0.6, alpha = 0.9) +
  scale_fill_manual(values = site_cols) +
  labs(
    title = plot_title,
    x = x_title,
    y = y_title
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# Display and Save
print(p_box)
ggsave(save_path, p_box, width = 8, height = 5, dpi = 300, bg = "white")


####################### ADONIS perMANOVA (October, 2025) ############################
# Randomly sample 200 rows per site
# Legacy code, we used 400 from the PCA above!!! 
set.seed(123) # for reproducibility
VIs_df_sampled <- VIs_df %>%
  group_by(Amoeba_ID) %>%
  slice_sample(n = 200) %>% 
  ungroup()

library(vegan)
# First, install remotes if you don't have it
#install.packages("remotes")

# Then install the GitHub version
#remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

# Load it
library(pairwiseAdonis)

# 1. Extract the index matrix (columns ARI1_5nm through WBI_5nm)
indices_matrix <- VI_subset %>% # Use VI_subset from ANOVA code above
  dplyr::select(ARI1_5nm:WBI_5nm)

# 2. Standardize (z-score) each index to mean = 0, sd = 1
#    Bray-Curtis is sensitive to absolute magnitudes but fine with negative values.
indices_scaled <- scale(indices_matrix)

# 3. Compute Bray–Curtis dissimilarity
D <- dist(indices_scaled, method = "euclidean")

# 4. Run global PERMANOVA
set.seed(123)
adonis_res <- adonis2(D ~ Amoeba_ID, data = VI_subset, permutations = 999)
cat("\n### Global PERMANOVA Results ###\n")
print(adonis_res)

# 5. Pairwise PERMANOVA with FDR correction
pairwise_res <- pairwiseAdonis::pairwise.adonis(
  x = indices_scaled,
  factors = VI_subset$Amoeba_ID,
  sim.method = "euclidean",
  p.adjust.m = "BH",
  perm = 999
)
cat("\n### Pairwise PERMANOVA (FDR corrected) ###\n")
print(pairwise_res)

# 6. Check multivariate dispersion (homogeneity of spread)
bd <- betadisper(D, VI_subset$Amoeba_ID)
cat("\n### Beta-dispersion ANOVA ###\n")
print(anova(bd))
cat("\n### Tukey HSD pairwise dispersion ###\n")
print(TukeyHSD(bd))

#### Save PERMANOVA ouputs! ####
out_dir <- file.path("./R_outputs/modelling/perMANOVA/amoeba_pm")
write.csv(pairwise_res, file.path(out_dir, "pairwise_PERMANOVA.csv"),
          row.names = FALSE)
saveRDS(pairwise_res, file.path(out_dir, "pairwise_PERMANOVA.rds"))

# 1) Global PERMANOVA (adonis2)
sink(file.path(out_dir, "global_PERMANOVA.txt"))
cat("### Global PERMANOVA (adonis2) ###\n\n")
print(adonis_res)
sink()
saveRDS(adonis_res, file.path(out_dir, "global_PERMANOVA.rds"))

# 3) Beta-dispersion (betadisper)
bd_anova <- anova(bd)
bd_tukey <- TukeyHSD(bd)

sink(file.path(out_dir, "betadisper_ANOVA.txt"))
cat("### Beta-dispersion ANOVA ###\n\n")
print(bd_anova)
sink()

# TukeyHSD output (list) → CSV for the $group table if present + RDS for full
if (is.list(bd_tukey) && !is.null(bd_tukey$group)) {
  tukey_df <- as.data.frame(bd_tukey$group)
  tukey_df$contrast <- rownames(tukey_df)
  tukey_df <- tukey_df[, c("contrast","diff","lwr","upr","p adj")]
  utils::write.csv(tukey_df,
                   file.path(out_dir, "betadisper_Tukey_group.csv"),
                   row.names = FALSE)
}
saveRDS(bd_tukey, file.path(out_dir, "betadisper_Tukey.rds"))

# 4) Useful metadata for reproducibility
sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()
writeLines(c(
  "Indices used (columns):",
  paste(colnames(indices_matrix), collapse = ", ")
), con = file.path(out_dir, "indices_columns.txt"))

cat("\nWrote outputs to:\n", normalizePath(out_dir), "\n")


#Ordination + drivers: NMDS (Euclidean) + envfit -
#set.seed(123)
#nmds <- metaMDS(indices_matrix, distance = "euclidean", k = 2, trymax = 50, autotransform = FALSE)
#cat("\nNMDS stress:\n"); print(nmds$stress)

# Fit VI vectors (drivers) onto NMDS
#fit <- envfit(nmds, indices_matrix, permutations = 999)
#cat("\nEnvfit results (r^2, p):\n"); print(fit)

# (Minimal base plot; swap to ggplotify later if desired)
# Option A: reset margins
#par(mar = c(4, 4, 2, 2))  

#plot(nmds, type = "n")  # sets up blank ordination space
#points(scores(nmds, display = "sites"), 
      # col = as.integer(VIs_df_sampled$Amoeba_ID), 
      # pch = 16, cex = 0.5)
#ordiellipse(nmds, VIs_df_sampled$Amoeba_ID, draw = "polygon", col = 1:10, alpha = 0.2, label = TRUE)
#plot(fit, p.max = 0.05, col = "black")  # only show vectors with p ≤ 0.05

# Extract site scores
#pts <- scores(nmds, display = "sites")
# Set up the blank NMDS plot with limited y-range
#plot(pts[,1], pts[,2],
##    xlab = "NMDS1", ylab = "NMDS2",
#     xlim = range(pts[,1]), ylim = c(-1, 1),
#     type = "n")

# Add points
#points(pts, col = as.integer(VIs_df_sampled$Amoeba_ID),
#       pch = 16, cex = 0.5)

# Add site ellipses
#ordiellipse(nmds, VIs_df_sampled$Amoeba_ID,
#            draw = "polygon", col = 1:10,
#            alpha = 0.2, label = TRUE)

# Add envfit vectors (only significant ones)
#plot(fit, p.max = 0.05, col = "black")


####################### ADONIS perMANOVA (April, 2025) ############################
# Randomly sample 200 rows per site
set.seed(123) # for reproducibility
VIs_df_sampled <- VIs_df %>%
  group_by(Amoeba_ID) %>%
  slice_sample(n = 200) %>%
  ungroup()

library(vegan)
library(pairwiseAdonis)

# Extract the vegetation indices you want to test
indices_matrix <- VIs_df_sampled %>%
  select(ARI1_1nm, ARI2_1nm, CRI1_1nm, CRI2_1nm, PSRI_1nm, PRI_1nm)

# Run perMANOVA (adonis2)
adonis_results <- adonis2(indices_matrix ~ Amoeba_ID, 
                          data = VIs_df_sampled, 
                          method = "bray", # "bray or "euclidean"
                          permutations = 500)

# Check results
print(adonis_results)

####################### PAIRWISE PERMANOVA (April, 2025) ############################
library(pairwiseAdonis)

# Perform pairwise perMANOVA comparisons
pairwise_results <- pairwise.adonis(indices_matrix, 
                                    factors = VIs_df_sampled$Amoeba_ID, 
                                    sim.method = "euclidean", # "bray or "euclidean"
                                    p.adjust.m = "fdr") # adjust for multiple comparisons

# View pairwise comparison results
pairwise_results
#Write to .csv
write.csv(adonis_results, "./R_outputs/modelling/perMANOVA/adonispma_results_amoebas_v1.csv")
write.csv(pairwise_results,"./R_outputs/modelling/perMANOVA/pairwise_results_amoebas_v1.csv")

############### BOX PLOT CREATION ##################################

###Read in canopy VIs and dendro attributes###
canopy_VIs<-VIs_df
# Scale the indices and combine into a singular stress index (HSI)
indices <- c("ARI1", "ARI2", "CRI1", "CRI2", "PSRI", "PRI")

# Scaling indices to a standard range (0-1)
canopy_VIs_scaled <- canopy_VIs
canopy_VIs_scaled[indices] <- lapply(canopy_VIs[indices], function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
})

# Compute HSI as equally-weighted average
canopy_VIs_scaled$HSI <- rowMeans(canopy_VIs_scaled[indices], na.rm = TRUE)

# Check resulting dataframe
head(canopy_VIs_scaled)

##### Box plot ####
library(tidyr)

# Define custom site order
site_order <- c("GI", "CE", "HI", "BI_North", "BI_South", "CC_North", "CC_South",
                "ET", "FP", "WP", "GW", "RI_IP", "RI_GB")

# Random sample of 1000 pixels per site
set.seed(124)
sampled_data <- canopy_VIs_scaled %>%
  group_by(Sample) %>%
  slice_sample(n = 1000) %>%
  ungroup() %>%
  mutate(Sample = factor(Sample, levels = site_order))  # apply site order once

### 1. PRI Boxplot
ggplot(sampled_data, aes(x = Sample, y = PRI, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4) +
  coord_cartesian(ylim = c(0.45, 0.6)) +
  labs(title = "PRI Distributions Across Coastal Spruce Sites",
       x = "Site",
       y = "PRI Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

### 2. PSRI Boxplot
ggplot(sampled_data, aes(x = Sample, y = PSRI, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4) +
  coord_cartesian(ylim = c(0.275, .425)) +
  labs(title = "PSRI Distributions Across Coastal Spruce Sites",
       x = "Site",
       y = "PSRI Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

### 3. CRI1 Boxplot
ggplot(sampled_data, aes(x = Sample, y = CRI1, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4) +
  coord_cartesian(ylim = c(0.085, .12)) +
  labs(title = "CRI1 Distributions Across Coastal Spruce Sites",
       x = "Site",
       y = "CRI1 Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

install.packages("patchwork")
library(patchwork)

# Create shared theme elements
shared_theme <- theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.text = element_blank())  # removes facet-style titles

# PRI
p1 <- ggplot(sampled_data, aes(x = Sample, y = PRI, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4) +
  coord_cartesian(ylim = c(0.5, 0.6)) +
  ylab("PRI Value") +
  xlab("") +
  shared_theme

# PSRI
p2 <- ggplot(sampled_data, aes(x = Sample, y = PSRI, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4) +
  coord_cartesian(ylim = c(0.2, 0.4)) +
  ylab("PSRI Value") +
  xlab("") +
  shared_theme

# CRI1
p3 <- ggplot(sampled_data, aes(x = Sample, y = CRI1, fill = Sample)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.4) +
  coord_cartesian(ylim = c(0, .2)) +
  ylab("CRI1 Value") +
  xlab("Site") +  # only this one gets the x-axis label
  shared_theme

# Combine vertically
p1 / p2 / p3

