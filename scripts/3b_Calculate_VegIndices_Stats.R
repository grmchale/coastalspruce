library(dplyr)

### CALCULATE SD, MEAN, AND QUANTILES FOR VEGETATION INDICES ###

# base path and target nm spacings
base_path <- "./R_outputs/speclib_chronologies/veg_indices/"
NMs <- c("1nm", "5nm", "10nm", "15nm")

for (NM in NMs) {
  # READ IN VEGETATION INDICES (per-pixel)
  VIs_df <- readRDS(paste0(base_path, "chronologies_VIs_", NM, ".rds"))
  
  # Assuming columns 3:end are vegetation indices
  VI_cols <- names(VIs_df)[3:ncol(VIs_df)]
  
  # CALCULATE SD, MEAN, AND QUANTILES FOR VEGETATION INDICES W/ ARIs
  VI_stats_df <- VIs_df %>%
    group_by(Sample) %>%
    summarize(
      across(
        all_of(VI_cols),
        list(
          Mean   = ~ mean(.x, na.rm = TRUE),
          SD     = ~ sd(.x, na.rm = TRUE),
          Q25    = ~ quantile(.x, 0.25, na.rm = TRUE),
          Median = ~ quantile(.x, 0.50, na.rm = TRUE),
          Q75    = ~ quantile(.x, 0.75, na.rm = TRUE)
        ),
        .names = paste0("{.col}_", NM, "_{.fn}")
      ),
      .groups = "drop"
    )
  
  # create objects named chronologies_VIstats_#nm in the environment
  assign(paste0("chronologies_VIstats_", NM), VI_stats_df)
}

# (optional) quick confirmation
cat("Created: ", paste0("chronologies_VIstats_", NMs, collapse = ", "), "\n")

# -----------Merge 1, 5, 10, 15 nm spectral libraries OR stats(?)---------------
library(dplyr)
library(purrr)

NM_values  <- c("1nm", "5nm", "10nm", "15nm")
df_names   <- paste0("chronologies_VIstats_", NM_values)

# grab the existing data frames from the environment
VIstats_list <- mget(df_names, inherits = TRUE)

# inner-join on Sample across all four
chronologies_VIstats_combined <- reduce(VIstats_list, ~ inner_join(.x, .y, by = "Sample"))

# quick peek
dim(chronologies_VIstats_combined)
head(names(chronologies_VIstats_combined), 20)

# Save the combined dataframe explicitly as a single RDS file
saveRDS(chronologies_VIstats_combined, "./R_outputs/speclib_chronologies/veg_indices/chronologies_VIstats_combined.rds")

###################### AMOEBAS JOINING SPECLIBS + VIS ###################################
# Read in spectal libraries
rds_path <- "./R_outputs/speclib_amoebas_final/dataframes"
canopy_image_spectra_df      <- readRDS(file.path(rds_path, "amoebas_speclib_1nm.rds"))
canopy_image_spectra_5nm_df  <- readRDS(file.path(rds_path, "amoebas_speclib_5nm.rds"))
canopy_image_spectra_10nm_df <- readRDS(file.path(rds_path, "amoebas_speclib_10nm.rds"))
canopy_image_spectra_15nm_df <- readRDS(file.path(rds_path, "amoebas_speclib_15nm.rds"))
# Read in vegetation indices
rds_path <- "G:/HyperspectralUAV/R_outputs/speclib_amoebas_final/veg_indices"
amoebas_VIs_1nm_df  <- readRDS(file.path(rds_path, "amoebas_VIs_1nm.rds"))
amoebas_VIs_5nm_df  <- readRDS(file.path(rds_path, "amoebas_VIs_5nm.rds"))
amoebas_VIs_10nm_df <- readRDS(file.path(rds_path, "amoebas_VIs_10nm.rds"))
amoebas_VIs_15nm_df <- readRDS(file.path(rds_path, "amoebas_VIs_15nm.rds"))

suppressPackageStartupMessages({ library(dplyr) })

# ---- helper: suffix non-key cols ----
suffix_nonkeys <- function(df, suffix, keys = c("Sample","sample_name","Amoeba_ID")) {
  keep <- intersect(keys, names(df))
  ren  <- setdiff(names(df), keep)
  names(df)[match(ren, names(df))] <- paste0(ren, suffix)
  df
}

# 1) Suffix columns by resampling
v1 <- suffix_nonkeys(amoebas_VIs_1nm_df,  "_1nm")
v5 <- suffix_nonkeys(amoebas_VIs_5nm_df,  "_5nm")
v10<- suffix_nonkeys(amoebas_VIs_10nm_df, "_10nm")
v15<- suffix_nonkeys(amoebas_VIs_15nm_df, "_15nm")

s1 <- suffix_nonkeys(canopy_image_spectra_df,      "_1nm")
s5 <- suffix_nonkeys(canopy_image_spectra_5nm_df,  "_5nm")
s10<- suffix_nonkeys(canopy_image_spectra_10nm_df, "_10nm")
s15<- suffix_nonkeys(canopy_image_spectra_15nm_df, "_15nm")

# 2) Keep Amoeba_ID only from the first (v1); drop from others to avoid .x/.y
drop_ai <- function(df) { if ("Amoeba_ID" %in% names(df)) df <- df %>% select(-Amoeba_ID); df }
v5 <- drop_ai(v5); v10 <- drop_ai(v10); v15 <- drop_ai(v15)
s1 <- drop_ai(s1); s5  <- drop_ai(s5);  s10 <- drop_ai(s10); s15 <- drop_ai(s15)

# 3) Join (inner join on keys so rows match across all)
keys <- c("Sample","sample_name")
joined <- v1 %>%
  inner_join(v5,  by = keys) %>%
  inner_join(v10, by = keys) %>%
  inner_join(v15, by = keys) %>%
  inner_join(s1,  by = keys) %>%
  inner_join(s5,  by = keys) %>%
  inner_join(s10, by = keys) %>%
  inner_join(s15, by = keys)

# 4) Column order: Sample, Amoeba_ID, sample_name, VIs (1nm→5nm→10nm→15nm), Spectra (1nm→5nm→10nm→15nm)
vi_1_cols  <- setdiff(names(v1),  c("Sample","sample_name","Amoeba_ID"))
vi_5_cols  <- setdiff(names(v5),  c("Sample","sample_name"))
vi_10_cols <- setdiff(names(v10), c("Sample","sample_name"))
vi_15_cols <- setdiff(names(v15), c("Sample","sample_name"))

sp_1_cols  <- setdiff(names(s1),  c("Sample","sample_name"))
sp_5_cols  <- setdiff(names(s5),  c("Sample","sample_name"))
sp_10_cols <- setdiff(names(s10), c("Sample","sample_name"))
sp_15_cols <- setdiff(names(s15), c("Sample","sample_name"))

final_order <- c("Sample", "Amoeba_ID", "sample_name",
                 vi_1_cols, vi_5_cols, vi_10_cols, vi_15_cols,
                 sp_1_cols, sp_5_cols, sp_10_cols, sp_15_cols)

# guard in case any columns are missing from joins
final_order <- final_order[final_order %in% names(joined)]
joined <- joined[, final_order]

# 5) Export
out_dir <- "G:/HyperspectralUAV/R_outputs/speclib_amoebas_final"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

csv_path <- file.path(out_dir, "amoebas_VIs_and_spectra_combined.csv")
rds_path <- file.path(out_dir, "amoebas_VIs_and_spectra_combined.rds")

write.csv(joined, csv_path, row.names = FALSE)
saveRDS(joined, rds_path)

message("Wrote:\n  ", csv_path, "\n  ", rds_path)











