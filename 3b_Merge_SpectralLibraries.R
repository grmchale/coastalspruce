library(purrr)
library(dplyr)

# ------------- MERGE SPECTRAL LIBRARIES ------------------------------------
NM_values <- c("1nm", "5nm", "10nm", "15nm")

# Read, clean, and rename columns before merging
VIstats_list <- lapply(NM_values, function(nm) {
  df <- read.csv(paste0("./R_outputs/speclib_dendrometers/dendrometer_canopy_bytreeid_", nm, ".csv"))
  
  # Remove "CC" and "Species" columns if they exist
  df <- df %>% select(-any_of(c("CC", "Species", "X")))
  
  # Rename all columns except "TreeID" by prefixing with nm value
  colnames(df) <- ifelse(colnames(df) == "TreeID", "TreeID", paste0(nm, "_", colnames(df)))
  
  return(df)
})

# Merge all dataframes using inner join on "TreeID"
VIstats_combined <- reduce(VIstats_list, ~inner_join(.x, .y, by = "TreeID"))

#Export the joined boys as .csv and dataframe
write.csv(VIstats_combined, "./R_outputs/speclib_dendrometers/dendrometer_canopy_bytreeid_MERGED.csv")

speclibraries_combined <- VIstats_combined
VIstats_merged <-read.csv("./R_outputs/speclib_dendrometers/veg_indices_stats/dendrometer_VIstats_MERGED.csv")

# Merge the dataframes using TreeID as the key
final_combined <- inner_join(speclibraries_combined, VIstats_merged, by = "TreeID")

write.csv(final_combined,"./R_outputs/speclib_dendrometers/dendrometer_canopy_bytreeid-speclibs_MERGED.csv")

#Just a quick scatterplot
ggplot(df_rf, aes(x = DPI_5nm_Q25, y = TWD)) +
  geom_point() +  # Scatter plot points
  labs(x = "DPI_5nm_Q25", y = "TWD", title = "Scatterplot of DPI_5nm_Q25 vs TWD") +
  theme_minimal()  # Clean theme

# ----------Structural and spectral characteristics???? HUHHHHHH???-------------
library(dplyr)
library(readr)
install.packages("Hmisc")
library(Hmisc)

# Load structural metrics (LiDAR)
lidar_metrics <- read_csv("G:/LiDAR/Metrics/Dendro_Crown_Metrics.csv")

# Load spectral metrics (Vegetation Indices)
spectral_metrics <- read_csv("G:/HyperspectralUAV/R_outputs/speclib_dendrometers/veg_indices/dendrometer_VIs_5nm_bytreeid.csv")

# Join datasets using TreeID
combined_data <- lidar_metrics %>%
  inner_join(spectral_metrics, by = "TreeID") %>%
  select(-c(58:67))

# Separate structural and spectral variables
structural_vars <- lidar_metrics %>%
  select(-TreeID) %>%
  select(where(is.numeric))

spectral_vars <- spectral_metrics %>%
  select(-TreeID) %>%
  select(where(is.numeric))

# Combine structural and spectral metrics (order matters!)
combined_numeric <- cbind(structural_vars, spectral_vars)

# Run Spearman correlations
cor_results <- rcorr(as.matrix(combined_numeric), type = "spearman")

# Get index positions:
structural_cols <- 1:ncol(structural_vars)
spectral_cols <- (ncol(structural_vars) + 1):ncol(combined_numeric)

# Now correctly subset the correlation matrix:
cor_matrix_structural_spectral <- cor_results$r[rownames(cor_results$r)[1:ncol(structural_vars)],
                                                colnames(cor_results$r)[(ncol(structural_vars)+1):ncol(combined_numeric)]]

# Extract p-values in the same manner:
p_values_structural_spectral <- cor_results$P[rownames(cor_results$P)[1:ncol(structural_vars)],
                                              colnames(cor_results$P)[(ncol(structural_vars)+1):ncol(combined_numeric)]]

# Check outputs:
print(cor_matrix_structural_spectral)
print(p_values_structural_spectral)

# --------------------Just a quick scatterplot for spectral vs structural------
library(ggplot2)

# Compute R² values manually
lm_dpi <- lm(DPI ~ isd, data = combined_data)
r2_dpi <- summary(lm_dpi)$r.squared

lm_ci <- lm(CI ~ isd, data = combined_data)
r2_ci <- summary(lm_ci)$r.squared

# Scatterplot for ISD vs DPI
plot1 <- ggplot(combined_data, aes(x = isd, y = DPI)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Linear trendline
  annotate("text", x = min(combined_data$isd), y = max(combined_data$DPI), 
           label = paste("R² =", round(r2_dpi, 3)), hjust = 0, size = 5) +
  labs(x = "ISD", y = "DPI", title = "ISD vs DPI") +
  theme_minimal()

# Scatterplot for ISD vs CI
plot2 <- ggplot(combined_data, aes(x = isd, y = CI)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +  # Linear trendline
  annotate("text", x = min(combined_data$isd), y = max(combined_data$CI), 
           label = paste("R² =", round(r2_ci, 3)), hjust = 0, size = 5) +
  labs(x = "ISD", y = "CI", title = "ISD vs CI") +
  theme_minimal()

# Display the plots side by side
library(gridExtra)
grid.arrange(plot1, plot2, ncol = 2)


#-------------- CORRELATION MATRIX!!--------------------------------------------

library(plotly)
library(reshape2)

# reshape correlation matrix into long format
cor_long <- melt(cor_matrix_structural_spectral, varnames = c("Structural_Metric", "Spectral_Metric"), value.name = "Correlation")

# Interactive heatmap
plot_ly(cor_long, x = ~Spectral_Metric, y = ~Structural_Metric, z = ~Correlation,
        type = "heatmap", colors = colorRamp(c("blue", "white", "red"))) %>%
  layout(title = "Interactive Correlation Heatmap: Structural vs. Spectral Metrics")

library(dplyr)
library(igraph)

# Assume cor_matrix_structural_spectral already exists:
library(reshape2)

# Convert correlation matrix to a dataframe with edges clearly specified
cor_long <- melt(cor_matrix_structural_spectral, 
                 varnames = c("Structural", "Spectral"),
                 value.name = "Correlation")

# Filter strong correlations only (absolute correlation > 0.5; adjust as needed)
cor_filtered <- cor_long %>% 
  filter(abs(Correlation) > 0.5)

# Prepare edge list clearly: FROM structural metric TO spectral metric
edges <- cor_filtered %>%
  select(Structural, Spectral, Correlation)

# Build the igraph object explicitly from edges:
library(igraph)

# Create a unique list of vertices
vertices <- data.frame(name = unique(c(edges$Structural, edges$Spectral)),
                       type = ifelse(unique(c(edges$Structural, edges$Spectral)) %in% cor_filtered$Structural,
                                     "Structural", "Spectral"))

# Create graph object
graph <- graph_from_data_frame(d = edges, vertices = vertices, directed = FALSE)

# Set node colors based on vertex type (structural vs spectral)
V(graph)$color <- ifelse(V(graph)$type == "Structural", "skyblue", "salmon")

# Set edge widths based on strength of correlation
E(graph)$width <- abs(edges$Correlation) * 3

# Plot clearly
plot(graph,
     layout = layout_with_fr,
     vertex.label.cex = 0.8,
     vertex.size = 15,
     edge.width = E(graph)$Correlation * 3, # Adjusted visually based on strength
     main = "Structural vs Spectral Correlation Network (>0.5)",
     vertex.label.color = "black")

# ---------------------- MAP MOISTURE INDEX AGAINST TWD ------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

drone_dendro<-read.csv("G:/Dendrometers/drone_dendro.csv")
VIstats_merged <-read.csv("./R_outputs/speclib_dendrometers/veg_indices_stats/dendrometer_VIstats_MERGED.csv")

# Left join the data frames by TreeID
merged_df <- drone_dendro %>%
  left_join(VIstats_merged, by = "TreeID")

# Function to pivot and create scatterplots for specified WBI columns
plot_wbi_scatter <- function(data, xvar, yvars, measure_label) {
  data_long <- data %>%
    pivot_longer(
      cols = all_of(yvars),
      names_to = "WBI_Spacing",
      values_to = "WBI_Value"
    )
  
  ggplot(data_long, aes_string(x = xvar, y = "WBI_Value")) +
    geom_point() +
    facet_wrap(~ WBI_Spacing, scales = "free_y") +
    labs(
      title = paste("TWD vs. WBI", measure_label, "at Different nm Spacings"),
      x = "TWD",
      y = paste("WBI", measure_label)
    ) +
    theme_minimal()
}

# 1) Create scatter plots for the WBI Median columns
plot_wbi_scatter(
  data          = merged_df,
  xvar          = "TWD",
  yvars         = c("WBI_1nm_Median", "WBI_5nm_Median", "WBI_10nm_Median", "WBI_15nm_Median"),
  measure_label = "Median"
)

# 2) Create scatter plots for the WBI Mean columns
plot_wbi_scatter(
  data          = merged_df,
  xvar          = "TWD",
  yvars         = c("WBI_1nm_Mean", "WBI_5nm_Mean", "WBI_10nm_Mean", "WBI_15nm_Mean"),
  measure_label = "Mean"
)

