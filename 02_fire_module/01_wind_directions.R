# --------------------------------------------------------------
# Fire Module Scipt 01: Wind data preparation for the fire module
# --------------------------------------------------------------
# Author: Marc Grünig
# Date: 2025-01-07
# 
# Description:
# This script processes wind data from ERA5 to generate mean wind 
# speed and direction over a 100km reference grid. It includes:
# 
# 1. Reprojecting wind rasters to match a reference grid (EPSG:3035).
# 2. Calculating mean U and V wind components across multiple layers.
# 3. Deriving wind speed and direction for each grid cell.
# 4. Aggregating results to compute mean and standard deviation 
#    of wind speed and direction across time steps.
# 5. Exporting the processed wind data to a CSV file.
# 
# The resulting outputs are used in SVD to define wind directions
# 
# --------------------------------------------------------------



# load libraries
library(raster)
library(terra)
library(sf)
library(rgdal)
library(data.table)


# define path
path <- "/.../"

# Load wind rasters and reference grid
wind_rasters <- rast(file.path(path, "/era5_wind.nc"))
ref_grid_100km <- rast(file.path(path, "/07_reference_grids/reference_grid_100km.tif"))

# Reproject wind rasters to match reference grid
wind_rasters <- terra::project(wind_rasters, y = "epsg:3035", method = "bilinear")


# Define helper functions
calculate_wind_speed <- function(u, v) {
  sqrt(u^2 + v^2)
}

calculate_wind_direction <- function(u, v) {
  (180 + atan2(u, v) * 180 / pi) %% 360
}


# Process U and V components
process_wind_component <- function(component_rasters, ref_grid) {
  # mean_component <- mean(component_rasters)
  component_rasters <- terra::project(component_rasters, ref_grid)
  component_rasters <- terra::mask(component_rasters, ref_grid)
  return(component_rasters)
}

u_winds <- wind_rasters[[grepl("u10", names(wind_rasters))]]
v_winds <- wind_rasters[[grepl("v10", names(wind_rasters))]]

u_mean <- process_wind_component(u_winds, ref_grid_100km)
v_mean <- process_wind_component(v_winds, ref_grid_100km)

# Calculate wind speed and direction
wind_speed <- calculate_wind_speed(u_mean, v_mean)
wind_direction <- calculate_wind_direction(u_mean, v_mean)

# Initialize data frames
direction_df <- as.data.frame(ref_grid_100km, xy = TRUE)
speed_df <- direction_df

# Loop through layers to calculate wind direction and speed for all time steps
for (i in 1:nlyr(u_mean)) {
  u <- u_mean[[i]]
  v <- v_mean[[i]]
  
  wind_dir <- calculate_wind_direction(u, v)
  wind_spe <- calculate_wind_speed(u, v)
  
  dir_df <- as.data.frame(wind_dir, xy = TRUE)
  spe_df <- as.data.frame(wind_spe, xy = TRUE)
  
  direction_df <- left_join(direction_df, dir_df, by = c("x", "y"))
  speed_df <- left_join(speed_df, spe_df, by = c("x", "y"))
}

# Calculate statistics for wind direction and speed
calculate_statistics <- function(df) {
  df %>%
    rowwise() %>%
    mutate(
      mean_value = mean(c_across(starts_with("u10"))),
      sd_value = sd(c_across(starts_with("u10")))
    ) %>%
    ungroup() %>%
    select(gridid, mean_value, sd_value)
}

direction_stats <- calculate_statistics(direction_df)
speed_stats <- calculate_statistics(speed_df)

# Combine direction and speed statistics
wind_df <- left_join(direction_stats, speed_stats, by = "gridid")
colnames(wind_df) <- c("gridid", "mean_dir", "sd_dir", "mean_speed", "sd_speed")

# save result
write.csv(wind_df, paste0(path, "/04_disturbance_modules/fire_module/wind_table.csv"), row.names = FALSE)
