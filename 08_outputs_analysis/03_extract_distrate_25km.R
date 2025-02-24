# ------------------------------------------------------------------------------
# Script Name: extract_distrate_25km
#
# Description: This scripts take the 10-year timesteps calculated
# in script 01 and aggregates the data to 25km hexagons.
# 
# Author: Marc Grünig
# Date: 7.2.2025
#
# Input Data:
#   - Disturbance simulation outputs (raster files)
#   - Europe scale forest mask and other GIS layers
#
#
# Output:
#   - CSV files with disturbance rates, areas, and frequencies for each year and RCP scenario
#   - Disturbance metrics including the percentage of disturbed area, disturbance rate, and frequency
#
# Processing Steps:
#   1. Load necessary geospatial data and simulation outputs.
#   2. Loop through each RCP (climate scenario) and year, extract disturbance data.
#   3. Compute disturbance rates and frequencies based on simulation results.
#   4. Save the calculated disturbance data in CSV format for further analysis.
#
# ------------------------------------------------------------------------------



### libraries ------------------------------------------------------------------
library(rnaturalearth)
library(rnaturalearthdata)
library(MetBrewer)
library(tidyterra)
library(gridExtra)
library(raster)
library(ggplot2)
library(sf)
library(terra)
library(stringr)



### load data ------------------------------------------------------------------


# define path
path <- "/.../"
path_results <- paste0(path, "/svd_simulations/results_eu/"

# load masks and hexagons
eu_shp <- vect(paste0(path,"/reference_grids/eu_mask.shp"))
hex_ecu <- shapefile(paste0(path, "/reference_grids/eu_mask_hexagons_25km.shp"))

sf_obj_forest_mask <- read_sf(paste0(path, "/reference_grids/hex_forest_mask_25km.gpkg"))
sf_obj_forest_mask_ids <- sf_obj_forest_mask %>% st_drop_geometry()


master_tab <- read.csv(paste0(path, "/svd_simulations/svd_simulations_ids.csv"), sep = ";")


### then all together --- 
files_wind <- list.files(paste0(path_results, "/wind/"), pattern = "*_10year.tif", full.names = T)
files_fire <- list.files(paste0(path_results, "/fire/"), pattern = "*_10year.tif", full.names = T)
files_bbtl <- list.files(paste0(path_results, "/bbtl/"), pattern = "*_10year.tif", full.names = T)

rcps <- c("rcp85", "rcp45", "rcp26", "historical")
gcms <- c("ncc", "ichec", "mpi")


# loop over rcps
for(c in 1:length(rcps)){
  
  print(c)
  
  # get files from the rcp
  wind_fls_rcps <- files_wind[grepl(paste0(rcps[c]), files_wind)]
  fire_fls_rcps <- files_fire[grepl(paste0(rcps[c]), files_fire)]
  bbtl_fls_rcps <- files_bbtl[grepl(paste0(rcps[c]), files_bbtl)]
  
  
  # create empty lists
  mean_plot_list <- list()
  sd_plot_list <- list()
  r <- 0
  mean_year_all <- c()
  all_maps_dat <- list()
  
  # loop over years
  for(y in seq(10, 80, 10)){
    
    if(y < 80){next}
    # counter and progress
    r <- r + 1
    print(y)
    
    # get files from the timestep
    fls_wind <- wind_fls_rcps[grepl(paste0("wind_", y, "_"), wind_fls_rcps)]
    fls_fire <- fire_fls_rcps[grepl(paste0("fire_", y, "_"), fire_fls_rcps)]
    fls_bbtl <- bbtl_fls_rcps[grepl(paste0("bbtl_", y, "_"), bbtl_fls_rcps)]
    fls_all <- c(fls_wind, fls_fire, fls_bbtl)
    rastas <- rast(fls_all)
    hex_srtm <- as.data.frame(exactextractr::exact_extract(rastas, hex_ecu, "sum", progress = F))

    rm(rastas)
    gc()
    
    
    hex_srtm_mean <- hex_srtm %>% 
      mutate(gridid = 1:nrow(.)) %>% 
      left_join(., sf_obj_forest_mask_ids, by = "gridid")
    
    long_df <- hex_srtm_mean %>%
      tidyr::pivot_longer(
        cols = starts_with("sum"),
        names_to = "variable",
        values_to = "value"
      ) %>%
      mutate(
        agent = str_extract(variable, "(?<=^sum\\.)[^_]+"),
        sim = str_extract(variable, "\\d+$")  # Extract the trailing digits
      ) %>% 
      mutate(dist_rate = as.numeric(100 / forest_area * value / 10),
             dist_freq = as.numeric(10 * forest_area / value)) %>% 
      mutate(year = 2020 + y,
             rcp = rcps[c],
             gcm = master_tab[sim, "GCM"],
             rep = master_tab[sim, "Rep"]) %>% 
      dplyr::select(gridid, agent, rcp, gcm, rep, year, sim,
                    dist_rate, dist_freq, dist_area = value) %>% 
      mutate(dist_rate = ifelse(is.infinite(dist_rate), NA, dist_rate),
             dist_rate = ifelse(is.infinite(dist_freq), NA, dist_rate))

    write_csv(long_df, paste0(path_results, "/dist_rates_25km/all_dist_rates_", y, "_", rcps[c], ".csv"))
    
    rm(long_df, hex_srtm_mean, hex_srtm)
    gc()
    
  } # close year loop
  
} # close rcp


### end ------------------------------------------------------------------------






