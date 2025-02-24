# ------------------------------------------------------------------------------
# Script Name: old_young_forests.R
# Description: This script processes the spatial simulation outputs to extract
# the share of young and old forests.
#
# Author: Marc Grünig
# Date: 7.2.2025
#
# Input Files:
#   - Raw simulation raster files for wind, fire, bark beetle outbreaks, and
#   forest management impacts.
#
# Output:
#   - Spatial data files (GPKG) summarizing young and old forests 
#     for each simulation run.
#
# Notes:
#   - We stored some example outputs in the biodiv folder as calculations are
#     are computationally intensive.
#
# ------------------------------------------------------------------------------



### libraries ------------------------------------------------------------------

library(tidyterra)
library(gridExtra)
library(raster)
library(ggplot2)
library(sf)
library(terra)
library(parallel)
library(exactextractr)
library(readr)
library(rnaturalearth)
library(rnaturalearthdata)
library(MetBrewer)


### settings ---------------------------------------------------------------------
path <- "/.../"
path_sims <- paste0(path, "/svd_simulations/raw/")
path_results <-  paste0(path, "/svd_simulations/results_eu/")

# terra options
tmpFiles(remove = T)


### load some data ---------------------------------------------------------------
eu_shp <- vect(paste0(path, "/reference_grids/eu_mask.shp"))
hex_ecu <- shapefile(paste0(path, "/reference_grids/eu_mask_hexagons_25km.shp"))

# plot(hex_ecu)
proj_laea <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
sf_obj_forest_mask <- st_read(paste0(path, "/reference_grids/hex_forest_mask_25km.gpkg"))
sf_obj_forest_mask_ids <- sf_obj_forest_mask %>% st_drop_geometry()


# load some data
master_tab <- read.csv(paste0(path, "/svd_simulations/svd_simulations_ids.csv"), sep = ";")
ref_grid <- rast(paste0(path, "reference_grids/reference_grid_10km.tif"))

# define agents to loop through
agents <- c("wind", "fire", "bbtl", "mgmt")

# define years
years <- seq(10, 80, 10)


### start processing ----------------------------------------------------------------
for(i in c(1:120)){
  
  print(i)
  rcp <- master_tab[i, "RCP"]
  gcm <- master_tab[i, "GCM"]
  rep <- master_tab[i, "Rep"]
  
  if(file.exists(paste0(path_results, "/biodiv/undist_vs_recently_dist_forest_", i, ".gpkg"))){next}
  
  
  out_list_sim <- list()
  mean_undist_sim <- list()
  
  for(y in 1:length(years)){
    
    year <- years[y]
    print(year)
    
    # total dist recently - get the 10y timestep file
    wind <- rast(paste0(path_results, "wind/wind_", year, "_", rcp, "_", gcm, "_", rep, "_10year.tif"))
    fire <- rast(paste0(path_results, "fire/fire_", year, "_", rcp, "_", gcm, "_", rep, "_10year.tif"))
    bbtl <- rast(paste0(path_results, "bbtl/bbtl_", year, "_", rcp, "_", gcm, "_", rep, "_10year.tif"))
    mgmt <- rast(paste0(path_results, "mgmt/mgmt_", year, "_", rcp, "_", gcm, "_", rep, "_10year.tif"))
    
    total_dist <- c(wind, fire, bbtl, mgmt)
    total_dist <- app(total_dist, "sum", cores = 36)
    total_dist[total_dist > 1] <- 1
    
    # total disturbed over whole simulation, get output at the end of simulation period
    if(rcp == "historical"){
      wind <- rast(paste0(path_lorien, "_historical/output_sim_eu_", i, "/ccwind/wind_year_", year, "_sim_eu_", i, ".tif"))
      fire <- rast(paste0(path_lorien, "_historical/output_sim_eu_", i, "/ccfire/firegrid_", year, "_sim_eu_", i, ".tif"))
      bbtl <- rast(paste0(path_lorien, "_historical/output_sim_eu_", i, "/bb_out/bbgrid_", year, "_sim_eu_", i, ".tif"))
      mgmt <- rast(paste0(path_lorien, "_historical/output_sim_eu_", i, "/mgmt_out/mgmtgrid_", year, "_sim_eu_", i, ".tif"))
    }else{
      wind <- rast(paste0(path_lorien, "/output_sim_eu_", i, "/ccwind/wind_year_", year, "_sim_eu_", i, ".tif"))
      fire <- rast(paste0(path_lorien, "/output_sim_eu_", i, "/ccfire/firegrid_", year, "_sim_eu_", i, ".tif"))
      bbtl <- rast(paste0(path_lorien, "/output_sim_eu_", i, "/bb_out/bbgrid_", year, "_sim_eu_", i, ".tif"))
      mgmt <- rast(paste0(path_lorien, "/output_sim_eu_", i, "/mgmt_out/mgmtgrid_", year, "_sim_eu_", i, ".tif"))
    }
    
    
    # add together
    total_undist <- c(wind, fire, bbtl, mgmt)
    total_undist <- app(total_undist, "sum", cores = 36)
    total_undist[total_undist > 1] <- 1
    

    # extract by hexagon
    hex_srtm_dist_recently <- as.data.frame(exactextractr::exact_extract(total_dist, hex_ecu, "sum", progress = F))
    colnames(hex_srtm_dist_recently) <- "dist_recently"
    
    hex_srtm_undist <- as.data.frame(exactextractr::exact_extract(total_undist, hex_ecu, "sum", progress = F))
    colnames(hex_srtm_undist) <- "disturbed"
    
    hex_srtm <- cbind(hex_srtm_dist_recently, hex_srtm_undist)
    
    hex_srtm_mean <- hex_srtm %>% 
      mutate(gridid = 1:nrow(.)) %>% 
      left_join(., sf_obj_forest_mask_ids, by = "gridid") %>%
      rowwise() %>%
      mutate(dist_prct = 100/forest_area * disturbed,
             dist_prct = ifelse(dist_prct < 0, 0, 
                                ifelse(dist_prct > 100, 100, dist_prct))) %>% 
      
      mutate(undisturbed = forest_area - disturbed) %>%
      mutate(undist_prct = 100/forest_area * undisturbed) %>% 
      mutate(undist_prct = ifelse(undist_prct < 0, 0, undist_prct)) %>% 
      
      mutate(dist_rec_prct = 100/forest_area * dist_recently) %>% 
      mutate(dist_rec_prct = ifelse(dist_rec_prct < 0, 0, dist_rec_prct)) %>% 
      
      mutate(year = 2020 + year) %>% 
      mutate(dist_prct = dist_prct - dist_rec_prct) %>% 
      dplyr::select(gridid, year, disturbed, dist_prct, undisturbed, undist_prct, dist_recently, dist_rec_prct)

    
    out_df <- hex_srtm_mean %>% 
      mutate(rcp = rcp,
             gcm = gcm, 
             rep = rep)
    
    
    out_list_sim[[y]] <- out_df

    rm(out_df)
    tmpFiles(remove = T)
    for (z in 1:5){gc(reset = T)} 
    
    
  }

  data_sim <- do.call(rbind, out_list_sim)
  write_sf(data_sim, paste0(path_results, "/biodiv/undist_vs_recently_dist_forest_", i, ".gpkg"))
  
  rm(data_sim)
  tmpFiles(remove = T)
  for (z in 1:5){gc(reset = T)} 
  
}


### end ------------------------------------------------------------------------