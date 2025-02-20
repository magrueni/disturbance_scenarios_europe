# --------------------------------------------------------------
# Initial landscape Script 01: Prepare Canopy height layer
# --------------------------------------------------------------
# Author: Marc Grünig
# Date: 2025-01-07
# 
# Description:
# This script processes the Lang et al. 2023 canopy height data 
# to 100m resolution (SVD resolution). Hereby we aggregated the 
# pixels as the 80th percentile of the canopy height, approximating
# dominant height.
#
# Notes: for the underlying data see:
# Lang et al., (2023). A high-resolution canopy height model of the Earth.
# 
# --------------------------------------------------------------


### libraries -----------------------------------------------------------------
library(terra)
library(raster)
library(tidyverse)
library(rgeos)
library(sf)
library(lookup)
library(data.table)


# path to folder
path_new <- "/.../"



### Building bouding boxes for the full canopy height dataset ----------------------------
min.long <- -12
min.lat <- 33
max.long <- 42
max.lat <- 72

seq.long <- seq(min.long, max.long, by = 3)
seq.lat <- seq(min.lat, max.lat, by = 3)

combination.min <- expand.grid(seq.long[-length(seq.long)], seq.lat[-length(seq.lat)])
combination.max <- expand.grid(seq.long[-1], seq.lat[-1])

full.combination <- tibble(min.long = combination.min[,1],
                           max.long = combination.max[,1],
                           min.lat = combination.min[,2],
                           max.lat = combination.max[,2])


full.combination <- as.data.frame(full.combination)

bbox <- full.combination %>%
  mutate(left.coord = paste0(ifelse(min.long < 0, "W", "E"), ifelse(abs(min.long) < 10, "00", "0"), round(abs(min.long), 0)),
         top.coord = paste0(ifelse(max.lat < 0, "S", "N"), round(abs(max.lat), 0) - 3)) 

bbox



### load data for each tile ---------------
for(i in 1:nrow(bbox)){
  
  print(paste0(i,  " / ", nrow(bbox)))
  
  name <- paste0(bbox[i, "top.coord"], bbox[i, "left.coord"])
  
  
  if(file.exists(paste0(path_new, "/ETH_GlobalCanopyHeight_10m_2020_", name, "_Map_aggregated_100m.tif"))){next}
  # get height
  
  if(file.exists(paste0(path_new, "/ETH_GlobalCanopyHeight_10m_2020_", name, "_Map.tif"))){
    height_rast <- rast(paste0(path_new, "//ETH_GlobalCanopyHeight_10m_2020_", name, "_Map.tif"))
  }else{
    next
  }
  
  # aggregate and calculate dominant height as 80th precentile
  height_rast_agg <- terra::aggregate(height_rast, fact = 10, cores = 24, fun = function(i,...) quantile(i, probs=0.8, na.rm=T)) 
 
  # double check
  if(origin(height_rast_agg)[1] != origin(height_rast)[1] & origin(height_rast_agg)[2] != origin(height_rast)[2]){
    print("Warning: origin not matching!")
  }
  if(ext(height_rast_agg) != ext(height_rast)){
    print("Warning: extent not matching!")
  }
  
  # save 
  writeRaster(height_rast_agg, paste0(path_new, "/ETH_GlobalCanopyHeight_10m_2020_", name, "_Map_aggregated_100m.tif"), overwrite = T)
  
  rm(height_rast_agg)
  tmpFiles(remove = T)
  gc()
  
  
  
}


tmpFiles(remove = T)


### end ------------------------------------------------------------------------