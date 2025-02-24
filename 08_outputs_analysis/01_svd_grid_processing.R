# --------------------------------------------------------------------------------------
# Script Name: 01_svd_grid_processing
# Description: This script processes the spatial simulation data outputs for all modules.
# first, the future simulations, in a second step the historical simulations.
# 
# The spatial outputs are at 100 x 100m scale across Europe, therefore large .tif files.
# These outputs are written at 10-year timestep by SVD.

# To process them, it is important to note, that the outputs for wind and fire disturbances
# show cumulated number of disturbances per grid cell. Meanwhile, bark beetle and 
# management outputs show years since last disturbance. Therefore, the processing is
# slightly different.
#
# 
# Author: Marc Grünig
# Date: 7.2.2025
#
# Dependencies:
#   - R packages: tidyverse, terra
#
# Output:
#   - Raster files representing 10-year timestep changes in wind, fire, bark beetle 
#     outbreaks, and forest management interventions.
#   - Summed impact raster files for specific years (40, 80).
#
# Notes:
#   - Computationally intensive, temporary files are periodically removed
#     to manage memory usage.
#
# --------------------------------------------------------------------------------------



### libraries -----------------------------------------------------------
library(tidyverse)
library(terra)

### settings ---------------------------------------------------------------
path <- "/.../"
path_public <- "/.../"
path_out <- paste0(path_public, "/svd_simulations/results_eu/")


master_tab <- read.csv(paste0(path_public, "/svd_simulations/svd_simulations_ids.csv"), sep = ";")


### start processing ----------------------------------------------------------

## wind --

# extract the wind per 10 year timestep
for(i in sims){
  
  print(i)
  
  rcp <- master_tab[i, "RCP"]
  gcm <- master_tab[i, "GCM"]
  rep <- master_tab[i, "Rep"]

  
  for(y in seq(10, 80, 10)){
    
    if(file.exists(paste0(path_out, "wind/wind_", y, "_", rcp, "_", gcm, "_", rep, "_10year.tif"))){next}
    if(y == 10){
      wind_rast <- rast(paste0(path, "/output_sim_eu_", i, "/ccwind/wind_year_", y, "_sim_eu_", i, ".tif"))
      wind_rast_focal <- wind_rast
    }else{
      
      if(!file.exists(paste0(path, "/output_sim_eu_", i, "/ccwind/wind_year_", y, "_sim_eu_", i, ".tif"))){
        print(paste0("sim ", i, " incomplete!"))
        next}
      r <- rast(paste0(path, "/output_sim_eu_", i, "/ccwind/wind_year_", y, "_sim_eu_", i, ".tif"))
      
      # in some years we want to know the total
      if(y %in% c(40, 80)){
        writeRaster(r, paste0(path_out, "wind/wind_", y, "_", rcp, "_", gcm, "_", rep, "_summed.tif"), overwrite = T)
      }
      
      wind_rast_focal <- r - wind_rast
      wind_rast <- r
      
    }

    writeRaster(wind_rast_focal, paste0(path_out, "wind/wind_", y, "_", rcp, "_", gcm, "_", rep, "_10year.tif"), overwrite = T)

    gc()
  }
  tmpFiles(remove = T)
  gc()
}




## fire --
sims <- c(1:90)

for(i in sims){
  
  print(i)

  rcp <- master_tab[i, "RCP"]
  gcm <- master_tab[i, "GCM"]
  rep <- master_tab[i, "Rep"]
  
  for(y in seq(10, 80, 10)){

    
    if(y == 10){
      fire_rast <- rast(paste0(path, "/output_sim_eu_", i, "/ccfire/firegrid_", y, "_sim_eu_", i, ".tif"))
      fire_rast_focal <- fire_rast
    }else{
      
      if(!file.exists(paste0(path, "/output_sim_eu_", i, "/ccfire/firegrid_", y, "_sim_eu_", i, ".tif"))){
        print(paste0("sim ", i, " incomplete!"))
        next}
      r <- rast(paste0(path, "/output_sim_eu_", i, "/ccfire/firegrid_", y, "_sim_eu_", i, ".tif"))
      
      
      if(y %in% c(40, 80)){
        writeRaster(r, paste0(path_out, "fire/fire_", y, "_", rcp, "_", gcm, "_", rep, "_summed.tif"), overwrite = T)
      }
      
      
      fire_rast_focal <- r - fire_rast
      fire_rast <- r
      
      gc()
    }
    
    writeRaster(fire_rast_focal, paste0(path_out, "fire/fire_", y, "_", rcp, "_", gcm, "_", rep, "_10year.tif"), overwrite = T)

    gc()
  }
  
  tmpFiles(remove = T)
  gc()
}



## bark beetle --
sims <- c(1:90)

for(i in sims){
  
  print(i)
  if(!file.exists(paste0(path, "/output_sim_eu_", i, "/bb_out/bbgrid_10_sim_eu_", i, ".tif"))){next}

  rcp <- master_tab[i, "RCP"]
  gcm <- master_tab[i, "GCM"]
  rep <- master_tab[i, "Rep"]
  
  if(file.exists(paste0(path_out, "bbtl/bbtl_80_", rcp, "_", gcm, "_", rep, "_summed.tif"))){next}
  
  
  for(y in seq(10, 80, 10)){

    if(y == 10){
      bbtl_rast <- rast(paste0(path, "/output_sim_eu_", i, "/bb_out/bbgrid_", y, "_sim_eu_", i, ".tif"))
      bbtl_rast[bbtl_rast > 0] <- 1

      bbtl_sum <- bbtl_rast
      
    }else{
      
      if(!file.exists(paste0(path, "/output_sim_eu_", i, "/bb_out/bbgrid_", y, "_sim_eu_", i, ".tif"))){
        print(paste0("sim ", i, " incomplete!"))
        next}
      # 
      bbtl_rast <- rast(paste0(path, "/output_sim_eu_", i, "/bb_out/bbgrid_", y, "_sim_eu_", i, ".tif"))
      bbtl_rast[bbtl_rast < (y-10)] <- 0
      bbtl_rast[bbtl_rast > 0] <- 1

      bbtl_sum <- c(bbtl_sum, bbtl_rast)
      
      if(y %in% c(40, 80)){
        bbtl_sum_out <- app(bbtl_sum, "sum", cores = 4)
        writeRaster(bbtl_sum_out, paste0(path_out, "bbtl/bbtl_", y, "_", rcp, "_", gcm, "_", rep, "_summed.tif"), overwrite = T)
        rm(bbtl_sum_out)
      }
      
      gc()
      
    }
    
    writeRaster(bbtl_rast, paste0(path_out, "bbtl/bbtl_", y, "_", rcp, "_", gcm, "_", rep, "_10year.tif"), overwrite = T)
    
    rm(bbtl_rast)
    for(x in 1:5){gc()}
  }
  
  rm(bbtl_sum)
  tmpFiles(remove = T)
  for(x in 1:5){gc()}
  
}




## management --
sims <- c(1:90)

for(i in sims){
  
  print(i)
  if(!file.exists(paste0(path, "/output_sim_eu_", i, "/mgmt_out/mgmtgrid_10_sim_eu_", i, ".tif"))){next}
  
  rcp <- master_tab[i, "RCP"]
  gcm <- master_tab[i, "GCM"]
  rep <- master_tab[i, "Rep"]
  if(file.exists(paste0(path_out, "mgmt/mgmt_80_", rcp, "_", gcm, "_", rep, "_summed.tif"))){next}
  
  for(y in seq(10, 80, 10)){

    if(y == 10){
      mgmt_rast <- rast(paste0(path, "/output_sim_eu_", i, "/mgmt_out/mgmtgrid_", y, "_sim_eu_", i, ".tif"))
      mgmt_rast[mgmt_rast > 0] <- 1

      mgmt_sum <- mgmt_rast
      
    }else{
      
      if(!file.exists(paste0(path, "/output_sim_eu_", i, "/mgmt_out/mgmtgrid_", y, "_sim_eu_", i, ".tif"))){
        print(paste0("sim ", i, " incomplete!"))
        next}

      mgmt_rast <- rast(paste0(path, "/output_sim_eu_", i, "/mgmt_out/mgmtgrid_", y, "_sim_eu_", i, ".tif"))
      mgmt_rast[mgmt_rast < (y-10)] <- 0
      mgmt_rast[mgmt_rast > 0] <- 1

      mgmt_sum <- c(mgmt_sum, mgmt_rast)

      if(y %in% c(40, 80)){
        mgmt_sum_out <- app(mgmt_sum, "sum", cores = 4)
        writeRaster(mgmt_sum_out, paste0(path_out, "mgmt/mgmt_", y, "_", rcp, "_", gcm, "_", rep, "_summed.tif"), overwrite = T)
        rm(mgmt_sum_out)
      }
      
      for(x in 1:5){gc()}
      
    }
    
    writeRaster(mgmt_rast, paste0(path_out, "mgmt/mgmt_", y, "_", rcp, "_", gcm, "_", rep, "_10year.tif"), overwrite = T)
    
    rm(mgmt_rast)
    for(x in 1:5){gc()}
  }
  
  rm(mgmt_sum)
  
  tmpFiles(remove = T)
  for(x in 1:5){gc()}
  
}





##################################################################################
### now historical simulations ---------------------------------------------------
##################################################################################

sims <- c(91:120)

## wind --

for(i in sims){
  
  print(i)
  if(!file.exists(paste0(path, "/output_sim_eu_", i, "/ccwind/wind_year_10_sim_eu_", i, ".tif"))){next}
  
  rcp <- master_tab[i, "RCP"]
  gcm <- master_tab[i, "GCM"]
  rep <- master_tab[i, "Rep"]
  
  if(file.exists(paste0(path_out, "wind/wind_80_", rcp, "_", gcm, "_", rep, "_summed.tif"))){next}
  
  for(y in seq(10, 80, 10)){
    
    if(y == 10){
      wind_rast <- rast(paste0(path, "/output_sim_eu_", i, "/ccwind/wind_year_", y, "_sim_eu_", i, ".tif"))
      wind_rast_focal <- wind_rast
    }else{
      
      if(!file.exists(paste0(path, "/output_sim_eu_", i, "/ccwind/wind_year_", y, "_sim_eu_", i, ".tif"))){
        print(paste0("sim ", i, " incomplete!"))
        next}
      r <- rast(paste0(path, "/output_sim_eu_", i, "/ccwind/wind_year_", y, "_sim_eu_", i, ".tif"))
      
      # in some years we want to know the total
      if(y %in% c(40, 80)){
        writeRaster(r, paste0(path_out, "wind/wind_", y, "_", rcp, "_", gcm, "_", rep, "_summed.tif"), overwrite = T)
      }
      
      wind_rast_focal <- r - wind_rast
      wind_rast <- r
      
    }

    writeRaster(wind_rast_focal, paste0(path_out, "wind/wind_", y, "_", rcp, "_", gcm, "_", rep, "_10year.tif"), overwrite = T)

    gc()
  }
  tmpFiles(remove = T)
  gc()
}



## fire --
for(i in sims){
  
  print(i)
  if(!file.exists(paste0(path, "/output_sim_eu_", i, "/ccfire/firegrid_10_sim_eu_", i, ".tif"))){next}
  
  rcp <- master_tab[i, "RCP"]
  gcm <- master_tab[i, "GCM"]
  rep <- master_tab[i, "Rep"]
  
  if(file.exists(paste0(path_out, "fire/fire_80_", rcp, "_", gcm, "_", rep, "_summed.tif"))){next}
  
  
  
  for(y in seq(10, 80, 10)){
    
    if(y == 10){
      fire_rast <- rast(paste0(path, "/output_sim_eu_", i, "/ccfire/firegrid_", y, "_sim_eu_", i, ".tif"))
      fire_rast_focal <- fire_rast
    }else{
      
      if(!file.exists(paste0(path, "/output_sim_eu_", i, "/ccfire/firegrid_", y, "_sim_eu_", i, ".tif"))){
        print(paste0("sim ", i, " incomplete!"))
        next}
      r <- rast(paste0(path, "/output_sim_eu_", i, "/ccfire/firegrid_", y, "_sim_eu_", i, ".tif"))
      
      
      if(y %in% c(40, 80)){
        writeRaster(r, paste0(path_out, "fire/fire_", y, "_", rcp, "_", gcm, "_", rep, "_summed.tif"), overwrite = T)
      }
      
      
      fire_rast_focal <- r - fire_rast
      fire_rast <- r
      
      gc()
    }
    
    writeRaster(fire_rast_focal, paste0(path_out, "fire/fire_", y, "_", rcp, "_", gcm, "_", rep, "_10year.tif"), overwrite = T)
    
    gc()
  }
  
  tmpFiles(remove = T)
  gc()
}




## bbtl --

sims <- c(91:120)

for(i in sims){
  
  print(i)
  if(!file.exists(paste0(path, "/output_sim_eu_", i, "/bb_out/bbgrid_10_sim_eu_", i, ".tif"))){next}
  
  rcp <- master_tab[i, "RCP"]
  gcm <- master_tab[i, "GCM"]
  rep <- master_tab[i, "Rep"]
  
  if(file.exists(paste0(path_out, "bbtl/bbtl_80_", rcp, "_", gcm, "_", rep, "_summed.tif"))){next}
  
  for(y in seq(10, 80, 10)){

    
    if(y == 10){
      bbtl_rast <- rast(paste0(path, "/output_sim_eu_", i, "/bb_out/bbgrid_", y, "_sim_eu_", i, ".tif"))
      bbtl_rast[bbtl_rast > 0] <- 1

      bbtl_sum <- bbtl_rast
      
    }else{
      
      if(!file.exists(paste0(path, "/output_sim_eu_", i, "/bb_out/bbgrid_", y, "_sim_eu_", i, ".tif"))){
        print(paste0("sim ", i, " incomplete!"))
        next}
       
      bbtl_rast <- rast(paste0(path, "/output_sim_eu_", i, "/bb_out/bbgrid_", y, "_sim_eu_", i, ".tif"))
      bbtl_rast[bbtl_rast < (y-10)] <- 0
      bbtl_rast[bbtl_rast > 0] <- 1
      
      bbtl_sum <- c(bbtl_sum, bbtl_rast)

      if(y %in% c(40, 80)){
        bbtl_sum_out <- app(bbtl_sum, "sum", cores = 4)
        writeRaster(bbtl_sum_out, paste0(path_out, "bbtl/bbtl_", y, "_", rcp, "_", gcm, "_", rep, "_summed.tif"), overwrite = T)
        rm(bbtl_sum_out)
      }
      
      
      gc()
      
    }
    
    writeRaster(bbtl_rast, paste0(path_out, "bbtl/bbtl_", y, "_", rcp, "_", gcm, "_", rep, "_10year.tif"), overwrite = T)
    
    rm(bbtl_rast)
    gc()
  }
  
  rm(bbtl_sum)
  tmpFiles(remove = T)
  gc()
  
}



## mgmt --

sims <- c(91:120)

for(i in sims){
  
  print(i)
  if(!file.exists(paste0(path, "/output_sim_eu_", i, "/mgmt_out/mgmtgrid_10_sim_eu_", i, ".tif"))){next}
  
  rcp <- master_tab[i, "RCP"]
  gcm <- master_tab[i, "GCM"]
  rep <- master_tab[i, "Rep"]
  
  if(file.exists(paste0(path_out, "mgmt/mgmt_80_", rcp, "_", gcm, "_", rep, "_summed.tif"))){next}
  
  for(y in seq(10, 80, 10)){
    
    # cat(print(paste0(" ", y, " ")))
    
    if(y == 10){
      mgmt_rast <- rast(paste0(path, "/output_sim_eu_", i, "/mgmt_out/mgmtgrid_", y, "_sim_eu_", i, ".tif"))
      mgmt_rast[mgmt_rast > 0] <- 1
      
      # mgmt_rast_focal <- mgmt_rast
      mgmt_sum <- mgmt_rast
      
    }else{
      
      if(!file.exists(paste0(path, "/output_sim_eu_", i, "/mgmt_out/mgmtgrid_", y, "_sim_eu_", i, ".tif"))){
        print(paste0("sim ", i, " incomplete!"))
        next}
      
      mgmt_rast <- rast(paste0(path, "/output_sim_eu_", i, "/mgmt_out/mgmtgrid_", y, "_sim_eu_", i, ".tif"))
      mgmt_rast[mgmt_rast < (y-10)] <- 0
      mgmt_rast[mgmt_rast > 0] <- 1
      
      mgmt_sum <- c(mgmt_sum, mgmt_rast)
      
      if(y %in% c(40, 80)){
        mgmt_sum_out <- app(mgmt_sum, "sum", cores = 4)
        writeRaster(mgmt_sum_out, paste0(path_out, "mgmt/mgmt_", y, "_", rcp, "_", gcm, "_", rep, "_summed.tif"), overwrite = T)
        rm(mgmt_sum_out)
      }
      
      
      gc()
      
    }
    
    writeRaster(mgmt_rast, paste0(path_out, "mgmt/mgmt_", y, "_", rcp, "_", gcm, "_", rep, "_10year.tif"), overwrite = T)
    
    rm(mgmt_rast)
    gc()
  }
  
  rm(mgmt_sum)
  
  tmpFiles(remove = T)
  for(g in 1:10){
    gc()}
  
}



