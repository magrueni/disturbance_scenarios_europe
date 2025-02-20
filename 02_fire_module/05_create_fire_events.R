# --------------------------------------------------------------
# Fire Module Scipt 05: Create series of fire events
# --------------------------------------------------------------
#
# This script combines the models from the previous steps to simulates fire events across historical and future 
# climate scenarios. It integrates in the first step the model for fire frequency per year. 
# In the second step, for each fire the location is sampled from the spatial probability map. 
# In the third step, depending on the location the fire size is modelled. 

# The output of this script is a fire event series where each fire for all simulation year is provided along its potential size.
# This is done for all different climate change scenarios as well as historical climate.
#
# Notes: the pre-calculated results of this script are located in ./fire_module/fire_event_series/
#-------------------------------------------------------------------------------




### libraries -----------------------------

library(raster)
library(terra)
library(sf)
library(dplyr)
library(tidyverse)
library(data.table)
library(VGAM)
library(parallel)
library(glmmTMB)
library(MetBrewer)


# define path
path <- "/.../"



### load data ----------------------------------------------------------------------------

# load some supporting data
forest_lookup <- readRDS(paste0(path, "fire_module/gis/forest_lookup_grids.rds"))
gridID_rast <- rast(paste0(path, "/reference_grids/gridID_10k.tif"))

# load probability map and plot
eu <- vect(paste0(path, "/gis/europe_lowres.shp"))
proba_rast <- rast(paste0(path, "/fire_module/models/probability_raster_baseline_cordex_biascorrected_10k.tif"))
cols <- MetBrewer::met.brewer("Homer2", 100, "continuous", direction = -1)
par(oma = c(2,2,2,2))
plot(eu, col = "black", box = F, axes = F)
plot(proba_rast, col = cols, box = F, axes = F, add = T)

# dev.print(png, paste0(path, "fire_module/figures/fire_probability_map.png"), width = 7, height = 7, units = "in", res = 300)


# load frequency model
mod_freq <- readRDS(paste0(path, "/fire_module/models/final_freq_lmer_100km.rds"))
mod_size <- readRDS(paste0(path, "/fire_module/models/final_size_lmer_100km.rds"))


# load other data
proj_leae <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs" 
ref_grid_100km <- rast(paste0(path, "/reference_grids/reference_grid_100km.tif"))
ref_grid_10km <- rast(paste0(path, "/reference_grids/reference_grid_10km.tif"))
df_100km <- as.data.frame(ref_grid_100km, xy = T)
df_10km <- as.data.frame(ref_grid_10km, xy = T)

#wind table
wind_table <- read.csv(paste0(path, "/fire_module/wind_table.csv"))



### processing ---------------------------------------------------------------------------
# run the pipeline to create series of fire events

# difine how many simulations
nsim <- c(1:10)


for(n in nsim){
  
  
  # for the three rcps and historical
  rcps <- c("historical", "rcp85", "rcp45", "rcp26") 
  
  
  for(c in rcps){

    # three GCMs per RCP
    gcms <- c("ichec", "mpi", "ncc")
    
    
    for(g in gcms){
      
      
      r <- 0
      output <- list()
      
      # for historical runs we sample years from 1986:2020. For future runs we use 2020:2100
      if(c == "historical"){
        
        set.seed(42)
        years <- sample(c(1986:2020), 82, replace = T)
        
      }else{
        
        years <- 2020:2100
        
      }
      
      # loop over the years
      
      # we provide here the example year 2026 (RCP2.6, ICHEC scenario) and the year 2093 (RCP8.5, ICHEC scenario) for code testing
      for(y in years){
        
        r <- r + 1
        
      
        ### in each year, we need to estimate the number of fires
        
        # load the climate data
        if(y %in% 1981:2020){
          clim <- read.csv(paste0(path, "/fire_module/climate/hist_vpd_100km.csv"), header = T)
        }else{
          clim <- fread(paste0(path, "/fire_module/climate/fut_data_vpd_100km_", g, "_", c, ".csv"), sep =",")
        }
        
        colnames(clim) <- gsub("X", "", colnames(clim))
        
        clim_year <- clim %>% dplyr::select(paste0(y), gridid) %>% drop_na() %>% 
          rename("vpd_summer_mean_rollmax" = paste0(y)) %>% 
          group_by(gridid) %>% 
          summarise(vpd_summer_mean_rollmax = mean(vpd_summer_mean_rollmax)) %>% 
          mutate(year = paste0(y)) %>% 
          dplyr::select(year, vpd_summer_mean_rollmax, gridid)

        
        
        ## 1. make prediction for the number of fires
        
        lambda <- predict(mod_freq, newdata = clim_year, type = "conditional",  allow.new.levels = T) 
        p <- predict(mod_freq, newdata = clim_year, type="zprob",  allow.new.levels = T)
        pred_new <- emdbook::rzinbinom(nrow(clim_year),
                                       mu=lambda, size=sigma(mod_freq),
                                       zprob=p)
        
        pred_tiles <- pred_new %>%
          as_tibble() %>%
          mutate(gridid = clim_year$gridid) %>%
          group_by(gridid) %>% 
          summarize(n_fires = round(sum(value, na.rm = T), 0))
        
        nfires <- sum(pred_tiles$n_fires)
        
        cat(paste0("In year ", y, " ", nfires, " fires occurred! \n "))
        
        
        ## 2. now sample the number of fires from the probability map
        
        # load raster depending on the year
        t <- ifelse(y %in% c(2011:2040), "2011-2040", 
                    ifelse(y %in% c(2041:2070), "2041-2070", "2071-2100"))
        
        if(y %in% c(1981:2020)){
          proba_rast <- rast(paste0(path, "/fire_module/models/probability_raster_baseline_cordex_biascorrected_10k.tif"))
          crs(proba_rast) <- proj_leae
        }else{
          proba_rast <- rast(paste0(path, "/fire_module/models/probability_rasters/probability_rast_fut_", t, "_", c, "_biascorrected_10k.tif"))
        }
        
        proba_rast <- terra::project(proba_rast, proj_leae)
        
        
        # sample locations from the fire probability map
        sampled_locations <- list()
        empty <- c()
        
        # loop over 100km tiles
        for (i in 1:nrow(pred_tiles)) {
          
          if(as.numeric(pred_tiles[i, "n_fires"]) == 0){next}
          
          # get the focal grid
          focal_grid <- df_100km[i, "gridid"]
          
          # sample points according to frequency prediction
          n_fires_tile <- as.numeric(pred_tiles[pred_tiles$gridid == focal_grid, "n_fires"])
          
          region <- ref_grid_100km
          region[region != focal_grid] <- NA
          region <- terra::project(region, proba_rast)
          region_mask <- terra::mask(proba_rast, region)
          proba_df <- as.data.frame(region_mask, xy = T)
          
          if(nrow(proba_df) == 0){
            
            # Find the neighboring cells
            neighbors <- as.data.frame(adjacent(ref_grid_100km, cells = focal_grid, pairs=TRUE))
            region <- ref_grid_100km
            region[!region %in% neighbors$to] <- NA
            region <- terra::project(region, proba_rast)
            region_mask <- terra::mask(proba_rast, region)
            proba_df <- as.data.frame(region_mask, xy = T)
            
          }
          
          if(nrow(proba_df) == 0){
            
            print(paste0("no grids ", n_fires_tile, " fires missing."))
            empty <- c(empty, focal_grid)
            next
            
          }
          
          
          # Sample points from the region based on the raster values
          sampled <- proba_df[sample(nrow(proba_df), size = n_fires_tile, replace = TRUE, prob = proba_df$prob),] %>% 
            mutate(gridid =  focal_grid)
          
          # Add the sampled points to the overall collection
          sampled_locations[[i]] <- sampled
        }
        
        
        # add grid ID and convert to list again
        sampled_locations <- do.call(rbind, sampled_locations)
        sampled_locations <- cbind(sampled_locations, gridID_10km = terra::extract(ref_grid_10km, sampled_locations[, c("x", "y")])[,2])
        sampled_locations_list <- split(sampled_locations, seq(nrow(sampled_locations)))
        
        # filter lookup table with grid ids in the sampled locations
        lookup_filtered <- forest_lookup[forest_lookup$gridID %in% sampled_locations$gridID_10km,]
        
        # sampled forest grid cell
        subsample_forestgrid <- function(x, lookup){
          
          grid <- as.numeric(x["gridID_10km"])
          forestgrids <- lookup[lookup$gridID == grid,]
          
          if(nrow(forestgrids) == 0){
            
            neighbors <- as.data.frame(adjacent(ref_grid_10km, cells = grid, pairs=TRUE))
            closest_index <- sample(neighbors$to, 1)
            grid <- lookup$gridID[closest_index]
            forestgrids <- lookup[lookup$gridID == grid,]
          }
          
          sampled_grid <- forestgrids[sample(nrow(forestgrids), size = 1, replace = FALSE),]
          
          x_new <- c(x, sampled_grid)
          return(x_new)
          
        }
        
        # sample 1 forestgrid cell within the 10km cell randomly
        sampled_locations_new <- parallel::mclapply(sampled_locations_list, subsample_forestgrid, lookup = lookup_filtered,
                                                    mc.cleanup = TRUE, mc.cores = 24)
        sampled_locations_df <- do.call(rbind, sampled_locations_new)
        colnames(sampled_locations_df) <- c("x_10k", "y_10k", "prob", "gridid_100km", "gridid_10km", "forestgridID", "x", "y", "gridID2")
        sampled_locations_df <- sampled_locations_df %>% data.table() %>% 
          dplyr::select(x, y, gridid_100km) %>% 
          mutate(x = as.numeric(x),
                 y = as.numeric(y),
                 gridid_100km = as.numeric(gridid_100km)) %>% 
          drop_na()
        
        df_full <- na.omit(sampled_locations_df)
        
        dim(df_full)
        if(nfires != dim(df_full)[1]){
          print(paste0(nfires - dim(df_full)[1], " fires lost!"))
        }
        
        ### finally get the max size for each fire -----------------------------------------------------------------------
        # get year before and after
        c2 <- ifelse(c == "historical" & y > 2005, "rcp45", c)
        
        if(c2 == "historical"){
          
          if(y == 1986){year_prev <- y}else{year_prev <- y - 1}
          if(y == 2005){year_aftr <- y}else{year_aftr <- y + 1}
          
        }else{
          
          if(y == 2006){year_prev <- y}else{year_prev <- y - 1}
          if(y == 2100){year_aftr <- y}else{year_aftr <- y + 1}
          
        }
        
        # load the data
        summer_vpd <- rast(paste0(path, "/fire_module/climate/vpd_summer_", y,"_", g, "_", c2, "_leae.tif"))
        
        # load data of previous year
        prev_summer_vpd <- rast(paste0(path, "/fire_module/climate/vpd_summer_", year_prev, "_", g, "_", c2, "_leae.tif"))
        
        # load data of following year
        aftr_summer_vpd <- rast(paste0(path, "/fire_module/climate/vpd_summer_", year_aftr, "_", g, "_", c2, "_leae.tif"))
        
        # calculate rolling max for each grid
        rolling_max_vpd_summer <- max(summer_vpd, prev_summer_vpd, aftr_summer_vpd)
        
        
        size_pred <- cbind(df_full, vpd_summer_mean_rollmax = terra::extract(rolling_max_vpd_summer, df_full[, c("x", "y")])[,2]) %>% 
          rename("gridid" =  gridid_100km)
        
        # make prediction
        pred <- predict(mod_size, newdata = size_pred, allow.new.levels = T, re.form = ~(1 + log10(vpd_summer_mean_rollmax) | gridid))
        pred <- 10^pred %>%
          as_tibble()
        
        df_fin <- cbind(fire_id = c(1:nrow(size_pred)), year = rep(r, nrow(size_pred)), size_pred, max_size = pred$value * 0.0001)
        df_fin <- df_fin %>% 
          left_join(wind_table, by = "gridid") %>%
          rowwise() %>%
          mutate(
            winddirection = rnorm(1, mean_dir, sd_dir),
            windspeed = rnorm(1, mean_speed, sd_speed)
          ) %>%
          ungroup()
        
        
        output[[r]] <- df_fin %>% dplyr::select(fire_id, year, x, y, max_size, winddirection, windspeed)
      }
      
      
      all_outputs <- do.call(rbind, output)
      
      
      if(c == "historical"){
        write.csv(all_outputs, paste0(path, "/fire_module/fire_event_series/future_fire_simulation_", g, "_", c, "_simulation_", n, "_long.csv"), row.names = FALSE)
      }else{
        write.csv(all_outputs, paste0(path, "/fire_module/fire_event_series/future_fire_simulation_", g, "_", c, "_simulation_", n, ".csv"), row.names = FALSE)
      }
      
    }
    
  }
  
}


### end ----------------------------------------------------