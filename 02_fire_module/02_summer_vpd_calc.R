
# --------------------------------------------------------------
# Fire Module Scipt 02: vapor pressure deficit calculation for Europe
# --------------------------------------------------------------
# Author: Marc Grünig
# Date: 2025-01-07
# 
# This script calculates the Vapor Pressure Deficit (VPD) from temperature and humidity data for different climate scenarios and future projections. 
# The script utilizes climate model data from Euro-CORDEX and CMIP5 GCMs, including historical and future periods, to calculate VPD values.
# The resulting VPD data is bias-corrected using ERA5 reanalysis data and then aggregated to a 100km resolution grid for further analysis. 
#
# The script processes and aggregates data from multiple GCMs and RCP scenarios, applies bias correction for VPD values, and stores the results 
# as raster files for spatial analysis. Additionally, a rolling maximum is applied to extract peak summer VPD values for each year.
# 
# --------------------------------------------------------------


# Packages ----------------------------------------------------------------

library(tidyverse)
library(raster)
library(lubridate)
library(sf)
library(patchwork)
library(terra)
library(ncdf4)
library(weathermetrics)


path <- "/.../"


# Studyregion -------------------------------------------------------------
studyregion <- read_sf(paste0(path, "/fire_module/climate/climategrid_epsg3035.gpkg"))
studyregion_latlng <- st_transform(studyregion, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

cntrs <- read_sf(paste0(path, "/gis/countries/europe.shp"))

# define projections
proj_wgs84 <- "+proj=longlat +datum=WGS84 +no_defs"
proj_nc <- "+proj=ob_tran +o_proj=longlat +o_lon_p=-162 +o_lat_p=39.25 +lon_0=180"
proj_leae <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs" 


# Functions ---------------------------------------------------------------
vpd_calc <- function(t, td, c1 = 0.611, c2 = 17.67, c3 = 243.5) {
  c1 * exp((c2 * (t)) / (t + c3)) - c1 * exp((c2 * (td)) / (td + c3))
}


# Load ERA5 soil moisture data --------------------------------------------
gcms <- c("ncc", "mpi", "ichec")
rcps <- c("rcp26", "rcp45", "rcp85")
rcms <- c("smhi_rca4")


# Loop through RCPS and GCMS 
for (rcp in rcps) {
  for (gcm in gcms) {
    print(gcm)
    
    # Set GCM details
    gcm_details <- list(
      ichec = list(name = "ICHEC-EC-EARTH", run = "r12i1p1", version = "v1"),
      ncc = list(name = "NCC-NorESM1-M", run = "r1i1p1", version = "v1"),
      mpi = list(name = "MPI-M-MPI-ESM-LR", run = "r1i1p1", version = "v1a")
    )
    gcm.long <- gcm_details[[gcm]]$name
    run <- gcm_details[[gcm]]$run
    version <- gcm_details[[gcm]]$version
    
    timesteps <- if (rcp == "historical") {
      c(1981, 1991, 2001)
    } else {
      c(2006, 2011, 2021, 2031, 2041, 2051, 2061, 2071, 2081, 2091)
    }
    
    # Loop over time steps ----------------------------------------------
    for (r in seq_along(timesteps)) {
      timestep <- timesteps[r]
      timestep_long <- if (rcp == "historical") {
        c("198101-199012", "199101-200012", "200101-200512")
      } else {
        c("200601-201012", "201101-202012", "202101-203012", "203101-204012", "204101-205012", "205101-206012", "206101-207012", "207101-208012", "208101-209012", "209101-210012")
      }
      
      # Load climate data ------------------------------------------------
      tas_ras <- brick(paste0(path, "/climatedata/euro-cordex/monthly/", rcp, "/tas_EUR-11_", gcm.long, "_", rcp, "_", run, "_SMHI-RCA4_", version, "_mon_", timestep_long[r], ".nc"))
      hurs_ras <- brick(paste0(path, "/climatedata/euro-cordex/monthly/", rcp, "/hurs_EUR-11_", gcm.long, "_", rcp, "_", run, "_SMHI-RCA4_", version, "_mon_", timestep_long[r], ".nc"))
      
      # Extract summer months (June, July, August)
      tas_ras_summer <- stack(tas_ras[[grep(".06.16", names(tas_ras))]], tas_ras[[grep(".07.16", names(tas_ras))]], tas_ras[[grep(".08.16", names(tas_ras))]])
      hurs_ras_summer <- stack(hurs_ras[[grep(".06.16", names(hurs_ras))]], hurs_ras[[grep(".07.16", names(hurs_ras))]], hurs_ras[[grep(".08.16", names(hurs_ras))]])
      
      # Process each year ------------------------------------------------
      end_year <- ifelse(timestep_long[r] %in% c("200601-201012", "200101-200512"), 4, 9)
      for (y in 0:end_year) {
        year <- timestep + y
        if (!any(grepl(paste0(year), names(tas_ras_summer)))) next
        
        # Process climate variables for the year
        summer_tas <- rast(mean(tas_ras_summer[[grepl(paste0(year), names(tas_ras_summer))]])) - 273.15
        crs(summer_tas) <- proj_wgs84
        summer_tas_proj <- terra::project(summer_tas, proj_nc)
        
        # Repeat for humidity ----------------
        summer_hurs <- rast(mean(hurs_ras_summer[[grepl(paste0(year), names(hurs_ras_summer))]]))
        crs(summer_hurs) <- proj_wgs84
        summer_hurs_proj <- terra::project(summer_hurs, proj_nc)
        
        # Dewpoint and VPD calculation
        td <- humidity.to.dewpoint(rh = summer_hurs_proj, t = summer_tas_proj, temperature.metric = "celsius")
        vpd_method2 <- vpd_calc(summer_tas_proj, td)
        
        # Save results ----------------------------------------------------
        terra::writeRaster(vpd_method2, paste0(path, "/fire_module/climate/vpd_summer_", year, "_", gcm, "_", rcp, "_leae.tif"), overwrite = TRUE)
      }
    }
  }
}




# Bias correction using ERA5 data --------------------------------------

# Define baseline years range
years <- 1986:2020

# Initialize stacks for diff, era, and cordex
diff_stack <- raster::stack()
era_stack <- raster::stack()
cordex_stack <- raster::stack()

# Loop through each year for bias correction
for (y in years) {
  # Load ERA5 and Cordex data
  era5 <- rast(paste0(path, "/fire_module/climate/era5_summer_vpd_", y, ".tif"))
  cord <- rast(paste0(path, "/fire_module/climate/vpd_summer_", y, "_ensemble_historical_leae.tif"))
  
  # Project ERA5 to match Cordex
  era5 <- terra::project(era5, cord)
  era5 <- terra::mask(era5, cord)
  
  # Stack the data
  era_stack <- stack(era_stack, raster(era5))
  cordex_stack <- stack(cordex_stack, raster(cord))
  
  # Calculate difference between ERA5 and Cordex
  dif <- era5 - cord
  diff_stack <- raster::stack(diff_stack, raster(dif))
}

# Calculate and plot bias
bias <- mean(diff_stack)
bias_leae <- terra::project(rast(bias), proj_leae)
plot(bias_leae)

# Calculate and plot mean of ERA5 and Cordex
mean_era <- mean(era_stack)
mean_era <- terra::project(rast(mean_era), proj_leae)

mean_cor <- mean(cordex_stack)
mean_cor <- terra::project(rast(mean_cor), proj_leae)

# check visually
# plot(mean_era, main = "ERA5")
# plot(mean_cor, main = "Cordex")

# Calculate and plot the difference
diff2 <- mean_era - mean_cor
plot(diff2)

# Subtract bias from the difference
diff_diff <- diff2 - bias_leae
plot(diff_diff)


# Bias correction for historical ensembles
for (y in years) {
  rast <- rast(paste0(path, "/fire_module/climate/vpd_summer_", y, "_ensemble_historical_leae.tif"))
  
  # Project, crop, and mask the raster
  rast <- terra::project(rast, bias_leae)
  rast <- terra::crop(rast, bias_leae)
  rast <- terra::mask(rast, bias_leae)
  
  # Apply bias correction
  rast_new <- rast + bias_leae
  
  # Save the corrected raster
  terra::writeRaster(rast_new, paste0(path, "/fire_module/climate/vpd_summer_", y, "_ensemble_historical_leae_v2_biascorrected.tif"), overwrite = TRUE)
}


# Bias correction for future data
gcms <- c("mpi", "ichec", "ncc")
rcps <- c("historical", "rcp26", "rcp45", "rcp85")

for (r in rcps) {
  for (g in gcms) {
   
    # Set years for each RCP
    years <- ifelse(r == "historical", 1986:2005, 2006:2100)
    
    for (y in years) {
      # Check if the file exists, if not, create the bias
      file_path <- paste0(path, "/fire_module/climate/vpd_summer_", y, "_", g, "_", r, "_leae.tif")
      rast <- rast(file_path)
      
      # Project, crop, and mask the raster, then apply bias correction
      rast <- terra::project(rast, bias_leae)
      rast <- terra::crop(rast, bias_leae)
      rast <- terra::mask(rast, bias_leae)
      rast_new <- rast + bias_leae
      
      # Save the corrected raster
      terra::writeRaster(rast_new, paste0(path, "/fire_module/climate/summer_vpd/vpd_summer_", y, "_", g, "_", r, "_leae_v2_biascorrected.tif"), overwrite = TRUE)
    }
  }
}



# extract VPD for the 100km resolution --------------------------------------------------------------------

# first for historical climate

r <- 0
vpd_stack <- ref_grid_100km  # Initialize the stack for VPD data

# Loop through each year from 1986 to 2020 to extract VPD data
for(y in 1986:2020) {
  print(y)
  
  r <- r + 1  # Increment the counter for each year
  
  # Get the previous and following years to calculate rolling max
  year_prev <- ifelse(y == 1986, y, y - 1)
  year_aftr <- ifelse(y == 2020, y, y + 1)
  
  # Load data for the current, previous, and following years
  summer_vpd <- rast(paste0(path, "/fire_module/climate/summer_vpd/vpd_summer_", y, "_ensemble_historical_leae_v2_biascorrected.tif"))
  prev_summer_vpd <- rast(paste0(path, "/fire_module/climate/summer_vpd/vpd_summer_", year_prev, "_ensemble_historical_leae_v2_biascorrected.tif"))
  aftr_summer_vpd <- rast(paste0(path, "/fire_module/climate/summer_vpd/vpd_summer_", year_aftr, "_ensemble_historical_leae_v2_biascorrected.tif"))
  
  # Calculate the rolling maximum of the VPD data across the three years
  rolling_max_vpd_summer <- max(summer_vpd, prev_summer_vpd, aftr_summer_vpd)
  
  # Resample to 10km grid and aggregate to 100km grid
  vpd_10km <- terra::resample(rolling_max_vpd_summer, ref_grid_10km)
  vpd_100km <- terra::aggregate(vpd_10km, fact = 10, fun = "mean", na.rm = T)
  vpd_100km <- terra::resample(vpd_100km, ref_grid_100km)
  
  # Apply focal operation to smooth the 100km grid
  vpd_100km <- focal(vpd_100km, c(3, 3), "mean", na.policy = "only")
  vpd_100km <- terra::mask(vpd_100km, ref_grid_100km)  # Mask out regions outside the reference grid
  names(vpd_100km) <- y  # Assign the year as the name
  
  # Add the processed data to the VPD stack
  vpd_stack <- c(vpd_stack, vpd_100km)
}

# Convert the VPD stack to a data frame and save it as a CSV
vpd_df <- as.data.frame(vpd_stack)
write.csv(vpd_df, paste0(path, "/fire_module/climate/hist_vpd_100km.csv"), row.names = F)


# then also for future climate

# List of GCMs (Global Climate Models) and RCPs (Representative Concentration Pathways)
gcms <- c("ichec", "mpi", "ncc")
rcps <- c("rcp85", "rcp45", "rcp26")

# Loop through each GCM and RCP to extract future VPD data
for(g in gcms) {
  for(rcp in rcps) {
    
    r <- 0  # Reset the counter for each GCM and RCP
    
    vpd_stack <- ref_grid_100km  # Initialize the stack for future VPD data
    
    for(y in 2020:2100) {
      r <- r + 1  # Increment the counter for each year
      
      # Get the previous and following years to calculate rolling max
      year_prev <- ifelse(y == 2020, y, y - 1)
      year_aftr <- ifelse(y == 2100, y, y + 1)
      
      # Load data for the current, previous, and following years for future scenarios
      summer_vpd <- rast(paste0(path, "/fire_module/climate/summer_vpd/vpd_summer_", y, "_", g, "_", rcp, "_leae_v2_biascorrected.tif"))
      prev_summer_vpd <- rast(paste0(path, "/fire_module/climate/summer_vpd/vpd_summer_", year_prev, "_", g, "_", rcp, "_leae_v2_biascorrected.tif"))
      aftr_summer_vpd <- rast(paste0(path, "/fire_module/climate/summer_vpd/vpd_summer_", year_aftr, "_", g, "_", rcp, "_leae_v2_biascorrected.tif"))
      
      # Calculate the rolling maximum of the VPD data
      rolling_max_vpd_summer <- max(summer_vpd, prev_summer_vpd, aftr_summer_vpd)
      
      # Resample and aggregate to 100km grid as done for historical data
      vpd_10km <- terra::resample(rolling_max_vpd_summer, ref_grid_10km)
      vpd_100km <- terra::aggregate(vpd_10km, fact = 10, fun = "mean", na.rm = T)
      vpd_100km <- terra::resample(vpd_100km, ref_grid_100km)
      
      # Apply focal operation to smooth the 100km grid
      vpd_100km <- focal(vpd_100km, c(3, 3), "mean", na.policy = "only")
      vpd_100km <- terra::mask(vpd_100km, ref_grid_100km)  # Mask out regions outside the reference grid
      names(vpd_100km) <- y  # Assign the year as the name
      
      # Add the processed data to the VPD stack
      vpd_stack <- c(vpd_stack, vpd_100km)
    }
    
    # Convert the VPD stack to a data frame and save it as a CSV
    vpd_df <- as.data.frame(vpd_stack)
    write.csv(vpd_df, paste0(path, "/fire_module/climate/fut_data_allgrids", g, "_", rcp, ".csv"), row.names = F)
    
  }
}




