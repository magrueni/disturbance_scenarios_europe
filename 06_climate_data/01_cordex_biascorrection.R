# --------------------------------------------------------------
# Climate Data 01: Biascorrection
# --------------------------------------------------------------
# Author: Marc Grünig
# Date: 2025-01-07
# 
# Description:
# In this script we calculate the biascorrection factor between 
# ERA-Interim data and the EURO-CORDEX data for each grid cell. This correction
# is then applied in the next script.
# 
# 
# Notes:
# - The script relies on data from ERA5 and Euro-Cordex which we cannot provide here.
# - The resulting dataframes are stored in the folder.
# 
# --------------------------------------------------------------

library(dplyr)
library(terra)
library(raster)
library(sf)
library(rgeos)
library(viridis)
library(ncdf4)

terraOptions(tempdir = "./tmp")

# define file paths
path <- "/.../"
dat_path <- paste0(path, "/climatedata/euro-cordex/historical/")
eval_path <- paste0(path, "/climatedata/euro-cordex/evaluation/")

# define CRS
proj_nc <- "+proj=ob_tran +o_proj=longlat +o_lon_p=-162 +o_lat_p=39.25 +lon_0=180 +ellps=WGS84"
proj <- "+proj=longlat +datum=WGS84 +no_defs"
proj_new <- "+proj=tmerc +lat_0=0 +lon_0=10.3333333333333 +k=1 +x_0=0 +y_0=-5000000 +ellps=bessel +units=m +no_defs +type=crs"

# get lists of all files
all_ncs_tas <- list.files(paste0(dat_path, "./"), patter = "tas_/*", full.names = TRUE)
all_ncs_tasmin <- list.files(paste0(dat_path, "./"), patter = "tasmin_/*", full.names = TRUE)
all_ncs_tasmax <- list.files(paste0(dat_path, "./"), patter = "tasmax_/*", full.names = TRUE)
all_ncs_pr <- list.files(paste0(dat_path, "./"), patter = "pr_/*", full.names = TRUE)
all_ncs_rh <- list.files(paste0(dat_path, "./"), patter = "hurs_/*", full.names = TRUE)
all_ncs_rad <- list.files(paste0(dat_path, "./"), patter = "rsds_/*", full.names = TRUE)

eval_ncs_tas <- list.files(paste0(eval_path, "./"), patter = "tas_/*", full.names = TRUE)
eval_ncs_tasmin <- list.files(paste0(eval_path, "./"), patter = "tasmin_/*", full.names = TRUE)
eval_ncs_tasmax <- list.files(paste0(eval_path, "./"), patter = "tasmax_/*", full.names = TRUE)
eval_ncs_pr <- list.files(paste0(eval_path, "./"), patter = "pr_/*", full.names = TRUE)
eval_ncs_rh <- list.files(paste0(eval_path, "./"), patter = "hurs_/*", full.names = TRUE)
eval_ncs_rad <- list.files(paste0(eval_path, "./"), patter = "rsds_/*", full.names = TRUE)

# functions
vpd_calc <- function(t, td, c1 = 0.611, c2 = 17.67, c3 = 243.5) {
  c1 * exp((c2 * t) / (t + c3)) - c1 * exp((c2 * td) / (td + c3))
}


# define timesteps
timesteps <- c("1981","1986","1991","1996","2001")
schaltjahr <- c(1984,1988,1992,1996,2000,2004)

# define gcms
gcms <- c("ICHEC-EC-EARTH", "MPI-M-MPI-ESM-LR", "NCC-NorESM1-M")


### start processing ----------------------------------------------------------

for(g in gcms){
  
  #create empty stacks for output
  tas_stack <- rast()
  tasmin_stack <- rast()
  tasmax_stack <- rast()
  pr_stack <- rast()
  rh_stack <- rast()
  rad_stack <- rast()
  vpd_stack <- rast()
  
  # get files of focal gcm
  dat_tas_layers_gcm <- all_ncs_tas[which(grepl(paste(g), all_ncs_tas))]
  dat_tasmin_layers_gcm <- all_ncs_tasmin[which(grepl(paste(g),all_ncs_tasmin))]
  dat_tasmax_layers_gcm <- all_ncs_tasmax[which(grepl(paste(g),all_ncs_tasmax))]
  dat_pr_layers_gcm <- all_ncs_pr[which(grepl(paste(g),all_ncs_pr))]
  dat_rh_layers_gcm <- all_ncs_rh[which(grepl(paste(g),all_ncs_rh))]
  dat_rad_layers_gcm <- all_ncs_rad[which(grepl(paste(g), all_ncs_rad))]
  
  # loop over timesteps
  for(i in timesteps){
    
    cat(i)
    
    # get files
    dat_tas_layers <- dat_tas_layers_gcm[which(grepl(paste(i),dat_tas_layers_gcm))]
    dat_tasmin_layers <- dat_tasmin_layers_gcm[which(grepl(paste(i),dat_tasmin_layers_gcm))]
    dat_tasmax_layers <- dat_tasmax_layers_gcm[which(grepl(paste(i),dat_tasmax_layers_gcm))]
    dat_pr_layers <- dat_pr_layers_gcm[which(grepl(paste(i),dat_pr_layers_gcm))]
    dat_rh_layers <- dat_rh_layers_gcm[which(grepl(paste(i),dat_rh_layers_gcm))]
    dat_rad_layers <- dat_rad_layers_gcm[which(grepl(paste(i),dat_rad_layers_gcm))]
    
    eval_tas_layers <- eval_ncs_tas[which(grepl(paste(i),eval_ncs_tas))]
    eval_tasmin_layers <- eval_ncs_tasmin[which(grepl(paste(i),eval_ncs_tasmin))]
    eval_tasmax_layers <- eval_ncs_tasmax[which(grepl(paste(i),eval_ncs_tasmax))]
    eval_pr_layers <- eval_ncs_pr[which(grepl(paste(i),eval_ncs_pr))]
    eval_rh_layers <- eval_ncs_rh[which(grepl(paste(i),eval_ncs_rh))]
    eval_rad_layers <- eval_ncs_rad[which(grepl(paste(i),eval_ncs_rad))]
    
    # convert to raster stacks
    dat_tas <- rast(dat_tas_layers)
    dat_tasmin <- rast(dat_tasmin_layers)
    dat_tasmax <- rast(dat_tasmax_layers)
    dat_pr <- rast(dat_pr_layers)
    dat_rh <- rast(dat_rh_layers)
    dat_rad <- rast(dat_rad_layers)
    
    eval_tas <- rast(eval_tas_layers)
    eval_tasmin <- rast(eval_tasmin_layers)
    eval_tasmax <- rast(eval_tasmax_layers)
    eval_pr <- rast(eval_pr_layers)
    eval_rh <- rast(eval_rh_layers)
    eval_rad <- rast(eval_rad_layers)
    
    years <- c(as.numeric(i),as.numeric(i)+1,as.numeric(i)+2,as.numeric(i)+3,as.numeric(i)+4)
    
    # calc annualmean temp
    days <- 0
    
    for(j in 1:5){
      
      print(years[j])
      
      # check whether year is a schaltjahr
      if(years[j] %in% schaltjahr){x <- 366}else{x <- 365}
      days_end <- days+x
      days_start <- days+1
      
      ## extract daily temps for the year
      dat_tas_extract <- app(dat_tas[[days_start:days_end]], fun = "mean")
      dat_tas_extract <- dat_tas_extract-273.15
      dat_tasmin_extract <- app(dat_tasmin[[days_start:days_end]], fun = "mean")
      dat_tasmax_extract <- app(dat_tasmax[[days_start:days_end]], fun = "mean")
      dat_pr_extract <- app(dat_pr[[days_start:days_end]], fun = "sum")
      dat_pr_extract <- dat_pr_extract*86400
      dat_rh_extract <- app(dat_rh[[days_start:days_end]], fun = "mean")
      dat_rad_extract <- app(dat_rad[[days_start:days_end]], fun = "sum")
      
      # calc vpd
      dat_td <- weathermetrics::humidity.to.dewpoint(rh = dat_rh_extract, t = dat_tas_extract, temperature.metric = "celsius")
      dat_vpd <- vpd_calc(t = dat_tas_extract, td = dat_td)
      
      #same for eval dat
      eval_tas_extract <- app(eval_tas[[days_start:days_end]], fun = "mean")
      eval_tas_extract <- eval_tas_extract-273.15
      eval_tasmin_extract <- app(eval_tasmin[[days_start:days_end]], fun = "mean")
      eval_tasmax_extract <- app(eval_tasmax[[days_start:days_end]], fun = "mean")
      eval_pr_extract <- app(eval_pr[[days_start:days_end]], fun = "sum")
      eval_pr_extract <- eval_pr_extract*86400
      eval_rh_extract <- app(eval_rh[[days_start:days_end]], fun = "mean")
      eval_rad_extract <- app(eval_rad[[days_start:days_end]], fun = "sum")
      
      td_eval <- weathermetrics::humidity.to.dewpoint(rh = eval_rh_extract, t = eval_tas_extract, temperature.metric = "celsius")
      eval_vpd <- vpd_calc(t = eval_tas_extract, td = td_eval)
      
      # get bias
      tas_out <- eval_tas_extract - dat_tas_extract 
      tasmin_out <- eval_tasmin_extract - dat_tasmin_extract 
      tasmax_out <- eval_tasmax_extract - dat_tasmax_extract
      pr_out <- eval_pr_extract/dat_pr_extract 
      rh_out <- eval_rh_extract - dat_rh_extract
      rad_out <- eval_rad_extract/dat_rad_extract 
      vpd_out <- eval_vpd/dat_vpd
      
      # add to output stack
      tas_stack <- c(tas_stack, tas_out)
      tasmin_stack <- c(tasmin_stack, tasmin_out)
      tasmax_stack <- c(tasmax_stack, tasmax_out)
      pr_stack <- c(pr_stack, pr_out)
      rh_stack <- c(rh_stack, rh_out)
      rad_stack <- c(rad_stack, rad_out)
      vpd_stack <- c(vpd_stack, vpd_out)
      
      
      days <- days_end
      
      # check
      if(j == 5){if(days != nlyr(dat_tas)){cat("ERROR: number of days not correct!", i)}}
      
    }
    
  } # close timestep
  
  # calculate mean over whole period
  tas_bias <- app(tas_stack, "mean")
  tasmin_bias <- app(tasmin_stack, "mean")
  tasmax_bias <- app(tasmax_stack, "mean")
  pr_bias <- app(pr_stack, "mean")
  rh_bias <- app(rh_stack, "mean")
  rad_bias <- app(rad_stack, "mean")
  vpd_bias <- app(vpd_stack, "mean")
  
  # save
  writeRaster(tas_bias, paste0(path, "/biascorrection/tas_hist_corr_", g, ".tif"), overwrite = TRUE)
  writeRaster(tasmin_bias, paste0(path, "/biascorrection/tasmin_hist_corr_", g, ".tif"), overwrite = TRUE)
  writeRaster(tasmax_bias, paste0(path, "/biascorrection/tasmax_hist_corr_", g, ".tif"), overwrite = TRUE)
  writeRaster(pr_bias, paste0(path, "/biascorrection/pr_hist_corr_", g, ".tif"), overwrite = TRUE)
  writeRaster(rh_bias, paste0(path, "/biascorrection/rh_hist_corr_", g, ".tif"), overwrite = TRUE)
  writeRaster(rad_bias, paste0(path, "/biascorrection/rad_hist_corr_", g, ".tif"), overwrite = TRUE)
  writeRaster(vpd_bias, paste0(path, "/biascorrection/vpd_hist_corr_", g, ".tif"), overwrite = TRUE)
  
  tmpFiles(remove=TRUE, current=T, orphan=T)
  gc()
  
}





### obtain tables with values for each grid cell -------------------------------

# reference system 
proj_nc <- "+proj=ob_tran +o_proj=longlat +o_lon_p=-162 +o_lat_p=39.25 +lon_0=180"
proj <- "+proj=longlat +datum=WGS84 +no_defs"
proj_forestmask <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

#europe poly
eu_mask_wgs <- rast("~/home/climate_data/EU_Mask_NetCDF_in_WGS.tif")!!!
CellsXY <- as.data.frame(eu_mask_wgs, xy = T)

# convert the points to the "normal" projection
points.orig <- cbind(id = 1:nrow(CellsXY), CellsXY[, c("x", "y")]) 
head(points.orig)
points <- points.orig
coordinates(points) <- ~x+y
points.v <- vect(points)

# transform grid cell coordinates
d <- sf::st_as_sf(points.v)
x_nc <- st_transform_proj(d, c("+proj=ob_tran +o_proj=longlat +o_lon_p=-162 +o_lat_p=39.25 +lon_0=180 +to_meter=0.01745329", proj))
st_crs(x_nc) = NA 
st_crs(x_nc) = proj_nc
nc_df <- data.frame(st_coordinates(x_nc[,2]))
head(nc_df)
nc_df_sp <- nc_df
coordinates(nc_df_sp) <- ~X+Y
points.v <- vect(nc_df_sp)

# test by backconverting
#d <- sf::st_as_sf(points.v)
#x <- st_transform_proj(d, c(st_crs(4326)$proj4string, "+proj=ob_tran +o_proj=longlat +o_lon_p=-162 +o_lat_p=39.25 +lon_0=180 +to_meter=0.01745329"))
#st_crs(x) = NA # ?! (!)
#st_crs() = 4326
#wgs_df <- data.frame(st_coordinates(x[,1]))
#head(wgs_df)

# define gcms
gcms <- c("MPI-M-MPI-ESM-LR", "ICHEC-EC-EARTH", "NCC-NorESM1-M")

fl_nme_public <- paste0(path, "/biascorrection/")

all_ncs_tas <- list.files(fl_nme_public, pattern = "tas_hist_corr/*")
all_ncs_tasmax <- list.files(fl_nme_public, patter = "tasmax_hist_corr/*")
all_ncs_tasmin <- list.files(fl_nme_public, patter = "tasmin_hist_corr/*")
all_ncs_pr <- list.files(fl_nme_public, patter = "pr_hist_corr/*")
all_ncs_rad <- list.files(fl_nme_public, patter = "rad_hist_corr/*")
all_ncs_vpd <- list.files(fl_nme_public, patter = "vpd_hist_corr/*")


# start loop
for(g in gcms){
  print(g)
  
  # read which files are from the selected GCM
  all_tas <- all_ncs_tas[which(grepl(paste(g), all_ncs_tas))]
  all_tasmax <- all_ncs_tasmax[which(grepl(paste(g), all_ncs_tasmax))]
  all_tasmin <- all_ncs_tasmin[which(grepl(paste(g), all_ncs_tasmin))]
  all_pr <- all_ncs_pr[which(grepl(paste(g), all_ncs_pr))]
  all_rad <- all_ncs_rad[which(grepl(paste(g), all_ncs_rad))]
  all_vpd <- all_ncs_vpd[which(grepl(paste(g), all_ncs_vpd))]
  
  all_tas_r <- rast(paste0(fl_nme_public, all_tas))
  all_tasmax_r <- rast(paste0(fl_nme_public, all_tasmax))
  all_tasmin_r <- rast(paste0(fl_nme_public, all_tasmin))
  all_pr_r <- rast(paste0(fl_nme_public, all_pr))
  all_rad_r <- rast(paste0(fl_nme_public, all_rad))
  all_vpd_r <- rast(paste0(fl_nme_public, all_vpd))
  
  # define the variables we want to extract data for
  vari <- c("tas", "tasmax", "tasmin", "pr", "rad", "vpd") #
  
  all_var <- c(all_tas_r, all_tasmax_r, all_tasmin_r, all_pr_r, all_rad_r, all_vpd_r)
  
  # extract the data  
  extract_df <- terra::extract(all_var, points.v)
  
  # organise the dataframe
  colnames(extract_df) <- c("id", vari)
  df.all <- cbind(extract_df[,1], points.orig[ ,c("x", "y")], nc_df[ , c("X", "Y")], extract_df[, c(2:ncol(extract_df))])
  
  # change rownames and colnames of dataframe
  rownames(df.all) <- c(1:nrow(df.all))
  colnames(df.all) <- c("ID", "wgs_x", "wgs_y", "rotpole_x", "rotpole_y", vari)
  df.all <- data.table(df.all)
  
  gc()
  tmpFiles(remove = T)
  
  
  clim_db_path <- "/.../climate_database/biascorr/"
  write_csv(df.all, paste0(clim_db_path, "/biascorr_", g, ".csv"))
  
  
}


### end ------------------------------------------------------------------------
