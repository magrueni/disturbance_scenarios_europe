# --------------------------------------------------------------
# Climate Data 02: Climate data preparation for SVD
# --------------------------------------------------------------
# Author: Marc Grünig
# Date: 2025-01-07
# 
# Description:
# In this script we prepare the climate data for SVD. We extract the daily data
# from our climate database obtained from Euro-Cordex data and store it.
# After this step, the daily climate data is compressed using a DNN (see methods for details).
# The compressed data is then used to drive climate in SVD.
# An example of the compressed data is stored in /svd/svd_example/climate/
# 
# Notes:
#  - We cannot share the climate database with daily data
#  - The results of this script are files with several GB size
#  
# --------------------------------------------------------------



library(terra)
library(dplyr)
library(raster)
library(compiler)
library(data.table)
library(parallel)
library(collapse)
library(sparklyr)
library(arrow)
library(collapse)
library(data.table)
library(Rcpp)
library(readr)



path <- "/.../"

#define options for tmp files storage
terraOptions(tempdir = "./tmp/raster")
gc()
tmpFiles(remove = T)



### functions ----------
# load reference grid an function to find grid ID
x_rast <- rast(paste0(path, "/reference_grids/reference_grid.tif"))
crs(x_rast) <- "+proj=longlat +datum=WGS84 +no_defs"
ref_tab <- read.csv(paste0(path, "reference_grids/reference_grid_tab.csv"))

get_point_id <- function(x){
  
  # extract the grid ID
  lon <- as.numeric(x[1])
  lat <- as.numeric(x[2])
  
  coords <- as.data.frame(cbind(lon, lat))
  
  point_extr <- terra::extract(x_rast, coords)[2]
  point_extr <- as.numeric(point_extr)
  
  if(is.na(point_extr)){
    closest_grid <- ref_tab[which.min(abs(lon - ref_tab$wgs_x) + abs(lat - ref_tab$wgs_y)), ]
    point_extr <- as.numeric(closest_grid[1])}
  
  return(point_extr)
  
}

### spark config ---------------------------------------------------------------
# Set memory allocation for whole local Spark instance
Sys.setenv("SPARK_MEM" = "13g")

# Set driver and executor memory allocations
config <- spark_config()
config$`spark.driver.memory` <- "24G"
config$`sparklyr.shell.driver-memory` <- "24G"
config$`spark.executor.memory` <- "24G"
config$`spark.driver.maxResultSize` <- "24G"
config$`spark.yarn.executor.memoryOverhead` <- "8G"
config$sparklyr.gateway.port = 8892
config$sparklyr.gateway.start.timeout = 180
config$sparklyr.gateway.connect.timeout = 1
config$`sparklyr.cores.local` <- 4
config$`spark.executor.cores` <- 4

options(java.parameters = "-Xmx8048m")


### data preparation ------
fl_pth <- "/.../"

###  read spark climate data for each simulation according the scenario and 
scenarios_fut <- c( "MPI-M-MPI-ESM-LR_rcp_8_5", "MPI-M-MPI-ESM-LR_rcp_4_5", "MPI-M-MPI-ESM-LR_rcp_2_6", 
                    "ICHEC-EC-EARTH_rcp_8_5", "ICHEC-EC-EARTH_rcp_4_5", "ICHEC-EC-EARTH_rcp_2_6",
                    "NCC-NorESM1-M_rcp_8_5", "NCC-NorESM1-M_rcp_4_5", "NCC-NorESM1-M_rcp_2_6")

scenarios_hist <- c("ICHEC-EC-EARTH_historical", "NCC-NorESM1-M_historical", "MPI-M-MPI-ESM-LR_historical")
scenarios <- c(scenarios_hist, scenarios_fut)


# load phenips module
Rcpp::sourceCpp(paste0(path, "/bb_module/phenips/phenips.cpp"))

# function to create aux data
aux_data <- function(clim, lai){
  
  pointid <- unlist(as.numeric(clim["point_id"]))
  lat <- r_df %>% fsubset(reference_grid == pointid) %>% dplyr::select(y) %>% unlist() %>% as.numeric()
  
  min_temp <- unlist(as.vector(clim[names(clim)[grepl("min_temp_", names(clim))]]))
  max_temp <- unlist(as.vector(clim[names(clim)[grepl("max_temp_", names(clim))]]))
  rad <- unlist(as.vector(clim[names(clim)[grepl("rad", names(clim))]]))
  vpd <- unlist(as.vector(clim[names(clim)[grepl("vpd", names(clim))]]))
  vpd <- mean(vpd[152:243])
  
  cmat <- as.matrix(as.data.frame(cbind(min_temp, max_temp, rad)))
  
  bbgen <- new(BarkBeetleGen)
  capture.output(bbgen$setup(latitude_deg = lat, lai = lai, do_log = FALSE))
  res <- bbgen$run(cmat)
  result <- c(pointid, res, vpd)
  return(result)
}


# get all grid cells
r_df <- as.data.frame(x_rast, xy = T)

gc()
tmpFiles(remove = T)

# and point ids from where clim data is extracted
point_ids <- unique(as.vector(unlist(r_df$reference_grid)))
length(point_ids)


for(s in scenarios){
  
  print(s)
  
  if(grepl("historical", s)){
    timesteps <- c("1981-1990", "1991-2000", "2001-2010")
  }else{
    timesteps <- c("2011-2020", "2021-2030", "2031-2040", "2041-2050", "2051-2060", "2061-2070", "2071-2080", "2081-2090", "2091-2100")
  }
  
  
  # load biascor df
  g <- strsplit(s, "_")[[1]][1]
  bias <- read_csv(paste0("/clim_data/biascorr_", g, "_v3.csv"))
  bias <- bias %>% dplyr::select(point_id = ID, tas_bias = tas, min_temp_bias = tasmin,
                                 max_temp_bias = tasmax, pr_bias = pr, rad_bias = rad, vpd_bias = vpd)
  # connect to climate database
  sc <- sparklyr::spark_connect(master = "local", config = config)
  spark_tbl_handle <- spark_read_parquet(sc, name="climate", path=paste0("/.../climate_database/spark_db_v3/", s, "/"), memory = FALSE)
  
  for(r in 1:length(timesteps)){
    
    timestep <- as.vector(c(as.numeric(substring(timesteps[r], 0, 4)): as.numeric(substring(timesteps[r], 6, 9))))
    
    timestep_list <- list()
    
    for(y in 1:length(timestep)){
      
      
      print(timestep[y])
      focal_year <- as.numeric(timestep[y])
      
      
      # need to open the open the other spark databases for the years 2006:2010
      if(grepl("historical", s) & focal_year > 2005){
        
        fut_scenarios <- scenarios_fut[grepl(substring(s, 0, nchar(s)-11), scenarios_fut)]
        
        clim_fut_batch <- list()
        
        for(sp_fut in 1:length(fut_scenarios)){
          
          spark_tbl_handle_fut <- spark_read_parquet(sc, name="climate", path=paste0("/mnt/dss/data/climate_database/spark_db_v3/", fut_scenarios[sp_fut], "/"), memory = FALSE)
          
          clim_tab_fut <- spark_tbl_handle_fut %>%
            filter(point_id %in% point_ids) %>%
            filter(year == focal_year) %>% 
            sparklyr::select(., c(point_id, year, month, day, tas, min_temp, max_temp, prec, rad, vpd)) %>%
            sparklyr::distinct() %>%
            collect() %>% 
            arrange(month, day)
          
          clim_tab_fut <- left_join(clim_tab_fut, bias, by = "point_id") %>% 
            mutate(tas = tas + tas_bias,
                   min_temp = min_temp + min_temp_bias,
                   max_temp = max_temp + max_temp_bias,
                   prec = prec * pr_bias,
                   rad = rad * rad_bias,
                   vpd = vpd * vpd_bias)
          
          clim_fut_batch[[sp_fut]] <- clim_tab_fut
          
          rm(spark_tbl_handle_fut)
          
        }
        
        clim_fut_batch_df <- data.table(do.call(bind_rows, clim_fut_batch))
        clim_tab <- collap(clim_fut_batch_df, ~ point_id + year + month + day, fmean, na.rm = TRUE)
        
      }else{
        
        clim_tab <- spark_tbl_handle %>%
          filter(point_id %in% point_ids) %>%
          filter(year == focal_year) %>% 
          sparklyr::select(., c(point_id, year, month, day, tas, min_temp, max_temp, prec, rad, vpd)) %>%
          sparklyr::distinct() %>%
          collect() %>% 
          arrange(month, day)
        
        clim_tab <- left_join(clim_tab, bias, by = "point_id") %>% 
          mutate(tas = tas + tas_bias,
                 min_temp = min_temp + min_temp_bias,
                 max_temp = max_temp + max_temp_bias,
                 prec = prec * pr_bias,
                 rad = rad * rad_bias,
                 vpd = vpd * vpd_bias)
        
      }
      
      batch <- clim_tab  %>% 
        tidyr::pivot_wider(names_from = c(month, day), values_from = c(tas, prec, rad, vpd, min_temp, max_temp)) %>% 
        dplyr::select(-tas_bias, -pr_bias, -min_temp_bias, -max_temp_bias, -rad_bias, -vpd_bias)
      colnames(batch) <- c("point_id", "year", paste0("tas_", 1:365), paste0("prec_", 1:365), paste0("rad_", 1:365),
                           paste0("vpd_", 1:365), paste0("min_temp_", 1:365), paste0("max_temp_", 1:365))
      
      batch_aux <- mclapply(asplit(batch, 1), aux_data, lai = 5, mc.cores = 12, mc.cleanup = T)
      batch_aux_df <- as.data.frame(do.call(rbind, batch_aux))
      colnames(batch_aux_df) <- c("point_id", "bbgen", "wintermin", "summervpd")
      batch <- batch %>%
        dplyr::select(point_id:vpd_365) %>% 
        left_join(., batch_aux_df, by = c("point_id"))

      timestep_list[[y]] <- batch

    }
    
    # but all years together
    batch_timestep <- do.call(bind_rows, timestep_list)
    
    # write out the batch for the timestep
    fname <- paste0(path, "/clim_data/daily_climate/all_grids_", timesteps[r], "_", s,".dat")
    write_feather(batch_timestep, fname, compression="lz4")
    rm(batch_timestep)
    gc()
    
  }
  
  # clear the sparklyr session
  sc %>% spark_session() %>% sparklyr::invoke("catalog") %>% sparklyr::invoke("dropTempView", "spark_tbl_handle")
  sc %>% spark_session() %>% sparklyr::invoke("catalog") %>% sparklyr::invoke("clearCache")
  sparklyr::spark_disconnect(sc)
  sparklyr::spark_disconnect_all()
  sc <- NULL
  spark_tbl_handle <- NULL
  gc()
  
}


### end -------------------

