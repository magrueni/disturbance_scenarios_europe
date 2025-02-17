# -----------------------------------------------------------------------------------------------------------
# Script for Extracting and Climate Data for Simulations
#
# Description:
# This script is designed to extract and process climate data associated with simulation metadata 
# and specific scenarios. It integrates Apache Spark for distributed data handling, utilizes 
# Terra for raster operations, and employs parallel processing for efficiency. The output is a 
# harmonized SQLite database with daily climate data for each simulation and scenario.
#
# Key Steps:
# 1. Configure Spark for distributed data handling with memory optimizations.
# 2. Load and process simulation metadata and climate data.
# 4. Extract bias-corrected daily climate variables for each simulation.
# 5. Store harmonized climate data into SQLite databases by scenario and user.
#
# Outputs:
# - Harmonized SQLite database containing simulation-specific climate data.
# - Daily climate variables include temperature, precipitation, radiation, and vapor pressure deficit.
#
# Dependencies:
# - R packages: `terra`, `sparklyr`, `dplyr`, `collapse`, `data.table`, `parallel`, `DBI`, `arrow`.
# - Input data: 
#   - Reference raster grid and climate grid metadata.
#   - Simulation metadata and timeseries data.
#   - Climate data stored in Apache Parquet format, organized by scenario.
#
# Notes:
# - Adjust Spark and memory configurations as needed for your environment.
# - Ensure input paths for climate data, metadata, and SQLite databases are correct.
# - Designed for high-performance environments; optimal runtime depends on system resources.
# - Again, we show this with an example dataset, as the whole climate database is too large
# 
# -----------------------------------------------------------------------------------------------------------


### script for extracting climate data for each simulation and compress it
### this script will be translated to python and is just here to set the stage
library(terra)
library(RColorBrewer)
library(sf)
library(dplyr)
library(DBI)
library(sparklyr)
library(arrow)
library(data.table)
library(collapse)
library(compiler)
library(parallel)
library(readr)



### spark config ---------------------------------------------------------------
# Set memory allocation for whole local Spark instance
# Sys.setenv("SPARK_MEM" = "26g")
# 
# config <- spark_config()
# config$`spark.driver.memory` <- "24g"
# config$`sparklyr.shell.driver-memory` <- "24g"
# config$`spark.executor.memory` <- "24g"
# config$`spark.driver.maxResultSize` <- "12g"
# config$`spark.yarn.executor.memoryOverhead` <- "4g"
# config$sparklyr.gateway.port = 8890
# config$sparklyr.gateway.start.timeout = 180
# config$sparklyr.gateway.connect.timeout = 180
# config$`sparklyr.cores.local` <- 8
# config$`spark.executor.cores` <- 8
# config$spark.executor.extraJavaOptions <- "-XX:MaxHeapFreeRatio=60 -XX:MinHeapFreeRatio=30"    # Set garbage collection settings
# config$spark.driver.extraJavaOptions <- "-XX:MaxHeapFreeRatio=60 -XX:MinHeapFreeRatio=30"    # Set garbage collection settings
# 
# options(java.parameters = "-Xmx8048m")
# 
# 


### functions --------------------
path <- "/.../"
x_rast <- rast("/reference_grid.tif")
crs(x_rast) <- "+proj=longlat +datum=WGS84 +no_defs"
ref_tab <- read.csv("/reference_grid_tab.csv")


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


# climate extraction function
extract_clim_func <- function(foc_year, clim_tab, sim_point, temp_offset, prec_offset, rad_bias, vpd_bias){

  clim_tab_sim <- clim_tab[point_id == sim_point & year == foc_year, .(tas, prec, rad, vpd)]
  clim_tab_sim[, c("tas", "prec", "rad", "vpd") := .(as.numeric(tas - temp_offset), as.numeric(prec / prec_offset),
                                                     as.numeric(rad * rad_bias), as.numeric(vpd * vpd_bias))]
  
  return(as.vector(c(Year_clim = foc_year, as.matrix(clim_tab_sim[,c("tas","prec","rad","vpd")]))))
  
}
extract_clim_func <- cmpfun(extract_clim_func)


# wrapper fun
wrapper_fun <- function(x){
  
  x <- as.vector(as.data.frame(x))
  sim_point <- get_point_id(as.data.frame(x[c("Lon", "Lat")]))
  
  simdata_new <- simdata[simulationID == as.numeric(x$ID)
                         & !is.na(Temp) & !is.na(Precip)]
  
  
  # get offset
  temp_offset <- as.numeric(x["temp_offset"])
  prec_offset <- as.numeric(x["prec_offset"])
  
  # get vpd and rad bias
  vpd_bias <- as.numeric(bias_df[bias_df$ID == sim_point, "vpd"])
  rad_bias <- as.numeric(bias_df[bias_df$ID == sim_point, "rad"])
  
  
  # load the timeseries
  focal_ts <- as.vector(unlist(
    timeseries_focal_df[timeseries_focal_df[, "uniqueID"] == as.numeric(x["uniqueID"]) & 
                          timeseries_focal_df[, "simulationID"] == as.numeric(x["ID"]), ][-c(1,2)]))
  
  simdata_new <- cbind(simdata_new, Year_clim = na.omit(focal_ts))
  
  unique_focal_ts <- unique(na.omit(focal_ts))
  unique_focal_ts_list <- split(unique_focal_ts, seq(length(unique_focal_ts)))  
  
  clim_list <- lapply(unique_focal_ts_list, FUN = extract_clim_func,
                      clim_tab = clim_tab, sim_point = sim_point,
                      temp_offset = temp_offset, prec_offset = prec_offset,
                      rad_bias = rad_bias, vpd_bias = vpd_bias)
  
  clims <- t(rbindlist(list(clim_list)))
  
  colnames(clims) <- c("Year_clim", paste0("tas_", 1:365), paste0("prec_", 1:365),
                       paste0("rad_", 1:365), paste0("vpd_", 1:365))
  sim_df <- simdata_new %>% left_join(., as.data.frame(clims), by = c("Year_clim"))
  
  return(sim_df)
  
}
wrapper_fun <- cmpfun(wrapper_fun)



### load metadata --------------------------------------------------------------------------

# example dataset
metadata <- read_csv(paste0(path_expl, "/data_pipeline/metadata_processed_expl.csv"))
samples <- read_csv(paste0(path_expl, "/data_pipeline/examples_expl.csv"))


### order metadata so we can training_db bunch of climate data from spark -------------------
# for instance order with target scenario in order to not switch between the scenario folders 
# and load all grid cells of 1000 simulations at the time
metadata <- metadata %>% arrange(scenario)
metadata <- metadata %>% filter(!is.na(scenario)) %>% data.table()


### loop over the simulations by groups of scenarios --------------------------------------
scenarios <- metadata %>% dplyr::select(scenario) %>% distinct() %>% unlist() %>% as.vector()

for(s in scenarios){
    
  
  # load clim bias correction file
  bias_df <- read_csv(paste0(path_expl, "/data_pipeline/clim_bias_expl.csv"))
  
    # subset the metadata by scenario
  meta_sub_scen <- metadata %>% collapse::fsubset(scenario == s)
    
  # get the coords
  coords <- meta_sub_scen[, c("Lon", "Lat")] %>%
      distinct() %>% data.table()
    
  my_point_id <- apply(coords, 1, get_point_id)
    
  users <- as.vector(unlist(unique(meta_sub_scen[, "uniqueID"])))
    
  ### read spark climate data for each simulation according the scenario and 
  # sc <- sparklyr::spark_connect(master = "local", config = config)
  # spark_tbl_handle <- spark_read_parquet(sc, name="climate",
  #                                          path=paste0("/mnt/dss/data/climate_database/spark_db_v3/", s, "/"),
  #                                          memory = FALSE)
  #   
  # clim_tab <- spark_tbl_handle %>%
  #     filter(point_id %in% my_point_id) %>% 
  #     select(., c(point_id, year, month, day, tas, prec, rad, vpd)) %>%
  #     collect %>%
  #     setDT()
    
    
  # clear the sparklyr session
  # sc %>% spark_session() %>% sparklyr::invoke("catalog") %>% 
  #   sparklyr::invoke("dropTempView", "spark_tbl_handle")
  # sc %>% spark_session() %>% sparklyr::invoke("catalog") %>%
  #   sparklyr::invoke("clearCache")
  # sparklyr::spark_disconnect(sc)
  # sparklyr::spark_disconnect_all()
  # sc <- NULL
  # spark_tbl_handle <- NULL
  # rm(sc)
  # gc()
  
  # get expl dataset instead
  clim_tab <- read_csv(paste0(path_expl, "/data_pipeline/climate_expl.csv")) %>% 
    setDT()
  
  # loop over unique IDs
  for(u in users){

      # separate the meta data of the user
      meta_sub_scen_user <- meta_sub_scen %>% collapse::fsubset(uniqueID == u)
      
      # load the table with the new timeseries
      timeseries_focal_df <- read_csv(paste0(path_expl, "/data_pipeline/time_series.csv"))
      
      # split to chunks if too large
      chnks <- split(c(1:nrow(meta_sub_scen_user)),
                     ceiling(seq_along(c(1:nrow(meta_sub_scen_user)))/3000))
      
      # loop over chunks
      for(chnk in 1:length(chnks)){
        
        # print progress  
        print(paste0(chnk," / ", length(chnks)))
        
        # take the metadata of the chunk
        meta_chnk <- meta_sub_scen_user[chnks[[chnk]][1]:chnks[[chnk]][length(chnks[[chnk]])], ]
        
        # convert it to a lst so we can mc lapply over it
        meta_chnk_list <- split(meta_chnk, seq(nrow(meta_chnk)))  
        
        # get all ids
        chnk_ids_sub <- as.vector(unlist(meta_chnk[,"ID"]))
        
        # load the simulation data of all those IDs
        # simulation_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(path, "/forest_simulation_db_v1.1.sqlite"))
        # tables_con <- dbListTables(simulation_db)
        # sim_tables <- tables_con[grepl(paste0("simulation"), tables_con)]
        # 
        # simdata <- simulation_db %>%
        #   tbl(sim_tables[grepl(u, sim_tables)]) %>%
        #   filter(simulationID %in% chnk_ids_sub) %>%
        #   collect() %>%
        #   arrange(simulationID)
        # simdata <- setDT(simdata)
        # dbDisconnect(simulation_db)
        
        # expl dataset instead
        simdata <- read_csv(paste0(path_expl, "/data_pipeline/simdata_expl.csv")) %>% 
          setDT()
        
        
        simdata_new_list <- parallel::mclapply(meta_chnk_list, FUN = wrapper_fun,
                                               mc.cores = 40, mc.cleanup = TRUE)
        #gc(verbose = F)
        
        sim_df_user <- rbindlist(simdata_new_list)
        sim_df_user <- sim_df_user %>% 
          dplyr::select(uniqueID:Proportion5, sp_comp = sp_state2, lai_class = lai_state,
                        dom_height = height_state, veg_state = svd_state, tas_1:vpd_365) %>%
          left_join(., meta_chnk %>%
                      dplyr::select(uniqueID, simulationID = ID,
                                    WHC:SoilDepth, AvailableNitrogen),
                    by = c("uniqueID", "simulationID") ) %>% 
          dplyr::select(uniqueID:veg_state, WHC:AvailableNitrogen, tas_1:vpd_365)
        
        # write to DB
        training_db <- DBI::dbConnect(
          RSQLite::SQLite(),
          paste0(path, "/data_portal/simulation_data/harmonized_db_", s, "_v1.1.sqlite"))
        
        
        if(chnk == 1){
          dbWriteTable(conn = training_db, name = paste0("harmonized_simulations_", u, "_", s),
                       value = sim_df_user, overwrite = T)
        }else{
          dbWriteTable(conn = training_db, name = paste0("harmonized_simulations_", u, "_", s),
                       value = sim_df_user, append = T)
        } 
        
        
        
        dbDisconnect(training_db)
        
        sim_df_user <- NULL
        simdata_new_list <- NULL
        simdata <- NULL
        meta_chnk_list <- NULL
        meta_chnk <- NULL
        training_db <- NULL
        simulation_db <- NULL
        rm(simdata, sim_df_user, meta_chnk, meta_chnk_list,
           simdata_new_list, training_db, simulation_db,
           chnk_ids_sub, tables_con, sim_tables)
        
        tmp_dir <- tempdir()
        temp_files <- list.files(tmp_dir, full.names = T, pattern = "^file")
        file.remove(temp_files)
        
        for (i in 1:10){gc(reset = T)} 
        
      } # close chnks
      
      
      rm(timeseries_focal_df, meta_sub_scen_user)
      for (i in 1:10){gc(reset = T)} 
      
    } # close user

  print(paste0("all user time: ",  users_time[3]))
  
  
  clim_tab <- NULL
  rm(clim_tab, meta_sub_scen)
  for (i in 1:10){gc(reset = T)} 
  

}



