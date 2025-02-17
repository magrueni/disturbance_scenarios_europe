# -----------------------------------------------------------------------------------------------------------
# Script for extracting climate time series
# Auhtor: Marc Grünig
#
# Description:
# This script processes climate data from climate data together with simulation data.
# It focuses on extracting temperature and precipitation trends, calculating their slopes 
# over time, as well as aggregating annual statistics (mean and sum) for each geographic location. 

# In a first step, we assign the best matching climate change scenario to each forest simulation
# In the second step, we create a timeseries based on the euro-cordex data that matches the simulation data
# This allows then to assign daily climate data to the simulation years.

# Dependencies:
# - Input data:
#     - A reference grid raster file (`reference_grid.tif`) and related metadata.
#     - Climate data stored in Apache Parquet format and managed in a distributed Spark database.
#     - Simulation metadata stored in an SQLite database.
#
# Notes: we provide example data as the climate data cannot be provided
#
# -----------------------------------------------------------------------------------------------------------



library(raster)
library(terra)
library(RColorBrewer)
library(sf)
library(dplyr)
library(DBI)
library(stars)
library(ggplot2)
library(data.table)
library(collapse)
library(parallel)
library(compiler)



## load functions
source("functions/get_point_id.R")

path <- "/.../"


# example datasets
# metadata <- read_csv(paste0(path, "/data_pipeline/metadata_expl.csv"))
# simdata <- read_csv(paste0(path, "/data_pipeline/simdata_expl.csv"))
# climate_slopes <- read_csv(paste0(path, "/data_pipeline/climate_slopes.csv"))




### load the simulation DB -----------------------------------
simulation_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(path, "/forest_simulation_db_v1.1.sqlite"))
tables_con <- dbListTables(simulation_db)

# load metadata
metadata <- dbReadTable(simulation_db, "metadata_all_cleaned_soil")

# load the climate slopes, extracted from the climate database in the previous script
climate_slopes <- dbReadTable(simulation_db, "climate_slopes_coordinates") %>% data.table()

# isolate the simulation data tables
sim_tables <- tables_con[grepl(paste0("simulation"), tables_con)]


# create rcp and gcm list to choose from
rcp_list <- list(c("historical"), c("rcp_2_6"), c("rcp_4_5"), c("rcp_8_5"), c("rcp_4_5", "rcp_8_5"), c("rcp_2_6", "rcp_4_5", "rcp_8_5"))
gcm_list <- list(c("ICHEC-EC-EARTH"), c("MPI-M-MPI-ESM-LR"), c("NCC-NorESM1-M"), c("ICHEC-EC-EARTH", "MPI-M-MPI-ESM-LR", "NCC-NorESM1-M")) 



### Create a function that compares the slope of the simulation with slopes from the different climate scenarios -----------------------------------------------
scenario_compare <- function(x, column_names){
  
  x <- as.vector(x)
  names(x) <- column_names
  sim_id <- as.numeric(x["ID"])
  
  # get the climate data from the simulation
  simdata <- simulation_data %>%
    filter(simulationID == as.numeric(sim_id)) %>% collect()
  
  # if empty return NAs
  if(nrow(simdata) == 0){return(c(NA, NA, NA, NA))}
  
  # get the temp and precip and calibrate linear models
  simdata <- simdata %>% collapse::fsubset(!is.na(Temp) & !is.na(Precip))
  sim_temp_mod <- collapse::flm(Temp ~ Year, data = simdata)
  sim_prec_mod <- collapse::flm(Precip ~ Year, data = simdata)
  
  # calculate averages
  sim_temp_avg <- mean(simdata$Temp)
  sim_prec_avg <- mean(simdata$Precip)
  
  # get climate for the gridcell
  sim_coords <- x[c("Lon", "Lat")]
  grid_id <- get_point_id(sim_coords)
  
  # check which rcp we have to look at
  rcp_indicator <- ifelse(grepl("aseline", x[c("Climate")]) | grepl("bserved", x[c("Climate")]), 1, 
                          ifelse(grepl("26", x[c("Climate")]) | grepl("2.6", x[c("Climate")]), 2, 
                                 ifelse(grepl("45", x[c("Climate")]) | grepl("4.5", x[c("Climate")]), 3, 
                                        ifelse(grepl("85", x[c("Climate")]) | grepl("8.5", x[c("Climate")]), 4,
                                               ifelse(grepl("6.0", x[c("Climate")]) | grepl("60", x[c("Climate")]), 5,
                                                      ifelse(grepl("SRES", x[c("Climate")]), 5, 6))))))
  
  rcps <- rcp_list[[rcp_indicator]]
  
  # check if there is information on the gcm in the metadata
  gcm_indicator <- ifelse(is.na(x[c("GCM")]), 4,
                          ifelse(grepl("ICHEC", x[c("GCM")]) | grepl("EARTH", x[c("GCM")]) | grepl("earth", x[c("GCM")]), 1, 
                                 ifelse(grepl("MPI", x[c("GCM")]) | grepl("mpi", x[c("GCM")]), 2, 
                                        ifelse(grepl("NCC", x[c("GCM")]) | grepl("NorESM1", x[c("GCM")]) | grepl("ncc", x[c("GCM")]), 3, 4))))
  
  gcms <- gcm_list[[gcm_indicator]]
  
  # create empty output DF
  diff_df <- NULL
  
  # loop over the potential gcms and rcps
  for(g in gcms){
    
    for(c in rcps){
      
      # for the coordinates of the simulation we get the precip and temp slope from the climate database for the focal scenario
      focal_gridcell <- climate_slopes %>% collapse::fsubset(point_id == grid_id) %>% collapse::fselect(-Lat, -Lon) %>% distinct() %>% as.data.frame()
      temp_slope_clim <- as.numeric(focal_gridcell[colnames(focal_gridcell) == paste0("temp_slope_", gsub("-", ".", g), "_", c)])
      prec_slope_clim <- as.numeric(focal_gridcell[colnames(focal_gridcell) == paste0("prec_slope_", gsub("-", ".", g), "_", c)])
      
      # calculate the deltas. temp is additive, precip multiplicative
      delta_t <- abs(as.numeric(sim_temp_mod[2,1]) - temp_slope_clim)
      # the alternative is to do it like for 
      delta_p <- abs(as.numeric(sim_prec_mod[2,1]) - prec_slope_clim)
      
      # add together with a arbitrary 0.1 weight for precip
      dev_metric <- delta_t + 0.1 * delta_p
      dev_metric <- data.table(dev_metric)
      colnames(dev_metric) <- paste0(g, "_", c)
      
      diff_df <- cbind(diff_df, dev_metric)
    }
  }
  
  # get the closest scenario
  min.val <- colnames(diff_df)[which.min(diff_df)]
  
  # calculate the offset from the averages
  temp_offset <- as.numeric(focal_gridcell[colnames(focal_gridcell) == paste0("temp_avg_", gsub("-", ".", min.val))]) - sim_temp_avg
  prec_offset <- as.numeric(focal_gridcell[colnames(focal_gridcell) == paste0("prec_avg_", gsub("-", ".", min.val))])/sim_prec_avg
  
  # return the values for scenario, temp offset and prec offset
  if(is.na(min.val)){return(c(sim_id, NA, temp_offset, prec_offset))}else{
    return(c(as.numeric(sim_id), as.character(min.val), as.numeric(temp_offset), as.numeric(prec_offset)))
  }
  
  gc()
}

# compile the function to make it faster
scenario_compare_cmp <- cmpfun(scenario_compare)


### process data ----------------------------------------------------------------------------------

users <- unique(metadata$uniqueID)
target_scen_all_list <- list()

for(u in 1:length(users)){
    
  print(users[u])
    
  #open database
  simulation_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(path, "/forest_simulation_db_v1.1.sqlite"))
  tables_con <- dbListTables(simulation_db)
  unique(tables_con)
  
  userID <- users[u]
    
  # subset meta data
  metadata_sub <- metadata %>% collapse::fsubset(uniqueID == userID)
    
  # if there is nothing, we skip
  if(nrow(metadata_sub) == 0){
    cat(paste0(userID))
    next}
    
  # convert to list to enable mclapply
  metadata_sub_list <- split(metadata_sub, seq(nrow(metadata_sub)))  
    
  # load the simulation data of the user
  # in case of expl dataset: simulation_data <- simdata %>%
  # dplyr::select(simulationID, Year, Temp, Precip) 
  
  simulation_data <- tbl(simulation_db, sim_tables[grepl(userID, sim_tables)]) %>%
      dplyr::select(simulationID, Year, Temp, Precip) %>% collect()

  # apply function in parallel
  target_scen <- parallel::mclapply(metadata_sub_list, FUN = scenario_compare_cmp, column_names = names(metadata_sub_list[[1]]), mc.cores = 16)
    
  # combine to a table 
  target_scen_df <- data.table(do.call(rbind, target_scen))
    
  colnames(target_scen_df) <- c("ID", "scenario", "temp_offset", "prec_offset")
  target_scen_df <- target_scen_df %>%
      mutate_at(c('ID', 'temp_offset', 'prec_offset'), as.numeric) %>% 
      mutate(uniqueID = userID)

  # add to list of users  
  target_scen_all_list[[u]] <- target_scen_df
    
  # close DB  
  dbDisconnect(simulation_db)
    
}# close users loop


# combine all tables
target_scen_all <- data.table(do.call(rbind, target_scen_all_list))
dim(target_scen_all)

# bind together to metadata and save in the SQLite DB
metadata_fit_scenario <- metadata %>% 
  left_join(., target_scen_all, by = c("ID", "uniqueID"))

metadata_fit_scenario <- metadata_fit_scenario %>% filter(!is.na(scenario))
head(metadata_fit_scenario)


dbWriteTable(conn = simulation_db, name = "metadata_with_fitting_scenario_v4", value = metadata_new, overwrite = T)
dbDisconnect(simulation_db)




### time series extraction --------------------------------------------------------------------

# next step is to extract the time series of years that matches best to the simulation series.
# this is done to extract daily climate data for each simulation year that matches the original climate data as good as possible.
# essentially, we compare each year of the simulation data with the climate database (euro-cordex) from the previously matched scenario.
# for this we use the temp and prec offset we calculated before. we check which year in the climate database has the closest values in mean temp and prec.
# from the 3 best matching years, we sample 1 to not get the same year in many cases.


# create time series function. This function compares the temp and prec from the 
# scenarios with the one from the simulations
time_series_creation <- function(x, temp, prec){
  
  dev_metric <- unlist(as.vector(abs(as.numeric(x["Temp"]) - temp) + (0.1 * as.numeric(x["Precip"])/prec)))
  
  min.year.val <- sample(names(base::sort(dev_metric)[1:3]), 1)
  min.year.val <- as.numeric(gsub("X", "", min.year.val, perl = TRUE))
  
  return(min.year.val)
  
}

time_series_creation_cmp <- cmpfun(time_series_creation)


# core function
# this function takes the climate data of each simulation and 
# uses then the function above to match it with the best fitting scenario
my_function <- function(meta_data, col_names){
  
  # obtain the sim id
  meta_data <- as.vector(meta_data)
  names(meta_data) <- col_names
  sim_id <- as.numeric(meta_data["ID"])
  
  # get the simulation data and extract temp and precip
  simdata_sub <- simulation_data %>%
    filter(simulationID == as.numeric(sim_id)) %>%
    collapse::fsubset(!is.na(Temp) & !is.na(Precip)) %>%
    collapse::fselect(Year, Temp, Precip)
  
  # get scenrio grid id
  sim_coords <- as.numeric(meta_data[c("Lon", "Lat")])
  grid_id <- get_point_id(sim_coords)
  
  # get the information of the scenario from the metadata
  target_scen <- as.character(meta_data["scenario"])
  target_scen <- gsub("-", "_", target_scen)
  
  temp_data_scen <- eval(parse(text = paste0("temp_", target_scen)))
  prec_data_scen <- eval(parse(text = paste0("prec_", target_scen)))
  
  # subset the temp data for the coordinate
  temp_data <- temp_data_scen %>% 
    collapse::fsubset(point_id == grid_id) 
  temp_data <- data.table(unique(temp_data[, 4:ncol(temp_data)]))
  temp_data <- temp_data - meta_data$temp_offset
  
  # subset the prec data 
  prec_data <- prec_data_scen %>% 
    collapse::fsubset(point_id == grid_id) 
  prec_data <- data.table(unique(prec_data[, 4:ncol(prec_data)]))
  prec_data <- prec_data / meta_data$prec_offset
  
  # apply matching function
  new_time_series <- apply(simdata_sub, 1, time_series_creation_cmp, temp = temp_data, prec = prec_data)
  
  return(new_time_series)
  
  
  # # testing
  # temp_data_test <- temp_data %>% pivot_longer(cols = everything(), names_to = "Year", values_to = "Temp") %>%
  #   select(Year, Temp) %>%
  #   mutate(Year = as.numeric(gsub("X", "", Year))) %>%
  #   mutate(Year = Year - 2006)
  # 
  # new <- new_time_series - 2006
  # new <- as.data.frame(new) %>%
  #   mutate(Year = new) %>%
  #   left_join(., temp_data_test, by = "Year")
  # 
  # ggplot() +
  #   geom_point(data=simdata_sub, aes(x=Year, y=Temp), color='green') +
  #   geom_point(data=temp_data_test, aes(x=Year, y=Temp), color='red') +
  #   geom_point(data=new, aes(x=Year, y=Temp), color='blue')
  # 
  # 
  # prec_data_test <- prec_data %>% pivot_longer(cols = everything(), names_to = "Year", values_to = "Prec") %>%
  #   select(Year, Prec) %>%
  #   mutate(Year = as.numeric(gsub("X", "", Year))) %>%
  #   mutate(Year = Year - 2006)
  # 
  # new <- new_time_series - 2006
  # new <- as.data.frame(new) %>%
  #   mutate(Year = new) %>%
  #   left_join(., prec_data_test, by = "Year")
  # 
  # ggplot() +
  #   geom_point(data=simdata_sub, aes(x=Year, y=Precip), color='green') +
  #   geom_point(data=prec_data_test, aes(x=Year, y=Prec), color='red') +
  #   geom_point(data=new, aes(x=Year, y=Prec), color='blue')
  # 
  # 
  
  
}

myFuncCmp <- cmpfun(my_function)


# load data from sqlite -----------------------------------------------------------------------------------

# example datasets
metadata <- read_csv(paste0(path, "/data_pipeline/metadata_processed_expl.csv"))
simdata <- read_csv(paste0(path, "/data_pipeline/simdata_expl.csv"))
climate_slopes <- read_csv(paste0(path, "/data_pipeline/climate_slopes.csv"))


# define all scenarios
gcms <- c("ICHEC-EC-EARTH", "MPI-M-MPI-ESM-LR", "NCC-NorESM1-M")
rcps <- c("rcp_2_6", "rcp_4_5", "rcp_8_5", "historical")


# load precalculated means for temp and precip - this is only one grid id for demonstration...
# the full dataset is too extensive
for(g in gcms){
  
  for( c in rcps){
    
    # combine gcms and rcps to scneario names
    scen_name <- paste0(g, "_", c)
    scen_name <- gsub("-", "_", scen_name)
    
    # read in climate means
    temp_data <- read_csv(paste0("/data/public/Projects/Resonate/publi/data_pipeline/clim_slopes/temp_", scen_name, ".csv"))
    assign(paste0("temp_", scen_name), temp_data)

    prec_data <- read_csv(paste0("/data/public/Projects/Resonate/publi/data_pipeline/clim_slopes/prec_", scen_name, ".csv"))
    assign(paste0("prec_", scen_name), prec_data)

  }
}


# get list of users we need to loop over from the harmonized database

# users <- unique(metadata$uniqueID)
# expl
users <- c("1005")



# loop over users
for(u in 1:length(users)){
    
    # connect to DB
    simulation_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(path, "/forest_simulation_db_v1.1.sqlite"))
    
    # monitor progress
    print(users[u])
    userID <- users[u]
    
    # subset meta data
    metadata_sub <- metadata %>% collapse::fsubset(uniqueID == userID) %>% filter(!is.na(scenario))
    
    # if empty metadata skip
    if(nrow(metadata_sub) == 0){
      cat(paste0(users[u]))
      next}
    
    # for the big user 
    if(nrow(metadata_sub) > 50000){
      
      chnks <- split(c(1:nrow(metadata_sub)), ceiling(seq_along(c(1:nrow(metadata_sub)))/50000))
      
      # load simulation data
      simulation_data <- tbl(simulation_db, sim_tables[grepl(userID, sim_tables)]) %>%
        dplyr::select(simulationID, Year, Temp, Precip) %>% collect()
      
      
      for(chn in 1:length(chnks)){
        
        print(paste0(chn, " / ", length(chnks)))
        
        metadata_sub_chnk <- metadata_sub[chnks[[chn]][1]:chnks[[chn]][length(chnks[[chn]])], ]
        
        # convert to list
        metadata_sub_list <- split(metadata_sub_chnk, seq(nrow(metadata_sub_chnk)))  
        
        # lapply over all <- entries and apply matching functions defined above
        set.seed(10)
        new_series <- parallel::mclapply(metadata_sub_list, FUN = myFuncCmp, col_names = names(metadata_sub_list[[1]]), mc.cores = 40)
        # rearrange 
        new_series_df <- data.table(t(sapply(new_series, "length<-", max(lengths(new_series)))))
        colnames(new_series_df) <- paste0("Year", seq_along(new_series_df))
        
        # add user ID
        time_series_df <- cbind(uniqueID = rep(userID, nrow(metadata_sub_chnk)), simulationID = metadata_sub_chnk$ID, new_series_df)
        
        # write out to database
        if(chn == 1){
          dbWriteTable(conn = simulation_db, name = paste0("user_", users[u], "_timeseries"), value = time_series_df, overwrite = T)
        }else{
          dbWriteTable(conn = simulation_db, name = paste0("user_", users[u], "_timeseries"), value = time_series_df, overwrite = F, append = T)
        }
        
      }
      
      # clean up
      rm(simulation_data)
      gc()
      
      # disconnect 
      dbDisconnect(simulation_db)
      
    }else{
      
      # convert to list
      metadata_sub_list <- split(metadata_sub, seq(nrow(metadata_sub)))  
      
      # load simulation data
      simulation_data <- tbl(simulation_db, sim_tables[grepl(userID, sim_tables)]) %>%
        dplyr::select(simulationID, Year, Temp, Precip) %>% collect()
      
      # lapply over all metadata entries and apply matching functions defined above
      set.seed(10)
      new_series <- parallel::mclapply(metadata_sub_list, FUN = myFuncCmp, col_names = names(metadata_sub_list[[1]]), mc.cores = 40)
      
      # rearrange 
      new_series_df <- data.table(t(sapply(new_series, "length<-", max(lengths(new_series)))))
      colnames(new_series_df) <- paste0("Year", seq_along(new_series_df))
      
      # add user ID
      time_series_df <- cbind(uniqueID = rep(userID, nrow(metadata_sub)), simulationID = metadata_sub$ID, new_series_df)
      
      # write out to database
      dbWriteTable(conn = simulation_db, name = paste0("user_", users[u], "_timeseries"), value = time_series_df, overwrite = T)
      
      # clean up
      rm(simulation_data)
      gc()
      
      # disconnect 
      dbDisconnect(simulation_db)
      
  }
}
  



