# =============================================================================
# Title: Training sample extraction and augmentation from Simulation Data
# Author: Marc Grünig

# 
# Description:
# This script extracts forest state transitions from a forest simulation 
# database. The extracted transitions serve as training data for a deep neural 
# network (DNN). The script includes functionalities for data augmentation to 
# enhance the robustness of the training dataset. Key tasks include:
# 
# 1. Connecting to the simulation SQLite database and extracting relevant 
#    simulation tables.
# 2. Using custom functions to process simulation data and identify vegetation 
#    transitions.
# 3. Creating augmentation samples by calculating residence times, target states, 
#    and historical vegetation states for each transition.
# 4. Compiling all processed data into a cohesive dataset for downstream analysis.
# 5. Writing the results back to the database for persistent storage.
# 
# Dependencies:
# - R Packages: DBI, tidyverse, ggplot2, parallel, doParallel, compiler, data.table
# - Custom Functions: create_examples_function.R (for transition and augmentation logic)
# - External Data: SQLite database (forest_simulation_db_v1.1.sqlite)
# 
# Outputs:
# - Processed examples of vegetation state transitions saved in a new database table 
#   (`examples_2m_augmentation`).
# - Checks and summaries for data consistency, including NA detection and state-change 
#   statistics.
# 
# Notes
# - we provide a example dataset to test the code: ./data_pipeline/metadata_expl.csv and ./data_pipeline/simdata_expl.csv
# - Ensure the simulation database is available at the specified path. It can be downloaded 
#   from: https://zenodo.org/records/12750180.
# - Adjust the `path` variable to point to the correct database location.
# - Ensure all required packages are installed and loaded before running the script.
# =============================================================================


# libraries
library(DBI)
library(tidyverse)
library(ggplot2)
library(parallel)
library(doParallel)
library(compiler)
library(data.table)


# define path to folder
path <- "/.../"


### load sim DB ----------------
simulation_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(path, "/forest_simulation_db_v1.1.sqlite"))

tables_con <- dbListTables(simulation_db)
simulations <- tables_con[grepl("simulation", tables_con)]

# in the case of the example data use:
# simdata <- read_csv(paste0(path, "/data_pipeline/simdata_expl.csv"))


### load functions -------------------
source("./functions/create_examples_function.R")
source("./functions/create_aug_examples_function.R")


### Run the functions for all simulations ---------------------
system.time({
  
  # loop over all simulation tables in the database
  examples_data <- as.list(1:length(simulations))
  
  r <- 0
  
  for(s in simulations){
    r <- r + 1
    print(s)
    
    # read the table
    tab <- dbReadTable(simulation_db, paste(s))
    tab <- as.data.table(tab)

    # in case of example data use:
    # tab <- as.data.table(simdata)
    
    sim_ids <- unique(tab$simulationID)
    examples_data[[r]] <- examples_data_fun(sim_ids)
    
  }
  
  

  
  # remove empty tables from list
  examples_data_new <- Filter(function(x) dim(x)[1] > 0, examples_data)
  examples_all_simulations <- do.call(bind_rows, examples_data_new)
  
  # save to SQLite DB
  dbWriteTable(conn = simulation_db, name = "examples_2m",
               value = examples_all_simulations, overwrite = T)
  
  dbDisconnect(simulation_db)
})

dim(expls_new)



# make some checks and obtain some numbers
dim(examples_all_simulations)
states <- unique(examples_all_simulations$svd_state)
states[grepl("NA", states)]
states <- unique(examples_all_simulations$target_state)
states[grepl("NA", states)]
length(states)

dim(examples_all_simulations[grepl("NA", examples_all_simulations$svd_state)])

# check whats the percentage of no change
no_change <- examples_all_simulations[examples_all_simulations$svd_state == examples_all_simulations$target_state,]
dim(no_change)



### now create augmentation training samples ----------------------------------------------

system.time({
  
  # loop over all simulation tables in the database
  examples_data <- as.list(1:length(simulations))
  
  r <- 0
  
  for(s in simulations){
    r <- r + 1
    print(s)
    # read the table
    tab <- dbReadTable(simulation_db, paste(s))
    tab <- as.data.table(tab)
    sim_ids <- unique(tab$simulationID)
    
    examples_data[[r]] <- examples_data_fun(sim_ids)
    
  }
  
  # remove empty tables from list
  examples_data_new <- Filter(function(x) dim(x)[1] > 0, examples_data)
  examples_all_simulations <- do.call(bind_rows, examples_data_new)
  
})

# save to SQLite DB
dbWriteTable(conn = simulation_db,
             name = "examples_2m_augmentation",
             value = examples_all_simulations, overwrite = T)

dbDisconnect(simulation_db)


