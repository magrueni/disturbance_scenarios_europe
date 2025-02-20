# --------------------------------------------------------------
# Initial landscape Scipt 02: Create vegetation stacks
# --------------------------------------------------------------
# Author: Marc Grünig
# Date: 2025-01-07
# 
# Description:
# The script uses two datasets, specifically Brus and Bonanella vegetation data,
# to create raster layers representing the spatial distribution of tree species in Europe.
# The workflow involves loading data, reprojecting, normalizing species proportions,
# and saving the results for further analysis. The resulting raster stacks show for each 
# pixel a species share for each species.


# Notes:
# - The resulting data is quite big (8Gb), but can be provided upon request
# -----------------------------------------------------------------



### load libraries
library(terra)
library(raster)
library(dplyr)
library(sf)
library(lookup)
library(data.table)
library(stringr)
library(tidyr)
library(data.table)
library(parallel)


# options settings
terraOptions(tempdir = "./tmp/raster")
gc()
tmpFiles(remove = T)
tempdir(check = FALSE)

path <- "/.../"
path_sp <- paste0(path, "/.../") # path to Bonnanella data


# define used projections
proj_wgs <- "+proj=longlat +datum=WGS84 +no_defs"
proj_leae <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"


# get vegetation from Brus maps
dom_veg <- list.files(paste0(path, "/EU_TreeMap_Brus_etal"), 
                      pattern = "_wgs.tif", full.names = TRUE)
dom_veg <- dom_veg[!grepl("Dominant", dom_veg)]
r <- rast(dom_veg)

# assign names
brus_names <- c("abal", "algl", "bepe", "broa", "cabe",
                "casa", "coni", "eugl", "fasy", "frex", 
                "lade", "piab", "pine", "pipi", "pisy",
                "potr", "psme", "oaks", "quro", "rops")


# define species
bonanella_species <- c("abies.alba", "castanea.sativa", "corylus.avellana",
                       "fagus.sylvatica", "olea.europaea", "picea.abies",
                       "pinus.halepensis", "pinus.nigra", "pinus.pinea", 
                       "pinus.sylvestris", "prunus.avium","quercus.ilex", 
                       "quercus.robur", "quercus.suber")

broa_lookup <- c("castanea.sativa", "prunus.avium",
                 "olea.europaea", "quercus.cerris", "quercus.ilex",
                 "quercus.robur", "quercus.suber")

coni_lookup <- c("abies.alba", "picea.abies", "pinus.halepensis",
                 "pinus.nigra",  "pinus.pinea", "pinus.sylvestris")

pine_lookup <- c("pinus.halepensis", "pinus.nigra",  "pinus.pinea",
                 "pinus.sylvestris")

oaks_lookup <- c("quercus.ilex", "quercus.robur", "quercus.suber")


# load bonanella maps
bon_rasts <- list.files(paste0(path_sp),
                        pattern = paste0("v0.3.tif"),
                        full.names = T)



### Building bouding boxes for the full canopy height dataset ---------------------------------------------------------
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


### start looping over the 3° tiles ---
for(i in 1:nrow(bbox)){
  
  # print progress  
  cat(paste0(i, " / ", nrow(bbox), " \n"))
  
  # get the canopy height raster as reference. If not available this means the tile is empty.
  name <- paste0(bbox[i, "top.coord"], bbox[i, "left.coord"])
  cat(paste0(name, " \n"))
  
  if(file.exists(paste0(path, "/ETH_GlobalCanopyHeight_10m_2020_", name, "_Map_aggregated_100m.tif"))){
    height_rast <- rast(paste0(path, "/ETH_GlobalCanopyHeight_10m_2020_", name, "_Map_aggregated_100m.tif"))
  }else{
    cat("height raster not available! \n")
    next
  }
  
  extent <- ext(height_rast)
  
  height_rast_laea <- terra::project(height_rast, proj_leae)
  veg_stack <- terra::crop(r, extent, snap = "out")
  veg_stack <- terra::project(veg_stack, height_rast, threads = TRUE)
  origin(veg_stack) <- origin(height_rast)
  
  # if empty value skip
  if(length(na.omit(values(veg_stack))) == 0){
    cat("vegetation stack has no values! \n")
    next}

  # gap filling
  veg_stack <- terra::focal(veg_stack, w = c(9, 9), fun = "mean", na.rm = T, na.policy = "only")
  veg_stack <- terra::focal(veg_stack, w = c(17, 17), fun = "mean", na.rm = T, na.policy = "only")
  veg_stack <- terra::mask(veg_stack, height_rast)
  r_df <- as.data.frame(veg_stack, xy = T, na.rm = F)
  colnames(r_df) <- c("x", "y", brus_names)
  
  
  r_df <- r_df %>% mutate_if(is.factor, as.numeric)
  r_df[is.na(r_df)] <- 0
  
  
  ### do the magic -------------------------------------------------------------------------------------------------------
  
  
  ## for broadleaved ----
  broa <- dom_veg[grepl("Broadleaves", dom_veg)]
  r_broa <- rast(broa)
  r_broa <- terra::crop(r_broa, height_rast)
  r_broa <- terra::project(r_broa, height_rast)
  r_broa <- terra::mask(r_broa, height_rast)
  
  # get bonanella maps for the broaferous
  broa_bon_rasts <- bon_rasts[grepl(paste(broa_lookup, collapse = "|"), bon_rasts)]
  broa_bon_rasts <- rast(broa_bon_rasts)
  names(broa_bon_rasts) <- broa_lookup
  
  # reproject to same epsg and res as height rast
  broa_bon_rasts <- terra::crop(broa_bon_rasts, height_rast_laea, snap = "out")
  broa_bon_rasts <- terra::project(broa_bon_rasts, height_rast, method = "average")
  broa_bon_rasts <- terra::mask(broa_bon_rasts, height_rast)
  
  # convert to dataframe
  broa_bon_df <- as.data.frame(broa_bon_rasts, xy = T, na.rm = F)
  broa_bon_df[is.na(broa_bon_df)] <- 0
  
  # normalize
  df_normalized <- prop.table(as.matrix(broa_bon_df[, c(3:ncol(broa_bon_df))]), margin = 1)
  df_normalized <- cbind(broa_bon_df[, c("x", "y")], df_normalized)
  df_normalized[is.na(df_normalized)] <- 0
  
  # join table from Brus map and Bonanella maps and multiply the Brus value with the proportions
  r_broa_df <- as.data.frame(r_broa, xy = T, na.rm = F)
  colnames(r_broa_df) <- c("x", "y", "broa")
  r_broa_df <- r_broa_df %>% mutate(broa = as.numeric(broa))
  r_broa_df[is.na(r_broa_df)] <- 0
  
  broa_df <- r_broa_df %>% 
    left_join(df_normalized, by = c("x", "y"))  %>%
    mutate(across(all_of(broa_lookup), ~ . * broa)) %>% 
    dplyr::select(-broa)
  
  broa_names <- strsplit(broa_lookup, "\\.")
  broa_names <- sapply(broa_names, function(x) {
    paste0(substr(x[[1]], 1, 2), substr(x[[2]], 1, 2))
  })
  
  colnames(broa_df) <- c("x", "y", broa_names)
  
  
  
  ## for conifers ----
  coni <- dom_veg[grepl("Conifers", dom_veg)]
  r_coni <- rast(coni)
  r_coni <- terra::crop(r_coni, height_rast)
  r_coni <- terra::project(r_coni, height_rast)
  r_coni <- terra::mask(r_coni, height_rast)
  
  # get bonanella maps for the coniferous
  coni_bon_rasts <- bon_rasts[grepl(paste(coni_lookup, collapse = "|"), bon_rasts)]
  coni_bon_rasts <- rast(coni_bon_rasts)
  names(coni_bon_rasts) <- coni_lookup
  
  # reproject to same epsg and res as height rast
  coni_bon_rasts <- terra::crop(coni_bon_rasts, height_rast_laea, snap = "out")
  coni_bon_rasts <- terra::project(coni_bon_rasts, height_rast, method = "average")
  coni_bon_rasts <- terra::mask(coni_bon_rasts, height_rast)
  
  # convert to dataframe
  coni_bon_df <- as.data.frame(coni_bon_rasts, xy = T, na.rm = F)
  coni_bon_df[is.na(coni_bon_df)] <- 0
  
  # normalize
  df_normalized <- prop.table(as.matrix(coni_bon_df[, c(3:ncol(coni_bon_df))]), margin = 1)
  df_normalized <- cbind(coni_bon_df[, c("x", "y")], df_normalized)
  df_normalized[is.na(df_normalized)] <- 0
  
  # join table from Brus map and Bonanella maps and multiply the Brus value with the proportions
  r_coni_df <- as.data.frame(r_coni, xy = T, na.rm = F)
  colnames(r_coni_df) <- c("x", "y", "coni")
  r_coni_df <- r_coni_df %>% mutate(coni = as.numeric(coni))
  r_coni_df[is.na(r_coni_df)] <- 0
  
  coni_df <- r_coni_df %>% 
    left_join(df_normalized, by = c("x", "y"))  %>%
    mutate(across(all_of(coni_lookup), ~ . * coni)) %>% 
    dplyr::select(-coni)
  
  coni_names <- strsplit(coni_lookup, "\\.")
  coni_names <- sapply(coni_names, function(x) {
    paste0(substr(x[[1]], 1, 2), substr(x[[2]], 1, 2))
  })
  
  colnames(coni_df) <- c("x", "y", coni_names)
  
  
  
  ## for Oaks ----
  oaks <- dom_veg[grepl("QuercusMisc", dom_veg)]
  r_oaks <- rast(oaks)
  r_oaks <- terra::crop(r_oaks, height_rast)
  r_oaks <- terra::project(r_oaks, height_rast)
  r_oaks <- terra::mask(r_oaks, height_rast)
  
  # get bonanella maps for the oaksferous
  oaks_bon_rasts <- bon_rasts[grepl(paste(oaks_lookup, collapse = "|"), bon_rasts)]
  oaks_bon_rasts <- rast(oaks_bon_rasts)
  names(oaks_bon_rasts) <- oaks_lookup
  
  # reproject to same epsg and res as height rast
  oaks_bon_rasts <- terra::crop(oaks_bon_rasts, height_rast_laea, snap = "out")
  oaks_bon_rasts <- terra::project(oaks_bon_rasts, height_rast, method = "average")
  oaks_bon_rasts <- terra::mask(oaks_bon_rasts, height_rast)
  
  # convert to dataframe
  oaks_bon_df <- as.data.frame(oaks_bon_rasts, xy = T, na.rm = F)
  oaks_bon_df[is.na(oaks_bon_df)] <- 0
  
  # normalize
  df_normalized <- prop.table(as.matrix(oaks_bon_df[, c(3:ncol(oaks_bon_df))]), margin = 1)
  df_normalized <- cbind(oaks_bon_df[, c("x", "y")], df_normalized)
  df_normalized[is.na(df_normalized)] <- 0
  
  # join table from Brus map and Bonanella maps and multiply the Brus value with the proportions
  r_oaks_df <- as.data.frame(r_oaks, xy = T, na.rm = F)
  colnames(r_oaks_df) <- c("x", "y", "oaks")
  r_oaks_df <- r_oaks_df %>% mutate(oaks = as.numeric(oaks))
  r_oaks_df[is.na(r_oaks_df)] <- 0
  
  oaks_df <- r_oaks_df %>% 
    left_join(df_normalized, by = c("x", "y"))  %>%
    mutate(across(all_of(oaks_lookup), ~ . * oaks)) %>% 
    dplyr::select(-oaks)
  
  oaks_names <- strsplit(oaks_lookup, "\\.")
  oaks_names <- sapply(oaks_names, function(x) {
    paste0(substr(x[[1]], 1, 2), substr(x[[2]], 1, 2))
  })
  
  colnames(oaks_df) <- c("x", "y", oaks_names)
  
  
  
  ## for Pines ----
  pine <- dom_veg[grepl("PinesMisc", dom_veg)]
  r_pine <- rast(pine)
  r_pine <- terra::crop(r_pine, height_rast)
  r_pine <- terra::project(r_pine, height_rast)
  r_pine <- terra::mask(r_pine, height_rast)
  
  # get bonanella maps for the pineferous
  pine_bon_rasts <- bon_rasts[grepl(paste(pine_lookup, collapse = "|"), bon_rasts)]
  pine_bon_rasts <- rast(pine_bon_rasts)
  names(pine_bon_rasts) <- pine_lookup
  
  # reproject to same epsg and res as height rast
  pine_bon_rasts <- terra::crop(pine_bon_rasts, height_rast_laea, snap = "out")
  pine_bon_rasts <- terra::project(pine_bon_rasts, height_rast, method = "average")
  pine_bon_rasts <- terra::mask(pine_bon_rasts, height_rast)
  
  # convert to dataframe
  pine_bon_df <- as.data.frame(pine_bon_rasts, xy = T, na.rm = F)
  pine_bon_df[is.na(pine_bon_df)] <- 0
  
  # normalize
  df_normalized <- prop.table(as.matrix(pine_bon_df[, c(3:ncol(pine_bon_df))]), margin = 1)
  df_normalized <- cbind(pine_bon_df[, c("x", "y")], df_normalized)
  df_normalized[is.na(df_normalized)] <- 0
  
  # join table from Brus map and Bonanella maps and multiply the Brus value with the proportions
  r_pine_df <- as.data.frame(r_pine, xy = T, na.rm = F)
  colnames(r_pine_df) <- c("x", "y", "pine")
  r_pine_df <- r_pine_df %>% mutate(pine = as.numeric(pine))
  r_pine_df[is.na(r_pine_df)] <- 0
  
  pine_df <- r_pine_df %>% 
    left_join(df_normalized, by = c("x", "y"))  %>%
    mutate(across(all_of(pine_lookup), ~ . * pine)) %>% 
    dplyr::select(-pine)
  
  pine_names <- strsplit(pine_lookup, "\\.")
  pine_names <- sapply(pine_names, function(x) {
    paste0(substr(x[[1]], 1, 2), substr(x[[2]], 1, 2))
  })
  
  colnames(pine_df) <- c("x", "y", pine_names)
  
  
  ## combine all layers
  r_final <- r_df %>% 
    left_join(broa_df, by = c("x", "y")) %>% 
    left_join(coni_df, by = c("x", "y")) %>% 
    left_join(oaks_df, by = c("x", "y")) %>% 
    left_join(pine_df, by = c("x", "y"))
  
  # Function to sum columns with the same first four characters
  sum_columns_with_same_prefix <- function(data) {
    prefixes <- unique(substr(names(data), 1, 4))
    
    for (prefix in prefixes) {
      matching_columns <- grep(paste0("^", prefix), names(data))
      if(length(matching_columns) > 1){
        data[prefix] <- rowSums(data[, matching_columns])
      }
    }
    return(data)
  }
  
  # Apply the function to the data frame
  result_df <- sum_columns_with_same_prefix(r_final)
  result_df <- result_df[, -grep("\\.", names(result_df))]
  result_df <- result_df %>% dplyr::select(-broa, -coni, -oaks, -pine)
  
  result_raster <- rast(result_df)
  crs(result_raster) <- proj_wgs
  par(mfrow = c(1,2))
  plot(result_raster[[5]])
  plot(result_raster[[16]])
  
  
  # double check
  if(origin(result_raster)[1] != origin(height_rast)[1] & origin(result_raster)[2] != origin(height_rast)[2]){
    print("Warning: origin not matching!")
  }
  if(ext(result_raster) != ext(height_rast)){
    print("Warning: extent not matching!")
  }
  if(round(res(result_raster)[1], 8) != round(res(height_rast)[1], 8) & round(res(result_raster)[2], 8) != round(res(height_rast)[2], 8)){
    print("Warning: resolution not matching!")
  }
  
  
  
  # save raster
  terra::writeRaster(x = result_raster, filename = paste0(path, "/initial_states/vegetation/species_stack_", name, ".tif"), overwrite = T)
  gc()
  
  
}


tmpFiles(remove = T)


### end --------------
