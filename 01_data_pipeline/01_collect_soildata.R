# ------------------------------------------------------------------------------------------------
# Script for Cleaning and Processing Soil Data
#
# Author: Marc Grünig
# 
# Description:
# This script performs a comprehensive workflow for cleaning, processing, and analyzing soil data 
# across Europe. It includes steps for loading and projecting geospatial datasets, masking datasets 
# to the European extent, filling missing data through spatial extraction, and calculating key 
# soil properties. The script also integrates soil information into simulation metadata.
#
# Key Components:
# 1. Set up paths and define geographic extents.
# 2. Load and process soil datasets, including texture, depth, and nitrogen availability.
# 3. Fill missing soil property values in simulation metadata using spatial extraction and calculations.
# 4. Write cleaned and enriched metadata back to a database.
#
# Outputs:
# - Cleaned and enriched soil metadata stored in a SQLite database.
# - Diagnostic histograms comparing original and cleaned metadata distributions.
#
# Dependencies:
# - R packages: `raster`, `terra`, `RColorBrewer`, `sf`, `dplyr`, `DBI`, `stars`, `ggplot2`, `exactextractr`.
# - Required input data:
#   - Soil and texture raster datasets.
#   - European extent shapefiles.
#   - Simulation metadata stored in a SQLite database.
#
# Notes:
# - Ensure all input datasets and paths are correctly specified in the `public_pth` variable.
# - Projections are standardized for consistency across all processing steps.
#
# -----------------------------------------------------------------------------------------------



### Load Required Libraries ---------------------------------------------------------------------
library(raster)
library(terra)
library(RColorBrewer)
library(sf)
library(dplyr)
library(DBI)
library(stars)
library(ggplot2)
library(exactextractr)


### Define Paths and Constants
public_pth <- "/.../"
eu_ext <- c(-10, 50, 35, 80)
proj_longlat <- "+proj=longlat +datum=WGS84 +no_defs"
proj_forest <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

# Load European boundary shapefiles
eu_poly <- vect(paste0(public_pth, "/06_06_gis/europe.shp"))
eu_poly_lr <- vect(paste0(public_pth, "/06_06_gis/europe_lowres.shp"))



### Load and Process Soil Data -----------------------------------------------------------------------------


# we cannot provide all of those datasets, but see here:
# https://esdac.jrc.ec.europa.eu/content/european-soil-database-derived-data

# Soil depth data
depth <- rast(paste0(public_pth, "/06_gis/STU/STU_EU_DEPTH_ROOTS.rst"))
crs(depth) <- proj_forest
depth <- depth %>% 
  terra::project(proj_forest) %>% 
  terra::mask(eu_poly_lr) %>% 
  terra::project(proj_longlat)
depth[depth == 0] <- NA
depth <- depth * 10 # Convert to mm

# Texture (clay, silt, sand)
load_texture <- function(layer_name) {
  layer <- rast(paste0(public_pth, "/06_gis/STU/", layer_name))
  crs(layer) <- proj_forest
  layer <- layer %>% 
    terra::project(proj_forest) %>% 
    terra::mask(eu_poly_lr) %>% 
    terra::project(proj_longlat)
  layer[layer == 0] <- NA
  return(layer)
}

clay_wgs <- load_texture("STU_EU_T_CLAY.rst")
silt_wgs <- load_texture("STU_EU_T_SILT.rst")
sand_wgs <- load_texture("STU_EU_T_SAND.rst")

# Water Holding Capacity (WHC)
twhc_s <- rast(paste0(public_pth, "/06_gis/STU/STU_EU_S_TAWC.rst"))
twhc_t <- rast(paste0(public_pth, "/06_gis/STU/STU_EU_T_TAWC.rst"))
whc <- twhc_s + twhc_t

crs(whc) <- proj_forest
whc <- whc %>% 
  terra::project(proj_forest) %>% 
  terra::mask(eu_poly_lr) %>% 
  terra::project(proj_longlat)
whc[whc == 0] <- NA

# Texture catergory layer
r_texture <- rast(paste0(public_pth, "/06_gis/texture_categories.tif"))
crs(r_texture) <- proj_forest
r_texture <- terra::project(r_texture, proj_forest, method = "mode")
r_texture[r_texture == 0] <- NA
r_texture <- terra::as.factor(r_texture)

# Nitrogen availability
nitro <- rast(paste0(public_pth, "/06_gis/predicted_N_available_1km_v2_interpolated.tif"))
nitro <- nitro %>% 
  terra::project(proj_longlat) %>% 
  terra::mask(eu_poly_lr)
nitro[nitro > 120] <- NA




### look at the simulations -----------------------------------
simulation_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(path, "/forest_simulation_db_v1.1.sqlite"))
tables_con <- dbListTables(simulation_db)
unique(tables_con)

metadata <- dbReadTable(simulation_db, "metadata_all")
# metadata <- read_csv(paste0(path_expls, "/data_pipeline/metadata_expl.csv"))

# some users used cm instead of mm
cm_users <- c(1003, 1021, 1023, 1032, 1033, 1034, 1035, 1036)
metadata <- metadata %>% 
  mutate(SoilDepth = ifelse(uniqueID %in% cm_users, SoilDepth * 10, SoilDepth))

# if there is no texture and no soildepth in the data we extract it from the maps
# sand silt clay depth, if one is NA, extract for all the values
if(nrow(metadata[is.na(metadata$TextureSand),]) > 0){
  
  metadata <- as.data.frame(cbind(metadata, 
                                  sand_extract = terra::extract(sand_wgs, metadata[, c("Lon", "Lat")])[,2],
                                  silt_extract = terra::extract(silt_wgs, metadata[, c("Lon", "Lat")])[,2],
                                  clay_extract = terra::extract(clay_wgs, metadata[, c("Lon", "Lat")])[,2],
                                  depth_extract = terra::extract(depth, metadata[, c("Lon", "Lat")])[,2]))
  
  metadata <- metadata %>% mutate(TextureSand = ifelse(is.na(TextureSand), sand_extract, TextureSand),
                                  TextureSilt = ifelse(is.na(TextureSilt), silt_extract, TextureSilt),
                                  TextureClay = ifelse(is.na(TextureClay), clay_extract, TextureClay),
                                  SoilDepth = ifelse(is.na(SoilDepth), depth_extract, SoilDepth))
  
  
  if(nrow(metadata[is.na(metadata$SoilDepth), ]) > 0){
    
    nas <- metadata[is.na(metadata$SoilDepth), ]
    coords <- nas[, c("Lon", "Lat")]
    coordinates(coords) <- ~Lon + Lat
    coords <- st_as_sf(coords)
    st_crs(coords) <- proj
    sp_buffer <- st_buffer(coords, dist = 1000) 
    
    # extract the data
    nas <- as.data.frame(cbind(nas, depth_buffer = exact_extract(depth, sp_buffer, fun = "mean", progress = F)))
    metadata[is.na(metadata$SoilDepth), "SoilDepth"] <- as.numeric(nas[, "depth_buffer"])
    
  }
  
  
  if(nrow(metadata[is.na(metadata$SoilDepth), ]) > 0){
    
    nas <- metadata[is.na(metadata$SoilDepth), ]
    coords <- nas[, c("Lon", "Lat")]
    coordinates(coords) <- ~Lon + Lat
    coords <- st_as_sf(coords)
    st_crs(coords) <- proj
    sp_buffer <- st_buffer(coords, dist = 20000) 
    
    # extract the data
    nas <- as.data.frame(cbind(nas, depth_buffer2 = exact_extract(depth, sp_buffer, fun = "mean", progress = F)))
    metadata[is.na(metadata$SoilDepth), "SoilDepth"] <- as.numeric(nas[, "depth_buffer2"])
    
  }
  
  
  # same for the other textures
  if(nrow(metadata[is.na(metadata$TextureSand), ]) > 0){
    
    nas <- metadata[is.na(metadata$TextureSand), ]
    coords <- nas[, c("Lon", "Lat")]
    coordinates(coords) <- ~Lon + Lat
    coords <- st_as_sf(coords)
    st_crs(coords) <- proj
    sp_buffer <- st_buffer(coords, dist = 1000) 
    
    # extract the data
    nas <- as.data.frame(cbind(nas, 
                               sand_buffer = exact_extract(sand_wgs, sp_buffer, fun = "mean", progress = F) ,
                               silt_buffer = exact_extract(silt_wgs, sp_buffer, fun = "mean", progress = F) ,
                               clay_buffer = exact_extract(clay_wgs, sp_buffer, fun = "mean", progress = F) ))
    
    metadata[is.na(metadata$TextureSand), "TextureSand"] <- as.numeric(nas[, "sand_buffer"])
    metadata[is.na(metadata$TextureSilt), "TextureSilt"] <- as.numeric(nas[, "silt_buffer"])
    metadata[is.na(metadata$TextureClay), "TextureClay"] <- as.numeric(nas[, "clay_buffer"])
    
  }
  
  if(nrow(metadata[is.na(metadata$TextureSand), ]) > 0){
    
    nas <- metadata[is.na(metadata$TextureSand), ]
    coords <- nas[, c("Lon", "Lat")]
    coordinates(coords) <- ~Lon + Lat
    coords <- st_as_sf(coords)
    st_crs(coords) <- proj
    sp_buffer <- st_buffer(coords, dist = 20000) 
    
    # extract the data
    nas <- as.data.frame(cbind(nas, 
                               sand_buffer2 = exact_extract(sand_wgs, sp_buffer, fun = "mean", progress = F) ,
                               silt_buffer2 = exact_extract(silt_wgs, sp_buffer, fun = "mean", progress = F) ,
                               clay_buffer2 = exact_extract(clay_wgs, sp_buffer, fun = "mean", progress = F) ))
    
    metadata[is.na(metadata$TextureSand), "TextureSand"] <- as.numeric(nas[, "sand_buffer2"])
    metadata[is.na(metadata$TextureSilt), "TextureSilt"] <- as.numeric(nas[, "silt_buffer2"])
    metadata[is.na(metadata$TextureClay), "TextureClay"] <- as.numeric(nas[, "clay_buffer2"])
    
  }
  
  
}


## WHC calculation
# Werner Rammer, April 2019

# Calculate water holding capacity following the approach of Schwalm & Ek (2004, Ecol. Mod.)
# pct_sand, pct_silt, pct_clay: soil texture fractions in %
# soil_depth: stone-free effective soil depth (mm)
# fc_threshold: field capacity (J kg-1 or kPa), default: 15kPa
# pwp_threshold: permanent wilting point (J kg-1 or kPa), default 1500kPa
# return: water holding capacity (between fc and pwp) in mm
calcWHC <- function(pct_sand, pct_silt, pct_clay, soil_depth, fc_threshold=15, pwp_threshold=4000) {
  
  theta_sat = 0.01 * (50.5 - 0.142*pct_sand - 0.037*pct_clay); # Eq. 78
  bt <- 11.43 - 0.103*pct_sand - 0.0687*pct_silt # Eq. 79
  rho_e <- -5.99 + 0.0544*pct_sand + 0.0451*pct_silt # Eq 80
  
  
  fc <- theta_sat * (rho_e / -fc_threshold)^(1/bt) * soil_depth # Eq 76
  pwp <- theta_sat * (rho_e / -pwp_threshold)^(1/bt) * soil_depth # Eq 77
  whc <- fc-pwp
  whc
}


# check which simulations have no data for WHC but a rating
# rating is from WHC 30mm to 200mm
if(nrow(metadata[is.na(metadata$WHC),]) > 0){
  
  metadata <- metadata %>% mutate(WHC = ifelse(is.na(WHC), SoilWaterRating*200, WHC))
  
  metadata <- metadata %>% mutate(WHC = ifelse(is.na(WHC), 
                                               calcWHC(pct_sand = TextureSand,
                                                       pct_silt = TextureSilt,
                                                       pct_clay = TextureClay,
                                                       soil_depth = SoilDepth),
                                               WHC))
  
  metadata <- as.data.frame(cbind(metadata, 
                                  whc_extract = terra::extract(r_whc_wgs, metadata[, c("Lon", "Lat")])[,2]))
  
  metadata <- metadata %>% mutate(WHC = ifelse(is.na(WHC) & is.na(SoilDepth), whc_extract, WHC))
  
  if(nrow(metadata[is.na(metadata$WHC), ]) > 0){
    
    nas <- metadata[is.na(metadata$WHC), ]
    coords <- nas[, c("Lon", "Lat")]
    coordinates(coords) <- ~Lon + Lat
    coords <- st_as_sf(coords)
    st_crs(coords) <- proj
    sp_buffer <- st_buffer(coords, dist = 1000) 
    
    # extract the data
    nas <- as.data.frame(cbind(nas, whc_buffer = exact_extract(r_whc_wgs, sp_buffer, fun = "mean", progress = F)))
    
    #nas <- nas %>% mutate(WHC = ifelse(is.na(WHC), whc_buffer, WHC))
    metadata[is.na(metadata$WHC), "WHC"] <- as.numeric(nas[, "whc_buffer"])
    
  }
}



# if available nitrogen is missing we take the fertiliy rating
# for some cases there is no rating, therefore we extract values from the map
if(nrow(metadata[is.na(metadata$AvailableNitrogen), ]) > 0){
  
  # transform the rating value to what we assigned them in the model data interface
  metadata <- metadata %>% mutate(AvailableNitrogen = ifelse(is.na(AvailableNitrogen), 
                                                             ifelse(FertilityRating < 0.5, FertilityRating*80 + 20*(1-FertilityRating), FertilityRating*100), 
                                                             AvailableNitrogen))
  metadata <- as.data.frame(cbind(metadata, 
                                  nitro_extract = terra::extract(nitro_wgs, metadata[, c("Lon", "Lat")])[,2]))
  
  metadata <- metadata %>% mutate(AvailableNitrogen = ifelse(is.na(AvailableNitrogen), nitro_extract, AvailableNitrogen))
  
  # there is a handful (12) simulations in two locations that are on the coast in northern germany and there get a NA from the extract function, we make correction by adding a 1km buffer around the locations
  if(nrow(metadata[is.na(metadata$AvailableNitrogen), ]) > 0){
    
    nas <- metadata[is.na(metadata$AvailableNitrogen), ]
    coords <- nas[, c("Lon", "Lat")]
    coordinates(coords) <- ~Lon + Lat
    coords <- st_as_sf(coords)
    st_crs(coords) <- proj
    sp_buffer <- st_buffer(coords, dist = 1000) 
    
    # extract the data
    nas <- as.data.frame(cbind(nas, nitro_buffer = exact_extract(nitro_wgs, sp_buffer, fun = "mean", progress = F)))
    metadata[is.na(metadata$AvailableNitrogen), "AvailableNitrogen"] <- as.numeric(nas[, "nitro_buffer"])
    
  }
  
  # there is a handful (12) simulations in two locations that are on the coast in northern germany and there get a NA from the extract function, we make correction by adding a 1km buffer around the locations
  if(nrow(metadata[is.na(metadata$AvailableNitrogen), ]) > 0){
    
    nas <- metadata[is.na(metadata$AvailableNitrogen), ]
    coords <- nas[, c("Lon", "Lat")]
    coordinates(coords) <- ~Lon + Lat
    coords <- st_as_sf(coords)
    st_crs(coords) <- proj
    sp_buffer <- st_buffer(coords, dist = 1500) 
    
    # extract the data
    nas <- as.data.frame(cbind(nas, nitro_buffer2 = exact_extract(nitro_wgs, sp_buffer, fun = "mean", progress = F)))
    metadata[is.na(metadata$AvailableNitrogen), "AvailableNitrogen"] <- as.numeric(nas[, "nitro_buffer2"])
    
  }
  
}

# check if there are NAs
# nrow(metadata[is.na(metadata$AvailableNitrogen) | is.na(metadata$WHC) | is.na(metadata$SoilDepth) | is.na(metadata$TextureSand), ])

# save
metadata <- metadata[, 1:19]
dbWriteTable(conn = simulation_db, name = "metadata_all_cleaned_soil", value = metadata, overwrite = T)
dbDisconnect(simulation_db)



### add soil conditions to examples ------------------------------------------------------------------------------

path <- "/.../"


simulation_db <- DBI::dbConnect(RSQLite::SQLite(), paste0(path, "/forest_simulation_db_v1.1.sqlite"))
tables_con <- dbListTables(simulation_db)
unique(tables_con)

metadata <- dbReadTable(simulation_db, "metadata_all_cleaned_soil" )
examples <- dbReadTable(simulation_db, "examples" )


# for example dataset:
# metadata <- read_csv(paste0(path_expls, "/data_pipeline/metadata_expl.csv"))
# examples <- read_csv(paste0(path_expls, "/data_pipeline/examples_expl.csv"))

examples_with_soil <- examples %>% left_join(metadata %>% dplyr::select(uniqueID, ID, WHC, TextureSand, SoilDepth, AvailableNitrogen) %>% 
                                                 mutate(simulationID = ID) %>% 
                                                 dplyr::select(-ID), by = c("uniqueID", "simulationID")) %>% 
    dplyr::select(-flag) %>% 
    mutate(uniqueID = as.integer(uniqueID),
           simulationID = as.integer(simulationID),
           Year = as.integer(Year), 
           residence_time = as.integer(residence_time),
           target_time = as.integer(target_time))


dbWriteTable(conn = simulation_db, name = "examples_2m_soil", value = examples_with_soil, overwrite = T)

dbDisconnect(simulation_db)


