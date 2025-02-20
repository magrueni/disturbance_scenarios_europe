# --------------------------------------------------------------
# Initial landscape Scipt 03: Create initial soil conditions layer
# --------------------------------------------------------------
# Author: Marc Grünig
# Date: 2025-01-07
# 
#
# Description:
# This script creates a soil layer for Europe based on different data sources.
# 
#
# Notes:
# - soil data was obtained from the ESDAC dataset
# - available nitrogen layer was created using a spatial modelling approach (see Grünig et al., 2024)
#
# --------------------------------------------------------------

library(raster)
library(terra)
library(RColorBrewer)
library(sf)
library(dplyr)
library(DBI)
library(stars)
library(ggplot2)
library(exactextractr)


path <- "/.../"



### gis load -------------------------------------------------------------------
eu.ext <- c(-10,50,35,80)
eu_poly <- vect(paste0(path, "/3pg/gis/europe.shp"))
eu_poly_lr <- vect(paste0(path, "/3pg/gis/europe_lowres.shp"))

proj_forest <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
proj <- "+proj=longlat +datum=WGS84 +no_defs"
proj_new <- "+proj=tmerc +lat_0=0 +lon_0=10.3333333333333 +k=1 +x_0=0 +y_0=-5000000 +ellps=bessel +units=m +no_defs +type=crs"
proj_laea <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"


### load data -------------------------------------------------------------------

# soil depth
depth <- rast(paste0(path, "/soil/esdac/STU/STU_EU_DEPTH_ROOTS.rst"))
crs(depth) <- proj_forest
depth <- terra::project(depth, proj_forest)
depth <- terra::mask(depth, eu_poly_lr)
depth <- terra::project(depth, proj_laea)
depth <- depth * 10
origin(depth) <- c(0, 0)

# sand content
sand_s <- rast(paste0(path, "/soil/esdac/STU/STU_EU_T_SAND.rst"))
crs(sand_s) <- proj_forest
sand_s <- terra::project(sand_s, proj_forest)
sand_s <- terra::mask(sand_s, eu_poly_lr)
sand_wgs <- terra::project(sand_s, proj_laea)
origin(sand_wgs) <- c(0, 0)

# water holding capacity
twhc_s <- rast(paste0(path, "/soil/esdac/STU/STU_EU_S_TAWC.rst"))
twhc_t <- rast(paste0(path, "/soil/esdac/STU/STU_EU_T_TAWC.rst"))
r_whc_wgs <- twhc_t + twhc_s

crs(r_whc_wgs) <- proj_forest
r_whc_wgs <- terra::project(r_whc_wgs, proj_forest)
r_whc_wgs <- terra::mask(r_whc_wgs, eu_poly_lr)
r_whc_wgs[r_whc_wgs == 0] <- NA
r_whc_wgs <- terra::project(r_whc_wgs, proj_laea)
r_whc_wgs <- terra::focal(r_whc_wgs, w = c(9, 9), fun = "mean", na.rm = T, na.policy = "only")
origin(r_whc_wgs) <- c(0, 0)

# load nitrogen
r_nitro <- rast(paste0(path, "/gis/predicted_N_available.tif"))
r_nitro <- terra::project(r_nitro, depth)
r_nitro <- terra::mask(r_nitro, depth)
r_nitro[r_nitro > 120] <- NA
nitro_wgs <- terra::project(r_nitro, proj_laea)
nitro_wgs <- terra::focal(nitro_wgs, w = c(9, 9), fun = "mean", na.rm = T, na.policy = "only")
nitro_wgs <- terra::focal(nitro_wgs, w = c(17, 17), fun = "mean", na.rm = T, na.policy = "only")
nitro_wgs <- terra::focal(nitro_wgs, w = c(25, 25), fun = "mean", na.rm = T, na.policy = "only")
nitro_wgs <- terra::focal(nitro_wgs, w = c(35, 35), fun = "mean", na.rm = T, na.policy = "only")
origin(nitro_wgs) <- c(0, 0)
rm(r_nitro)


# create a new soil raster
soil_rast <- depth
values(soil_rast) <- 1:ncell(soil_rast)
soil_rast <- terra::mask(soil_rast, depth)

# climate raster
x_rast <- rast(paste0(path, "/reference_grids/reference_grid.tif"))
crs(x_rast) <- "+proj=longlat +datum=WGS84 +no_defs"
climate_rast <- terra::project(climate_rast, proj_laea, "near")
origin(climate_rast) <- c(0, 0)


# for each grid cell of the soil raster extract the values for all layers
CellsXY <- as.data.frame(soil_rast, xy = T)
data_df <- as.data.frame(cbind(CellsXY[, c("x", "y")],
                               climateId = terra::extract(climate_rast, CellsXY[, c("x", "y")])[,2],
                               depth = terra::extract(depth, CellsXY[, c("x", "y")])[,2],
                               sand = terra::extract(sand_wgs, CellsXY[, c("x", "y")])[,2],
                               whc = terra::extract(r_whc_wgs, CellsXY[, c("x", "y")])[,2],
                               nitro = terra::extract(nitro_wgs, CellsXY[, c("x", "y")])[,2]
))


data_df <- cbind(id = c(1:nrow(data_df)),
                 data_df[, c("x", "y", "climateId", "depth",
                             "sand", "whc", "nitro")])

colnames(data_df) <- c("id", "x", "y", "climateId",
                       "soil_depth", "sand_pct", "whc",
                       "available_nitrogen")

# 
env_rast <- rast(data_df[, c("x", "y", "id")])
env_rast[is.na(env_rast)] <- 0
crs(env_rast) <- proj_laea
env_rast <- terra::project(env_rast, "epsg:3035", method = "near")

# # forest mask
forest_mask <- rast(psate0(path, "/gis/forest_mask_new_laea_masked.tif"))

# change to forest mask resolution
env_rast <- terra::project(env_rast, forest_mask, "near")
forest_mask <- terra::mask(forest_mask, env_rast)
env_rast <- terra::mask(env_rast, forest_mask)
env_rast[env_rast == 0] <- NA


env_df <- as.data.frame(env_rast, xy = T)


writeRaster(env_rast,
            paste0(path, "/initial_states/initial_soil_rast.tif"),
            overwrite = T)


env_rast <- rast(paste0(path, "/initial_states/initial_soil_rast.tif"))
env_df <- as.data.frame(env_rast, xy = T, na.rm = T)


data_df <- data_df[, c("id", "climateId", "soil_depth",
                       "sand_pct", "whc", "available_nitrogen")]


data_df_new <- left_join(env_df, data_df, by = c("id")) %>% 
  dplyr::select(id, climateId, soil_depth, sand_pct, whc, available_nitrogen) %>% 
  mutate(soil_depth = soil_depth * 10) %>% 
  distinct()

head(data_df_new)
nrow(data_df_new[is.na(data_df_new$climateId),])

write.csv(data_df_new,
          paste0(path, "/initial_states/initial_soil.csv"),
          row.names = F)



### load initial vegetation file and aggregate to that resolution -------------------------------
init_soil_100 <- rast(paste0(path, "/initial_states/initial_soil_rast.tif"))

env_df <- as.data.frame(init_soil_100, xy = T)


init_veg_100 <- rast(paste0(path, "/initial_states/init_veg.tif"))

# convert soil to the resolution of the veg if not correct
init_soil_100 <- terra::resample(env_rast, init_veg_100, method="mode")
init_soil_100 <- terra::crop(init_soil_100, init_veg_100)  
init_soil_100 <- crop(extend(init_soil_100, init_veg_100), init_veg_100)

# the grid cells without veg are put to NAs and then masked in the soil raster
init_soil_100 <- terra::mask(init_soil_100, init_veg_100)

writeRaster(init_soil_100, 
            paste0(path, "/initial_states/init_soil.tif"),
            overwrite = T,
            datatype="INT4S")

writeRaster(init_veg_100, 
            paste0(path, "/initial_states/init_veg_v2.tif"),
            datatype="INT2S", overwrite=T)  


# if some IDs are missing in the CSV file:
init_soil_100 <- rast(paste0(path, "/initial_states/init_soil.tif"))
soil_vals <- unique(init_soil_100)
init_soil_dat <- read.csv(paste0(path, "/initial_states/initial_soil.csv"))


length(unique(init_soil_dat$id))
length(unique(soil_vals$id))

# check for missing IDs
misses <- soil_vals$id[ !(soil_vals$id %in% init_soil_dat$id) ]
length(misses)

gc()
tmpFiles(remove = T)


### end -------------------------------------------------------------------------
