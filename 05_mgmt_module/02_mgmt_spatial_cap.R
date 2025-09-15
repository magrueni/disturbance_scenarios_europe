# --------------------------------------------------------------
# Management Module Script 02: Management cap
# --------------------------------------------------------------
# Author: Marc Grünig
# Date: 2025-01-07
# 
# Description:
# This script analyzes forest disturbance data across a 100km grid in Europe,
# focusing on the mean annual harvested areas. It processes spatial data,
# integrates disturbance information with grid cells, and calculates key statistics 
# such as the total disturbed area per year and the mean harvested area per grid cell.
# The output maps highlight areas of highest mean harvest activity, that serves
# as cap for the management within SVD
#
# Disturbance patches can be downloaded here: https://zenodo.org/records/8202241
# --------------------------------------------------------------
  
  
 
### libraries ------------------------------------------------------------------
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)

# paths
path <- "/.../"



### Data Preprocessing ---------------------------------------------------------

# load spatial data
ref_grid_100km <- rast(paste0(path, "/07_reference_grids/reference_grid_100km.tif"))
gridid <- ref_grid_100km
gridid_vect <- as.polygons(ref_grid_100km)
gridid_df <- as.data.frame(gridid, xy = T)

proj_laea <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"


# load patches that were assigned to management
d <- list.files(
  paste0(path, "/04_disturbance_modules/mgmt_module/patches"), 
  pattern = "patches", 
  full.names = TRUE
) %>%
  map_df(~{
    country <- str_sub(basename(.), 9, nchar(basename(.)) - 4)
    data <- read_csv(.)
    data <- mutate(data, 
                   country = country,
                   patch = as.integer(patch),
                   year = as.integer(year),
                   x = as.double(x),
                   y = as.double(y),
                   area = as.double(area))
    return(data)
  })

head(d)


years <- c(1986:2020)

# make spatial
d_sf <- d %>% 
  st_as_sf(
    .,
    coords = c("x", "y"),
    crs = st_crs(ref_grid_100km)
  )


# get dist per gridcells 100km per year
k <- 0
dat_annual <- vector("list", length(years))

# loop over years 
for (y in years) {
  
  print(y)
  
  k <- k + 1
  
  gridid_sel <- terra::extract(gridid, d_sf %>% filter(year == y)) %>%
    .$gridid
  gridid_sel <- as.data.frame(cbind(d_sf %>% filter(year == y), gridid = gridid_sel)) %>% 
    dplyr::select(patch, year, area, gridid)
  
  dat_annual[[k]] <- gridid_sel
  
}

dat_annual_summary <- dat_annual %>%
  bind_rows()


# set graphic
par(mfrow = c(1,2))

# plot the management over the years
plot_df <- dat_annual_summary %>% group_by(year) %>% summarize(sum_area = sum(area * 0.0001))
ggplot(plot_df, aes(x = year, y = sum_area))+
  geom_line()


# check the mean disturbance activity year for each pixel
mean_dist <- dat_annual_summary %>%
  group_by(gridid, year) %>%
  summarize(area_tot = sum(area, na.rm = T)) %>%
  ungroup()  %>%
  group_by(gridid) %>%
  mutate(area_tot = mean(area_tot, na.rm = T)) %>% 
  ungroup()
mean_dist <- mean_dist %>% mutate(area_tot = area_tot * 0.0001)

# convert to raster
grid_new <- left_join(gridid_df, mean_dist, by = c("gridid"))
new_rast_mean <- rast(grid_new[, c("x", "y", "area_tot")])
plot(new_rast_mean, main = "maximum harvest area")
crs(new_rast_mean) <- "epsg:3035"

# save
terra::writeRaster(new_rast_mean, paste0(path, "/04_disturbance_modules/mgmt_module/mean_harvest_100km.tif"), overwrite = T, datatype = "FLT4S")


