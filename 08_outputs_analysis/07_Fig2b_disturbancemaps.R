# ------------------------------------------------------------------------------
# Script Name: Plotting disturbance rates
# Description: This script combines spatial simulation output and creates maps
# of disturbance hotspots. 
#
# Author: Marc Grünig
# Date: 7.2.2025
#
# Input Files:
#   - outputs from script 3 disturbance rates at hexagon level
#
# Output:
#   - Figure 2b, Figure S2 - S10
# ------------------------------------------------------------------------------


### libraries ------------------------------------------------------------------
library(rnaturalearth)
library(rnaturalearthdata)
library(MetBrewer)
library(RColorBrewer)
library(tidyterra)
library(gridExtra)
library(raster)
library(ggplot2)
library(sf)
library(terra)
library(stringr)
library(readr)
library(sfheaders)
library(dplyr)
library(tidyverse)



### load data ------------------------------------------------------------------

# define path
path <- "/.../"
path_results <- paste0(path, "/svd_simulations/results_eu/")

# load shp and hexagon
eu_shp <- vect(paste0(path, "/reference_grids/eu_mask.shp"))
hex_ecu <- shapefile(paste0(path, "/reference_grids/eu_mask_hexagons_25km.shp"))


# load forest mask
sf_obj_forest_mask <- read_sf(paste0(path, "/reference_grids/hex_forest_mask_25km.gpkg"))
sf_obj_forest_mask_ids <- sf_obj_forest_mask %>% st_drop_geometry()


# load simulation overview to assign sims to rcps
master_tab <- read.csv(paste0(path, "/svd_simulations/svd_simulations_ids.csv"), sep = ";")


# define rcps
rcps <- c("rcp85", "rcp45", "rcp26")


# create output and counter
r <- 0
all_maps_dat <- list()


# loop over rcps
for(c in 1:length(rcps)){
  
  # progress
  print(c)
  
  # loop over years
  for(y in seq(10, 80, 10)){
    
    r <- r + 1
    
    dat <- read_csv(paste0(path_results, "/dist_rates_25km/all_dist_rates_", y, "_", rcps[c], ".csv"))
    
    all_maps_dat[[r]] <- dat 
    
  }
  
}


# combine outputs
all_dat <- do.call(rbind, all_maps_dat)  


# get the size of one hexagon
hex_size <- max(sf_obj_forest_mask$forest_area, na.rm = T)


# get europe map for plotting
proj_wgs <- "+proj=longlat +datum=WGS84 +no_defs"
proj_leae <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
eu_shp <- terra::project(eu_shp, proj_wgs)
eu_shp_sf <- st_as_sf(eu_shp)
eu_shp_sf_leae <- st_transform(eu_shp_sf, proj_leae)
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, proj_leae)
europe <- st_crop(world, eu_shp_sf_leae)
europe <- europe[!europe$sovereignt %in%
                   c("Tunisia", "Russia", "Algeria", "Ukraine", "Belarus", "Moldova",
                     "Turkey", "Cyprus", "Northern Cyprus", "Morocco", "Iceland"), ]





### plotting ---------------------------------------------------------------------------------------

## first RCP8.5 -------------------------------------------------------------------------------------

# filter the data
dat_85 <- all_dat %>% filter(rcp == "rcp85")

# prepare data and calculate dist rate per hexagon
dat_85_all <- dat_85 %>% dplyr::select(-sim) %>% 
  group_by(gridid, rcp, gcm, rep, year) %>% 
  summarise(dist_area = sum(dist_area, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>% 
  mutate(dist_rate = as.numeric(100 / forest_area * dist_area / 10),
         dist_freq = as.numeric(10 * forest_area / dist_area)) %>% 
  group_by(gridid, year) %>% 
  summarize(mean_dist_rate = mean(dist_rate, na.rm = T),
            sd_dist_rate = sd(dist_rate, na.rm = T),
            mean_dist_freq = mean(dist_freq, na.rm = T),
            sd_dist_freq = sd(dist_freq, na.rm = T),
            n = n()) %>% 
  ungroup() %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>%
  st_as_sf(.)



## then the mean across the whole period
dat_85_all_plot <- dat_85_all %>%
  mutate(mean_dist_rate = ifelse(is.nan(mean_dist_rate), NA, mean_dist_rate), 
         mean_dist_freq = ifelse(is.nan(mean_dist_freq), NA, mean_dist_freq)) %>% 
  group_by(gridid) %>% 
  summarize(mean_dist_rate = mean(mean_dist_rate, na.rm = T),
            land_area = mean(land_area, na.rm = T),
            forest_area = mean(forest_area, na.rm = T)) %>% 
  mutate(forestShare = (1/land_area * forest_area)) %>% 
  mutate(mean_dist_rate = ifelse(mean_dist_rate > 0.5, 0.5, mean_dist_rate))



centroids <- st_centroid(dat_85_all_plot)
centroid_coords <- st_coordinates(centroids)

# Add the centroid coordinates to the original data frame
sf_with_centroids <- dat_85_all_plot %>%
  mutate(
    x = centroid_coords[, 1],  # X coordinate of centroid
    y = centroid_coords[, 2]   # Y coordinate of centroid
  )

dat_85_all_plot_mean <- sf_with_centroids

cols <- rev(met.brewer("OKeeffe2", 100, "continuous"))
DistMap <- ggplot() +
  geom_point(data=dat_85_all_plot_mean,
             aes(x=x, y=y, color=mean_dist_rate), size=0)+
  geom_point(data=dat_85_all_plot_mean,
             aes(x=x, y=y, color=mean_dist_rate, size=forestShare), show.legend = TRUE)+
  geom_spatvector(data = eu_shp, fill = "black", col = "black", alpha = 0.1) +
  geom_spatvector(data = europe, fill = "transparent", col = "black") +
  theme_classic() +
  scale_size_continuous(range = c(-0.5, 1.5)) + 
  scale_color_gradientn(colours = rev(cols), na.value = "lightgrey", limits = c(0, 0.5), breaks = c(0, 0.5)) +
  labs(color = bquote("disturbance \n rate [% yr"^-1*"]"),        size = bquote("forest share")) + # Change color legend title here
  theme(axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, vjust = 3),
        legend.position = "right") +
  xlab("") + ylab("")

DistMap
ggsave(DistMap, filename = paste0(path_results, "/figures/Fig2b.png"), width = 12, height = 9)
write_csv(dat_85_all_plot_mean, paste0(path_results, "/figures/plot_data/figure2b_data_rcp85.csv"))


# then the delta ------------
all_bas <- dat_85_all %>% filter(year == 2030) %>% st_drop_geometry() %>%
  mutate(mean_dist_rate = ifelse(mean_dist_rate == 0, 0.0001, mean_dist_rate))

all_fut <- dat_85_all %>% filter(year == 2100) %>% st_drop_geometry() %>%
  mutate(mean_dist_rate = ifelse(mean_dist_rate == 0, 0.0001, mean_dist_rate))

diff_df <- left_join(all_bas, all_fut, by = c("gridid")) %>%
  mutate(mean_diff = mean_dist_rate.y - mean_dist_rate.x,
         mean_diff_rel = mean_dist_rate.y / mean_dist_rate.x - 1,
         mean_diff_prct = (mean_dist_rate.y - mean_dist_rate.x) / mean_dist_rate.x * 100) %>%
  mutate(forestShare = (1/land_area.x * forest_area.x)) %>%
  mutate(across(c(mean_diff, mean_diff_rel, mean_diff_prct), ~ ifelse(is.finite(.), ., NA))) %>%
  left_join(., sf_obj_forest_mask, by = "gridid") %>%
  st_as_sf(.)


100/length(diff_df$mean_diff)*nrow(diff_df[diff_df$mean_diff > 0,])
100/length(diff_df$mean_diff)*nrow(diff_df[diff_df$mean_diff < 0,])


centroids <- st_centroid(diff_df)
centroid_coords <- st_coordinates(centroids)

# Add the centroid coordinates to the original data frame
sf_with_centroids <- diff_df %>%
  mutate(
    x = centroid_coords[, 1],  # X coordinate of centroid
    y = centroid_coords[, 2]   # Y coordinate of centroid
  )

dat_85_all_plot_diff <- sf_with_centroids

cols <- met.brewer("OKeeffe1", 100, "continuous")
DistMap <- ggplot() +
  geom_point(data=dat_85_all_plot_diff,
             aes(x=x, y=y, color=mean_diff), size=0)+
  geom_point(data=dat_85_all_plot_diff,
             aes(x=x, y=y, color=mean_diff, size=forestShare), show.legend = TRUE)+
  geom_spatvector(data = eu_shp, fill = "black", col = "black", alpha = 0.1) +
  geom_spatvector(data = europe, fill = "transparent", col = "black") +
  theme_classic() +
  scale_size_continuous(range = c(-0.5, 1.5)) +
  scale_color_gradientn(colours = rev(cols), na.value = "lightgrey", limits = c(-0.5, 0.5), breaks = c(-0.5, 0.5)) +
  labs(color = bquote("disturbance \n rate [% year"^-1*"]"),        size = bquote("forest share")) + # Change color legend title here
  theme(axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, vjust = 3),
        legend.position = "right") +
  xlab("") + ylab("")

DistMap
ggsave(DistMap, filename = paste0(path_results, "/figures/FigS4.png"), width = 12, height = 9)
write_csv(dat_85_all_plot_mean, paste0(path_results, "/figures/plot_data/figS4_data_map.csv"))



### RCP4.5 ---------------------------------------------------------------------

# filter the data
dat_45 <- all_dat %>% filter(rcp == "rcp45")

# prepare data and calculate dist rate per hexagon
dat_45_all <- dat_45 %>% dplyr::select(-sim) %>% 
  group_by(gridid, rcp, gcm, rep, year) %>% 
  summarise(dist_area = sum(dist_area, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>% 
  mutate(dist_rate = as.numeric(100 / forest_area * dist_area / 10),
         dist_freq = as.numeric(10 * forest_area / dist_area)) %>% 
  group_by(gridid, year) %>% 
  summarize(mean_dist_rate = mean(dist_rate, na.rm = T),
            sd_dist_rate = sd(dist_rate, na.rm = T),
            mean_dist_freq = mean(dist_freq, na.rm = T),
            sd_dist_freq = sd(dist_freq, na.rm = T),
            n = n()) %>% 
  ungroup() %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>%
  st_as_sf(.)

## then the mean across the whole period
dat_45_all_plot <- dat_45_all %>%
  mutate(mean_dist_rate = ifelse(is.nan(mean_dist_rate), NA, mean_dist_rate), 
         mean_dist_freq = ifelse(is.nan(mean_dist_freq), NA, mean_dist_freq)) %>% 
  group_by(gridid) %>% 
  summarize(mean_dist_rate = mean(mean_dist_rate, na.rm = T),
            land_area = mean(land_area, na.rm = T),
            forest_area = mean(forest_area, na.rm = T)) %>% 
  mutate(forestShare = (1/land_area * forest_area)) %>% 
  mutate(mean_dist_rate = ifelse(mean_dist_rate > 0.5, 0.5, mean_dist_rate))



centroids <- st_centroid(dat_45_all_plot)
centroid_coords <- st_coordinates(centroids)

# Add the centroid coordinates to the original data frame
sf_with_centroids <- dat_45_all_plot %>%
  mutate(
    x = centroid_coords[, 1],  # X coordinate of centroid
    y = centroid_coords[, 2]   # Y coordinate of centroid
  )

dat_45_all_plot_mean <- sf_with_centroids


cols <- rev(met.brewer("OKeeffe2", 100, "continuous"))
DistMap <- ggplot() +
  geom_point(data=dat_45_all_plot_mean,
             aes(x=x, y=y, color=mean_dist_rate), size=0)+
  geom_point(data=dat_45_all_plot_mean,
             aes(x=x, y=y, color=mean_dist_rate, size=forestShare), show.legend = TRUE)+
  geom_spatvector(data = eu_shp, fill = "black", col = "black", alpha = 0.1) +
  geom_spatvector(data = europe, fill = "transparent", col = "black") +
  theme_classic() +
  scale_size_continuous(range = c(-0.5, 1.5)) + 
  scale_color_gradientn(colours = rev(cols), na.value = "lightgrey", limits = c(0, 0.5), breaks = c(0, 0.5)) +
  labs(color = bquote("disturbance \n rate [% year"^-1*"]"),        size = bquote("forest share")) + # Change color legend title here
  theme(axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, vjust = 3),
        legend.position = "right") +
  xlab("") + ylab("")

DistMap
ggsave(DistMap, filename = paste0(path_results, "/figures/FigS2.png"), width = 12, height = 9)
write_csv(dat_45_all_plot_mean, paste0(path_results, "/figures/plot_data/figS2_data_map.csv"))



## then the delta ------------
all_bas <- dat_45_all %>% filter(year == 2030) %>% st_drop_geometry() %>%
  mutate(mean_dist_rate = ifelse(mean_dist_rate == 0, 0.0001, mean_dist_rate))

all_fut <- dat_45_all %>% filter(year == 2100) %>% st_drop_geometry() %>%
  mutate(mean_dist_rate = ifelse(mean_dist_rate == 0, 0.0001, mean_dist_rate))

diff_df <- left_join(all_bas, all_fut, by = c("gridid")) %>%
  mutate(mean_diff = mean_dist_rate.y - mean_dist_rate.x,
         mean_diff_rel = mean_dist_rate.y / mean_dist_rate.x - 1,
         mean_diff_prct = (mean_dist_rate.y - mean_dist_rate.x) / mean_dist_rate.x * 100) %>%
  mutate(forestShare = (1/land_area.x * forest_area.x)) %>%
  mutate(across(c(mean_diff, mean_diff_rel, mean_diff_prct), ~ ifelse(is.finite(.), ., NA))) %>%
  left_join(., sf_obj_forest_mask, by = "gridid") %>%
  st_as_sf(.)


centroids <- st_centroid(diff_df)
centroid_coords <- st_coordinates(centroids)

# Add the centroid coordinates to the original data frame
sf_with_centroids <- diff_df %>%
  mutate(
    x = centroid_coords[, 1],  # X coordinate of centroid
    y = centroid_coords[, 2]   # Y coordinate of centroid
  )

dat_45_all_plot_diff <- sf_with_centroids

cols <- met.brewer("OKeeffe1", 100, "continuous")
DistMap <- ggplot() +
  geom_point(data=dat_45_all_plot_diff,
             aes(x=x, y=y, color=mean_diff), size=0)+
  geom_point(data=dat_45_all_plot_diff,
             aes(x=x, y=y, color=mean_diff, size=forestShare), show.legend = TRUE)+
  geom_spatvector(data = eu_shp, fill = "black", col = "black", alpha = 0.1) +
  geom_spatvector(data = europe, fill = "transparent", col = "black") +
  theme_classic() +
  scale_size_continuous(range = c(-0.5, 1.5)) +
  scale_color_gradientn(colours = rev(cols), na.value = "lightgrey", limits = c(-0.5, 0.5), breaks = c(-0.5, 0.5)) +
  labs(color = bquote("disturbance \n rate [% year"^-1*"]"),        size = bquote("forest share")) + # Change color legend title here
  theme(axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, vjust = 3),
        legend.position = "right") +
  xlab("") + ylab("")

DistMap
ggsave(DistMap, filename = paste0(path_results, "/figures/FigS5.png"), width = 12, height = 9)
write_csv(dat_45_all_plot_mean, paste0(path_results, "/figures/plot_data/figS5_data_map.csv"))



## RCP2.6 --------------------------------------------------------------------------------------------------

# filter the data
dat_26 <- all_dat %>% filter(rcp == "rcp26")

# prepare data and calculate dist rate per hexagon
dat_26_all <- dat_26 %>% dplyr::select(-sim) %>% 
  group_by(gridid, rcp, gcm, rep, year) %>% 
  summarise(dist_area = sum(dist_area, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>% 
  mutate(dist_rate = as.numeric(100 / forest_area * dist_area / 10),
         dist_freq = as.numeric(10 * forest_area / dist_area)) %>% 
  group_by(gridid, year) %>% 
  summarize(mean_dist_rate = mean(dist_rate, na.rm = T),
            sd_dist_rate = sd(dist_rate, na.rm = T),
            mean_dist_freq = mean(dist_freq, na.rm = T),
            sd_dist_freq = sd(dist_freq, na.rm = T),
            n = n()) %>% 
  ungroup() %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>%
  st_as_sf(.)


## then the mean across the whole period
dat_26_all_plot <- dat_26_all %>%
  mutate(mean_dist_rate = ifelse(is.nan(mean_dist_rate), NA, mean_dist_rate), 
         mean_dist_freq = ifelse(is.nan(mean_dist_freq), NA, mean_dist_freq)) %>% 
  group_by(gridid) %>% 
  summarize(mean_dist_rate = mean(mean_dist_rate, na.rm = T),
            land_area = mean(land_area, na.rm = T),
            forest_area = mean(forest_area, na.rm = T)) %>% 
  mutate(forestShare = (1/land_area * forest_area)) %>% 
  mutate(mean_dist_rate = ifelse(mean_dist_rate > 0.5, 0.5, mean_dist_rate))



centroids <- st_centroid(dat_26_all_plot)
centroid_coords <- st_coordinates(centroids)

# Add the centroid coordinates to the original data frame
sf_with_centroids <- dat_26_all_plot %>%
  mutate(
    x = centroid_coords[, 1],  # X coordinate of centroid
    y = centroid_coords[, 2]   # Y coordinate of centroid
  )

dat_26_all_plot_mean <- sf_with_centroids


cols <- rev(met.brewer("OKeeffe2", 100, "continuous"))
DistMap <- ggplot() +
  geom_point(data=dat_26_all_plot_mean,
             aes(x=x, y=y, color=mean_dist_rate), size=0)+
  geom_point(data=dat_26_all_plot_mean,
             aes(x=x, y=y, color=mean_dist_rate, size=forestShare), show.legend = TRUE)+
  geom_spatvector(data = eu_shp, fill = "black", col = "black", alpha = 0.1) +
  geom_spatvector(data = europe, fill = "transparent", col = "black") +
  theme_classic() +
  scale_size_continuous(range = c(-0.5, 1.5)) + 
  scale_color_gradientn(colours = rev(cols), na.value = "lightgrey", limits = c(0, 0.5), breaks = c(0, 0.5)) +
  labs(color = bquote("disturbance \n rate [% year"^-1*"]"),        size = bquote("forest share")) + # Change color legend title here
  theme(axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, vjust = 3),
        legend.position = "right") +
  xlab("") + ylab("")

DistMap
ggsave(DistMap, filename = paste0(path_results, "/figures/FigS3.png"), width = 12, height = 9)
write_csv(dat_26_all_plot_mean, paste0(path_results, "/figures/plot_data/figS3_data_map.csv"))



## then the delta ------------
all_bas <- dat_26_all %>% filter(year == 2030) %>% st_drop_geometry() %>%
  mutate(mean_dist_rate = ifelse(mean_dist_rate == 0, 0.0001, mean_dist_rate))

all_fut <- dat_26_all %>% filter(year == 2100) %>% st_drop_geometry() %>%
  mutate(mean_dist_rate = ifelse(mean_dist_rate == 0, 0.0001, mean_dist_rate))

diff_df <- left_join(all_bas, all_fut, by = c("gridid")) %>%
  mutate(mean_diff = mean_dist_rate.y - mean_dist_rate.x,
         mean_diff_rel = mean_dist_rate.y / mean_dist_rate.x - 1,
         mean_diff_prct = (mean_dist_rate.y - mean_dist_rate.x) / mean_dist_rate.x * 100) %>%
  mutate(forestShare = (1/land_area.x * forest_area.x)) %>%
  mutate(across(c(mean_diff, mean_diff_rel, mean_diff_prct), ~ ifelse(is.finite(.), ., NA))) %>%
  left_join(., sf_obj_forest_mask, by = "gridid") %>%
  st_as_sf(.)


centroids <- st_centroid(diff_df)
centroid_coords <- st_coordinates(centroids)

# Add the centroid coordinates to the original data frame
sf_with_centroids <- diff_df %>%
  mutate(
    x = centroid_coords[, 1],  # X coordinate of centroid
    y = centroid_coords[, 2]   # Y coordinate of centroid
  )

dat_26_all_plot_diff <- sf_with_centroids

cols <- met.brewer("OKeeffe1", 100, "continuous")
DistMap <- ggplot() +
  geom_point(data=dat_26_all_plot_diff,
             aes(x=x, y=y, color=mean_diff), size=0)+
  geom_point(data=dat_26_all_plot_diff,
             aes(x=x, y=y, color=mean_diff, size=forestShare), show.legend = TRUE)+
  geom_spatvector(data = eu_shp, fill = "black", col = "black", alpha = 0.1) +
  geom_spatvector(data = europe, fill = "transparent", col = "black") +
  theme_classic() +
  scale_size_continuous(range = c(-0.5, 1.5)) +
  scale_color_gradientn(colours = rev(cols), na.value = "lightgrey", limits = c(-0.5, 0.5), breaks = c(-0.5, 0.5)) +
  labs(color = bquote("disturbance \n rate [% year"^-1*"]"),        size = bquote("forest share")) + # Change color legend title here
  theme(axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, vjust = 3),
        legend.position = "right") +
  xlab("") + ylab("")

DistMap
ggsave(DistMap, filename = paste0(path_results, "/figures/FigS6.png"), width = 12, height = 9)
write_csv(dat_26_all_plot_mean, paste0(path_results, "/figures/plot_data/figuS6_data_map.csv"))



### maps for the agents --------------------------------------------------------


# fire
cols <- brewer.pal(9, "YlOrBr")

# filter the data
dat_85 <- all_dat %>% filter(rcp == "rcp85") %>% filter(agent == "firegrid")

# prepare data and calculate dist rate per hexagon
dat_85_all <- dat_85 %>% dplyr::select(-sim) %>% 
  group_by(gridid, rcp, gcm, rep, year) %>% 
  summarise(dist_area = sum(dist_area, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>% 
  mutate(dist_rate = as.numeric(100 / forest_area * dist_area / 10),
         dist_freq = as.numeric(10 * forest_area / dist_area)) %>% 
  group_by(gridid, year) %>% 
  summarize(mean_dist_rate = mean(dist_rate, na.rm = T),
            sd_dist_rate = sd(dist_rate, na.rm = T),
            mean_dist_freq = mean(dist_freq, na.rm = T),
            sd_dist_freq = sd(dist_freq, na.rm = T),
            n = n()) %>% 
  ungroup() %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>%
  st_as_sf(.)


pct <- round(as.numeric(quantile(dat_85_all$mean_dist_rate, 0.9, na.rm = T)), 2)
print(pct)
max(dat_85_all$mean_dist_rate, na.rm = T)

## then the mean across the whole period
dat_85_all_plot <- dat_85_all %>%
  mutate(mean_dist_rate = ifelse(is.nan(mean_dist_rate), NA, mean_dist_rate), 
         mean_dist_freq = ifelse(is.nan(mean_dist_freq), NA, mean_dist_freq)) %>% 
  group_by(gridid) %>% 
  summarize(mean_dist_rate = mean(mean_dist_rate, na.rm = T),
            land_area = mean(land_area, na.rm = T),
            forest_area = mean(forest_area, na.rm = T)) %>% 
  mutate(forestShare = (1/land_area * forest_area)) %>% 
  mutate(mean_dist_rate = ifelse(forestShare < 0.1, 0, mean_dist_rate)) %>% 
  mutate(mean_dist_rate = ifelse(mean_dist_rate > pct, pct, mean_dist_rate))



centroids <- st_centroid(dat_85_all_plot)
centroid_coords <- st_coordinates(centroids)

# Add the centroid coordinates to the original data frame
sf_with_centroids <- dat_85_all_plot %>%
  mutate(
    x = centroid_coords[, 1],  # X coordinate of centroid
    y = centroid_coords[, 2]   # Y coordinate of centroid
  )

dat_85_all_plot_mean <- sf_with_centroids

DistMap <- ggplot() +
  geom_point(data=dat_85_all_plot_mean,
             aes(x=x, y=y, color=mean_dist_rate), size=0)+
  geom_point(data=dat_85_all_plot_mean,
             aes(x=x, y=y, color=mean_dist_rate, size=forestShare), show.legend = TRUE)+
  geom_spatvector(data = eu_shp, fill = "black", col = "black", alpha = 0.1) +
  geom_spatvector(data = europe, fill = "transparent", col = "black") +
  theme_classic() +
  scale_size_continuous(range = c(-0.5, 1.5)) + 
  scale_color_gradientn(colours = cols, na.value = "lightgrey", limits = c(0, pct), breaks = c(0, pct)) +
  labs(color = bquote("disturbance \n rate [% yr"^-1*"]"),        size = bquote("forest share")) + # Change color legend title here
  theme(axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, vjust = 3),
        legend.position = "right") +
  xlab("") + ylab("")

DistMap

print(pct)
ggsave(DistMap, filename = paste0(path_results, "/figures/FigS7.png"), width = 12, height = 9)
write_csv(dat_85_all_plot_mean, paste0(path_results, "/figures/plot_data/figS7_data_map.csv"))




# wind
cols <- brewer.pal(9, "BuPu")

# filter the data
dat_85 <- all_dat %>% filter(rcp == "rcp85") %>% filter(agent == "wind")

# prepare data and calculate dist rate per hexagon
dat_85_all <- dat_85 %>% dplyr::select(-sim) %>% 
  group_by(gridid, rcp, gcm, rep, year) %>% 
  summarise(dist_area = sum(dist_area, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>% 
  mutate(dist_rate = as.numeric(100 / forest_area * dist_area / 10),
         dist_freq = as.numeric(10 * forest_area / dist_area)) %>% 
  group_by(gridid, year) %>% 
  summarize(mean_dist_rate = mean(dist_rate, na.rm = T),
            sd_dist_rate = sd(dist_rate, na.rm = T),
            mean_dist_freq = mean(dist_freq, na.rm = T),
            sd_dist_freq = sd(dist_freq, na.rm = T),
            n = n()) %>% 
  ungroup() %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>%
  st_as_sf(.)

pct <- round(as.numeric(quantile(dat_85_all$mean_dist_rate, 0.9, na.rm = T)), 2)
print(pct)



## then the mean across the whole period
dat_85_all_plot <- dat_85_all %>%
  mutate(mean_dist_rate = ifelse(is.nan(mean_dist_rate), NA, mean_dist_rate), 
         mean_dist_freq = ifelse(is.nan(mean_dist_freq), NA, mean_dist_freq)) %>% 
  group_by(gridid) %>% 
  summarize(mean_dist_rate = mean(mean_dist_rate, na.rm = T),
            land_area = mean(land_area, na.rm = T),
            forest_area = mean(forest_area, na.rm = T)) %>% 
  mutate(forestShare = (1/land_area * forest_area)) %>% 
  mutate(mean_dist_rate = ifelse(forestShare < 0.1, 0, mean_dist_rate)) %>% 
  mutate(mean_dist_rate = ifelse(mean_dist_rate > pct, pct, mean_dist_rate))



centroids <- st_centroid(dat_85_all_plot)
centroid_coords <- st_coordinates(centroids)

# Add the centroid coordinates to the original data frame
sf_with_centroids <- dat_85_all_plot %>%
  mutate(
    x = centroid_coords[, 1],  # X coordinate of centroid
    y = centroid_coords[, 2]   # Y coordinate of centroid
  )

dat_85_all_plot_mean <- sf_with_centroids

DistMap <- ggplot() +
  geom_point(data=dat_85_all_plot_mean,
             aes(x=x, y=y, color=mean_dist_rate), size=0)+
  geom_point(data=dat_85_all_plot_mean,
             aes(x=x, y=y, color=mean_dist_rate, size=forestShare), show.legend = TRUE)+
  geom_spatvector(data = eu_shp, fill = "black", col = "black", alpha = 0.1) +
  geom_spatvector(data = europe, fill = "transparent", col = "black") +
  theme_classic() +
  scale_size_continuous(range = c(-0.5, 1.5)) + 
  scale_color_gradientn(colours = cols, na.value = "lightgrey", limits = c(0, pct), breaks = c(0, pct)) +
  labs(color = bquote("disturbance \n rate [% yr"^-1*"]"),        size = bquote("forest share")) + # Change color legend title here
  theme(axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, vjust = 3),
        legend.position = "right") +
  xlab("") + ylab("")

DistMap

ggsave(DistMap, filename = paste0(path_results, "/figures/FigS8.png"), width = 12, height = 9)
write_csv(dat_85_all_plot_mean, paste0(path_results, "/figures/plot_data/figS8_data_map.csv"))




# bark beetle
cols <- brewer.pal(9, "YlGn")

# filter the data
dat_85 <- all_dat %>% filter(rcp == "rcp85") %>% filter(agent == "bbgrid")

# prepare data and calculate dist rate per hexagon
dat_85_all <- dat_85 %>% dplyr::select(-sim) %>% 
  group_by(gridid, rcp, gcm, rep, year) %>% 
  summarise(dist_area = sum(dist_area, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>% 
  mutate(dist_rate = as.numeric(100 / forest_area * dist_area / 10),
         dist_freq = as.numeric(10 * forest_area / dist_area)) %>% 
  group_by(gridid, year) %>% 
  summarize(mean_dist_rate = mean(dist_rate, na.rm = T),
            sd_dist_rate = sd(dist_rate, na.rm = T),
            mean_dist_freq = mean(dist_freq, na.rm = T),
            sd_dist_freq = sd(dist_freq, na.rm = T),
            n = n()) %>% 
  ungroup() %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>%
  st_as_sf(.)

pct <- round(as.numeric(quantile(dat_85_all$mean_dist_rate, 0.97, na.rm = T)), 2)
print(pct)


## then the mean across the whole period
dat_85_all_plot <- dat_85_all %>%
  mutate(mean_dist_rate = ifelse(is.nan(mean_dist_rate), NA, mean_dist_rate), 
         mean_dist_freq = ifelse(is.nan(mean_dist_freq), NA, mean_dist_freq)) %>% 
  group_by(gridid) %>% 
  summarize(mean_dist_rate = mean(mean_dist_rate, na.rm = T),
            land_area = mean(land_area, na.rm = T),
            forest_area = mean(forest_area, na.rm = T)) %>% 
  mutate(forestShare = (1/land_area * forest_area)) %>% 
  mutate(mean_dist_rate = ifelse(forestShare < 0.1, 0, mean_dist_rate)) %>% 
  mutate(mean_dist_rate = ifelse(mean_dist_rate > pct, pct, mean_dist_rate))



centroids <- st_centroid(dat_85_all_plot)
centroid_coords <- st_coordinates(centroids)

# Add the centroid coordinates to the original data frame
sf_with_centroids <- dat_85_all_plot %>%
  mutate(
    x = centroid_coords[, 1],  # X coordinate of centroid
    y = centroid_coords[, 2]   # Y coordinate of centroid
  )

dat_85_all_plot_mean <- sf_with_centroids

DistMap <- ggplot() +
  geom_point(data=dat_85_all_plot_mean,
             aes(x=x, y=y, color=mean_dist_rate), size=0)+
  geom_point(data=dat_85_all_plot_mean,
             aes(x=x, y=y, color=mean_dist_rate, size=forestShare), show.legend = TRUE)+
  geom_spatvector(data = eu_shp, fill = "black", col = "black", alpha = 0.1) +
  geom_spatvector(data = europe, fill = "transparent", col = "black") +
  theme_classic() +
  scale_size_continuous(range = c(-0.5, 1.5)) + 
  scale_color_gradientn(colours = cols, na.value = "lightgrey", limits = c(0, pct), breaks = c(0, pct)) +
  labs(color = bquote("disturbance \n rate [% yr"^-1*"]"),        size = bquote("forest share")) + # Change color legend title here
  theme(axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, vjust = 3),
        legend.position = "right") +
  xlab("") + ylab("")

DistMap

ggsave(DistMap, filename = paste0(path_results, "/figures/FigS9.png"), width = 12, height = 9)
write_csv(dat_85_all_plot_mean, paste0(path_results, "/figures/plot_data/figS9_data_map.csv"))


### end ------------------------------------------------------------------------
