# --------------------------------------------------------------
# Bark Beetle Module Scipt 01: Background outbreak probability
# --------------------------------------------------------------
# Author: Marc Grünig
# Date: 2025-01-07
# 
# Description:
# This script calculates a background bark beetle outbreak probability map
# which is used in SVD. Additionally, this script calibrates a model that 
# defines how strong the background outbreak probability is modified by VPD.
# 
# 
# Notes:
# - The script relies on data from the European disturbance maps (https://zenodo.org/records/4570157)
# 
# --------------------------------------------------------------






### libraries --------------------------------------------
library(terra)
library(sf)
library(tidyverse)
library(tidyterra)
library(patchwork)


path <- "/.../"


### load data --------------------------------------------
ref_grid_100km <- rast(paste0(path, "/07_reference_grids/reference_grid_100km.tif"))
ref_grid_10km <- rast(paste0(path, "/07_reference_grids/reference_grid_10km.tif"))
ref_grid_1km <- rast(paste0(path, "/07_reference_grids/reference_grid_1km.tif"))


gridid <- ref_grid_100km
gridid_vect <- as.polygons(ref_grid_100km)
gridid_df <- as.data.frame(gridid, xy = T)


# load wind/bb patches - this comes from European disturbance maps
# note that the disturbance maps have one category for wind and bark beetle
# we remove thus the extreme wind years where we assume that wind contriuted over proportionally (see below)
d <- list.files(
  paste0(path, "/bbtl_module/patches"), 
  pattern = "patches", 
  full.names = TRUE
) %>%
  map(read_csv) %>%
  map(., ~ mutate(.,
                  patch = as.integer(patch),
                  year = as.integer(year),
                  x = as.double(x),
                  y = as.double(y),
                  area = as.double(area))) %>%
  bind_rows()

# get years and get rid of extreme wind years
years <- unique(d$year)
years <- years %>% sort()
years <- years[!years %in% c(1990, 2000, 2007, 2016)]

# make spatial
d_sf <- d %>% 
  filter(year %in% years) %>%
  # mutate(x = unlist(map(d$geometry,1)),
  #        y = unlist(map(d$geometry,2))) %>%
  st_as_sf(
    .,
    coords = c("x", "y"),
    crs = st_crs(ref_grid_100km)
  )


# get dist per gridcells 100km
k <- 0
dat_annual <- vector("list", length(years))

for (y in years) {
  
  print(y)
  
  k <- k + 1
  
  gridid_sel <- terra::extract(gridid, d_sf %>% filter(year == y)) %>%
    .$gridid
  # biomeid <- terra::extract(ecoreg, d_sf %>% filter(year == y)) %>%
  #   .$BIOME
  gridid_sel <- as.data.frame(cbind(d_sf %>% filter(year == y), gridid = gridid_sel)) %>% 
    dplyr::select(patch, year, area, gridid)
  
  dat_annual[[k]] <- gridid_sel
  
}

dat_annual_summary <- dat_annual %>%
  bind_rows()

df <- dat_annual_summary %>% mutate(area_ha = area * 0.0001)
mean(df$area)
summary(df)

par(mfrow = c(2,2))



# check the mean disturbance activity year for each pixel - This results in the disturbance rate per 100km grid cell
mean_bb <- dat_annual_summary %>% 
  group_by(gridid, year) %>%
  summarize(area_tot = sum(area)) %>% 
  ungroup()  %>%
  group_by(gridid) %>%
  #slice(which.max(area_tot)) %>%
  summarize(area_tot = mean(area_tot)) %>% 
  ungroup()

mean_bb <- mean_bb %>% mutate(area_tot = area_tot / 1000000, # to km2
                              area_prop = area_tot / 10000, # get the proportion within 100x100km grid
                              #area_prop = area_prop * 0.3, # forest area
                              area_prop = area_prop * 0.5) # deduct the wind

summary(mean_bb)



# create raster from it
grid_new <- left_join(gridid_df, mean_bb, by = c("gridid")) 
new_rast <- rast(grid_new[, c("x", "y", "area_prop")])
new_rast[is.na(new_rast)] <- 0
crs(new_rast) <- crs(ref_grid_100km)
new_rast <- terra::project(new_rast, ref_grid_100km)
plot(new_rast, main = "mean disturbance")
summary(grid_new)


### get forest mask - Now this needs to be corrected with the forest mask 
eu_shp <- vect(paste0(path, "/reference_grids/eu_mask.shp"))
proj_laea <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

# load forest mask
forest_mask <- rast(paste0(path, "/gis/forest_mask_new_laea_masked.tif"))
forest_mask <- terra::crop(forest_mask, eu_shp)
forest_mask <- terra::mask(forest_mask, eu_shp)
forest_mask[forest_mask > 0] <- 1
forest_mask[is.na(forest_mask)] <- 0

# aggregate to 100km grid
forest_share <- terra::project(forest_mask, ref_grid_100km, method = "bilinear")
plot(forest_share)

# then multiply with the disturbance layer
new_rast2 <- new_rast * forest_share
plot(new_rast2)


### Then within the forest, the spruce share differs
# spruce shares were calculated from the initial forest layer
spruce_share <- rast(paste0(path, "/bbtl_module/spruce_share.tif"))
plot(spruce_share)
spruce_share[spruce_share < 0.1] <- 0
reduced_proba <- new_rast2 / spruce_share
plot(reduced_proba)

reduced_proba[is.infinite(reduced_proba)] <- 0
plot(reduced_proba)

### look at the statistics
df_grid <- as.data.frame(reduced_proba)
summary(df_grid)


### Finally we can mask it with the spruce distribution range
# get vegetation from Brus maps
dom_veg <- list.files(paste0(path, "/gis/EU_TreeMap_Brus_etal"), 
                      pattern = "_wgs.tif", full.names = TRUE)
dom_veg <- dom_veg[grepl("Picea", dom_veg)]
shp <- rast(dom_veg)
shp <- terra::project(shp, ref_grid_100km)
shp <- terra::focal(shp, na.policy = "only", w = c(3,3), fun = mean)
shp <- terra::mask(shp, ref_grid_100km)

shp[shp < 1] <- 0
shp[shp > 0] <- 0.00001
shp[shp == 0] <- NA
shp <- terra::project(shp, new_rast)
plot(shp)

reduced_proba <- terra::mask(reduced_proba, shp)
df_grid <- as.data.frame(reduced_proba)
summary(df_grid)



# get vegetation from Brus maps
shp <- rast(dom_veg)
shp <- terra::project(shp, ref_grid_100km)
shp <- terra::focal(shp, na.policy = "only", w = c(3,3), fun = mean)
shp <- terra::mask(shp, ref_grid_100km)

# assign base proba
shp[shp < 1] <- 0
shp[shp > 0] <- 0.00001
plot(shp)

shp <- terra::project(shp, new_rast)
plot(shp)

# save
terra::writeRaster(shp, paste0(path, "/bbtl_module/background_equalproba_100km.tif"), overwrite = T, datatype = "FLT4S")



### now the increase in probability with vpd ------------------------------------
# first, we extract the VPD conditions for each barkbeetle patch from 1986:2020
# the results are stored in the folder
# 
# d_clim <- list()
# 
# for(y in 1:length(years)){
#   
#   focal_year <- as.numeric(years[y])
#   d_sub <- d %>% filter(year == focal_year)
#   
#   x <- rast(paste0(path, "/Data/climatedata/euro-cordex/yearly/summer_vpd/vpd_summer_", focal_year, "_ensemble_historical_leae_v2_biascorrected.tif"))
#   vpd_df <- as.data.frame(x, xy = T)
#   vpd_df$cell <- as.numeric(rownames(vpd_df))
#   
#   new_df <- cbind(d_sub, terra::extract(x, d_sub[, c("x", "y")]))
#   d_sf <- new_df %>% 
#     st_as_sf(
#       .,
#       coords = c("x", "y"),
#       crs = st_crs(x)
#     )
#   
#   x_extract <- d_sf %>% 
#     terra::extract(
#       x = x,
#       y = ., 
#       cells = TRUE,
#       ID = FALSE
#     ) %>%
#     mutate(
#       year = d_sf$year,
#       area = d_sf$area
#     )
#   
#   disturbances <- x_extract %>%
#     group_by(cell) %>%
#     summarize(vpd = mean(mean, na.rm = T),
#               n_total = n(),
#               n_year = length(unique(year)),
#               area_sum = sum(area),
#               area_prop = area_sum / 10000^2)
#   
#   nodisturbance <- vpd_df[, c("x", "y", "cell", "mean")] %>%
#     rename("vpd" = "mean") %>% 
#     filter(!(cell %in% disturbances$cell)) %>%
#     mutate(
#       cell = as.numeric(cell),
#       n_total = 0,
#       n_year = 0,
#       area_sum = 0,
#       area_prop = 0
#     ) %>%
#     select(-x, -y)
#   
#   dat <- list(disturbances, nodisturbance) %>%
#     bind_rows() %>%
#     ungroup() %>%
#     mutate(
#       binary = ifelse(n_total > 0, 1, 0)
#     ) %>% 
#     drop_na()
#   
#   dat <- dat %>%
#     left_join(vpd_df, by = c("cell")) %>% 
#     mutate(year = focal_year)
#   
#   d_clim[[y]] <- dat
#   
# }
# 
# d_clim <- do.call(rbind, d_clim)

# long_term <- d_clim %>%
#   dplyr::select(cell, vpd, n_total) %>% 
#   group_by(cell) %>% 
#   summarize(mean_vpd = mean(vpd, na.rm = T),
#             mean_disturb = mean(n_total)) %>% 
#   left_join(., vpd_df, by = c("cell")) %>% 
#   dplyr::select(x, y, mean_disturb, mean_vpd)
# 
# write_csv(long_term, paste0(path, "bbtl_module/bb_patches_vpd.csv"))

long_term <- read_csv(paste0(path, "bbtl_module/bb_patches_vpd.csv"))

# convert to raster
disturb_rast <- rast(long_term[, c("x", "y", "mean_disturb")])
crs(disturb_rast) <- proj_laea
vpd_rast <- rast(long_term[, c("x", "y", "mean_vpd")])
plot(vpd_rast)
cut_vpd <- vpd_rast

# remove everything with really low vpd assuming that it was winter storm
cut_vpd[cut_vpd < 0.5] <- NA
plot(cut_vpd)

# use dom veg to remove 
shp <- rast(dom_veg)
shp <- terra::project(shp, proj_laea)
shp[shp < 5] <- NA
disturb_rast <- rast(long_term[, c("x", "y", "mean_disturb")])
crs(disturb_rast) <- proj_laea
shp <- terra::project(shp, disturb_rast)
disturb_rast <- terra::mask(disturb_rast, shp)
plot(disturb_rast)

long_term_masked <- na.omit(cbind(long_term[, c("x", "y")], 
                                  mean_disturb = terra::extract(disturb_rast, long_term[, c("x", "y")])[,2], 
                                  mean_vpd = terra::extract(vpd_rast, long_term[, c("x", "y")])[,2] ))

# here we select only the years 2017:2020
years <- unique(d$year)
new_list <- list()
years <- years %>% sort()
years <- c(2017, 2018, 2019, 2020)

# then we obtain the patches and add the vpd to create a presence - absences dataframe 
for(y in 1:length(years)){
  
  focal_year <- as.numeric(years[y])
  df_year <- d_clim %>% filter(year == focal_year)
  df_year <- df_year %>%
    left_join(., long_term_masked, by = c("x", "y")) %>% 
    drop_na() %>% 
    mutate(diff = n_total / mean_disturb) %>% 
    mutate(diff_vpd = vpd - mean_vpd) 

  pres <- df_year %>%
    #filter(diff >= 1) %>% 
    filter(binary == 1) %>%
    dplyr::select(vpd, diff_vpd) %>% 
    mutate(P = 1)
  
  abs <- df_year %>% 
    filter(binary == 0) %>% 
    #filter(diff < 1) %>% 
    dplyr::select(vpd, diff_vpd) %>% 
    #sample_n(nrow(pres)) %>%  # get same number as presences
    mutate(P = 0)
  
  df_occ <- rbind(pres, abs)
  
  
  new_list[[y]] <- df_occ
  
}
new_list_df <- do.call(rbind, new_list)


# finally we fit a model
p1 <- ggplot(data = new_list_df,  aes(x = vpd, y = P)) +
  geom_point() +
  geom_smooth(method = "glm", , se = T, method.args = list(family = "binomial")) +
  ggtitle("Disturbance probability ~ summer VPD") 
p1

fit_vpd <- glm(P ~ vpd,
               data = new_list_df_filt,
               family = "binomial")


### visualize how the function modifies the background probability
# test set of vpds
new_data <- data.frame(vpd = seq(0, 2, 0.1))
predicted_probabilities <- predict(fit_vpd, newdata = as.data.frame(new_data), type = "response")

# get coefs from model
new_vpd <- 0.0
intercept <- as.numeric(coef(fit_vpd)[1])
coef_vpd <- as.numeric(coef(fit_vpd)[2])

# Adjust the intercept for vpd = 0
adjusted_intercept <- intercept - coef_vpd * mean(new_list_df_filt$vpd)

# Calculate the log-odds
log_odds <- adjusted_intercept + coef_vpd * new_data

# Apply the logistic function to get the probability
intercept_param <- 1
background_proba <- 0.0005

predicted_probabilities <- as.vector((intercept_param + 1 / (1 + exp(-log_odds))) * background_proba)
print(predicted_probabilities)

new_df <- cbind(new_data, probs = predicted_probabilities)
colnames(new_df) <- c("vpd", "probs")

ggplot(new_df, aes(x = vpd, y = probs)) +
  geom_line()

