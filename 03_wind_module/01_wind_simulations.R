#-------------------------------------------------------------------------------
# Script Name: Storm Occurrence Modeling and Simulation
# Description: This script models storm occurrences based on historical wind data 
#              and simulates storm events for future scenarios. It performs data 
#              preprocessing, builds predictive models for storm probabilities, 
#              and generates simulations for storm number, area, and severity.
#              Outputs include storm occurrence predictions, storm size distributions,
#              and simulated storm event datasets.
#
# Author: Marc Grünig & Cornelius Senf
#
# Input Data:
# - CSV files: Patch data, storm occurrence data, storm return interval rasters.
# - Raster files: Storm data, including wind speed return intervals, mean, and max windspeed.
#
# Output:
# - CSV files: predicted storm events with location, storm size and severity.
#
# Note:
# - In the folder "wind_module/wind_event_series/" we provide one precalculated 
#   example.
# - underlying disturbance data can be obtained here: https://zenodo.org/records/4607230
#
#-------------------------------------------------------------------------------



### libraries ------------------------------------------------------------------

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(mgcViz)


path <- "/.../"


#### Data Preprocessing ####

### Load patch data and combine with storm data
# this data is obtained from the European disturbance maps (see https://zenodo.org/records/4607230 for download) 

# predictions <- "/DisturbanceMappingEurope/attribution/results/predictions/"
# cntrs <- list.files(predictions, ".csv$") %>% 
#   str_sub(., 12, nchar(.) - 4)
# list.files(predictions, ".csv$", full.names = TRUE) %>%
#   map(read_csv) %>%
#   map(., ~ filter(., prediction_breakage_prob > 0.5)) %>%
#   map(., ~ dplyr::select(., patch, year, x, y, area)) %>%
#   map2(.x = ., .y = cntrs, ~ write_csv(.x, paste0(path, "/wind_module/patches/wind_patches_", .y, ".csv")))

d <- list.files(
  paste0(path, "/wind_module/patches/"), 
  pattern = "wind_patches", 
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


# this data is from copernicus storm database
stormdata <- rast(
  c(
    paste0(path, "/wind_module/copernicus_storms/windspeed_return_intervals_epsg3035_europe_10km.tif"),
    paste0(path, "/wind_module/copernicus_storms/average_windspeed_epsg3035_europe_10km.tif"),
    paste0(path, "/wind_module/copernicus_storms/maximum_windspeed_epsg3035_europe_10km.tif")
  )
)

years <- 1987:2016

stormdata_df <- stormdata %>%
  as.data.frame(
    xy = TRUE,
    cells = TRUE
  ) %>%
  rename(
    id = cell,
    return_5 = lyr.1,
    return_10 = lyr.2,
    return_50 = lyr.3,
    return_100 = lyr.4
  )

d_sf <- d %>% 
  filter(year %in% years) %>%
  st_as_sf(
    .,
    coords = c("x", "y"),
    crs = st_crs(stormdata)
  )

d_stormdata <- d_sf %>% 
  terra::extract(
    x = stormdata,
    y = ., 
    cells = TRUE,
    ID = FALSE
  ) %>%
  mutate(
    year = d_sf$year,
    area = d_sf$area
  )

names(d_stormdata) <- c("return_5", "return_10", "return_50", "return_100", "mean", "max", "id", "year", "area")
d_stormdata <- as_tibble(
  d_stormdata
)

stormdata_disturbances <- d_stormdata %>%
  group_by(id, return_5, return_10, return_50, return_100, mean, max) %>%
  summarize(n_total = n(),
            n_year = length(unique(year)),
            area_sum = sum(area),
            area_prop = area_sum / 10000^2)

stormdata_nodisturbance <- stormdata_df %>%
  filter(!(id %in% stormdata_disturbances$id)) %>%
  mutate(
    n_total = 0,
    n_year = 0,
    area_sum = 0,
    area_prop = 0
  ) %>%
  select(-x, -y)

stormdata_dat <- list(stormdata_disturbances, stormdata_nodisturbance) %>%
  bind_rows() %>%
  ungroup() %>%
  filter(
    if_any(return_5:return_100, ~!is.na(.))
  ) %>%
  mutate(
    binary = ifelse(n_total > 0, 1, 0)
  )

stormdata_dat <- stormdata_dat %>%
  left_join(stormdata_df)

ggplot(
  data = stormdata_dat
) +
  geom_point(
    aes(
      x = return_5,
      y = max
    )
  )




#### Storm occurence Model ####

fit_prop_random <- mgcv::gam(binary ~ 1,
                             data = stormdata_dat,
                             family = "binomial")

fit_prop_wind_nosa_5 <- mgcv::gam(binary ~ s((return_5)),
                                  data = stormdata_dat,
                                  family = "binomial")

fit_prop_wind_nosa_5_log <- mgcv::gam(binary ~ s(log10(return_5)),
                                      data = stormdata_dat,
                                      family = "binomial")

fit_prop_wind_nosa_10 <- mgcv::gam(binary ~ s((return_10)),
                                   data = stormdata_dat,
                                   family = "binomial")

fit_prop_wind_nosa_10_log <- mgcv::gam(binary ~ s(log10(return_10)),
                                       data = stormdata_dat,
                                       family = "binomial")

fit_prop_wind_nosa_50 <- mgcv::gam(binary ~ s((return_50)),
                                   data = stormdata_dat,
                                   family = "binomial")

fit_prop_wind_nosa_50_log <- mgcv::gam(binary ~ s(log10(return_50)),
                                       data = stormdata_dat,
                                       family = "binomial")


fit_prop_wind_nosa_5_log_max <- mgcv::gam(binary ~ s(log10(return_5)) + s(max),
                                          data = stormdata_dat,
                                          family = "binomial")

AIC(
  fit_prop_random, 
  fit_prop_wind_nosa_5,
  fit_prop_wind_nosa_10,
  fit_prop_wind_nosa_50,
  fit_prop_wind_nosa_5_log,
  fit_prop_wind_nosa_10_log,
  fit_prop_wind_nosa_50_log,
  fit_prop_wind_nosa_5_log_max
)


dev.off()

plot(fit_prop_wind_nosa_5_log_max)

summary(fit_prop_wind_nosa_5_log_max)

simulation <- DHARMa::simulateResiduals(
  fittedModel = fit_prop_wind_nosa_5_log_max, 
  plot = FALSE
)

plot(simulation)

stormdata_predict <- stormdata_df

stormdata_predict$pred <- predict(fit_prop_wind_nosa_5_log_max, newdata = stormdata_predict)

stormdata_predict$prob <- boot::inv.logit(stormdata_predict$pred) / length(years)

stormdata_map <- stormdata_predict %>%
  st_as_sf(
    x = .,
    coords = c("x", "y"),
    crs = st_crs(stormdata)
  )

ggplot() +
  geom_sf(data = stormdata_map, aes(col = prob)) +
  scale_color_viridis_c() +
  labs(fill = "Annual probability")

stormdata_map <- tidyterra::as_spatraster(stormdata_predict, xycols = 2:3)

plot(stormdata_map, 9)

writeRaster(
  stormdata_map,
  paste0(path, "/wind_module/results/storm_probability_prediction.tif"),
  overwrite = TRUE
)



#### Group disturbances to model storm area and number of storms ####

gridid <- subset(stormdata, 1)
values(gridid) <- 1:ncell(gridid)
names(gridid) <- "gridid"
gridid_vect <- as.polygons(gridid)

k <- 0
stormdata_annual <- vector("list", length(years))

for (y in years) {
  
  print(y)
  
  k <- k + 1
  
  gridid_sel <- terra::extract(gridid, d_sf %>% filter(year == y)) %>%
    .$gridid %>%
    unique(.)
  
  stormdata_annual[[k]] <- gridid_vect %>%
    as_sf() %>%
    filter(gridid %in% gridid_sel) %>%
    summarize() %>%
    st_cast(., "POLYGON") %>%
    mutate(id = 1:n(),
           year = y)
  
}


stormdata_annual_summary <- stormdata_annual %>%
  bind_rows() %>%
  mutate(
    area = as.double(st_area(.)),
    ncells = area / 10000^2
  ) %>%
  st_drop_geometry()

write_csv(
  stormdata_annual_summary,
  paste0(path, "/wind_module/results/stormsizes.csv")
)

ggplot(stormdata_annual_summary,
       aes(x = ncells)) +
  geom_density(adjust = 6, linewidth = 1) +
  labs(x = "Storm footprint area (# of 100 sqkm cells)",
       y = "Density") +
  scale_x_log10()



### Storm severity

stormdata_dat_only_storms <- stormdata_dat %>%
  filter(binary == 1)

write_csv(
  stormdata_dat_only_storms,
  paste0(path, "/wind_module/results/stormseverity.csv")
)



ggplot(
  data = stormdata_dat_only_storms
) +
  geom_density(
    aes(
      x = area_prop * 100
    ),
    adjust = 3,
    linewidth = 1
  ) +
  labs(x = "Storm severity (% of footprint disturbed)",
       y = "Density") +
  scale_x_log10()


### Number of storms

number_of_storms <- stormdata_annual_summary %>%
  group_by(year) %>%
  summarise(n = n()) %>%
  ungroup()

ggplot(data = number_of_storms, 
       aes(x = n)) +
  geom_histogram() +
  labs(x = "Number of storms",
       y = "Density")





#### Create simulations #### ------------------------------------------------------------------------------------------
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)

# load probability raster
stormdata_map <- rast(paste0(path, "wind_module/results/storm_probability_prediction.tif"))
prob_raster <- subset(stormdata_map, 9)

# stormdata tables created above
stormdata_annual_summary <- read_csv(paste0(path, "wind_module/results/stormsizes.csv"))
number_of_storms <- stormdata_annual_summary %>%
  group_by(year) %>%
  summarise(n = n()) %>%
  ungroup()

stormdata_dat_only_storms <- read_csv(paste0(path, "wind_module/results/stormseverity.csv"))


n_sim_y <- 100
n_sim_d <- 90

sim_years <- 1:n_sim_y

for (d in 1:n_sim_d){
  
  print(d)
  
  storm_samples <- vector("list", n_sim_y)
  
  for (y in 1:n_sim_y) {
    
    number_storms_tmp <- sample(
      number_of_storms$n, 
      size = 1
    )
    
    year_tmp <- rep(sim_years[y], number_storms_tmp)
    
    size_storms_tmp <- sample(
      stormdata_annual_summary$ncells,
      size = number_storms_tmp
    )
    
    severity_storm_tmp <- sample(
      stormdata_dat_only_storms$area_prop,
      size = number_storms_tmp
    )
    
    storm_samples[[y]] <- spatSample(
      x = prob_raster, 
      size = number_storms_tmp, 
      method = "weight",
      xy = TRUE,
      cell = TRUE
    ) %>%
      select(
        -prob
      ) %>%
      mutate(
        number_of_cells = size_storms_tmp,
        proportion_of_cell = severity_storm_tmp,
        year = year_tmp
      ) %>%
      rename(
        cell_position = cell
      )
    
  }
  
  storm_samples <- storm_samples %>%
    bind_rows()
  
  write_csv(
    storm_samples,
    paste0(path, "/wind_module/wind_event_series/simulated_storms_draw", d, ".csv"
    )
  )
} 



