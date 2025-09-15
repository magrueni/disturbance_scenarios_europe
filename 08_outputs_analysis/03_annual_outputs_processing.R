# --------------------------------------------------------------------------------------
# Script Name: annual_outputs_processing.R
# Description: This script processes and visualizes the annual outputs from the
#              disturbance simulation for fire, wind, bark beetle outbreaks, 
#              and forest management.
#              The data is extracted from multiple simulation runs and historical datasets, 
#              cleaned, aggregated. 
#
# Author: Marc Grünig
# Date: 7.2.2025
#
#
# Input:
#   - Historical data from fire, wind, and bark beetle disturbance events.
#   - Climate simulation data (fire, wind, bark beetle) from multiple RCP scenarios.
#   - Simulation paths and CSV/GPKG files containing model outputs.
#
# Output:
#   - CSV files summarizing disturbance impacts (fire, wind, bark beetle) over time.
#
# Notes:
# - Data from the disturbance mapping cannot be provided here but downloaded here:
# https://zenodo.org/records/4570157
# - Data from Patacca et al., 2023 can be downloaded here: https://zenodo.org/records/7319179
#
# ------------------------------------------------------------------------------



### load libraries -------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(dplyr)
library(MetBrewer)
library(sf)
library(stars)
library(terra)


# define paths
path <- "/.../"
path_public <- "/.../"
path_results <- paste0(path, "/10_results/")
path_sims <- paste0(path, "/09_svd_simulations/annual_outputs/")

### load historical data -------------------------------------------------------
cntrs <- list.files(paste0(path_public, "/Projects/DisturbanceMappingEurope/attribution/results/complexes"), ".gpkg$") %>% 
  str_sub(., 14, nchar(.) - 5)
cntrs <- cntrs[grepl("*eps150", cntrs)] %>% str_sub(., 0, nchar(.) - 7)

k <- 0
out <- vector("list", length(cntrs))

cntrs.fullnam <- list.files(paste0(path, "/firesize/data/complexes"), ".gpkg", full.names = TRUE)
cntrs.fullnam <- cntrs.fullnam[grepl("*eps150", cntrs.fullnam)]
cntrs.fullnam <- cntrs.fullnam[!grepl("ukraine|belarus|moldova", cntrs.fullnam)]

for (i in cntrs.fullnam) {
  
  k <- k + 1
  out[[k]] <- read_sf(i) 
  
}

a <- out %>%
  bind_rows(.id = "clust")


# read the grid
grid <- read_sf(paste0(path, "/firemodelling/clim_data/climategrid_epsg3035_090122.gpkg"))
grid_raster <- st_rasterize(grid)
grid_r <- rast(grid_raster)

d_sf <- a %>% 
  st_as_sf(., coords = c("x", "y"), crs = st_crs(grid))

fire_vect <- vect(d_sf)


a_fire <- as.data.frame(fire_vect) %>% 
  dplyr::select(year, complex_area_burnt_m2) %>% 
  group_by(year) %>% 
  summarise(n_fires = n(),
            area_tot  = sum(complex_area_burnt_m2*0.0001)) %>% 
  mutate(dat = "real")

# load historical data for wind and barkbeetle
cntrs <- list.files(paste0(path_public, "/Projects/DisturbanceMappingEurope/attribution/results/update/predictions/"), ".csv") %>% 
  str_sub(., 12, nchar(.) - 4)
cntrs <- cntrs[!grepl("ukraine|belarus|moldova", cntrs)]

k <- 0
out <- vector("list", length(cntrs))
cntrs.fullnam <- list.files(paste0(path_public, "/Projects/DisturbanceMappingEurope/attribution/results/update/predictions/"), ".csv", full.names = TRUE)
cntrs.fullnam <- cntrs.fullnam[!grepl("ukraine|belarus|moldova", cntrs.fullnam)]
cntrs.fullnam <- cntrs.fullnam[grepl("biotic", cntrs.fullnam)]

for (i in cntrs.fullnam) {
  
  country_name <- i %>% 
    str_sub(., 117, nchar(.) - 4)
  k <- k + 1
  x <- read_csv(i) %>% 
    mutate(country = country_name)
  out[[k]] <- x
  
}

a <- out %>%
  bind_rows(.id = "clust")


a_wind <- a %>%
  filter(prediction_class == "barkbeetle_windthrow") %>%
  dplyr::select(country, patch, year, area_m2 = area) %>%
  group_by(year) %>% 
  summarise(n_storms = n(),
            area_tot = sum(area_m2) * 0.0001) %>% 
  mutate(dat = "real")

a_bb <- a %>%
  filter(prediction_class == "barkbeetle_windthrow") %>%
  dplyr::select(country, patch, year, x, y, area_m2 = area) %>%
  group_by(year) %>% 
  summarise(n_storms = n(),
            area_tot = sum(area_m2) * 0.0001) %>% 
  mutate(dat = "real") 

# load data from Patacca et al., 2023 and calculate percentages to split up the remote sensing data
dat <- read.csv(paste0(path, "/disturbance_timeseries_database.csv"), sep = ";")


dat_filt <- dat %>% 
  filter(Year %in% c(1986:2020)) %>% 
  filter(disturbance_driver == "Bark Beetles") %>% 
  dplyr::select(year = Year, Reported_m3) %>% 
  group_by(year) %>% 
  summarize(Reported_m3 = sum(Reported_m3, na.rm = T))

dat_filt <- rbind(dat_filt, c(2020, mean(dat_filt$Reported_m3, na.rm = T))) # for 2020 take average ratio


dat_summed <- dat %>% 
  filter(Year %in% c(1986:2020)) %>% 
  filter(disturbance_driver %in% c("Bark Beetles", "Wind")) %>% 
  dplyr::select(year = Year, Reported_m3) %>% 
  group_by(year) %>% 
  summarize(tot = sum(Reported_m3, na.rm = T))

dat_summed <- rbind(dat_summed, c(2020, mean(dat_summed$tot, na.rm = T))) # for 2020 take average ratio

dat_filt_new <- dat_filt %>% 
  left_join(., dat_summed, by = c("year")) %>% 
  mutate(pct = 100/tot * Reported_m3)

a_bb_new <- a_bb %>% 
  left_join(dat_filt_new, by = "year") %>% 
  mutate(area_tot = area_tot * pct/100)

a_wind_new <- a_wind %>% 
  left_join(dat_filt_new, by = "year") %>% 
  mutate(area_tot = area_tot * (100-pct)/100)

a_mgmt <- a %>%
  filter(prediction_class == "harvest") %>%
  dplyr::select(country, patch, year, x, y, area_m2 = area) %>%
  group_by(year) %>% 
  summarise(n_storms = n(),
            area_tot = sum(area_m2) * 0.0001) %>% 
  mutate(dat = "real")




### load simulation data -------------------------------------------------------


# fire --
fire_dat <- list()
for(i in c(1:90)){
  
  dat <- read_csv(paste0(path_sims, "/output_sim_eu_", i, "/fire_run_sim_eu_", i, ".csv"))
  dat <- dat %>% dplyr::select(year, realized_size) %>% 
    mutate(Year = as.numeric(year) + 2019) %>% 
    mutate(sim = paste(i),
           scen = ifelse(i %in% seq(1, 90, 3), "rcp85", 
                         ifelse(i %in% seq(2, 90, 3), "rcp45", "rcp26"))) %>% 
    mutate(scen = ifelse(scen == "rcp85", "RCP8.5",
                         ifelse(scen == "rcp45", "RCP4.5", "RCP2.6")))
  
  fire_dat[[i]] <- dat
}

fire_dat <- do.call(rbind, fire_dat)

# group 
fire_dat_plot <- fire_dat %>%
  filter(Year > 2020, 
         Year <= 2100)

write_csv(fire_dat_plot, paste0(path_results, "/fire_results_fut_final.csv"))

fire_dat_plot <- fire_dat_plot %>% 
  group_by(sim, year, scen) %>% 
  summarize(burned_area = sum(realized_size)) %>% 
  mutate(period = ifelse(year <= 20, "2021-2040",
                         ifelse(year <= 40, "2041-2060",
                                ifelse(year <= 60, "2061-2080", "2081-2100"))))


# bring historical data to same format
a_fire_plot <- a_fire %>% 
  mutate(sim = "0",
         period = "1986-2020",
         year = year - 1985) %>% 
  dplyr::select(sim, year, burned_area = area_tot, period)

hist_fire <- rbind(a_fire_plot %>% mutate(scen = "RCP2.6"),
                   a_fire_plot %>% mutate(scen = "RCP4.5"),
                   a_fire_plot %>% mutate(scen = "RCP8.5"))


write_csv(a_fire_plot, paste0(path_results, "/fire_results_hist_final.csv"))

# calculate historical mean for plottin
hist_mean <- mean(a_fire_plot$burned_area)

# combine historical and simulation data
fire_dat_plot_hist <- rbind(fire_dat_plot, hist_fire)
write_csv(fire_dat_plot_hist, paste0(path_results, "/fire_results_final.csv"))






# wind ---
wind_dat <- list()
for(i in c(1:90)){
  
  dat <- read_csv(paste0(path_sims, "/output_sim_eu_", i, "/wind_europe_sim_eu_", i, ".csv"))
  dat <- dat %>% dplyr::select(year, cells_affected) %>% 
    mutate(Year = as.numeric(year) + 2019) %>% 
    mutate(sim = paste(i),
           scen = ifelse(i %in% seq(1, 90, 3), "rcp85", 
                         ifelse(i %in% seq(2, 90, 3), "rcp45", "rcp26"))) %>% 
    mutate(scen = ifelse(scen == "rcp85", "RCP8.5",
                         ifelse(scen == "rcp45", "RCP4.5", "RCP2.6")))
  
  wind_dat[[i]] <- dat
}

wind_dat <- do.call(rbind, wind_dat)

wind_dat_plot <- wind_dat %>%
  filter(Year > 2020, 
         Year <= 2100)

write_csv(wind_dat_plot, paste0(path_results, "/wind_results_fut_final.csv"))


wind_dat_plot <- wind_dat_plot %>% 
  group_by(sim, year, scen) %>% 
  summarize(dist_area = sum(cells_affected)) %>% 
  mutate(period = ifelse(year <= 20, "2021-2040",
                         ifelse(year <= 40, "2041-2060",
                                ifelse(year <= 60, "2061-2080", "2081-2100"))))


a_wind_new_plot <- a_wind_new %>% 
  mutate(sim = "0",
         period = "1986-2020",
         year = year - 1985) %>% 
  dplyr::select(sim, year, dist_area = area_tot, period)

write_csv(a_wind_new_plot, paste0(path_results, "/wind_results_hist_final.csv"))


hist_wind <- rbind(a_wind_new_plot %>% mutate(scen = "RCP2.6"),
                   a_wind_new_plot %>% mutate(scen = "RCP4.5"),
                   a_wind_new_plot %>% mutate(scen = "RCP8.5"))

hist_mean <- mean(a_wind_new_plot$dist_area)

wind_dat_plot_hist <- rbind(wind_dat_plot, hist_wind)
write_csv(wind_dat_plot_hist, paste0(path_results, "/wind_results_final.csv"))



### bark beetle  ---
bbtl_dat <- list()
for(i in c(1:90)){
  
  dat <- read_csv(paste0(path_sims, "/output_sim_eu_", i, "/barkbeetle.csv"))
  dat <- dat %>% dplyr::select(year, n_impact) %>% 
    mutate(Year = as.numeric(year) + 2019) %>% 
    mutate(sim = paste(i),
           scen = ifelse(i %in% seq(1, 90, 3), "rcp85", 
                         ifelse(i %in% seq(2, 90, 3), "rcp45", "rcp26"))) %>% 
    mutate(scen = ifelse(scen == "rcp85", "RCP8.5",
                         ifelse(scen == "rcp45", "RCP4.5", "RCP2.6")))
  
  bbtl_dat[[i]] <- dat
}

bbtl_dat <- do.call(rbind, bbtl_dat)

bbtl_dat_plot <- bbtl_dat %>%
  filter(Year > 2020, 
         Year <= 2100)


write_csv(bbtl_dat_plot, paste0(path_results, "/bbtl_results_fut_final.csv"))


bbtl_dat_plot <- bbtl_dat_plot %>% 
  group_by(sim, year, scen) %>% 
  summarize(dist_area = sum(n_impact)) %>% 
  mutate(period = ifelse(year <= 20, "2021-2040",
                         ifelse(year <= 40, "2041-2060",
                                ifelse(year <= 60, "2061-2080", "2081-2100"))))


a_bbtl_new_plot <- a_bb_new %>% 
  mutate(sim = "0",
         period = "1986-2020",
         year = year - 1985) %>% 
  dplyr::select(sim, year, dist_area = area_tot, period)

write_csv(a_bbtl_new_plot, paste0(path_results, "/bbtl_results_hist_final.csv"))


hist_bbtl <- rbind(a_bbtl_new_plot %>% mutate(scen = "RCP2.6"),
                   a_bbtl_new_plot %>% mutate(scen = "RCP4.5"),
                   a_bbtl_new_plot %>% mutate(scen = "RCP8.5"))

hist_mean <- mean(a_bbtl_new_plot$dist_area)

bbtl_dat_plot_hist <- rbind(bbtl_dat_plot, hist_bbtl)
write_csv(bbtl_dat_plot_hist, paste0(path_results, "/bbtl_results_final.csv"))



### bring the agents together ----------------------------------
# 
# veg_init <- rast(paste0(path, "/initial_states/init_veg_v2.tif"))
# veg_init[!is.na(veg_init)] <- 1
# total_forest_area <- sum(values(veg_init), na.rm = T)

total_forest_area <- 186739717
historical_area <- 207334385 # different because different resolution (30m)

bbtl_dat_plot_hist <- read_csv(paste0(path_results, "/bbtl_results_final.csv"))
fire_dat_plot_hist <- read_csv(paste0(path_results, "/fire_results_final.csv"))
wind_dat_plot_hist <- read_csv(paste0(path_results, "/wind_results_final.csv"))

bbtl_dat_plot <- bbtl_dat_plot_hist %>% mutate(agent = paste("Barkbeetle"))
fire_dat_plot <- fire_dat_plot_hist %>% mutate(agent = paste("Fire")) %>% rename("dist_area" = "burned_area") 
wind_dat_plot <- wind_dat_plot_hist %>% mutate(agent = paste("Wind"))
plot_dat_agents <- rbind(wind_dat_plot, fire_dat_plot, bbtl_dat_plot)

write_csv(plot_dat_agents, paste0(path_results, "/all_agents_distarea_final.csv"))


### numbres all --- 

# disturbed area

numbers_change <- plot_dat_agents %>%
  group_by(year, sim, scen, period) %>% 
  summarize(dist_area = sum(dist_area)) %>% 
  group_by(year, scen, period) %>% 
  mutate(dist_area = mean(dist_area))

hist_dat <- numbers_change %>%
  filter(period %in% c("1986-2000", "2001-2020")) %>% 
  mutate(dist_rate = 100/historical_area*dist_area,
         dist_freq = historical_area/dist_area) %>% 
  mutate(Year = year + 1985)

fut_dat <- numbers_change %>%
  filter(!period %in% c("1986-2000", "2001-2020")) %>%  
  mutate(dist_rate = 100/total_forest_area*dist_area,
         dist_freq = total_forest_area/dist_area) %>% 
  mutate(Year = year + 2020) %>% 
  filter(Year <= 2100)


df_combined_sub <- rbind(hist_dat, fut_dat)

numbers_change <- df_combined_sub %>% drop_na() %>% 
  group_by(period, scen) %>% 
  summarise(mean_dist_rate = mean(dist_rate),
            sd = sd(dist_rate), 
            n = n(),
            SE = sd/sqrt(n),
            mean_dist_freq = mean(dist_freq),
            sd_freq = sd(dist_freq),
            SE_freq = sd_freq/sqrt(n)) %>% 
  mutate(ci_low = mean_dist_rate - 1.96 * SE,
         ci_high = mean_dist_rate + 1.96 * SE,
         ci_low_freq = mean_dist_freq - 1.96 * SE_freq,
         ci_high_freq = mean_dist_freq + 1.96 * SE_freq)

df_combined_sub_hist_21th <- as.vector(unlist(numbers_change[numbers_change$period == "1986-2020", "mean_dist_rate"]))[1]
df_combined_sub_hist_21th_freq <- as.vector(unlist(numbers_change[numbers_change$period == "1986-2020", "mean_dist_freq"]))[1]

numbers_change <- numbers_change %>% 
  mutate(dist_rate_change = 100/df_combined_sub_hist_21th*mean_dist_rate,
         dist_freq_change = 100/df_combined_sub_hist_21th_freq*mean_dist_freq)

View(numbers_change)

write_csv(numbers_change, paste0(path_results, "/disturbance_change_final.csv"))


###  end ---------------------------------------------------------------------------
