# --------------------------------------------------------------
# Fire Module Scipt 03: statistical model calibration for fire size and frequency
# --------------------------------------------------------------
# Author: Marc Grünig
# Date: 2025-01-07
# 
# This script calibrates statistical models for number of fires and fire size depending on
# Vapor Pressure Deficit (VPD). The script uses previously calculated VPD values for 100km grid cells across Europe.
# Further, fire complexes are defined from Senf & Seidl 2021, and published in Grünig et al., 2023 (see below for link).

# The resulting models can be applied to derive a number of fires per year, depending on the aridity conditions,
# as well as the specific fire size, depending on aridity and location of a fire. These two models are two of the three core elements 
# to obtain a series of fire events across Europe which are used in the SVD fire module.

# Notes:
# - The underlying data for the fire complexes can be found here https://zenodo.org/records/7386862

# --------------------------------------------------------------


### packages ----------------------------------------------------
library(raster)
library(terra)
library(sf)
library(rgdal)
library(tidyverse)
library(sf)
library(gridExtra)
library(stars)
library(lme4)
library(dplyr)
library(glmmTMB)
library(data.table)
library(MuMIn)
library(DHARMa)

### settings -------------------------------
path <- "/.../"

# Raster package
options(dplyr.summarise.inform = FALSE)
tmpFiles(remove=TRUE, current=T, orphan=T)
removeTmpFiles()


### load fire complexes -----------------------------------------------
cntrs <- list.files(paste0(path, "/fire_module/complexes"), ".gpkg$") %>% 
  str_sub(., 14, nchar(.) - 5)
cntrs <- cntrs[grepl("*eps150", cntrs)] %>% str_sub(., 0, nchar(.) - 7)

k <- 0
out <- vector("list", length(cntrs))

cntrs.fullnam <- list.files(paste0(path, "/fire_module/complexes"), ".gpkg$", full.names = TRUE)
cntrs.fullnam <- cntrs.fullnam[grepl("*eps150", cntrs.fullnam)]

for (i in cntrs.fullnam) {
  
  k <- k + 1
  
  complexes <- read_sf(i) %>%
    mutate(area = st_area(.)) %>%
    mutate(area = as.integer(area))
  
  centroids <- complexes %>% 
    st_centroid() %>% 
    st_geometry()
  
  complexes <- complexes %>% st_drop_geometry()
  
  out[[k]] <- cbind(complexes, centroids)
  
}

nfire <- out %>% map(nrow) %>% unlist()

a <- out[which(nfire > 0)] %>%
  set_names(cntrs[which(nfire > 0)]) %>%
  bind_rows()


### Group disturbances to model area and number of storms per year -----------------------------------------------------
ref_grid_100km <- rast(paste0(path, "/reference_grids/reference_grid_100km.tif"))
ref_grid_10km <- rast(paste0(path, "/reference_grids/reference_grid_10km.tif"))

gridid <- ref_grid_100km
gridid_vect <- as.polygons(ref_grid_100km)


years <- 1986:2020
d_sf <- a %>% 
  filter(year %in% years) %>%
  mutate(x = unlist(map(a$geometry,1)),
         y = unlist(map(a$geometry,2))) %>%
  st_as_sf(
    .,
    coords = c("x", "y"),
    crs = st_crs(ref_grid_100km)
  )

k <- 0
dat_annual <- vector("list", length(years))

for (y in years) {
  
  print(y)
  
  k <- k + 1
  
  gridid_sel <- terra::extract(gridid, d_sf %>% filter(year == y)) %>%
    .$gridid
  gridid_sel <- as.data.frame(cbind(d_sf %>% filter(year == y), gridid = gridid_sel)) %>% 
    dplyr::select(clust, year, complex_size_m2, gridid)
  
  dat_annual[[k]] <- gridid_sel
  
}

dat_annual_summary <- dat_annual %>%
  bind_rows()


n_fires <- dat_annual_summary %>%
  group_by(year) %>% 
  summarise(n = n()) %>%
  ungroup()


ggplot(data = n_fires, 
       aes(x = year, y = n)) +
  geom_bar(stat = "identity") +
  labs(x = "fires",
       y = "year")

dat_new <- dat_annual_summary %>%
  group_by(gridid, year) %>% 
  summarize(complex_size_m2_max = max(complex_size_m2),
            n_fires = n())



### combine with vpd data -----------------------------------------------------------------
vpd_hist <- read.csv(paste0(path, "/fire_module/climate/hist_vpd_100km.csv"))
colnames(vpd_hist) <- c("gridid", 1986:2020)
vpd_hist <- vpd_hist %>% 
  pivot_longer(`1986`:`2020`, names_to = "year", values_to = "vpd_summer") %>% 
  group_by(gridid, year) %>% 
  summarize(vpd_summer = mean(vpd_summer)) %>% 
  mutate(year = as.numeric(year))

full_df <- dat_new %>% 
  left_join(., vpd_hist, by = c("gridid", "year"))

dat_model_size <- full_df %>%
  dplyr::select(year, gridid, complex_size_m2_max, n_fires, vpd_summer) %>%
  dplyr::rename("vpd_summer_mean_rollmax" = vpd_summer) %>% 
  drop_na()

gridid_df <- as.data.frame(ref_grid_100km)
gridid_df <- expand.grid(unlist(as.vector(gridid_df)), years)
colnames(gridid_df) <- c("gridid", "year")

full_df <- gridid_df %>% 
  left_join(., dat_new, by = c("gridid", "year")) %>% 
  left_join(., vpd_hist, by = c("gridid", "year"))

full_df[is.na(full_df$n_fires), "n_fires"] <- 0
full_df[is.na(full_df$complex_size_m2_max), "complex_size_m2_max"] <- 0

dat_model_freq <- full_df %>%
  dplyr::select(year, gridid, complex_size_m2_max, n_fires, vpd_summer) %>%
  mutate(vpd_summer = ifelse(vpd_summer == Inf | vpd_summer == -Inf, NA, vpd_summer)) %>% 
  mutate(complex_size_m2_max = ifelse(complex_size_m2_max == Inf | complex_size_m2_max == -Inf, NA, complex_size_m2_max)) %>% 
  mutate(n_fires = ifelse(n_fires == Inf | n_fires == -Inf, NA, n_fires)) %>% 
  dplyr::rename("vpd_summer_mean_rollmax" = vpd_summer) %>% 
  drop_na()


# Check for infinite values in specific columns
columns_with_infinite <- sapply(dat_model_freq, function(x) any(is.infinite(x)))
if (any(columns_with_infinite)) {
  print("There are infinite values in the following columns:")
  print(names(columns_with_infinite)[columns_with_infinite])
}


###  modelling -----------------------------------------------------------------------------


###  modelling -----------------------------------------------------------------------------

# calibrate model
mod_freq <- glmmTMB((n_fires) ~ (vpd_summer_mean_rollmax) + (1 + (vpd_summer_mean_rollmax) | year),
                    data = dat_model_freq,  family = "nbinom2", ziformula= ~ 1 + (vpd_summer_mean_rollmax))

# visualize the data
ggplot(dat_model_freq , 
       aes(x = vpd_summer_mean_rollmax, 
           y = log10(n_fires))) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  theme(legend.position = "none") +
  facet_wrap(~year)

# check r2
r.squaredGLMM(mod_freq)


# some model checking
# simulationOutput <- simulateResiduals(fittedModel = mod_freq)
# testDispersion(simulationOutput)
# plot(simulationOutput)
# res <- simulateResiduals(mod_freq, plot = T)
# testZeroInflation(res)

# save the model
saveRDS(mod_freq, paste0(path, "/fire_module/models/final_freq_lmer_100km.rds"))


# use the model to predict number of fires for all simulation years to test the prediction ability
fires_dat <- list()
for(n in 1:30){
  
  fires <- c()
  for(i in 1986:2020){
    
    print(i)
    test_dat <- dat_model_freq %>%
      filter(year == i) %>%
      dplyr::select(year, vpd_summer_mean_rollmax)
    
    
    # according to https://stats.stackexchange.com/questions/189005/simulate-from-a-zero-inflated-poisson-distribution
    lambda <- predict(mod_freq, newdata = test_dat, type = "conditional",  allow.new.levels = T) # allow.new.levels = T
    p <- predict(mod_freq, newdata = test_dat, type="zprob",  allow.new.levels = T)
    # pred_new <- VGAM::rzipois(nrow(test_dat), lambda = lambda, pstr0 = p)
    pred_new <- emdbook::rzinbinom(nrow(test_dat),
                                   mu=lambda, size=sigma(mod_freq),
                                   zprob=p)
    print(sum(pred_new))
    fires <- c(fires, sum(pred_new))
  }
  
  
  ## future preds
  fires_fut <- c()
  for(y in 2021:2100){
    
    print(y)
    clim <- fread(paste0(path, "/fire_module/climate/fut_data_vpd_100km_ncc_rcp85.csv"), sep =",")
    colnames(clim) <- gsub("X", "", colnames(clim))
    
    test_dat <- clim %>% dplyr::select(paste0(y), gridid) %>% drop_na() %>% 
      dplyr::rename("vpd_summer_mean_rollmax" = paste0(y)) %>% 
      group_by(gridid) %>% 
      summarise(vpd_summer_mean_rollmax = mean(vpd_summer_mean_rollmax)) %>% 
      ungroup() %>% 
      mutate(year = paste0(y)) %>% 
      dplyr::select(year, vpd_summer_mean_rollmax)
    
    
    # according to https://stats.stackexchange.com/questions/189005/simulate-from-a-zero-inflated-poisson-distribution
    lambda <- predict(mod_freq, newdata = test_dat, type = "conditional",  allow.new.levels = T) # allow.new.levels = T
    p <- predict(mod_freq, newdata = test_dat, type="zprob",  allow.new.levels = T)
    pred_new <- emdbook::rzinbinom(nrow(test_dat),
                                   mu=lambda, size=sigma(mod_freq),
                                   zprob=p)
    print(sum(pred_new))
    fires_fut <- c(fires_fut, sum(pred_new))
  }
  
  fires_dat[[n]] <- as.data.frame(cbind(year = 1986:2100, n_fires = c(fires, fires_fut), sim_nr = n))
}


fires_dat <- do.call(bind_rows, fires_dat)

dat_new <- dat_model_freq %>% 
  group_by(year) %>% 
  summarise(n_fires = sum(n_fires),
            mean_size = mean(complex_size_m2_max),
            total_area = sum(complex_size_m2_max, na.rm = T)) %>% 
  dplyr::select(year, n_fires) %>% 
  mutate(sim_nr = 0)

plot_dat <- fires_dat %>% bind_rows(dat_new) %>%
  rowwise() %>% 
  mutate(dat = ifelse(sim_nr == 0, "real", "sim"))

ggplot(plot_dat, aes(x = year, y = n_fires, color = dat, group = sim_nr)) + 
  geom_line(aes(alpha = 0.2)) +
  geom_line(data = plot_dat %>% filter(dat == "real"), aes(alpha = 1, color = "real")) +
  scale_alpha_identity() +
  scale_color_manual(values = c("real" = "red", "sim" = "lightblue")) # Add color



### modelling fire size ----------------------------------------------------------------------------------

# this part follows the model that was published in Grünig et al., 2023  (https://doi.org/10.1111/gcb.16547)
mod_size <- lmer(log10(complex_size_m2_max) ~ log10(vpd_summer_mean_rollmax) + (1 + log10(vpd_summer_mean_rollmax) | gridid),
                  data = dat_model_size)

# save model
saveRDS(mod_size,  paste0(path, "/fire_module/models/final_size_lmer_100km.rds"))

### -----------------------------------------------------------------------------------
