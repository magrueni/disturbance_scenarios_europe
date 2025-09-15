# --------------------------------------------------------------------------------------
# Script Name: module_evaluation
# Description: This scripts processes the simulations under historical climate conditions
#              and compares the outputs with historical disturbance observations
#              from remote sensing data.
#
# 
# Author: Marc Grünig
# Date: 7.2.2025
#
# Output:
#   - Figures and statistics to evaluate the disturbance modules 
#   - Figures S30 - S34
#
# Notes:
# - underlying disturbance data can be obtained here: https://zenodo.org/records/4607230
# --------------------------------------------------------------------------------------


### libraries ------------------------------------------------------------------
library(ggplot2)
library(tidyverse)
library(dplyr)
library(MetBrewer)
library(sf)
library(stars)
library(terra)
library(ggtext)
library(ggh4x)


# paths
path <- "/.../"
path_public <- "/.../"

# define colors
cols <- rev(met.brewer(name = "Gauguin", n = 3, type="discrete"))



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
grid <- read_sf(paste0(path, "/firesize/data/climate/climategrid_epsg3035_090122.gpkg"))
grid_raster <- st_rasterize(grid)
grid_r <- rast(grid_raster)

d_sf <- a %>% 
  st_as_sf(., coords = c("x", "y"), crs = st_crs(grid))

fire_vect <- vect(d_sf)


a_fire <- as.data.frame(fire_vect) %>% 
  dplyr::select(year, complex_size_m2, complex_area_burnt_m2) %>% 
  group_by(year) %>% 
  summarise(n_fires = n(),
            area_tot  = sum(complex_area_burnt_m2*0.0001),
            area_tot_planned  = sum(complex_size_m2*0.0001)) %>% 
  mutate(dat = "real")

fire_df <- as.data.frame(fire_vect) %>% 
  dplyr::select(year, complex_size_m2, complex_area_burnt_m2) %>% 
  group_by(year) %>% 
  summarise(n_fires = n(),
            burned_area_realized = sum(complex_area_burnt_m2*0.0001),
            burned_area_planned = sum(complex_size_m2*0.0001)) %>% 
  mutate(sim_run = 0) %>% 
  mutate(ratio = burned_area_realized/burned_area_planned)

hist(fire_df$ratio)
mean(fire_df$burned_area_realized)
mean(fire_df$burned_area_planned)



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


a_fire <- a %>%
  filter(prediction_class == "fire") %>%
  dplyr::select(country, patch, year, area_m2 = area) %>%
  group_by(year) %>% 
  summarise(n_storms = n(),
            area_tot = sum(area_m2) * 0.0001) %>% 
  mutate(dat = "real")

mean(a_fire$area_tot, na.rm = T)



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
  #filter(country %in% c("germany", "sweden", "finland")) %>% 
  group_by(year) %>% 
  summarise(n_storms = n(),
            area_tot = sum(area_m2) * 0.0001) %>% 
  mutate(dat = "real") 

# load Patacca data and calculate percentages to split up the remote sensing data
dat <- read.csv(paste0(path, "/disturbance_timeseries_database.csv"), sep = ";")

dat_filt <- dat %>% 
  filter(Year %in% c(1986:2020)) %>% 
  filter(disturbance_driver == "Bark Beetles") %>% 
  dplyr::select(year = Year, Reported_m3) %>% 
  group_by(year) %>% 
  summarize(Reported_m3 = sum(Reported_m3, na.rm = T))

dat_summed <- dat %>% 
  filter(Year %in% c(1986:2020)) %>% 
  filter(disturbance_driver %in% c("Bark Beetles", "Wind")) %>% 
  dplyr::select(year = Year, Reported_m3) %>% 
  group_by(year) %>% 
  summarize(tot = sum(Reported_m3, na.rm = T))

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
  group_by(year, country) %>% 
  summarise(n_storms = n(),
            area_tot = sum(area_m2) * 0.0001) %>% 
  mutate(dat = "real")




### now load simulation data ---------------------------------------------------

path_sims <- paste0(path, "/09_svd_simulations/raw/")


### fire ---
fire_dat <- list()
for(i in c(91:120)){
  
  if(!file.exists(paste0(path_sims, "/output_sim_eu_", i, "/ccfire/fire_run_sim_eu_", i, ".csv"))){next}
  dat <- read_csv(paste0(path_sims, "/output_sim_eu_", i, "/ccfire/fire_run_sim_eu_", i, ".csv"))
  dat <- dat %>% dplyr::select(year, realized_size) %>% 
    mutate(sim = paste(i),
           scen = ifelse(i %in% seq(1, 90, 3), "rcp85", 
                         ifelse(i %in% seq(2, 90, 3), "rcp45",
                                ifelse(i %in% seq(3, 90, 3), "rcp26",
                                       "Historical simulations"))))
  fire_dat[[i]] <- dat
}

fire_dat <- do.call(rbind, fire_dat)

fire_dat_plot <- fire_dat %>% group_by(sim, year, scen) %>% 
  summarize(dist_area = sum(realized_size)) %>% 
  mutate(period = "1981-2010")


a_fire_new_plot <- a_fire %>% 
  mutate(period = ifelse(year <= 2000, "1986-2000", "2001-2020")) %>% 
  mutate(sim = "0",
         year = year - 1985) %>% 
  dplyr::select(sim, year, dist_area = area_tot, period) %>% 
  mutate(scen = "Historical data")


hist_fire <- a_fire_new_plot
hist_mean <- mean(a_fire_new_plot$dist_area, na.rm = T)

fire_dat_plot_hist <- rbind(fire_dat_plot, hist_fire)

test <- fire_dat_plot_hist
mean(test$dist_area, na.rm = T)

fire_dat_plot_hist %>% group_by(scen) %>% 
  summarize(mean = mean(dist_area, na.rm = T))
fire_dat_plot_hist <- fire_dat_plot_hist %>% 
  mutate(period = ifelse(period != "1981-2010", "1986-2020", "1981-2010"))

# save results
write_csv(wind_dat_plot_hist, paste0(path, "/11_figures/figure_data/figS31_fire_eval.csv"))


# Create a density plot for burned area grouped by sim, colored and filled by scen
p2_fire <- ggplot(fire_dat_plot_hist, aes(x = as.factor(period), y = (dist_area), fill = scen, color = scen)) +
  geom_boxplot(alpha = 0.2, linewidth = 1, outlier.shape = NA) +
  geom_boxplot(data = subset(fire_dat_plot_hist, period == "2000-2020"), aes(fill = "Historical", color = "Historical"), alpha = 0.5, linewidth = 1, outlier.shape = NA) +
  geom_jitter(data = subset(fire_dat_plot_hist, period != "2000-2020"), position = position_jitter(width = 0.2, height = 0), alpha = 0.5, aes(color = scen, fill = scen)) +
  geom_jitter(data = subset(fire_dat_plot_hist, period == "2000-2020"), aes(fill = "Historical", color = "Historical"), alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point", col = "red", show.legend = FALSE) +  stat_summary(fun = "mean", geom = "point", col = "red", show.legend = FALSE) +
  labs(title = "Annual burned area",
       x = "Period",
       y = "Annual burned area [ha]") +
  scale_color_manual(values = c("Historical simulations" = "navy",
                                "Historical data" = "lightblue"),
                     name = "Scenario",
                     breaks = c("Historical simulations", "Historical data"),
                     labels = c("Historical simulations", "Historical data")) +
  scale_fill_manual(values = c("Historical simulations" = "navy",
                               "Historical data" = "lightblue"),
                    name = "Scenario",
                    breaks = c("Historical simulations", "Historical data"),
                    labels = c("Historical simulations", "Historical data")) +
  theme_bw() +
  geom_point(aes(shape = "mean"), alpha = 0, col = "red")+  # <-- added this line of code and next
  guides(shape=guide_legend(title=NULL, override.aes = list(alpha = 1)))


p2_fire
ggsave(p2_fire, filename = paste0(path_results, "/11_figures/FigS31_fire.png"), width = 8, height = 6)

# numbers
dat_real <- subset(fire_dat_plot_hist, period != "1981-2010")
mean(dat_real$dist_area)
median(dat_real$dist_area)
sd(dat_real$dist_area)

quantile(dat_real$dist_area, 0.05)
quantile(dat_real$dist_area, 0.95)


dat_simu <- subset(fire_dat_plot_hist, period == "1981-2010")
mean(dat_simu$dist_area)
median(dat_simu$dist_area)
sd(dat_simu$dist_area)
quantile(dat_simu$dist_area, 0.05)
quantile(dat_simu$dist_area, 0.95)



### wind ---
wind_dat <- list()
for(i in c(91:120)){
  
  if(!file.exists(paste0(path_sims, "/output_sim_eu_", i, "/ccwind/wind_europe_sim_eu_", i, ".csv"))){next}
  dat <- read_csv(paste0(path_sims, "/output_sim_eu_", i, "/ccwind/wind_europe_sim_eu_", i, ".csv"))
  dat <- dat %>% dplyr::select(year, cells_affected) %>% 
    mutate(sim = paste(i),
           scen = ifelse(i %in% seq(1, 90, 3), "rcp85", 
                         ifelse(i %in% seq(2, 90, 3), "rcp45",
                                ifelse(i %in% seq(3, 90, 3), "rcp26",
                                       "Historical simulations"))))
  wind_dat[[i]] <- dat
}

wind_dat <- do.call(rbind, wind_dat)

wind_dat_plot <- wind_dat %>% group_by(sim, year, scen) %>% 
  summarize(dist_area = sum(cells_affected)) %>% 
  mutate(period = "1981-2010")


a_wind_new_plot <- a_wind_new %>% 
  # filter(year > 1999) %>% 
  mutate(period = ifelse(year <= 2000, "1986-2000", "2001-2020")) %>% 
  mutate(sim = "0",
         year = year - 1985) %>% 
  dplyr::select(sim, year, dist_area = area_tot, period) %>% 
  mutate(scen = "Historical data")



hist_wind <- a_wind_new_plot
hist_mean <- mean(a_wind_new_plot$dist_area, na.rm = T)

wind_dat_plot_hist <- rbind(wind_dat_plot, hist_wind)

test <- wind_dat_plot_hist
mean(test$dist_area, na.rm = T)

wind_dat_plot_hist %>% group_by(period, scen) %>% 
  summarize(mean = mean(dist_area, na.rm = T))


second <- a_wind_new_plot %>% filter(period == "2001-2020")
mean(second$dist_area, na.rm = T)
wind_dat_plot_hist <- rbind(wind_dat_plot, hist_wind)

test <- wind_dat_plot_hist
mean(test$dist_area, na.rm = T)

wind_dat_plot_hist %>% group_by(period, scen) %>% 
  summarize(mean = mean(dist_area, na.rm = T))
write_csv(wind_dat_plot_hist, paste0(path, "/11_figures/figure_data/figS32_wind_eval.csv"))


# Create a density plot for burned area grouped by sim, colored and filled by scen
p2_wind <- ggplot(wind_dat_plot_hist, aes(x = as.factor(period), y = (dist_area), fill = scen, color = scen)) +
  geom_boxplot(alpha = 0.2, linewidth = 1, outlier.shape = NA) +
  geom_boxplot(data = subset(wind_dat_plot_hist, period == "2000-2020"), aes(fill = "Historical", color = "Historical"), alpha = 0.5, linewidth = 1, outlier.shape = NA) +
  geom_jitter(data = subset(wind_dat_plot_hist, period != "2000-2020"), position = position_jitter(width = 0.2, height = 0), alpha = 0.5, aes(color = scen, fill = scen)) +
  geom_jitter(data = subset(wind_dat_plot_hist, period == "2000-2020"), aes(fill = "Historical", color = "Historical"), alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point", col = "red", show.legend = FALSE) +  stat_summary(fun = "mean", geom = "point", col = "red", show.legend = FALSE) +
  labs(title = "Annual windthrow area",
       x = "Period",
       y = "Annual windthrow area [ha]") +
  scale_color_manual(values = c("Historical simulations" = "navy",
                                "Historical data" = "lightblue"),
                     name = "Scenario",
                     breaks = c("Historical simulations", "Historical data"),
                     labels = c("Historical simulations", "Historical data")) +
  scale_fill_manual(values = c("Historical simulations" = "navy",
                               "Historical data" = "lightblue"),
                    name = "Scenario",
                    breaks = c("Historical simulations", "Historical data"),
                    labels = c("Historical simulations", "Historical data")) +
  theme_bw() +
  geom_point(aes(shape = "mean"), alpha = 0, col = "red")+  # <-- added this line of code and next
  guides(shape=guide_legend(title=NULL, override.aes = list(alpha = 1)))


p2_wind
ggsave(p2_wind, filename = paste0(path_results, "/11_figures/FigS32_wind.png"), width = 8, height = 6)


# numbers
dat_real <- subset(wind_dat_plot_hist, period == "1986-2000")
mean(dat_real$dist_area, na.rm = T)
median(dat_real$dist_area, na.rm = T)
sd(dat_real$dist_area, na.rm = T)

quantile(dat_real$dist_area, 0.05, na.rm = T)
quantile(dat_real$dist_area, 0.95, na.rm = T)


dat_real <- subset(wind_dat_plot_hist, period == "2001-2020")
mean(dat_real$dist_area, na.rm = T)
median(dat_real$dist_area, na.rm = T)
sd(dat_real$dist_area, na.rm = T)

quantile(dat_real$dist_area, 0.05, na.rm = T)
quantile(dat_real$dist_area, 0.95, na.rm = T)



dat_simu <- subset(wind_dat_plot_hist, period == "1981-2010")
mean(dat_simu$dist_area, na.rm = T)
median(dat_simu$dist_area, na.rm = T)
sd(dat_simu$dist_area, na.rm = T)
quantile(dat_simu$dist_area, 0.05, na.rm = T)
quantile(dat_simu$dist_area, 0.95, na.rm = T)


### bark beetle  ---
bbtl_dat <- list()
for(i in c(91:105)){
  
  if(!file.exists(paste0(path, "/output_sim_eu_", i, "/barkbeetle.csv"))){next}
  dat <- read_csv(paste0(path, "/output_sim_eu_", i, "/barkbeetle.csv"))
  if(nrow(dat) == 0){next}
  dat <- dat %>% dplyr::select(year, n_impact) %>% 
    mutate(sim = paste(i),
           scen = ifelse(i %in% seq(1, 90, 3), "rcp85", 
                         ifelse(i %in% seq(2, 90, 3), "rcp45",
                                ifelse(i %in% seq(3, 90, 3), "rcp26",
                                       "Historical simulations"))))
  
  bbtl_dat[[i]] <- dat
}


bbtl_dat <- do.call(rbind, bbtl_dat)

bbtl_dat_plot <- bbtl_dat %>% group_by(sim, year, scen) %>% 
  summarize(dist_area = sum(n_impact)) %>% 
  mutate(period = "1981-2010")


a_bbtl_new_plot <- a_bb_new %>% 
  mutate(period = ifelse(year <= 2000, "1986-2000", "2001-2020")) %>% 
  mutate(sim = "0",
         year = year - 1985) %>% 
  dplyr::select(sim, year, dist_area = area_tot, period) %>% 
  mutate(scen = "Historical data")


hist_bbtl <- a_bbtl_new_plot
hist_mean <- mean(a_bbtl_new_plot$dist_area, na.rm = T)

second <- a_bbtl_new_plot %>% filter(period == "2001-2020")
mean(second$dist_area, na.rm = T)
bbtl_dat_plot_hist <- rbind(bbtl_dat_plot, hist_bbtl)

test <- bbtl_dat_plot_hist
mean(test$dist_area, na.rm = T)

bbtl_dat_plot_hist %>% group_by(period, scen) %>% 
  summarize(mean = mean(dist_area, na.rm = T))


write_csv(bbtl_dat_plot_hist, paste0(path, "/11_figures/figure_data/figS33_bbtl_results.csv"))



# Create a density plot for burned area grouped by sim, colored and filled by scen
p2_bbtl <- ggplot(bbtl_dat_plot_hist, aes(x = as.factor(period), y = (dist_area), fill = scen, color = scen)) +
  geom_boxplot(alpha = 0.2, linewidth = 1, outlier.shape = NA) +
  geom_boxplot(data = subset(bbtl_dat_plot_hist, period == "2000-2020"), aes(fill = "Historical", color = "Historical"), alpha = 0.5, linewidth = 1, outlier.shape = NA) +
  geom_jitter(data = subset(bbtl_dat_plot_hist, period != "2000-2020"), position = position_jitter(width = 0.2, height = 0), alpha = 0.5, aes(color = scen, fill = scen)) +
  geom_jitter(data = subset(bbtl_dat_plot_hist, period == "2000-2020"), aes(fill = "Historical", color = "Historical"), alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point", col = "red", show.legend = FALSE) +  stat_summary(fun = "mean", geom = "point", col = "red", show.legend = FALSE) +
  labs(title = "Annual bark beetle disturbed area",
       x = "Period",
       y = "Annual disturbed area [ha]") +
  scale_color_manual(values = c("Historical simulations" = "navy",
                                "Historical data" = "lightblue"),
                     name = "Scenario",
                     breaks = c("Historical simulations", "Historical data"),
                     labels = c("Historical simulations", "Historical data")) +
  scale_fill_manual(values = c("Historical simulations" = "navy",
                               "Historical data" = "lightblue"),
                    name = "Scenario",
                    breaks = c("Historical simulations", "Historical data"),
                    labels = c("Historical simulations", "Historical data")) +
  theme_bw() +
  geom_point(aes(shape = "mean"), alpha = 0, col = "red")+  # <-- added this line of code and next
  guides(shape=guide_legend(title=NULL, override.aes = list(alpha = 1)))+
  theme(text = element_text(size = 12))

p2_bbtl
ggsave(p2_bbtl, filename = paste0(path_results, "/11_figures/S33_bbtl.png"), width = 8, height = 6)

# numbers
dat_real <- subset(bbtl_dat_plot_hist, period == "1986-2000")
mean(dat_real$dist_area, na.rm = T)
median(dat_real$dist_area, na.rm = T)
sd(dat_real$dist_area, na.rm = T)

quantile(dat_real$dist_area, 0.05, na.rm = T)
quantile(dat_real$dist_area, 0.95, na.rm = T)


dat_real <- subset(bbtl_dat_plot_hist, period == "2001-2020")
mean(dat_real$dist_area, na.rm = T)
median(dat_real$dist_area, na.rm = T)
sd(dat_real$dist_area, na.rm = T)

quantile(dat_real$dist_area, 0.05, na.rm = T)
quantile(dat_real$dist_area, 0.95, na.rm = T)



dat_simu <- subset(bbtl_dat_plot_hist, period == "1981-2010")
mean(dat_simu$dist_area, na.rm = T)
median(dat_simu$dist_area, na.rm = T)
sd(dat_simu$dist_area, na.rm = T)
quantile(dat_simu$dist_area, 0.05, na.rm = T)
quantile(dat_simu$dist_area, 0.95, na.rm = T)



###management ---------------------------
mgmt_dat <- list()
for(i in c(91:105)){
  
  if(!file.exists(paste0(path, "/output_sim_eu_", i, "/management_sim_eu_", i, ".csv"))){next}
  dat <- read_csv(paste0(path, "/output_sim_eu_", i, "/management_sim_eu_", i, ".csv"))
  if(nrow(dat) == 0){next}
  dat <- dat %>% 
    group_by(year) %>% 
    summarize(dist_area = sum(n, na.rm = T)) %>% 
    mutate(sim = paste(i),
           scen = ifelse(i %in% seq(1, 90, 3), "rcp85", 
                         ifelse(i %in% seq(2, 90, 3), "rcp45",
                                ifelse(i %in% seq(3, 90, 3), "rcp26",
                                       "Historical simulations"))))
  
  mgmt_dat[[i]] <- dat
}


mgmt_dat <- do.call(rbind, mgmt_dat)

mgmt_dat_plot <- mgmt_dat %>% group_by(sim, year, scen) %>% 
  summarize(dist_area = sum(dist_area)) %>% 
  mutate(period = "1981-2010")


a_mgmt_new_plot <- a_mgmt %>% 
  group_by(year) %>% 
  summarise(area_tot = sum(area_tot, na.rm = T)) %>% 
  # filter(year > 1999) %>% 
  mutate(sim = "0",
         period = "1986-2020",
         year = year - 1985) %>% 
  dplyr::select(sim, year, dist_area = area_tot, period) %>% 
  mutate(scen = "Historical data")


hist_mgmt <- a_mgmt_new_plot
hist_mean <- mean(a_mgmt_new_plot$dist_area, na.rm = T)

mgmt_dat_plot_hist <- rbind(mgmt_dat_plot, hist_mgmt)

test <- mgmt_dat_plot_hist
mean(test$dist_area, na.rm = T)

mgmt_dat_plot_hist %>% group_by(period, scen) %>% 
  summarize(mean = mean(dist_area, na.rm = T))


write_csv(bbtl_dat_plot_hist, paste0(path, "/11_figures/figure_data/figS34_bbtl_results.csv"))



# Create a density plot for burned area grouped by sim, colored and filled by scen
p2_bbtl <- ggplot(mgmt_dat_plot_hist, aes(x = as.factor(scen), y = (dist_area), fill = scen, color = scen)) +
  geom_boxplot(alpha = 0.2, linewidth = 1, outlier.shape = NA) +
  geom_boxplot(data = subset(mgmt_dat_plot_hist, period == "2000-2020"), aes(fill = "Historical", color = "Historical"), alpha = 0.5, linewidth = 1, outlier.shape = NA) +
  geom_jitter(data = subset(mgmt_dat_plot_hist, period != "2000-2020"), position = position_jitter(width = 0.2, height = 0), alpha = 0.5, aes(color = scen, fill = scen)) +
  geom_jitter(data = subset(mgmt_dat_plot_hist, period == "2000-2020"), aes(fill = "Historical", color = "Historical"), alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point", col = "red", show.legend = FALSE) +  stat_summary(fun = "mean", geom = "point", col = "red", show.legend = FALSE) +
  labs(title = "Annual harvest area",
       x = "Period",
       y = "Annual harvested area [ha]") +
  scale_color_manual(values = c("Historical simulations" = "navy",
                                "Historical data" = "lightblue"),
                     name = "Scenario",
                     breaks = c("Historical simulations", "Historical data"),
                     labels = c("Historical simulations", "Historical data")) +
  scale_fill_manual(values = c("Historical simulations" = "navy",
                               "Historical data" = "lightblue"),
                    name = "Scenario",
                    breaks = c("Historical simulations", "Historical data"),
                    labels = c("Historical simulations", "Historical data")) +
  theme_bw() +
  geom_point(aes(shape = "mean"), alpha = 0, col = "red")+  # <-- added this line of code and next
  guides(shape=guide_legend(title=NULL, override.aes = list(alpha = 1)))+
  theme(text = element_text(size = 12))

p2_bbtl
ggsave(p2_bbtl, filename = paste0(path_results, "/11_figures/FigS34_mgmt.png"), width = 8, height = 6)

# numbers
dat_real <- subset(mgmt_dat_plot_hist, period != "1981-2010")
mean(dat_real$dist_area, na.rm = T)
median(dat_real$dist_area, na.rm = T)
sd(dat_real$dist_area, na.rm = T)

dat_simu <- subset(mgmt_dat_plot_hist, period == "1981-2010")
mean(dat_simu$dist_area, na.rm = T)
median(dat_simu$dist_area, na.rm = T)
sd(dat_simu$dist_area, na.rm = T)

