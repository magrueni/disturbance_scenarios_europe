# --------------------------------------------------------------------------------------
# Script Name: Extreme events analysis
# Description: This scripts analyses the annual disturbance outputs from SVD simulations
# to calculate the frequency of extreme events
#
# Author: Marc Grünig
# Date: 7.2.2025
#
# Dependencies:
#   - R packages: dplyr, readr
#
# Input Files:
#   - processed annual outputs from script 04_annual_output_processing.R
#
# Output:
#   - Numbers on the extreme event frequency
#   - Data used for plotting figure 3b,c,d.
#
# ------------------------------------------------------------------------------




### libraries ------------------------------------------------------------------
library(dplyr)
library(readr)

path <- "/.../"

### load data ------------------------------------------------------------------

# master tab to get the gcm and rcp information
master_tab <- read.csv(paste0(path, "/svd_simulations/svd_simulations_ids.csv"), sep = ";")

# processed simulation outputs
fut_dat <- read_csv(paste0(path, "/svd_simulations/results_eu/fut_dist_annual_final.csv"))
hist_dat <- read_csv(paste0(path, "/svd_simulations/results_eu/hist_dist_annual_final.csv"))


### data processing ------------------------------------------------------------

### for fire ---

# load data and get the scenario information
fire_dat_plot_hist <- read_csv(paste0(path, "/svd_simulations/results_eu/fire_results_final.csv"))
fire_dat_plot_hist <- fire_dat_plot_hist %>%
  mutate(scen = ifelse(sim == "0", "Historical", scen),
         sim = ifelse(sim == "0", "99", sim)) %>% 
  distinct() %>% 
  rowwise() %>% 
  mutate(gcm = master_tab[as.numeric(sim), "GCM"],
         rep = master_tab[as.numeric(sim), "Rep"]) %>% 
  drop_na()

# then calculate extremes 
result_fire_historical <- fire_dat_plot_hist %>%
  filter(scen == "Historical") %>% 
  group_by(scen, gcm, rep) %>%
  summarize(mean_area = mean(burned_area , na.rm = TRUE),
            q90 = quantile(burned_area , 0.90, na.rm = TRUE),
            q99 = quantile(burned_area , 0.99, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(scen) %>% 
  summarize(mean_area = mean(mean_area),
            sd_mean = sd(mean_area),
            mean_q90 = mean(q90),
            sd_q90 = sd(q90),
            mean_q99 = mean(q99),
            sd_q99 = sd(q99))

result_fire <- fire_dat_plot_hist %>%
  filter(scen != "Historical") %>% 
  mutate(period = ifelse(year > 40, "p2", "p1")) %>% 
  group_by(scen, gcm, rep) %>%
  summarize(mean_area = mean(burned_area , na.rm = TRUE),
            q90 = quantile(burned_area , 0.90, na.rm = TRUE),
            q99 = quantile(burned_area , 0.99, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(scen) %>% 
  summarize(mean_area = mean(mean_area),
            sd_mean = sd(mean_area),
            mean_q90 = mean(q90),
            sd_q90 = sd(q90),
            mean_q99 = mean(q99),
            sd_q99 = sd(q99))


# check how much the q90 increases
result_fire %>% 
  mutate(percentual_difference = if_else(scen == "Historical", 0,
                                         ((mean_q90 - as.numeric(result_fire_historical[result_fire_historical$scen == "Historical", "mean_q90"])) /
                                            as.numeric(result_fire_historical[result_fire_historical$scen == "Historical", "mean_q90"])) * 100))


# create a summary table to see how often historical q90 was crossed
summary_table_fire <- fire_dat_plot_hist %>%
  filter(!period %in% c("1986-2020")) %>% 
  group_by(scen, gcm, rep) %>%
  mutate(n_tot = n()) %>% 
  filter(burned_area > as.numeric(result_fire_historical[result_fire_historical$scen == "Historical", "mean_q90"])) %>% 
  group_by(scen, gcm, rep) %>% 
  summarize(n_tot = mean(n_tot),
            n = n(),
            n_share = n/n_tot) %>% 
  group_by(scen) %>% 
  summarise(n_tot = sum(n),
            mean_n = mean(n),
            mean_n_share = mean(n_share),
            sd_n_share = sd(n_share)) %>% 
  mutate(mean_years = 1/mean_n_share,
         sd_years = mean_years*sd_n_share,
         se_years = sd_years / sqrt(n_tot),
         ci_year = se_years * 1.96)
summary_table_fire


summary_table_fire <- fire_dat_plot_hist %>%
  group_by(scen, gcm, rep, period) %>%
  mutate(n_tot = n()) %>% 
  filter(burned_area > as.numeric(result_fire_historical[result_fire_historical$scen == "Historical", "mean_q90"])) %>% 
  group_by(scen, gcm, rep, period) %>% 
  summarize(n_tot = mean(n_tot),
            n = n(),
            n_share = n/n_tot) %>% 
  group_by(scen, period) %>% 
  summarise(n_tot = sum(n),
            mean_n = mean(n),
            mean_n_share = mean(n_share),
            sd_n_share = sd(n_share)) %>% 
  mutate(mean_years = 1/mean_n_share,
         sd_years = mean_years*sd_n_share,
         se_years = sd_years / sqrt(n_tot),
         ci_year = se_years * 1.96)
summary_table_fire




### for wind ---
wind_dat_plot_hist <- read_csv(paste0(path, "/svd_simulations/results_eu/wind_results_final.csv"))
wind_dat_plot_hist <- wind_dat_plot_hist %>% 
  filter(!(sim == 0 & year == 35)) # remove 2020 from historical data as there was no splitting up to agents
wind_dat_plot_hist <- wind_dat_plot_hist %>% 
  mutate(scen = ifelse(sim == 0, "Historical", scen),
         sim = ifelse(sim == 0, 99, sim)) %>% 
  distinct() %>% 
  rowwise() %>% 
  mutate(gcm = master_tab[master_tab$ID == sim, "GCM"],
         rep = master_tab[master_tab$ID == sim, "Rep"])


# calculate extremes
result_wind <- wind_dat_plot_hist %>%
  group_by(scen, gcm, rep) %>%
  summarize(mean_area = mean(dist_area , na.rm = TRUE),
            q90 = quantile(dist_area , 0.90, na.rm = TRUE),
            q99 = quantile(dist_area , 0.99, na.rm = TRUE)) %>%
  group_by(scen) %>%
  summarize(mean_area = mean(mean_area),
            sd_mean = sd(mean_area),
            mean_q90 = mean(q90),
            sd_q90 = sd(q90),
            mean_q99 = mean(q99),
            sd_q99 = sd(q99)) %>%
  mutate(percentual_difference = if_else(scen == "Historical", 0, ((mean_q90 - first(mean_q90)) / first(mean_q90)) * 100))


result_wind


summary_table_wind <- wind_dat_plot_hist %>%
  group_by(scen, gcm, rep) %>% 
  mutate(n_tot = n()) %>% 
  filter(dist_area > as.numeric(result_wind[result_wind$scen == "Historical", "mean_q90"])) %>% 
  group_by(scen, gcm, rep) %>% 
  summarize(n_tot = mean(n_tot),
            n = n(),
            n_share = n/n_tot) %>% 
  group_by(scen) %>% 
  summarise(n_tot = sum(n),
            mean_n = mean(n),
            mean_n_share = mean(n_share),
            sd_n_share = sd(n_share)) %>% 
  mutate(mean_years = 1/mean_n_share,
         sd_years = mean_years*sd_n_share,
         se_years = sd_years / sqrt(n_tot))


summary_table_wind <- wind_dat_plot_hist %>%
  group_by(scen, gcm, rep, period) %>%
  mutate(n_tot = n()) %>% 
  filter(dist_area > as.numeric(result_wind[result_wind$scen == "Historical", "mean_q90"])) %>% 
  group_by(scen, gcm, rep, period) %>% 
  summarize(n_tot = mean(n_tot),
            n = n(),
            n_share = n/n_tot) %>% 
  group_by(scen, period) %>% 
  summarise(n_tot = sum(n),
            mean_n = mean(n),
            mean_n_share = mean(n_share),
            sd_n_share = sd(n_share)) %>% 
  mutate(mean_years = 1/mean_n_share,
         sd_years = mean_years*sd_n_share,
         se_years = sd_years / sqrt(n_tot))
summary_table_wind



### for BB ---
bbtl_dat_plot_hist <- read_csv(paste0(path, "/svd_simulations/results_eu/bbtl_results_final.csv"))
bbtl_dat_plot_hist <- bbtl_dat_plot_hist %>% 
  filter(!(sim == 0 & year == 35))# remove 2020 from historical data as there was no splitting up to agents

bbtl_dat_plot_hist <- bbtl_dat_plot_hist %>% 
  mutate(scen = ifelse(sim == 0, "Historical", scen),
         sim = ifelse(sim == 0, 99, sim)) %>% 
  distinct() %>% 
  rowwise() %>% 
  mutate(gcm = master_tab[master_tab$ID == sim, "GCM"],
         rep = master_tab[master_tab$ID == sim, "Rep"]) 



result_bbtl <- bbtl_dat_plot_hist %>%
  filter(scen == "Historical") %>% 
  mutate(period = "1986-2020") %>% 
  group_by(scen, gcm, rep, period) %>%
  summarize(mean_area = mean(dist_area , na.rm = TRUE),
            q95 = quantile(dist_area , 0.90, na.rm = TRUE),
            q99 = quantile(dist_area , 0.99, na.rm = TRUE)) %>% 
  group_by(scen) %>% 
  summarize(mean_area = mean(mean_area),
            sd_mean = sd(mean_area),
            mean_q95 = mean(q95),
            sd_q95 = sd(q95),
            mean_q99 = mean(q99),
            sd_q99 = sd(q99))



result_bbtl %>% 
  mutate(percentual_difference = if_else(scen == "Historical", 0,
                                         ((mean_q95 - as.numeric(result_bbtl[result_bbtl$scen == "Historical", "mean_q95"])) /
                                            as.numeric(result_bbtl[result_bbtl$scen == "Historical", "mean_q95"])) * 100))


bbtl_dat_plot_hist %>%
  group_by(scen, gcm, rep) %>% 
  mutate(n_tot = n()) %>% 
  filter(dist_area > as.numeric(result_bbtl[result_bbtl$scen == "Historical", "mean_q95"])) %>% 
  group_by(scen, gcm, rep) %>% 
  summarize(n_tot = mean(n_tot),
            n = n(),
            n_share = n/n_tot) %>% 
  group_by(scen) %>% 
  summarise(mean_n = mean(n),
         mean_n_share = mean(n_share),
         sd_n_share = sd(n_share)) %>% 
  mutate(mean_years = 1/mean_n_share,
         sd_years = mean_years*sd_n_share)


summary_table_bbtl <- bbtl_dat_plot_hist %>%
  group_by(scen, gcm, rep, period) %>% 
  mutate(n_tot = n()) %>% 
  filter(dist_area > as.numeric(result_bbtl[result_bbtl$scen == "Historical", "mean_q95"])) %>% 
  group_by(scen, gcm, rep, period) %>% 
  summarize(n_tot = mean(n_tot),
            n = n(),
            n_share = n/n_tot) %>% 
  group_by(scen, period) %>% 
  summarise(n_tot = sum(n),
            mean_n = mean(n),
            mean_n_share = mean(n_share),
            sd_n_share = sd(n_share)) %>% 
  mutate(mean_years = 1/mean_n_share,
         sd_years = mean_years*sd_n_share,
         se_years = sd_years / sqrt(n_tot),
         ci_year = se_years * 1.96)

summary_table_bbtl


summary_table_fire <- bbtl_dat_plot_hist %>%
  filter(!period %in% c("1986-2020")) %>% 
  group_by(scen, gcm, rep) %>%
  mutate(n_tot = n()) %>% 
  filter(dist_area > as.numeric(result_bbtl[result_bbtl$scen == "Historical", "mean_q95"])) %>% 
  group_by(scen, gcm, rep) %>% 
  summarize(n_tot = mean(n_tot),
            n = n(),
            n_share = n/n_tot) %>% 
  group_by(scen) %>% 
  summarise(n_tot = sum(n),
            mean_n = mean(n),
            mean_n_share = mean(n_share),
            sd_n_share = sd(n_share)) %>% 
  mutate(mean_years = 1/mean_n_share,
         sd_years = mean_years*sd_n_share,
         se_years = sd_years / sqrt(n_tot),
         ci_year = se_years * 1.96)
summary_table_fire



### end ------------------------------------
