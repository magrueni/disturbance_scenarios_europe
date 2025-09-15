# --------------------------------------------------------------------------------------
# Script Name: Fig3bcd_extreme_events
# Description: This script plots the data from the script 04_annual_outputs_processing.R
#              as extreme events of the three disturbance agents
#
# Author: Marc Grünig
# Date: 7.2.2025
#
# ------------------------------------------------------------------------------


### libraries ------------------------------------------------------------------
library(dplyr)
library(readr)
library(ggplot2)
library(scales)


path <- "/.../"

master_tab <- read.csv(paste0(path, "/09_svd_simulations/svd_simulations_ids.csv"), sep = ";")



### for fire ---

# load data and get the scenario information
fire_dat_plot_hist <- read_csv(paste0(path, "/10_results/fire_results_final.csv"))
fire_dat_plot_hist <- fire_dat_plot_hist %>% 
  mutate(scen = ifelse(sim == "0", "Historical", scen),
         sim = ifelse(sim == "0", "99", sim)) %>% 
  distinct() %>% 
  rowwise() %>% 
  mutate(gcm = master_tab[as.numeric(sim), "GCM"],
         rep = master_tab[as.numeric(sim), "Rep"]) %>% 
  group_by(year, scen) %>% 
  summarize(burned_area = mean(burned_area, na.rm = T)) 


write_csv(fire_dat_plot_hist, paste0(path, "/11_figures/figure_data/fig3b_fire_extremes.csv"))

# get historical data and calculate range from Q10 to Q90
hist_mean <- fire_dat_plot_hist %>% 
  filter(scen == "Historical") %>% 
  group_by(scen) %>% 
  summarise(hist_mean = mean(burned_area)) %>% 
  dplyr::select(hist_mean) %>% 
  as.numeric()

hist_q10 <- fire_dat_plot_hist %>% 
  filter(scen == "Historical") %>% 
  group_by(scen) %>% 
  summarise(hist_q90 = quantile(burned_area, 0.10)) %>% 
  dplyr::select(hist_q90) %>% 
  as.numeric()


hist_q90 <- fire_dat_plot_hist %>% 
  filter(scen == "Historical") %>% 
  group_by(scen) %>% 
  summarise(hist_q90 = quantile(burned_area, 0.90)) %>% 
  dplyr::select(hist_q90) %>% 
  as.numeric()



# Define the colors for each scenario and the vertical lines
cols <- c("Historical" = "lightgrey", "RCP2.6" = "#ed9b49", "RCP4.5" = "#8d1c06", "RCP8.5" = "#3c0d03",
          "Historical mean" = "blue", "Historical q90" = "red")

# Create the density plot for burned area grouped by sim, colored and filled by scen
p1_fire <- ggplot(fire_dat_plot_hist %>% filter(scen != "Historical")) +
  
  labs(x =bquote("Annual area disturbed [ha yr"^-1*"]"), y = "Density") +
  scale_color_manual(values = cols, name = "Scenario") +
  scale_fill_manual(values = cols, name = "Scenario") +
  theme_classic() +
  geom_rect(aes(xmin = round(hist_q10, 0), xmax = round(hist_q90, 0), ymin = 0, ymax = Inf),
            fill = "lightgrey", alpha = 0.1) +  # Add a grey rectangle between historical mean and q90
  geom_density(aes(x = burned_area, fill = scen, color = scen, group = as.factor(scen))
               , alpha = 0.3, linewidth = 0.5) +
  scale_x_continuous(breaks = c(0, round(hist_q10, -4), round(hist_q90, -3), round(as.numeric(quantile(fire_dat_plot_hist$burned_area, 1, na.rm = T)), -4)),
                     labels = scales::comma,
                     limits = c(0, max(fire_dat_plot_hist$burned_area))) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(alpha = 1)),  
    color = guide_legend(order = 1, override.aes = list(alpha = 1))  
  ) +
  theme(legend.position = "none",
        axis.text = element_text(size = 10),
        plot.margin = margin(1, 1, 1, 1, "cm")) 

p1_fire



ggsave(p1_fire, filename = paste0(path, "/11_figures/Fig3b_extremes_fires.png"),
       width = 5, height = 4)



### for wind ---
wind_dat_plot_hist <- read_csv(paste0(path, "/09_svd_simulations/results_eu/wind_results_final.csv"))
wind_dat_plot_hist <- wind_dat_plot_hist %>% 
  filter(!(sim == 0 & year == 35))# remove 2020 from historical data as there was no splitting up to agents
wind_dat_plot_hist <- wind_dat_plot_hist %>% 
  mutate(scen = ifelse(sim == 0, "Historical", scen),
         sim = ifelse(sim == 0, 99, sim)) %>% 
  distinct() %>% 
  rowwise() %>% 
  mutate(gcm = master_tab[master_tab$ID == sim, "GCM"],
         rep = master_tab[master_tab$ID == sim, "Rep"]) %>% 
  group_by(year, scen) %>% 
  summarize(dist_area = mean(dist_area, na.rm = T)) 

write_csv(wind_dat_plot_hist, paste0(path, "/11_figures/figure_data/fig3c_wind_extremes.csv"))


hist_mean <- wind_dat_plot_hist %>% 
  filter(scen == "Historical") %>% 
  group_by(scen) %>% 
  summarise(hist_mean = mean(dist_area, na.rm = T)) %>% 
  dplyr::select(hist_mean) %>% 
  as.numeric()

hist_q10 <- wind_dat_plot_hist %>% 
  filter(scen == "Historical") %>% 
  group_by(scen) %>% 
  summarise(hist_q90 = quantile(dist_area, 0.05, na.rm = T)) %>% 
  dplyr::select(hist_q90) %>% 
  as.numeric()

hist_q90 <- wind_dat_plot_hist %>% 
  filter(scen == "Historical") %>% 
  group_by(scen) %>% 
  summarise(hist_q90 = quantile(dist_area, 0.95, na.rm = T)) %>% 
  dplyr::select(hist_q90) %>% 
  as.numeric()


# Define the colors for each scenario and the vertical lines
cols <- c("Historical" = "darkgrey", "RCP2.6" = "#ed9b49", "RCP4.5" = "#8d1c06", "RCP8.5" = "#3c0d03",
          "Historical mean" = "blue", "Historical q90" = "red")

# Create the density plot for burned area grouped by sim, colored and filled by scen
p1_wind <- ggplot(wind_dat_plot_hist %>% filter(scen != "Historical")) +
  
  labs(x = bquote("Annual area disturbed [ha yr"^-1*"]"), y = "Density") +
  scale_color_manual(values = cols, name = "Scenario") +
  scale_fill_manual(values = cols, name = "Scenario") +
  theme_classic() +
  geom_rect(aes(xmin = round(hist_q10, 0), xmax = round(hist_q90, 0), ymin = 0, ymax = Inf),
            fill = "lightgrey", alpha = 0.1) +  # Add a grey rectangle between historical mean and q90
  geom_density(aes(x = dist_area, fill = scen, color = scen, group = as.factor(scen))
               , alpha = 0.3, linewidth = 0.5) +
  scale_x_continuous(breaks = c(round(hist_q10, -4), round(hist_q90, -3),
                                round(as.numeric(quantile(wind_dat_plot_hist$dist_area, 1, na.rm = T)), -4)),
                     labels = scales::comma,
                     limits = c(0, hist_q90+ 30000)) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(alpha = 1)),  
    color = guide_legend(order = 1, override.aes = list(alpha = 1))  
  ) +
  theme(legend.position = "none",
        axis.text = element_text(size = 10),
        plot.margin = margin(1, 1, 1, 1, "cm")) 


p1_wind



ggsave(p1_wind, filename = paste0(path, "/11_figures/Fig3c_extremes_wind.png"),
       width = 5, height = 4)





### for bbtl ---
bbtl_dat_plot_hist <- read_csv(paste0(path, "/09_svd_simulations/results_eu/bbtl_results_final.csv"))
bbtl_dat_plot_hist <- bbtl_dat_plot_hist %>% 
  filter(!(sim == 0 & year == 35))# remove 2020 from historical data as there was no splitting up to agents
bbtl_dat_plot_hist <- bbtl_dat_plot_hist %>% 
  mutate(scen = ifelse(sim == 0, "Historical", scen),
         sim = ifelse(sim == 0, 99, sim)) %>% 
  distinct() %>% 
  rowwise() %>% 
  mutate(gcm = master_tab[master_tab$ID == sim, "GCM"],
         rep = master_tab[master_tab$ID == sim, "Rep"])%>% 
  group_by(year, scen) %>% 
  summarize(dist_area = mean(dist_area, na.rm = T)) 


write_csv(bbtl_dat_plot_hist, paste0(path, "/11_figures/figure_data/fig3d_bbtl_extremes.csv"))


hist_mean <- bbtl_dat_plot_hist %>% 
  filter(scen == "Historical") %>% 
  group_by(scen) %>% 
  summarise(hist_mean = mean(dist_area, na.rm = T)) %>% 
  dplyr::select(hist_mean) %>% 
  as.numeric()

hist_q10 <- bbtl_dat_plot_hist %>% 
  filter(scen == "Historical") %>% 
  group_by(scen) %>% 
  summarise(hist_q90 = quantile(dist_area, 0.10, na.rm = T)) %>% 
  dplyr::select(hist_q90) %>% 
  as.numeric()

hist_q90 <- bbtl_dat_plot_hist %>% 
  filter(scen == "Historical") %>% 
  group_by(scen) %>% 
  summarise(hist_q90 = quantile(dist_area, 0.90, na.rm = T)) %>% 
  dplyr::select(hist_q90) %>% 
  as.numeric()


# Define the colors for each scenario and the vertical lines
cols <- c("Historical" = "darkgrey", "RCP2.6" = "#ed9b49", "RCP4.5" = "#8d1c06", "RCP8.5" = "#3c0d03",
          "Historical mean" = "blue", "Historical q90" = "red")

# Create the density plot for burned area grouped by sim, colored and filled by scen
p1_bbtl <- ggplot(bbtl_dat_plot_hist %>% filter(scen != "Historical")) +
  
  labs(x = bquote("Annual area disturbed [ha yr"^-1*"]"), y = "Density") +
  scale_color_manual(values = cols, name = "Scenario") +
  scale_fill_manual(values = cols, name = "Scenario") +
  theme_classic() +
  geom_rect(aes(xmin = round(hist_q10, 0), xmax = round(hist_q90, 0), ymin = 0, ymax = Inf),
            fill = "lightgrey", alpha = 0.1, ) +  # Add a grey rectangle between historical mean and q90
  geom_density(aes(x = dist_area, fill = scen, color = scen, group = as.factor(scen))
               , alpha = 0.3, linewidth = 0.5) +
  scale_x_continuous(breaks = c(round(hist_q10, -4), round(hist_q90, -3),
                                round(as.numeric(quantile(bbtl_dat_plot_hist$dist_area, 0.99, na.rm = T)), -4)),
                     labels = scales::comma,
                     limits = c(0, quantile(bbtl_dat_plot_hist$dist_area, 0.991,  na.rm = T))) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(alpha = 1)),  
    color = guide_legend(order = 1, override.aes = list(alpha = 1))  
  ) +
  theme(legend.position = "none",
        axis.text = element_text(size = 10),
        plot.margin = margin(0.5, 1, 0.5, 1, "cm"))  

p1_bbtl

ggsave(p1_bbtl, filename = paste0(path, "/11_figures/Fig3d_extremes_bbtl.png"),
       width = 5, height = 4)


## legend
# Define the colors for each scenario and the vertical lines
cols <- c("Historical" = "lightgrey", "RCP2.6" = "#ed9b49", "RCP4.5" = "#8d1c06", "RCP8.5" = "#3c0d03",
          "Historical mean" = "blue", "Historical 90th percentile" = "red")

# Create the density plot for burned area grouped by sim, colored and filled by scen
legend <- ggplot(bbtl_dat_plot_hist) +
  geom_density(aes(x = (dist_area), fill = scen, color = scen, group = as.factor(scen), y = ..count..),
               alpha = 0.3, linewidth = 0.8) +
  labs(x = "Annual disturbed area [ha]", y = "Density") +
  scale_color_manual(values = cols, name = "Scenario") +
  scale_fill_manual(values = cols, name = "Scenario") +
  theme_classic() +
  geom_vline(aes(xintercept = (hist_mean), linetype = "Historical mean", color = "Historical mean"), size = 0.8) +
  geom_vline(aes(xintercept = (hist_q90), linetype = "Historical 90th percentile", color = "Historical 90th percentile"), size = 0.8) +
  scale_linetype_manual(name = "", values = c("Historical mean" = "dashed", "Historical 90th percentile" = "dashed")) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(alpha = 1)),
    linetype = guide_legend(order = 2, override.aes = list(color = c("red", "blue")))
  )


legend

legend_extremes <- cowplot::get_legend(legend)

# ggsave(legend_extremes, filename = paste0(path, "/09_svd_simulations/results_eu/11_figures/extremes_legend_heat.png"), width = 7.5, height = 7.5)


### end ------------------------------------------------------------------------

