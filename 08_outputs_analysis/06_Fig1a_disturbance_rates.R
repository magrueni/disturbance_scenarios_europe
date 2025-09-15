# ------------------------------------------------------------------------------
# Script Name: fig1_disturbance_rates.R
# Description: This script processes the annual outputs of all agents and creates
# plot. In particular figure 1a, figure S1 and figures S16-S18 are created
#
# Author: Marc Grünig
# Date: 7.2.2025
#
# Input Files:
#   - Annual outputs from script 4
#
# Output:
#   - figure 1a, figure S1, S16 - S18
# ------------------------------------------------------------------------------


### libraries & settings -------------------------------------------------------

library(terra)
library(tidyverse)
library(sf)
library(stars)

# paths
path <- "/.../"


### load fire data -------------------------------------------------------------
df_plot_fire <- read_csv(paste0(path, "/10_results/fire_results_fut_final.csv"))
df_plot_fire <- df_plot_fire %>% 
  group_by(Year, sim, scen) %>% 
  summarize(burned_area = sum(realized_size, na.rm = T))

# Calculate the mean of n_impact for each year
mean_by_year <- df_plot_fire %>%
  group_by(Year, scen) %>%
  summarise(burned_area_mean = mean(burned_area))


# combine with historical data
hist_fire <- read_csv(paste0(path, "/10_results/fire_results_hist_final.csv"))
hist_fire_plot <- hist_fire %>% 
  mutate(scen = "hist",
         sim = 0,
         Year = year + 1985) %>% 
  dplyr::select(Year, scen, burned_area_mean = burned_area)

df_plot_all_fire <- rbind(hist_fire_plot, mean_by_year %>% 
                            dplyr::select(Year, scen, burned_area_mean = burned_area_mean))



write_csv(df_plot_fire, paste0(path_public, "/11_figures/figure_data/figure_data/FigS16_fut.csv"))
write_csv(df_plot_all_fire, paste0(path_public, "/11_figures/figure_data/figure_data/FigS16_hist.csv"))


# plot 
p1_fire <- ggplot() +
  
  geom_line(data = df_plot_fire, aes(y = burned_area, x = Year, group = sim, color = scen), alpha = 0.1, linewidth = 0.5) +
  geom_line(data = df_plot_all_fire, aes(y = burned_area_mean, x = Year, color = scen), linewidth = 1.5) +
  scale_color_manual(values = c("RCP2.6" = "#ed9b49",
                                "RCP4.5" = "#8d1c06",
                                "RCP8.5" = "#3c0d03",
                                "Historical" = "grey"),
                     name = "Scenario",
                     breaks = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical"),
                     labels = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical")) +
  scale_fill_manual(values = c("RCP2.6" = "#ed9b49",
                               "RCP4.5" = "#8d1c06",
                               "RCP8.5" = "#3c0d03",
                               "Historical" = "grey"),
                    name = "Scenario",
                    breaks = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical"),
                    labels = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical")) +
  scale_x_continuous(breaks = c(1986, 2001, 2020, 2050, 2100)) +
  theme_classic() +
  labs(title = "",
       x = "Year",
       y = "disturbed area [ha]")+
  guides(col="none") +
  theme(text = element_text(size = 16),
        axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 1))

p1_fire
ggsave(p1_fire, filename = paste0(path, "/11_figures/FigS16.png"), width = 10, height = 8)

# 
# # Using the cowplot package
# a_fire_plot_egend <- ggplot() +
#   geom_line(data = df_plot_all_fire, aes(y = burned_area_mean, x = Year, color = scen), linewidth = 0.5) +
#   geom_line(data = df_plot_fire, aes(y = burned_area, x = Year, group = run, color = scen), alpha = 0.3, linewidth = 0.5) +
#   scale_color_manual(values = c("RCP2.6" = "#ed9b49",
#                                 "RCP4.5" = "#8d1c06",
#                                 "RCP8.5" = "#3c0d03",
#                                 "Historical" = "grey"),
#                      name = "Scenario",
#                      breaks = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical"),
#                      labels = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical")) +
#   scale_fill_manual(values = c("RCP2.6" = "#ed9b49",
#                                "RCP4.5" = "#8d1c06",
#                                "RCP8.5" = "#3c0d03",
#                                "Historical" = "grey"),
#                     name = "Scenario",
#                     breaks = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical"),
#                     labels = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical")) +
#   scale_y_continuous(breaks = seq(0, 1000000, 300000)) +
#   theme_classic() +
#   labs(#title = "Barkbeetle area",
#     x = "Years",
#     y = "disturbed area [ha]") +
#   theme(text = element_text(size = 12),
#         axis.title.x = element_text(vjust = -1),
#         axis.title.y = element_text(vjust = 1))
# legend <- cowplot::get_legend(a_fire_plot_egend)
# 
# library(grid)
# grid.newpage()
# grid.draw(legend)




### wind data ------------------------------------------------------------------
df_plot_wind <- read_csv(paste0(path, "/10_results/wind_results_fut_final.csv"))
df_plot_wind <- df_plot_wind %>% 
  group_by(Year, sim, scen) %>% 
  summarize(wind_area = sum(cells_affected  , na.rm = T))

# Calculate the mean of n_impact for each year
mean_by_year <- df_plot_wind %>%
  group_by(Year, scen) %>%
  summarise(wind_area_mean = mean(wind_area))


# combine with historical data
hist_wind <- read_csv(paste0(path, "/10_results/wind_results_hist_final.csv"))
hist_wind_plot <- hist_wind %>% 
  mutate(scen = "hist",
         sim = 0,
         Year = year + 1985) %>% 
  dplyr::select(Year, scen, wind_area_mean = dist_area )

df_plot_all_wind <- rbind(hist_wind_plot, mean_by_year %>% 
                            dplyr::select(Year, scen, wind_area_mean = wind_area_mean))%>% 
  filter(Year != 2020) # exclude because we don't have the proportion between wind and bbtl


write_csv(df_plot_wind, paste0(path_public, "/11_figures/figure_data/FigS18_fut.csv"))
write_csv(df_plot_all_wind, paste0(path_public, "/11_figures/figure_data/figure_data/FigS18_hist.csv"))


# Plot the data and mean lines
p1_wind <- ggplot() +
  geom_line(data = df_plot_wind, aes(y = wind_area, x = Year, group = sim, color = scen), alpha = 0.1, linewidth = 0.5) +
  geom_line(data = df_plot_all_wind, aes(y = wind_area_mean, x = Year, color = scen), linewidth = 1.5) +
  scale_color_manual(values = c("RCP2.6" = "#ed9b49",
                                "RCP4.5" = "#8d1c06",
                                "RCP8.5" = "#3c0d03",
                                "Historical" = "grey"),
                     name = "Scenario",
                     breaks = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical"),
                     labels = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical")) +
  scale_fill_manual(values = c("RCP2.6" = "#ed9b49",
                               "RCP4.5" = "#8d1c06",
                               "RCP8.5" = "#3c0d03",
                               "Historical" = "grey"),
                    name = "Scenario",
                    breaks = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical"),
                    labels = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical")) +
  theme_classic() +
  labs(title = "",
       x = "Years",
       y = "disturbed area [ha]") +
  guides(col="none") +
  theme(text = element_text(size = 16))  # Adjust the size as needed
p1_wind
ggsave(p1_wind, filename = paste0(path, "/11_figures/FigS18.png"), width = 10, height = 8)





### bark beetle ----------------------------------------------------------------
df_plot_bbtl <- read_csv(paste0(path, "/10_results/bbtl_results_fut_final.csv"))
df_plot_bbtl <- df_plot_bbtl %>% 
  group_by(Year, sim, scen) %>% 
  summarize(bbtl_area = sum(n_impact  , na.rm = T))

# Calculate the mean of n_impact for each year
mean_by_year <- df_plot_bbtl %>%
  group_by(Year, scen) %>%
  summarise(bbtl_area_mean = mean(bbtl_area))


# combine with historical data
hist_bbtl <- read_csv(paste0(path, "/10_results/bbtl_results_hist_final.csv"))
hist_bbtl_plot <- hist_bbtl %>% 
  mutate(scen = "hist",
         sim = 0,
         Year = year + 1985) %>% 
  dplyr::select(Year, scen, bbtl_area_mean = dist_area )

df_plot_all_bbtl <- rbind(hist_bbtl_plot, mean_by_year %>% 
                            dplyr::select(Year, scen, bbtl_area_mean = bbtl_area_mean)) %>% 
  filter(Year != 2020) # exclude because we don't have the proportion between wind and bbtl


write_csv(df_plot_wind, paste0(path_public, "/11_figures/figure_data/FigS17_fut.csv"))
write_csv(df_plot_all_wind, paste0(path_public, "/11_figures/figure_data/figure_data/FigS17_hist.csv"))


# plot 
p1_bbtl <- ggplot() +
  
  geom_line(data = df_plot_bbtl, aes(y = bbtl_area, x = Year, group = sim, color = scen), alpha = 0.1, linewidth = 1) +
  geom_line(data = df_plot_all_bbtl, aes(y = bbtl_area_mean, x = Year, color = scen), linewidth = 1.5) +
  scale_color_manual(values = c("RCP2.6" = "#ed9b49",
                                "RCP4.5" = "#8d1c06",
                                "RCP8.5" = "#3c0d03",
                                "Historical" = "grey"),
                     name = "Scenario",
                     breaks = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical"),
                     labels = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical")) +
  scale_fill_manual(values = c("RCP2.6" = "#ed9b49",
                               "RCP4.5" = "#8d1c06",
                               "RCP8.5" = "#3c0d03",
                               "Historical" = "grey"),
                    name = "Scenario",
                    breaks = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical"),
                    labels = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical")) +
  theme_classic() +
  labs(title = "",
       x = "Year",
       y = "disturbed area [ha]")+
  guides(col="none") +
  theme(text = element_text(size = 16),
        axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 1))

p1_bbtl
ggsave(p1_bbtl, filename = paste0(path, "/11_figures/FigS17.png"), width = 10, height = 8)





### combine everything ---------------------------------------------------------

df_plot_all_fut <- rbind(df_plot_fire %>% 
                           rename("dist_area" = "burned_area") %>% 
                           mutate(agent = "Fire"),
                         df_plot_wind %>% 
                           mutate(sim =as.numeric(sim)) %>% 
                           rename("dist_area" = "wind_area") %>% 
                           mutate(agent = "Wind"), 
                         df_plot_bbtl %>% 
                           rename("dist_area" = "bbtl_area") %>% 
                           mutate(agent = "Bark beetle"))


df_plot_all_hist <- rbind(df_plot_all_fire %>%
                            filter(Year <= 2020) %>% 
                            rename("dist_area" = "burned_area_mean") %>%
                            mutate(agent = "Fire") %>%
                            mutate(scen = "Historical"),
                          df_plot_all_wind %>%
                            filter(Year <= 2020) %>% 
                            rename("dist_area" = "wind_area_mean") %>%
                            mutate(agent = "Wind") %>%
                            mutate(scen = "Historical"),
                          df_plot_all_bbtl %>%
                            filter(Year <= 2020) %>% 
                            rename("dist_area" = "bbtl_area_mean") %>%
                            mutate(agent = "Bark beetle")) %>%
  mutate(sim = as.numeric(0))

# combine historical and future data
df_combined <- bind_rows(df_plot_all_hist, df_plot_all_fut)

df_combined_tot <- df_combined %>%
  group_by(Year, scen, sim) %>% 
  summarize(dist_area = sum(dist_area)) %>% 
  mutate(agent = "All")

df_combined <- bind_rows(df_combined_tot) %>% 
  mutate(agent = ifelse(agent == "Barkbeete", "Bark beetle", agent)) %>% 
  mutate(gcm = ifelse(sim %in% c(1, 2, 3, 10, 11, 12, 19, 20, 21, 28, 29, 30, 37, 38, 39), "ichec", 
                      ifelse(sim %in% c(4, 5, 6, 13, 14, 15, 22, 23, 24, 31, 32, 33, 40, 41, 42), "mpi", "ncc"))) %>% 
  mutate(rep = ceiling(as.numeric(sim)/9))

hist_dat <- df_combined %>%
  filter(scen == "Historical") 

fut_dat <- df_combined %>% 
  filter(scen != "Historical")


write_csv(fut_dat, paste0(path, "/10_results/fut_dist_annual_final.csv"))
write_csv(hist_dat, paste0(path, "/10_results/hist_dist_annual_final.csv"))


### figure 1a -------------------------------------------------------------------

# calculate the forest area
# init_veg <- rast(paste0(path, "/initial_states/init_veg_v2.tif"))
# init_veg[!is.na(init_veg)] <- 1
# total_forest_area <- sum(values(init_veg), na.rm = T)

# precalculated - difference to historical because remote sensing data is 30m resolution
total_forest_area <- 186739717
historical_area <- 207334385


# # load data
# fut_dat <- read_csv(paste0(path, "/10_results/fut_dist_annual_final.csv"))
# hist_dat <- read_csv(paste0(path, "/10_results/hist_dist_annual_final.csv"))


df_combined_sub <- rbind(hist_dat %>% filter(Year <= 2020), fut_dat %>% filter(Year > 2020))
df_combined_sub <- df_combined_sub %>% filter(agent == "All")
df_combined_sub_hist <- df_combined_sub %>%
  filter(sim %in% c(0)) %>% 
  filter(scen == "Historical") %>% 
  mutate(dist_rate = 100/historical_area*dist_area,
         dist_freq = historical_area/dist_area)

df_plot <- fut_dat %>% filter(Year > 2020) %>% 
  mutate(dist_rate = 100/total_forest_area*dist_area,
         dist_freq = total_forest_area/dist_area) %>% 
  group_by(Year, scen) %>% 
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


# change historical values to a box --
df_boxplot <- df_combined_sub_hist %>%
  filter(scen == "Historical") %>% 
  mutate(Year = ifelse(Year >= 1986 & Year <= 2000, 2000,
                       ifelse(Year >= 2001 & Year <= 2020, 2010, NA))) %>% 
  mutate(scen_new = ifelse(Year == 2000, "hist1",
                           ifelse(Year == 2010, "hist2", "hist3"))) %>% 
  drop_na()




write_csv(df_boxplot, paste0(path, "/11_figures/figure_data/Fig1a_FigS1_hist.csv"))
write_csv(df_plot, paste0(path, "/11_figures/figure_data/Fig1a_FigS1_fut.csv"))



# Define the custom labels
custom_labels <- c("1990", "2010" , seq(2020, max(df_plot$Year), by = 10))
all_plot <- ggplot() +
  # Add the boxplot at 1990 and 2010 with the same fill color for both
  geom_boxplot(data = df_boxplot, aes(x = Year, y = dist_rate, fill = "Historical", color = scen_new), 
               width = 5, alpha = 0.5, outlier.shape = NA) +
  # Add the lines and ribbon for future scenarios starting from 2020
  geom_line(data = df_plot, aes(y = mean_dist_rate, x = Year, color = scen), alpha = 1, linewidth = 0.5) +
  geom_ribbon(data = df_plot, aes(ymin = ci_low, ymax = ci_high, x = Year, fill = scen), alpha = 0.2) +
  scale_color_manual(values = c("RCP2.6" = "#ed9b49",
                                "RCP4.5" = "#8d1c06",
                                "RCP8.5" = "#3c0d03"),
                     name = "Scenario",
                     breaks = c("RCP8.5", "RCP4.5", "RCP2.6"),
                     labels = c("RCP8.5", "RCP4.5", "RCP2.6"),
                     guide = guide_legend(override.aes = list(fill = NA))) +
  scale_fill_manual(values = c("RCP2.6" = "#ed9b49",
                               "RCP4.5" = "#8d1c06",
                               "RCP8.5" = "#3c0d03",
                               "Historical" = "grey"),
                    name = "Scenario",
                    breaks = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical"),
                    labels = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical"),
                    guide = guide_legend(override.aes = list(linetype = 1))) +
  theme_classic() +
  labs(x = "Year",
       y = bquote("Disturbance rate [% yr"^-1*"]")) +
  theme(text = element_text(size = 12),
        axis.title.x = element_text(vjust = -2),
        axis.title.y = element_text(vjust = 2),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_x_continuous(breaks = c(2000, 2010, seq(2020, max(df_plot$Year), by = 10)),
                     labels = custom_labels, limits = c(1995, max(df_plot$Year)))

p_fin <- all_plot + theme(legend.position = "none")
p_fin
ggsave(p_fin, filename = paste0(path, "/11_figures/Fig1a.png"), width = 6, height = 5)




### figure S1 ------------------------------------------------------------------

# Define the custom labels
custom_labels <- c("1990", "2010" , seq(2020, max(df_plot$Year), by = 10))
all_plot <- ggplot() +
  # Add the boxplot at 1990 and 2010 with the same fill color for both
  geom_boxplot(data = df_boxplot, aes(x = Year, y = dist_freq, fill = "Historical", color = scen_new), 
               width = 5, alpha = 0.5, outlier.shape = NA) +
  # Add the lines and ribbon for future scenarios starting from 2020
  geom_line(data = subset(df_plot, Year >= 2020), aes(y = mean_dist_freq, x = Year, color = scen), alpha = 1, linewidth = 0.5) +
  geom_ribbon(data = subset(df_plot, Year >= 2020), aes(ymin = ci_low_freq, ymax = ci_high_freq, x = Year, fill = scen), alpha = 0.2) +
  scale_color_manual(values = c("RCP2.6" = "#ed9b49",
                                "RCP4.5" = "#8d1c06",
                                "RCP8.5" = "#3c0d03"),
                     name = "Scenario",
                     breaks = c("RCP8.5", "RCP4.5", "RCP2.6"),
                     labels = c("RCP8.5", "RCP4.5", "RCP2.6"),
                     guide = guide_legend(override.aes = list(fill = NA))) +
  scale_fill_manual(values = c("RCP2.6" = "#ed9b49",
                               "RCP4.5" = "#8d1c06",
                               "RCP8.5" = "#3c0d03",
                               "Historical" = "grey"),
                    name = "Scenario",
                    breaks = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical"),
                    labels = c("RCP8.5", "RCP4.5", "RCP2.6", "Historical"),
                    guide = guide_legend(override.aes = list(linetype = 1))) +
  theme_classic() +
  labs(x = "Year",
       y = bquote("Rotation period [year]")) +
  theme(text = element_text(size = 12),
        axis.title.x = element_text(vjust = -2),
        axis.title.y = element_text(vjust = 2),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_x_continuous(breaks = c(2000, 2010, seq(2020, max(df_plot$Year), by = 10)),
                     labels = custom_labels, limits = c(1995, max(df_plot$Year)))

p_fin <- all_plot + theme(legend.position = "none")
p_fin

ggsave(p_fin, filename = paste0(path, "/11_figures/FigS1.png"), width = 6, height = 5)




### end ------------------------------------------------------------------------
