# ------------------------------------------------------------------------------
# Script Name: Plotting disturbance rates
# Description: This script combines annual disturbance outputs with climate data
# and creates a heatmap for Figure 1b, and heatmaps for the individual agents 
#
# Author: Marc Grünig
# Date: 7.2.2025
#
# Input Files:
#   - Annual outputs from script 4
#
# Output:
#   - Figure 1b
#   - Figure S2
#   - Figure S3
#   - Figure S4
# ------------------------------------------------------------------------------


### libraries ------------------------------------------------------------------
library(terra)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tidyterra)
library(arrow)
library(collapse)
library(sf)
library(cowplot)
library(grid)
library(gridExtra) 
library(patchwork)
library(readr)



### load and prepare data ------------------------------------------------------
path <- "/.../"

gc()
tmpFiles(remove = T)


### load climate data ----------------------------------------------------------
timestep <- "1981-2010"
scenarios_hist <- c("ICHEC-EC-EARTH_historical", "NCC-NorESM1-M_historical", "MPI-M-MPI-ESM-LR_historical")

# get data from the historical scenarios
temp <- list()
for(s in 1:length(scenarios_hist)){
  
  fls <- read_csv(paste0(path, "/05_clim_data/annual_climate/svd_clim_", scenarios_hist[[s]], "_v3_aux_biascorr.csv"))
  temp_df <- fls %>% 
    filter(year >= 2001) %>% 
    rowwise() %>% 
    mutate(summerP = sum(prec_month6*30,
                         prec_month7*31,
                         prec_month8*31),
           ANP = ANP*365) %>%
    dplyr::select(year, MAT, ANP, summerP, summervpd) %>% 
    mutate(scen = strsplit(scenarios_hist[s], "_")[[1]][1])
  temp[[s]] <- temp_df
  
}


# calculate the historical means
temp_hist <- do.call(rbind, temp)

temp_hist <- temp_hist %>% 
  #filter(year >= 2000) %>% 
  group_by(year, scen) %>% 
  summarize(MAT_hist = mean(MAT),
            ANP_hist = mean(ANP),
            summerP_hist = mean(summerP),
            summervpd_hist = mean(summervpd)) %>% 
  mutate(gcm = ifelse(grepl("ICHEC", scen), "ichec",
                      ifelse(grepl("NCC", scen), "ncc", "mpi"))) %>% 
  dplyr::select(-scen)


delta_hist <- temp_hist %>%
  #filter(year >= 1981 & year < 2006) %>%
  group_by(gcm) %>%
  summarize(mean_Mat_hist = mean(MAT_hist),
            mean_ANP_hist = mean(ANP_hist),
            mean_summervpd_hist = mean(summervpd_hist))


# load future data
scenarios_fut <- c("ICHEC-EC-EARTH_rcp_2_6", "NCC-NorESM1-M_rcp_2_6", "MPI-M-MPI-ESM-LR_rcp_2_6",
                   "ICHEC-EC-EARTH_rcp_4_5", "NCC-NorESM1-M_rcp_4_5", "MPI-M-MPI-ESM-LR_rcp_4_5", 
                   "ICHEC-EC-EARTH_rcp_8_5", "NCC-NorESM1-M_rcp_8_5", "MPI-M-MPI-ESM-LR_rcp_8_5")


temp <- list()
for(s in 1:length(scenarios_fut)){
  
  fls <- read_csv(paste0(path, "/05_clim_data/annual_climate/svd_clim_", scenarios_fut[[s]], "_v3_aux_biascorr.csv"))
  temp_df <- fls %>%
    rowwise() %>% 
    mutate(summerP = sum(prec_month6*30,
                         prec_month7*31,
                         prec_month8*31),
           ANP = ANP*365) %>%
    dplyr::select(year, climateId, MAT, ANP, summerP, summervpd) %>% 
    group_by(year) %>% 
    summarize(MAT = mean(MAT),
              ANP = mean(ANP),
              summerP = mean(summerP),
              summervpd = mean(summervpd)) %>% 
    mutate(scen = paste(scenarios_fut[[s]]))
  temp[[s]] <- temp_df
  
}

temp_fut <- do.call(rbind, temp)


# add the 2011 to 2020 period to the historical data
temp_hist2 <- temp_fut %>% 
  filter(year < 2021) %>% 
  mutate(gcm = ifelse(grepl("ICHEC", scen), "ichec",
                      ifelse(grepl("NCC", scen), "ncc", "mpi"))) %>% 
  group_by(year, gcm) %>% 
  summarize(MAT_hist = mean(MAT),
            ANP_hist = mean(ANP), 
            summerP_hist = mean(summerP),
            summervpd_hist = mean(summervpd))

temp_hist <- rbind(temp_hist, temp_hist2) %>% 
  group_by(gcm) %>% 
  summarize(MAT_hist = mean(MAT_hist),
            ANP_hist = mean(ANP_hist),
            summerP_hist = mean(summerP_hist),
            summervpd_hist = mean(summervpd_hist))


# calculate the difference in climate for each year compared to baseline
temp_fut_df <- temp_fut %>% 
  filter(year >= 2021) %>% 
  mutate(gcm = ifelse(grepl("ICHEC", scen), "ichec",
                      ifelse(grepl("NCC", scen), "ncc", "mpi"))) %>% 
  mutate(scen = ifelse(grepl("rcp_2_6", scen), "RCP2.6",
                       ifelse(grepl("rcp_4_5", scen), "RCP4.5", "RCP8.5"))) %>% 
  left_join(., temp_hist, by = "gcm") %>% 
  mutate(MAT_diff = MAT - MAT_hist,
         ANP_diff = ANP - ANP_hist,
         summerP_diff = summerP - summerP_hist,
         summervpd_diff = summervpd - summervpd_hist) 


temp_fut_df %>% 
  mutate(period = case_when(
    year >= 2021 & year <= 2040 ~ "2021-2040",
    year >= 2041 & year <= 2060 ~ "2041-2060",
    year >= 2061 & year <= 2080 ~ "2061-2080",
    year >= 2081 & year <= 2100 ~ "2081-2100"
  )) %>% 
  group_by(scen, period) %>% 
  summarise(delta_MAT = mean(MAT_diff),
            delta_VPD = mean(summervpd_diff))

# calculate delta T and delta P
# delta_fut <- temp_fut_df %>% 
#   filter(year >= 2076) %>% 
#   group_by(scen, gcm) %>% 
#   summarize(mean_Mat = mean(MAT),
#             mean_ANP = mean(ANP))
# 
# delta_clim <- delta_fut %>% 
#   left_join(., delta_hist, by = "gcm") %>% 
#   mutate(delta_T = mean_Mat - mean_Mat_hist,
#          delta_P = mean_ANP - mean_ANP_hist)




### load the disturbed area ---------------------------------------------------

# precalculated - difference to historical because remote sensing data is 30m resolution
total_forest_area <- 186739717
historical_area <- 207334385


# load the collected disturbance data
fut_dat <- read_csv(paste0(path, "/10_results/fut_dist_annual_final.csv"))
df_plot <- fut_dat %>% filter(Year > 2020) %>% 
  mutate(dist_rate = 100/total_forest_area*dist_area,
         dist_freq = total_forest_area/dist_area) 

df_plot_new <- df_plot %>% rename("year" = "Year") %>% 
  left_join(., temp_fut_df, by = c("year", "scen", "gcm")) %>% 
  drop_na()


# Calculate breaks for bins based on integer rounding
mat_diff_breaks <- seq(floor(min(df_plot_new$MAT_diff, na.rm = TRUE)), ceiling(max(df_plot_new$MAT_diff, na.rm = TRUE)), by = 0.5)

# Calculate breaks for ANP_diff (bins with intervals of 50 units)
anp_diff_min <- floor(min(df_plot_new$summervpd_diff, na.rm = TRUE) / 0.05) * 0.05
anp_diff_max <- ceiling(max(df_plot_new$summervpd_diff, na.rm = TRUE) / 0.05) * 0.05
anp_diff_breaks <- seq(anp_diff_min, anp_diff_max, by = 0.05)

# Bin MAT_diff and ANP_diff using case_when
# Assign MAT_diff_bin and ANP_diff_bin directly from breaks
df_plot_fin <- df_plot_new %>%
  mutate(
    MAT_diff_bin = as.factor(mat_diff_breaks[findInterval(MAT_diff, mat_diff_breaks)]),
    ANP_diff_bin = as.factor(anp_diff_breaks[findInterval(summervpd_diff, anp_diff_breaks)])
  )


# Calculate mean dist_rate for each combination of MAT_diff_bin and ANP_diff_bin
df_plot_heatmap <- df_plot_fin %>%
  group_by(MAT_diff_bin, ANP_diff_bin) %>%
  mutate(dist_rate = mean(dist_rate, na.rm = TRUE)) %>%
  ungroup() %>%
  tidyr::unite("combined", MAT_diff_bin, ANP_diff_bin, remove = FALSE) %>%
  mutate(unique_id = as.integer(as.factor(combined))) %>%
  group_by(unique_id) %>% 
  mutate(n = n()) %>% filter(n > 0)


write_csv(df_plot_heatmap, paste0(path, "/11_figures/figure_data/fig1b_heatmap.csv"))

# Create the base heatmap with scatter plot
base_plot <- ggplot() +
  geom_tile(data = df_plot_heatmap, aes(x = MAT_diff_bin, y = ANP_diff_bin, fill = dist_rate)) +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_x_discrete(breaks = c(-2, 0, 2, 4, 6)) +  # Adjust with your desired breaks
  scale_y_discrete(breaks = c(-0.2, 0, 0.2, 0.4, 0.6)) +  # Adjust with your desired breaks
  labs(x = "Δ Temperature [°C]",
       y = "Δ Summer VPD [kPa]",
       fill = "Disturbance rate") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 12))  # Remove legend for heatmap

base_plot


# Define custom colors for scen
custom_colors <- c("RCP2.6" = "#ed9b49", "RCP4.5" = "#8d1c06", "RCP8.5" = "#3c0d03")

# Create marginal density plots
x_density <- ggplot(df_plot_new, aes(x = MAT_diff, fill = scen, color = scen)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_void() +
  theme(legend.position = "none")  # Remove legend for x_density plot

y_density <- ggplot(df_plot_new, aes(x = summervpd_diff, fill = scen, color = scen)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  coord_flip() +
  theme_void() +
  labs(color = "Scenario") +
  theme(legend.position = "none")  # Remove legend for y_density plot

# Combine the plots using patchwork
combined_plot <- x_density + plot_spacer() + base_plot + y_density +
  plot_layout(ncol = 2, widths = c(4, 1), heights = c(1, 4))

# Print the combined plot (without legends)
print(combined_plot)
ggsave(combined_plot, filename = paste0(path, "/11_figures/Fig1b.pdf"), width = 5, height = 5)

# Create legends plot
y_density <- ggplot(df_plot_new, aes(x = ANP_diff, fill = scen, color = scen)) +
  geom_density(alpha = 0.5) +
  # scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors, guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
  coord_flip() +
  theme_void() +
  labs(fill = "Scenario") +
  theme(legend.position = "right")  # Remove legend for y_density plot
legend_rcp <- cowplot::get_legend(y_density)
# ggsave(legend_rcp, filename = paste0(path, "/11_figures/heatmap_legend.png"), width = 5, height = 5)


base_plot <- ggplot() +
  geom_tile(data = df_plot_heatmap, aes(x = MAT_diff_bin, y = ANP_diff_bin, fill = dist_rate)) +
  scale_fill_gradient(low = "blue", high = "red",  breaks = c(0.1, 0.2, 0.3, 0.4), 
                      guide = guide_colourbar(title.position = "top", title.hjust = 0.5)) + #, breaks = c(0.2, 0.4)) +
  labs(x = "Δ Temperature",
       y = "Δ Precipitation",
       fill =  bquote("Disturbance \n rate [% yr"^-1*"]")) +
  theme_classic() +
  theme(legend.position = "bottom",
        text = element_text(size = 12))  # Remove legend for heatmap
legend_heat <- cowplot::get_legend(base_plot)
plot(legend_heat)
# ggsave(legend_heat, filename = paste0(path, "/11_figures/heatmap_legend_heat.png"), width = 5, height = 5)



### numbers ---
df_plot_heatmap %>% 
  group_by(MAT_diff_bin, ANP_diff_bin) %>% 
  summarize(dist_rate = mean(dist_rate)) %>% 
  View()




### create plots for the individual agents ----------------------------------------------------------------------

## first fire

# load the collected disturbance data
fut_dat_fire <- read_csv(paste0(path, "/10_results/fire_results_final.csv"))
df_plot <- fut_dat_fire %>% filter(year > 1) %>% 
  mutate(year = year + 2019) %>% 
  mutate(dist_rate = 100/total_forest_area*burned_area,
         dist_freq = total_forest_area/burned_area) %>% 
  mutate(gcm = ifelse(sim %in% c(1, 2, 3, 10, 11, 12, 19, 20, 21, 28, 29, 30, 37, 38, 39), "ichec", 
                      ifelse(sim %in% c(4, 5, 6, 13, 14, 15, 22, 23, 24, 31, 32, 33, 40, 41, 42), "mpi", "ncc")))

df_plot_new <- df_plot %>% 
  left_join(., temp_fut_df, by = c("year", "scen", "gcm")) %>% 
  drop_na()

# Calculate breaks for bins based on integer rounding
mat_diff_breaks <- seq(floor(min(df_plot_new$MAT_diff, na.rm = TRUE)), ceiling(max(df_plot_new$MAT_diff, na.rm = TRUE)), by = 0.5)

# Calculate breaks for ANP_diff (bins with intervals of 50 units)
anp_diff_min <- floor(min(df_plot_new$summervpd_diff, na.rm = TRUE) / 0.05) * 0.05
anp_diff_max <- ceiling(max(df_plot_new$summervpd_diff, na.rm = TRUE) / 0.05) * 0.05
anp_diff_breaks <- seq(anp_diff_min, anp_diff_max, by = 0.05)

# Bin MAT_diff and ANP_diff using case_when
# Assign MAT_diff_bin and ANP_diff_bin directly from breaks
df_plot_fin <- df_plot_new %>%
  mutate(
    MAT_diff_bin = as.factor(mat_diff_breaks[findInterval(MAT_diff, mat_diff_breaks)]),
    ANP_diff_bin = as.factor(anp_diff_breaks[findInterval(summervpd_diff, anp_diff_breaks)])
  )


# Calculate mean dist_rate for each combination of MAT_diff_bin and ANP_diff_bin
df_plot_heatmap <- df_plot_fin %>%
  group_by(MAT_diff_bin, ANP_diff_bin) %>%
  mutate(dist_rate = mean(dist_rate, na.rm = TRUE)) %>%
  ungroup() %>%
  tidyr::unite("combined", MAT_diff_bin, ANP_diff_bin, remove = FALSE) %>%
  mutate(unique_id = as.integer(as.factor(combined))) %>%
  group_by(unique_id) %>% 
  mutate(n = n()) %>% filter(n > 0)


write_csv(df_plot_heatmap, paste0(path, "/11_figures/figure_data/FigS2.csv"))

# Create the base heatmap with scatter plot
base_plot <- ggplot() +
  geom_tile(data = df_plot_heatmap, aes(x = MAT_diff_bin, y = ANP_diff_bin, fill = dist_rate)) +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_x_discrete(breaks = c(-2, 0, 2, 4, 6)) +  # Adjust with your desired breaks
  scale_y_discrete(breaks = c(-0.2, 0, 0.2, 0.4, 0.6)) +  # Adjust with your desired breaks
  labs(x = "Δ Temperature [°C]",
       y = "Δ Summer VPD [kPa]",
       fill = "Disturbance rate") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 14))  # Remove legend for heatmap

base_plot


# Define custom colors for scen
custom_colors <- c("RCP2.6" = "#ed9b49", "RCP4.5" = "#8d1c06", "RCP8.5" = "#3c0d03")

# Create marginal density plots
x_density <- ggplot(df_plot_new, aes(x = MAT_diff, fill = scen, color = scen)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_void() +
  theme(legend.position = "none")  # Remove legend for x_density plot

y_density <- ggplot(df_plot_new, aes(x = summervpd_diff, fill = scen, color = scen)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  coord_flip() +
  theme_void() +
  labs(color = "Scenario") +
  theme(legend.position = "none")  # Remove legend for y_density plot

# Combine the plots using patchwork
combined_plot <- x_density + plot_spacer() + base_plot + y_density +
  plot_layout(ncol = 2, widths = c(4, 1), heights = c(1, 4))

# Print the combined plot (without legends)
print(combined_plot)

ggsave(combined_plot, filename = paste0(path, "/11_figures/FigS2.png"), width = 5, height = 5)



## now bark beetle

# load the collected disturbance data
fut_dat_fire <- read_csv(paste0(path, "/10_results/bbtl_results_final.csv"))
df_plot <- fut_dat_fire %>% filter(year > 1) %>% 
  mutate(year = year + 2019) %>% 
  mutate(dist_rate = 100/total_forest_area*dist_area ,
         dist_freq = total_forest_area/dist_area ) %>% 
  mutate(gcm = ifelse(sim %in% c(1, 2, 3, 10, 11, 12, 19, 20, 21, 28, 29, 30, 37, 38, 39), "ichec", 
                      ifelse(sim %in% c(4, 5, 6, 13, 14, 15, 22, 23, 24, 31, 32, 33, 40, 41, 42), "mpi", "ncc")))

df_plot_new <- df_plot %>% 
  left_join(., temp_fut_df, by = c("year", "scen", "gcm")) %>% 
  drop_na()


# Calculate breaks for bins based on integer rounding
mat_diff_breaks <- seq(floor(min(df_plot_new$MAT_diff, na.rm = TRUE)), ceiling(max(df_plot_new$MAT_diff, na.rm = TRUE)), by = 0.5)

# Calculate breaks for ANP_diff (bins with intervals of 50 units)
anp_diff_min <- floor(min(df_plot_new$summervpd_diff, na.rm = TRUE) / 0.05) * 0.05
anp_diff_max <- ceiling(max(df_plot_new$summervpd_diff, na.rm = TRUE) / 0.05) * 0.05
anp_diff_breaks <- seq(anp_diff_min, anp_diff_max, by = 0.05)

# Bin MAT_diff and ANP_diff using case_when
# Assign MAT_diff_bin and ANP_diff_bin directly from breaks
df_plot_fin <- df_plot_new %>%
  mutate(
    MAT_diff_bin = as.factor(mat_diff_breaks[findInterval(MAT_diff, mat_diff_breaks)]),
    ANP_diff_bin = as.factor(anp_diff_breaks[findInterval(summervpd_diff, anp_diff_breaks)])
  )


# Calculate mean dist_rate for each combination of MAT_diff_bin and ANP_diff_bin
df_plot_heatmap <- df_plot_fin %>%
  group_by(MAT_diff_bin, ANP_diff_bin) %>%
  mutate(dist_rate = mean(dist_rate, na.rm = TRUE)) %>%
  ungroup() %>%
  tidyr::unite("combined", MAT_diff_bin, ANP_diff_bin, remove = FALSE) %>%
  mutate(unique_id = as.integer(as.factor(combined))) %>%
  group_by(unique_id) %>% 
  mutate(n = n()) %>% filter(n > 0)


write_csv(df_plot_heatmap, paste0(path, "/11_figures/figure_data/FigS3.csv"))

# Create the base heatmap with scatter plot
base_plot <- ggplot() +
  geom_tile(data = df_plot_heatmap, aes(x = MAT_diff_bin, y = ANP_diff_bin, fill = dist_rate)) +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_x_discrete(breaks = c(-2, 0, 2, 4, 6)) +  # Adjust with your desired breaks
  scale_y_discrete(breaks = c(-0.2, 0, 0.2, 0.4, 0.6)) +  # Adjust with your desired breaks
  labs(x = "Δ Temperature [°C]",
       y = "Δ Summer VPD [kPa]",
       fill = "Disturbance rate") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 14))  # Remove legend for heatmap

base_plot

# Define custom colors for scen
custom_colors <- c("RCP2.6" = "#ed9b49", "RCP4.5" = "#8d1c06", "RCP8.5" = "#3c0d03")

# Create marginal density plots
x_density <- ggplot(df_plot_new, aes(x = MAT_diff, fill = scen, color = scen)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_void() +
  theme(legend.position = "none")  # Remove legend for x_density plot

y_density <- ggplot(df_plot_new, aes(x = summervpd_diff, fill = scen, color = scen)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  coord_flip() +
  theme_void() +
  labs(color = "Scenario") +
  theme(legend.position = "none")  # Remove legend for y_density plot

# Combine the plots using patchwork
combined_plot <- x_density + plot_spacer() + base_plot + y_density +
  plot_layout(ncol = 2, widths = c(4, 1), heights = c(1, 4))

# Print the combined plot (without legends)
print(combined_plot)

ggsave(combined_plot, filename = paste0(path, "/11_figures/FigS3.png"), width = 5, height = 5)



## finally wind

# load the collected disturbance data
fut_dat_fire <- read_csv(paste0(path, "/10_results/wind_results_final.csv"))
df_plot <- fut_dat_fire %>% filter(year > 1) %>% 
  mutate(year = year + 2019) %>% 
  mutate(dist_rate = 100/total_forest_area*dist_area ,
         dist_freq = total_forest_area/dist_area ) %>% 
  mutate(gcm = ifelse(sim %in% c(1, 2, 3, 10, 11, 12, 19, 20, 21, 28, 29, 30, 37, 38, 39), "ichec", 
                      ifelse(sim %in% c(4, 5, 6, 13, 14, 15, 22, 23, 24, 31, 32, 33, 40, 41, 42), "mpi", "ncc")))

df_plot_new <- df_plot %>% 
  left_join(., temp_fut_df, by = c("year", "scen", "gcm")) %>% 
  drop_na()


# Calculate breaks for bins based on integer rounding
mat_diff_breaks <- seq(floor(min(df_plot_new$MAT_diff, na.rm = TRUE)), ceiling(max(df_plot_new$MAT_diff, na.rm = TRUE)), by = 0.5)

# Calculate breaks for ANP_diff (bins with intervals of 50 units)
anp_diff_min <- floor(min(df_plot_new$summervpd_diff, na.rm = TRUE) / 0.05) * 0.05
anp_diff_max <- ceiling(max(df_plot_new$summervpd_diff, na.rm = TRUE) / 0.05) * 0.05
anp_diff_breaks <- seq(anp_diff_min, anp_diff_max, by = 0.05)

# Bin MAT_diff and ANP_diff using case_when
# Assign MAT_diff_bin and ANP_diff_bin directly from breaks
df_plot_fin <- df_plot_new %>%
  mutate(
    MAT_diff_bin = as.factor(mat_diff_breaks[findInterval(MAT_diff, mat_diff_breaks)]),
    ANP_diff_bin = as.factor(anp_diff_breaks[findInterval(summervpd_diff, anp_diff_breaks)])
  )


# Calculate mean dist_rate for each combination of MAT_diff_bin and ANP_diff_bin
df_plot_heatmap <- df_plot_fin %>%
  group_by(MAT_diff_bin, ANP_diff_bin) %>%
  mutate(dist_rate = mean(dist_rate, na.rm = TRUE)) %>%
  ungroup() %>%
  tidyr::unite("combined", MAT_diff_bin, ANP_diff_bin, remove = FALSE) %>%
  mutate(unique_id = as.integer(as.factor(combined))) %>%
  group_by(unique_id) %>% 
  mutate(n = n()) %>% filter(n > 0)


write_csv(df_plot_heatmap, paste0(path, "/11_figures/figure_data/FigS4.csv"))

# Create the base heatmap with scatter plot
base_plot <- ggplot() +
  geom_tile(data = df_plot_heatmap, aes(x = MAT_diff_bin, y = ANP_diff_bin, fill = dist_rate)) +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_x_discrete(breaks = c(-2, 0, 2, 4, 6)) +  # Adjust with your desired breaks
  scale_y_discrete(breaks = c(-0.2, 0, 0.2, 0.4, 0.6)) +  # Adjust with your desired breaks
  labs(x = "Δ Temperature [°C]",
       y = "Δ Summer VPD [kPa]",
       fill = "Disturbance rate") +
  theme_classic() +
  theme(legend.position = "none",
        text = element_text(size = 14))  # Remove legend for heatmap

base_plot


# Define custom colors for scen
custom_colors <- c("RCP2.6" = "#ed9b49", "RCP4.5" = "#8d1c06", "RCP8.5" = "#3c0d03")

# Create marginal density plots
x_density <- ggplot(df_plot_new, aes(x = MAT_diff, fill = scen, color = scen)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  theme_void() +
  theme(legend.position = "none")  # Remove legend for x_density plot

y_density <- ggplot(df_plot_new, aes(x = summervpd_diff, fill = scen, color = scen)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  coord_flip() +
  theme_void() +
  labs(color = "Scenario") +
  theme(legend.position = "none")  # Remove legend for y_density plot

# Combine the plots using patchwork
combined_plot <- x_density + plot_spacer() + base_plot + y_density +
  plot_layout(ncol = 2, widths = c(4, 1), heights = c(1, 4))

# Print the combined plot (without legends)
print(combined_plot)

ggsave(combined_plot, filename = paste0(path, "/11_figures/FigS4.png"), width = 5, height = 5)



### end --------------------------------------------------------------------------
