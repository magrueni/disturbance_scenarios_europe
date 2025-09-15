# ------------------------------------------------------------------------------
# Script Name: Fig4_structural_trends.R
#
# Description: This script processes the outputs from script old_young_forests.R
#              and plot the trends per biome
#
# Author: Marc Grünig
# Date: 7.2.2025
#
#
# Output:
#   - figure 4a & 4b, S19 - 23
# ------------------------------------------------------------------------------


### libraries ------------------------------------------------------------------
library(terra)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(sf)
library(readr)
library(tidyverse)
library(MetBrewer)


path <- "/.../"


### load climate data ----------------------------------------------------------
timestep <- "1981-2010"
scenarios_hist <- c("ICHEC-EC-EARTH_historical", "NCC-NorESM1-M_historical", "MPI-M-MPI-ESM-LR_historical")

temp <- list()
for(s in 1:length(scenarios_hist)){
  
  fls <- read_csv(paste0(path, "/05_clim_data/annual_climate/svd_clim_", scenarios_hist[[s]], "_v3_aux_biascorr.csv"))
  temp_df <- fls %>% 
    group_by(climateId) %>% 
    summarize(MAT = mean(MAT),
              ANP = mean(ANP) * 365,
              summervpd = mean(summervpd)) %>% 
    dplyr::select(climateId, MAT, ANP, summervpd)
  temp[[s]] <- temp_df
  
}


temp <- do.call(rbind, temp)
MAT_hist <- mean(temp$MAT)
ANP_hist <- mean(temp$ANP)
summervpd_hist <- mean(temp$summervpd)


# future data
scenarios_fut <- c("ICHEC-EC-EARTH_rcp_2_6", "NCC-NorESM1-M_rcp_2_6", "MPI-M-MPI-ESM-LR_rcp_2_6",
                   "ICHEC-EC-EARTH_rcp_4_5", "NCC-NorESM1-M_rcp_4_5", "MPI-M-MPI-ESM-LR_rcp_4_5", 
                   "ICHEC-EC-EARTH_rcp_8_5", "NCC-NorESM1-M_rcp_8_5", "MPI-M-MPI-ESM-LR_rcp_8_5")

temp <- list()
for(s in 1:length(scenarios_fut)){
  
  fls <- read_csv(paste0(path, "/05_clim_data/annual_climate/svd_clim_", scenarios_fut[[s]], "_v3_aux_biascorr.csv"))
  temp_df <- fls %>%
    mutate(ANP = ANP*365) %>%
    dplyr::select(year, climateId, MAT, ANP, summervpd) %>% 
    group_by(year) %>% 
    summarize(MAT = mean(MAT),
              ANP = mean(ANP),
              summervpd = mean(summervpd)) %>% 
    mutate(scen = paste(scenarios_fut[[s]]))
  temp[[s]] <- temp_df
  
}

temp_fut <- do.call(rbind, temp)

temp_fut_df <- temp_fut %>% 
  mutate(MAT_diff = MAT - MAT_hist,
         ANP_diff = ANP - ANP_hist,
         summervpd_diff = summervpd - summervpd_hist) %>% 
  mutate(gcm = ifelse(grepl("ICHEC", scen), "ichec",
                      ifelse(grepl("NCC", scen), "ncc", "mpi"))) %>% 
  mutate(scen = ifelse(grepl("rcp_2_6", scen), "RCP2.6",
                       ifelse(grepl("rcp_4_5", scen), "RCP4.5", "RCP8.5")))

temp_fut_df_10 <- temp_fut_df %>% 
  mutate(time_step = (year %/% 10) * 10) %>% 
  group_by(scen, gcm, time_step) %>% 
  summarize(MAT = mean(MAT, na.rm = T),
            MAT_diff = mean(MAT_diff),
            ANP = mean(ANP, na.rm = T),
            ANP_diff = mean(ANP_diff),
            VPD = mean(summervpd, na.rm = T),
            VPD_diff = mean(summervpd_diff),)



### load biodiv data -----------------------------------------------------------

# in case the step above was done already
data_all <- list()

for(i in c(1:120)){
  
  x <-  read_sf(paste0(path, "/09_svd_simulations/results_eu/biodiv/undist_vs_recently_dist_forest_", i, ".gpkg"))
    
  data_all[[i]] <- x %>% st_drop_geometry()
  
}

data_all_df <- do.call(bind_rows, data_all)

# get some numbers
data_all_df %>% 
  group_by(year, rcp) %>% 
  summarise(old = mean(undist_prct, na.rm = T),
            young = mean(dist_rec_prct, na.rm = T)) %>% 
  filter(year == 2100)



# load biomes
hex_ecoreg <- st_read(paste0(path, "/07_reference_grids/biomes_hex.gpkg"))

df_plot <- data_all_df %>% 
  left_join(., hex_ecoreg, by = c("gridid")) %>% 
  group_by(year, rcp, gcm, biome) %>% 
  summarize(undist_prct = mean(undist_prct, na.rm = T),
            dist_rec_prct = mean(dist_rec_prct, na.rm = T)) %>% 
  dplyr::select(year, rcp, gcm, biome, dist_rec_prct) %>% 
  mutate(rcp = ifelse(rcp == "historical", "Historical", rcp))

# get some numbers
df_plot %>% 
  group_by(year, rcp, biome) %>% 
  summarise(young = mean(dist_rec_prct, na.rm = T)) %>% 
  filter(biome == "Mediterranean")



df_plot_dist_rec <- df_plot %>% 
  pivot_wider(names_from = c(rcp, biome), values_from = dist_rec_prct) %>%
  mutate(across(starts_with("rcp26"), ~ . - get(sub("rcp26", "Historical", cur_column())), .names = "difference_{col}"),
         across(starts_with("rcp45"), ~ . - get(sub("rcp45", "Historical", cur_column())), .names = "difference_{col}"),
         across(starts_with("rcp85"), ~ . - get(sub("rcp85", "Historical", cur_column())), .names = "difference_{col}"),
         across(starts_with("rcp26"), ~ (. - get(sub("rcp26", "Historical", cur_column()))) / get(sub("rcp26", "Historical", cur_column())) * 100, .names = "relative_change_{col}"),
         across(starts_with("rcp45"), ~ (. - get(sub("rcp45", "Historical", cur_column()))) / get(sub("rcp45", "Historical", cur_column())) * 100, .names = "relative_change_{col}"),
         across(starts_with("rcp85"), ~ (. - get(sub("rcp85", "Historical", cur_column()))) / get(sub("rcp85", "Historical", cur_column())) * 100, .names = "relative_change_{col}")) %>%
  pivot_longer(cols = starts_with("relative_change"), names_to = "rcp_biome", values_to = "diff") %>%
  mutate(rcp_biome = str_remove(rcp_biome, "relative_change_")) %>%
  separate(rcp_biome, into = c("rcp", "biome"), sep = "_(?=[^_]+$)", extra = "merge") %>%
  mutate(rcp = case_when(
    rcp == "rcp26" ~ "RCP2.6",
    rcp == "rcp45" ~ "RCP4.5",
    rcp == "rcp85" ~ "RCP8.5"
  )) %>%
  dplyr::select(year, gcm, rcp, biome, rec_dist_diff = diff)


df_plot_new <- df_plot_dist_rec %>% rename("time_step" = "year",
                                      "scen" = "rcp") %>% 
  left_join(., temp_fut_df_10, by = c("time_step", "scen", "gcm")) %>% 
  drop_na() %>% 
  filter(biome != "Temperate Grasslands")

df_plot_europe <- df_plot_new %>% 
  group_by(time_step, gcm, scen) %>% 
  summarize(rec_dist_diff = mean(rec_dist_diff),
            MAT = mean(MAT),
            MAT_diff = mean(MAT_diff),
            ANP = mean(ANP),
            ANP_diff = mean(ANP_diff),
            VPD = mean(VPD),
            VPD_diff = mean(VPD_diff)) %>% 
  mutate(biome = "Europe")

# extract some numbers
df_plot_europe %>% 
  group_by(time_step, scen) %>% 
  summarize(young_forests = mean(rec_dist_diff)) %>% 
  filter(time_step %in% c(2100))

df_plot_new <- rbind(df_plot_new, df_plot_europe)
df_plot_new <- df_plot_new %>% 
  mutate(biome = case_when(
    biome == "Mediterranean" ~ "Mediterranean",
    biome == "Temperate Broadleaf" ~ "Temperate broadleaved",
    biome == "Temperate Coniferous" ~ "Temperate coniferous",
    biome == "Boreal Forests" ~ "Boreal",
    biome == "Tundra" ~ "Tundra",
    biome == "Europe" ~ "Europe",
    TRUE ~ biome   # keep original if no match
  ))

df_plot_new$biome <- factor(df_plot_new$biome,
                            levels = c("Mediterranean", "Temperate broadleaved", "Temperate coniferous",
                                       "Boreal", "Tundra", "Europe"))


df_plot_new %>% filter(biome == "Mediterranean") %>% 
  filter(MAT_diff >= 4) %>% 
  group_by(biome) %>% 
  summarize(diff = mean(rec_dist_diff))


write_csv(df_plot_new, paste0(path, "/11_figures/figure_data/Fig4a_FigS23a.csv"))


## fig 4a
cols <- c(met.brewer("Isfahan1", 5), "#0065bd", "lightgrey")
p1 <- ggplot(df_plot_new %>% filter(biome != "Temperate Grasslands"),
             aes(x = MAT_diff, y = rec_dist_diff, color = biome, fill = biome, group = biome)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = cols) +
  geom_smooth(method = "gam", se = TRUE, alpha = 0.1) +
  labs(
    x = "Δ Temperature [°C]",
    y = "Δ Young forests [%]",
    color = "Biome",
    fill = "Biome"
  ) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(
    text = element_text(size = 10),
    legend.position = c(0.05, 0.98),
    legend.justification = c("left", "top")
  )

p1

ggsave(p1, filename = paste0(path, "/11_figures/Fig4a.png"), width = 5, height = 5)

p2 <- ggplot(df_plot_new %>% filter(biome != "Temperate Grasslands"),
             aes(x = VPD_diff, y = rec_dist_diff, color = biome, fill = biome, group = biome)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = cols) +
  geom_smooth(method = "gam", se = TRUE, alpha = 0.1) +
  labs(
    x = "Δ Summer VPD [kPa]",
    y = "Δ Young forests [%]",
    color = "Biome",
    fill = "Biome"
  ) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(
    text = element_text(size = 10),
    legend.position = c(0.05, 1.02),  # Adjust position as needed
    legend.justification = c("left", "top")
  )

p2

ggsave(p2, filename = paste0(path, "/11_figures/FigS23a.png"), width = 5, height = 5)


# timeseries
cols <- c(met.brewer("Isfahan1", 5), "#0065bd", "lightgrey")
df_rev <- df_plot_new %>%
  filter(biome != "Temperate Grasslands")

p3 <- ggplot(df_rev) +
  geom_boxplot(aes(x = as.factor(time_step), y = rec_dist_diff, 
                   color = scen, fill = scen), alpha = 0.5)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") + # Add horizontal line at y = 0
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
  labs(
    x = "Year",
    y = "Δ Young forests [%]"
  )+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+
  facet_wrap(~biome)

p3

ggsave(p3, filename = paste0(path, "/11_figures/FigS20.png"), width = 7, height = 5)
write_csv(df_rev, paste0(path, "/11_figures/figure_data/FigS20.csv"))



### now plot the old forests ----------------------------------------------------
data_all_df_end <- data_all_df %>% filter(year == 2100)

hex_ecoreg <- st_read(paste0(path, "/reference_grids/biomes_hex.gpkg"))

data_all_df_end <- data_all_df_end %>% 
  left_join(., hex_ecoreg, by = c("gridid"))


df_plot <- data_all_df_end %>% 
  group_by(year, rcp, gcm, rep, biome) %>% 
  summarize(undist_prct = mean(undist_prct, na.rm = T)) %>% 
  dplyr::select(year, rcp, gcm, rep, biome, undist_prct)

df_plot <- df_plot %>% 
  mutate(biome = case_when(
    biome == "Mediterranean" ~ "Mediterranean",
    biome == "Temperate Broadleaf" ~ "Temperate broadleaved",
    biome == "Temperate Coniferous" ~ "Temperate coniferous",
    biome == "Boreal Forests" ~ "Boreal",
    biome == "Tundra" ~ "Tundra",
    biome == "Europe" ~ "Europe",
    TRUE ~ biome   # keep original if no match
  ))


# Combine data and pivot wider with rcp and biome
df_plot_dist_rec <- df_plot %>% mutate(rcp = ifelse(rcp == "historical", "Historical", rcp)) %>% 
  pivot_wider(names_from = c(rcp, biome), values_from = undist_prct)

# Extract the unique biomes and RCP scenarios from the data
biomes <- unique(df_plot$biome)
rcps <- c("rcp26", "rcp45", "rcp85")

# Loop through each RCP and biome combination to calculate differences and relative changes
for (rcp in rcps) {
  for (biome in biomes) {
    df_plot_dist_rec <- df_plot_dist_rec %>%
      mutate(
        !!paste0("difference_", rcp, "_", biome) := !!sym(paste0(rcp, "_", biome)) - !!sym(paste0("Historical_", biome)),
        !!paste0("relative_change_", rcp, "_", biome) := ( !!sym(paste0(rcp, "_", biome)) - !!sym(paste0("Historical_", biome)) ) / !!sym(paste0("Historical_", biome)) * 100
      )
  }
}

# Select only the relevant columns
df_plot_dist_rec <- df_plot_dist_rec %>% ungroup() %>% 
  dplyr::select(year, gcm, rep, starts_with("difference_"), starts_with("relative_change_"))

# Pivot longer to have a single column for the relative changes
df_plot_dist_rec <- df_plot_dist_rec %>%
  pivot_longer(cols = starts_with("relative_change_"), names_to = "rcp_biome", values_to = "diff") %>%
  separate(rcp_biome, into = c("type", "rcp", "biome"), sep = "_", extra = "merge") %>%
  mutate(rcp = case_when(
    rcp == "rcp26" ~ "RCP2.6",
    rcp == "rcp45" ~ "RCP4.5",
    rcp == "rcp85" ~ "RCP8.5"
  )) %>%
  dplyr::select(year, gcm, rep, biome, rcp, undist_prct_diff = diff) %>% 
  mutate(rcp = ifelse(grepl("rcp26", biome), "RCP2.6", 
                      ifelse(grepl("rcp45", biome), "RCP4.5", "RCP8.5"))) %>% 
  mutate(Biome = ifelse(grepl("Boreal", biome), "B",
                        ifelse(grepl("Mediterranean", biome), "M",
                               ifelse(grepl("Temperate broadleaved", biome), "TB",
                                      ifelse(grepl("Tundra", biome), "T",
                                             ifelse(grepl("Temperate coniferous", biome), "TC", "TG")))))) 


df_plot_dist_rec %>% 
  group_by(Biome, rcp, year) %>% 
  summarize(undist_mean = mean(undist_prct_diff, na.rm = T),
            undist_sd = sd(undist_prct_diff, na.rm = T))

mean_diff_biomes <- df_plot_dist_rec %>% filter(rcp == "RCP8.5") %>% 
  group_by(Biome) %>%
  summarise(mean_undist = mean(undist_prct_diff, na.rm = TRUE)) %>%
  arrange(mean_undist) %>%
  pull(Biome)

df_plot_dist_rec <- df_plot_dist_rec %>% 
  mutate(Biome = factor(Biome, levels = mean_diff_biomes)) %>% 
  filter(Biome != "TG")


write_csv(df_plot_dist_rec, paste0(path, "/11_figures/figure_data/FigS19.csv"))


p2 <- ggplot(df_plot_dist_rec) +
  geom_boxplot(aes(x = Biome, y = undist_prct_diff, 
                   color = rcp, fill = rcp), alpha = 0.5)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") + 
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
  ylab("Δ Old forests [%]") +
  theme_classic()


p2
ggsave(p2, filename = paste0(path, "/11_figures/FigS19.png"), width = 5, height = 5)


# get the numbers
df_plot_dist_rec %>% 
  group_by(Biome, rcp) %>% 
  summarize(mean = mean(undist_prct_diff, na.rm = T))




### line plot for undisturbed ---
df_plot <- data_all_df %>% 
  left_join(., hex_ecoreg, by = c("gridid")) %>% 
  group_by(year, rcp, gcm, biome) %>% 
  summarize(undist_prct = mean(undist_prct, na.rm = T),
            dist_rec_prct = mean(dist_rec_prct, na.rm = T)) %>% 
  dplyr::select(year, rcp, gcm, biome, undist_prct) %>% 
  mutate(rcp = ifelse(rcp == "historical", "Historical", rcp))

df_plot_dist_rec <- df_plot%>% 
  pivot_wider(names_from = c(rcp, biome), values_from = undist_prct) %>%
  mutate(across(starts_with("rcp26"), ~ . - get(sub("rcp26", "Historical", cur_column())), .names = "difference_{col}"),
         across(starts_with("rcp45"), ~ . - get(sub("rcp45", "Historical", cur_column())), .names = "difference_{col}"),
         across(starts_with("rcp85"), ~ . - get(sub("rcp85", "Historical", cur_column())), .names = "difference_{col}"),
         across(starts_with("rcp26"), ~ (. - get(sub("rcp26", "Historical", cur_column()))) / get(sub("rcp26", "Historical", cur_column())) * 100, .names = "relative_change_{col}"),
         across(starts_with("rcp45"), ~ (. - get(sub("rcp45", "Historical", cur_column()))) / get(sub("rcp45", "Historical", cur_column())) * 100, .names = "relative_change_{col}"),
         across(starts_with("rcp85"), ~ (. - get(sub("rcp85", "Historical", cur_column()))) / get(sub("rcp85", "Historical", cur_column())) * 100, .names = "relative_change_{col}")) %>%
  pivot_longer(cols = starts_with("relative_change"), names_to = "rcp_biome", values_to = "diff") %>%
  mutate(rcp_biome = str_remove(rcp_biome, "relative_change_")) %>%
  separate(rcp_biome, into = c("rcp", "biome"), sep = "_(?=[^_]+$)", extra = "merge") %>%
  mutate(rcp = case_when(
    rcp == "rcp26" ~ "RCP2.6",
    rcp == "rcp45" ~ "RCP4.5",
    rcp == "rcp85" ~ "RCP8.5"
  )) %>%
  dplyr::select(year, gcm, rcp, biome, rec_dist_diff = diff)



# df_plot_fin <- left_join(df_plot_dist_rec, df_plot_undist, by = c("year", "gcm", "rcp"))
df_plot_fin <- df_plot_dist_rec

df_plot_new <- df_plot_fin %>% rename("time_step" = "year",
                                      "scen" = "rcp") %>% 
  left_join(., temp_fut_df_10, by = c("time_step", "scen", "gcm")) %>% 
  drop_na() %>% 
  filter(biome != "Temperate Grasslands")


df_plot_new %>% 
  group_by(biome, time_step, scen) %>% 
  summarize(young_forests = mean(rec_dist_diff)) %>% 
  filter(time_step %in% c(2100))


df_plot_europe <- df_plot_new %>% 
  group_by(time_step, gcm, scen) %>% 
  summarize(rec_dist_diff = mean(rec_dist_diff),
            MAT = mean(MAT),
            MAT_diff = mean(MAT_diff),
            ANP = mean(ANP),
            ANP_diff = mean(ANP_diff),
            VPD = mean(VPD),
            VPD_diff = mean(VPD_diff)) %>% 
  mutate(biome = "Europe")

# extract some numbers
df_plot_europe %>% 
  group_by(time_step, scen) %>% 
  summarize(young_forests = mean(rec_dist_diff)) %>% 
  filter(time_step %in% c(2100))



df_plot_new <- rbind(df_plot_new, df_plot_europe)

df_plot_new$biome <- factor(df_plot_new$biome,
                            levels = c("Mediterranean", "Temperate Broadleaf", "Temperate Coniferous",
                                       "Boreal Forests", "Tundra", "Europe"))

write_csv(df_plot_new, paste0(path, "/11_figures/figure_data/Fig4b_FigS23b.csv"))



cols <- c(met.brewer("Isfahan1", 5), "#0065bd", "lightgrey")
p1 <- ggplot(df_plot_new %>% filter(biome != "Temperate Grasslands"),
             aes(x = MAT_diff, y = rec_dist_diff, color = biome, fill = biome, group = biome)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = cols) +
  geom_smooth(method = "gam", se = TRUE, alpha = 0.1) +
  labs(
    x = "Δ Temperature [°C]",
    y = "Δ Old forests [%]",
    color = "Biome",
    fill = "Biome"
  ) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(
    text = element_text(size = 10),
    legend.position = "none",  # Adjust position as needed
    legend.justification = c("left", "top")
  )

p1

ggsave(p1, filename = paste0(path, "/11_figures/Fig4b.png"), width = 5, height = 5)
ggsave(p1, filename = paste0(path, "/11_figures/Fig4b.pdf"), width = 5, height = 5)

p2 <- ggplot(df_plot_new %>% filter(biome != "Temperate Grasslands"),
             aes(x = VPD_diff, y = rec_dist_diff, color = biome, fill = biome, group = biome)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = cols) +
  geom_smooth(method = "gam", se = TRUE, alpha = 0.1) +
  labs(
    x = "Δ summer VPD [kPa]",
    y = "Δ Old forests [%]",
    color = "Biome",
    fill = "Biome"
  ) +
  scale_fill_manual(values = cols) +
  theme_classic() +
  theme(
    text = element_text(size = 10),
    legend.position = "none",  # Adjust position as needed
    legend.justification = c("left", "top")
  )

p2

ggsave(p2, filename = paste0(path, "/11_figures/FigS23b.png"), width = 5, height = 5)




cols <- c(met.brewer("Isfahan1", 5), "#0065bd", "lightgrey")
df_rev <- df_plot_new %>%
  filter(biome != "Temperate Grasslands")

p3 <- ggplot(df_rev) +
  geom_boxplot(aes(x = as.factor(time_step), y = rec_dist_diff, 
                   color = scen, fill = scen), alpha = 0.5)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") + # Add horizontal line at y = 0
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
  labs(
    x = "Year",
    y = "Δ Old forests [%]"
  )+
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+
  facet_wrap(~biome) +
  ylim(c(-10, 5))

p3

ggsave(p3, filename = paste0(path, "/11_figures/FigS21.png"), width = 7, height = 5)
write_csv(df_rev, paste0(path, "/11_figures/figure_data/FigS21.csv"))





### end ------------------------------------------------------------------------
