# ------------------------------------------------------------------------------
# Script Name: Uncertainty analysis
# Author: Marc Grünig
# Date: 8.9.2025
#
# Input Files:
#   - Processed annual disturbance data
#
# Output:
#   - Figure S24
# ------------------------------------------------------------------------------


### libraries & settings -------------------------------------------------------

library(terra)
library(tidyverse)
library(sf)
library(stars)
library(ggplot2)
library(dplyr)


# paths
path <- "/.../"


### load data ------------------------------------------------------------------
df_plot_bbtl <- read_csv(paste0(path, "/10_results/bbtl_results_final.csv"))
df_plot_bbtl <- df_plot_bbtl %>%
  filter(sim != 0) %>% 
  mutate(Year = year + 2019) %>% 
  group_by(Year, sim, scen) %>% 
  summarize(dist_area = sum(dist_area, na.rm = T)) %>% 
  mutate(agent = "bbtl")

df_plot_fire <- read_csv(paste0(path, "/10_results/fire_results_final.csv"))
df_plot_fire <- df_plot_fire %>% 
  filter(sim != 0) %>% 
  mutate(Year = year + 2019) %>% 
  group_by(Year, sim, scen) %>% 
  summarize(dist_area = sum(burned_area, na.rm = T)) %>% 
  mutate(agent = "fire")


df_plot_wind <- read_csv(paste0(path, "/10_results/wind_results_final.csv"))
df_plot_wind <- df_plot_wind %>% 
  filter(sim != 0) %>% 
  mutate(Year = year + 2019) %>% 
  group_by(Year, sim, scen) %>% 
  summarize(dist_area = sum(dist_area, na.rm = T)) %>% 
  mutate(agent = "wind")


df_plot_all <- rbind(df_plot_fire, df_plot_wind, df_plot_bbtl)


id_map <- tibble(
  sim = 1:90,
  RCP = rep(c("rcp85", "rcp45", "rcp26", "rcp85", "rcp45", "rcp26", "rcp85", "rcp45", "rcp26"), times = 10),
  GCM = rep(c("ichec", "ichec", "ichec", "mpi", "mpi", "mpi", "ncc", "ncc", "ncc"), times = 10)
)


df <- df_plot_all %>%
  left_join(id_map, by = "sim")

df <- df %>%
  mutate(time_period = cut(Year,
                           breaks = c(2021, 2040, 2060, 2080, 2101),
                           labels = c("2021–2040", "2041–2060", "2061–2080", "2081–2100"),
                           right = FALSE))



## violin plot
my_cols <- c("RCP2.6" = "#ed9b49", "RCP4.5" = "#8d1c06", "RCP8.5" = "#3c0d03")

p1_all <- ggplot(df, aes(x = scen, y = dist_area, fill = scen)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  scale_fill_manual(values = my_cols) +
  labs(x = "", y = "area disturbed [ha]", fill = "") +
  theme_minimal(base_size = 14)
p1_all
ggsave(p1_all, filename = paste0(path, "/11_figures/figures/FigS24a.png"), width = 7.5, height = 7.5)
write_csv(df, paste0(path, "11_figures/figure_data/FigS24a.csv"))




# add sim id to see how much variance comes just from the years
# variance decomposition function on raw data
get_variance_explained_raw <- function(data) {
  data <- data %>%
    mutate(RCP = factor(RCP),
           GCM = factor(GCM),
           sim = factor(sim),
           Year = factor(Year))
  
  # Run ANOVA on the raw 'burned_area' variable (not averaged)
  model <- aov(dist_area ~ RCP + GCM + Year + sim, data = data)
  
  ss <- summary(model)[[1]][, "Sum Sq"]
  total_ss <- sum(ss)
  
  tibble(
    RCP = ss[1] / total_ss,
    GCM = ss[2] / total_ss,
    Year = ss[3] / total_ss,
    # Year = ss[1] / total_ss,
    Replicate = ss[4] / total_ss,
    Residuals = ss[5] / total_ss
  ) %>%
    pivot_longer(everything(), names_to = "Source", values_to = "Proportion")
}

df %>% group_by(RCP, GCM, Year, time_period, sim) %>% 
  summarize(n = n())

decomp_df <- df %>%
  group_by(time_period) %>%
  group_modify(~ get_variance_explained_raw(.x)) %>%
  ungroup()


# scaled - explained variance
decomp_df_scaled <- decomp_df %>%
  # drop residuals
  filter(Source != "Residuals") %>%
  # rescale proportions to sum to 1 within each time_period
  group_by(time_period) %>%
  mutate(Proportion = Proportion / sum(Proportion)) %>%
  ungroup()

decomp_df %>% filter(Source == "Residuals")

# Plot without residuals, scaled to 100%
p3 <- ggplot(decomp_df_scaled, aes(x = time_period, y = Proportion, fill = Source)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    x = "Time Period",
    y = "Proportion of explained variance",
    fill = "Source"
  ) +
  theme_minimal(base_size = 14)
p3
ggsave(p3, filename = paste0(path, "/11_figures/figures/FigS24b.png"), width = 7.5, height = 7.5)
write_csv(decomp_df_scaled, paste0(path, "11_figures/figure_data/FigS24b.csv"))

### end ------------------------------------------------------------------------
