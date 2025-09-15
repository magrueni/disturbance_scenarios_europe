# --------------------------------------------------------------------------------------
# Script Name: Fig3a_stacked_bars.R
# Description: This script takes the data from the script 04_annual_outputs_processing.R
#              to plot the disturbed area per disturbance agent.
#
# Author: Marc Grünig
# Date: 7.2.2025
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
library(ggh4x)

# define paths
path <- "/.../"

### load data ------------------------------------------------------------------


plot_dat_agents_fut <- read_csv(paste0(path, "/10_results/all_agents_distarea_fut.csv"))
plot_dat_agents_hist <- read_csv(paste0(path, "/10_results/all_agents_distarea_hist.csv"))


# remove 2020 - As no data were available for the year 2020, disturbed area for bark beetle and wind were omitted for 2020 in this analysis. 
plot_dat_agents <- rbind(plot_dat_agents_fut, plot_dat_agents_hist) 
plot_dat_agents <- plot_dat_agents %>% filter(!(year == 35 & period == "2001-2020" & agent %in% c("Wind", "Barkbeetle")))


plot_dat <- plot_dat_agents %>% 
  group_by(sim, year, scen, period, agent) %>% 
  summarize(dist_area = sum(dist_area)) %>% 
  drop_na()

plot_dat %>% filter(scen == "RCP8.5") %>% 
  group_by(sim, scen, period, agent) %>%
  summarize(sum_area = sum(dist_area)) %>% arrange(sum_area)

mean_values <- plot_dat %>% 
  group_by(scen, year, agent, period) %>% 
  summarise(dist_area = mean(dist_area)) %>% 
  group_by(scen, period, agent) %>%
  summarize(mean_area = mean(dist_area, na.rm = T),
            n = n(),
            se = sd(dist_area, na.rm = T) / sqrt(n()),
            sd_area = sd(dist_area, na.rm = T)) %>%
  mutate(
    ci = 1.96 * se,
    ci_lower = mean_area - ci,
    ci_upper = mean_area + ci
  )

# Reorder the agents as specified
mean_values$agent <- factor(mean_values$agent, levels = c("Wind", "Barkbeetle", "Fire"))
mean_values <- mean_values %>% arrange(agent)

mean_values <- mean_values %>%
  arrange(scen, period, agent) %>%
  group_by(scen, period) %>%
  mutate(
    cum_mean_area = cumsum(mean_area),
    cum_ci_lower = cumsum(ci_lower),
    cum_ci_upper = cumsum(ci_upper),
    cum_mean_area_prev = lag(cum_mean_area, default = 0),
    cum_ci_lower_prev = lag(cum_ci_lower, default = 0),
    cum_ci_upper_prev = lag(cum_ci_upper, default = 0)
  ) %>%
  mutate(
    cum_mean_area_sd = cumsum(mean_area),
    cum_sd_area_sd = cumsum(sd_area),
    cum_mean_area_prev_sd = lag(cum_mean_area_sd, default = 0),
    cum_sd_area_prev_sd = lag(cum_sd_area_sd, default = 0)
  ) %>% 
  mutate(agent = ifelse(agent == "Barkbeetle", "Bark beetle", paste(agent)))


write_csv(mean_values, paste0(path, "/11_figures/figure_data/Fig3a.csv"))


### Plotting -------------------------------------------------------------------
met <- met.brewer("Cassatt2", 9)[7]
cols <- c("darkblue", met, "darkred")


# inspired from https://www.cedricscherer.com/2019/08/05/a-ggplot2-tutorial-for-beautiful-plotting-in-r/
stacked_bars <- ggplot(mean_values, aes(x = factor(period), y = mean_area, fill = factor(agent),
                                        color = factor(agent), group = scen)) +
  geom_col(alpha = 0.5) +
  geom_errorbar(aes(ymin = cum_mean_area_prev + (mean_area - se), 
                    ymax = cum_mean_area_prev + (mean_area + se)), width = 0.2) +
  geom_col(data = subset(mean_values, period == "2001-2020" | period == "1986-2000"),
           aes(fill = "Historical", color = "Historical"), alpha = 0.5, linewidth = 1) +
  scale_color_manual(values = c("Fire" = "darkred",
                                "Bark beetle" = met,
                                "Wind" = "darkblue",
                                "Historical" = "grey"),
                     name = "Agent",
                     breaks = c("Fire", "Bark beetle", "Wind", "Historical"),
                     labels = c("Fire", "Bark beetle", "Wind", "Historical")) +
  scale_fill_manual(values = c("Fire" = "darkred",
                               "Bark beetle" = met,
                               "Wind" = "darkblue",
                               "Historical" = "grey"),
                    name = "Agent",
                    breaks = c("Fire", "Bark beetle", "Wind", "Historical"),
                    labels = c("Fire", "Bark beetle", "Wind", "Historical")) +
  labs(title = "",
       x = "Time period",
       y = bquote("Annual area disturbed [ha yr"^-1*"]")) +
  theme_classic() +
  theme(legend.position = "right") +
  facet_wrap2(~scen, strip = strip_themed(clip = "inherit",
                                          background_x = elem_list_rect(fill = c("#ed9b49", "#8d1c06", "#3c0d03")))) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        axis.title.x = element_text(size = 12, margin = margin(t = 20)),              # Adjust x-axis title size
        axis.title.y = element_text(size = 12),              # Adjust y-axis title size
        axis.text = element_text(size = 10),                 # Adjust axis text size
        legend.text = element_text(size = 10),               # Adjust legend text size
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold", color = "white")
  )


stacked_bars
ggsave(stacked_bars, filename = paste0(path, "/11_figures/Fig3a.pdf"), width = 8, height = 8)


### end ------------------------------------------------------------------------