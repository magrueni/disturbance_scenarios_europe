# ------------------------------------------------------------------------------
# Script Name: Feedbacks and Interactions
# Author: Marc Grünig
# Date: 8.9.2025
#
# Description: This script compares outputs from simulation without bark
#              beetle wind interactions with original simulations. And
#              static vegatation simulations to originals simulations.
#

# Input Files:
#   - raw simulation data from additional experiment simulations
#
# Output:
#   - Figure S35 & S36
# ------------------------------------------------------------------------------




### libraries ------------------------------------------------------------------
library(terra)
library(tidyverse)
library(sf)
library(stars)
library(readr)



# paths 
path_lorien <- "/.../"
path <- "/.../"
path_int <- "/.../"



### interaction experiment -----------------------------------------------------

# load original simulations
df_out <- list()

for(i in 1:3){
  
  if(!file.exists(paste0(path_lorien, "/output_sim_eu_", i, "/barkbeetle.csv"))){next}
  df <- read_csv(paste0(path_lorien, "/output_sim_eu_", i, "/barkbeetle.csv"))
  df <- df %>% group_by(year) %>% summarize(bb_area = sum(n_impact))
  df <- df %>% mutate(Year = as.numeric(year) + 2019) %>% 
    mutate(run = paste(i),
           scen = ifelse(i %in% seq(1, 90, 3), "rcp85", 
                         ifelse(i %in% seq(2, 90, 3), "rcp45", "rcp26")))
  df_impact <- df %>% dplyr::select(Year, bb_area, run, scen) %>% arrange(Year) %>% drop_na()
  df_out[[i]] <- df_impact
}

df_plot_bb <- do.call(rbind, df_out)

df_plot_bb <- df_plot_bb %>%
  mutate(scen = ifelse(scen == "rcp85", "RCP8.5",
                       ifelse(scen == "rcp45", "RCP4.5", "RCP2.6"))) %>% 
  mutate(source = "orig")


# load interaction experiment simulations
for(i in 1:3){
  
  if(!file.exists(paste0(path_int, "/output_sim_eu_", i, "/barkbeetle.csv"))){next}
  df <- read_csv(paste0(path_int, "/output_sim_eu_", i, "/barkbeetle.csv"))
  df <- df %>% group_by(year) %>% summarize(bb_area = sum(n_impact ))
  df <- df %>% mutate(Year = as.numeric(year) + 2019) %>% 
    mutate(run = paste(i),
           scen = ifelse(i %in% seq(1, 90, 3), "rcp85", 
                         ifelse(i %in% seq(2, 90, 3), "rcp45", "rcp26")))
  df_impact <- df %>% dplyr::select(Year, bb_area, run, scen) %>% arrange(Year) %>% drop_na()
  df_out[[i]] <- df_impact
}

df_plot_bbtl_test <- do.call(rbind, df_out)
df_plot_bbtl_test <- df_plot_bbtl_test %>% filter(Year <= 2100) %>%
  mutate(scen = ifelse(scen == "rcp85", "RCP8.5",
                       ifelse(scen == "rcp45", "RCP4.5", "RCP2.6"))) %>% 
  mutate(source = "test")



df_plot_all_bb <- rbind(df_plot_bb, df_plot_bbtl_test)

# create boxplot
numbers <- df_plot_all_bb %>% 
  pivot_wider(., names_from = "source", values_from = "bb_area") %>% 
  mutate(ratio = test/orig)
mean(numbers$ratio, na.rm = T)

# get the percentage that is explained by the interactions
mean(numbers$orig)/mean(numbers$test, na.rm = T)

# Sample data (assuming df_plot_all_bb is already loaded)
df_plot_all_bb <- df_plot_all_bb %>% mutate(source = ifelse(source == "orig", "with interaction", "without interaction")) %>% 
  rename("disturbed_area" = "bb_area")

write_csv(df_plot_all_bb, paste0(path, "/11_figures/figure_data/FigS36.csv"))

int_plot <- ggplot(df_plot_all_bb, aes(x = scen, y = disturbed_area, fill = source)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Effect of wind-barkbeetle interactions",
       x = "",
       y = "disturbed area [ha]") +
  scale_fill_brewer(palette = "Set2")  +
  theme(text = element_text(size = 12),
        axis.title.x = element_text(vjust = -1),
        axis.title.y = element_text(vjust = 1))

int_plot
ggsave(int_plot, filename = paste0(path, "/11_figures/FigS36.png"), width = 6, height = 6)






### feedback experiment --------------------------------------------------------

master_tab <- read.csv(paste0(path, "/09_svd_simulations/svd_simulations_ids.csv"), sep = ";")

# original file
orig_fire <- rast(paste0(path, "/09_svd_simulations/raw/fire/fire_80_rcp85_mpi_2_summed.tif"))
orig_wind <- rast(paste0(path, "/09_svd_simulations/raw/wind/wind_80_rcp85_mpi_2_summed.tif"))
orig_bbtl <- rast(paste0(path, "/09_svd_simulations/raw/bbtl/bbtl_80_rcp85_mpi_2_summed.tif"))
orig <- app(c(orig_fire, orig_wind, orig_bbtl), "sum") 

# clean
tmpFiles(remove = T)
gc()

# test data
test_fire <- rast(paste0(path, "/09_svd_simulations/output_final_v9_nochange/ccfire/firegrid_80_sim_eu_15.tif"))
test_wind <- rast(paste0(path, "/09_svd_simulations/output_final_v9_nochange/ccwind/wind_year_80_sim_eu_15.tif"))

test_bbtl <- rast()
for(i in c(10, 20, 30, 40, 50, 60, 70, 80)){
  
  test <- rast(paste0(path, "/09_svd_simulations/output_final_v9_nochange/bb_out/bbgrid_", i, "_sim_eu_15.tif"))
  test[test < (as.numeric(i)-10)] <- 0
  test[test > 0] <- 1
  plot(test)
  test_bbtl <- c(test_bbtl, test)
  gc()
}

test_bbtl <- app(test_bbtl, "sum")
test <- app(c(test_fire, test_wind, test_bbtl), "sum") 

## convert to dataframes
orig_df <- as.data.frame(orig, xy  = T)
colnames(orig_df)[3] <- "disturbed"

test_df <- as.data.frame(test, xy  = T)
colnames(test_df)[3] <- "disturbed"


nrow(orig_df[orig_df$disturbed > 0,])
nrow(orig_df[orig_df$disturbed > 1,])

nrow(test_df[test_df$disturbed > 0,])
nrow(test_df[test_df$disturbed > 1,])

proba_test <- nrow(test_df[test_df$disturbed > 1,])/nrow(test_df)
proba_orig <- nrow(orig_df[orig_df$disturbed > 1,])/nrow(orig_df)

100/proba_orig*proba_test

proba_test <- nrow(test_df[test_df$disturbed > 2,])/nrow(test_df)
proba_orig <- nrow(orig_df[orig_df$disturbed > 2,])/nrow(orig_df)

100/proba_orig*proba_test



nrow(orig_df[orig_df$disturbed > 1,])/nrow(test_df[test_df$disturbed > 1,])

orig_df_hist <- orig_df %>%
  filter(disturbed > 1) %>% 
  filter(disturbed < 20)
test_df_hist <- test_df %>% 
  filter(disturbed > 1) %>% 
  filter(disturbed < 20)


par(mfrow = c(1,2))
hist(orig_df_hist[, "disturbed"], main = "Original")
hist(test_df_hist[, "disturbed"], main = "Statische Veg")


# Histogram for Original
p1 <- ggplot(orig_df_hist, aes(x = disturbed)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  theme_minimal() +
  ggtitle("Dynamic vegetation") +
  xlab("Disturbance frequency") +
  ylab("Number of grid cells") +
  scale_x_continuous(breaks = seq(2, 20, by = 2), limits = c(1, 20)) +
  scale_y_continuous(breaks = seq(0, 20000000, by = 5000000), limits = c(0, 6000000))
p1
# Histogram for static veg
p2 <- ggplot(test_df_hist, aes(x = disturbed)) +
  geom_histogram(binwidth = 1, fill = "red", color = "black", alpha = 0.7) +
  theme_minimal() +
  ggtitle("Static vegetation") +
  xlab("Disturbance frequency") +
  ylab("Number of grid cells") +
  scale_x_continuous(breaks = seq(2, 20, by = 2), limits = c(1, 20)) +
  scale_y_continuous(breaks = seq(0, 20000000, by = 5000000), limits = c(0, 6000000))

p2
# Arrange the plots side by side
library(gridExtra)
plot_feedback <- grid.arrange(p1, p2, ncol = 2)

ggsave(plot_feedback, filename = paste0(path, "/11_figures/FigS35.png"), width = 6, height = 6)

out_df <- rbind(orig_df, test_df)
write_csv(out_df, paste0(path, "/11_figures/figure_data/FigS35.csv"))


### end ------------------------------------------------------------------------
