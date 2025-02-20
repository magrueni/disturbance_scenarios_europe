# --------------------------------------------------------------
# Management Module Script 01: Management probabilities based on Increment threshold
# --------------------------------------------------------------
# Author: Marc Grünig
# Date: 2025-01-07
# 
# Description:
# This script processes a forest state lookup table to compute management-related parameters based on tree height.
# It calculates a height increment threshold and scales a management probability (pManagement) with height,
# setting it to zero for trees shorter than 20 meters. 
# The results are visualized with histograms and scatter plots,
# and the processed data is saved as a CSV file for use in forest management simulations.
#
# --------------------------------------------------------------



### libraries ----------------------------
library(tidyverse)
library(dplyr)
library(DBI)
library(terra)
library(stringr)


# define path
path <- "/.../"


# load state lookup
dnn_lookup <- read.csv(paste0(path, "/dnn/dnn_states_lookup.csv")) 

# we put an height increment threshold of 0.5 % of the height. 
states <- dnn_lookup %>% 
  rowwise() %>% 
  mutate(height = as.numeric(strsplit(state, "_")[[1]][3]) + 1) %>%  # middle value of the class
  mutate(increment_threshold = height * 0.005 * height/30)

# put 0 in missing state
states[1, c(3,4)] <- 0

# scale the pManagement with the height and set to 0 if below 20m
states <- states %>% mutate(pManagement = 0.33*height/30) %>% 
  dplyr::rename("heightIncrementThreshold" = "increment_threshold") %>% 
  rowwise() %>% 
  mutate(pManagement = ifelse(height < 20, 0, pManagement))

# View(states)
hist(states$heightIncrementThreshold)
hist(states$pManagement)

plot(states$height, states$pManagement)
plot(states$height, states$heightIncrementThreshold)


states <- states %>% 
  dplyr::select(stateID, heightIncrementThreshold, pManagement)
write.csv(states, paste0(path, "/mgmt_module/management.csv"), row.names = F)

