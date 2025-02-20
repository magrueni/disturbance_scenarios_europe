# --------------------------------------------------------------
# Initial landscape Scipt 06: DNN states for SVD
# --------------------------------------------------------------
# Author: Marc Grünig
# Date: 2025-01-07
# 
#
# Description:
# This script creates tables with information on all forests states for SVD
#
# --------------------------------------------------------------



### libraries ------------------------------------------------------------------
library(readr)
library(tidyverse)

path <- "/.../"


dnn_lookup <- read_csv(paste0(path, "/initial_states/states_lookup_pruned.csv"))
colnames(dnn_lookup) <- c("code", "stateId")
dim(dnn_lookup)

# write result
write.csv(dnn_lookup, paste0(path, "/dnn/states.csv"))


# create states europe file with information on structure comp, etc.
states_europe <- dnn_lookup %>% dplyr::select(stateId, unique.str = code) %>% 
  rowwise() %>% 
  mutate(composition = strsplit(unique.str, "_")[[1]][1],
         structure = strsplit(unique.str, "_")[[1]][3],
         fct = strsplit(unique.str, "_")[[1]][2],
         n = 0) %>% 
  mutate(composition = gsub("(.{4})", "\\1 ", composition, perl = TRUE))
states_europe <- cbind(states_europe, type = rep(NA, nrow(states_europe)))

states_europe[1, ] <- c(0, "INVALID", "nope", 0, 0, 0, "nope")

# check
# states_europe %>% filter(composition == "miss")
# states_europe %>% filter(stateId == 3058)

# save result
write_csv(states_europe, paste0(path, "/dnn/states_europe.csv"))
