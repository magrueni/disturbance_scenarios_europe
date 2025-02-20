# --------------------------------------------------------------
# Fire Module Scipt 06: Define fire state transitions
# --------------------------------------------------------------
#
# This script defines the transitions for the forest states after a fire
# Essentially, the height and LAI dimension of the forest state is reduced to 0 (1 is the lowest class for LAI)
# However, the new state does not necesseraly occur in the forest state space. 
# Therefore we use a mapping function to find the closest forest state.
# 
# Note: we save the state transitions also for wind and mgmt module, as potential transitions are the same.
#
#-------------------------------------------------------------------------------




### libraries ----------------
library(terra)
library(sf)
library(tidyverse)
library(dplyr)
library(MetBrewer)
library(DBI)


### settings -------------------------------
path <- "/.../"

# Raster package
terraOptions(memfrac=0.5, tempdir = "./tmp/")
tmpFiles(remove=TRUE, current=T, orphan=T)


### load functions --------------------

# load a function which assigns any forest state (from the millions possible combinations) to a state that occurs in the training data.
source(paste0("./functions/states_mapping_function.R"))



### load data -----------------------

dnn_lookup <- read_csv(paste0(path, "/dnn/dnn_states_lookup.csv"))
states <- as.data.frame(dnn_lookup) %>% 
  mutate(stateID = as.numeric(stateID))
dnn_lookup <- as.data.frame(rbind(c("non_state", -1), dnn_lookup)) %>% 
  mutate(stateID = as.numeric(stateID))

### define for each state the state transition after fire --------------------------------------------------------------

# but for the moment, I just create a lookup that sets the stand to the same veg, LAI = 1 and height = 0_2
state_transition_mat <- states %>% 
  rename("state_old" = "state") %>% 
  rowwise() %>% 
  mutate(state = strsplit(state_old, "_")[[1]][1]) %>% 
  mutate(state = paste0(state, "_1_0_2")) %>% 
  mutate(state = ifelse(state == "missing_1_0_2", "missing", state)) %>% 
  mutate(state = states_mapping_function(state, dnn_lookup, max_height_diff_allowed = 4)) %>% 
  left_join(., states, by = c("state")) %>% 
  mutate(key = 1,
         p = 1) %>% 
  dplyr::select(stateID.x, key, stateID.y, p)

colnames(state_transition_mat) <- c("stateId", "key", "targetId", "p")
state_transition_mat <- state_transition_mat %>% 
  mutate(targetId = ifelse(is.na(targetId), 0, targetId))



state_transition_mat2 <- state_transition_mat %>% 
  mutate(key = 0)
state_transition_mat <- rbind(state_transition_mat, state_transition_mat2)



# write them out for the other dist agents.
write.csv(state_transition_mat, paste0(path, "/fire_module/fire_transitions.csv"), row.names = F)


# same transitions are used for the wind module and the mgmt module
write.csv(state_transition_mat, paste0(path, "/mgmt_module/mgmt_transitions.csv"), row.names = F)
write.csv(state_transition_mat, paste0(path, "/wind_module/wind_transitions.csv"), row.names = F)





### states burn probability ------------------------------------------

# we can define the burn probability. In this case all states are the same
firestates <- state_transition_mat %>% 
  dplyr::select(stateId) %>% 
  mutate(pSeverity = 1,
         pBurn = 0.75)

write.csv(firestates, paste0(path, "/fire_module/firestates.csv"), row.names = F)


## end ---------------------------------------------------------------------