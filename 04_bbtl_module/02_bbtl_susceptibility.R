# --------------------------------------------------------------
# Bark Beetle Module Scipt 02: Forest state susceptibility to BB
# --------------------------------------------------------------
# Author: Marc Grünig
# Date: 2025-01-07
# 
# Description:
# For each forest state, we want to assign the susceptibility to BB
# damage based on height and spruce share. Further, we create a transition matrix
# that defines the forest state after disturbance for each forest state.
# Here we remove the Spruce from mixed stand instead of only setting the height and LAI to 0
# 
# Notes:
# - The formulas are from Neherer et al. (Acute Drought Is an Important Driver of Bark Beetle Infestation in Austrian Norway Spruce Stands)
# and Seidl et al. 2005
# --------------------------------------------------------------


### libraries --------
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyverse)
library(DBI)
library(terra)


# define folder path
path <- "/.../"



### define the susceptibility based on stand age and spruce share

# -	Spruce share: bas = 1-0.92*exp(-4.9*x^2.4); x = spruce share
# -	Age: 1-0.8*exp(-4.2*(x/100)^3.3) x = stand age 
# -	Ks (Kronenschluss):  if(x>70, 0.5, 1-0.021*x-0.0008*x^2+3*10^-5*x^3-2.25*10^-7*x^4) -> we can omit
# -	Dd: 1-exp(-150*smi^2.5)
# -	0.3*bas + 0.25*age + 0.15*ks + 0.3*dd
# -	Without ks it will be: 0.35*bas + 0.3*age + 0.35*dd

spruce_share <- 0.1
stand_age <- 40

bas <- 1-0.92*exp(-4.9*spruce_share^2.4)
age <- 1-0.8*exp(-4.2*(stand_age/100)^3.3)
smi <- 0.15
dd <- 1-exp(-150*smi^2.5) # smi =?

susceptibility <- 0.35*bas + 0.3*age + 0.35*dd
susceptibility

susc_bb <- function(x, smi){
  
  spruce_share <- as.numeric(x["spruce_share"])
  if(spruce_share == 0){
    return(0)
  }else{
    stand_age <- as.numeric(x["age"])
    bas <- 1-0.92*exp(-4.9*spruce_share^2.4)
    age <- 1-0.8*exp(-4.2*(stand_age/100)^3.3)
    dd <- 1-exp(-150*smi^2.5) # smi =?
    
    susceptibility <- 0.35*bas + 0.3*age + 0.35*dd
    return(susceptibility)
  }
}


x <- seq(1, 150, 1)
plot(1-0.8*exp(-4.2*(x/100)^3.3) )



# get yield tables to calibrate a model that predicts age from height
# we cannot share this data, the calibrated model is saved in the folder

# Ertragstafeln <- read_excel(paste0(path, "/mgmt_module/Ertragstafeln.xls"))
# tafel <- Ertragstafeln %>% filter(ERTRAGSTAFEL_KURZ=="FI13")
# tafel <- tafel %>% group_by(ET_TABELLE_BONITAET) %>%
#   mutate(delta_h = (ET_TAB_OBERHOEHE - lag(ET_TAB_OBERHOEHE))/((ET_TABELLE_ALTER - lag(ET_TABELLE_ALTER))) ) %>% 
#   mutate(r_years = 2/delta_h)
# 
# tafel <- tafel %>%
#   mutate(h_cls = round(ET_TAB_OBERHOEHE / 2)*2) %>%
#   filter(h_cls > 14) %>% 
#   filter(ET_TABELLE_BONITAET > 5 & ET_TABELLE_BONITAET < 10)
# 
# ggplot(tafel, aes(y = ET_TABELLE_ALTER, x = ET_TAB_OBERHOEHE, col = ET_TABELLE_BONITAET)) +
#   geom_point() +
#   geom_smooth(method = lm, formula = y ~ splines::bs(x, 2), se = FALSE)
# 
# 
# tafel <- tafel[, c("ET_TABELLE_ALTER", "ET_TAB_OBERHOEHE")]
# colnames(tafel) <- c("age", "height")
# mod <- glm(age ~ poly(height, 2), data = tafel)
# newdat <- as.data.frame(seq(1, 50, 1))
# colnames(newdat) <- "height"
# predict(mod, newdata = newdat)
# 
# saveRDS(mod, paste0(path, "/bbtl_module/age_estimator.rds"))




### do the matrix for the different states!
# now that we get the correct information on the probability damage, we can load the svd state matrix and calculate the susceptibility for each state!

# load states lookup table
states <- read_csv(paste0(path, "/dnn/dnn_states_lookup.csv"))

# state matrix
states_meta <- read_csv(paste0(path, '/dnn/states_europe_meta.csv'))
states_meta_spruce <- states_meta %>% dplyr::select(stateId, unique.str, piab) %>% 
  rename("state" = "unique.str")

# load the age estimation model
mod <- readRDS(paste0(path, "/bbtl_module/age_estimator.rds"))

# lookup for each 
states_check <- states %>% 
  rename("stateId" = "stateID") %>% 
  left_join(., states_meta_spruce, by = c("state", "stateId")) %>% 
  rename("spruce_share" = "piab") %>% 
  rowwise() %>% 
  mutate(height = as.numeric(strsplit(state, "_")[[1]][3])) %>% 
  mutate(height = ifelse(height < 0, 0, height))
states_check$age <- predict(mod, states_check)
states_check[1, c("spruce_share", "height", "age")] <- 0 

states_check <- states_check %>% mutate(age = ifelse(age < 0, 0, age))

# suceptbility for trees smaller than 10m is set to 0
states_check$bb_susceptibility <- apply(states_check , 1, susc_bb, smi = 2.5)
states_check <- states_check %>% mutate(bb_susceptibility = ifelse(height < 10, 0, bb_susceptibility))
hist(states_check$bb_susceptibility)

# visualize
ggplot(states_check %>% 
         filter(grepl("PIAB|piab", state)) %>%
         mutate(bb_susceptibility = ifelse(bb_susceptibility > 0.9, 0.9, bb_susceptibility)), aes(y = bb_susceptibility, x = height)) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x), se = FALSE)

# rename to SVD format
states_check <- states_check %>%
  dplyr::select(state, stateId, pBarkBeetleDamage = bb_susceptibility)

# damage probability larger than 0.9 is set to 0.9
states_check <- states_check %>% mutate(pBarkBeetleDamage = ifelse(pBarkBeetleDamage > 0.9, 0.9, pBarkBeetleDamage))
states_check %>% filter(grepl("PIAB", state)) %>% View()

write_csv(states_check, paste0(path, "/bbtl_module/states_bb.csv"))    


### define transition matrix after bb damage
# states with PIAB only are reduced in height and LAI.
# mixed states only loose piab share and are reduced in LAI

### define for each state the state transition after fire --------------------------------------------------------------
source("./functions/states_mapping_function.R")

# get the lookup for the states and extract PIAB states
states <- read_csv(paste0(path, "/dnn/dnn_states_lookup.csv"))
dnn_lookup <- rbind(c("non_state", -1), states)
dnn_lookup_piab <- dnn_lookup %>% filter(grepl("PIAB", state)) 
dnn_lookup_no_piab <- dnn_lookup %>% filter(!grepl("piab|PIAB", state)) %>% bind_rows(., dnn_lookup_piab[1,]) %>% mutate(stateID = as.numeric(stateID)) %>% as.data.frame()

# but for the moment, just create a lookup that sets the stand to the same veg, LAI = 1 and height = 0_2
state_transition_mat <- states %>% 
  rename("state_old" = "state") %>% 
  rowwise() %>% 
  mutate(sp = strsplit(state_old, "_")[[1]][1]) %>% 
  mutate(lai = strsplit(state_old, "_")[[1]][2]) %>% 
  mutate(height = gsub(paste0(sp, "_", lai, "_"), "", state_old)) %>% 
  mutate(sp_new = gsub("PIAB|piab", "", sp)) %>% 
  mutate(sp_new = ifelse(sp_new == "", sp, sp_new)) %>% 
  mutate(lai_new = as.numeric(ifelse(grepl("PIAB|piab", sp_new), "1", lai))) %>%
  mutate(height_new = ifelse(grepl("PIAB|piab", sp_new), paste0("0_2"), height)) %>% 
  mutate(state_new = paste0(sp_new, "_", lai_new, "_", height_new)) %>% 
  mutate(state_new = ifelse(grepl("missing", state_new), "missing_1_0_2", state_new)) %>% 
  filter(state_new != "missing_1_0_2") %>% 
  rowwise() %>% 
  mutate(state = states_mapping_function(state_new, dnn_lookup_no_piab, max_height_diff_allowed = 4)) %>% 
  left_join(., states, by = c("state")) %>% 
  mutate(key = 1, p = 1) 

state_transition_mat %>% filter(grepl("PIAB", state_old)) %>% View()


state_transition_mat <- state_transition_mat %>% 
  dplyr::select(stateID.x, key, stateID.y, p)

state_transition_mat <- rbind(state_transition_mat, c(0,1,0,1))
colnames(state_transition_mat) <- c("stateId", "key", "targetId", "p")
state_transition_mat %>% filter(is.na(targetId))

state_transition_mat <- state_transition_mat %>% 
  mutate(targetId = ifelse(is.na(targetId), as.numeric(dnn_lookup_piab[1,2]), targetId))
state_transition_mat2 <- state_transition_mat %>% 
  mutate(key = 0)
state_transition_mat <- rbind(state_transition_mat, state_transition_mat2)

write_csv(state_transition_mat, paste0(path, "/bbtl_module/transition_matrix_bb.csv"))    


### end ---------------
