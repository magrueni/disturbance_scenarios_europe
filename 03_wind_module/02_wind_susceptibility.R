# --------------------------------------------------------
# Script to assign wind susceptibility to forest states
#
#
# This script calculates the probability of wind damage 
# for different tree species based on height and diameter.
# It applies susceptibility coefficients from Schmidt et al. 2010 
# and processes forest states to assign wind susceptibility scores.
# --------------------------------------------------------



### libraries --------------------------------
library(ggplot2)
library(dplyr)

path <- "/.../"

### load functions ---------------------------

# We use the coefficients from Schmidt et al. 2010
coefficents <- as.data.frame(rbind(
  c('abal', -8.46, -0.505, -2.449 ),
  c('acps', -8.78, -0.287, -1.770 ),
  c('alin', -8.78, -0.287, -1.770 ),
  c('alvi', -8.78, -0.287, -1.770 ),
  c('bepe', -8.78, -0.287, -1.770 ),
  c('cabe', -8.78, -0.287, -1.770 ),
  c('fasy', -13.04,  -0.998, -3.994 ),
  c('lade', -8.59, -1.625,  -3.525 ),
  c('piab', -12.27,  -1.775, -5.128 ),
  c('pice', -8.59, -1.625, -3.525 ),
  c('pisy', -8.59, -1.625, -3.525 ),
  c('poni', -8.78, -0.287,  -1.770 ),
  c('potr', -8.78, -0.287,  -1.770 ),
  c('qupe', -13.04, -0.998,  -3.994 ),
  c('tico', -8.78, -0.287,  -1.770 )))

colnames(coefficents) <- c("species", "intercept", "dbh", "height")
coefficents <- coefficents %>% 
  mutate(intercept = as.numeric(intercept),
         dbh = as.numeric(dbh),
         height = as.numeric(height))


# susceptibility functions according to iLand script:
fun1 <- function(b1, a1, y1, DBH, height){
  g <- b1 + log(DBH^(a1)/height^(y1)) + 3
  return(g)
}

fun2 <- function(g){
  p <-  exp(g)/(exp(g) + 1)
  return(p)
}


# function to get the probability based on the height and coefficients
proba_fun <- function(species, height, coefficents){
  
  DBH <- 100*height/80
  
  if(height == 0){p <- 0}else{
    
    b1 <- coefficents[coefficents$species == species, "intercept"]
    a1 <- coefficents[coefficents$species == species, "dbh"]
    y1 <- coefficents[coefficents$species == species, "height"]
    
    
    g <- fun1(b1, a1, y1, DBH, height)
    p <- fun2(g)
  }
  
  return(p)
}



# demonstrate how susceptibility changes with different heights
heights <- seq(1, 50, 1)
all_sp <- list()

# for(i in 1:length(coefficents$species)){
for(i in c(9)){ # example for FASY
  
  probas_sp <- c()
  
  for(h in 1:length(heights)){
    
    height <- heights[h]
    DBH <- 100*height/80
    
    b1 <- coefficents[i, "intercept"]
    a1 <- coefficents[i, "dbh"]
    y1 <- coefficents[i, "height"]
    
    
    g <- fun1(b1, a1, y1, DBH, height)
    p <- fun2(g)
    
    probas_sp <- c(probas_sp, p)
    
  }
  
  
  df <- as.data.frame(cbind(heights, probas_sp, species = paste(coefficents[i, "species"])))
  all_sp[[i]] <- df
  
  ggplot(df, aes(x = heights, y = probas_sp)) +
    geom_point() + 
    geom_smooth(method = "gam")
  
}

all_sp <- do.call(rbind, all_sp)
all_sp <- all_sp %>% 
  mutate(heights = as.numeric(heights),
         probas_sp = as.numeric(probas_sp))

ggplot(all_sp, aes(x = heights, y = probas_sp, colour = species)) +
  geom_point() + 
  geom_smooth(method = "gam") +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) + 
  theme_bw() 




### processing -------------------------------------------------------

# load states lookup table
dnn_lookup <- read_csv(paste0(path, "/dnn/dnn_states_lookup.csv"))
states <- as.data.frame(dnn_lookup) %>% 
  mutate(stateID = as.numeric(stateID))
dnn_lookup <- as.data.frame(rbind(c("non_state", -1), dnn_lookup)) %>% 
  mutate(stateID = as.numeric(stateID))


#let's check for each state if there is a species in the wind susceptibility data
coefficents <- as.data.frame(rbind(
  c('fasypiab', -12.655, -1.3865, -4.561),
  c('abalfasy', -10.75, -0.7515, -3.2215),
  c('abal', -8.46, -0.505, -2.449 ),
  c('acps', -8.78, -0.287, -1.770 ),
  c('alin', -8.78, -0.287, -1.770 ),
  c('alvi', -8.78, -0.287, -1.770 ),
  c('bepe', -8.78, -0.287, -1.770 ),
  c('cabe', -8.78, -0.287, -1.770 ),
  c('fasy', -13.04,  -0.998, -3.994 ),
  c('lade', -8.59, -1.625,  -3.525 ),
  c('piab', -12.27,  -1.775, -5.128 ),
  c('pice', -8.59, -1.625, -3.525 ),
  c('pisy', -8.59, -1.625, -3.525 ),
  c('poni', -8.78, -0.287,  -1.770 ),
  c('potr', -8.78, -0.287,  -1.770 ),
  c('qupe', -13.04, -0.998,  -3.994 ),
  c('tico', -8.78, -0.287,  -1.770 ),
  c('other_broa', -8.78, -0.287, -1.770 ),
  c('other_coni', -8.59, -1.625, -3.525 ),
  c('other_oaks', -13.04, -0.998,  -3.994)))

colnames(coefficents) <- c("species", "intercept", "dbh", "height")
coefficents <- coefficents %>% 
  mutate(intercept = as.numeric(intercept),
         dbh = as.numeric(dbh),
         height = as.numeric(height))


other_broa <- c("acca", "acmo", "acne", "acop", "algl", "amov", "arun", 
                "buse", "casa", "ceau", "cesi", "coav", "cosa", "crmo", 
                "ecgl", "fica", "frex", "ilaq",  "lano", "masy", "moal",
                "oleu", "phla", "plhi", "prav",  "pypy", "rhal", "bepu",
                "rops", "saat", "sani", "soar", "tipl", "ulgl", "ulmi", "chhu", "taba", "acpl", "myco")

other_coni <- c("celi", "cuse", "juco", "juox", "juph", "jure", "juth", 
                "piha", "pini", "pipi", "pipn", "pira", "piun", "psme")

other_oaks <- c("quca", "qufa",
                "qufr", "quil", "qupu", "qupy", "quro", "qusu")




# extract the dominant species, if not available assign to "other"
states_check <- states %>% 
  rowwise() %>% 
  mutate(sp = strsplit(state, "_")[[1]][1]) %>% 
  mutate(sp_dom = ifelse(str_detect(sp, "[[:upper:]]"), 1, 0))

states_dom <- states_check[states_check$sp_dom == 1,]
states_dom <- states_dom %>% 
  mutate(sp_new = tolower(substring(sp, 0, 4))) %>% 
  mutate(sp_new = ifelse(sp_new %in% other_broa, "other_broa",
                         ifelse(sp_new %in% other_coni, "other_coni",
                                ifelse(sp_new %in% other_oaks, "other_oaks", sp_new))))


# assign species to states without dominant species
states_mix <- states_check[states_check$sp_dom == 0,]
states_mix <- states_mix[-(1),] # remove "missing" state
states_mix <- states_mix %>% 
  mutate(sp_new = ifelse(sp %in% coefficents$species, sp, substring(sp, 0, 4))) %>% 
  mutate(sp_new = ifelse(sp_new %in% other_broa, "other_broa",
                         ifelse(sp_new %in% other_coni, "other_coni",
                                ifelse(sp_new %in% other_oaks, "other_oaks", sp_new))))





states_new <- rbind(states_dom, states_mix)
states_new <- states_new %>% 
  rowwise() %>% 
  mutate(height = as.numeric(strsplit(state, "_")[[1]][3])) %>% 
  mutate(pDamage = proba_fun(sp_new, height, coefficents))


# test
states_new %>% filter(grepl("PIAB_3", state))


# quick check
ggplot(states_new, aes(x = height, y = pDamage)) + geom_point()


## save output
fin_df <- states_new %>% dplyr::select(stateID, pDamage) %>% dplyr::rename("stateId" = "stateID") %>% arrange(stateId)
fin_df <- rbind(c(0,0), fin_df)
write_csv(fin_df, paste0(path, "/wind_module/results/states_wind.csv"))    


