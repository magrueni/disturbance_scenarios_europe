# --------------------------------------------------------------
# Initial landscape Scipt 05: Evaluation of the initial landscape
# --------------------------------------------------------------
# Author: Marc Grünig
# Date: 2025-01-07
# 
#
# Description:
# This script analyzes and compares estimated vegetation data against real-world 
# datasets and other reference datasets.
# It specifically evaluates vegetation patterns, species distributions, 
# and forest composition across countries, with a focus on Germany and European regions.
# 
# In specific, we first spatially evaluate the distribution of the tree species
# with a tree distribution map from Blickensdörfer et al.
# In a second step, we look at the share of broadleaved,
# coniferous and mixed forests for all countries in the study area.
#
#
# Notes:
# - We can't share the evaluation dataset.
# but see https://www.sciencedirect.com/science/article/pii/S0034425724000804#da0005
#
# --------------------------------------------------------------



### Libraries ------------------------------------------------------------------

library(terra)
library(dplyr)
library(purrr)
library(data.table)
library(stringr)
library(parallel)
library(tidyr)
library(ggplot2)
library(readr)


### load data ------------------------------------------------------------------

# paths
path <-"/.../"


# dnn states lookup table
dnn_lookup <- read_csv(paste0(path, "/dnn/dnn_states_lookup.csv"))
dnn_lookup <- dnn_lookup %>% mutate(stateID = as.numeric(stateID))

# load init veg landscape 
init_veg <- rast(paste0(path, "/03_initial_forest_states/init_veg_v2.tif"))


### compare to Blickensdorfer et al., 2023 data from Remote sensing for Germany ----------------------------------------

ger <- rast(paste0(path, "/species_class_sum_INT1U.tif")) # Blickensdörfer map
ger_agg <- aggregate(ger, fact = 10, fun = "modal")

# ID,species
# 2,Birch
# 3,Beech
# 4,Douglas fir
# 5,Oak
# 6,Alder
# 8,Spruce
# 9,Pine
# 10,Larch
# 14,Fir
# 16,ODH
# 17,ODL

ger_shp <- vect(paste0(path, "/06_gis/countries/germany.shp"))

sp <- c("bepe", "fasy", "psme", "quro", "algl", "piab", "pisy", "lade", "abal")
sp_name <- c("Quercus robur", "Fagus")
num <- c(2, 3, 4, 5, 6, 8, 9, 10, 14)
area_map <- rep(NA, length(sp))
area_calc <- rep(NA, length(sp))
sp_df <- as.data.frame(cbind(sp, num, area_map, area_calc))

# load init veg landscape 
origin(init_veg) <- origin(ger_agg)
init_veg_crop <- terra::crop(init_veg, ger_agg)
init_veg_crop <- terra::mask(init_veg_crop, ger_agg)


for(i in 1:nrow(sp_df)){
  
  
  # isolate species from species map
  ger_sp <- ger_agg
  
  sp_num <- as.numeric(sp_df[i, "num"])
  sp_code <- sp_df[i, "sp"]
  
  sp_name <- sp_df[i, "sp"]
  
  ger_sp[ger_sp != sp_num] <- NA
  ger_sp[ger_sp == sp_num] <- 1
  
  # get the number of ha
  sp_df[i, "area_map"] <- sum(values(ger_sp), na.rm = T)
  
  # get states for dominant and other occuring
  dom_states <- dnn_lookup %>% filter(grepl(toupper(sp_code), state)) %>% 
    dplyr::select(stateID) %>% as.vector() %>% unlist()
  other <- dnn_lookup %>% filter(grepl(sp_code, state)) %>% 
    dplyr::select(stateID) %>% as.vector() %>% unlist()
  
  dom_states <- unlist(as.vector(dom_states))
  other <- unlist(as.vector(other))
  
  # dominant cells get 1
  dom <- init_veg_crop
  dom[!dom %in% dom_states] <- 0
  dom[dom != 0] <- 2
  
  # others get 0.2
  oth <- init_veg_crop
  oth[!oth %in% other] <- 0
  oth[oth %in% other] <- 1
  oth[oth != 0] <- 1
  
  # put together
  veg <- dom + oth
  veg[veg == 0] <- NA
  
  # get number of estimated ha
  sp_df[i, "area_calc"] <- (sum(values(dom), na.rm = T) * 0.8)  +  (sum(values(oth), na.rm = T) * 0.2)
  
  # plot
  
  veg[veg > 0] <- 1
  
  
  dev.off()
  par(mfrow =  c(1,2))
  
  
  
  #draw germany around
  plot(ger_shp, main = paste0("Blickensdöfer et al. 2024"), axes = F, box = F)
  plot(ger_sp, col = "red", add = T, axes = FALSE, box = F, legend = F)
  
  plot(ger_shp, main = paste0("Initial forest landscape"), axes = F, box = F)
  plot(veg, col = "blue", add = T, axes = FALSE, box = F, legend = F)
  
  dev.print(png, paste0(path, "/03_initial_forest_states/evaluation/", sp_code, "_germany.png"), res = 300, width = 2500)
  
}

print(sp_df)
write.csv(sp_df, paste0(path, "/03_initial_forest_states/evaluation/species_area_germany.csv"))



### country level number of species: mixed vs. pure stands -----------------------------------------------------------

# calculate per country the percentage of area with only 1 species vs the percentage with multiple species
# and also per continent
vals <- na.omit(as.data.frame(init_veg))

# convert vals to data.table
dt_vals <- data.table(vals)
colnames(dt_vals) <- "stateID"

# join dt_vals with dnn_lookup on stateID
dnn_lookup <- dnn_lookup %>% 
  rowwise() %>% 
  mutate(species = str_split(state, "_")[[1]][1]) %>% 
  mutate(number = ifelse(nchar(species) > 4, 2, 1))

dt_result <- dt_vals %>% left_join(dnn_lookup, by = c("stateID"))

nrow(dt_result[dt_result$number == 1, ])/ nrow(dt_result)


# mask for some countries
countries <- list.files(paste0(path, "/06_gis/countries"), pattern = ".shp", full.names = TRUE)
countries <- countries[!grepl("europe", countries)]

countries_check <- c("norway", "sweden", "finland", "ireland", "portugal",
                     "spain", "france", "belgium", "netherlands", "denmark", 
                     "switzerland", "austria", "estonia", "latvia", "lithuania",
                     "poland", "czechia", "slovakia", "croatia")

countries <- countries[grepl(paste(countries_check, collapse = "|"), countries)]


# forest2000 data on percentage of mixed forests
pcts <- as.numeric(c(45, 10, 30, 25, 50,
                     15, 20, 25, 10, 30, 
                     10, 45, 10, 40, 10, 
                     45, 10, 10, 15))
pcts_calc <- rep(NA, length(countries))
output_df <- as.data.frame(cbind(countries, pcts, pcts_calc)) 



for(c in 1:length(countries_check)) {
  
  country <- countries[grepl(paste(countries_check[c]), countries)]
  country <- substr(country, 64, nchar(country)-4)
  country_shp <- vect(countries[grepl(paste(countries_check[c]), countries)])
  
  veg_country <- terra::crop(init_veg, country_shp)
  veg_country <- terra::mask(veg_country, country_shp)
  
  vals <- na.omit(as.data.frame(veg_country))
  dt_vals <- data.table(vals)
  colnames(dt_vals) <- "stateID"
  
  dt_result <- dt_vals %>% left_join(dnn_lookup, by = c("stateID"))
  
  output_df[c, "pcts_calc"] <- round(nrow(dt_result[dt_result$number == 1, ])/ nrow(dt_result), 2) * 100
  output_df[c, "countries"] <- country
  
}

df <- output_df %>% mutate(pcts_calc = round(as.numeric(pcts_calc), 2))
write.csv(df, paste0(path, "/03_initial_forest_states/evaluation/mixed_vs_pure_country_level.csv"))
df <- read.csv(paste0(path, "/03_initial_forest_states/evaluation/mixed_vs_pure_country_level.csv"))

df <- df %>% dplyr::select(countries, pcts, pcts_calc) %>% 
  mutate(pcts = as.numeric(pcts)) %>% 
  dplyr::rename("real_dat" = pcts,
                "estimated" = pcts_calc) %>% 
  pivot_longer(cols = c(real_dat, estimated), names_to = "dat", values_to = "value")



p <- ggplot(df, aes(x = countries, y = value, fill = dat)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "% mixed species stands") +
  xlab("Countries") +
  ylab("%") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw() + 
  theme(axis.text = element_text(size = 16),
        title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave(p, filename = paste0(path, "/03_initial_forest_states/evaluation/mixed_species_stand_per_country.png"), width = 8, height = 6)



### cconiferous vs decidous -------------------------------------------------------------------------------------------


# load copernicus forest mask data
forest_mask <- rast(paste0(path, "/06_gis/forest_mask_new_laea_masked.tif"))
forest_mask <- aggregate(forest_mask, fact = 10, fun = "modal", na.rm = T)
init_veg <- rast(paste0(path, "/03_initial_forest_states/init_veg_v2.tif"))
init_veg <- aggregate(init_veg, fact = 10, fun = "modal", na.rm = T)

forest_crop <- terra::crop(extend(forest_mask, init_veg), init_veg)
forest_crop <- terra::project(forest_crop, init_veg)
forest_crop <- terra::mask(forest_crop, init_veg)

rm(forest_mask)

forest_broa <- forest_crop
forest_broa[forest_broa != 1] <- NA

forest_coni <- forest_crop
forest_coni[forest_coni != 2] <- NA

forest_mix <- forest_crop
forest_mix[forest_mix != 3] <- NA

broa_sum <- sum(values(forest_broa), na.rm = T)
coni_sum <- sum(values(forest_coni), na.rm = T)
mix_sum <- sum(values(forest_mix), na.rm = T)

rm(forest_crop)


# classification of the initial vegetation layer to those three classes

broads <- c("fasy", "bepe", "quro", "quil", "acpl", "acps", "frex", "potr", "qupe", "saat", "soar", "tico", "tipl",
            "ulgl", "acca", "acmo", "juco", "ulmi", "rops", "prav", "phla", "qupu", "qusu", "coav")
conis <- c("abal", "lade", "piab", "piun", "pini", "pipi", "piha")

dnn_lookup_new <- dnn_lookup %>% 
  mutate(stateID = as.numeric(stateID)) %>% 
  rowwise() %>% 
  mutate(species = str_split(state, "_")[[1]][1]) %>% 
  mutate(type = if_else(grepl(paste0(c(conis, toupper(conis)), collapse = "|"), species) & grepl(paste0(c(broads, toupper(broads)), collapse = "|"), species), 3,
                        if_else(grepl(paste0(c(conis, toupper(conis)), collapse = "|"), species) & !grepl(paste0(c(broads, toupper(broads)), collapse = "|"), species), 2, 1)))


# check which countries are not well represented to get an idea about the regions
# mask for some countries
countries <- list.files(paste0(path, "/06_gis/countries"), pattern = ".shp", full.names = TRUE)
countries <- countries[!grepl("europe", countries)]

# forest2000 data.
output_df <- list()
countries_names <- c()

for(c in 1:length(countries)) {
  
  
  country <- substr(countries[c], 58, nchar(countries[c])-4)
  print(country)
  
  country_shp <- vect(countries[c])
  countries_names <- c(countries_names, country)
  
  
  # mask the rasters with the country
  forest_broa_c <- terra::mask(forest_broa, country_shp)
  forest_coni_c <- terra::mask(forest_coni, country_shp)
  forest_mix_c <- terra::mask(forest_mix, country_shp)
  
  broa_sum_c <- sum(values(forest_broa_c), na.rm = T)
  coni_sum_c <- sum(values(forest_coni_c), na.rm = T)
  mix_sum_c <- sum(values(forest_mix_c), na.rm = T)
  sum_cop <- broa_sum_c + coni_sum_c + mix_sum_c
  
  rm(forest_broa_c, forest_coni_c, forest_mix_c)
  gc()
  
  # same for estimated vegetation
  init_veg_c <- terra::mask(init_veg, country_shp)
  init_veg_df <- as.data.frame(init_veg_c, xy = T)
  colnames(init_veg_df) <- c("x", "y", "stateID")
  init_veg_df_class <- init_veg_df %>% 
    left_join(dnn_lookup_new, by = "stateID")
  
  
  broa_sum_est_c <- nrow(init_veg_df_class[init_veg_df_class$type == 1,])
  coni_sum_est_c <- nrow(init_veg_df_class[init_veg_df_class$type == 2,])
  mix_sum_est_c <- nrow(init_veg_df_class[init_veg_df_class$type == 3,])
  sum_est <- broa_sum_est_c + coni_sum_est_c + mix_sum_est_c
  
  rm(init_veg_c, init_veg_df)
  
  out <- c(broa_sum_c, coni_sum_c, mix_sum_c, sum_cop, broa_sum_est_c, coni_sum_est_c, mix_sum_est_c, sum_est)
  
  output_df[[c]] <- out
  
  
  gc()
  
}

gc()
tmpFiles(remove = T)


# check the values
init_veg_df <- as.data.frame(init_veg, xy = T)
colnames(init_veg_df) <- c("x", "y", "stateID")
init_veg_df_class <- init_veg_df %>% 
  left_join(dnn_lookup_new, by = "stateID")

broa_sum_est <- nrow(init_veg_df_class[init_veg_df_class$type == 1,])
coni_sum_est <- nrow(init_veg_df_class[init_veg_df_class$type == 2,])
mix_sum_est <- nrow(init_veg_df_class[init_veg_df_class$type == 3,])
sum_est <- broa_sum_est + coni_sum_est + mix_sum_est
sum_cop <- broa_sum + coni_sum + mix_sum


output_df[[length(countries) + 1]] <- c(broa_sum, coni_sum, mix_sum, sum_cop, broa_sum_est, coni_sum_est, mix_sum_est, sum_est)
output_df <- do.call(rbind, output_df)
rownames(output_df) <- c(countries_names, "eu")
colnames(output_df) <- c("broa_sum_c", "coni_sum_c", "mix_sum_c", "sum_cop",
                         "broa_sum_est_c", "coni_sum_est_c", "mix_sum_est_c", "sum_est")

write.csv(output_df, paste0(path, "/03_initial_forest_states/evaluation/forest_type_comparison.csv"))


# plotting
df <- read.csv(paste0(path, "/03_initial_forest_states/evaluation/forest_type_comparison.csv"))


df_mixed <- df %>% rename("countries" = "X") %>% dplyr::select(countries, mix_sum_c , mix_sum_est_c, sum_cop, sum_est ) %>% 
  mutate(real_dat = 100/sum_cop * mix_sum_c,
         estimated = 100/sum_est * mix_sum_est_c ) %>% 
  pivot_longer(cols = c(real_dat, estimated), names_to = "dat", values_to = "value") %>% 
  mutate(value = value)


p_mixed <- ggplot(df_mixed %>% filter(countries != "eu"), aes(x = countries, y = value, fill = dat)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Area mixed species stands") +
  xlab("Countries") +
  ylab("area [ha]") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw() + 
  theme(axis.text = element_text(size = 16),
        title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_mixed
ggsave(p_mixed, filename = paste0(path, "/03_initial_forest_states/evaluation/mixed_species_stand_per_country_copernicus.png"), width = 8, height = 6)

countries <- df_mixed %>% filter(countries != "eu")
cor(countries$mix_sum_c, countries$mix_sum_est_c)


df <- read.csv(paste0(path, "/03_initial_forest_states/evaluation/forest_type_comparison.csv"))
df_broa <- df %>% rename("countries" = "X") %>% dplyr::select(countries, broa_sum_c, broa_sum_est_c , sum_cop, sum_est  ) %>% 
  mutate(real_dat = 100/sum_cop * broa_sum_c,
         estimated = 100/sum_est * broa_sum_est_c ) %>% 
  pivot_longer(cols = c(real_dat, estimated), names_to = "dat", values_to = "value") 

p_broa <- ggplot(df_broa %>% filter(countries != "eu"), aes(x = countries, y = value, fill = dat)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Area broadleaved species stands") +
  xlab("Countries") +
  ylab("area [ha]") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw() + 
  theme(axis.text = element_text(size = 16),
        title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_broa
ggsave(p_broa, filename = paste0(path, "/03_initial_forest_states/evaluation/broadleaved_species_stand_per_country_copernicus.png"), width = 8, height = 6)

countries <- df_broa %>% filter(countries != "eu")
cor(countries$broa_sum_c, countries$broa_sum_est_c)

df <- read.csv(paste0(path, "/03_initial_forest_states/evaluation/forest_type_comparison.csv"))
df_coni <- df %>% rename("countries" = "X") %>% dplyr::select(countries, coni_sum_c , coni_sum_est_c, sum_cop, sum_est  ) %>% 
  mutate(real_dat = 100/sum_cop * coni_sum_c,
         estimated = 100/sum_est * coni_sum_est_c ) %>% 
  pivot_longer(cols = c(real_dat, estimated), names_to = "dat", values_to = "value")

p_broa <- ggplot(df_coni %>% filter(countries != "eu"), aes(x = countries, y = value, fill = dat)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Area coniferous species stands") +
  xlab("Countries") +
  ylab("area [ha]") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw() + 
  theme(axis.text = element_text(size = 16),
        title = element_text(size = 18),
        legend.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1))
p_broa
ggsave(p_broa, filename = paste0(path, "/03_initial_forest_states/evaluation/coniferous_species_stand_per_country_copernicus.png"), width = 8, height = 6)

countries <- df_coni %>% filter(countries != "eu")
cor(countries$coni_sum_c, countries$coni_sum_est_c)


### put together
df <- df_mixed %>% dplyr::select(countries, dat, value_mixed = value) %>% 
  left_join(., df_coni %>% 
              dplyr::select(countries, dat, value_coni = value), 
            by = c("countries", "dat")) %>% 
  left_join(., df_broa %>% 
              dplyr::select(countries, dat, value_broa = value), 
            by = c("countries", "dat")) %>% 
  mutate(sum = value_mixed + value_coni + value_broa)

head(df)
df <- df %>% filter(!countries %in% c("moldova", "andorra"))


# Pivot the data from wide to long format
df_long <- tidyr::pivot_longer(df, cols = c(value_mixed, value_coni, value_broa), names_to = "variable", values_to = "value")
df_long <- df_long %>% rowwise() %>% 
  mutate(countries = str_to_title(countries))
df_long$countries <- forcats::fct_rev(df_long$countries)
df_long <- df_long %>% rowwise() %>% 
  mutate(variable = ifelse(variable == "value_mixed", "mixed",
                           ifelse(variable == "value_coni", "coniferous",
                                  "broadleaved"))) %>% 
  mutate(dat = ifelse(dat == "estimated", "Inital forest landscape", "Copernicus data"))

df_long <- df_long %>% filter(countries != "Eu")


# Create the grouped stacked bar plot
library(MetBrewer)
copern <- ggplot(df_long, aes(x = countries, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.7) +
  facet_grid(.~dat, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = met.brewer("Isfahan1", 8)[c(1,3, 6)]) +
  labs(x = "Countries", y = "Proportion", fill = "Forest type") +
  theme_minimal() +
  coord_flip()

copern
ggsave(copern,
       filename = paste0(path, "/03_initial_forest_states/evaluation/species_stand_per_country_copernicus_all.png"),
       width = 6, height = 6)


### end ------------------------------------------------------------------------