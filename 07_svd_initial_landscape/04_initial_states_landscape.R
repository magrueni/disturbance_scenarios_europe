# --------------------------------------------------------------
# Initial landscape Scipt 04: Create the initial landscape
# --------------------------------------------------------------
# Author: Marc Grünig
# Date: 2025-01-07
# 
#
# Description:
# This script creates a forest landscape for Europe that should reflect current forest conditions.
# This was done by combining different data sources to reflect conditions in specific SVD states.
# 
#
# Notes:
# - The height data comes from Lang et al. 2023
# - LAI data was obtained through google earth engine from MODIS
# - For species composition, we use 2 different datasets (Brus et al., 2012 and Bonannella et al. 2022)
#   these two datasets complement each other to allow the sampling from tree species leading to a 
#   representative forest landscape for the whole continent of Europe.
# - The datasets created here are quite large, instead of providing the underlying data
#   we provide the end product.
#
# --------------------------------------------------------------


### load libraries -------------------------------------------------------------
library(terra)
library(raster)
library(dplyr)
library(sf)
library(lookup)
library(data.table)
library(stringr)
library(tidyr)
library(data.table)
library(parallel)
library(DBI)
library(tidyverse)


path <- "/.../"


# define used projections
proj_leae <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
proj_wgs <- "+proj=longlat +datum=WGS84 +no_defs"


### Building bounding boxes for the full canopy height dataset -----------------
min.long <- -12
min.lat <- 33
max.long <- 42
max.lat <- 72

seq.long <- seq(min.long, max.long, by = 3)
seq.lat <- seq(min.lat, max.lat, by = 3)

combination.min <- expand.grid(seq.long[-length(seq.long)], seq.lat[-length(seq.lat)])
combination.max <- expand.grid(seq.long[-1], seq.lat[-1])

full.combination <- tibble(min.long = combination.min[,1],
                           max.long = combination.max[,1],
                           min.lat = combination.min[,2],
                           max.lat = combination.max[,2])


full.combination <- as.data.frame(full.combination)

bbox <- full.combination %>%
  mutate(left.coord = paste0(ifelse(min.long < 0, "W", "E"), ifelse(abs(min.long) < 10, "00", "0"), round(abs(min.long), 0)),
         top.coord = paste0(ifelse(max.lat < 0, "S", "N"), round(abs(max.lat), 0) - 3)) # minus 3 wüu canopy height vo unge agschribe isch

bbox


### load data for forest mask, vegetation and lai  -----------------------------

# forest mask - obtained from copernicus landcover data (forest type)
forest_mask <- rast(paste0(path_new, "/06_gis/forest_mask_wgs.tif"))

# lai raster - obtained from MODIS through google earth engine
lai_eu <- rast(paste0(path, "/06_gis/LAI_annual_mean_2020_2021.tif"))
lai_eu <- lai_eu/10
lai_eu <- terra::project(lai_eu, proj_wgs)

# vegetation rasters - obtained from Bonanella et al.2022
veg_sdms <- list.files(paste0(path, "/vegetation_data/"), pattern = "_v0.2.tif")
veg_sdms <- veg_sdms[grepl("_p_", veg_sdms)]

# dominant vegetation - obtained from Brus et al. 2012
dom_veg <- list.files(paste0(path, "/EU_TreeMap_Brus_etal/"),
                      pattern = ".tif")

dom_veg <- dom_veg[!grepl("Dominant", dom_veg)]
dom_veg <- dom_veg[!grepl("_wgs", dom_veg)]

# dom veg names to svd state names dict
classification <- substr(dom_veg, 0, nchar(dom_veg)-4)

sp <- c("abal", "algl", "bepe", "broa", "cabe",
        "casa", "coni", "eugl", "fasy", "frex", 
        "lade", "piab", "pine", "pipi", "pisy",
        "potr", "psme", "oaks", "quro", "rops")

# create a lookup table
lookup_domveg <- data.table(cbind(classification, sp))
lookup("AbiesSpp", lookup_domveg$classification, lookup_domveg$sp)


### load functions -------------------------------------------------------------------------------------------

# translate the numeric code to species code function
translator <- function(x){
  
  p2 <- str_split(x, "_")
  p3 <- na.omit(sapply(p2, lookup, lookup_df$classification, lookup_df$sp))
  output <- paste(unlist(p3), collapse = "")
  
  return(output)
  
}

# string reordering function to bring alphabetical order to species codes
reorder_strings <- function(x){
  
  str <- as.character(x)
  
  # check if there are capital letters
  condition <- str_detect(str,"[[:upper:]]")
  
  sorted_str <- sort(unique(strsplit(str, "(?<=.{4})", perl = TRUE)[[1]]))
  
  # if there are capital letters, bring them to the front, else sort alphabetically
  output <- ifelse(condition, paste0(sorted_str[which(str_detect(sorted_str,"[[:upper:]]"))], 
                                     paste0(sorted_str[which(!str_detect(sorted_str,"[[:upper:]]"))], collapse = ""),
                                     collapse = ""), 
                   paste0(sort(unique(strsplit(str, "(?<=.{4})", perl = TRUE)[[1]])), collapse = ""))
  
  return(output)
  
}


# toupper the first 4 characters 
first4up <- function(x) {
  substr(x, 1, 4) <- toupper(substr(x, 1, 4))
  x
}


# combine sp names for lookup table
combine_sp_names <- function(x){
  
  str <- unique(as.character(unlist(x)))
  output <- paste(str, collapse = "")
  return(output)
  
}


# translate the numeric code to species code function
translator_advanced <- function(x){
  
  x <- as.data.frame(x)
  sp_code <- as.character(x["translated_value",])
  code_nr <- as.numeric(lookup_new[.(sp_code), .(code_nr)])

  return(code_nr)
}


# load cpp sampler
# the sampler samples a number of species based on the diversity (Shannon Index) of the Brus et al data
Rcpp::sourceCpp(paste0(path_new, "/03_initial_forest_states/sampler/vegetation_sampler.cpp"))
nspecies_thresholds <- c(1.5, 2.5, 3.5, 100)


### define lookup tables -------------------------------------------------------

# basically, we define the species and create all possible combinations of states
# and with this create a lookup table into SVD states that we have in the database

# first a lookup species for the brus et al., dataset
sp <- c("abal", "algl", "bepe", "cabe", "coav",
        "oleu", "prav", "quil", "qusu", "piha",
        "pini", "pipn",
        "casa", "coni", "eugl", "fasy", "frex", 
        "lade", "piab", "pine", "pipi", "pisy",
        "potr", "psme", "oaks", "quro", "rops",
        "")

sp2 <- c(toupper(sp))


# second, create a lookup table with all possible combinations of species
sp_all <- unique(c(sp, sp2))
all_sp <- do.call(expand.grid, rep(list(sp_all), 4))
all_sp[1,]
all_sp <- apply(all_sp, 1, combine_sp_names)
all_sp[all_sp == ""] <- NA
all_sp <- paste0(substring(all_sp, 1, 4), tolower(substring(all_sp, 5, 8)), tolower(substring(all_sp, 9, 12)))
all_sp <- unique(all_sp)
all_sp_reord <- sapply(unique((all_sp)), reorder_strings)
all_sp_reord <- na.omit(unique(all_sp_reord))
length(all_sp_reord)

# create matrix with admixed species
lowersp <- do.call(expand.grid, rep(list(sp), 4))
lowersp <- apply(lowersp, 1, combine_sp_names)
sp_all <- unique(c(sp, sp2))
all_sp <- do.call(expand.grid, rep(list(sp), 2))
all_sp <- apply(all_sp, 1, combine_sp_names)
all_sp <- unlist(lapply(all_sp, function(i) {
  paste(sp2, i, sep = "")
}))


# combine to have all possible composition combos
all_sp <- unique(c(all_sp, lowersp))
all_sp[all_sp == ""] <- NA


#all_sp <- paste0(substring(all_sp, 1, 4), tolower(substring(all_sp, 5, 8)), tolower(substring(all_sp, 9, 12)))
all_sp <- unique(all_sp)
all_sp_reord <- sapply(unique((all_sp)), reorder_strings)
all_sp_reord <- na.omit(unique(all_sp_reord))
length(all_sp_reord)


# then define all possible LAIs and heights
lai <- c("_1_", "_2_", "_3_")
heights <- paste0(seq(0, 48, 2), "_", seq(2, 50, 2))
combinations <- expand.grid(lai, heights, stringsAsFactors = FALSE)
result <- apply(combinations, 1, paste0, collapse = "")


# combined lai-height combinations with vegetation combos
combinations_all <- expand.grid(all_sp_reord, result, stringsAsFactors = FALSE)
result_fin <- apply(combinations_all, 1, paste0, collapse = "")

# define lookup. Here we add the missing state 0
lookup <- data.table(as.data.frame(cbind(species_code = result_fin, code_nr = c(1:length(result_fin)))) %>%
                       rbind(c("missing", 0), .))


# check dimensions
dim(unique(lookup))


# translate the states to the SVD states
# load the mapping function and the lookup table that are used for the DNN.
source("./functions/states_mapping_function.R")


# now load the dnn states and lookup from the database
dnn_lookup <- read_csv(paste0(path_new, "/03_initial_forest_states/states_lookup_pruned.csv"))
colnames(dnn_lookup) <- c("state", "stateID")


# add missing state and non_state. Non-state is when the state is not translatable to svd states
dnn_lookup <- rbind(c("non_state", -1), dnn_lookup)
dnn_lookup <- dnn_lookup %>% distinct() %>% mutate(stateID = as.numeric(stateID)) %>% as.data.frame()
dim(dnn_lookup)

# just a test
states_mapping_function("PIABpisy_1_20_22", dnn_lookup = dnn_lookup, max_height_diff_allowed = 2)

# do it for all states
svd_translated <- unlist(mclapply(split(lookup$species_code, factor(1:length(lookup$species_code))),
                                  FUN = states_mapping_function,
                                  max_height_diff_allowed = 2,
                                  dnn_lookup = dnn_lookup, mc.cores = 12, mc.cleanup = TRUE))

# bind together
lookup <- cbind(lookup, state = svd_translated)


# check how many states were not mapped
nrow(lookup[lookup$state == "non_state", ])
lookup %>% filter(state == "non_state")


# assign SVD stateId to the states 
lookup_num <- lookup %>% 
  dplyr::select(species_code, code_nr, state) %>% 
  left_join(., dnn_lookup, by = c("state")) %>%
  mutate(stateID = as.numeric(stateID))


# define final lookup table
lookup_new <- lookup_num %>% dplyr::select(species_code, stateID)
colnames(lookup_new) <- c("species_code", "code_nr")
hist(lookup_num$stateID, breaks = 100)
head(lookup_new)
lookup_new <- lookup_new %>% distinct()
setkey(lookup_new, species_code)

# check
lookup_new %>% filter(code_nr == 415)


### start processing the landscape ----------------------------------------------------

# define some objects for output
all_tiles_df <- list()

total_forest_cells <- 0
total_calc_cells <- 0

# start looping over the 3° tiles ---
for(i in 1:nrow(bbox)){
  
  # print progress  
  cat(paste0(i, " / ", nrow(bbox), " \n"))
    
  # get the canopy height raster as reference. If not available this means the tile is empty.
  name <- paste0(bbox[i, "top.coord"], bbox[i, "left.coord"])

  if(file.exists(paste0("/ETH_GlobalCanopyHeight_10m_2020_", name, "_Map_aggregated_100m.tif"))){
    height_rast <- rast(paste0("/ETH_GlobalCanopyHeight_10m_2020_", name, "_Map_aggregated_100m.tif"))
    }else{next}

  # apply forest mask! From corine landcover data to each layer
  forest_mask <- rast(paste0(path_new, "/06_gis/forest_mask_wgs.tif"))
  forest_mask_tile <- terra::crop(forest_mask, ext(height_rast))
  forest_mask_tile <- terra::resample(forest_mask_tile, height_rast)
  forest_mask_tile[forest_mask_tile == 0 | forest_mask_tile == 255] <- NA
  height_rast_masked <- terra::mask(height_rast, forest_mask_tile)

  # if empty value skip
  if(!is.finite(max(values(height_rast_masked), na.rm = T))){next}
  
  # get LAI and project to canopy height layer
  lai_sub <- rast(paste0(paste0(path_new, "/06_gis/LAI_annual_mean_2020_2021_projwgs.tif")))
  lai_sub <- terra::crop(lai_sub, height_rast_masked)
  lai_sub <- terra::project(lai_sub, height_rast_masked)
  lai_sub <- terra::mask(lai_sub, height_rast_masked)
  recl_mat <- matrix(c(0, 2, 1, 2.00000001, 4, 2, 4.00000001, 10, 3), byrow = T, nrow = 3)
  lai_class <- classify(lai_sub, recl_mat)
  
  # use focal to close gaps in the data
  lai_class <- focal(lai_class, w = c(9,9), fun = "modal", na.rm = T, na.policy = "only")
  lai_class <- focal(lai_class, w = c(17,17), fun = "modal", na.rm = T, na.policy = "only")
  lai_class <- terra::mask(lai_class, height_rast_masked)
  rm(lai_sub)

  # get vegetation STACK
  if(file.exists(paste0(path, "/03_initial_forest_states/vegetation/species_stack_", name, ".tif"))){
    veg_stack <- rast(paste0(path, "/03_initial_forest_states/vegetation/species_stack_", name, ".tif"))
  }else{next}
 
  # get the values but skip if no values are available 
  veg_sum <- app(veg_stack, sum)
  if(sum(values(veg_sum)) == 0){next}
  stack_vals <- values(veg_stack)
  snames <- colnames(stack_vals)
  if(!is.finite(max(stack_vals[[1]], na.rm = T))){next}
  
  # create 1km version of the veg layer (with count pixels)
  rast_1km <- terra::aggregate(height_rast_masked, fact=10, fun=function(x, ...) {sum(!is.na(x))})
  rast_1km[is.na(rast_1km)] <- 0
  plot(rast_1km)
  
  # get value of the stack for center point of 1km pixel and sample N values (sampleComp function)
  # but first fill corresponding cells in 100m layer, calculate stats per km.
  rast_id <- rast_1km
  rast_id <- terra::crop(rast_id, veg_stack)
  rast_id <- setValues(rast_id, 1:ncell(rast_id))
  rast_id <- terra::resample(rast_id, height_rast_masked, method = "near")
  out_df1 <- as.data.frame(rast_id, xy = T, na.rm = F)
  colnames(out_df1) <- c("x", "y", "coarse_id")
  
  # get a dataframe with the raster cells that need to be sampled
  rast_id <- terra::mask(rast_id, height_rast_masked)
  out_df2 <- as.data.frame(rast_id, xy = T, na.rm = F)
  colnames(out_df2) <- c("x", "y", "val")
  out_df2 <- out_df2 %>% 
    mutate(fine_id = 1:nrow(out_df2))
  out_df <- out_df1 %>% left_join(out_df2, by = c("x", "y"))
  
  # empty list
  out_list <- list()
  # go through the cells and sample species from the vegetation stack
  for (cl in 1:ncell(rast_1km)) {
    
    # skip if no cells
    if (rast_1km[cl] == 0)
      next;
    
    # get values from stack:
    cellvals <- stack_vals[cellFromXY(veg_stack, xyFromCell(rast_1km, cl)),]
    if (sum(is.na(cellvals)) > 0)
      next;
    # outside
    
    # sample composition states:
    df <- as.data.frame(
      sampleComp(cellvals, 
                 min_threshold = 0.1, 
                 temperature = 0.5, 
                 thresholds = nspecies_thresholds, 
                 threshold_random_range = 0.2,
                 dom_threshold = 0.40,
                 dom_range = 0.10,
                 do_debug = F)
    )
    
    
    df <- df[1:as.numeric(rast_1km[cl]),] ### we use just so many cells as we have 100m px
    states <- getStateId(as.matrix(df), snames )
    out_list[[cl]] <- states
    
  }
  
  out_vect <- unlist(as.vector(out_list))
  if(length(out_vect) != nrow(na.omit(out_df))){
    cat("check the sampler \n ")
  }
  
  out_vect <- gsub(" ", "", out_vect)
  
  # assign the veg_states to the grid cells
  veg_df <- out_df %>% 
    filter(!is.na(val)) %>% 
    arrange(coarse_id) %>% 
    mutate(veg_state = out_vect)
  
  veg_df <- out_df2 %>% left_join(veg_df, by = c("x", "y", "fine_id", "val"))
 
  
  
  ### put together the states -------------------------------------------------
 
   # get height, lai and species dataframes
  height_df <- as.data.frame(height_rast_masked, xy = T)
  lai_df <- as.data.frame(lai_class, xy = T)
  lai_df$focal_modal <- round(lai_df$focal_modal, 0)
  
  combined_df <- left_join(height_df, lai_df, by = c("x", "y")) %>% 
    left_join(., veg_df[, c("x", "y", "veg_state")], by = c("x", "y"))
  
  colnames(combined_df) <- c("x", "y", "height_state", "lai_state", "veg_state")
  
  # mutate to svd classes
  combined_df <- combined_df %>% 
    mutate(height_state = case_when(height_state < 2 ~ "0_2",
                                    height_state < 4 ~ "2_4",
                                    height_state < 6 ~ "4_6",
                                    height_state < 8 ~ "6_8",
                                    height_state < 10 ~ "8_10",
                                    height_state < 12 ~ "10_12",
                                    height_state < 14 ~ "12_14",
                                    height_state < 16 ~ "14_16",
                                    height_state < 18 ~ "16_18",
                                    height_state < 20 ~ "18_20",
                                    height_state < 22 ~ "20_22",
                                    height_state < 24 ~ "22_24",
                                    height_state < 26 ~ "24_26",
                                    height_state < 28 ~ "26_28",
                                    height_state < 30 ~ "28_30",
                                    height_state < 32 ~ "30_32",
                                    height_state < 34 ~ "32_34",
                                    height_state < 36 ~ "34_36",
                                    height_state < 38 ~ "36_38",
                                    height_state < 40 ~ "38_40",
                                    height_state < 42 ~ "40_42",
                                    height_state < 44 ~ "42_44",
                                    height_state < 46 ~ "44_46",
                                    height_state < 48 ~ "46_48",
                                    height_state < 50 ~ "48_50",
                                    height_state >= 50 ~ "50_")) 
  
  # convert to a single string
  combined_df <- combined_df %>% 
    mutate(svd_state = paste0(veg_state, "_", lai_state, "_", height_state)) %>% 
    mutate(svd_state = ifelse(grepl("NA", svd_state), NA, svd_state)) %>% 
    mutate(svd_state = ifelse(grepl("missing", svd_state), "missing", svd_state)) %>% 
    mutate(id = row_number())
  
  # we put NAs to SVD_state -9999 for the moment
  combined_df_isna <- combined_df[is.na(combined_df$svd_state), ] %>% 
    mutate(svd_num = -9999)
  
  # for the others we use the lookup table to get a svd state number
  combined_df_notna <- combined_df[!is.na(combined_df$svd_state), ]
  combined_df_notna <- combined_df_notna %>% 
    mutate(svd_num = as.numeric(paste0(lookup(svd_state, lookup_new$species_code, lookup_new$code_nr)))) %>% 
    rowwise() %>% 
    mutate(svd_num = ifelse(is.na(svd_num), -9999, svd_num)) %>% 
    mutate(svd_num = ifelse(grepl("EUGL", svd_state), -7777, svd_num)) 
  
  # combine again
  combined_df_new <- rbind(combined_df_notna, combined_df_isna)
  plot_df_height <- left_join(height_df[, c("x", "y")], combined_df_new[, c("x", "y", "height_state")], by = c("x", "y")) 
  plot_df_height <- plot_df_height %>% 
    rowwise() %>% 
    mutate(height_state = as.numeric(strsplit(height_state, "_")[[1]][1]))
  r_height <- rast(plot_df_height[, c("x", "y", "height_state")])
  plot(r_height)
  rm(combined_df_notna, combined_df_isna)
  
  
  # convert the dataframes to rasters
  plot_df <- left_join(height_df[, c("x", "y")], combined_df_new[, c("x", "y", "svd_num")], by = c("x", "y")) 
  
  # eucalyptus we put to -7777 to filter out later
  eugls <- plot_df %>% filter(svd_num == -7777) %>% nrow()
  eugls <- 100 / nrow(plot_df) * eugls
  # record also the number of non-states
  non_states <- plot_df %>% filter(svd_num == -1) %>% nrow()
  non_states <- 100 / nrow(plot_df) * non_states
  
  # convert to raster
  r_df <- rast(plot_df[, c("x", "y", "svd_num")])
  
  # NAs and non-states are put to NA
  r_df[r_df == -9999] <- NA
  r_df[r_df == -1] <- NA
  
  # use focal to close gaps from non-states
  r_df <- terra::focal(r_df, w = c(5,5), "modal", na.policy = "only", na.rm = T)
  r_df <- terra::focal(r_df, w = c(7,7), "modal", na.policy = "only", na.rm = T)

  # get the height mask
  height_rast_masked2 <- crop(height_rast_masked, r_df)
  crs(r_df) <- crs(height_rast_masked2)
  r_df <- terra::mask(r_df, height_rast_masked2)
  
  # kick out EUGL
  r_df[r_df == -7777] <- NA 
  
  # convert back to dataframe and save
  plot_df <- as.data.frame(r_df, xy = T)
  colnames(plot_df) <- c("x", "y", "svd_state_id")
  #write.table(plot_df, paste0(path, "/03_initial_forest_states/initialstate_", name, ".csv"), sep = ";")
  
  # print progress and the number of NAs, EUGLs and Non_states
  total_forest_cells <- total_forest_cells + nrow(na.omit(height_df))
  total_calc_cells <- total_calc_cells + nrow(na.omit(plot_df))
  
  print(round(100 / nrow(na.omit(height_df)) * nrow(na.omit(plot_df))), 2)
  print(paste0(round(eugls, 2), " % EUGL"))
  print(paste0(round(non_states, 2), " % non_states"))
  par(mfrow = c(1,2))
  plot(r_df)
  plot(height_rast_masked2)
 
  # save raster
  writeRaster(r_df, paste0(path, "/03_initial_forest_states/tiles_rasters/", name, ".tif"), overwrite = T)
  
  # clean up
  tmpFiles(remove = T)
  gc()
  
  rm(height_df, lai_df, veg_df, plot_df, height_rast_masked, 
     r_df, height_rast_masked2,
     combined_df, combined_df_new)
  
  print(all_time)
}




### aggregate all rasters tiles --------------------------------------------------------------------


# put together all the tiles
# needs to be done in multiple steps for computational reasons
for(p in 1:5){
  
  if(p == 1){start <- 1; end <- 44}
  if(p == 2){start <- 45; end <- 90}
  if(p == 3){start <- 91; end <- 130}
  if(p == 4){start <- 131; end <- 170}
  if(p == 5){start <- 171; end <- 234}
  
  mosaic <- 0
  
  for(i in start:end){
    
    # print progress  
    cat(paste0(i, " / ", nrow(bbox), " \n"))
    name <- paste0(bbox[i, "top.coord"], bbox[i, "left.coord"])
    
    if(file.exists(paste0(path, "/03_initial_forest_states/tiles_rasters/", name, ".tif"))){
      r <- rast(paste0(path, "/03_initial_forest_states/tiles_rasters/", name, ".tif"))
      r <- terra::project(r, crs(height_rast), res = res(height_rast), method = "mode")
      
    }else{next}
    
    mosaic <- mosaic + 1
    if(mosaic == 1){out_raster <- r}else{
      out_raster <- terra::mosaic(out_raster, r, fun = "min")

    }
    
  }
  
  plot(out_raster)
  writeRaster(out_raster, paste0(path, "/03_initial_forest_states/initial_ls_part", p, ".tif"), overwrite = T)
  
  tmpFiles(remove = T) 
  gc()
  
}


# load the three parts and combine them
part1 <- rast(paste0(path, "/03_initial_forest_states/initial_ls_part1.tif"))
part2 <- rast(paste0(path, "/03_initial_forest_states/initial_ls_part2.tif"))
part3 <- rast(paste0(path, "/03_initial_forest_states/initial_ls_part3.tif"))
part4 <- rast(paste0(path, "/03_initial_forest_states/initial_ls_part4.tif"))
part5 <- rast(paste0(path, "/03_initial_forest_states/initial_ls_part5.tif"))

# mosaic to one landscape layer
veg_merged <- terra::mosaic(part1, part2, part3, part4, part5, fun = "min")
plot(veg_merged)

# save result
writeRaster(veg_merged, paste0(path, "/03_initial_forest_states/init_veg_wgs.tif"),
            overwrite = T)

# project to laea projection
proj_laea <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
veg_merged <- terra::project(veg_merged, proj_laea, res=100, method = "mode")

# save result
writeRaster(veg_merged, paste0(path, "/03_initial_forest_states/init_veg.tif"),
            datatype="INT2S", overwrite = T)



### create rasters for SVD ---------------------------------------------------------

#load soil raster - see script 03_initial_soilconditions.R
init_soil_100 <- rast(paste0(path, "/03_initial_forest_states/initial_soil_rast.tif"))

#load init veg file
init_veg_100 <- rast(paste0(path, "/03_initial_forest_states/init_veg.tif"))

# make sure soil raster and veg layer match!
init_soil_100 <- terra::crop(init_soil_100, init_veg_100)  
init_veg_100 <- terra::crop(init_veg_100, init_soil_100)  

# the grid cells without veg are put to NAs and then masked in the soil raster
init_soil_100 <- terra::mask(init_soil_100, init_veg_100)

# save results
writeRaster(init_soil_100, 
            paste0(path, "/03_initial_forest_states/init_soil.tif"),
            overwrite = T,
            datatype="INT4S")

writeRaster(init_veg_100, 
            paste0(path, "/03_initial_forest_states/init_veg_v2.tif"),
            datatype="INT2S", overwrite=T)  

gc()
tmpFiles(remove = T)




### create residence time grid ---------------------------------------------------------------------
init_restime <- rast(paste0(path, "/03_initial_forest_states/init_veg_v2.tif"))

# fill up with random values between 1 and 10
init_restime[][!is.na(init_restime[])] <- sample(1:10, sum(!is.na(init_restime[])), replace=T)
init_restime <- terra::mask(init_restime, init_soil_100)
plot(init_restime)
writeRaster(init_restime,
            paste0(path, "/03_initial_forest_states/init_restime.tif"),
            datatype="INT2U", overwrite = T)
  

### end --------------------------------------------------------------------------------------------



