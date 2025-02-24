# --------------------------------------------------------------------------------------
# Script Name: Fig2_biomes_dist_change
# Description: This script processes and visualizes disturbance rates at the biome
#              level and creates Fig2a (boxplot of disturbance rate change 
#              per biome. Fig2a shows the RCP8.5 data but all scenarios are plotted here.
#
#
# Author: Marc Grünig
# Date: 7.2.2025
#
# Input:
#   - aggregated disturbance data (25km) from script 3
#
# Output:
#   - boxplots for disturbance change per biome
#
# Notes:
# - Biomes data was obtained from Olson et al.
#
# ------------------------------------------------------------------------------


### libraries ------------------------------------------------------------------
library(rnaturalearth)
library(rnaturalearthdata)
library(MetBrewer)
library(tidyterra)
library(gridExtra)
library(raster)
library(ggplot2)
library(sf)
library(terra)
library(stringr)
library(readr)



# define path
path <- "/.../"
path_results <- paste0(path, "/04_work/svd_simulations/results_eu/")
                       

### load data ------------------------------------------------------------------

# load shp and hexagon
eu_shp <- vect(paste0(path, "/reference_grids/eu_mask.shp"))
hex_ecu <- shapefile(paste0(path, "/reference_grids/eu_mask_hexagons_25km.shp"))


# load forest mask
sf_obj_forest_mask <- read_sf(paste0(path, "/reference_grids/hex_forest_mask_25km.gpkg"))
sf_obj_forest_mask_ids <- sf_obj_forest_mask %>% st_drop_geometry()


# load simulation overview to assign sims to rcps
master_tab <- read.csv(paste0(path, "/svd_simulations/svd_simulations_ids.csv"), sep = ";")


# define rcps
rcps <- c("rcp85", "rcp45", "rcp26")


# create output and counter
r <- 0
all_maps_dat <- list()


# loop over rcps to load the data
for(c in 1:length(rcps)){
  
  # progress
  print(c)
  
  # loop over years
  for(y in seq(10, 80, 10)){
    
    r <- r + 1
    
    dat <- read_csv(paste0(path_results, "/dist_rates_25km/all_dist_rates_", y, "_", rcps[c], ".csv"))
    
    all_maps_dat[[r]] <- dat 
    
  }
  
}


# combine outputs
all_dat <- do.call(rbind, all_maps_dat)  
write_csv(all_dat, paste0(path, "/figures/plot_data/biomes_data_all.csv"))


# get the size of one hexagon
hex_size <- max(sf_obj_forest_mask$forest_area, na.rm = T)


# get europe map for plotting
proj_wgs <- "+proj=longlat +datum=WGS84 +no_defs"
proj_leae <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
eu_shp <- terra::project(eu_shp, proj_wgs)
eu_shp_sf <- st_as_sf(eu_shp)
eu_shp_sf_leae <- st_transform(eu_shp_sf, proj_leae)
world <- ne_countries(scale = "medium", returnclass = "sf")
world <- st_transform(world, proj_leae)
europe <- st_crop(world, eu_shp_sf_leae)
europe <- europe[!europe$sovereignt %in% c("Tunisia", "Russia", "Algeria", "Ukraine", "Belarus", "Moldova",
                                           "Turkey", "Cyprus", "Northern Cyprus", "Morocco", "Iceland"), ]

hex_ecu <- shapefile(paste0(path, "/reference_grids/eu_mask_hexagons_25km.shp"))
# eu_shp <- vect(eu_shp)

# load forest mask
sf_obj_forest_mask <- st_read(paste0(path, "/reference_grids/hex_forest_mask_25km.gpkg"))

# load biomes
ecoregions <- st_read(paste0(path, "/gis/ecoregions/terrestrial_ecoregions_olson.shp"))
ecoregions <- st_transform(ecoregions, proj_leae)
ecoregions <- vect(ecoregions)

# convert to raster
forest_mask_for_biomes <- rast(paste0(path, "/gis/forest_mask_new_laea_masked_nas.tif"))
forest_mask_agg <- terra::aggregate(forest_mask_for_biomes, fact = 10)
ecoreg_rast <- rasterize(ecoregions, forest_mask_agg, field = "BIOME")
reclass_matrix <- cbind(c(12, 5, 6, 11, 4, 8), # old values
                        c(1, 3, 4, 5, 2, 6))  # new class values

# Apply the reclassification
discrete_raster <- classify(ecoreg_rast, reclass_matrix)
ecoreg_rast_masked <- terra::mask(discrete_raster, europe)

# plot ecoregions
cols <- c(met.brewer("Isfahan1", 5), "lightgrey")
plot(ecoreg_rast_masked, col = cols, box = F, axes = F, legend = F)
dev.print(png, filename = paste0(path, "/figures/biomes_map.png"), width = 1500, height = 1500, res = 300)


# assign biome to each hexagon
hex_ecu <- shapefile(paste0(path, "/reference_grids/eu_mask_hexagons_25km.shp"))
hex_ecoreg <- exactextractr::exact_extract(ecoreg_rast, hex_ecu, fun = "mode")
hex_ecoreg <- cbind(sf_obj_forest_mask, hex_ecoreg)
hex_ecoreg <- hex_ecoreg %>% dplyr::select(gridid, biome = hex_ecoreg)
biomes <- c(12, 5, 8, 6, 11, 4)
names_biomes <- c("Mediterranean", "Temperate Coniferous", "Temperate Grasslands",
                  "Boreal Forests", "Tundra", "Temperate Broadleaf")
biomes_df <- data.frame(cbind(as.numeric(biomes), names_biomes))
hex_ecoreg <- cbind(hex_ecoreg, Biome = biomes_df$names_biomes[match(hex_ecoreg$biom, biomes_df$V1)]) 
hex_ecoreg <- hex_ecoreg %>% dplyr::select(gridid, biome = Biome)



### plotting ---------------------------------------------------------------------------------------

## first RCP8.5 ---

# filter the data
dat_85 <- all_dat %>% filter(rcp == "rcp85")

# prepare data and calculate dist rate per hexagon
dat_85_all <- dat_85 %>% 
  group_by(gridid, rcp, gcm, rep, year, sim) %>% 
  summarise(dist_area = sum(dist_area, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>% 
  mutate(dist_rate = as.numeric(100 / forest_area * dist_area / 10),
         dist_freq = as.numeric(10 * forest_area / dist_area)) %>% 
  mutate(dist_rate = ifelse(forest_area < land_area/20, NA, dist_rate),
         dist_freq = ifelse(forest_area < land_area/20, NA, dist_freq)) %>%
  st_as_sf(.)

all_bas <- dat_85_all %>% filter(year == 2030) %>% st_drop_geometry()
all_fut <- dat_85_all %>% filter(year == 2100) %>% st_drop_geometry()

diff_df <- left_join(all_bas, all_fut, by = c("gridid", "rcp", "gcm", "sim", "rep")) %>% 
  mutate(mean_diff = dist_rate.y - dist_rate.x) %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>% 
  st_as_sf(.)

length(na.omit(diff_df$mean_diff))

hist(diff_df$mean_diff)


diff_df_out <- diff_df %>% group_by(gridid) %>% 
  summarise(mean_diff = mean(mean_diff)) %>% 
  drop_na()

diff_df_out %>% st_drop_geometry() %>% 
  summarise(
    negative_percent = mean(mean_diff < 0) * 100,
    positive_percent = mean(mean_diff > 0) * 100,
    zero_percent = mean(mean_diff == 0) * 100
  )


# add biomes
diff_df_out_biomes <- left_join(diff_df_out, hex_ecoreg %>% st_drop_geometry(), by = c("gridid"))

diff_df_out_biomes %>% 
  st_drop_geometry() %>% 
  group_by(biome) %>% 
  summarise(
    negative_percent = mean(mean_diff < 0) * 100,
    positive_percent = mean(mean_diff > 0) * 100,
    zero_percent = mean(mean_diff == 0) * 100
  )



# for the plot - calculate how much % of the gridcells are posy 
diff_df_plot <- diff_df %>%
  left_join(., hex_ecoreg %>% st_drop_geometry(), by = c("gridid")) %>% 
  group_by(biome, gcm, rep, sim) %>% 
  summarise(
    negative_percent = mean(mean_diff < 0, na.rm = T) * 100,
    positive_percent = mean(mean_diff > 0, na.rm = T) * 100,
    zero_percent = mean(mean_diff == 0, na.rm = T) * 100
  ) %>% st_drop_geometry() %>% 
  drop_na() %>% 
  filter(biome != "Temperate Grasslands") %>% 
  rename("Biome" = "biome") 

# bar plot
level_order <- diff_df_plot %>%
  group_by(Biome) %>%
  summarise(avg = mean(as.numeric(positive_percent ))) %>%
  arrange(-avg) %>% dplyr::select(Biome) %>% unique() %>% unlist() %>% as.vector()

diff_df_plot$Biome <- factor(diff_df_plot$Biome, levels = level_order)

cols <- c(met.brewer("Isfahan1", 5)[1],
          met.brewer("Isfahan1", 5)[3],
          met.brewer("Isfahan1", 5)[2],
          met.brewer("Isfahan1", 5)[4],
          met.brewer("Isfahan1", 5)[5])


# Plotting with ggplot, flipped axis and grey dashed line at 0
boxplot1 <- ggplot(diff_df_plot,
                   aes(x = Biome, y = positive_percent , fill = Biome)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  labs(x = "", y = "Area with disturbance increase [%]") +
  scale_fill_manual(values = cols) +
  theme_classic() +
  coord_flip(ylim = c(40, 80)) +
  
  theme(legend.position = "",
        text = element_text(size = 14),
        axis.title.x = element_text(vjust = -2),
        axis.title.y = element_text(angle=0, vjust=0.5),
        plot.margin = margin(1, 1, 1, 0, "cm")) +
  guides(fill = guide_legend(reverse = TRUE)) 

boxplot1

ggsave(boxplot1,
       filename = paste0(path, "/figures/Fig2a.png"),
       width = 7, height = 8)

write_csv(diff_df_plot, paste0(path, "/figures/plot_data/figure_data/fig2a_biomes.csv"))



### same for RCP4.5 ----

# filter the data
dat_45 <- all_dat %>% filter(rcp == "rcp45")

# prepare data and calculate dist rate per hexagon
dat_45_all <- dat_45 %>% 
  group_by(gridid, rcp, gcm, rep, year, sim) %>% 
  summarise(dist_area = sum(dist_area, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>% 
  mutate(dist_rate = as.numeric(100 / forest_area * dist_area / 10),
         dist_freq = as.numeric(10 * forest_area / dist_area)) %>% 
  mutate(dist_rate = ifelse(forest_area < land_area/20, NA, dist_rate),
         dist_freq = ifelse(forest_area < land_area/20, NA, dist_freq)) %>%
  st_as_sf(.)

all_bas <- dat_45_all %>% filter(year == 2030) %>% st_drop_geometry()
all_fut <- dat_45_all %>% filter(year == 2100) %>% st_drop_geometry()

diff_df <- left_join(all_bas, all_fut, by = c("gridid", "rcp", "gcm", "sim", "rep")) %>% 
  mutate(mean_diff = dist_rate.y - dist_rate.x) %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>% 
  st_as_sf(.)

length(na.omit(diff_df$mean_diff))

hist(diff_df$mean_diff)


diff_df_out <- diff_df %>% group_by(gridid) %>% 
  summarise(mean_diff = mean(mean_diff)) %>% 
  drop_na()

diff_df_out %>% st_drop_geometry() %>% 
  summarise(
    negative_percent = mean(mean_diff < 0) * 100,
    positive_percent = mean(mean_diff > 0) * 100,
    zero_percent = mean(mean_diff == 0) * 100
  )


# add biomes
diff_df_out_biomes <- left_join(diff_df_out, hex_ecoreg %>% st_drop_geometry(), by = c("gridid"))

diff_df_out_biomes %>% 
  st_drop_geometry() %>% 
  group_by(biome) %>% 
  summarise(
    negative_percent = mean(mean_diff < 0) * 100,
    positive_percent = mean(mean_diff > 0) * 100,
    zero_percent = mean(mean_diff == 0) * 100
  )



# for the plot - calculate how much % of the gridcells are posy 
diff_df_plot <- diff_df %>%
  left_join(., hex_ecoreg %>% st_drop_geometry(), by = c("gridid")) %>% 
  group_by(biome, gcm, rep, sim) %>% 
  summarise(
    negative_percent = mean(mean_diff < 0, na.rm = T) * 100,
    positive_percent = mean(mean_diff > 0, na.rm = T) * 100,
    zero_percent = mean(mean_diff == 0, na.rm = T) * 100
  ) %>% st_drop_geometry() %>% 
  drop_na() %>% 
  filter(biome != "Temperate Grasslands") %>% 
  rename("Biome" = "biome") %>% 
  mutate(Biome = ifelse(grepl("Boreal Forests", Biome), "BF",
                        ifelse(grepl("Mediterranean", Biome), "M",
                               ifelse(grepl("Temperate Broadleaf", Biome), "TB",
                                      ifelse(grepl("Tundra", Biome), "T",
                                             ifelse(grepl("Temperate Coniferous", Biome), "TC", "TG")))))) 


# bar plot
level_order <- diff_df_plot %>%
  group_by(Biome) %>%
  summarise(avg = mean(as.numeric(positive_percent ))) %>%
  arrange(-avg) %>% dplyr::select(Biome) %>% unique() %>% unlist() %>% as.vector()

diff_df_plot$Biome <- factor(diff_df_plot$Biome, levels = level_order)

cols <- c(met.brewer("Isfahan1", 5)[1],
          met.brewer("Isfahan1", 5)[3],
          met.brewer("Isfahan1", 5)[2],
          met.brewer("Isfahan1", 5)[4],
          met.brewer("Isfahan1", 5)[5])



# Plotting with ggplot, flipped axis and grey dashed line at 0
boxplot2 <- ggplot(diff_df_plot,
                   aes(x = Biome, y = positive_percent , fill = Biome)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +  # Add dashed line at y = 0
  labs(x = "Biome", y = "Area with disturbance increase [%]") +
  scale_fill_manual(values = cols) +
  theme_classic() +
  coord_flip(ylim = c(40, 80)) +
  
  theme(legend.position = "",
        text = element_text(size = 16),
        axis.title.x = element_text(vjust = -2),
        axis.title.y = element_text(vjust = 2),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  guides(fill = guide_legend(reverse = TRUE))  # Reverse legend order
boxplot2



### same for RCP2.6

# filter the data
dat_26 <- all_dat %>% filter(rcp == "rcp26")

# prepare data and calculate dist rate per hexagon
dat_26_all <- dat_26 %>% 
  group_by(gridid, rcp, gcm, rep, year, sim) %>% 
  summarise(dist_area = sum(dist_area, na.rm = T)) %>% 
  ungroup() %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>% 
  mutate(dist_rate = as.numeric(100 / forest_area * dist_area / 10),
         dist_freq = as.numeric(10 * forest_area / dist_area)) %>% 
  mutate(dist_rate = ifelse(forest_area < land_area/20, NA, dist_rate),
         dist_freq = ifelse(forest_area < land_area/20, NA, dist_freq)) %>%
  st_as_sf(.)

all_bas <- dat_26_all %>% filter(year == 2030) %>% st_drop_geometry()

all_fut <- dat_26_all %>% filter(year == 2100) %>% st_drop_geometry()

diff_df <- left_join(all_bas, all_fut, by = c("gridid", "rcp", "gcm", "sim", "rep")) %>% 
  mutate(mean_diff = dist_rate.y - dist_rate.x) %>% 
  left_join(., sf_obj_forest_mask, by = "gridid") %>% 
  st_as_sf(.)

length(na.omit(diff_df$mean_diff))

hist(diff_df$mean_diff)


diff_df_out <- diff_df %>% group_by(gridid) %>% 
  summarise(mean_diff = mean(mean_diff)) %>% 
  drop_na()

diff_df_out %>% st_drop_geometry() %>% 
  summarise(
    negative_percent = mean(mean_diff < 0) * 100,
    positive_percent = mean(mean_diff > 0) * 100,
    zero_percent = mean(mean_diff == 0) * 100
  )


# add biomes
diff_df_out_biomes <- left_join(diff_df_out, hex_ecoreg %>% st_drop_geometry(), by = c("gridid"))

diff_df_out_biomes %>% 
  st_drop_geometry() %>% 
  group_by(biome) %>% 
  summarise(
    negative_percent = mean(mean_diff < 0) * 100,
    positive_percent = mean(mean_diff > 0) * 100,
    zero_percent = mean(mean_diff == 0) * 100
  )



# for the plot - calculate how much % of the gridcells are posy 
diff_df_plot <- diff_df %>%
  left_join(., hex_ecoreg %>% st_drop_geometry(), by = c("gridid")) %>% 
  group_by(biome, gcm, rep, sim) %>% 
  summarise(
    negative_percent = mean(mean_diff < 0, na.rm = T) * 100,
    positive_percent = mean(mean_diff > 0, na.rm = T) * 100,
    zero_percent = mean(mean_diff == 0, na.rm = T) * 100
  ) %>% st_drop_geometry() %>% 
  drop_na() %>% 
  filter(biome != "Temperate Grasslands") %>% 
  rename("Biome" = "biome") %>% 
  mutate(Biome = ifelse(grepl("Boreal Forests", Biome), "BF",
                        ifelse(grepl("Mediterranean", Biome), "M",
                               ifelse(grepl("Temperate Broadleaf", Biome), "TB",
                                      ifelse(grepl("Tundra", Biome), "T",
                                             ifelse(grepl("Temperate Coniferous", Biome), "TC", "TG")))))) 


# bar plot
level_order <- diff_df_plot %>%
  group_by(Biome) %>%
  summarise(avg = mean(as.numeric(positive_percent ))) %>%
  arrange(-avg) %>% dplyr::select(Biome) %>% unique() %>% unlist() %>% as.vector()

diff_df_plot$Biome <- factor(diff_df_plot$Biome, levels = level_order)

cols <- c(met.brewer("Isfahan1", 5)[1],
          met.brewer("Isfahan1", 5)[3],
          met.brewer("Isfahan1", 5)[2],
          met.brewer("Isfahan1", 5)[4],
          met.brewer("Isfahan1", 5)[5])



# Plotting with ggplot, flipped axis and grey dashed line at 0
boxplot3 <- ggplot(diff_df_plot,
                   aes(x = Biome, y = positive_percent , fill = Biome)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +  # Add dashed line at y = 0
  labs(x = "Biome", y = "Area with disturbance increase [%]") +
  scale_fill_manual(values = cols) +
  theme_classic() +
  coord_flip(ylim = c(30, 80)) +
  
  theme(legend.position = "",
        text = element_text(size = 16),
        axis.title.x = element_text(vjust = -2),
        axis.title.y = element_text(vjust = 2),
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  guides(fill = guide_legend(reverse = TRUE))  # Reverse legend order
boxplot3
        

### end ------------------------------------------------------------------------
