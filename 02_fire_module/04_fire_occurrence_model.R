# --------------------------------------------------------------
# Fire Module Scipt 04: Fire Occurrence Modeling and Prediction
# --------------------------------------------------------------
#
# This script performs fire occurrence modeling and prediction using environmental predictors in a spatial GAM.
# It includes steps for loading and preparing data, calibrating a predictive model, and projecting fire probabilities under 
# historical and future climate scenarios.
# 
# Calibrating this model, is the third core element for predicting fire events (see script 05_create_fire_events)
#
# Notes: the pre-calculated results of this script are located in ./fire_module/models/
#-------------------------------------------------------------------------------



### libraries --------------------------
library(terra)
library(sf)
library(tidyverse)
library(stars)
library(PresenceAbsence)
library(grid)
library(gridExtra)
library(MuMIn)
library(lme4)
library(DHARMa)
library(sp)
library(data.table)
library(mgcv)
library(corrplot)



# define folder path
path <- "/.../"


### load fire complexes -----------------------------------------------
cntrs <- list.files(paste0(path, "/fire_module/complexes"), ".gpkg$") %>% 
  str_sub(., 14, nchar(.) - 5)
cntrs <- cntrs[grepl("*eps150", cntrs)] %>% str_sub(., 0, nchar(.) - 7)

k <- 0
out <- vector("list", length(cntrs))

cntrs.fullnam <- list.files(paste0(path, "/fire_module/complexes"), ".gpkg", full.names = TRUE)
cntrs.fullnam <- cntrs.fullnam[grepl("*eps150", cntrs.fullnam)]

for (i in cntrs.fullnam) {
  
  k <- k + 1
  out[[k]] <- read_sf(i) 
  
}

a <- out %>%
  bind_rows(.id = "clust")


# read the grid
grid <- read_sf(paste0(path, "/fire_module/climate/climategrid_epsg3035.gpkg"))
grid_raster <- st_rasterize(grid)
grid_r <- rast(grid_raster)

d_sf <- a %>% 
  st_as_sf(., coords = c("x", "y"), crs = st_crs(grid))



### load all the predictor variables --------------------------------------------

#  define CRS
crs_leae <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs" 

# refernece grid
ref_grid_100km <- rast(paste0(path, "/07_reference_grids/reference_grid_100km.tif"))
ref_grid_10km <- rast(paste0(path, "/07_reference_grids/reference_grid_10km.tif"))
ref_grid_1km <- rast(paste0(path, "/07_reference_grids/reference_grid_1km.tif"))

# mask for Europe
eu_mask_wgs <- rast(paste0(path, "/fire_module/gis/EU_Mask.tif"))
eu_mask <- terra::project(eu_mask_wgs, crs_leae)
eu_mask <- resample(eu_mask, dem)
rm(eu_mask_wgs)


# digital elevation model from Copernicus
dem <- rast(paste0(path, "/fire_module/gis/europe_dem_10km.tif"))

# calculate slope from the DEM
slope <- terra::terrain(dem, "slope", unit = "degrees")
slope <- terra::resample(slope, ref_grid_10km, method = "bilinear")

# calculate aspect
aspect <- terrain(dem, "aspect", unit = "degrees")
# Re-classify the aspect into four main directions: N, E, S, and W
recl_mat <- matrix(c(0, 67, 1, 67, 113, 2, 113, 248, 3, 248, 293, 4, 293, 360, 1), ncol = 3, byrow = T)
aspect_reclassified <- classify(aspect, recl_mat)

# Aggregate the aspect raster from 100m to 1km resolution and calculate the proportion of each aspect
proportion_aspect <- aggregate(aspect_reclassified, fact = 10, fun = function(x) {
  x <- na.omit(x)
  proportion <- length(x[x == 3]) / length(x)
  return(proportion)
})


# climate data from Euro-Cordex
temp_base <- rast(paste0(path, "fire_module/gis/mean_tas_historical.tif"))
seas_base <- rast(paste0(path, "fire_module/gis/mean_seas_historical.tif"))
prec_base <- rast(paste0(path, "fire_module/gis/mean_prec_historical.tif"))


# LAI from Modis
lai <- rast(paste0(path, "/fire_module/gis/LAI_layer.tif"))


# pop density from HYDE dataset
pop_dens <- rast(paste0(path, "/fire_module/gis/pop_dense_layer.tif"))


# distance to infrastructure derived from Copernicus landcover data
infrastructure <- rast(paste0(path, "/fire_module/gis/landcover_distance.tif"))


# distance to water derived from Copernicus landcover data
water <- rast(paste0(path, "/fire_module/gis/water_distance.tif"))


# summer lightning density from https://doi.pangaea.de/10.1594/PANGAEA.904253
lightning <- rast(paste0(path, "/fire_module/gis/lightning_density.tif"))




### create a 10km layer to get a presence absence layer --------------------------------------------
fire_vect <- vect(d_sf)
absence_rast <- ref_grid_1km
absence_rast[!is.na(absence_rast)] <- 0
absence_rast <- aggregate(absence_rast, fact = 10, fun = "mean", na.rm = T)
plot(absence_rast)

absence_rast_new <- terra::mask(absence_rast, fire_vect, touches = TRUE, inverse = T)
absence_rast_new[!is.na(absence_rast_new)] <- 0
absence_rast_new[is.na(absence_rast_new)] <- 1

fire_rast <- terra::mask(absence_rast_new, absence_rast)

plot(fire_rast)
plot(fire_vect, add = T, col = "red")



### extract data from all rasters ------------------------------------------------


cellsXY <- as.data.frame(fire_rast, xy = T)

# extract env data for all gridcells
env <- as.data.frame(cbind(cellsXY,
                           gridid_10km = terra::extract(ref_grid_10km, cellsXY[, c("x", "y")])[,2], 
                           pop_dens = terra::extract(pop_dens, cellsXY[, c("x", "y")])[,2],
                           dist_road = terra::extract(infrastructure, cellsXY[, c("x", "y")])[,2],
                           water = terra::extract(water, cellsXY[, c("x", "y")])[,2],
                           lightning = terra::extract(lightning, cellsXY[, c("x", "y")])[,2],
                           
                           dem = terra::extract(dem, cellsXY[, c("x", "y")])[,2], 
                           slope = terra::extract(slope, cellsXY[, c("x", "y")])[,2], 
                           aspect = terra::extract(proportion_aspect, cellsXY[, c("x", "y")])[,2],
                           
                           temp_base = terra::extract(temp_base, cellsXY[, c("x", "y")])[,2],
                           prec_base = terra::extract(prec_base, cellsXY[, c("x", "y")])[,2],
                           seas_base = terra::extract(seas_base, cellsXY[, c("x", "y")])[,2]))


# write.csv(env, paste0(path, "/fire_module/models/presabs_extract_complexes.csv"))
gc()
tmpFiles(remove=TRUE, current=T, orphan=T)



### modelling --------------------------------------------------------------------------------------

dat <- data.table::fread(paste0(path, "/fire_module/models/presabs_extract_complexes.csv"))
dat <- dat %>% 
  dplyr::rename("P" = "gridid") %>% 
  drop_na()


# set seed
set.seed(42)

# sample the same amount of absences as we have presences
presences <- dat[dat$P == 1, ]
absences <- dat[dat$P == 0, ] %>% sample_n(nrow(presences)) 

dat <- rbind(presences, absences)

# scale the variables
df_preds_all_scaled <- dat  %>% 
  dplyr::select(P, x, y, pop_dens, dist_road, water, lightning, dem, slope, 
                aspect, temp_base, prec_base, seas_base) %>% 
  mutate(prec_base = log(prec_base))

df_preds_all_scaled <- na.omit(df_preds_all_scaled)



### check correlations in predictor 
# 
# corr_matrix <- cor(df_preds_all_scaled[, 5:ncol(df_preds_all_scaled)],
#                    use = "complete.obs")
# corrplot(corr_matrix, method = "color")


### calibrate fire occurence model --------------------------------------------------

mod_full <- mgcv::gam(P ~ s(pop_dens) + s(dist_road) + s(lightning) + s(water) +
                        s(dem) + s(slope) + s(aspect) +
                        s(temp_base) + s(prec_base) + s(seas_base),
                      data = df_preds_all_scaled, family = "binomial")



plot.gam(mod_full, pages=1)



### testing ---------------------------------------------------------------------
# library(DHARMa)
# testDispersion(mod_full)
# simulationOutput <- simulateResiduals(fittedModel = mod_full, plot = T)
# residuals(simulationOutput)
# residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
# testZeroInflation(simulationOutput)
# plot(simulationOutput)
# 
# plotResiduals(simulationOutput, form = df_preds_all_scaled$temp_base)
# testQuantiles(simulationOutput)
# simulationOutput$scaledResiduals
# plot(simulationOutput, quantreg = T)



### evaluation -----------------------------------------------------------------


# predict
Prob <- as.vector(predict(mod_full, type = "response"))

# get kappa threshold
kappa.th <- optimal.thresholds(cbind(c(1:length(Prob)), df_preds_all_scaled[,"P"], Prob), opt.methods = "MaxKappa")

# binary for confusion matrix
Pred <- if_else(Prob > as.numeric(kappa.th[2]), 1, 0)
ConfusionMatrix <- table(Pred, pull(df_preds_all_scaled, P)) #`pull` results in a vector

# AUC
Pred <- ROCR::prediction(Prob, as.vector(pull(df_preds_all_scaled, P)))
AUC <- ROCR::performance(Pred, measure = "auc")
AUC <- AUC@y.values[[1]]
AUC

# r2 
summary_gam <- summary(mod_full)
r2 <- summary_gam$r.sq

eval <- cbind(r2, AUC, kappa.th[2])
print(eval)


summary.gam(mod_full)



### predict model --------------------------------------------------------------

# create newdata
newdata <- dat %>% 
  mutate(prec_base = log(prec_base)) 

# make prediction
prob <- predict(mod_full, newdata = newdata, type="response")
prob_df <- cbind(newdata[, c("x", "y")], prob)

# create raster from prediction
prob_rast <- rast(prob_df)

# divide through number of years because prediction is for 1 year
prob_rast <- prob_rast / 35
prob_rast[prob_rast < 0] <- 0
plot(prob_rast)
crs(prob_rast) <- crs_leae


# save
# writeRaster(prob_rast, paste0(path, "/fire_module/models/probability_raster_baseline_cordex_biascorrected_10k.tif"), overwrite = T)



### to future climate conditions ------------------------------------------------

years <- c("2011-2040", "2041-2070", "2071-2100")

rcps <- c("rcp85", "rcp45", "rcp26")
for(c in rcps){
  
  # read future climate
  for(y in years){
    
    temp_fut <- rast(paste0(path, "/fire_module/gis/mean_tas_", y, "_", c, "_biascorrected.tif"))
    seas_fut <- rast(paste0(path, "/fire_module/gis/mean_seas_", y, "_", c, "_biascorrected.tif"))
    prec_fut <- rast(paste0(path, "/fire_module/gis/mean_prec_", y, "_", c, "_biascorrected.tif"))

    temp_base <- terra::resample(temp_base, ref_grid_10km, method = "bilinear")
    seas_base <- terra::resample(seas_base, ref_grid_10km, method = "bilinear")
    prec_base <- terra::resample(prec_base, ref_grid_10km, method = "bilinear")
    
    env_fut <- as.data.frame(cbind(cellsXY,
                                   pop_dens = terra::extract(pop_dens, cellsXY[, c("x", "y")])[,2],
                                   dist_road = terra::extract(infrastructure, cellsXY[, c("x", "y")])[,2],
                                   water = terra::extract(water, cellsXY[, c("x", "y")])[,2],
                                   lightning = terra::extract(lightning, cellsXY[, c("x", "y")])[,2],
                                   
                                   dem = terra::extract(dem, cellsXY[, c("x", "y")])[,2], 
                                   slope = terra::extract(slope, cellsXY[, c("x", "y")])[,2], 
                                   aspect = terra::extract(proportion_aspect, cellsXY[, c("x", "y")])[,2],
                                   
                                   temp_base = terra::extract(temp_fut, cellsXY[, c("x", "y")])[,2],
                                   prec_base = terra::extract(prec_fut, cellsXY[, c("x", "y")])[,2],
                                   seas_base = terra::extract(seas_fut, cellsXY[, c("x", "y")])[,2]))
    
    env_fut <- env_fut %>% 
      mutate(prec_base = log(prec_base)) 
    
    # make prediction
    prob <- predict(mod_full, newdata = env_fut, type="response")
    prob_df <- cbind(env[, c("x", "y")], prob)
    
    # create raster from prediction
    prob_rast_fut <- rast(prob_df)
    
    # divide through number of years because prediction is for 1 year
    prob_rast_fut <- prob_rast_fut / 35
    prob_rast_fut[prob_rast_fut] <- 0
    plot(prob_rast_fut)
    crs(prob_rast_fut) <- crs_leae
    
    # save
    # writeRaster(prob_rast_fut, paste0(path, "/fire_module/models/probability_rast_fut_", y, "_", c, "_biascorrected_10k.tif"), overwrite = T)
    
  }
  
}


