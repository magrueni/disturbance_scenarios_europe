#' get reference grid id from coordinates
#' 
#' @description For a given set of coordinates, this function gets the reference grid ID for the climate dataset
#' @param x (file) coordinates
#' 
#' 
#' @author Marc Grünig
#' Last modified: 16/03/2023.
#' 
#' 
#' 
#' 

library(terra)


path <- "/data/public/"

x_rast <- rast(paste0(path, "/Projects/Resonate/clim_data/reference_grid.tif"))
crs(x_rast) <- "+proj=longlat +datum=WGS84 +no_defs"
ref_tab <- read.csv(paste0(path, "/Projects/Resonate/clim_data/reference_grid_tab.csv"))


get_point_id <- function(x){
  
  # extract the grid ID
  lon <- as.numeric(x[1])
  lat <- as.numeric(x[2])
  
  coords <- as.data.frame(cbind(lon, lat))
  
  point_extr <- terra::extract(x_rast, coords)[2]
  point_extr <- as.numeric(point_extr)
  
  if(is.na(point_extr)){
    closest_grid <- ref_tab[which.min(abs(lon - ref_tab$wgs_x) + abs(lat - ref_tab$wgs_y)), ]
    point_extr <- as.numeric(closest_grid[1])}
  
  return(point_extr)
  
}
