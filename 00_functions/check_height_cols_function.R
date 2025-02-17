#' Check if all columns are there
#' 
#' @description This function checks whether those three height columns are in the simulation data
#' 
#' @param x (file) simulation data
#' 
#' 
#' @author Marc Grünig
#' Last modified: 16/03/2023.
#' 
#' 
#' 

check_height_cols <- function(x){
  
  if(!("MinHeight" %in% colnames(x))){
    x <- x %>% mutate(MinHeight = rep(NA, nrow(x)))
  }
  
  if(!("MaxHeight" %in% colnames(x))){
    x <- x %>% mutate(MaxHeight = rep(NA, nrow(x)))
  }
  
  if(!("MeanHeight" %in% colnames(x))){
    x <- x %>% mutate(MeanHeight = rep(NA, nrow(x)))
  }
  
  return(x)
}