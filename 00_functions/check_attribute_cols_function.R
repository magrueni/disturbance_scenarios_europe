#' Check if all columns are there
#' 
#' @description This function checks if those three attribute columns are in the simulation data
#' 
#' @param x (file) simulation data
#' 
#' 
#' @author Marc Grünig
#' Last modified: 16/03/2023.
#' 
#' 
#' 

check_attribute_cols <- function(x){
  
  if(!("Volume" %in% colnames(x))){
    x <- x %>% mutate(Volume = rep(NA, nrow(x)))
  }
  
  if(!("MeanDiameter" %in% colnames(x))){
    x <- x %>% mutate(MeanDiameter = rep(NA, nrow(x)))
  }
  
  if(!("StemNumber" %in% colnames(x))){
    x <- x %>% mutate(StemNumber = rep(NA, nrow(x)))
  }
  
  return(x)
}

