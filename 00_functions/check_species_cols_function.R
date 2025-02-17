#' Check if all columns are there
#' 
#' @description This function checks if the columns for the 5 species proportionas are in the simulation data
#' 
#' @param x (file) simulation data
#' 
#' 
#' @author Marc Grünig
#' Last modified: 16/03/2023.
#' 
#' 
#' 
#' 
check_species_cols <- function(x){
  
  if(!("Species2" %in% colnames(x))){
    x <- x %>% mutate(Species2 = rep(NA, nrow(x)))
  }
  
  if(!("Proportion2" %in% colnames(x))){
    x <- x %>% mutate(Proportion2 = rep(NA, nrow(x)))
  }
  
  if(!("Species3" %in% colnames(x))){
    x <- x %>% mutate(Species3 = rep(NA, nrow(x)))
  }
  
  if(!("Proportion3" %in% colnames(x))){
    x <- x %>% mutate(Proportion3 = rep(NA, nrow(x)))
  }
  
  if(!("Species4" %in% colnames(x))){
    x <- x %>% mutate(Species4 = rep(NA, nrow(x)))
  }
  
  if(!("Proportion4" %in% colnames(x))){
    x <- x %>% mutate(Proportion4 = rep(NA, nrow(x)))
  }
  
  if(!("Species5" %in% colnames(x))){
    x <- x %>% mutate(Species5 = rep(NA, nrow(x)))
  }
  
  if(!("Proportion5" %in% colnames(x))){
    x <- x %>% mutate(Proportion5 = rep(NA, nrow(x)))
  }
  
  return(x)
}
