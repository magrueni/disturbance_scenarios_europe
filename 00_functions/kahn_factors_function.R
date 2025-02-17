#' Check if all columns are there
#' 
#' @description This function calculated the dominant height from Maximum Height for those species we have the Kahn-factors.
#' @param x (file) simulation data
#' 
#' 
#' @author Marc Grünig
#' Last modified: 16/03/2023.
#' 
#' 
#' 
#' 

library(compiler)

kahn_factors_function <- function(x){
  
  x$dom_species <- strsplit(as.character(x$sp_state2), "(?<=.{4})", perl = TRUE)[1]
  
  x <- x %>% mutate(Hdom = ifelse(dom_species %in% c("PIAB", "ABAL", "LADE"), MaxHeight / 1.138,
                                  ifelse(dom_species %in% c("PISY"), MaxHeight / 1.18, 
                                         ifelse(dom_species %in% c("FASY"), MaxHeight / 1.132, 
                                                ifelse(dom_species %in% c("QURO", "QUCE", "QUCO", "QUFR", "QUIL", "QUPE", "QUPU", "QUSU", "QUTR"),
                                                       MaxHeight / 1.184, MaxHeight / 1.14)))))
  

  
  
  # if not, we calc the dominant height as before 
  x <- x %>% 
    mutate(Hdom = Hdom*Proportion1/100 + (MaxHeight*0.8)*(1-Proportion1/100))
  
  
  x <- x[, -c("dom_species")]
  
  return(x)
  
}

kahn_factors_function_cmp <- cmpfun(kahn_factors_function)