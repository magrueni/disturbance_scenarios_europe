#' Check if all columns are there
#' 
#' @description This function checks if there are wierd values in the simulations
#' # for some simulations I found that there is one single entry of another species with 100% proportion
#' # this needs to be sorted out....
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
#' 
#' 
#' 



check_species_condition <- function(df) {
  result <- df %>%
    group_by(Species1) %>%
    summarise(num_occurrences = n(), total_proportion = mean(Proportion1)) %>%
    filter(num_occurrences < 10, total_proportion == 100)
  
  return(nrow(result) > 0)
}
