#' Check if all columns are there
#' 
#' @description This function checks if there are capital letters in the string indicating that the species is dominant.
#' If this is the case, the Capital letters (the dominant species) get move to the beginning of the string.
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
library(compiler)
library(stringr)

reorder_strings <- function(x) {
  output <- mclapply(1:nrow(x), function(i) {
    if (str_detect(x$sp_state2[i], "[[:upper:]]")) {
      x$sp_state2[i]
    } else {
      substrings <- unlist(strsplit(x$sp_state2[i], "(?<=.{4})", perl = TRUE))
      sorted_substrings <- sort(substrings)
      paste0(sorted_substrings, collapse = "")
    }
  }, mc.cores = 12)
  
  x$sp_state2 <- unlist(output)
  return(x)
}


reorder_strings <- cmpfun(reorder_strings)

