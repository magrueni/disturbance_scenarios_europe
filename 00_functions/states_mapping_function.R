library(tidyverse)


### find height function ----------------------

#' find_height
#' 
#' @description Identifies the state hight from a lookup table whose heights are within a specified range of the target height.
#' The function uses an iterative process to expand the allowable height difference until matches are found or the maximum allowed difference is reached.
#' 
#' @param species_match (vector) Indices referring to species in `lookup_states` that match the target species composition.
#' @param lookup_states (vector) Vector of strings representing state information, formatted as `species_LAI_height`.
#' @param max_height_diff_allowed (numeric) Maximum allowable height difference.
#' @param state_height (numeric) Height of the target state.
#' 
#' @return A vector of indices in `lookup_states` that meet the height difference criteria.
#' 
#' @author Your Name
#' Last modified: 10/01/2025.
#' 


find_height <- function(species_match, lookup_states, max_height_diff_allowed, state_height){
  
  max_height_diff_allowed <- max_height_diff_allowed
  height_diff_allowed <- -2
  height_matches <- NULL
  
  while(length(height_matches) == 0 & height_diff_allowed <= max_height_diff_allowed){
    
    height_diff_allowed <- height_diff_allowed + 2
    
    # Check if the height is within 2 meters of the state
    height_matches <- lapply(1:length(species_match), function(x) {
      match_height <- as.numeric(strsplit(lookup_states[species_match[x]], "_")[[1]][3])
      if (abs(match_height - state_height) <= height_diff_allowed) {
        species_match[x]
      } else {
        NA
      }
    })
    
    height_matches <- unlist(as.vector(height_matches[!is.na(height_matches)]))
    
  }
  
  return(as.vector(na.omit(height_matches)))
}# close function



### find LAI function ------------
#' find_lai
#' 
#' @description Identifies the LAI class from a set of height matches whose Leaf Area Index (LAI) values match the target LAI within a specified tolerance.
#' If no match is found, a random height match is selected.
#' 
#' @param height_matches (vector) Indices from `lookup_states` that match the height criteria.
#' @param lookup_states (vector) Vector of strings representing state information, formatted as `species_LAI_height`.
#' @param state_lai (numeric) LAI of the target state.
#' 
#' @return The best-matching state as a string from `lookup_states`.
#' 
#' @author Your Name
#' Last modified: 10/01/2025.
#' 



find_lai <- function(height_matches, lookup_states, state_lai){
  
  lai_diff_allowed <- -1
  lai_matches <- NULL
  
  while(length(lai_matches) == 0 & lai_diff_allowed <= 2){
    
    lai_diff_allowed <- lai_diff_allowed + 1
    
    lai_matches <- lapply(1:length(height_matches), function(x) {
      match_lai <- as.numeric(strsplit(lookup_states[height_matches[x]], "_")[[1]][2])
      if(is.na(match_lai)){match_lai <- 10} # if no LAI - write ten so we are for sure above the allowed difference
      if (abs(match_lai - state_lai) == lai_diff_allowed) {
        height_matches[x]
      } else {
        NA
      }
    })
  }
  
  lai_matches <- unlist(as.vector(lai_matches[!is.na(lai_matches)]))
  
  # if we find a state with lai matching we take that one
  if(length(lai_matches) > 0){
    final_match <- dplyr::sample_n(as.data.frame(lai_matches), 1) %>% unlist()
    new_state <- lookup_states[final_match]
    return(new_state)
    
    # otherwise we take another of the matching heights  
  }else{
    final_match <- dplyr::sample_n(as.data.frame(height_matches), 1) %>% unlist()
    new_state <- lookup_states[final_match]
    return(new_state)
    
  }
}# close function



### find state function ----------------------


#' find_state
#' 
#' @description Determines the most suitable state based on species composition, height, and LAI values.
#' Handles cases for different numbers of species (up to four) and ensures the best match is returned using helper functions.
#' 
#' @param sp (character) Species composition of the target state.
#' @param lookup_states (vector) Vector of strings representing state information, formatted as `species_LAI_height`.
#' @param max_height_diff_allowed (numeric) Maximum allowable height difference.
#' @param state_lai (numeric) LAI of the target state.
#' @param state_height (numeric) Height of the target state.
#' 
#' @return A string representing the most suitable state or `"non_state"` if no match is found.
#' 
#' Last modified: 10/01/2025.
#' 


find_state <- function(sp, lookup_states, max_height_diff_allowed, state_lai, state_height){
  
  max_height_diff_allowed <- max_height_diff_allowed
  sp_dom <- ifelse(nchar(sp) > 4, substr(sp, 0, 4), sp)
  sp2 <- ifelse(nchar(sp) > 4, substr(sp, 5, 8), "")
  sp3 <- ifelse(nchar(sp) > 8, substr(sp, 9, 12), "")
  sp4 <- ifelse(nchar(sp) > 12, substr(sp, 13, 16), "")
  
  state_height <- state_height
  state_lai <- state_lai
  
  # find species composition matching
  # species_match <- lookup_states[grep(sp, lookup_states)]
  species_match <- grep(paste0("^", sp, "_"), lookup_states)
  if(length(species_match) > 0){
    heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
    if(length(heights) == 0){species_match <- NULL}else{
      new_state <- find_lai(heights, lookup_states, state_lai)
    }
  }
  
  
  # If there is no state but we have a dominant species, we check if there are some other states with this species as dominant and no other speices
  if(length(species_match) == 0 & str_detect(sp_dom, "[[:upper:]]")){
    species_match <- grep(paste0(toupper(sp_dom), "_"), lookup_states)
    if(length(species_match) > 0){
      heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
      if(length(heights) == 0){species_match <- NULL}else{
        new_state <- find_lai(heights, lookup_states, state_lai)
      }
    }
  }
  

  
  # if still nothing and there is 
  # if its only 1 species
  if(nchar(sp) == 4){
    
    if(length(species_match) == 0){
      species_match <- grep(paste0("^", tolower(sp_dom), "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    if(length(species_match) == 0){
      species_match <- grep(paste0(toupper(sp_dom), "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
  }
  
  # if two species
  if(nchar(sp) == 8){
    
    if(length(species_match) == 0){
      species_match <- grep(paste0("^", tolower(sp_dom), sp2, "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    if(length(species_match) == 0){
      species_match <- grep(paste0(toupper(sp_dom), "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    if(length(species_match) == 0){
      species_match <- grep(paste0("^", sp2, "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    if(length(species_match) == 0){
      species_match <- grep(paste0(toupper(sp2), "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    if(length(species_match) == 0){
      species_match <- grep(paste0("^", tolower(sp_dom), "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    
    
  }
  
  # if three species
  if(nchar(sp) == 12){
    
    if(length(species_match) == 0){
      species_match <- grep(paste0("^", tolower(sp_dom), sp2, sp3, "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    if(length(species_match) == 0){
      species_match <- grep(paste0("^", tolower(sp_dom), sp3, "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    if(length(species_match) == 0 & !(str_detect(sp_dom, "[[:upper:]]"))){
      species_match <- grep(paste0("^", sp2, sp3, "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    if(length(species_match) == 0){
      species_match <- grep(paste0(toupper(sp_dom), "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    if(length(species_match) == 0){
      species_match <- grep(paste0("^", tolower(sp_dom), "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
  }
  
  # if four species
  if(nchar(sp) == 16){
    
    if(length(species_match) == 0){
      species_match <- grep(paste0("^", sp_dom, sp2, sp3, sp4), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    if(length(species_match) == 0){
      species_match <- grep(paste0("^", sp_dom, sp2, sp3, "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    if(length(species_match) == 0){
      species_match <- grep(paste0("^", sp_dom, sp2, sp4, "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    if(length(species_match) == 0){
      species_match <- grep(paste0("^", sp_dom, sp3, sp4, "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    if(length(species_match) == 0 & !(str_detect(sp_dom, "[[:upper:]]"))){
      species_match <- grep(paste0("^", sp2, sp4, "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    if(length(species_match) == 0 & !(str_detect(sp_dom, "[[:upper:]]"))){
      species_match <- grep(paste0("^", sp3, sp4, "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    if(length(species_match) == 0){
      species_match <- grep(paste0(toupper(sp_dom), "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
    if(length(species_match) == 0){
      species_match <- grep(paste0("^", tolower(sp_dom), "_"), lookup_states)
      if(length(species_match) > 0){
        heights <- find_height(species_match, lookup_states, max_height_diff_allowed, state_height)
        if(length(heights) == 0){species_match <- NULL}else{
          new_state <- find_lai(heights, lookup_states, state_lai)
        }
      }
    }
    
  }
  
  if(length(species_match) == 0){
    return("non_state")
  }else{
    return(new_state)
  }
  
} # close function




### mapping function ------------------------------------
#' states_mapping_function
#' 
#' @description Maps a given state to the most appropriate state from a lookup table, using criteria for species composition, height, and LAI.
#' Utilizes helper functions `find_state`, `find_height`, and `find_lai` for refining the match.
#' 
#' @param state (character) The state to be mapped, formatted as `species_LAI_height`.
#' @param dnn_lookup (data.frame) Data frame where the first column contains state strings formatted as `species_LAI_height`.
#' @param max_height_diff_allowed (numeric) Maximum allowable height difference for matching.
#' 
#' @return The most appropriate state as a string or `NA` if no valid match is found.
#' 
#' Last modified: 10/01/2025.
#' 




states_mapping_function <- function(state, dnn_lookup, max_height_diff_allowed){
  
  max_height_diff_allowed <- max_height_diff_allowed
  lookup_states <- dnn_lookup[,1]
  lookup_states_dom <- substr(lookup_states, 0, 4)
  
  if(grepl("NA", state)){return(NA)}
  
  
  if(state %in% lookup_states){
    
    return(state)
    
  }else{
    
    # extract the details
    sp <- strsplit(state, "_")[[1]][1]
    state_lai <- as.numeric(strsplit(state, "_")[[1]][2])
    state_height <-  as.numeric(strsplit(state, "_")[[1]][3])
    
    
    new_state <- find_state(sp, lookup_states, max_height_diff_allowed, state_lai, state_height)
    return(new_state)
    
  }
  
}# close function
