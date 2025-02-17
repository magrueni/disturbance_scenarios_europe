#' Check if all columns are there
#'
#' @description This function creates examples from the simulation dataset.
#' For each year that passes in the simulation, the residence time increases by one year if the vegetation state remains the same as the previous year; otherwise, it resets to 1.
#' We then identify whether the vegetation state will change within the next ten years. If a state change occurs, we estimate the timing (TARGET TIME) of the change relative to the current year.
#'
#' The resulting dataframe includes, for each timestep, the current state, the previous year's state, the next state, and the target time (i.e., the number of years from the current year when the change will happen).
#'
#' We also flag the years where a state change occurs, defined as when the current state differs from the previous year's state. Additionally, we flag years where no change occurs within a 10-year period (no change examples).
#'
#' Finally, for each flagged example, we retrieve the state history, which includes the last three states and their corresponding residence times.
#'
#' The wrapper function processes the simulation data by splitting it into individual simulations. This ensures that simulations are handled separately to avoid mixing data. After processing, it combines the results of multiple simulations into a single dataset.
#'
#' @param x (file) Simulation data.
#'
#' @author Marc Grünig
#' Last modified: 16/03/2023.
#' 
#' 
#' 
#



### function
create_aug_examples <- function(x){
  
  # create empty vectors
  residence_time <- c()
  target_state <- c()
  previous_year_state <- c()
  target_time <- c()
  
  # for each row in the simulation dataset we check when the state changes
  res <- 1
  
  for(n in 1:nrow(x)){
    
    # the residence time increases when the vegetation state in the year before was the same. otherwise it's 1.
    res <- ifelse(n == 1, 1,
                  ifelse(as.character(x[n, "svd_state"]) == as.character(x[n-1, "svd_state"]),
                         res + 1, 1))
    
    # the next step is always +9 years unless the year is at the end
    next_step <- ifelse(n < nrow(x) - 9, n + 9, nrow(x))
    
    # obtain current state and the state at the next step
    current_state <- as.character(x[n, "svd_state"])
    next_states <- x[n:next_step, ]
    
    # residence time
    residence_time <- c(residence_time, res)
    
    # if the next state is not the same as the current, we check where is the timestep where it changes
    change_in <- min(which(next_states[, "svd_state"] != current_state)) - 1
    
    # if it's the same state, we put 10 years
    change_in <- ifelse(is.infinite(change_in), 10, change_in)
    
    # if the end is not covering the residence time, omit the data
    change_in <- ifelse(n > (nrow(x)-9) & change_in > (nrow(x) - n), NA, change_in)
    
    # add change horizon to the vector
    target_time <- c(target_time, change_in)
    
    # extract the previous year
    prev_year <- ifelse(n > 1, n-1, n)
    previous_year_state <- as.vector(c(previous_year_state,
                                       x[prev_year, "svd_state"]))
    
    # extract the target state
    target <- ifelse(change_in == 10, current_state,
                     next_states[min(which(next_states[, "svd_state"] != current_state)),
                                 "svd_state"] )
    
    target_state <- c(target_state, target)
    
  }
  
  # combine
  x$previous_year_state <- previous_year_state
  x$residence_time <- residence_time
  x$target_state <- target_state
  x$target_time <- target_time
  
  
  # flag years with state changes
  flag <- ifelse(previous_year_state != x$svd_state, 1, 0)
  
  # flag also the years where nothing happened for 10 years
  for(f in 1:length(flag)){
    
    if(f < 10){next}
    flag[f] <- ifelse(sum(flag[(f-9):f]) == 0, 1, flag[f])
    
  }
  
  # now obtain all years where the state changes
  x$flag <- flag
  
  x <- na.omit(x)
  examples <- x[x$flag == 1,]
  
  hist_df <- as.data.frame(matrix(ncol = 6, nrow = nrow(examples)))
  colnames(hist_df) <- c("hist_state1", "hist_state2",
                         "hist_state3", "hist_time1",
                         "hist_time2", "hist_time3")
  
  #
  for(e in 1:nrow(examples)){
    
    if(e == 1){
      hist_state1 <- "missing" #as.character(examples[e, c("svd_state")])
      hist_state2 <- "missing" #hist_state1 
      hist_state3 <- "missing" #hist_state1
      
      hist_time1 <- 0 #as.character(examples[e, c("residence_time")])
      hist_time2 <- 0 #hist_time1 
      hist_time3 <- 0 #hist_time1
    }
    
    if(e == 2){
      hist_state1  <- as.character(examples[e-1, c("svd_state")])
      hist_state2  <- "missing" #hist_state1 
      hist_state3  <- "missing" #hist_state1
      
      hist_time1 <- as.character(examples[e-1, c("residence_time")])
      hist_time2 <- 0 #hist_time1 
      hist_time3 <- 0 #hist_time1
    }
    
    if(e == 3){
      hist_state1 <- as.character(examples[e-1, c("svd_state")])
      hist_state2 <- as.character(examples[e-2, c("svd_state")])
      hist_state3 <- "missing" #hist_state2
      
      hist_time1 <- as.character(examples[e-1, c("residence_time")])
      hist_time2 <- as.character(examples[e-2, c("residence_time")])
      hist_time3 <- 0 #hist_time2
    }
    
    if(e > 3){
      hist_state1  <- as.character(examples[e-1, c("svd_state")])
      hist_state2  <- as.character(examples[e-2, c("svd_state")])
      hist_state3  <- as.character(examples[e-3, c("svd_state")])
      
      hist_time1  <- as.numeric(examples[e-1, c("residence_time")])
      hist_time2  <- as.numeric(examples[e-2, c("residence_time")]) 
      hist_time3  <- as.numeric(examples[e-3, c("residence_time")]) 
    }
    
    row <- as.data.frame(t(cbind(unlist(c(hist_state1, hist_state2,
                                          hist_state3, hist_time1,
                                          hist_time2, hist_time3)))))
    
    colnames(row) <- c("hist_state1", "hist_state2",
                       "hist_state3", "hist_time1", 
                       "hist_time2", "hist_time3")
    
    hist_df <- rbind(hist_df, row)
    
  }
  
  hist_df <- hist_df %>% mutate(hist_time1 = as.numeric(hist_time1),
                                hist_time2 = as.numeric(hist_time2),
                                hist_time3 = as.numeric(hist_time3)) %>% drop_na()
  
  # if there are no examples, the history is not needed and we return the empty dataframe
  if(nrow(examples) == 0){examples_fin <- examples}else{
    examples_fin <- cbind(examples, hist_df)
    
    examples_fin <- examples_fin %>% filter(target_time != 10) %>% filter(target_time != 1)
    if(nrow(examples_fin) == 0){examples_fin <- examples_fin}else{
      
      for(i in 1:nrow(examples_fin)){
        random_year <- sample(1:(examples_fin[i, "target_time"]-1), 1)
        examples_fin[i, "residence_time"] <- examples_fin[i, "residence_time"] + random_year
        examples_fin[i, "target_time"] <- examples_fin[i, "target_time"] - random_year
      }
    }
    
    
  }
  
  return(examples_fin)
  
}

create_aug_examples <- cmpfun(create_aug_examples)



# wrapper function to obtain the simulation data and loop over them
examples_aug_data_fun <- function(sim_ids_sub){
  
  expl_data_list <- lapply(1:length(sim_ids_sub), function(i){
    
    tab_sim <- tab[simulationID == sim_ids_sub[i]]
    tab_sim_sub <- as.data.frame(tab_sim[, c("uniqueID", "simulationID", 
                                             "Year", "svd_state")])
    
    if(nrow(tab_sim_sub) == 0){return()}
    data.table(create_aug_examples(tab_sim_sub))
    
  })
  
  # remove empty tables
  expl_data_list <- Filter(function(x) dim(x)[1] > 0, expl_data_list)
  
  # bind together
  do.call(bind_rows, expl_data_list)
  
}

# compile function
examples_aug_data_fun <- cmpfun(examples_aug_data_fun)
