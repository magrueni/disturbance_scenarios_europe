#' Check if all columns are there
#' 
#' @description In this function we calculate the dominant height for species we don't have the Kahn-factors.
#' Essentially, this function is based on the yield table. There we estimate a dominant height from mean and max height and fit a linear mixed model to the data with the actuall
#' dominant height. We use the model to calculate dominant height here.
#' 
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


height_max_mean_function <- function(x){
  
  #species <- strsplit(as.character(x$sp_state2), "(?<=.{4})", perl = TRUE)[[1]]
  x <- as.data.frame(x)
  x_split <- strsplit(as.character(x["sp_state2"]), "(?<=.{4})", perl = TRUE)[[1]]
  
  # if more than 1 species, get the proportion of the species for weighting the slope
  for(i in 1:length(x_split)){
      
      columns_with_string <- apply(x, 2, function(col) any(col == tolower(x_split[i])))
      result_columns <- na.omit(colnames(x)[columns_with_string])
      prop <- ifelse(result_columns == "Species1", as.numeric(x["Proportion1"])/100,
                         ifelse(result_columns == "Species2", as.numeric(x["Proportion2"])/100, 
                                ifelse(result_columns == "Species3", as.numeric(x["Proportion3"])/100,
                                       ifelse(result_columns == "Species4", as.numeric(x["Proportion4"])/100, as.numeric(x["Proportion5"])/100))))
      assign(paste0("prop_sp", i), prop)
  }
  
  
  # get the total of the species in case it doesn't add up to 100
  tot <- c()
  for(i in 1:length(x_split)){
    tot <- c(tot, get(paste0("prop_sp", i)))
  }
  tot <- sum(tot)
  
  
  # proportion of species in here
  for(i in 1:length(x_split)){
    
    prop <- get(paste0("prop_sp", i))/tot
    assign(paste0("prop_sp", i), prop)
    
  }

  # if the species doesnt have a specific model, get the global model
  if(length(x_split) == 0){sp <- "global"}else{
    sp <- ifelse(tolower(x_split) %in% c("bepe", "frex", "potr", "cabe", "poca", "rops", "algl"), "laub", tolower(x_split)) 
    sp <- ifelse(!sp %in% c("piab", "abal", "lade", "pisy", "fasy", "quro", "laub"), "global", sp)
  }

  
  # if there are two species, take weighted model coeffictients
  if(length(sp) > 1){
    
    intercept_sp <- c()
    slope_sp <- c()
    for(i in 1:length(sp)){
     
      coefs_sp <- all_coefs_domheightmodel[all_coefs_domheightmodel$species == sp[i], ]
      prop_sp <- get(paste0("prop_sp", i))
      intercept_sp <- c(intercept_sp, coefs_sp[, "intercept"] * prop_sp)
      slope_sp <- c(slope_sp, coefs_sp[, "slope"] * prop_sp)
      
    }
    
    intercept_sp <- sum(intercept_sp)
    slope_sp <- sum(slope_sp)
    
  }else{
    
    coefs_sp <- all_coefs_domheightmodel[all_coefs_domheightmodel$species %in% sp, ]
    intercept_sp <- coefs_sp[, "intercept"]
    slope_sp <- coefs_sp[, "slope"]
    
  }
  
  MeanHeight <- as.numeric(x["MeanHeight"])
  MaxHeight <- as.numeric(x["MaxHeight"])
  
  if(is.na(MeanHeight)){
    MinHeight <- as.numeric(x["MinHeight"])
    MeanHeight <- (MaxHeight + MinHeight)/2
  }
  
  hdom_est <- intercept_sp + ((as.numeric(MaxHeight) + as.numeric(MeanHeight))/2) * slope_sp
  
  # everything that is below a factor of 0.8 we cut to 0.8
  hdom_est <- ifelse(hdom_est/MaxHeight < 0.8, MaxHeight * 0.8, hdom_est)
  hdom_est <- ifelse(hdom_est > MaxHeight, MaxHeight, hdom_est)
  
  return(as.numeric(hdom_est))
}

height_max_mean_function <- cmpfun(height_max_mean_function)




### MEAN HEIGHT ONLY FUNCTION ------------------------------------------------------------------------------------------
height_mean_only_function <- function(x){
  
  # x <- as.vector(x)
  x <- as.data.frame(x)
  x_split <- strsplit(as.character(x["sp_state2"]), "(?<=.{4})", perl = TRUE)[[1]]
   
  # if more than 1 species, get the proportion of the species for weighting the slope
  for(i in 1:length(x_split)){
       
     columns_with_string <- apply(x, 2, function(col) any(col == tolower(x_split[i])))
     result_columns <- na.omit(colnames(x)[columns_with_string])
     prop <- ifelse(result_columns == "Species1", as.numeric(x["Proportion1"])/100,
                    ifelse(result_columns == "Species2", as.numeric(x["Proportion2"])/100, 
                           ifelse(result_columns == "Species3", as.numeric(x["Proportion3"])/100,
                                  ifelse(result_columns == "Species4", as.numeric(x["Proportion4"])/100, as.numeric(x["Proportion5"])/100))))
     assign(paste0("prop_sp", i), prop)
   }
  
  # get the total of the species in case it doesn't add up to 100
  tot <- c()
  for(i in 1:length(x_split)){
     tot <- c(tot, get(paste0("prop_sp", i)))
   }
  tot <- sum(tot)
   
   
  # proportion of species in here
  for(i in 1:length(x_split)){
     
     prop <- get(paste0("prop_sp", i))/tot
     assign(paste0("prop_sp", i), prop)
     
   }
  
   
  if(length(x_split) == 0){
    sp <- "global"
  }else{
      sp <- ifelse(tolower(x_split) %in% all_coefs_maxheightmodel[, "species"], tolower(x_split), "global")
    }

  
  # if there are two species, take average of model coeffictients
  if(length(sp) > 1){
    
    intercept_sp <- c()
    slope_sp <- c()
    for(i in 1:length(sp)){
      
      coefs_sp <- all_coefs_maxheightmodel[all_coefs_maxheightmodel$species == tolower(sp[i]), ]
      prop_sp <- get(paste0("prop_sp", i))
      intercept_sp <- c(intercept_sp, coefs_sp[, "intercept"] * prop_sp)
      slope_sp <- c(slope_sp, coefs_sp[, "slope"] * prop_sp)
      
    }
    
    intercept_sp <- sum(intercept_sp)
    slope_sp <- sum(slope_sp)
    
  }else{
    
    coefs_sp <- all_coefs_maxheightmodel %>% filter(species %in% sp)
    intercept_sp <- mean(coefs_sp[, "intercept"])
    slope_sp <- mean(coefs_sp[, "slope"])
    
  }
  
  MeanHeight <- as.numeric(x["MeanHeight"])
  H_max <- intercept_sp + as.numeric(MeanHeight) * slope_sp
  #H_max <- ifelse(H_max/MeanHeight < 1.25, MeanHeight * 1.25, H_max)
  H_max <- ifelse(H_max < MeanHeight, MeanHeight, H_max)
  
  
  # if there is a dominant species, we do the kahn factors transformation
  x_split <- strsplit(as.character(x["sp_state2"]), "(?<=.{4})", perl = TRUE)[[1]]
  if(str_detect(x_split[1],"[[:upper:]]")){

    hdom_est1 <- ifelse(x_split[1] %in% c("PIAB", "ABAL", "LADE"), H_max / 1.138,
                       ifelse(x_split[1] %in% c("PISY"), H_max / 1.18,
                              ifelse(x_split[1] %in% c("FASY"), H_max / 1.132,
                                     ifelse(x_split[1] %in% c("QURO", "QUCE", "QUCO", "QUFR", "QUIL", "QUPE", "QUPU", "QUSU", "QUTR"),
                                            H_max / 1.184, H_max / 1.14))))


    x_split <- x_split[-1]
  }
    
  # if not, we calc the dominant height as before 
  if(length(x_split) == 0){sp <- "global"}else{
      sp <- ifelse(tolower(x_split) %in% c("bepe", "frex", "potr", "cabe", "poca", "rops", "algl"), "laub", tolower(x_split)) 
      sp <- ifelse(!sp %in% c("piab", "abal", "lade", "pisy", "fasy", "quro", "laub"), "global", sp)
  }
  
  # if there are two species, take average of model coeffictients
  if(length(sp) > 1){
      
    intercept_sp <- c()
    slope_sp <- c()
    for(i in 1:length(sp)){
        
      coefs_sp <- all_coefs_domheightmodel[all_coefs_domheightmodel$species == sp[i], ]
      prop_sp <- get(paste0("prop_sp", i))
      intercept_sp <- c(intercept_sp, coefs_sp[, "intercept"] * prop_sp)
      slope_sp <- c(slope_sp, coefs_sp[, "slope"] * prop_sp)
        
    }
      
    intercept_sp <- sum(intercept_sp)
    slope_sp <- sum(slope_sp)
      
  }else{
      
    coefs_sp <- all_coefs_domheightmodel %>% filter(species %in% sp)
    intercept_sp <- mean(coefs_sp[, "intercept"])
    slope_sp <- mean(coefs_sp[, "slope"])
      
  }
    
  hdom_est2 <- intercept_sp + ((as.numeric(H_max) + as.numeric(MeanHeight))/2) * slope_sp
  
  x_split <- strsplit(as.character(x["sp_state2"]), "(?<=.{4})", perl = TRUE)[[1]]
  if(str_detect(x_split[1],"[[:upper:]]")){
    hdom_est <- (hdom_est1 * prop_sp1 + hdom_est2 * (1-prop_sp1))
  }else{
    hdom_est <- hdom_est2
  }
  
  # everything that is below a factor of 0.8 we cut to 0.8
  hdom_est <- ifelse(hdom_est/H_max < 0.8, H_max * 0.8, hdom_est)
  hdom_est <- ifelse(hdom_est > H_max, H_max, hdom_est)

  return(as.numeric(hdom_est))
  
}

height_mean_only_function_cmp <- cmpfun(height_mean_only_function)
