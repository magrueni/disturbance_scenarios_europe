# --------------------------------------------------------------
# Bark Beetle Module Scipt 03: Bark beetle spreading kernels
# --------------------------------------------------------------
# Author: Werner Rammer
# Date: 2025-01-07
# 
# Description:
# The script models the spread of bark beetle (BB) populations in a forest using a probabilistic kernel-based approach.
# It defines the BB spread probabilities across a spatial grid, simulates the spatial dynamics of infestation under different generational setups (bbgen),
# and evaluates the effect of varying a reproductive factor (bb_multiplier).
# The output includes raster maps of BB spread and summary metrics for further analysis.
#
# --------------------------------------------------------------


### libraries ---
library(raster)

path <- "/.../"


### processing ---
# prepare a cumulative distribution of the bark beetle kernel
x <- 1:257
probs <- exp(-((x/4.5)^2)/4/40.5)  
probs <- probs / sum(probs)
cprobs <- cumsum(probs)
cprobidx <- rep(0, 100000)
pidx <- 1
for (i in 1:100000) {
  if (i / 100000 >= cprobs[pidx]) {
    pidx <- pidx + 1
    print(paste(i, pidx))
  }
  cprobidx[i] <- pidx
}
plot(cprobidx)
table(cprobidx)

plot( probs, type="l")
sum( probs[1:50] )
sum( probs[51:150] )
sum( probs[151:257] )

sum(exp(-((x/4.5)^2)/4/40.5))


bbgens <- list( c(1),  # 1 gen
                c(99), # 1 + 1
                c(1,2), # 2
                c(1,100), # 2 + 1
                c(1,2,3)) # 4

runBBExperiment <- function( bbgen, bb_multiplier, n_start=10000 ) {
  print(paste("running spread for", paste(bbgen, collapse="-"), ", factor k:", bb_multiplier))
  field <- raster(xmn = -500, xmx=600, ymn=-500, ymx=600, res=100, vals=0)
  n_pack <- n_start
  
  for (rounds in bbgen) {
    if (rounds == 1 | rounds==99) {
      # start with n_pack random points
      posx <- runif(n_pack, min=0, max=100 )
      posy <- runif(n_pack, min=0, max=100 )
    } else {
      # for additional rounds the number of beetles increases
      posx <- nposx
      posy <- nposy
    }
    eff_multiplier <- bb_multiplier
    if (rounds>=99) {
      # sisterbrood, add halve of the offspring
      eff_multiplier <- bb_multiplier * 1.5
    }
    
    n_pack <- round(n_pack * eff_multiplier) # number of offsping that spreads
    dir <- runif(n_pack, min=0, max=2*pi) # sample random direction
    dists <- sample(cprobidx, size=n_pack, replace=T) # sample distance from kernel
    nposx <- posx + dists*sin(dir) # landing coordinates (starting coords posx/posy are reused)
    nposy <- posy + dists*cos(dir) 
    cidx <- cellFromXY(field, cbind(nposx, nposy))
    # count how often beetles land in each 100m cell
    tcid <- table(cidx)
    cidx <- as.integer(names(tcid))
    # update the playing field
    field[][cidx] <- field[][cidx] + tcid
    #plot(nposy ~ nposx, main=rounds, pch=".")
  }
  
  # return the raster
  field[][which.max(field[])] <- 0 # remove source pixel
  field[] <- field[] / n_start # scale to number of tries
  
  field
  
}

# helper function to create a data frame from the raster
createDf_from_raster <- function(r, gen, k) {
  
  df <- as.data.frame(t(r[]))
  names(df) <- (  as.data.frame(coordinates(r)) %>% mutate(name=paste(x,y,sep="_")) )$name
  df <- cbind(gen=gen, k=k, df)
  df
}


### Run full experiment with multiple values for bb_multiplier (aka k)
result <- data.frame()

bb_multiplier <- 4

for (bb_multiplier in c(4,6,8,10,12,14,16)) {
  print(paste("Running full experiment with k=", bb_multiplier))
  fields <- stack()
  for (bbgen in bbgens) {
    if (bb_multiplier < 10) nstart <- 100000 else nstart <- 10000
    field <- runBBExperiment(bbgen, bb_multiplier, nstart)
    
    plot(log10(field))
    
    #sum(field[])
    #field[][ field[]>0 ]
    #field[][ field[]>1 ] <- 1
    #print( sum(field[]>0) )
    
    fields <- addLayer(fields, field)
    
  }
  
  names(fields) <- c("g1", "g1.5", "g2", "g2.5", "g3")
  for (x in 1:nlayers(fields)) {
    print( sum(fields[[x]][])  )
    result <- rbind(result, createDf_from_raster(fields[[x]], names(fields)[x], bb_multiplier) )
  }
  
}

names(result) <- paste0("v", names(result))
names(result) <- str_replace(names(result), "-", "m")
result$vgen <- as.numeric( str_replace(result$vgen, "g", "") )
## some checks

result_agg <- result %>% mutate(total = rowSums(across(vm450_550:v550_m450))) %>% dplyr::select(vgen, vk, total)
# limit to 1 (each cell counts max 1)
result_agg <- result %>% mutate(total = rowSums(across(vm450_550:v550_m450, ~ pmin(.x,1) ))) %>% dplyr::select(vgen, vk, total)

ggplot(result_agg, aes(x=vk, y=total, color=vgen)) + geom_line(lwd=1) + xlab("factor k")  


# save resulting data file
names(result[1:2]) <- c("gen", "k")
write_csv(result, paste0(path, "/bbtl_module/bb_kernels.csv"))
names(result)
dim(result)
summary(result)


##### Multiple random draws (kernel values > 1)
# kernel values > 1 are interpreted as "multiple" tries (with given susceptibility ps)

ps <- 0.4
kv <- 6
1 - (1 - ps)^kv

pmaps <- data.frame()
for (ps in (1:20)/20) {
  kv <- (1:100)/10
  est <- 1 - (1 - ps)^kv 
  pmaps <- rbind( pmaps, data.frame(ps = ps, kernel = kv, est = est) ) 
}
ggplot(pmaps, aes(x = kernel, y=est, color=ps, group=ps)) + geom_line() 

for ( vk in unique(result$vk))
  print(apply( result[ result$vk == vk,3:123],1, max ))
