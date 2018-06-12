### This Code does not run simulations and is only for the analysis of the model from agout-brazilnut.R

###=======================================================================
### Parameters
###=======================================================================
highHarvestFecundity <- 0.85 # Multiplier for fecundity rate for Adult trees under HIGH harvest
highHarvestSurvival <- 0.9 # Multiplier for survival rate of Adult trees under HIGH harvest
agoutiGrowth <- 1.1 # Growth rate of Agoutis in logistic model

lowHunting <- 0.1 # Percentage of agoutis hunted during LOW hunting
highHunting <- 0.25 # Percentage of agoutis hunted during HIGH hunting

seedlingCapacity <- 5000 # Carrying Capacities. Tree capacities don't matter as much.
saplingCapacity <- 500
adultCapacity <- 100
agoutiCapacity <- 5200
perc <- 0.9
ylimit <- 1 # For plotting

seedlingInit <- perc * seedlingCapacity#5000 # Initial Populations
saplingInit <- perc * saplingCapacity#500
adultInit <- perc * adultCapacity#100
agoutiInit <- perc * agoutiCapacity#5000

seedlingInit <- 5000 # Initial Populations
saplingInit <- 500
adultInit <- 100
agoutiInit <- 5000

m <- 0.05    # m is the desired proportion at which sigmoid(m) = m . Ideally it is small (~0.01-0.05).
agouti_to_PlantSteepness <- -(log(1-m)-log(m))/((m-0.5)*agoutiCapacity) # Steepness needed for sigmoid(m) = m
plant_to_AgoutiSteepness <- -(log(1-m)-log(m))/((m-0.5)*adultCapacity)  # Steepness needed for sigmoid(m) = m
#This formula above is derived from logistic function with "x = m*CAP" , "x0 = .5*CAP" , "y = m" , and solving for k. (CAP = carrying capacity)

time_end <- 500 # Length of simulation in years

# For linear functional form
m <- 1/agoutiCapacity
b <- 0
#=========================================================================



###========================================================================
### The Original 17-Stage Matrix from Zuidema and high-harvest multiplier
###========================================================================
plant_S_mat <- matrix( 0, nrow = 17, ncol = 17)
diag(plant_S_mat) <- c(0.455, 0.587, 0.78, 0.821, 0.941, 0.938, 0.961, 0.946, 0.94, 0.937, 0.936, 0.966, 0.968, 0.971, 0.965, 0.967, 0.985)
plant_S_mat[cbind(2:17,1:16)] <- matrix(c(0.091, 0.147, 0.134, 0.167, 0.044, 0.047, 0.034, 0.049, 0.055, 0.058, 0.059, 0.029, 0.027, 0.024, 0.020, 0.018))
plant_S_mat[1,12:17] <- matrix( c(12.3, 14.6, 16.9, 19.3, 22.3, 26.6) )

high_harv <- matrix(1, nrow = 17, ncol = 17)
high_harv[1,12:17] <- highHarvestFecundity    # Multiplier for fecundity rate for Adult trees
high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees
#===========================================================================



###====================================================================
### Plant and animal dynamics under harvest
###====================================================================
# Plant
plant_mat_low <- plant_S_mat
plant_mat_high <- plant_S_mat * high_harv

sigmoid <- function(k, x0, x) 
{
  1/(1+exp(-k*(x-x0))) #k: steepness #x0 = midpoint
} # plant_animal <- matrix(c(1,1,sigmoid(1, K_animal*0.5, N),1,1,1,1,1,1), nrow=3, byrow=T)
# Animal

linear <- function(m, x, b)
{
  y <- m*x + b
  return(y)
}

LogisticGrowthHunt <- function(R, N, K, H, p) 
{ # p is how the plant affects carrying capacity of agoutis (from 0 to 1)
  Nnext <- R*N*(1-N/(K*(p))) - H*N + N
  return(Nnext)
} # some proportion of N_t are harvested
###=====================================================================



###====================================================================
### Harvest markov chain
###====================================================================
library(markovchain)
# Specifying the markov chain
statesNames = c("low","high")
mcHarvest <- new("markovchain", states = statesNames, 
                 transitionMatrix = matrix(data = c(0.2, 0.8, 0.8, 0.2), byrow = TRUE, 
                                           nrow = 2, dimnames=list(statesNames,statesNames)), name="Harvest")
# Simulating a discrete time process for harvest
set.seed(100)
harvest_seq <- rmarkovchain(n=time_end, object = mcHarvest, t0="low")
head(harvest_seq)
###====================================================================




###===========================================================================
### Eigenvalues vs. Population of Agoutis
###===========================================================================
eigvals_vs_agouti <- function()
{
  
  domEigenvals = c()
  num <- 1
  interval <- 10 # 1720 is point at which eig ~= 1
  for (i in seq(0,5200,interval))
  {
    plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
    plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, i)
    #  plant_animal_mat[1,12:17] <- linear(m, i, b)
    plant_matrix <- plant_animal_mat * plant_mat_low
    
    eigenvals <- eigen(plant_matrix)$values
    
    dominant <- 0
    for (i in eigenvals)
    {
      dominant <- max(sqrt(Re(i*Conj(i))), dominant)
    }
    domEigenvals[num] <- dominant
    
    num <- num + 1
  }
  
  plot(seq(0,5200,interval), domEigenvals)
  # 1714 is point at which eig ~= 1, at m = 0.0001
  # 0.4265611 times the carrying capacity (p = 0.42806) is where dN/dt ~= 0 when N = 1714. N is pop of agouti
  # Since p is bounded between 0.9 and 1, it is impossible for  system to die
  
  # 103 is the point at which eig ~= 1, at m = 0.05
  # 0.02563348 times the carrying capacity (of agoutis) is inflection point
  # Since p is bounded between 0.9 and 1, it is impossible for system to die
  
  for (i in 1:length(domEigenvals))
  {
    if (domEigenvals[i] >= 1 && domEigenvals[i] <= 1.0001)
    {
      print(i)
    }
  }
}

#eigvals_vs_agouti()

###===========================================================================



###===========================================================================
### Eigenvalues vs. Percent Harvest
###===========================================================================
# This is only for harvest relating to the fecundity of the trees.
eigvals_vs_harvest <- function()
{

  domEigenvals = c()
  num <- 1
  interval <- .25
  
  for (i in seq(0,100,interval))
  {
    
    high_harv <- matrix(1, nrow = 17, ncol = 17)
    high_harv[1,12:17] <- 1-i/100    # Multiplier for fecundity rate for Adult trees
    high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees
    plant_mat_high <- plant_S_mat * high_harv
    
    eigenvals <- eigen(plant_mat_high)$values
    
    dominant <- 0
    for (i in eigenvals)
    {
      dominant <- max(sqrt(Re(i*Conj(i))), dominant)
    }
    domEigenvals[num] <- dominant
    
    num <- num + 1
  }

  plot(seq(0,100,interval), domEigenvals, xlab ='Percent Harvest', ylab = 'Dominant Eigenvalue')
}

#eigvals_vs_harvest()

###===========================================================================



###===========================================================================
### Sensitivity Matrix
###===========================================================================
sensitivity_matrix <- function(matrix1)
{
  eigvals <- eigen(matrix1)$values
  eigvecs <- eigen(matrix1)$vectors
  W <- eigvecs
  
  sensMat <- matrix(0, nrow = nrow(matrix1), ncol = ncol(matrix1))
  diag(sensMat) <- c(eigvals)
  
  count <- 0
  for (i in eigvals)
  {
    count <- count + 1
    
    if (dominant <= sqrt(Re(i*Conj(i))))
    {
      dominant <- max(sqrt(Re(i*Conj(i))), dominant)
      finalCount <- count
    }
  }
  
  V <- Conj(solve(W)) # solve finds inverse matrix
  w <- matrix(W[,finalCount])
  v <- matrix(Re(V[finalCount,]))
  
  sensMat <- Re(v %*% t(w))
  
  elasMat <- (sensMat * matrix1) / dominant
  
  stableAgeDistr <- ((matrix(Re(eigvecs[,finalCount]) * -1)) / norm(matrix(Re(eigvecs[,finalCount]) * -1))) # Stable age distribution
  
  return(sensMat)
}

#sensMat <- sensitivity_matrix(plant_S_mat)

###===========================================================================




###===========================================================================
###===========================================================================
###===========================================================================
### Sensitivity Analysis
###===========================================================================
###===========================================================================
###===========================================================================



###====================================================================
### Sensitivity Analysis of Percent Harvested (Under purely High Harvests)
###====================================================================
sens_percent_harvest <- function()
{
  plant_mat <- matrix(0, nrow = 17)
  
  plant_mat[1:4] <- seedlingInit/4   #Setting initial population of seedlings
  plant_mat[5:11] <- saplingInit/7   #Setting initial population of saplings
  plant_mat[12:17] <- adultInit/6  #Setting initial population of adult trees
  agouti_vec <- c(agoutiInit) # Initializing the vector containing agouti pop at each timestep
  
  plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
  
  plot(1, xlab="", ylab="Population size/Max size", col="brown", ylim=c(0,100), type="l",xlim=c(1,time_end),xaxs="i")
  
  for (j in seq(0,100,10)) 
  {
    
    plant_mat <- matrix(0, nrow = 17)
    
    plant_mat[1:4] <- seedlingInit/4   #Setting initial population of seedlings
    plant_mat[5:11] <- saplingInit/7   #Setting initial population of saplings
    plant_mat[12:17] <- adultInit/6  #Setting initial population of adult trees
    agouti_vec <- c(agoutiInit) # Initializing the vector containing agouti pop at each timestep
    
    plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
    
    high_harv <- matrix(1, nrow = 17, ncol = 17)
    high_harv[1,12:17] <- 1-j/100    # Multiplier for fecundity rate for Adult trees
    plant_mat_low <- plant_S_mat
    plant_mat_high <- plant_S_mat * high_harv
    
    for (i in 1:time_end) 
    {
      h_i <- harvest_seq[i]
      
      if (h_i == "low") 
      {
        #pmat <- plant_mat_low
        #h_off <- lowHunting
        pmat <- plant_mat_high
        h_off <- highHunting
      } 
      
      else 
      {
        pmat <- plant_mat_high
        h_off <- highHunting
      }
      
      p <- sigmoid(plant_to_AgoutiSteepness, 50, sum(plant_mat[12:17]))*.1 + 0.9 # bounded between 0.9 and 1.0.... k was 0.1
      agouti_vec[(i+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(i)],agoutiCapacity,h_off, p)
      plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
      plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(i+1)]) # k was 0.0025
      #  plant_animal_mat[1,12:17] <- linear(m, agouti_vec[(i+1)], b) # A different functional form
      plant_mat <- matrix( c((plant_animal_mat * pmat) %*% plant_mat))
      
      #Summing the stages into 3 categories for better plotting
      plant_mat_sum <- c( sum(plant_mat[1:4]), sum(plant_mat[5:11]), sum(plant_mat[12:17])) 
      plant_all <- cbind(plant_all, plant_mat_sum)
    }
    
    lines(plant_all[3,]/adultCapacity, col='forestgreen')
  }
  
  par(mar=c(5,4,1,1),oma=c(0,0,0,0))
  #plot(1, xlab="", ylab="Population size/Max size", col="brown", ylim=c(0,10), type="l",xlim=c(1,time_end),xaxs="i")
  #lines(plant_all[1,]/seedlingCapacity, col="forestgreen")
  #lines(plant_all[2,]/saplingCapacity, col="turquoise3")
  #lines(plant_all[3,]/adultCapacity, col="orange")
  #legend("bottomleft", c("Animal","Adult","Sapling","Seedling"),col=c("brown","orange","turquoise3","forestgreen"),lty=c(1,1,1,1), bty="n",ncol=2)
  #mtext("Time step", 1, line=1.85, at=25, col="black")
  #axis(1,1:time_end,labels=toupper( substr(harvest_seq,1,1) ),line=2,col=NA,col.ticks=NA,col.axis="black", cex.axis=0.65)
  #mtext("Harvest:",1,line=3,at=-2.5,col="black")
}

#sens_percent_harvest()

###===========================================================================



###====================================================================
### Sensitivity Analysis of Harvest Survival Probability of Adults (Under purely High Harvests)
###====================================================================
sens_surv_prob <- function()
{
  plant_mat <- matrix(0, nrow = 17)
  
  plant_mat[1:4] <- seedlingInit/4   #Setting initial population of seedlings
  plant_mat[5:11] <- saplingInit/7   #Setting initial population of saplings
  plant_mat[12:17] <- adultInit/6  #Setting initial population of adult trees
  agouti_vec <- c(agoutiInit) # Initializing the vector containing agouti pop at each timestep
  
  plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
  
  plot(1, xlab="", ylab="Population size/Max size", col="brown", ylim=c(0,4), type="l",xlim=c(1,time_end),xaxs="i")
  
  for (j in seq(0,100,10)) 
  {
    
    plant_mat <- matrix(0, nrow = 17)
    
    plant_mat[1:4] <- seedlingInit/4   #Setting initial population of seedlings
    plant_mat[5:11] <- saplingInit/7   #Setting initial population of saplings
    plant_mat[12:17] <- adultInit/6  #Setting initial population of adult trees
    agouti_vec <- c(agoutiInit) # Initializing the vector containing agouti pop at each timestep
    
    plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
    
    high_harv <- matrix(1, nrow = 17, ncol = 17)
    high_harv[1,12:17] <- highHarvestFecundity    # Multiplier for fecundity rate for Adult trees
    high_harv[cbind(12:17,12:17)] <- 1-j/100 # Multiplier for survival rate of Adult trees
    plant_mat_low <- plant_S_mat
    plant_mat_high <- plant_S_mat * high_harv
    
    for (i in 1:time_end) 
    {
      h_i <- harvest_seq[i]
      
      if (h_i == "low") 
      {
        #pmat <- plant_mat_low
        #h_off <- lowHunting
        pmat <- plant_mat_high
        h_off <- highHunting
      } 
      
      else 
      {
        pmat <- plant_mat_high
        h_off <- highHunting
      }
      
      p <- sigmoid(plant_to_AgoutiSteepness, 50, sum(plant_mat[12:17]))*.1 + 0.9 # bounded between 0.9 and 1.0.... k was 0.1
      agouti_vec[(i+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(i)],agoutiCapacity,h_off, p)
      plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
      plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(i+1)]) # k was 0.0025
      #  plant_animal_mat[1,12:17] <- linear(m, agouti_vec[(i+1)], b) # A different functional form
      plant_mat <- matrix( c((plant_animal_mat * pmat) %*% plant_mat))
      
      #Summing the stages into 3 categories for better plotting
      plant_mat_sum <- c( sum(plant_mat[1:4]), sum(plant_mat[5:11]), sum(plant_mat[12:17])) 
      plant_all <- cbind(plant_all, plant_mat_sum)
    }
    
    lines(plant_all[3,]/adultCapacity, col='forestgreen')
  }
  
  par(mar=c(5,4,1,1),oma=c(0,0,0,0))
  #plot(1, xlab="", ylab="Population size/Max size", col="brown", ylim=c(0,10), type="l",xlim=c(1,time_end),xaxs="i")
  #lines(plant_all[1,]/seedlingCapacity, col="forestgreen")
  #lines(plant_all[2,]/saplingCapacity, col="turquoise3")
  #lines(plant_all[3,]/adultCapacity, col="orange")
  #legend("bottomleft", c("Animal","Adult","Sapling","Seedling"),col=c("brown","orange","turquoise3","forestgreen"),lty=c(1,1,1,1), bty="n",ncol=2)
  #mtext("Time step", 1, line=1.85, at=25, col="black")
  #axis(1,1:time_end,labels=toupper( substr(harvest_seq,1,1) ),line=2,col=NA,col.ticks=NA,col.axis="black", cex.axis=0.65)
  #mtext("Harvest:",1,line=3,at=-2.5,col="black")
}

#sens_surv_prob()

###===========================================================================



###====================================================================
### Sensitivity Analysis of Agouti Hunting rate (Under purely High Harvests)
###====================================================================
sens_agouti_hunt <- function()
{
  plant_mat <- matrix(0, nrow = 17)
  
  plant_mat[1:4] <- seedlingInit/4   #Setting initial population of seedlings
  plant_mat[5:11] <- saplingInit/7   #Setting initial population of saplings
  plant_mat[12:17] <- adultInit/6  #Setting initial population of adult trees
  agouti_vec <- c(agoutiInit) # Initializing the vector containing agouti pop at each timestep
  
  plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
  
  plot(1, xlab="", ylab="Population size/Max size", col="brown", ylim=c(0,ylimit), type="l",xlim=c(1,time_end),xaxs="i")
  
  for (j in seq(0,100,10)) 
  {
    
    highHunting <- j/100
    
    plant_mat <- matrix(0, nrow = 17)
    
    plant_mat[1:4] <- seedlingInit/4   #Setting initial population of seedlings
    plant_mat[5:11] <- saplingInit/7   #Setting initial population of saplings
    plant_mat[12:17] <- adultInit/6  #Setting initial population of adult trees
    agouti_vec <- c(agoutiInit) # Initializing the vector containing agouti pop at each timestep
    
    plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
    
    for (i in 1:time_end) 
    {
      h_i <- harvest_seq[i]
      
      if (h_i == "low") 
      {
        #pmat <- plant_mat_low
        #h_off <- lowHunting
        pmat <- plant_mat_high
        h_off <- highHunting
      } 
      
      else 
      {
        pmat <- plant_mat_high
        h_off <- highHunting
      }
      
      p <- sigmoid(plant_to_AgoutiSteepness, 50, sum(plant_mat[12:17]))*.1 + 0.9 # bounded between 0.9 and 1.0.... k was 0.1
      agouti_vec[(i+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(i)],agoutiCapacity,h_off, p)
      plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
      plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(i+1)]) # k was 0.0025
      #  plant_animal_mat[1,12:17] <- linear(m, agouti_vec[(i+1)], b) # A different functional form
      plant_mat <- matrix( c((plant_animal_mat * pmat) %*% plant_mat))
      
      #Summing the stages into 3 categories for better plotting
      plant_mat_sum <- c( sum(plant_mat[1:4]), sum(plant_mat[5:11]), sum(plant_mat[12:17])) 
      plant_all <- cbind(plant_all, plant_mat_sum)
    }
    
    lines(agouti_vec/agoutiCapacity, col='orange') # Just looking at agouti and adult tree population
    lines(plant_all[3,]/adultCapacity, col='forestgreen')
#    lines(plant_all[2,]/saplingCapacity, col='turquoise3')
#    lines(plant_all[1,]/seedlingCapacity, col='brown')
    
  }
  
  par(mar=c(5,4,1,1),oma=c(0,0,0,0))
  #plot(1, xlab="", ylab="Population size/Max size", col="brown", ylim=c(0,10), type="l",xlim=c(1,time_end),xaxs="i")
  #lines(plant_all[1,]/seedlingCapacity, col="forestgreen")
  #lines(plant_all[2,]/saplingCapacity, col="turquoise3")
  #lines(plant_all[3,]/adultCapacity, col="orange")
  #legend("bottomleft", c("Animal","Adult","Sapling","Seedling"),col=c("brown","orange","turquoise3","forestgreen"),lty=c(1,1,1,1), bty="n",ncol=2)
  #mtext("Time step", 1, line=1.85, at=25, col="black")
  #axis(1,1:time_end,labels=toupper( substr(harvest_seq,1,1) ),line=2,col=NA,col.ticks=NA,col.axis="black", cex.axis=0.65)
  #mtext("Harvest:",1,line=3,at=-2.5,col="black")
}

#sens_agouti_hunt()

###===========================================================================



###====================================================================
### Sensitivity Analysis of Agouti Growth rate (Under purely High Harvests)
###====================================================================
sens_agouti_growth <- function()
{
  plant_mat <- matrix(0, nrow = 17)
  
  plant_mat[1:4] <- seedlingInit/4   #Setting initial population of seedlings
  plant_mat[5:11] <- saplingInit/7   #Setting initial population of saplings
  plant_mat[12:17] <- adultInit/6  #Setting initial population of adult trees
  agouti_vec <- c(agoutiInit) # Initializing the vector containing agouti pop at each timestep
  
  plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
  
  plot(1, xlab="", ylab="Population size/Max size", col="brown", ylim=c(0,4), type="l",xlim=c(1,time_end),xaxs="i")
  
  for (j in seq(0,100,10)) 
  {
    
    agoutiGrowth <- j/100
    
    plant_mat <- matrix(0, nrow = 17)
    
    plant_mat[1:4] <- seedlingInit/4   #Setting initial population of seedlings
    plant_mat[5:11] <- saplingInit/7   #Setting initial population of saplings
    plant_mat[12:17] <- adultInit/6  #Setting initial population of adult trees
    agouti_vec <- c(agoutiInit) # Initializing the vector containing agouti pop at each timestep
    
    plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
    
    for (i in 1:time_end) 
    {
      h_i <- harvest_seq[i]
      
      if (h_i == "low") 
      {
        #pmat <- plant_mat_low
        #h_off <- lowHunting
        pmat <- plant_mat_high
        h_off <- highHunting
      } 
      
      else 
      {
        pmat <- plant_mat_high
        h_off <- highHunting
      }
      
      p <- sigmoid(plant_to_AgoutiSteepness, 50, sum(plant_mat[12:17]))*.1 + 0.9 # bounded between 0.9 and 1.0.... k was 0.1
      agouti_vec[(i+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(i)],agoutiCapacity,h_off, p)
      plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
      plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(i+1)]) # k was 0.0025
      #  plant_animal_mat[1,12:17] <- linear(m, agouti_vec[(i+1)], b) # A different functional form
      plant_mat <- matrix( c((plant_animal_mat * pmat) %*% plant_mat))
      
      #Summing the stages into 3 categories for better plotting
      plant_mat_sum <- c( sum(plant_mat[1:4]), sum(plant_mat[5:11]), sum(plant_mat[12:17])) 
      plant_all <- cbind(plant_all, plant_mat_sum)
    }
    
    lines(agouti_vec/agoutiCapacity, col='orange')
    lines(plant_all[3,]/adultCapacity, col='forestgreen')
  }
  
  par(mar=c(5,4,1,1),oma=c(0,0,0,0))
  #plot(1, xlab="", ylab="Population size/Max size", col="brown", ylim=c(0,10), type="l",xlim=c(1,time_end),xaxs="i")
  #lines(plant_all[1,]/seedlingCapacity, col="forestgreen")
  #lines(plant_all[2,]/saplingCapacity, col="turquoise3")
  #lines(plant_all[3,]/adultCapacity, col="orange")
  #legend("bottomleft", c("Animal","Adult","Sapling","Seedling"),col=c("brown","orange","turquoise3","forestgreen"),lty=c(1,1,1,1), bty="n",ncol=2)
  #mtext("Time step", 1, line=1.85, at=25, col="black")
  #axis(1,1:time_end,labels=toupper( substr(harvest_seq,1,1) ),line=2,col=NA,col.ticks=NA,col.axis="black", cex.axis=0.65)
  #mtext("Harvest:",1,line=3,at=-2.5,col="black")
}

#sens_agouti_growth()

###===========================================================================




###===========================================================================
### PopBio Analysis
###===========================================================================

# Stochastic growth rate
lambda_sim <- function()
{
  brazilNut <- list(low=plant_mat_low, high=plant_mat_high)
  
  ### Part 1: Inspect the leading eigenvalue for each individual matrix
  lambdas <- lapply(brazilNut, lambda) # "List apply": similar in spirit to python's mapply; applies a function over each element in a list or vector
  lambdas # Note that 1985 and 1987 are "bad" years
  
  sgr1 <- stoch.growth.rate(brazilNut, prob=rep(0.50, 2)) # equal probabilities of all years
  sgr1 # note that $approx returns Tuljapurkar's approximation of log-lambda and $sim gives you a simulated estimate
  # To extract the true lambdas:
  exp(sgr1$approx)
  exp(sgr1$sim)
  
  return(exp(sgr1$sim))
}

growth_rate1 <- lambda_sim() #Does not take into account the agoutis' effect of plants (assumes max agouti population)



# Improved Stochastic Growth Rate. (Takes agoutis into account)
maxt <- 50000
brazilNut <- list(low=plant_mat_low, high=plant_mat_high)

stoch_growth <- function()
{
  r <- numeric(maxt)
  
  plant_mat <- matrix(0, nrow = 17)
  
  plant_mat[1:4] <- seedlingInit/4   #Setting initial population of seedlings
  plant_mat[5:11] <- saplingInit/7   #Setting initial population of saplings
  plant_mat[12:17] <- adultInit/6  #Setting initial population of adult trees
  agouti_vec <- c(agoutiInit) # Initializing the vector containing agouti pop at each timestep
  
  plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
  
  plant_mat <- plant_mat / sum(plant_mat)
  
  statesNames = c("low","high")
  mcHarvest <- new("markovchain", states = statesNames, 
                   transitionMatrix = matrix(data = c(0.2, 0.8, 0.8, 0.2), byrow = TRUE, 
                                             nrow = 2, dimnames=list(statesNames,statesNames)), name="Harvest")
  # Simulating a discrete time process for harvest
  set.seed(100)
  harvest_seq <- rmarkovchain(n=maxt, object = mcHarvest, t0="low")
  
  for (i in 1:maxt)
  {
    h_i <- harvest_seq[i]
    
    if (h_i == "low") 
    {
      pmat <- plant_mat_low
      h_off <- lowHunting
    } 
    
    else 
    {
      pmat <- plant_mat_high
      h_off <- highHunting
    }
    
    p <- sigmoid(plant_to_AgoutiSteepness, 50, sum(plant_mat[12:17]))*.1 + 0.9 # bounded between 0.9 and 1.0.... k was 0.1
    agouti_vec[(i+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(i)],agoutiCapacity,h_off, p)
    plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
    plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(i+1)]) # k was 0.0025
    #  plant_animal_mat[1,12:17] <- linear(m, agouti_vec[(i+1)], b) # A different functional form
    plant_mat <- matrix( c((plant_animal_mat * pmat) %*% plant_mat))
    
    #Summing the stages into 3 categories for better plotting
    plant_mat_sum <- c( sum(plant_mat[1:4]), sum(plant_mat[5:11]), sum(plant_mat[12:17])) 
    plant_all <- cbind(plant_all, plant_mat_sum)
    
    N <- sum(plant_mat)
    r[i] <- log(N)
    plant_mat <- plant_mat / N
  }
  
  loglambsim <- mean(r)
  
  return(loglambsim)
}

growth_rate <- exp(stoch_growth())
print(growth_rate)

###===========================================================================



###===========================================================================
### Population Viability Analysis
###===========================================================================

VertebratePVA <- function(reps) {
  
  nrun <- function() 
  {
    
    statesNames = c("low","high")
    mcHarvest <- new("markovchain", states = statesNames, 
                     transitionMatrix = matrix(data = c(0.2, 0.8, 0.8, 0.2), byrow = TRUE, 
                                               nrow = 2, dimnames=list(statesNames,statesNames)), name="Harvest")
#    set.seed(100)
    harvest_seq <- rmarkovchain(n=time_end, object = mcHarvest, t0="low")
    
    plant_mat <- matrix(0, nrow = 17)
    plant_mat[1:4] <- seedlingInit/4   #Setting initial population of seedlings
    plant_mat[5:11] <- saplingInit/7   #Setting initial population of saplings
    plant_mat[12:17] <- adultInit/6  #Setting initial population of adult trees
    agouti_vec <- c(agoutiInit) # Initializing the vector containing agouti pop at each timestep
    
    plant_all <- matrix( c(seedlingInit, saplingInit, adultInit) ) # This will contain the summed plant populations at ALL timesteps
    
    for (i in 1:time_end) 
    {
      h_i <- harvest_seq[i]
      
      if (h_i == "low") 
      {
        pmat <- plant_mat_low
        h_off <- lowHunting
      } 
      
      else 
      {
        pmat <- plant_mat_high
        h_off <- highHunting
      }
      
      p <- sigmoid(plant_to_AgoutiSteepness, 50, sum(plant_mat[12:17]))*.1 + 0.9 # bounded between 0.9 and 1.0.... k was 0.1
      agouti_vec[(i+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(i)],agoutiCapacity,h_off, p)
      plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
      plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(i+1)]) # k was 0.0025
      #  plant_animal_mat[1,12:17] <- linear(m, agouti_vec[(i+1)], b) # A different functional form
      plant_mat <- matrix( c((plant_animal_mat * pmat) %*% plant_mat))
      
      #Summing the stages into 3 categories for better plotting
      plant_mat_sum <- c( sum(plant_mat[1:4]), sum(plant_mat[5:11]), sum(plant_mat[12:17])) 
      plant_all <- cbind(plant_all, plant_mat_sum)
      
    }
    return(plant_all[3,])  # PVA for adult trees
#   return(agouti_vec)      #PVA for agoutis
  }
  
  PopMat <- replicate(reps,
                      nrun())
  
  
  return(list(Nmat=PopMat, # matrix of population sizes over time
              Nend=PopMat[time_end,] # final population size
              ))
}

PlotCloud <- function(simdata){
  # From Kevin Shoemaker: http://naes.unr.edu/shoemaker/teaching/NRES-470/LECTURE12.html#step_4:_simulate
  
  # Visualize population trends
  plot(c(1:nrow(simdata)),simdata[,1],col=gray(0.7),type="l",ylim=c(0,max(simdata)),xlab="Years",ylab="Abundance")
  
  for(r in 2:ncol(simdata)){
    lines(c(1:nrow(simdata)),simdata[,r],col=gray(0.7),type="l")
  }
}

### Running the functions
#PVAruns <- VertebratePVA(100)

### Visualize the different runs
#PlotCloud(PVAruns$Nmat)




###===========================================================================



###===========================================================================

#Perturbation Analysis ideas:
#   plot agouti population sequence on top of each other given different harvest rates
#   plot tree population sequence on top of each other give different harvest rates
#   Sensitivity of each matrix element as a matrix heatmap



