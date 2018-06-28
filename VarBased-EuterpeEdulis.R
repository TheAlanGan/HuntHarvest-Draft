# Variance-Based Sensitivity Analysis


###=======================================================================
### Parameters
###=======================================================================
#install.packages('sensitivity')
#install.packages('markovchain')

library(sensitivity)
library(boot)
library(lhs)

# highHarvestFecundity <- 0.85 # Multiplier for fecundity rate for Adult trees under HIGH harvest
# highHarvestSurvival <- 0.9 # Multiplier for survival rate of Adult trees under HIGH harvest
# animalGrowth <- 1.1 # Growth rate of animals in logistic model
# 
# lowHunting <- 0.1 # Percentage of animals hunted during LOW hunting
# highHunting <- 0.25 # Percentage of animals hunted during HIGH hunting

seedlingCapacity <- 24272 # Carrying Capacities. Tree capacities don't matter as much.
saplingCapacity <- 260
adultCapacity <- 492
animalCapacity <- 1
perc <- 0.9
ylimit <- 1 # For plotting

seedlingInit <- perc * seedlingCapacity#5000 # Initial Populations
saplingInit <- perc * saplingCapacity#500
adultInit <- perc * adultCapacity#100
animalInit <- perc * animalCapacity#5000

# m <- 0.05    # m is the desired proportion at which sigmoid(m) = m . Ideally it is small (~0.01-0.05).
# animal_to_PlantSteepness <- -(log(1-m)-log(m))/((m-0.5)*animalCapacity) # Steepness needed for sigmoid(m) = m
# plant_to_animalSteepness <- -(log(1-m)-log(m))/((m-0.5)*adultCapacity)  # Steepness needed for sigmoid(m) = m
#This formula above is derived from logistic function with "x = m*CAP" , "x0 = .5*CAP" , "y = m" , and solving for k. (CAP = carrying capacity)

time_end <- 500 # Length of simulation in years
maxt <- 500 # Length of simulation for stoch_growth function.

#brazilNut <- list(low=plant_mat_low, high=plant_mat_high)

# high_harv <- matrix(1, nrow = 17, ncol = 17)

##=======The Original 17-Stage Matrix from Zuidema and high-harvest multiplier
# plant_S_mat <- matrix( 0, nrow = 17, ncol = 17)
# diag(plant_S_mat) <- c(0.455, 0.587, 0.78, 0.821, 0.941, 0.938, 0.961, 0.946, 0.94, 0.937, 0.936, 0.966, 0.968, 0.971, 0.965, 0.967, 0.985)
# plant_S_mat[cbind(2:17,1:16)] <- matrix(c(0.091, 0.147, 0.134, 0.167, 0.044, 0.047, 0.034, 0.049, 0.055, 0.058, 0.059, 0.029, 0.027, 0.024, 0.020, 0.018))
# plant_S_mat[1,12:17] <- matrix( c(12.3, 14.6, 16.9, 19.3, 22.3, 26.6) )
# high_harv[1,12:17] <- highHarvestFecundity 
# high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees
# plant_mat_low <- plant_S_mat
# plant_mat_high <- plant_S_mat * high_harv


#===========================================================================
#============FUNCTIONS======================================================
#===========================================================================

sigmoid <- function(k, x0, x) 
{
  1/(1+exp(-k*(x-x0))) #k: steepness #x0 = midpoint
} 


linear <- function(m, x, b)
{
  y <- m*x + b
  return(y)
}


LogisticGrowthHunt<- function(R, N, K, H, p) 
{ # p is how the plant affects carrying capacity of animals (from 0 to 1)
  Nnext <- R*N*(1-N/(K*(p))) - H*N + N
  return(Nnext)
}


LogisticGrowthHuntRK<- function(R, N, K, H, p,s) 
{ # p is how the plant affects carrying capacity of animals (from 0 to 1)
  Nnext <- R*(s)*N*(1-N/((K*(p)))) - H*N + N
  return(Nnext)
} 


# Specifying the markov chain
library('markovchain')
markovChain<- function(tEnd){
  
  statesNames = c("low","high")
  mcHarvest <- new("markovchain", states = statesNames, 
                   transitionMatrix = matrix(data = c(0.2, 0.8, 0.8, 0.2), byrow = TRUE, 
                                             nrow = 2, dimnames=list(statesNames,statesNames)), name="Harvest")
  # Simulating a discrete time process for harvest
  #  set.seed(100)
  harvest_seq <- markovchain::rmarkovchain(n=tEnd, object = mcHarvest, t0="low")
  return(harvest_seq)
}

#harvest_seq <- markovChain(time_end)




###=====================================================================

# Parameter list is as follows:
# 1.  m (determines steepness of sigmoid) ... [0.03, 0.1]
# 2.  delta ... [0, 1]
# 3.  high harvest fecundity ... [0, 9] ... get better interval?
# 4.  high harvest adult survival ... [0, 9] ... get better interval?
# 5.  low hunting rate (herbivores) ... [0, 0.25]
# 6.  high hunting rate (herbivores) ... [0.25, 1]
# 7.  rMax (herbivore) ... [0.2, 1.1]
# 8.  Adult Carrying Capacity ... [50,150]
# 9.  
# 10. 
# 11. 
# 12. 

stoch_growth_sobol <- function(X) # X is matrix of parameters. Columns are each parameter. Rows are different parameter sets.
{
  # Function taking in matrix of parameters
  # Returns (nx1) matrix of stochastic growth rates
  
  numSamples <- nrow(X)
  stochGrowthRates <- numeric(numSamples)
  popAfterTime <- numeric(numSamples)
  
  for (i in 1:numSamples)
  {
    # Unpacking parameters from matrix
    m <- X[i, 1] * (0.1 - 0.03) + 0.03 # Mapping [0,1] to [0.03,0.1]
    delta <- X[i, 2] # No need for mapping since delta is in [0,1] already
    highHarvestFecundity <- X[i, 3] * (0.9 - 0) + 0 # No mapping?
    highHarvestSurvival <- X[i, 4] * (0.9 - 0) + 0 # Mapping to [0,0.9]
    lowHunting <- X[i, 5] * (0.25 - 0) + 0 # Mapping to [0,0.25]
    highHunting <- X[i, 6] * (1 - 0.25) + 0.25 # Mapping to [0.25,1]
    animalGrowth <- X[i, 7] * (1.1 - 0.2) + 0.2 # Mapping to [0.2,1.1]
    adultCapacity <- X[i, 8] * (150 - 50) + 50 # Mapping to [50,150]
    
    ### Using parameters ###
    # Setting low harvest matrix
    plant_S_mat <- matrix( 0, nrow = 4, ncol = 4)
    diag(plant_S_mat) <- c(0.5380,0.6769,0.5806,1)
    plant_S_mat[cbind(2:4,1:3)] <- matrix(c(0.0062,0.2615,0.413))
    plant_S_mat[1,4] <- matrix(c(54))
    
    # Setting high harvest matrix
    high_harv <- matrix(1, nrow = 4, ncol = 4)    
    high_harv[1,4] <- highHarvestFecundity 
    high_harv[4,4] <- highHarvestSurvival # Multiplier for survival rate of Adult trees
    
    plant_mat_low <- plant_S_mat
    plant_mat_high <- plant_S_mat * high_harv
    
    animal_to_PlantSteepness <- -(log(1-m)-log(m))/((m-0.5)*animalCapacity) # Steepness needed for sigmoid(m) = m
    plant_to_animalSteepness <- -(log(1-m)-log(m))/((m-0.5)*adultCapacity) 
    
    r <- numeric(maxt) # maxt is global constant
    plant_mat <- matrix(0, nrow = 4)
    plant_mat[1] <- seedlingInit   # Initial populations are constants
    plant_mat[2:3] <- saplingInit/2   # So they are not parameters
    plant_mat[4] <- adultInit 
    animal_vec <- c(animalInit)
    
    harvest_seq <- markovChain(maxt)
    
    # if (i == 19)
    # {
    #   print(plant_mat_low)
    # }
    
    N <- 0 # Pop of plants after time maxt
    
    # Running the simulation to find 'Stochastic Growth Rate'
    for (j in 1:maxt)
    {
      h_j <- harvest_seq[j]
      
      if (h_j == "low") 
      {
        pmat <- plant_mat_low
        h_off <- lowHunting
      } 
      
      else 
      {
        pmat <- plant_mat_high
        h_off <- highHunting
      }
      
      prevN <- sum(plant_mat)
      
      p <- sigmoid(plant_to_animalSteepness, adultCapacity/2, sum(plant_mat[4]))*delta + (1-delta)
      animal_vec[(j+1)] <- LogisticGrowthHunt(animalGrowth, animal_vec[(j)],animalCapacity,h_off, p)
      plant_animal_mat <- matrix(1, nrow = 4, ncol = 4)
      plant_animal_mat[1,4] <- sigmoid(animal_to_PlantSteepness, animalCapacity/2, animal_vec[(j+1)]) # k was 0.0025
      #  plant_animal_mat[1,12:17] <- linear(m, animal_vec[(j+1)], b) # A different functional form
      plant_mat <- matrix( c((plant_animal_mat * pmat) %*% plant_mat))
      
      N <- sum(plant_mat) 
      r[i] <- log(N / prevN) # Calculating Growth Rate
    }
    
    #    stochGrowthRates[i] <- exp(mean(r)) # Collect growth rates into column vector
    popAfterTime[i] <- N
    
  }
  #  return(stochGrowthRates)
  print(popAfterTime)
  return(popAfterTime)
  
}

###=====================================================================

# Regular Sobol
#n <- 100
#X1 <- data.frame(matrix(runif(8 * n), nrow = n))
#X2 <- data.frame(matrix(runif(8 * n), nrow = n))

#result <- sobol(model = stoch_growth_sobol, X1 = X1, X2 = X2, order = 1, nboot = 100)
#print(result)

#f <- stoch_growth_sobol(X1)


# LHS Sobol
lhs <- sobolroalhs(model = stoch_growth_sobol, factors = 8, N = 100, p = 1, order = 1, nboot = 100)
print(lhs)
plot(lhs)
