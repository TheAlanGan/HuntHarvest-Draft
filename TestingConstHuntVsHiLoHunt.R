# Monte-Carlo-type approach to testing the difference 

# Averaging Elasticities over randomly generated parameters
# using Latin hypercube sampling

library(lhs)
library(markovchain)

highHarvestFecundity <- 0.85 # Multiplier for fecundity rate for Adult trees under HIGH harvest
highHarvestSurvival <- 0.9 # Multiplier for survival rate of Adult trees under HIGH harvest
agoutiGrowth <- 1.1 # Growth rate of Agoutis in logistic model

lowHunting <- 0.1 # Percentage of agoutis hunted during LOW hunting
highHunting <- 0.25 # Percentage of agoutis hunted during HIGH hunting

delta <- 0.1
deltaPlant <- 1

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
#agoutiInit <- 5000


m <- 0.05    # m is the desired proportion at which sigmoid(m) = m . Ideally it is small (~0.01-0.05).
agouti_to_PlantSteepness <- -(log(1-m)-log(m))/((m-0.5)*agoutiCapacity) # Steepness needed for sigmoid(m) = m
plant_to_AgoutiSteepness <- -(log(1-m)-log(m))/((m-0.5)*adultCapacity)  # Steepness needed for sigmoid(m) = m
#This formula above is derived from logistic function with "x = m*CAP" , "x0 = .5*CAP" , "y = m" , and solving for k. (CAP = carrying capacity)

#time_end <- 8000 # Length of simulation in years
maxt <- 8000 #controls the markovChain

high_harv <- matrix(1, nrow = 17, ncol = 17)

##The Original 17-Stage Matrix from Zuidema and high-harvest multiplier---------------------------------
plant_S_mat <- matrix( 0, nrow = 17, ncol = 17)
diag(plant_S_mat) <- c(0.455, 0.587, 0.78, 0.821, 0.941, 0.938, 0.961, 0.946, 0.94, 0.937, 0.936, 0.966, 0.968, 0.971, 0.965, 0.967, 0.985)
plant_S_mat[cbind(2:17,1:16)] <- matrix(c(0.091, 0.147, 0.134, 0.167, 0.044, 0.047, 0.034, 0.049, 0.055, 0.058, 0.059, 0.029, 0.027, 0.024, 0.020, 0.018))
plant_S_mat[1,12:17] <- matrix( c(12.3, 14.6, 16.9, 19.3, 22.3, 26.6) )
high_harv[1,12:17] <- highHarvestFecundity 
high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees
plant_mat_low <- plant_S_mat
plant_mat_high <- plant_S_mat * high_harv

###===========================================================================
### Basic Functions
###===========================================================================
sigmoid <- function(k, x0, x) 
{
  1/(1+exp(-k*(x-x0))) #k: steepness #x0 = midpoint
} 


sigmoidMat <- function(k, x0, X) 
{
  1/(1+exp(-k*(X-x0))) #k: steepness #x0 = midpoint
}


markovChain<- function(maxt, i){
  
  statesNames = c("low","high")
  mcHarvest <- new("markovchain", states = statesNames, 
                   transitionMatrix = matrix(data = c(1-i, i, 1-i, i), byrow = TRUE, 
                                             nrow = 2, dimnames=list(statesNames,statesNames)), name="Harvest")
  # Simulating a discrete time process for harvest
  set.seed(100)
  harvest_seq <- markovchain::rmarkovchain(n=maxt, object = mcHarvest, t0="low")
  return(harvest_seq)
}


###===========================================================================
### Sensitivity Matrix
###===========================================================================
sensitivity_matrix <- function(A)
{
  # Computes the Sensitivity and Elasticity Matrices by method specified by Caswell 2001
  #   Chapter 9.
  #
  # Inputs:
  #   A : (n,n) matrix
  #
  # Returns:
  #   A list with 3 items:
  #     sensMat : (n,n) matrix
  #                 i,j values are sensitivity values for each i,j element of A
  #
  #     elasMat : (n,n) matrix
  #                 i,j values are elasticity values for each i,j element of A
  #
  #     stableAge : (n) vector
  #                 The stable Age distribution for A
  
  eigvals <- eigen(A)$values
  eigvecs <- eigen(A)$vectors
  W <- eigvecs
  
  sensMat <- matrix(0, nrow = nrow(A), ncol = ncol(A))
  diag(sensMat) <- c(eigvals)
  
  dominant <- 0
  for (i in eigvals)
  {
    dominant <- max(sqrt(Re(i*Conj(i))), dominant)
  }
  finalCount <- which(eigvals==dominant) # The index of dominant eigenvalue.
  
  V <- Conj(solve(W)) # solve() finds inverse matrix
  w <- matrix(W[,finalCount])
  v <- matrix(Re(V[finalCount,]))
  
  sensMat <- Re(v %*% t(w))
  
  elasMat <- (sensMat * A) / dominant
  
  return(elasMat)
}
###===========================================================================


###===========================================================================
### Elasticity Latin Hypercube Averaging
###===========================================================================
# Parameter multipliers are in following order:
#   1.  Low Hunting --- [0, 1]
#   2.  High Hunting --- [0, 1]
#   3.  High Hunting Proportion --- [0, 1]

test_Hunting <- function(X) # X is LHS matrix (each row is a parameter set multiplier)
{
  plantMat <- plant_S_mat

  numSamples <- nrow(X)

  stochGrowthHiLo <- numeric(numSamples)
  finalPopHiLo <- numeric(numSamples)
  
  stochGrowthConst <- numeric(numSamples)
  finalPopConst <- numeric(numSamples)
  
  for (i in 1:numSamples)
  {
    # Unpacking parameters from matrix
    lowHunting <- X[i, 1]
    highHunting <- X[i, 2]
    huntProp <- X[i, 3] # The (rough) proportion of high hunting
    
    harvest_seq <- markovChain(maxt, huntProp)
    
    constHunt <- highHunting * (huntProp) + lowHunting * (1 - huntProp) # The weighted average of hunting
    #constHunt <- highHunting * table(harvest_seq)[1]/maxt + lowHunting * (1 - table(harvest_seq)[2]/maxt) # A diff (better?) way to calculate
    
    ### Using parameters ###
    # Setting low harvest matrix
    plant_S_mat <- matrix( 0, nrow = 17, ncol = 17)
    diag(plant_S_mat) <- c(0.455, 0.587, 0.78, 0.821, 0.941, 0.938, 0.961, 0.946, 0.94, 0.937, 0.936, 0.966, 0.968, 0.971, 0.965, 0.967, 0.985)
    plant_S_mat[cbind(2:17,1:16)] <- matrix(c(0.091, 0.147, 0.134, 0.167, 0.044, 0.047, 0.034, 0.049, 0.055, 0.058, 0.059, 0.029, 0.027, 0.024, 0.020, 0.018))
    plant_S_mat[1,12:17] <- matrix( c(12.3, 14.6, 16.9, 19.3, 22.3, 26.6) )
    
    # Setting high harvest matrix
    high_harv <- matrix(1, nrow = 17, ncol = 17)    
    high_harv[1,12:17] <- highHarvestFecundity 
    high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees
    
    plant_mat_low <- plant_S_mat
    plant_mat_high <- plant_S_mat * high_harv
    
    agouti_to_PlantSteepness <- -(log(1-m)-log(m))/((m-0.5)*agoutiCapacity) # Steepness needed for sigmoid(m) = m
    plant_to_AgoutiSteepness <- -(log(1-m)-log(m))/((m-0.5)*adultCapacity) 
    
  
      
    nHiLo <- 0 # Pop of plants after time maxt
    rHiLo <- numeric(maxt) # Initializing Array for growth rate values. maxt is global constant
    plant_mat <- matrix(0, nrow = 17)
    plant_mat[1:4] <- seedlingInit/4   # Initial populations are constants
    plant_mat[5:11] <- saplingInit/7   # So they are not parameters
    plant_mat[12:17] <- adultInit/6 
    agouti_vec <- c(agoutiInit) * 0.5

    # Running the simulation to find 'Stochastic Growth Rate' of Hi/Lo Hunting
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
      
      p <- sigmoid(plant_to_AgoutiSteepness, adultCapacity/2, sum(plant_mat[12:17]))*delta + (1-delta)
      agouti_vec[(j+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(j)],agoutiCapacity,h_off, p)
      
      if (agouti_vec[j+1] < 0)
      {
        agouti_vec[j+1] <- 0
      }
      
      plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
      plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(j+1)])*deltaPlant + (1- deltaPlant)
      #  plant_animal_mat[1,12:17] <- linear(m, agouti_vec[(j+1)], b) # A different functional form
      plant_mat <- matrix( c((plant_animal_mat * pmat) %*% plant_mat))
      
      nHiLo <- sum(plant_mat) 
      rHiLo[i] <- log(nHiLo / prevN) # Calculating Growth Rate
    }
    
    
    nConst <- 0 # Pop plants after maxt time
    rConst <- numeric(maxt)
    plant_mat <- matrix(0, nrow = 17)
    plant_mat[1:4] <- seedlingInit/4   # Initial populations are constants
    plant_mat[5:11] <- saplingInit/7   # So they are not parameters
    plant_mat[12:17] <- adultInit/6 
    agouti_vec <- c(agoutiInit) * 0.5
    
    # Running Sim to find Stochastic Growth Rate of CONSTANT Hunting.
    for (j in 1:maxt)
    {
      h_j <- harvest_seq[j]
      
      if (h_j == "low") 
      {
        pmat <- plant_mat_low
        h_off <- constHunt
      } 
      
      else 
      {
        pmat <- plant_mat_high
        h_off <- constHunt
      }
      
      prevN <- sum(plant_mat)
      
      p <- sigmoid(plant_to_AgoutiSteepness, adultCapacity/2, sum(plant_mat[12:17]))*delta + (1-delta)
      agouti_vec[(j+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(j)],agoutiCapacity,h_off, p)
      
      if (agouti_vec[j+1] < 0)
      {
        agouti_vec[j+1] <- 0
      }
      
      plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
      plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(j+1)])*deltaPlant + (1- deltaPlant)
      #  plant_animal_mat[1,12:17] <- linear(m, agouti_vec[(j+1)], b) # A different functional form
      plant_mat <- matrix( c((plant_animal_mat * pmat) %*% plant_mat))
      
      nConst <- sum(plant_mat) 
      rConst[i] <- log(nConst / prevN) # Calculating Growth Rate
    }
    
    stochGrowthHiLo[i] <- exp(mean(rHiLo)) # Collect growth rates into column vector
    finalPopHiLo[i] <- nHiLo # Plant population after time. (not animals)
    
    stochGrowthConst[i] <- exp(mean(rConst))
    finalPopConst[i] <- nConst
    
  }

  diffStoch <- abs(stochGrowthHiLo - stochGrowthConst)
  diffFinalPop <- abs(finalPopHiLo - finalPopConst)

  return(list(stoch=diffStoch, pop=diffFinalPop))
}


###===========================================================================

lhSample <- t(randomLHS(4, 2000)) # Doing the Latin Hypercube Sampling. (Needs to be transposed)
a <- test_Hunting(lhSample)
plot((a$stoch))
summary(a$stoch)

#Smpls2000Yrs8000
