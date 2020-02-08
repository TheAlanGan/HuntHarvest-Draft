####################################################################################################
########## Modeling Hunting and Harvesting Dynamics between Plant-Disperser Pair ###################
## Kevin De Angeli | kevindeangeli@utk.edu
## Eeman Abbasi    | eabbasi@sas.upenn.edu 
## Alan Gan        | a@utk.edu
## Last updated: Nov 4 2019
##==================================================================================================
### Script description: 
# Variance-Based Sensitivity Analysis
#===================================================================================================
###=================================================================================================
### Parameters
###=================================================================================================
#install.packages('sensitivity')
#install.packages('markovchain')

library(sensitivity)
library(boot)
library(lhs)

seedlingCapacity <- 5000 # Carrying Capacities. Tree capacities don't matter as much.
saplingCapacity <- 500
adultCapacity <- 100
agoutiCapacity <- 5200
perc <- 0.9

seedlingInit <- perc * seedlingCapacity#5000 # Initial Populations
saplingInit <- perc * saplingCapacity#500
adultInit <- perc * adultCapacity#100
agoutiInit <- perc * agoutiCapacity#5000

maxt <- 500 # Length of simulation for stoch_growth function.

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
{ # p is how the plant affects carrying capacity of agoutis (from 0 to 1)
  Nnext <- (R*N*(1-N/(K*(p))) + N) * (1-H)
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
  harvest_seq <- markovchain::rmarkovchain(n=tEnd, object = mcHarvest, t0="low")
  return(harvest_seq)
}

###=====================================================================

# Parameter list is as follows:
# 1.  m (determines steepness of sigmoid) ... [0.03, 0.1]
# 2.  delta ... [0, 1]
# 3.  high harvest fecundity ... [0, .9] ... 
# 4.  high harvest adult survival ... [0, .9] ... 
# 5.  hunting rate (herbivores) ... [0.037, 0.57]
# 6.  rMax (herbivore) ... [0.67, 1.1]
# 7.  Adult Plant Carrying Capacity ... [20, 150]
# 8.  delta_plant ... [0, 1]

# Function 'plant_sobol'. To be used with Sobol global variance based sensitivity analysis (fast99)
#
# Arguments:
#  - X : (n x 8) matrix. Columns (8) are each parameter. Rows (n) are different parameter sets. 
#
# Output:
#  - plantPopAfterTime : (n x 1) vector. Each row is the total plant population after time 'maxt' for each parameter set from
#                       input matrix X.
#
plant_sobol <- function(X)
{
  # Function taking in matrix of parameters
  # Returns (nx1) matrix of plant population
  
  numSamples <- nrow(X)
  plantPopAfterTime <- numeric(numSamples)
  
  for (i in 1:numSamples)
  {
    
    # Unpacking parameters from matrix
    m <- X[i, 1] * (0.1 - 0.03) + 0.03                # Mapping [0,1] to [0.03,0.1]
    delta<- X[i, 2]                                   # No need for mapping since delta is in [0,1] already
    highHarvestFecundity <- X[i, 3] * (0.9- 0) + 0    # Mapping to [0.85,1]
    highHarvestSurvival <- X[i, 4] * (0.9 - 0) + 0    # Mapping to [0.9,1]
    constHunt <- X[i, 5]*(0.57- 0.037)+ 0.037         #
    agoutiGrowth <- X[i, 6] * (1.1 - 0.67) + 0.67     # Mapping to [0.67,1.1]
    adultCapacity <- X[i, 7] * (150 - 20) + 20        # Mapping to [50,150]
    deltaPlant <- X[i, 8]                             # No need for mapping
    
    ### Using parameters ###
    # Setting low harvest matrix
    plant_transition_mat <- matrix( 0, nrow = 17, ncol = 17)
    diag(plant_transition_mat) <- c(0.455, 0.587, 0.78, 0.821, 0.941, 0.938, 0.961, 0.946, 0.94, 0.937, 0.936, 0.966, 0.968, 0.971, 0.965, 0.967, 0.985)
    plant_transition_mat[cbind(2:17,1:16)] <- matrix(c(0.091, 0.147, 0.134, 0.167, 0.044, 0.047, 0.034, 0.049, 0.055, 0.058, 0.059, 0.029, 0.027, 0.024, 0.020, 0.018))
    plant_transition_mat[1,12:17] <- matrix( c(12.3, 14.6, 16.9, 19.3, 22.3, 26.6) )
    
    # Setting high harvest matrix
    high_harv <- matrix(1, nrow = 17, ncol = 17)    
    high_harv[1,12:17] <- highHarvestFecundity 
    high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees
    
    plant_t0_low <- plant_transition_mat
    plant_t0_high <- plant_transition_mat * high_harv
    
    agouti_to_PlantSteepness <- -(log(1-m)-log(m))/((m-0.5)*agoutiCapacity) # Steepness needed for sigmoid(m) = m
    plant_to_AgoutiSteepness <- -(log(1-m)-log(m))/((m-0.5)*adultCapacity) 
    
    r <- numeric(maxt) # maxt is global constant
    plant_t0 <- matrix(0, nrow = 17)
    plant_t0[1:4] <- seedlingInit/4   # Initial populations are constants
    plant_t0[5:11] <- saplingInit/7   # So they are not parameters
    plant_t0[12:17] <- adultInit/6 
    agouti_vec <- c(agoutiInit) * 0.5
    
    harvest_seq <- markovChain(maxt)
    
    
    N <- 0 # Pop of plants after time maxt
    
    # Running the simulation to find 'plant abundance'
    for (j in 1:maxt)
    {
      h_j <- harvest_seq[j]
      
      if (h_j == "low") 
      {
        pmat <- plant_t0_low
        h_off <- constHunt
        
      } 
      
      else 
      {
        pmat <- plant_t0_high
        h_off <- constHunt
      }
      
      prevN <- sum(plant_t0)
      
      p <- sigmoid(plant_to_AgoutiSteepness, adultCapacity/2, sum(plant_t0[12:17]))*delta + (1-delta)
      agouti_vec[(j+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(j)],agoutiCapacity,h_off, p)
      
      if (agouti_vec[j+1] < 0)
      {
        agouti_vec[j+1] <- 0
      }
      
      plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
      plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(j+1)])*deltaPlant + (1- deltaPlant)
      plant_t0 <- matrix( c((plant_animal_mat * pmat) %*% plant_t0))
      
      N <- sum(plant_t0) 
      r[i] <- log(N / prevN) # Calculating Growth Rate
    }
    
    plantPopAfterTime[i] <- N # Plant population after time
  }
  return(plantPopAfterTime)
  
}

# Function 'animal_sobol'. To be used with Sobol global variance based sensitivity analysis (fast99)
#
# Arguments:
#  - X : (n x 8) matrix. Columns (8) are each parameter. Rows (n) are different parameter sets. 
#
# Output:
#  - popAfterTime : (n x 1) vector. Each row is the animal population after time 'maxt' for each parameter set from
#                       input matrix X.
#

animal_sobol <- function(X) # X is matrix of parameters. Columns are each parameter. Rows are different parameter sets.
{
  # Function taking in matrix of parameters
  # Returns (nx1) matrix of seed-disperser abundance 
  
  numSamples <- nrow(X)
  popAfterTime <- numeric(numSamples)
  
  for (i in 1:numSamples)
  {
    # Unpacking parameters from matrix
    m <- X[i, 1] * (0.1 - 0.03) + 0.03              # Mapping [0,1] to [0.03,0.1]
    deltaDisp <- X[i, 2]                            # Disperser delta is in [0,1] already
    highHarvestFecundity <- X[i, 3] * (0.9- 0) + 0  # Mapping to [0,0.9]
    highHarvestSurvival <- X[i, 4] * (0.9 - 0) + 0  # Mapping to [0,0.9]
    constHunt <- X[i, 5]*(0.57- 0.037)+ 0.037
    agoutiGrowth <- X[i, 6] * (1.1 - 0.67) + 0.67   # Mapping to [0.67,1.1]
    adultCapacity <- X[i, 7] * (150 - 20) + 20      # Mapping to [20,150]
    deltaPlant <- X[i, 8] 
    
    
    ### Using parameters ###
    # Setting low harvest matrix
    plant_transition_mat <- matrix( 0, nrow = 17, ncol = 17)
    diag(plant_transition_mat) <- c(0.455, 0.587, 0.78, 0.821, 0.941, 0.938, 0.961, 0.946, 0.94, 0.937, 0.936, 0.966, 0.968, 0.971, 0.965, 0.967, 0.985)
    plant_transition_mat[cbind(2:17,1:16)] <- matrix(c(0.091, 0.147, 0.134, 0.167, 0.044, 0.047, 0.034, 0.049, 0.055, 0.058, 0.059, 0.029, 0.027, 0.024, 0.020, 0.018))
    plant_transition_mat[1,12:17] <- matrix( c(12.3, 14.6, 16.9, 19.3, 22.3, 26.6) )
    
    # Setting high harvest matrix
    high_harv <- matrix(1, nrow = 17, ncol = 17)
    high_harv[1,12:17] <- highHarvestFecundity
    high_harv[cbind(12:17,12:17)] <- highHarvestSurvival # Multiplier for survival rate of Adult trees
    
    plant_t0_low <- plant_transition_mat
    plant_t0_high <- plant_transition_mat * high_harv
    
    agouti_to_PlantSteepness <- -(log(1-m)-log(m))/((m-0.5)*agoutiCapacity) # Steepness needed for sigmoid(m) = m
    plant_to_AgoutiSteepness <- -(log(1-m)-log(m))/((m-0.5)*adultCapacity)
    
    r <- numeric(maxt) # maxt is global constant, 
    plant_t0 <- matrix(0, nrow = 17)
    plant_t0[1:4] <- seedlingInit/4   # Initial populations are constants
    plant_t0[5:11] <- saplingInit/7   # So they are not parameters
    plant_t0[12:17] <- adultInit/6
    agouti_vec <- c(agoutiInit) * 0.5
    
    harvest_seq <- markovChain(maxt)
    
    
    N <- 0 # Pop of animals after time maxt
    
    for (j in 1:maxt)
    {
      h_j <- harvest_seq[j]
      
      if (h_j == "low")
      {
        pmat <- plant_t0_low
        h_off <- constHunt
        
      }
      
      else
      {
        pmat <- plant_t0_high
        h_off <- constHunt
      }
      
      
      p <- sigmoid(plant_to_AgoutiSteepness, adultCapacity/2, sum(plant_t0[12:17]))*deltaDisp + (1-deltaDisp)
      agouti_vec[(j+1)] <- LogisticGrowthHunt(agoutiGrowth, agouti_vec[(j)],agoutiCapacity,h_off, p)
      
      if (agouti_vec[j+1] < 0)
      {
        agouti_vec[j+1] <- 0
      }
      
      plant_animal_mat <- matrix(1, nrow = 17, ncol = 17)
      plant_animal_mat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, agoutiCapacity/2, agouti_vec[(j+1)])*deltaPlant + (1- deltaPlant)
      plant_t0 <- matrix( c((plant_animal_mat * pmat) %*% plant_t0))
      
      N <- agouti_vec[j+1]
    }
    
    popAfterTime[i] <- N # disperser population 
    
  }
  return(popAfterTime)
  
}




###=====================================================================

#Fourier Amplitude Sens Test (Saltelli)
# Estimating the Sobol indices using the 'fast99' function from 'sensitivity' package.

fastPlant <- fast99(model = plant_sobol ,factors = 8, n = 10000, q.arg = list(min=0.0, max=1.0))
print(fastPlant)
plot(fastPlant)

fastAnimal <- fast99(model = animal_sobol, factors = 8, n = 10000, q.arg = list(min=0.0, max=1.0))
print(fastAnimal)
plot(fastAnimal)


# Plotting 1st order indices
indicesFirstPlant <- (fastPlant$D1 / fastPlant$V)
indicesFirstAnimal <- (fastAnimal$D1 / fastAnimal$V)

data <- data.frame(indicesFirstPlant, indicesFirstAnimal)
factList <- c('m', expression(delta["p->d"]), expression(G[t]), expression(S[t]), expression(hat(R)), expression(r[max]), expression(K["a"]), expression(delta["d->p"]))
counts <- table(indicesFirstPlant, indicesFirstAnimal)
barplot(t(as.matrix(data)), names.arg = factList, ylim = c(0,0.6), legend.text = c('Plant', 'Animal'), xlab = 'Parameters', ylab = 'First-Order Indices', col = c('yellowgreen', 'brown'), beside = TRUE)
par(bg='white')
dev.copy(png, 'SobolNoBackground.png')
dev.off()


# Plotting total-order indices
indicesTotalPlant <- (1 - fastPlant$Dt / fastPlant$V)
indicesTotalAnimal <- (1 - fastAnimal$Dt / fastAnimal$V)

data <- data.frame(indicesTotalPlant, indicesTotalAnimal)
factList <- c('m', expression(delta["p->d"]), expression(G[t]), expression(S[t]), expression(hat(R)), expression(r[max]), expression(K["a"]), expression(delta["d->p"]))
counts <- table(indicesTotalPlant, indicesTotalAnimal)
barplot(t(as.matrix(data)), names.arg = factList, ylim = c(0,1), legend.text = c('Plant', 'Animal'), xlab = 'Parameters', ylab= "Total-Order Indices", col = c('yellowgreen', 'brown'), beside= TRUE)
par(bg='white')
dev.copy(png, 'SobolNoBackground.png')
dev.off()
