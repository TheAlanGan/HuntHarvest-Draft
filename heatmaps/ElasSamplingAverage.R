# Averaging Elasticities over randomly generated parameters
# using Latin hypercube sampling

library(lhs)

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
#agoutiInit <- 5000


m <- 0.05    # m is the desired proportion at which sigmoid(m) = m . Ideally it is small (~0.01-0.05).
agouti_to_PlantSteepness <- -(log(1-m)-log(m))/((m-0.5)*agoutiCapacity) # Steepness needed for sigmoid(m) = m
plant_to_AgoutiSteepness <- -(log(1-m)-log(m))/((m-0.5)*adultCapacity)  # Steepness needed for sigmoid(m) = m
#This formula above is derived from logistic function with "x = m*CAP" , "x0 = .5*CAP" , "y = m" , and solving for k. (CAP = carrying capacity)

time_end <- 100 # Length of simulation in years
maxt <- 100 #controls the markovChain

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


sigmoid <- function(k, x0, x) 
{
  1/(1+exp(-k*(x-x0))) #k: steepness #x0 = midpoint
} 

sigmoidMat <- function(k, x0, X) 
{
  1/(1+exp(-k*(X-x0))) #k: steepness #x0 = midpoint
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
#   1.  Adult Survival under Harvest --- [0.85, 0.99]
#   2.  Germination under Harvest --- [10, 25] ?
#   3.  Animal proportion of carrying capacity --- [0, 1]
#   4.  Steepness of sigmoid --- [0.03, 0.1]

# Saved for later:
#   5.  Seedling Survival Probabilities --- []
#   6.  Sapling Survival Probabilities  --- []
#   7.  Adult Survival Probabilities --- []
#   8.  Seedling Transition Probabilities --- []
#   9.  Sapling Transition Probabilties --- []
#   10. Adult Transition Probabilities --- []

elas_lhs <- function(X) # X is LHS matrix (each row is a parameter set multiplier)
{
  plantMat <- plant_S_mat
  plantMat <- plantMat
  
  AdultSurvElas <- c()
  SaplingSurvElas <- c()
  SeedlingSurvElas <- c()
  
  GerminationElas <- c()
  
  SeedlingTransElas <- c()
  SaplingTransElas <- c()
  AdultTransElas <- c()
  
  for (i in 1:nrow(X))
  {
    plantMat[cbind(12:17,12:17)] <- X[i, 1] * (.99 - .85) + .85 # Adult Survival... # Mapping [0,1] to [0.85,0.99]
    plantMat[1,12:17] <- X[i, 2] * (25 - 10) + 10 # Adult Germination... # Mapping [0,1] to [10,25]
    m <- X[i, 4] * (0.1-0.03) + 0.03 # Steepness of Sigmoid... # Mapping [0,1] to [0.01,0.05]
    
    agouti_to_PlantSteepness <- -(log(1-m)-log(m))/(m-0.5) # Steepness needed for sigmoid(m) = m

    plantMat[1,12:17] <- sigmoid(agouti_to_PlantSteepness, 1/2, X[i, 3]) # Sigmoid using X[,3]
    
    elas <- sensitivity_matrix(plantMat)

    # Getting the proper elasticity values for each parameter set    
    SeedlingSurvElas[i] <- mean(elas[cbind(1:4,1:4)])    
    SaplingSurvElas[i] <- mean(elas[cbind(5:11,5:11)])
    AdultSurvElas[i] <- mean(elas[cbind(12:17,12:17)])
    
    GerminationElas[i] <- mean(elas[1,12:17])
    
    SeedlingTransElas[i] <- mean(elas[cbind(2:5,1:4)])
    SaplingTransElas[i] <- mean(elas[cbind(6:12,5:11)])
    AdultTransElas[i] <- mean(elas[cbind(13:17,12:16)])
  }
  
  # Averaging elasticity values over ALL parameter sets
  avgSeedlingSurvElas <- mean(SeedlingSurvElas)    
  avgSaplingSurvElas <- mean(SaplingSurvElas)
  avgAdultSurvElas <- mean(AdultSurvElas)
  
  avgGerminationElas <- mean(GerminationElas)
  
  avgSeedlingTransElas <- mean(SeedlingTransElas)
  avgSaplingTransElas <- mean(SaplingTransElas)
  avgAdultTransElas <- mean(AdultTransElas)

  return(c(avgSeedlingSurvElas, avgSaplingSurvElas, avgAdultSurvElas, avgGerminationElas, avgSeedlingTransElas, avgSaplingTransElas, avgAdultTransElas))
}


# This function uses less of a 'for loop'
elas_lhs_better <- function(X)
{
  numSamples <- nrow(X)
  
  plantMat <- array(0, dim = c(numSamples, 17, 17)) # First index is which sample
  
  AdultSurvElas <- matrix(numeric(numSamples))
  SaplingSurvElas <- matrix(numeric(numSamples))
  SeedlingSurvElas <- matrix(numeric(numSamples))
  
  GerminationElas <- matrix(numeric(numSamples))
  
  SeedlingTransElas <- matrix(numeric(numSamples))
  SaplingTransElas <- matrix(numeric(numSamples))
  AdultTransElas <- matrix(numeric(numSamples))
  
  agouti_to_PlantSteepness <- matrix(numeric(numSamples))
  
  
  plantMat[, cbind(12:17, 12:17)] <- X[, 1] * (.99 - .85) + .85 # Adult Survival... # Mapping [0,1] to [0.85,0.99]
  plantMat[, 1, 12:17] <- X[, 2] * (25 - 10) + 10 # Adult Germination... # Mapping [0,1] to [10,25]
  m <- X[, 4] * (0.1-0.03) + 0.03 # Steepness of Sigmoid... # Mapping [0,1] to [0.01,0.05]
    
  agouti_to_PlantSteepness <- -(log(1-m)-log(m))/(m-0.5) # Steepness needed for sigmoid(m) = m
    
  plantMat[, 1, 12:17] <- sigmoidMat(agouti_to_PlantSteepness, 1/2, X[, 3]) # Sigmoid using X[,3]

  elas <- matrix(numeric(numSamples))
  
  for (i in 1:numSamples)
  {
    elas[i] <- sensitivity_matrix(plantMat[i,,])
  }

  # Getting the proper elasticity values for each parameter set
  SeedlingSurvElas <- mean(elas[,cbind(1:4,1:4)])    
  SaplingSurvElas <- mean(elas[,cbind(5:11,5:11)])
  AdultSurvElas <- mean(elas[,cbind(12:17,12:17)])
    
  GerminationElas <- mean(elas[,1,12:17])
    
  SeedlingTransElas <- mean(elas[,cbind(2:5,1:4)])
  SaplingTransElas <- mean(elas[,cbind(6:12,5:11)])
  AdultTransElas <- mean(elas[,cbind(13:17,12:16)])
  
  # Averaging elasticity values over all parameter sets
  avgSeedlingSurvElas <- mean(SeedlingSurvElas)    
  avgSaplingSurvElas <- mean(SaplingSurvElas)
  avgAdultSurvElas <- mean(AdultSurvElas)
  
  avgGerminationElas <- mean(GerminationElas)
  
  avgSeedlingTransElas <- mean(SeedlingTransElas)
  avgSaplingTransElas <- mean(SaplingTransElas)
  avgAdultTransElas <- mean(AdultTransElas)
  
  return(c(avgSeedlingSurvElas, avgSaplingSurvElas, avgAdultSurvElas, avgGerminationElas, avgSeedlingTransElas, avgSaplingTransElas, avgAdultTransElas))
}
###===========================================================================
set.seed(100)
lhSample <- t(randomLHS(4, 100, TRUE))
a <- elas_lhs(lhSample)
plot(a)

