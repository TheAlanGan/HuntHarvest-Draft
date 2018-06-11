
###=======================================================================
### Parameters
###=======================================================================
highHarvestFecundity <- 0.85 # Multiplier for fecundity rate for Adult trees under HIGH harvest
highHarvestSurvival <- 0.9 # Multiplier for survival rate of Adult trees under HIGH harvest
agoutiGrowth <- 1.1 # Growth rate of Agoutis in logistic model

lowHunting <- 0.1 # Percentage of agoutis hunted during LOW hunting
highHunting <- 0.25 # Percentage of agoutis hunted during HIGH hunting

seedlingCapacity <- 5000 # Carrying Capacities
saplingCapacity <- 500
adultCapacity <- 100
agoutiCapacity <- 5200

seedlingInit <- 5000 # Initial Populations
saplingInit <- 500
adultInit <- 100
agoutiInit <- 5200

m <- 0.05    # m is the desired proportion at which sigmoid(m) = m . Ideally it is small (~0.01-0.05).
agouti_to_PlantSteepness <- -(log(1-m)-log(m))/((m-0.5)*agoutiCapacity) # Steepness needed for sigmoid(m) = m
plant_to_AgoutiSteepness <- -(log(1-m)-log(m))/((m-0.5)*adultCapacity)  # Steepness needed for sigmoid(m) = m
#This formula above is derived from logistic function with "x = m*CAP" , "x0 = .5*CAP" , "y = m" , and solving for k. (CAP = carrying capacity)

time_end <- 400 # Length of simulation in years

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




###===========================================================================
### Eigenvalues vs. Population of Agoutis
###===========================================================================
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
###===========================================================================



###===========================================================================
### Eigenvalues vs. Percent Harvest
###===========================================================================
# This is only for harvest relating to the fecundity of the trees.
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

###===========================================================================


#Perturbation Analysis ideas:
#   plot agouti population sequence on top of each other given different harvest rates
#   plot tree population sequence on top of each other give different harvest rates
#   



